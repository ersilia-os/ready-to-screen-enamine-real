from google.oauth2.service_account import Credentials
from googleapiclient.http import MediaFileUpload
from scipy.sparse import csr_matrix, save_npz
from rdkit.Chem import rdFingerprintGenerator
from googleapiclient.discovery import build
import pyarrow.parquet as pq
from rdkit import Chem
from tqdm import tqdm
import pyarrow as pa
import numpy as np
import httplib2
import os

def clip_sparse(vect, nBits=2048):
    """
    Convert a sparse RDKit fingerprint vector to a dense int8 list.

    Parameters
    ----------
    vect : rdkit.DataStructs.cDataStructs.ExplicitBitVect or similar
        Sparse fingerprint vector with nonzero elements.
    nBits : int, optional, default=2048
        Length of the dense output vector.

    Returns
    -------
    list of int
        Dense list representation of the fingerprint, where values are
        clipped to the maximum representable int8 value.
    """
    l = [0] * nBits
    for i, v in vect.GetNonzeroElements().items():
        l[i] = v if v < np.iinfo(np.int8).max else np.iinfo(np.int8).max
    return l

def smiles_to_ecfp(smiles, radius=2, nBits=2048):
    """
    Convert a list of SMILES strings into Extended Connectivity Fingerprints (ECFP).

    Parameters
    ----------
    smiles : list of str
        List of SMILES strings to convert.
    radius : int, optional, default=2
        Radius of the Morgan fingerprint (number of bond hops).
    nBits : int, optional, default=2048
        Length (number of bits) of the fingerprint vector.

    Returns
    -------
    OUTPUT_SMILES : list of str
        List of valid SMILES strings successfully converted.
    X : numpy.ndarray
        2D array of ECFP bit vectors (dtype=int8).

    Notes
    -----
    Invalid or unparsable SMILES are skipped.
    """
    mfpgen = rdFingerprintGenerator.GetMorganGenerator(radius=radius, fpSize=nBits)
    OUTPUT_SMILES, X = [], []
    for smi in tqdm(smiles):
        try:
            mol = Chem.MolFromSmiles(smi)
            ecfp = mfpgen.GetCountFingerprint(mol)
            ecfp = clip_sparse(ecfp)
            X.append(ecfp)
            OUTPUT_SMILES.append(smi)
        except:
            pass  
    return OUTPUT_SMILES, np.array(X, dtype=np.int8)

def save_fingerprints(fingerprints, filename):
    """
    Save molecular fingerprints to a compressed sparse NPZ file.

    Parameters
    ----------
    fingerprints : numpy.ndarray
        2D array of fingerprint vectors to save.
    filename : str
        Output file path for the compressed NPZ file.

    Notes
    -----
    Fingerprints are converted to uint8 and stored as a CSR (Compressed Sparse Row) matrix
    to minimize disk space usage.
    """
    sparse_fp = csr_matrix(fingerprints.astype(np.uint8))
    save_npz(filename, sparse_fp)

def upload(file_path, service_file, folder_id="1KqylmcjiLWX9LLPTr2UnjnGlTCpNvvGA"):
    """
    Upload a local file to a specified Google Drive folder using a service account.

    Parameters
    ----------
    file_path : str
        Path to the local file to upload.
    folder_id : str
        ID of the destination Google Drive folder.
    service_file : str
        Path to the Google Cloud service account JSON credentials file.

    Returns
    -------
    dict
        Dictionary containing the uploaded file’s metadata:
        - 'id': The unique ID of the uploaded file.
        - 'webViewLink': The URL to view the uploaded file in Google Drive.

    Notes
    -----
    This function authenticates with Google Drive using a service account,
    builds the Drive API client, and uploads the file using a resumable upload.
    The HTTP timeout is extended to 600 seconds to support large file transfers.
    """
    creds = Credentials.from_service_account_file(service_file, scopes=["https://www.googleapis.com/auth/drive.file"])
    drive_service = build("drive", "v3", credentials=creds)
    file_metadata = {"name": os.path.basename(file_path), "parents": [folder_id]}
    media = MediaFileUpload(file_path, resumable=True)
    httplib2.DEFAULT_TIMEOUT = 600
    # _remove_existing_in_folder(drive_service, os.path.basename(file_path), folder_id)
    drive_service.files().create(
        body=file_metadata,
        media_body=media,
        fields="id, webViewLink",
        supportsAllDrives=True).execute()

def save_ids(output_smiles, filename):
    """
    Save a list of [ID, SMILES] pairs to a Parquet file.

    Parameters
    ----------
    output_smiles : list of [str, str]
        List containing molecule identifiers and corresponding SMILES strings.
    filename : str
        Output path for the Parquet file.

    Notes
    -----
    The file is written using Zstandard compression (level 19) for efficient storage.
    Both columns are stored as UTF-8 strings.
    """
    table = pa.table({
        "ID": pa.array([i[0] for i in output_smiles], type=pa.large_string()),
        "SMILES": pa.array([i[1] for i in output_smiles], type=pa.large_string())})
    pq.write_table(table, filename, compression="zstd", compression_level=19)
    
# def _remove_existing_in_folder(service, name, folder_id):
#     resp = service.files().list(
#         q=f"name='{name}' and '{folder_id}' in parents and trashed=false",
#         fields="files(id)", supportsAllDrives=True).execute()
#     for f in resp.get("files", []):
#         service.files().update(fileId=f["id"], body={"trashed": True},
#                                supportsAllDrives=True).execute()