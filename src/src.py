from google.oauth2.service_account import Credentials
from googleapiclient.http import MediaFileUpload
from scipy.sparse import csr_matrix, save_npz
from rdkit.Chem import rdFingerprintGenerator
from googleapiclient.discovery import build
from rdkit import Chem
from tqdm import tqdm
import numpy as np
import pandas as pd
import httplib2
import os
import csv
import io
import zstandard as zstd

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
    MAX_I8 = 127
    for i, v in vect.GetNonzeroElements().items():
        l[i] = v if v < MAX_I8 else MAX_I8
    return np.array(l, dtype=np.int8)

def smiles_to_ecfp(smiles, radius=3, nBits=2048):
    """
    Convert a list of SMILES strings into Extended Connectivity Fingerprints (ECFP).

    Parameters
    ----------
    smiles : list of [SMILES, ID], (str, str)
        List of SMILES strings to convert.
    radius : int, optional, default=3
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
    for smi, _id in tqdm(smiles):
        try:
            mol = Chem.MolFromSmiles(smi)
            if mol is not None:
                ecfp = mfpgen.GetCountFingerprint(mol)
                ecfp = clip_sparse(ecfp)
                X.append(ecfp)
                OUTPUT_SMILES.append([smi, _id])
                del ecfp
            del mol
        except:
            pass
    assert len(OUTPUT_SMILES) == len(X), "Row mismatch between X and OUTPUT_SMILES"   
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

def upload(file_path, service_file, folder_id):
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
    Save a list of [SMILES, ID] pairs to a compressed CSV (.csv.zst) file.

    Parameters
    ----------
    output_smiles : list of [str, str]
        List containing SMILES strings and corresponding molecule identifiers.
    filename : str
        Output path for the compressed CSV file.

    Notes
    -----
    The file is written using Zstandard compression (level 10) for efficient storage.
    Rows are streamed directly to the output file to minimize memory usage.
    """
    with open(filename, "wb") as f, zstd.ZstdCompressor(level=10).stream_writer(f) as zf:
        writer = csv.writer(io.TextIOWrapper(zf, encoding="utf-8", newline=""))
        writer.writerow(["SMILES", "ID"])
        writer.writerows(output_smiles)
    
# def _remove_existing_in_folder(service, name, folder_id):
#     resp = service.files().list(
#         q=f"name='{name}' and '{folder_id}' in parents and trashed=false",
#         fields="files(id)", supportsAllDrives=True).execute()
#     for f in resp.get("files", []):
#         service.files().update(fileId=f["id"], body={"trashed": True},
#                                supportsAllDrives=True).execute()

def download(filename, out_path, service_file, folder_id):
    """
    Download a file from a Google Drive folder.

    Parameters
    ----------
    filename : str
        Name of the file to download.
    out_path : str
        Local path to save the file.
    service_file : str
        Path to the service account JSON credentials file.
    folder_id : str
        Google Drive folder ID.
    """
    creds = Credentials.from_service_account_file(service_file, scopes=["https://www.googleapis.com/auth/drive.readonly"])
    service = build("drive", "v3", credentials=creds)
    httplib2.DEFAULT_TIMEOUT = 600

    query = f"name='{filename}' and '{folder_id}' in parents and trashed=false"
    results = service.files().list(q=query, fields="files(id)", supportsAllDrives=True, includeItemsFromAllDrives=True).execute()
    files = results.get("files", [])
    if not files:
        raise FileNotFoundError(f"'{filename}' not found in folder {folder_id}.")

    file_id = files[0]["id"]
    request = service.files().get_media(fileId=file_id, supportsAllDrives=True)
    with open(out_path, "wb") as f:
        f.write(request.execute())


def list_files(service_file, folder_id):
    """
    List files in a Google Drive folder.

    Parameters
    ----------
    service_file : str
        Path to the service account JSON file.
    folder_id : str
        ID of the Google Drive folder.

    Returns
    -------
    set[str]
        Set of file names in the folder.
    """
    creds = Credentials.from_service_account_file(service_file, scopes=["https://www.googleapis.com/auth/drive.readonly"])
    service = build("drive", "v3", credentials=creds)
    query = f"'{folder_id}' in parents and trashed=false"

    names, token = set(), None
    while True:
        results = service.files().list(
            q=query,
            fields="nextPageToken, files(name)",
            pageSize=1000,
            pageToken=token,
            supportsAllDrives=True,
            includeItemsFromAllDrives=True
        ).execute()
        names.update(f["name"] for f in results.get("files", []))
        token = results.get("nextPageToken")
        if not token:
            break
    return names