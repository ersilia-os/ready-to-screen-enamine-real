from google.oauth2.service_account import Credentials
from googleapiclient.http import MediaFileUpload
# from scipy.sparse import csr_matrix, save_npz
from googleapiclient.errors import HttpError
from googleapiclient.http import MediaIoBaseDownload
from google_auth_httplib2 import AuthorizedHttp
from rdkit.Chem import rdFingerprintGenerator
from googleapiclient.discovery import build
from rdkit import Chem
from tqdm import tqdm
import numpy as np
import pandas as pd
import httplib2
import os
import sys
import time
import io

def clip_sparse(vect, nBits=2048):
    """
    Convert a sparse RDKit count fingerprint to a dense int8 NumPy array.

    Parameters
    ----------
    vect : rdkit.DataStructs.cDataStructs.SparseIntVect or similar
        Sparse count fingerprint exposing `GetNonzeroElements()`.
    nBits : int, optional, default=2048
        Length of the dense output vector.

    Returns
    -------
    numpy.ndarray
        Dense fingerprint vector of shape `(nBits,)` with dtype `int8`,
        where values are clipped to the maximum representable int8 value.
    """
    MAX_I8 = 127
    arr = np.zeros(nBits, dtype=np.int8)
    for i, v in vect.GetNonzeroElements().items():
        arr[i] = min(v, MAX_I8)
    return arr

def smiles_to_ecfp(smiles, radius=3, nBits=2048):
    """
    Convert a list of (SMILES, ID) into ECFP count fingerprints.

    Parameters
    ----------
    smiles : list of (str, str)
        List of (SMILES, ID) pairs.
    radius : int, optional, default=3
        Radius of the Morgan fingerprint.
    nBits : int, optional, default=2048
        Length of the fingerprint vector.

    Returns
    -------
    OUTPUT_SMILES : list of [str, str]
        [SMILES, ID] pairs.
    X : numpy.ndarray
        2D array of ECFP count fingerprints (dtype=int8).
    """
     
    mfpgen = rdFingerprintGenerator.GetMorganGenerator(radius=radius, fpSize=nBits)
    OUTPUT_SMILES, X = [], []
    for smi, _id in tqdm(smiles):
        try:
            mol = Chem.MolFromSmiles(smi)
            if mol is not None:
                ecfp = mfpgen.GetCountFingerprint(mol)
                ecfp = clip_sparse(ecfp, nBits=nBits)
                X.append(ecfp)
                OUTPUT_SMILES.append([smi, _id])
        except Exception:
            continue
    assert len(OUTPUT_SMILES) == len(X), "Row mismatch between X and OUTPUT_SMILES"   
    return OUTPUT_SMILES, np.array(X, dtype=np.int8)

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
    http = httplib2.Http(timeout=600)
    authed_http = AuthorizedHttp(creds, http=http)
    drive_service = build("drive", "v3", http=authed_http)
    file_metadata = {"name": os.path.basename(file_path), "parents": [folder_id]}
    media = MediaFileUpload(file_path, resumable=False)
    drive_service.files().create(
        body=file_metadata,
        media_body=media,
        fields="id, webViewLink",
        supportsAllDrives=True).execute()

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
    http = httplib2.Http(timeout=600)
    authed_http = AuthorizedHttp(creds, http=http)
    service = build("drive", "v3", http=authed_http)

    query = f"name='{filename}' and '{folder_id}' in parents and trashed=false"
    results = service.files().list(q=query, fields="files(id)", supportsAllDrives=True, includeItemsFromAllDrives=True).execute()
    files = results.get("files", [])
    if not files:
        raise FileNotFoundError(f"'{filename}' not found in folder {folder_id}.")

    file_id = files[0]["id"]
    request = service.files().get_media(fileId=file_id, supportsAllDrives=True)
    # with open(out_path, "wb") as f:
    #     f.write(request.execute())
    with io.FileIO(out_path, "wb") as fh:
        downloader = MediaIoBaseDownload(fh, request,chunksize=10 * 1024 * 1024)  # 100 MB chunks
        done = False
        retries = 0
        while not done:
            try:
                status, done = downloader.next_chunk()
                if status:
                    sys.stderr.write(f"Download {int(status.progress() * 100)}%\n")
                    sys.stderr.flush()
            except (HttpError, OSError):
                retries += 1
                if retries >= 10:
                    raise
                time.sleep(10)

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

# def save_fingerprints(fingerprints, filename):
#     """
#     Save molecular fingerprints to a compressed sparse NPZ file.

#     Parameters
#     ----------
#     fingerprints : numpy.ndarray
#         2D array of fingerprint vectors to save.
#     filename : str
#         Output file path for the compressed NPZ file.

#     Notes
#     -----
#     Fingerprints are converted to uint8 and stored as a CSR (Compressed Sparse Row) matrix
#     to minimize disk space usage.
#     """
#     sparse_fp = csr_matrix(fingerprints.astype(np.uint8))
#     save_npz(filename, sparse_fp)

# def save_ids(output_smiles, filename):
#     """
#     Save a list of [SMILES, ID] pairs to a compressed CSV (.csv.zst) file.

#     Parameters
#     ----------
#     output_smiles : list of [str, str]
#         List containing SMILES strings and corresponding molecule identifiers.
#     filename : str
#         Output path for the compressed CSV file.

#     Notes
#     -----
#     The file is written using Zstandard compression (level 10) for efficient storage.
#     Rows are streamed directly to the output file to minimize memory usage.
#     """
#     with open(filename, "wb") as f, zstd.ZstdCompressor(level=10).stream_writer(f) as zf:
#         writer = csv.writer(io.TextIOWrapper(zf, encoding="utf-8", newline=""))
#         writer.writerow(["SMILES", "ID"])
#         writer.writerows(output_smiles)
    
# def _remove_existing_in_folder(service, name, folder_id):
#     resp = service.files().list(
#         q=f"name='{name}' and '{folder_id}' in parents and trashed=false",
#         fields="files(id)", supportsAllDrives=True).execute()
#     for f in resp.get("files", []):
#         service.files().update(fileId=f["id"], body={"trashed": True},
#                                supportsAllDrives=True).execute()