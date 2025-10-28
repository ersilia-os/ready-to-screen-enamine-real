from google.oauth2.service_account import Credentials
from googleapiclient.http import MediaFileUpload
from scipy.sparse import csr_matrix, save_npz
from rdkit.Chem import rdFingerprintGenerator
from googleapiclient.discovery import build
from rdkit import Chem
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
    for smi in smiles:
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

def upload(file_path, folder_id, service_file):
    creds = Credentials.from_service_account_file(service_file, scopes=["https://www.googleapis.com/auth/drive.file"])
    drive_service = build("drive", "v3", credentials=creds)
    file_metadata = {"name": os.path.basename(file_path), "parents": [folder_id]}
    media = MediaFileUpload(file_path, resumable=True)
    httplib2.DEFAULT_TIMEOUT = 600
    uploaded = drive_service.files().create(
        body=file_metadata,
        media_body=media,
        fields="id, webViewLink",
        supportsAllDrives=True
    ).execute()

