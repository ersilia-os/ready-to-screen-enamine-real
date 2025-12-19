import pandas as pd
import numpy as np
import zipfile
import pickle
import time
import sys
import os

sys.path.append("/aloy/home/acomajuncosa/Ersilia/ready-to-screen-enamine-real/src")
from src import download, smiles_to_ecfp, upload, list_files

root = os.path.dirname(os.path.abspath(__file__))
os.chdir(root)

# Define variables
data_dir = os.path.join(root, "..", "data")
tmp_dir = os.path.join(root, "..", "tmp")
FOLDER_ID_LIBRARY = "1bWrCvi5FXodxQ2S88nYLecHDjk5Jer8Y"
FOLDER_ID_ECFP = "1FBELagBf9hlKVgvkaZ8YF60jKRAmsHPo"
PATH_TO_SERVICE = os.path.join(data_dir, "service.json")

# Get file
IND = int(sys.argv[1])
FILENAME_ZIP = pickle.load(open(os.path.join(data_dir, "splits.pkl"), "rb"))[IND]
FILENAME = FILENAME_ZIP.replace('.zip', '')

# Download
sys.stderr.write(f"Processing file {FILENAME_ZIP}\n")
sys.stderr.write("Downloading...\n")
download(FILENAME_ZIP, os.path.join(tmp_dir, FILENAME_ZIP), PATH_TO_SERVICE, FOLDER_ID_LIBRARY)

# Extract files from ZIP
with zipfile.ZipFile(os.path.join(tmp_dir, FILENAME_ZIP), "r") as zf:
    zf.extractall(tmp_dir)

# Get SMILES
sys.stderr.write("Getting SMILES...\n")
SMILES = pd.read_csv(os.path.join(tmp_dir, FILENAME), sep='\t').values.tolist()
sys.stderr.write(f"Number of original compounds: {len(SMILES)}\n")
sys.stderr.flush()

# Calculate ECFP4s
sys.stderr.write("Calculating ECFP6s...\n")
OUTPUT_SMILES, X = smiles_to_ecfp(SMILES)
sys.stderr.write(f"Len OUTPUT_SMILES:{len(OUTPUT_SMILES)}...\n")
sys.stderr.write(f"Shape X:{X.shape}...\n")
sys.stderr.flush()
assert len(OUTPUT_SMILES) == len(X), "Row mismatch between X and OUTPUT_SMILES"

# Create files
sys.stderr.write("Creating output files...\n")
SMILES_FILE = FILENAME.replace(".tsv", "_SMILES_IDs.tsv")
X_FILE = FILENAME.replace(".tsv", "_X.npz")
# SMILES and IDs
sys.stderr.write("SMILES and IDs...\n")
with open(os.path.join(tmp_dir, SMILES_FILE), "w") as out:
    out.write("smiles\tid\n")
    for smi, _id in OUTPUT_SMILES:
        out.write(f"{smi}\t{_id}\n")
# ECPF6s
sys.stderr.write("ECFP6s...\n")
np.savez_compressed (os.path.join(tmp_dir, X_FILE), X=X)

# Compress SMILES and IDs file
sys.stderr.write("Compressing SMILES and IDs...\n")
sys.stderr.flush()
with zipfile.ZipFile(os.path.join(tmp_dir, SMILES_FILE.replace('.tsv', '.tsv.zip')), "w", zipfile.ZIP_DEFLATED) as zf:
    zf.write(os.path.join(tmp_dir, SMILES_FILE), arcname=SMILES_FILE)

# Check that files are not uploaded yet
sys.stderr.write("Checking Google Drive before uploading...\n")
files_in_drive = list_files(PATH_TO_SERVICE, FOLDER_ID_ECFP)
if SMILES_FILE.replace('.tsv', '.tsv.zip') in files_in_drive:
    raise FileExistsError(f"{SMILES_FILE.replace('.tsv', '.tsv.zip')} already exists in Google Drive folder {FOLDER_ID_ECFP}")
if X_FILE in files_in_drive:
    raise FileExistsError(f"{X_FILE} already exists in Google Drive folder {FOLDER_ID_ECFP}")

# Upload files
sys.stderr.write("Uploading files...\n")
upload(os.path.join(tmp_dir, SMILES_FILE.replace('.tsv', '.tsv.zip')), PATH_TO_SERVICE, FOLDER_ID_ECFP)
upload(os.path.join(tmp_dir, X_FILE), PATH_TO_SERVICE, FOLDER_ID_ECFP)
time.sleep(10)

# Check that files are uploaded already
sys.stderr.write("Checking Google Drive after uploading...\n")
sys.stderr.flush()
files_in_drive = list_files(PATH_TO_SERVICE, FOLDER_ID_ECFP)
if SMILES_FILE.replace('.tsv', '.tsv.zip') not in files_in_drive:
    raise FileExistsError(f"{SMILES_FILE.replace('.tsv', '.tsv.zip')} does not exist in Google Drive folder {FOLDER_ID_ECFP}")
if X_FILE not in files_in_drive:
    raise FileExistsError(f"{X_FILE} does not exist in Google Drive folder {FOLDER_ID_ECFP}")

# Remove local files
sys.stderr.write("Removing local files...\n")
sys.stderr.flush()
os.remove(os.path.join(tmp_dir, FILENAME))
os.remove(os.path.join(tmp_dir, FILENAME_ZIP))
os.remove(os.path.join(tmp_dir, X_FILE))
os.remove(os.path.join(tmp_dir, SMILES_FILE))
os.remove(os.path.join(tmp_dir, SMILES_FILE.replace('.tsv', '.tsv.zip')))


sys.stderr.write("Job finished!\n")
sys.stderr.flush()