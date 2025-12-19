import pandas as pd
import pickle
import tqdm
import sys
import os
sys.path.append("/aloy/home/acomajuncosa/Ersilia/ready-to-screen-enamine-real/src")
from src import download, smiles_to_ecfp, save_fingerprints, save_ids, upload, list_files

root = os.path.dirname(os.path.abspath(__file__))
os.chdir(root)

# Get file
file = int(sys.argv[1])
file = pickle.load(open(os.path.join("..", "files.pkl"), "rb"))[file]

# Download
sys.stderr.write(f"Processing file {file}\n")
sys.stderr.write("Downloading...\n")
download(file, os.path.join("..", "tmp", file), "../service.json", "1P3DSLRpD0GL8th2B78b7OauEXp3I2ZY0")

# Get SMILES
sys.stderr.write("Getting SMILES...\n")
SMILES = pd.read_csv(os.path.join("..", "tmp", file)).values.tolist()

# Calculate ECFP4s
sys.stderr.write("Calculating ECFP6s...\n")
OUTPUT_SMILES, X = smiles_to_ecfp(SMILES)
sys.stderr.write(f"Len OUTPUT_SMILES:{len(OUTPUT_SMILES)}...\n")
sys.stderr.write(f"Shape ECFP6s:{X.shape}...\n")

# Create files
sys.stderr.write("Creating output files...\n")
X_file = os.path.join("..", "tmp", file.replace(".csv.zst", "") + "_X.npz")
SMILES_file = os.path.join("..", "tmp", file.replace(".csv.zst", "") + "_SMILES_IDs.csv.zst")
save_ids(OUTPUT_SMILES, SMILES_file)
del OUTPUT_SMILES
save_fingerprints(X, open(X_file, "wb"))
del X

# # Upload files
# upload_file = True
# FOLDER_ID = "19x9yAUySBXgrHBE3gjGjHsomcLQuLQmL"
# while upload_file:
#     upload(X_file, "../service.json", folder_id=FOLDER_ID)
#     upload(SMILES_file, "../service.json", folder_id=FOLDER_ID)
#     existing_files = list_files("../service.json", folder_id=FOLDER_ID)
#     if os.path.basename(X_file) in existing_files and os.path.basename(SMILES_file) in existing_files:
#         upload_file = False
#     else:
#         sys.stderr.write("Upload failed, retrying...\n")
#         sys.stderr.write(f"Existing files: {existing_files}\n"
#                          f"Expected files: {os.path.basename(X_file)}, {os.path.basename(SMILES_file)}\n")

# # Remove files
# os.remove(X_file)
# os.remove(SMILES_file)
# os.remove(os.path.join("..", "tmp", file))