import pandas as pd
import pickle
import tqdm
import sys
import os
sys.path.append("/aloy/home/acomajuncosa/Ersilia/ready-to-screen-enamine-real/src")
from src import download, smiles_to_ecfp, save_fingerprints, save_ids, upload

root = os.path.dirname(os.path.abspath(__file__))
os.chdir(root)

# Get file
file = int(sys.argv[1])
file = pickle.load(open(os.path.join("..", "files.pkl"), "rb"))[file]

# Download
sys.stderr.write("Downloading...\n")
download(file, os.path.join("..", "tmp", file), "../service.json", "1P3DSLRpD0GL8th2B78b7OauEXp3I2ZY0")

# Get SMILES
sys.stderr.write("Getting SMILES...\n")
SMILES = pd.read_csv(os.path.join("..", "tmp", file)).values.tolist()

# Calculate ECFP4s
sys.stderr.write("Calculating ECFP4s...\n")
OUTPUT_SMILES, X = smiles_to_ecfp(SMILES)

# Create files
sys.stderr.write("Creating output files...\n")
X_file = os.path.join("..", "tmp", file.replace(".csv.zst", "") + "_X.npz")
SMILES_file = os.path.join("..", "tmp", file.replace(".csv.zst", "") + "_SMILES_IDs.csv.zst")
save_ids(OUTPUT_SMILES, SMILES_file)
del OUTPUT_SMILES
save_fingerprints(X, open(X_file, "wb"))
del X

# Upload files
upload(X_file, "../service.json", folder_id="1bzGf7sVJ3xcZwzTYz3JTGzkasMXebzIN")
upload(SMILES_file, "../service.json", folder_id="1bzGf7sVJ3xcZwzTYz3JTGzkasMXebzIN")

# Remove files
os.remove(X_file)
os.remove(SMILES_file)
os.remove(os.path.join("..", "tmp", file))