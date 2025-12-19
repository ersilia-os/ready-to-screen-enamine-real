import zipfile
import time
import os
import bz2
import sys

t = time.time()

# root = os.path.dirname(os.path.abspath(__file__))
root = '.'
os.chdir(root)
sys.path.append(os.path.join(root, "..", "src"))
from src import upload, list_files
tmp_dir = os.path.join(root, "..", "tmp")
data_dir = os.path.join(root, "..", "data")
os.makedirs(tmp_dir, exist_ok=True)

# Define some variables
CHUNK_SIZE = 10_000_000

NAMES = ['NaturalProducts', "Sample", "LeadLike"]

# Mapping name to file
name_to_file = {"NaturalProducts": "Enamine_REAL_natural_products_like_cxsmiles.cxsmiles.bz2",
                "Sample": "2025.02_Enamine_REAL_DB_1B.cxsmiles.bz2",
                "LeadLike": "Enamine_REAL_lead-like_cxsmiles.cxsmiles.bz2"}

# Define variables
FOLDER_ID = "1bWrCvi5FXodxQ2S88nYLecHDjk5Jer8Y"
PATH_TO_SERVICE = os.path.join(data_dir, "service.json")

# For each library
for NAME in NAMES:

    # Parse the original SMILES file
    data_dir = "/aloy/scratch/acomajuncosa/Ersilia/Enamine_libs"
    with bz2.open(os.path.join(data_dir, name_to_file[NAME]), "rt") as f:

        # Read the header
        header = f.readline().rstrip("\n").split("\t")

        # Take 'smiles' and 'id' indexes
        smiles_idx = header.index("smiles")
        id_idx = header.index("id")

        sys.stderr.write(f"Parsing file {name_to_file[NAME]}...\n")
        sys.stderr.flush()

        LINE_COUNT = 0
        CHUNK_COUNT = 0
        out = None

        # For each line
        for line in f:

            if LINE_COUNT % CHUNK_SIZE == 0:
                
                # A file is already open
                if out is not None:

                    # Close file
                    out.close()

                    # Zip file
                    with zipfile.ZipFile(os.path.join(tmp_dir, FILENAME_ZIP), "w", zipfile.ZIP_DEFLATED) as zf:
                        zf.write(os.path.join(tmp_dir, FILENAME), arcname=FILENAME)

                    # Check that file is not uploaded yet
                    files_in_drive = list_files(PATH_TO_SERVICE, FOLDER_ID)
                    if FILENAME_ZIP in files_in_drive:
                        raise FileExistsError(f"{FILENAME_ZIP} already exists in Google Drive folder {FOLDER_ID}")

                    # Upload file
                    upload(os.path.join(tmp_dir, FILENAME_ZIP), PATH_TO_SERVICE, folder_id=FOLDER_ID)
                    time.sleep(5)

                    # Check that file is uploaded
                    files_in_drive = list_files(PATH_TO_SERVICE, FOLDER_ID)
                    if FILENAME_ZIP not in files_in_drive:
                        raise FileExistsError(f"{FILENAME_ZIP} does not exist in Google Drive folder {FOLDER_ID}!")

                    # Remove local files
                    os.remove(os.path.join(tmp_dir, FILENAME))
                    os.remove(os.path.join(tmp_dir, FILENAME_ZIP))
                                    
                # Create new file
                FILENAME = f"Enamine_REAL_{NAME}_{str(CHUNK_COUNT).zfill(3)}.tsv"
                FILENAME_ZIP = FILENAME + ".zip"
                out = open(os.path.join(tmp_dir, FILENAME), "w")
                out.write("smiles\tid\n")
                sys.stderr.write(f"  {NAME} →→→→ chunk {CHUNK_COUNT}\n")
                sys.stderr.flush()
                CHUNK_COUNT += 1

            # Read line
            fields = line.rstrip("\n").split("\t")
            out.write(f"{fields[smiles_idx]}\t{fields[id_idx]}\n")
            LINE_COUNT += 1

        # A file is already open
        if out is not None:

            # Close file
            out.close()

            # Zip file
            with zipfile.ZipFile(os.path.join(tmp_dir, FILENAME_ZIP), "w", zipfile.ZIP_DEFLATED) as zf:
                zf.write(os.path.join(tmp_dir, FILENAME), arcname=FILENAME)

            # Check that file is not uploaded yet
            files_in_drive = list_files(PATH_TO_SERVICE, FOLDER_ID)
            if FILENAME_ZIP in files_in_drive:
                raise FileExistsError(f"{FILENAME_ZIP} already exists in Google Drive folder {FOLDER_ID}")

            # Upload file
            upload(os.path.join(tmp_dir, FILENAME_ZIP), PATH_TO_SERVICE, folder_id=FOLDER_ID)
            time.sleep(5)

            # Check that file is uploaded
            files_in_drive = list_files(PATH_TO_SERVICE, FOLDER_ID)
            if FILENAME_ZIP not in files_in_drive:
                raise FileExistsError(f"{FILENAME_ZIP} does not exist in Google Drive folder {FOLDER_ID}!")

            # Remove local files
            os.remove(os.path.join(tmp_dir, FILENAME))
            os.remove(os.path.join(tmp_dir, FILENAME_ZIP))