import pandas as pd
import os
import sys
sys.path.append("../src")
from src import upload

root = '.'

########################
### NATURAL PRODUCTS ###
########################

# print("Splitting REAL Natural Products!")
# PATH_TO_DATA = "/home/acomajuncosa/Downloads/Enamine_REAL_natural_products_like_cxsmiles.cxsmiles.bz2"

# for c, chunk in enumerate(pd.read_csv(PATH_TO_DATA, compression="infer", chunksize=10_000_000, dtype_backend="pyarrow", sep='\t', usecols=["smiles", "id"])):

#     # Set filename
#     filename = f"Enamine_REAL_NaturalProducts_{str(c).zfill(3)}.csv.zst"
    
#     # Save chunk
#     chunk.to_csv(os.path.join(root, "..", "tmp", filename), index=False, compression={"method": "zstd", "level": 10})

#     # Upload to Google Drive and remove file
#     upload(os.path.join(root, "..", "tmp", filename), "../service.json", folder_id="1P3DSLRpD0GL8th2B78b7OauEXp3I2ZY0")
#     os.remove(os.path.join(root, "..", "tmp", filename))

########################
##### REAL SAMPLE ######
########################

# print("Splitting REAL Sample!")
# PATH_TO_DATA = "/home/acomajuncosa/Downloads/2025.02_Enamine_REAL_DB_1B.cxsmiles.bz2"

# for c, chunk in enumerate(pd.read_csv(PATH_TO_DATA, compression="infer", chunksize=10_000_000, dtype_backend="pyarrow", sep='\t', usecols=["smiles", "id"])):

#     # Set filename
#     filename = f"Enamine_REAL_Sample_{str(c).zfill(3)}.csv.zst"
    
#     # Save chunk
#     chunk.to_csv(os.path.join(root, "..", "tmp", filename), index=False, compression={"method": "zstd", "level": 10})

#     # Upload to Google Drive and remove file
#     upload(os.path.join(root, "..", "tmp", filename), "../service.json", folder_id="1P3DSLRpD0GL8th2B78b7OauEXp3I2ZY0")
#     os.remove(os.path.join(root, "..", "tmp", filename))


########################
###### LEAD-LIKE #######
########################

print("Splitting REAL Lead-Like!")
PATH_TO_DATA = "/home/acomajuncosa/Downloads/Enamine_REAL_lead-like_cxsmiles.cxsmiles.bz2"

for c, chunk in enumerate(pd.read_csv(PATH_TO_DATA, compression="infer", chunksize=10_000_000, dtype_backend="pyarrow", sep='\t', usecols=["smiles", "id"])):

    # Set filename
    filename = f"Enamine_REAL_LeadLike_{str(c).zfill(3)}.csv.zst"
    
    # Save chunk
    chunk.to_csv(os.path.join(root, "..", "tmp", filename), index=False, compression={"method": "zstd", "level": 10})

    # Upload to Google Drive and remove file
    upload(os.path.join(root, "..", "tmp", filename), "../service.json", folder_id="1P3DSLRpD0GL8th2B78b7OauEXp3I2ZY0")
    os.remove(os.path.join(root, "..", "tmp", filename))