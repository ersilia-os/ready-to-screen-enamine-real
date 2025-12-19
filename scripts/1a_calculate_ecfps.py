import numpy as np
import pickle
import os

root = os.path.dirname(os.path.abspath(__file__))
os.chdir(root)
DATA_DIR = os.path.join(root, "..", "data")

CHUNK_SIZE = 10_000_000
CHUNKS_NP = 52
CHUNKS_SAMPLE = 104
CHUNKS_LEAD_LIKE = 838

SPLITS = []
SPLITS.extend([f"Enamine_REAL_NaturalProducts_{str(r).zfill(3)}.tsv.zip" for r in range(CHUNKS_NP)])
SPLITS.extend([f"Enamine_REAL_Sample_{str(r).zfill(3)}.tsv.zip" for r in range(CHUNKS_SAMPLE)])
SPLITS.extend([f"Enamine_REAL_LeadLike_{str(r).zfill(3)}.tsv.zip" for r in range(CHUNKS_LEAD_LIKE)])

print(f"Number of splits: {len(SPLITS)}")

# Store list
pickle.dump(SPLITS, open(os.path.join(DATA_DIR, "splits.pkl"), "wb"))