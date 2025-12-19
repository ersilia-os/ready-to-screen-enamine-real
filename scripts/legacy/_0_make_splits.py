# import pickle
# import os

# root = os.path.dirname(os.path.abspath(__file__))
# os.chdir(root)
# DATA_DIR = os.path.join(root, "..", "data")

# CHUNK_SIZE = 10_000_000
# CHUNKS_NP = 52
# CHUNKS_SAMPLE = 104
# CHUNKS_LEAD_LIKE = 838

# SPLITS = []
# SPLITS.extend(["NaturalProducts", str(i).zfill(3), i * CHUNK_SIZE + 1, i * CHUNK_SIZE + CHUNK_SIZE + 1] for i in range(CHUNKS_NP))
# SPLITS.extend(["Sample", str(i).zfill(3), i * CHUNK_SIZE + 1, i * CHUNK_SIZE + CHUNK_SIZE + 1] for i in range(CHUNKS_SAMPLE))
# SPLITS.extend(["LeadLike", str(i).zfill(3), i * CHUNK_SIZE + 1, i * CHUNK_SIZE + CHUNK_SIZE + 1] for i in range(CHUNKS_LEAD_LIKE))

# # Store list
# pickle.dump(SPLITS, open(os.path.join(DATA_DIR, "splits.pkl"), "wb"))