# Ready-to-screen Enamine REAL

This repository provides the code to prepare the Ersilia subset of the Enamine REAL database, calculating ECFP6 count fingerprints (radius 3, 2048 bits, RDKit version 2025.9.1) at the IRB cluster and storing them at Ersilia's Google Drive.

# Ersilia subset

The Ersilia subset comprises ~10B compounds and includes:

- Enamine REAL Sample 1.04B (1,035,352,518 compounds)
- Enamine REAL Lead-Like 8.37B (8,371,778,942 compounds)
- Enamine REAL Natural-like Product 0.52B (517,797,846 compounds)

All three libraries are publicly available and can be downloaded directly from the [Enamine READ Database website](https://enamine.net/compound-collections/real-compounds/real-database-subsets?highlight=WyJlbmFtaW5lIiwicmVhbCIsInN1YnNldHMiXQ==). The original files are extremely heavy (31, 113 and 7 GB, respectively). For reproducibility and robustness, we downloaded them programmatically from the Enamine FTP server with the following commands:

```
wget -c --progress=dot:mega --tries=0 --timeout=60 "https://ftp.enamine.net/download/REAL/2025.02_Enamine_REAL_DB_1B.cxsmiles.bz2"
wget -c --progress=dot:mega --tries=0 --timeout=60 "https://ftp.enamine.net/download/REAL/Enamine_REAL_lead-like_cxsmiles.cxsmiles.bz2"
wget -c --progress=dot:mega --tries=0 --timeout=60 "https://ftp.enamine.net/download/REAL/Enamine_REAL_natural_products_like_cxsmiles.cxsmiles.bz2"
```

The downloaded files were stored in their original format in the `data/` directory, which is git-ignored. 

# Pipeline

Code in this repository does the following:

1. Splits all molecules (~10B) into 994 chunks of 10M compounds each, including SMILES and Enamine IDs (TSV format). Generated files are zip-compressed and directly uploaded into Ersilia's Google Drive. Running time is about 24h in local. 

2. This step is run at the IRB cluster. 




## Usage

Clone this repository and create a Conda environment to install package requirements:

```bash
git clone https://github.com/ersilia-os/ready-to-screen-enamine-real
cd ready-to-screen-enamine-real

conda create -n enamine python=3.12
conda activate enamine
pip install -r requirements.txt
```