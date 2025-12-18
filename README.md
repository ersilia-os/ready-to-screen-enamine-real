# Ready-to-screen Enamine REAL

This repository provides the code to prepare the Ersilia subset of the Enamine REAL database, calculating ECFP6 count fingerprints (radius 3, 2048 bits, RDKit version 2025.9.1) at the IRB cluster and storing them at Ersilia's Google Drive.

# Ersilia subset

The Ersilia subset comprises ~9.9B compounds and includes:

- Enamine REAL Sample 1.04B ()
- Enamine REAL Lead-Like 8.37B ()
- Enamine REAL Natural-like Product 0.52B (517,797,846 compounds)

All three libraries are publicly available and can be downloaded directly from the [Enamine READ Database website](https://enamine.net/compound-collections/real-compounds/real-database-subsets?highlight=WyJlbmFtaW5lIiwicmVhbCIsInN1YnNldHMiXQ==). The original files are extremelly heavy (31, 113 and 7 GB, respectively). For reproducibility and robustness, we downloaded them programatically from the Enamine FTP server with the following commands:

```
wget -c --progress=dot:mega --tries=0 --timeout=60 "https://ftp.enamine.net/download/REAL/2025.02_Enamine_REAL_DB_1B.cxsmiles.bz2"
wget -c --progress=dot:mega --tries=0 --timeout=60 "https://ftp.enamine.net/download/REAL/Enamine_REAL_lead-like_cxsmiles.cxsmiles.bz2"
wget -c --progress=dot:mega --tries=0 --timeout=60 "https://ftp.enamine.net/download/REAL/Enamine_REAL_natural_products_like_cxsmiles.cxsmiles.bz2"
```

# Pipeline






## Usage

Clone this repository and create a Conda environment to install package requirements:

```bash
git clone https://github.com/ersilia-os/ready-to-screen-enamine-real
cd ready-to-screen-enamine-real

conda create -n enamine python=3.12
conda activate enamine
pip install -r requirements.txt
```