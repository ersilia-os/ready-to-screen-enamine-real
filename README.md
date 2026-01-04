# Ready-to-screen Enamine REAL

This repository provides the code to prepare the Ersilia subset of the Enamine REAL database, calculating ECFP6 count fingerprints (radius 3, 2048 bits, RDKit version 2025.9.1) at the IRB cluster and storing them at Ersilia's Google Drive.

# Ersilia subset

The Ersilia subset comprises ~10B compounds and includes:

- Enamine REAL Sample 1.04B (1,035,352,518 compounds)
- Enamine REAL Lead-Like 8.37B (8,371,778,942 compounds)
- Enamine REAL Natural-like Product 0.52B (517,797,846 compounds)

All three libraries are publicly available and can be downloaded directly from the [Enamine REAL Database website](https://enamine.net/compound-collections/real-compounds/real-database-subsets?highlight=WyJlbmFtaW5lIiwicmVhbCIsInN1YnNldHMiXQ==). The original files are extremely heavy (31, 113 and 7 GB, respectively). For reproducibility and robustness, we downloaded them programmatically from the Enamine FTP server with the following commands:

```
wget -c --progress=dot:mega --tries=0 --timeout=60 "https://ftp.enamine.net/download/REAL/2025.02_Enamine_REAL_DB_1B.cxsmiles.bz2"
wget -c --progress=dot:mega --tries=0 --timeout=60 "https://ftp.enamine.net/download/REAL/Enamine_REAL_lead-like_cxsmiles.cxsmiles.bz2"
wget -c --progress=dot:mega --tries=0 --timeout=60 "https://ftp.enamine.net/download/REAL/Enamine_REAL_natural_products_like_cxsmiles.cxsmiles.bz2"
```

The downloaded files were stored in their original format in the `data/` directory, which is git-ignored. 

# Pipeline

Code in this repository does the following:

1. Splits all molecules (~10B) into 994 chunks of 10M compounds each, including SMILES and Enamine IDs (TSV format). Generated files are zip-compressed and directly uploaded into [Ersilia's Google Drive](https://drive.google.com/drive/folders/1bWrCvi5FXodxQ2S88nYLecHDjk5Jer8Y), under the name `{split_name}.tsv.zip`. Running time is about 24h in local. 

2. This step is run in parallel at the IRB cluster. For each split (chunk of 10M compounds), SMILES are downloaded and processed to calcuate ECFP6s using RDKit. Python functions to download and upload data, calculate ECFP6s and list files on remote repositories are located in `src/src.py`. Two main outputs are generated in this step: a `{split_name}_SMILES_IDs.tsv.zip` (including SMILES and compound Enamine IDs) and a `{split_name}_X.npz` file (including their corresponding ECFP6s). Notice that although the number of rows in both files must be identical, this number is not necessarily the same as the input number of compounds (10M), as invalid or failed compounds will be ommitted and dropped from final results. Running time depends a lot on parallelization capacity: each split needs 48GB of RAM and takes 1h30min, roughly. 

3. Several checks have been implemented to assess the consistency of final results. All checks are found in `ipynb` format in the `notebooks` folder:

- `check_errors.ipynb`: Checking that there are no errors/warnings in job logs.
- `check_files.ipynb`: Checking that files in Google Drive ("libraries" and "ecfps") are exactly the ones expected - no more, no less, no duplicates.
- `check_consistency.ipynb`: Checking that outcome results are consistent. For each split, (i) download both SMILES and ECFP6s, (ii) check that the number of rows between both files are identical, (iii) manually sample 10k rows (both SMILES and ECFP6s) from the huge 10M matrices, (iv) manually calculate ECFP6s for the 10k SMILES and (v) check that these newly calculated 10k ECFP6s are identical to the ones extracted from the huge matrix. 


## Usage

Clone this repository and create a Conda environment to install package requirements:

```bash
git clone https://github.com/ersilia-os/ready-to-screen-enamine-real
cd ready-to-screen-enamine-real

conda create -n enamine python=3.12
conda activate enamine
pip install -r requirements.txt
```
