#!/bin/bash

# Initialize conda
eval "$(conda shell.bash hook)"
conda activate /aloy/home/acomajuncosa/anaconda3/envs/enamine

# Get alpha from first argument
alpha=$1

# Determine the directory where run.sh is located
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"

# Call main.py using an absolute path
python "$SCRIPT_DIR/1b_calculate_ecfps.py" "$alpha"