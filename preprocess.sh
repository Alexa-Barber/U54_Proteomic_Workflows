#!/bin/bash
#SBATCH --job-name=preprocess
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=24G
#SBATCH --time=02:00:00
#SBATCH --output=preprocess.%A.out
#SBATCH --error=preprocess.%A.err

# Change to working directory
cd working_directory

# Load necessary modules
module load anaconda
module load python/3.9

# Specify input file and output directory
INPUT="Spectronaut_Input.tsv"
OUTDIR="preprocessed_files"
mkdir -p "$OUTDIR"

# Step 1: remove rows with F.ExcludedFromQuantification == True
# Step 2: split by R.FileName
awk -F'\t' '
NR==1 {
    # save header line
    header=$0
    # find which column is R.FileName and which is F.ExcludedFromQuantification
    for (i=1; i<=NF; i++) {
        if ($i == "R.FileName") filecol=i
        if ($i == "F.ExcludedFromQuantification") excludecol=i
    }
    next
}
{
    # only keep if not excluded
    if ($excludecol != "True") {
        file=$filecol
        out=file".tsv"
        if (!(file in seen)) {
            print header > "'$OUTDIR'/"out
            seen[file]=1
        }
        print $0 >> "'$OUTDIR'/"out
    }
}' "$INPUT"
