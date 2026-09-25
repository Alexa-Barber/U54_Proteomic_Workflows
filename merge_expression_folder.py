#!/usr/bin/env python3
import os
import sys
import pandas as pd
from glob import glob

def merge_expression_folder(input_folder, outdir,
                            file_pattern="*_SSPA_Human_Protein_Expression.txt",
                            output_file="Nina_Final_Merged_Dataset_Human_Protein_Expression.tsv"):

    # Find all matching files
    input_files = sorted(glob(os.path.join(input_folder, file_pattern)))
    n_files = len(input_files)
    
    if not input_files:

        print(f"[ERROR] No files found in {input_folder} matching pattern {file_pattern}")

        return

    merged_df = None
    sample_ids = []


    for idx, file_path in enumerate(input_files, start=1):
        df = pd.read_csv(file_path, sep="\t", low_memory=False)


        # Clean sample ID from filename
        sample_id = os.path.splitext(os.path.basename(file_path))[0]
        sample_id = sample_id.replace("_SSPA_Human_Protein_Expression", "")
        sample_ids.append(sample_id)

        # Keep necessary columns
        df = df[['Protein_ID', 'Gene', 'Peptide_Count', 'Adjusted_Expression']].copy()

        # Convert counts and intensities to integer
        df['Peptide_Count'] = df['Peptide_Count'].astype(int)
        df['Adjusted_Expression'] = df['Adjusted_Expression'].round(0).astype(int)

        # Rename columns

        df = df.rename(columns={
            'Peptide_Count': f'{sample_id}_Peptides',
            'Adjusted_Expression': f'{sample_id}_Intensity'

        })

        # Merge

        if merged_df is None:
            
            # Keep Gene from the first sample
            merged_df = df

        else:
            # Drop Gene from subsequent samples
            df = df

            merged_df = pd.merge(merged_df, df, on='Protein_ID', how='outer')

        # Print progress

        print(f"[INFO] Processed {idx}/{n_files} samples: {sample_id}")



    # Reorder columns: Protein_ID, Gene, all peptide columns, all intensity columns

    peptide_cols = [f'{sid}_Peptides' for sid in sample_ids]
    intensity_cols = [f'{sid}_Intensity' for sid in sample_ids]
    merged_df = merged_df[['Protein_ID'] + peptide_cols + intensity_cols]

    # Missing values remain as NaN

    # Write output
    os.makedirs(outdir, exist_ok=True)
    out_path = os.path.join(outdir, output_file)
    merged_df.to_csv(out_path, sep="\t", index=False)

    print(f"[INFO] Merged expression table written to: {out_path}")
    print(f"[INFO] Total samples processed: {n_files}")
    print(f"[INFO] Total proteins in merged table: {merged_df.shape[0]}")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python merge_expression_folder.py <input_folder> <output_dir>")
        sys.exit(1)

    input_folder = sys.argv[1]
    outdir = sys.argv[2]

    merge_expression_folder(input_folder, outdir)

