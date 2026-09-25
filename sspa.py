#!/usr/bin/env python3
import os
import sys
import pandas as pd

# -------------------------------
# Step 1: Preprocess the file
# -------------------------------
def filePreTreat(df):
    return df

# -------------------------------
# Step 2: Sort peptides into HU, MU, SH
# -------------------------------
def sortPeptides(protein_name):
    name = protein_name.upper()

    if "HUMAN" in name:
        return "HU"
    elif "MOUSE" in name:
        return "MU"
    elif "SHARED" in name:
        return "SH"
    else:
        return "UNKNOWN"

# -------------------------------
# Step 3: Calculate total R
# -------------------------------
def calculateR(df_HU, df_MU):
    hu_sum = df_HU["F.PeakArea"].fillna(0).sum()
    mu_sum = df_MU["F.PeakArea"].fillna(0).sum()
    return hu_sum / mu_sum if mu_sum != 0 else 0
	
# -------------------------------
# Step 4: Get unique human proteins
# -------------------------------
def getHumanProteins(df):
    # Identify proteins that have at least one human or ambiguous peptide
    protein_set = set()
    for idx, row in df.iterrows():
        if row["flag"] in ["HU", "SH"]:
            for gid in row["PG.ProteinGroups"].split(";"):
                protein_set.add(gid)
    return protein_set

# -------------------------------
# Step 5: S-SPA calculations using human percentage
# -------------------------------
def getHumanProteinExpressionSSPA(df):
    df_expression_all = pd.DataFrame()
    df_ratios_all = pd.DataFrame(columns=["Protein_ID", "Gene", "Protein_Human_Pct", "Ratio_Type"])

    # Calculate overall mouse percentage
    total_HU = df[df["flag"]=="HU"]["F.PeakArea"].fillna(0).sum()
    total_MU = df[df["flag"]=="MU"]["F.PeakArea"].fillna(0).sum()
    mouse_pct_overall = total_MU / (total_HU + total_MU) if (total_HU + total_MU) > 0 else 0
    human_pct_overall = 1 - mouse_pct_overall

    # Determine unique protein IDs based on first field of PG.ProteinNames
    df["Protein_ID"] = df["PG.ProteinNames"].apply(lambda x: x.split("_")[0])
    human_proteins = df.loc[df["flag"].isin(["HU", "SH"]), "Protein_ID"].unique()

    for protein in human_proteins:
        dfp = df[df["Protein_ID"] == protein]

        hu_intensity = dfp[dfp["flag"]=="HU"]["F.PeakArea"].fillna(0).sum()
        mu_intensity = dfp[dfp["flag"]=="MU"]["F.PeakArea"].fillna(0).sum()
        sh_intensity = dfp[dfp["flag"]=="SH"]["F.PeakArea"].fillna(0).sum()
        peptide_count_mu = len(dfp[dfp["flag"]=="MU"])
        peptide_count_hu = len(dfp[dfp["flag"]=="HU"])
        peptide_count_sh = len(dfp[dfp["flag"]=="SH"])
        peptide_count = peptide_count_hu + peptide_count_sh

        # Determine human percentage for this protein
        if (hu_intensity > 0 &  mu_intensity > 0 & sh_intensity > 0:
            human_pct_protein = hu_intensity / (hu_intensity + mu_intensity)
            ratio_type = "Protein-specific"
        else:
            human_pct_protein = human_pct_overall
            ratio_type = "Overall"

        # Adjusted protein expression
        adjusted_expression = hu_intensity + sh_intensity * human_pct_protein

        # Gene name (take first non-null gene)
        gene = dfp["PG.Genes"].dropna().iloc[0] if len(dfp["PG.Genes"].dropna()) > 0 else ""

        # Append to expression DataFrame
        df_expression_all = pd.concat([
            df_expression_all,
            pd.DataFrame({
                "Protein_ID": [protein],
                "PG.ProteinGroups": [";".join(dfp["PG.ProteinGroups"].unique())],
                "Gene": [gene.upper()],
                "Peptide_Count": [peptide_count],
                "Adjusted_Expression": [adjusted_expression]
            })
        ], ignore_index=True)

        # Append to ratios DataFrame
        df_ratios_all = pd.concat([
            df_ratios_all,
            pd.DataFrame({
                "Protein_ID": [protein],
                "PG.ProteinGroups": [";".join(dfp["PG.ProteinGroups"].unique())],
                "Gene": [gene],
                "Protein_Human_Pct": [human_pct_protein],
                "Ratio_Type": [ratio_type]
            })
        ], ignore_index=True)

    return df_expression_all, df_ratios_all
	
# -------------------------------
# Step 6: Pipeline wrapper (modified)
# -------------------------------
def S_SPA_pipeline(input_file, outdir):
    # Read input
    df = pd.read_csv(input_file, sep="\t", low_memory=False)
    df = filePreTreat(df)
    df["flag"] = df["PG.ProteinNames"].apply(sortPeptides)

    os.makedirs(outdir, exist_ok=True)

    # S-SPA calculations
    df_expression, df_ratios = getHumanProteinExpressionSSPA(df)

    # Extract sample name from input file (without path and extension)
    sample_name = os.path.splitext(os.path.basename(input_file))[0]

    # Write outputs with sample-specific filenames
    df_expression.to_csv(os.path.join(outdir, f"{sample_name}_SSPA_Human_Protein_Expression.txt"),
                         sep="\t", index=False)
    df_ratios.to_csv(os.path.join(outdir, f"{sample_name}_SSPA_Human_Mouse_Ratios.txt"),
                     sep="\t", index=False)


# -------------------------------
# Main (modified for array job)
# -------------------------------
if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python sspa.py <input_file> <output_dir>")
        sys.exit(1)

    input_file = sys.argv[1]
    outdir = sys.argv[2]

    S_SPA_pipeline(input_file, outdir)

    n_total = df_expression_all["Protein_ID"].nunique()

    sample_name = os.path.splitext(os.path.basename(input_file))[0]
    summary_path = os.path.join(outdir, "Overall_Human_Percentage_Table.tsv")
    lock_path = summary_path + ".lock"

    # ---------------------------
    # Append this sample’s summary to the shared table
    # ---------------------------
    new_data = pd.DataFrame({
        "Sample_Name": [sample_name],
        "Total_Proteins": [n_total],
        "HU_Intensity": [round(hu_intensity)],
        "MU_Intensity": [round(mu_intensity)],
        "Human_Percentage": [round(human_pct_overall * 100, 2)]
    })

    os.makedirs(outdir, exist_ok=True)

    with FileLock(lock_path, timeout=300):
        if os.path.exists(summary_path):
            existing = pd.read_csv(summary_path, sep="\t")
            # Drop any old record for this sample (if rerunning)
            existing = existing[existing["Sample_Name"] != sample_name]
            combined = pd.concat([existing, new_data], ignore_index=True)
        else:
            combined = new_data
        combined.to_csv(summary_path, sep="\t", index=False)

    print(f"\n[INFO] {sample_name} Using overall human percentage: "
          f"{human_pct_overall*100:.2f}% ({hu_intensity:.0f} HU / {mu_intensity:.0f} MU, "
          f"{n_total} total proteins)")
