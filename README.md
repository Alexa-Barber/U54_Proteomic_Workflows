# Full Proteomic Workflow
This method uses a modified verison of the Li-Lab-SJTU/pdxSPA repository that has been adapted for data-independent acquisition (DIA) proteomics and search results from Spectronaut. 

## Step 1: Create the combined fasta file: "Create_Fasta_110725.ipynb
The fasta file is generated from start to finish in the jupyter notebook "Create_Fasta_110725.ipynb". The human and mouse proteomes are downloaded from UniProt, digested in-silico using Trypsin-P/Lys-C rules, and filtered to maintain proteotypic peptides within each species. Proteotypic peptides that are found in both species are labeled as "shared peptides". The peptides are then aggregated by the protein root (Protein_ID) and split into three entries: one for human-unique peptides (_HUMAN), one for mouse-unique peptides (_MOUSE), and one for shared peptides (_SHARED). 

## Step 2: Spectronaut search with combined fasta file
The output from the spectronaut search will need to include the following columns: 
- R.Filename
- PG.ProteinNames
- PG.ProteinGroups
- PG.Genes
- F.ExcludedFromQuantification
- F.PeakArea

## Step 3: Filter and split peptides by sample: "preprocess.sh"
The tab-separated output from Spectronaut is the input for this script. First, peptides that would not have been used for protein quantification in Spectronaut are removed from the dataset by excluding columns where F.ExcludedFromQuantification = True. Second, peptides are separated into individual files by sample (R.FileName) for downstream processing and parallelization. 

## Step 4: Modified S-SPA method: "sspa.py"
This script estimates species compisiton of each sample using the ratio of the sum of human-unqiue peptides to the sum of the mouse-unique peptides. Human proteins are quantified by the sum of their unique human peptides and a corrected sum of the shared peptides. Shared peptides are allocated based on the ratio of the sum of human-unique peptides to the sum of mouse-unique peptides for that specific protein when all three peptide classes are defined. If species unique peptides are missing, the shared peptides are allocated based on the overall ratio of human-unique peptides to mouse-unique peptides. 

## Step 5: Merge expression data across samples: "merge_expression_folder.py"
This step can be run quickly on the command-line and creates a table that lists the number of peptides detected and the adjusted expression intensity for each protein across samples. 

## Downstream Analyses
"Limited_Filtering_Imputation_Script.ipynb" is used to remove proteins with fewer than 2 detected peptide across samples, correct for overall intensity imbalance, and impute missing values for proteins that were detected in 2/3 technical replicates within a sample. 
"MaxSum_Normalization_Clustering.ipynb" includes a stricter filtering in which proteins detected in fewer than 50% of samples are removed from the analysis. The remaining missing values are imputed using K-nearest neighbors imputation to allow for dimentionality reduction analyses such as principle component analysis (PCA) and partial-least square- discriminant analysis (PLS-DA). 
