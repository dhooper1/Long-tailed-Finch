import pandas as pd
from scipy.stats import pearsonr
import sys
import numpy as np

# Check if input files were provided
if len(sys.argv) != 5:
    print("Usage: python imputation_accuracy_DS_write2out.py truth_genotypes.tsv.gz imputed_genotypes.tsv.gz <INFO_SCORE_threshold> <output.tsv>")
    sys.exit(1)

# Assign input file paths from command line arguments
truth_file = sys.argv[1]
imputed_file = sys.argv[2]
info_threshold = float(sys.argv[3])
output_file = sys.argv[4]

# Load the genotype data, including allele frequency from the imputed file
truth_df = pd.read_csv(truth_file, sep="\t", names=["CHROM", "POS", "GT"], compression='gzip')
imputed_df = pd.read_csv(imputed_file, sep="\t", names=["CHROM", "POS", "INFO_SCORE", "AF", "GT", "DS"], compression='gzip')

# Fold AF to range 0.0 to 0.5
imputed_df['AF'] = imputed_df['AF'].apply(lambda x: min(x, 1 - x))

# Filter imputed data based on INFO_SCORE threshold
imputed_df = imputed_df[(imputed_df['INFO_SCORE'] >= info_threshold)]

# Merge data on chromosome and position
merged_df = pd.merge(truth_df, imputed_df, on=["CHROM", "POS"], suffixes=('_truth', '_imputed'))

# Function to convert genotype to numeric encoding (e.g., 0, 1, or 2 for biallelic)
def genotype_to_numeric(gt):
    if gt == "0/0": return 0
    elif gt == "0/1" or gt == "1/0": return 1
    elif gt == "1/1": return 2
    else: return None

# Apply conversion to both genotype columns
merged_df['GT_truth'] = merged_df['GT_truth'].apply(genotype_to_numeric)
#merged_df['GT_imputed'] = merged_df['GT_imputed'].apply(genotype_to_numeric)

# Drop any rows with missing values
merged_df.dropna(inplace=True)

# Output the total number of shared positions/sites
num_shared_sites = len(merged_df)
print("Number of shared positions/sites:", num_shared_sites)

# Define custom allele frequency bins
bins = [0.0, 0.001, 0.002, 0.005, 0.01, 0.05, 0.1, 0.2, 0.5]

# List to store results for each bin
results = []

# Calculate and print r^2 for each allele frequency bin
for i in range(len(bins) - 1):
    # Define current bin range
    bin_start = bins[i]
    bin_end = bins[i + 1]

    # Filter for rows where allele frequency is within the current bin
    bin_df = merged_df[(merged_df['AF'] >= bin_start) & (merged_df['AF'] < bin_end)]
    num_sites_in_bin = len(bin_df)

    # Check for sufficient variance in truth and imputed data
    if num_sites_in_bin > 0:
        if bin_df['GT_truth'].nunique() > 1 and bin_df['DS'].nunique() > 1:
            # Calculate r^2 for the bin
            r, _ = pearsonr(bin_df['GT_truth'], bin_df['DS'])
            r_squared = round(r ** 2, 3)
        else:
            r_squared = None  # Indicate insufficient variance
    else:
        r_squared = None  # Indicate no sites in bin

    # Append results for this bin
    results.append({"bin_start": bin_start, "bin_end": bin_end, "r2": r_squared, "num_sites": num_sites_in_bin})
    # Print each result to the screen
    print(f"AF bin: {bin_start}-{bin_end}, r^2: {r_squared}, number of sites: {num_sites_in_bin}")

# Calculate r^2 for all sites
if merged_df['GT_truth'].nunique() > 1 and merged_df['DS'].nunique() > 1:
    r_all, _ = pearsonr(merged_df['GT_truth'], merged_df['DS'])
    r_squared_all = round(r_all ** 2, 3)
else:
    r_squared_all = None

# Add the result for all sites to the results list and print it to the screen
results.append({"bin_start": 0.0, "bin_end": 1.0, "r2": r_squared_all, "num_sites": num_shared_sites})
print(f"AF bin: 0.0-1.0, r^2: {r_squared_all}, number of sites: {num_shared_sites}")

# Convert results to a DataFrame and save as TSV
results_df = pd.DataFrame(results)
results_df.to_csv(output_file, sep="\t", index=False)

print(f"Results saved to {output_file}")

