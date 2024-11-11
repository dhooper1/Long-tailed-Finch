import pandas as pd
import sys

# Check if input files were provided
if len(sys.argv) != 5:
    print("Usage: python genotype_discordance.py truth_genotypes.tsv.gz imputed_genotypes.tsv.gz <INFO_SCORE_threshold> <min_AF>")
    sys.exit(1)

# Assign input file paths from command line arguments
truth_file = sys.argv[1]
imputed_file = sys.argv[2]
info_threshold = float(sys.argv[3])
min_AF = float(sys.argv[4])

# Load the genotype data
truth_df = pd.read_csv(truth_file, sep="\t", names=["CHROM", "POS", "GT"], compression='gzip')
imputed_df = pd.read_csv(imputed_file, sep="\t", names=["CHROM", "POS", "INFO_SCORE", "AF", "GT", "DS"], compression='gzip')

# Filter imputed data based on INFO_SCORE and AF thresholds
imputed_df = imputed_df[(imputed_df['INFO_SCORE'] >= info_threshold) & (imputed_df['AF'] >= min_AF)]

# Merge data on chromosome and position
merged_df = pd.merge(truth_df, imputed_df, on=["CHROM", "POS"], suffixes=('_truth', '_imputed'))

# Function to convert genotype to a numeric code (e.g., 0, 1, or 2 for biallelic genotypes)
def genotype_to_numeric(gt):
    if gt == "0/0": return 0
    elif gt == "0/1" or gt == "1/0": return 1
    elif gt == "1/1": return 2
    else: return None

# Apply conversion to both genotype columns
merged_df['GT_truth'] = merged_df['GT_truth'].apply(genotype_to_numeric)
merged_df['GT_imputed'] = merged_df['GT_imputed'].apply(genotype_to_numeric)

# Drop any rows with missing values
merged_df.dropna(inplace=True)

# Calculate discordance across all sites
total_sites = len(merged_df)
discordant_all = (merged_df['GT_truth'] != merged_df['GT_imputed']).sum()
discordance_all = round((discordant_all / total_sites) * 100, 3)

# Calculate discordance for each genotype category
categories = {
    "all_sites": merged_df,
    "homozygous_major": merged_df[merged_df['GT_truth'] == 0],
    "heterozygous": merged_df[merged_df['GT_truth'] == 1],
    "homozygous_minor": merged_df[merged_df['GT_truth'] == 2],
}

# Calculate and print discordance for each category
print("Genotype Discordance Rates (%)")
print("--------------------------------")
for category, df in categories.items():
    if not df.empty:
        discordant_count = (df['GT_truth'] != df['GT_imputed']).sum()
        discordance_rate = round((discordant_count / len(df)) * 100, 3)
        print(f"{category.replace('_', ' ').title()}: {discordance_rate}% (out of {len(df)} sites)")
    else:
        print(f"{category.replace('_', ' ').title()}: No sites in this category.")

