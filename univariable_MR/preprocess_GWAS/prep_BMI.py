import pandas as pd

# Load and format the GWAS data 
gwas_file = '/Users/sanjeedahs/Desktop/NEWMR/AF/30061737-GCST006414-EFO_0000275.h.tsv'

try:
    # Read the TSV file, skipping the first line if it's a duplicate header
    df = pd.read_csv(gwas_file, sep='\t', low_memory=False, skiprows=1)
    print("File loaded successfully.")
except FileNotFoundError:
    print("The specified GWAS file was not found.")
    exit()
except pd.errors.ParserError:
    print("Error parsing the GWAS file. Please check the file format.")
    exit()

# Select and rename the relevant columns
try:
    df = df.rename(columns={
        'hm_variant_id': 'SNP',
        'hm_chrom': 'chr',
        'hm_pos': 'bp_loc',
        'hm_beta': 'beta',
        'standard_error': 'se',
        'hm_effect_allele': 'effect_allele',
        'hm_other_allele': 'other_allele',
        'hm_effect_allele_frequency': 'eaf',
        'p_value': 'pval'
    })

    # Keep only the necessary columns
    df = df[['SNP', 'chr', 'bp_loc', 'beta', 'se', 'effect_allele', 'other_allele', 'eaf', 'pval']]
    print("Columns selected and renamed.")
except KeyError as e:
    print(f"Error in selecting or renaming columns: {e}")
    exit()

# Ensure that the alleles columns are strings
df['effect_allele'] = df['effect_allele'].astype(str)
df['other_allele'] = df['other_allele'].astype(str)

# Remove rows with multi-base alleles
df = df[(df['effect_allele'].str.len() == 1) & (df['other_allele'].str.len() == 1)]

# Remove rows with NA values
df = df.dropna(subset=['SNP', 'chr', 'bp_loc', 'beta', 'se', 'effect_allele', 'other_allele', 'eaf', 'pval'])

# Save the filtered data to a new file 
filtered_gwas_file = '/Users/sanjeedahs/Downloads/NEWMR/AF/AF.tsv'
df.to_csv(filtered_gwas_file, sep='\t', index=False)

