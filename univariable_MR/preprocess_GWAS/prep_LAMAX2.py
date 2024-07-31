import pandas as pd

# Step 1: Load and format the GWAS data
gwas_file = '/Users/sanjeedahs/Desktop/NEWMR/LAMAX2/rntrn_ilamax.bgen.stats.betastd.tsv'
df = pd.read_csv(gwas_file, sep='\t')

# Select and rename the relevant columns
df = df.rename(columns={
    'SNP': 'SNP',
    'CHR': 'chr',
    'BP': 'bp_loc',
    'BETA': 'beta',
    'SE': 'se',
    'ALLELE1': 'effect_allele',
    'ALLELE0': 'other_allele',
    'A1FREQ': 'eaf',
    'P_LINREG': 'pval',
    'N': 'samplesize'
})

# Keep only the necessary columns
df = df[['SNP', 'chr', 'bp_loc', 'beta', 'se', 'effect_allele', 'other_allele', 'eaf', 'pval', 'samplesize']]

# Step 2: Filter for significant p-values
significant_snps = df[df['pval'] < 1e-5].copy()

# Step 3: Calculate the F-statistic and filter for F-statistic > 10
significant_snps['fstat'] = (significant_snps['beta'] ** 2) / (significant_snps['se'] ** 2)
filtered_snps = significant_snps[significant_snps['fstat'] > 10]

# Step 4: Remove rows with multi-base alleles
filtered_snps = filtered_snps[filtered_snps['effect_allele'].apply(len) == 1]
filtered_snps = filtered_snps[filtered_snps['other_allele'].apply(len) == 1]

# Step 5: Remove rows with NA values
filtered_snps = filtered_snps.dropna(subset=['SNP', 'chr', 'bp_loc', 'beta', 'se', 'effect_allele', 'other_allele', 'eaf', 'pval', 'samplesize', 'fstat'])

# Step 6: Save the filtered data to a new file (optional)
filtered_gwas_file = '/Users/sanjeedahs/Desktop/NEWMR/LAMAX2/LAMAX2.tsv'
filtered_snps.to_csv(filtered_gwas_file, sep='\t', index=False)

# Display the first few rows of the filtered DataFrame
print(filtered_snps.head())