import pandas as pd
from pyliftover import LiftOver

# Initialize the LiftOver object for GRCh37 to GRCh38 conversion
lo = LiftOver('hg19', 'hg38')

# Load your GWAS data
input_file = '/Users/sanjeedahs/Desktop/NEWMR/LAMAX2/LAMAX2.tsv'
output_file = '/Users/sanjeedahs/Desktop/NEWMR/LAMAX2/LAMAX2_38.tsv'

gwas_data = pd.read_csv(input_file, sep='\t')

# Function to perform the liftover and handle missing conversions
def lift_over_coordinates(chromosome, bp_loc):
    lifted = lo.convert_coordinate('chr' + str(chromosome), bp_loc)
    if lifted:
        return lifted[0][1]
    else:
        return None

# Apply the liftover to the bp_loc column and replace the original column
gwas_data['bp_loc'] = gwas_data.apply(lambda row: lift_over_coordinates(row['chr'], row['bp_loc']), axis=1)

# Save the updated GWAS data to a new file
gwas_data.to_csv(output_file, sep='\t', index=False)

print(f"Lifted over GWAS data saved to {output_file}")
