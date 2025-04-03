#!/Users/Daniel/opt/anaconda3/envs/spyder_env/bin/python
import pandas as pd
import pyranges as pr
import csv
from collections import defaultdict
import matplotlib.pyplot as plt
import numpy as np
from tqdm import tqdm

def process_provirus(provirus_data, exons_ranges):
    try:
        provirus_id, fragments = provirus_data
        fragment_df = pd.DataFrame(fragments, columns=['Chromosome', 'Start', 'End'])
        provirus_range = pr.PyRanges(fragment_df)
        overlap = provirus_range.intersect(exons_ranges)
        return overlap.df
    except Exception as e:
        print(f"Error processing provirus {provirus_id}: {e}")
        

#Reading provirus key
prv_key = pd.read_csv('/Users/Daniel/Documents/Med_School/G1-G4/Coffin_Lab/Thesis/Thesis_Project/RNA_Seq_Processing/TCGA_Analysis/Annotations/083123_hg38_p14_HERVs_HS_LS_CytoBand.tsv',sep='\t', header=None)

# Create a function to parse the gene_id
def parse_gene_id(gene_id_str):
    return gene_id_str.split('"')[1]

# Apply this function to the 8th column of the DataFrame
prv_key[8] = prv_key[8].apply(parse_gene_id)

# Create a dictionary where the keys are the provirus identifiers and the values are lists of fragments
print("GENERATING PROVIRUS IDENTIFIER DICTIONARY")
provirus_dict = prv_key.groupby(8)[[0, 3, 4]].apply(lambda x: list(map(tuple, x.values))).to_dict()

# Reading expression file
print("READING EXPRESSION FILE")
expr_file = '/Users/Daniel/Documents/Med_School/G1-G4/Coffin_Lab/Thesis/Thesis_Project/RNA_Seq_Processing/TCGA_Analysis/top_proviruses/090623_results_GBM_vs_Controls_top_proviruses_withSense.csv'
df_expr = pd.read_csv(expr_file, index_col=0)


# Reading GENCODE GTF
print("READING GENCODE GTF")
gencode_file = '/Users/Daniel/Documents/Med_School/G1-G4/Coffin_Lab/Thesis/Thesis_Project/RNA_Seq_Processing/TCGA_Analysis/Annotations/gencode.v44.annotation.gtf'
gencode_ranges = pr.read_gtf(gencode_file)

#Get rid of nonsense-mediated-decay products
gencode_ranges = gencode_ranges[gencode_ranges.transcript_type != 'nonsense_mediated_decay']

# Extract exon intervals from the GENCODE data
exons_ranges = gencode_ranges[gencode_ranges.Feature == 'exon']

# Maintain a dictionary that maps each transcript ID to its highest exon number

print("GENERATING MAX EXONS")
transcript_max_exon_num = exons_ranges.df['exon_number'].astype(int).groupby(exons_ranges.df['transcript_id']).max()

print("GENERATING DataFrame TO HOLD ALL PROVIRUSES")
# First, create a DataFrame with all proviruses
provirus_df = pd.concat((pd.DataFrame(fragments, columns=['Chromosome', 'Start', 'End']).assign(Provirus=provirus_id) for provirus_id, fragments in provirus_dict.items()), ignore_index=True)

print("Generating PyRanges object")
provirus_range = pr.PyRanges(provirus_df)

# Extract transcript intervals from GENCODE data 
transcript_ranges = gencode_ranges[gencode_ranges.Feature == 'transcript']

# Then intersect/join
print("RUNNING INTERSECTION WITH ALL FEATURES")
overlaps = provirus_range.join(gencode_ranges)
print("DONE")


# Convert PyRanges to DataFrames
df_overlaps = overlaps.df

# Reset the index of df_expr to make it a column for the merge
df_expr.reset_index(inplace=True)

# Rename the index column to match the 'Provirus' column in df_overlaps
df_expr.rename(columns={'index': 'Provirus'}, inplace=True)

# Merge the dataframes
merged_df = pd.merge(df_expr, df_overlaps, on='Provirus')


gene_ids = set()
gene_orientation = defaultdict(list)
orientation_count = defaultdict(lambda: defaultdict(int))  # For keeping track of the different combinations

# Specify the categories of transcripts to include
included_genetypes = {'lncRNA', 'protein_coding'}


# Preprocessing steps
merged_df['gene_id'] = merged_df['gene_id'].str.split('.').str[0]  # Remove .# suffix from column 17
gene_orientation_dict = dict(zip(merged_df['gene_id'], merged_df['Strand']))
#merged_df['gene_orientation'] = merged_df['gene_id'].map(gene_orientation_dict)

merged_df['exon_number'] = merged_df['exon_number'].astype('Int64')
merged_df['transcript_id'] = merged_df['transcript_id'].astype(str)


# Creating provirus_info DataFrame
provirus_info_df = merged_df.iloc[:, 0:6]
provirus_info_dict = provirus_info_df.groupby('Provirus').apply(lambda x: x.values.tolist()).to_dict()
Proviruses = set(merged_df['Provirus'])

# Determining exon positions
merged_df['is_intronic'] = merged_df['exon_number'].isna()
merged_df['exon_number'] = merged_df['exon_number'].fillna(False)
merged_df['is_first_exon'] = merged_df['exon_number'] == 1
merged_df['is_last_exon'] = merged_df['exon_number'] == merged_df['transcript_id'].map(transcript_max_exon_num)
merged_df['exon_pos'] = 'middle_exon'

merged_df.loc[merged_df['is_intronic'], 'exon_pos'] = 'intron'
merged_df.loc[merged_df['is_first_exon'] & merged_df['is_last_exon'] & ~merged_df['is_intronic'], 'exon_pos'] = 'single_exon'
merged_df.loc[merged_df['is_first_exon'] & ~merged_df['is_last_exon'] & ~merged_df['is_intronic'], 'exon_pos'] = 'first_exon'
merged_df.loc[~merged_df['is_first_exon'] & merged_df['is_last_exon'] & ~merged_df['is_intronic'], 'exon_pos'] = 'last_exon'

# Assigning expression types
merged_df['exp_type'] = 'underexpressed'
merged_df.loc[merged_df.iloc[:, 2].astype(float) > 1.5, 'exp_type'] = 'overexpressed'

def determine_orientation(row):
    if row['sense'] == row['Strand']:
        return 'same'
    else:
        return 'opposite'

# Compute orientations for each row
merged_df['relative_orientation'] = merged_df.apply(lambda row: determine_orientation(row), axis=1)

# List of columns to keep
columns_to_keep = [
    'Provirus','Chromosome','Start','End','Feature','sense', 'Strand', 'gene_id', 'gene_type', 'gene_name',
    'transcript_id', 'transcript_type', 'transcript_name', 'transcript_support_level',
    'exon_number', 'exon_id', 'is_first_exon', 'is_last_exon', 'exon_pos', 'exp_type', 'relative_orientation'
]

# Select only the columns to keep from merged_df
merged_df = merged_df[columns_to_keep]

# Remove fully duplicate rows
merged_df = merged_df.drop_duplicates()
 
# Count the number of occurrences of each 'Provirus'
provirus_counts = merged_df['Provirus'].value_counts()

# Compute the group counts
group_counts = merged_df.groupby(['Provirus', 'exp_type', 'gene_type', 'exon_pos', 'relative_orientation']).size().reset_index(name='counts')

# Add a new column to group_counts for the total counts of each 'Provirus'
group_counts['total_counts'] = group_counts['Provirus'].map(provirus_counts)

# Calculate the fractional counts
group_counts['fractional_count'] = 1 / group_counts['total_counts']

# Now we can merge 'group_counts' with the original DataFrame
merged_df = pd.merge(merged_df, group_counts, on=['Provirus', 'exp_type', 'gene_type', 'exon_pos', 'relative_orientation'], how='left')

# Check the result
#print(merged_df)

# List of transcript types to keep
#transcript_types_to_keep = ['protein_coding', 'lncRNA', 'retained_intron']
transcript_types_to_keep = ['protein_coding', 'lncRNA']

# Only keep rows where 'transcript_type' is in transcript_types_to_keep
merged_df = merged_df[merged_df['transcript_type'].isin(transcript_types_to_keep)]

# Print to csv
merged_df.to_csv('overlaps_DEprvs_genes.csv', index=False)
