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
provirus_dict = prv_key.groupby(8)[[0, 3, 4, 6]].apply(lambda x: list(map(tuple, x.values))).to_dict()

# Reading GENCODE GTF
print("READING GENCODE GTF")
gencode_file = '/Users/Daniel/Documents/Med_School/G1-G4/Coffin_Lab/Thesis/Thesis_Project/RNA_Seq_Processing/TCGA_Analysis/Annotations/gencode.v44.annotation.gtf'
gencode_ranges = pr.read_gtf(gencode_file)

#Get rid of nonsense-mediated-decay products
gencode_ranges = gencode_ranges[gencode_ranges.transcript_type != 'nonsense_mediated_decay']

# Extract exon intervals from the GENCODE data
exons_ranges = gencode_ranges[gencode_ranges.Feature == 'exon']

# Maintain a dictionary that maps each transcript ID to its highest exon number
transcript_max_exon_num = defaultdict(int)

print("GENERATING MAX EXONS")
transcript_max_exon_num = exons_ranges.df['exon_number'].astype(int).groupby(exons_ranges.df['transcript_id']).max()

print("GENERATING DataFrame TO HOLD ALL PROVIRUSES")
# First, create a DataFrame with all proviruses
provirus_df = pd.concat((pd.DataFrame(fragments, columns=['Chromosome', 'Start', 'End', 'sense']).assign(Provirus=provirus_id) for provirus_id, fragments in provirus_dict.items()), ignore_index=True)

print("Generating PyRanges object")
provirus_range = pr.PyRanges(provirus_df)

# Then intersect/join
#print("RUNNING INTERSECTION")
#overlaps = provirus_range.join(exons_ranges)
#print("DONE")

# Then intersect/join
print("RUNNING INTERSECTION")
overlaps = provirus_range.join(gencode_ranges)
print("DONE")

# Convert PyRanges to DataFrame
df_overlaps = overlaps.df

# Save DataFrame to CSV
#df_overlaps.to_csv('overlaps_AllPrvs.csv', index=False)

gene_ids = set()
gene_orientation = defaultdict(list)
orientation_count = defaultdict(lambda: defaultdict(int))  # For keeping track of the different combinations

# Specify the categories of transcripts to include
included_genetypes = {'lncRNA', 'protein_coding'}


# Create a dictionary to store the orientations of the proviruses overlapping each exon
exon_orientations = defaultdict(list)

# Preprocessing steps
df_overlaps['gene_id'] = df_overlaps['gene_id'].str.split('.').str[0]  # Remove .# suffix from column 17
gene_orientation_dict = dict(zip(df_overlaps['gene_id'], df_overlaps['Strand']))

df_overlaps['exon_number'] = df_overlaps['exon_number'].astype('Int64')
df_overlaps['transcript_id'] = df_overlaps['transcript_id'].astype(str)


# Creating provirus_info DataFrame
provirus_info_df = df_overlaps.iloc[:, 0:6]
provirus_info_dict = provirus_info_df.groupby('Provirus').apply(lambda x: x.values.tolist()).to_dict()
Proviruses = set(df_overlaps['Provirus'])

# Determining exon positions
df_overlaps['is_intronic'] = df_overlaps['exon_number'].isna()
df_overlaps['exon_number'] = df_overlaps['exon_number'].fillna(False)
df_overlaps['is_first_exon'] = df_overlaps['exon_number'] == 1
df_overlaps['is_last_exon'] = df_overlaps['exon_number'] == df_overlaps['transcript_id'].map(transcript_max_exon_num)
df_overlaps['exon_pos'] = 'middle_exon'

df_overlaps.loc[df_overlaps['is_intronic'], 'exon_pos'] = 'intron'
df_overlaps.loc[df_overlaps['is_first_exon'] & df_overlaps['is_last_exon'] & ~df_overlaps['is_intronic'], 'exon_pos'] = 'single_exon'
df_overlaps.loc[df_overlaps['is_first_exon'] & ~df_overlaps['is_last_exon'] & ~df_overlaps['is_intronic'], 'exon_pos'] = 'first_exon'
df_overlaps.loc[~df_overlaps['is_first_exon'] & df_overlaps['is_last_exon'] & ~df_overlaps['is_intronic'], 'exon_pos'] = 'last_exon'


def determine_orientation(row):
    if row['sense'] == row['Strand']:
        return 'same'
    else:
        return 'opposite'

# Compute orientations for each row
df_overlaps['relative_orientation'] = df_overlaps.apply(lambda row: determine_orientation(row), axis=1)

# List of columns to keep
columns_to_keep = [
    'Chromosome','Start','End', 'Provirus', 'sense', 'Strand', 'gene_id', 'gene_type', 'gene_name',
    'transcript_id', 'transcript_type', 'transcript_name', 'transcript_support_level',
    'exon_number', 'exon_id', 'is_first_exon', 'is_last_exon', 'exon_pos', 'relative_orientation'
]

# Select only the columns to keep from merged_df
merged_df = df_overlaps[columns_to_keep]

# Remove fully duplicate rows
merged_df = merged_df.drop_duplicates()

# Count the number of occurrences of each 'Provirus'
provirus_counts = merged_df['Provirus'].value_counts()

# Compute the group counts
group_counts = merged_df.groupby(['Provirus', 'gene_type', 'exon_pos', 'relative_orientation']).size().reset_index(name='counts')

# Add a new column to group_counts for the total counts of each 'Provirus'
group_counts['total_counts'] = group_counts['Provirus'].map(provirus_counts)

# Calculate the fractional counts
group_counts['fractional_count'] = 1 / group_counts['total_counts']

# Now we can merge 'group_counts' with the original DataFrame
merged_df = pd.merge(merged_df, group_counts, on=['Provirus', 'gene_type', 'exon_pos', 'relative_orientation'], how='left')

# Check the result
#print(merged_df)

# List of transcript types to keep
#transcript_types_to_keep = ['protein_coding', 'lncRNA', 'retained_intron']
transcript_types_to_keep = ['protein_coding', 'lncRNA']

# Only keep rows where 'transcript_type' is in transcript_types_to_keep
merged_df = merged_df[merged_df['transcript_type'].isin(transcript_types_to_keep)]

# Print to csv
merged_df.to_csv('merged_df_overlaps_ALLprvs.csv', index=False)

result = merged_df[merged_df['gene_type'].isin(included_genetypes)].groupby(['gene_type', 'exon_pos', 'relative_orientation'])['fractional_count'].sum()

result_df = result.reset_index(name='count')
   
print(result_df)

#df_plot = result_df.set_index(['exp_type']).sort_index()
df_plot = result_df
    
# Pivot the data: index='relative_orientation', columns='gene_type', values='count'
pivot_data = df_plot.pivot_table(index=['relative_orientation', 'exon_pos'], columns='gene_type', values='count', aggfunc='sum', fill_value=0)

# Plot the pivoted data as a stacked bar plot
ax = pivot_data.plot(kind='bar', stacked=True, title="All HERV-Gene Overlaps")
ax.set_ylabel('Count')
ax.set_xlabel('Exon Overlap Category')
title = ax.set_title('All Exonized HERVs')
title.set_size(24)

# Set the ylabel
ylabel = ax.set_ylabel('Count')
ylabel.set_size(20)

# Set the xlabel
xlabel = ax.set_xlabel('Exon Overlap Category')
xlabel.set_size(20)

# Increase the size of the axis labels and legend labels
ax.tick_params(axis='both', which='major', labelsize=20)
ax.legend(prop={'size': 20})

# Replace the word "exon" in each x-axis category label
new_labels = [label.get_text().replace('exon', '') for label in ax.get_xticklabels()]
ax.set_xticklabels(new_labels)

# Change the location of the legend to the upper left
legend = ax.legend(prop={'size': 20}, loc='upper left')

# Change the label for "protein_coding" in the legend to "mRNA"
for text in legend.get_texts():
    if text.get_text() == 'protein_coding':
        text.set_text('mRNA')
        
# Replace any underscores with spaces
new_labels = [label.replace('_', ' ') for label in new_labels]
new_labels = [label.replace('(', '') for label in new_labels]
new_labels = [label.replace(')', '') for label in new_labels]
new_labels = [label.replace('opposite', '-') for label in new_labels]
new_labels = [label.replace('same', '+') for label in new_labels]
new_labels = [label.replace(',', '|') for label in new_labels]
new_labels = [label.replace('first', 'First') for label in new_labels]
new_labels = [label.replace('middle', 'Middle') for label in new_labels]
new_labels = [label.replace('last', 'Last') for label in new_labels]
new_labels = [label.replace('single', 'Single') for label in new_labels]
ax.set_xticklabels(new_labels)

plt.show()
