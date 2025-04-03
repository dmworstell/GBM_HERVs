#!/usr/bin/env python3
import os
import numpy as np
import pandas as pd
from pybedtools import BedTool

os.chdir('/Users/Daniel/Documents/Med_School/G1-G4/Coffin_Lab/Thesis/lncRNA_Project/RNA_Seq_Processing/TCGA_Analysis/alternative_polyA/DaPars2_results')

# Load the ALL HERVs - genes overlap file
annotation_path = "/Users/Daniel/Documents/Med_School/G1-G4/Coffin_Lab/Thesis/lncRNA_Project/RNA_Seq_Processing/TCGA_Analysis/Provirus_Gene_OLs/overlaps_ALLprvs_exons.csv"


annotation_df = pd.read_csv(annotation_path, sep=",", header=0)


# Read the DaPars2 output
dapars_df = pd.read_csv('DaPars2_output_merged_filtered.tsv', sep='\t')


dapars_df[['chr', 'coords']] = dapars_df['Loci'].str.split(':', 1, expand=True)
dapars_df[['start', 'end']] = dapars_df['coords'].str.split('-', 1, expand=True)

# Split the 'Gene' column into separate fields
dapars_df[['Transcript_ID', 'Gene_ID', 'Chr', 'Strand']] = dapars_df['Gene'].str.split('|', 3, expand=True)
dapars_df.rename(columns={"Strand": "Gene_Strand"}, inplace=True)

# Convert 'start' and 'end' to integer
dapars_df['start'] = dapars_df['start'].astype(int)
dapars_df['end'] = dapars_df['end'].astype(int)
dapars_df.rename(columns={"start": "3-prime UTR start", "end": "3-prime UTR end"}, inplace=True)


# Create a BED formatted DataFrame
dapars_df_bed = dapars_df[['chr', '3-prime UTR start', '3-prime UTR end']].copy()


# Write the BED file
dapars_df_bed.to_csv('DaPars2_output.bed', sep='\t', header=False, index=False)

#Isolate dapars APA UTRs that overlap proviral fragments
annotation_df_bed = annotation_df[['Chromosome','Start','End']].copy()
annotation_df_bed.to_csv('prv_fragments.bed', sep='\t', header=False, index=False)

prv_fragments_bed = BedTool('prv_fragments.bed')
dapars_bed = BedTool('DaPars2_output.bed')
intersected_bed = dapars_bed.intersect(prv_fragments_bed, wa=True, wb=True, v=True)


# Load intersected into a DataFrame
intersected_df = pd.read_table(intersected_bed.fn, header=None)

# Convert 'start' and 'end' to integer in intersected_df
intersected_df[1] = intersected_df[1].astype(int)
intersected_df[2] = intersected_df[2].astype(int)

# Specify the columns you want to keep from each dataframe
dapars_columns_to_keep = ['Gene', 'fit_value', 'Predicted_Proximal_APA', 'Loci', 'Control_median_PDUI','Disease_median_PDUI', 'T-Test Raw p_value', 'T-Test Adjusted p_value', 'Kruskal-Wallis p_value', 'chr', 'coords', '3-prime UTR start', '3-prime UTR end', 'Transcript_ID', 'Gene_ID', 'Chr', 'Gene_Strand']

# Drop unwanted columns
dapars_df = dapars_df[dapars_columns_to_keep]
                     
# Assign column names to the dataframes

#annotation_df.columns = ['Chromosome','Frag_Start','Frag_End','Provirus_ID','Provirus_Sense','Gene_Strand','gene_id','gene_type','gene_name','transcript_id','transcript_type','transcript_name','transcript_support','exon_number','exon_id','is_first_exon','is_last_exon','Terminal_exon_overlap','relative_orientation','counts','total_counts','fractional_counts']
#annotation_columns_to_keep = ['Chromosome', 'Frag_Start', 'Frag_End', 'Provirus_ID', 'Provirus_Sense', 'Gene_Strand','transcript_type','gene_name','exon_number','Terminal_exon_overlap']


# Perform the merge, but only keep specified columns from each dataframe
merged_all = intersected_df.merge(dapars_df[dapars_columns_to_keep], left_on=[0, 1, 2], right_on=['chr', '3-prime UTR start', '3-prime UTR end'], how='left')

# Cleaning up the column names

merged_all.columns = ['3-prime UTR chr', '3-prime UTR start', '3-prime UTR end', 'Gene', 'fit_value', 'Predicted_Proximal_APA', 'Loci', 'Control_median_PDUI', 'Disease_median_PDUI', 'T-Test Raw p_value', 'T-Test Adjusted p_value', 'Kruskal-Wallis p_value', 'DaPars_chr', 'DaPars_coords', 'DaPars_3-prime UTR start', 'DaPars_3-prime UTR end', 'Transcript_ID', 'Gene_name', 'DaPars_Chr', 'Gene_Strand']


# Drop the now unnecessary info columns
merged_all.drop(['DaPars_Chr','DaPars_coords','Loci','3-prime UTR chr','3-prime UTR start','3-prime UTR end','Gene'], axis=1, inplace=True)

columns_order = ['Gene_name', 'DaPars_chr', 'Gene_Strand','DaPars_3-prime UTR start', 'DaPars_3-prime UTR end', 'fit_value', 'Predicted_Proximal_APA', 'Control_median_PDUI', 'Disease_median_PDUI', 'T-Test Raw p_value', 'T-Test Adjusted p_value', 'Kruskal-Wallis p_value', 'Transcript_ID']

merged_all = merged_all.reindex(columns=columns_order)

# Save the result
merged_all.to_csv('Intersect_NOprvs_Dapars', sep='\t', index=False)


# Create a column that represents the difference between Control_median_PDUI and Disease_median_PDUI
merged_all['PDUI_difference'] = (merged_all['Control_median_PDUI'] - merged_all['Disease_median_PDUI'])

# Explicitly create a new DataFrame instead of a view
simplified_df = merged_all.sort_values(['PDUI_difference', 'T-Test Raw p_value'], ascending=[False, True]).drop_duplicates(['Gene_name'], keep='first').copy()

# Save to a TSV file
simplified_df.to_csv("simplified_sorted_NOprvs.tsv", sep="\t", index=False)

