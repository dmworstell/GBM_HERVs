#!/usr/bin/env python3
import os
import numpy as np
import pandas as pd
from pybedtools import BedTool

os.chdir('/Users/Daniel/Documents/Med_School/G1-G4/Coffin_Lab/Thesis/lncRNA_Project/RNA_Seq_Processing/TCGA_Analysis/alternative_polyA/DaPars2_results')

# Load the ALL HERVs - genes overlap file
annotation_path = "/Users/Daniel/Documents/Med_School/G1-G4/Coffin_Lab/Thesis/lncRNA_Project/RNA_Seq_Processing/TCGA_Analysis/Provirus_Gene_OLs/overlaps_ALLprvs_genes.csv"

annotation_df = pd.read_csv(annotation_path, sep=",", header=0)
del(annotation_path)

# Read in a GTF to get full gene coordinates
Gencode_path = '/Users/Daniel/Documents/Med_School/G1-G4/Coffin_Lab/Thesis/lncRNA_Project/RNA_Seq_Processing/TCGA_Analysis/Annotations/gencode.v44.chr_patch_hapl_scaff.annotation.gtf'
Gencode_df = pd.read_table(Gencode_path, sep="\t", header=None, comment='#')
del(Gencode_path)
Gencode_df['Transcript_ID'] = Gencode_df[8].str.extract('transcript_id "([^"]*)"')

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


dapars_df = dapars_df.merge(Gencode_df[[3, 4, 'Transcript_ID']], how='left', on='Transcript_ID')
dapars_df.rename(columns={3: 'start_coordinate', 4: 'end_coordinate'}, inplace=True)
del(Gencode_df)

# Remove rows where 'start_coordinate' or 'end_coordinate' is NaN
dapars_df.dropna(subset=['start_coordinate', 'end_coordinate'], inplace=True)

# Create a BED formatted DataFrame
dapars_df_bed = dapars_df[['chr', 'start_coordinate', 'end_coordinate']].copy()
# Convert 'start' and 'end' to integer

dapars_df_bed['start_coordinate'] = dapars_df_bed['start_coordinate'].astype(int)
dapars_df_bed['end_coordinate'] = dapars_df_bed['end_coordinate'].astype(int)

# Write the BED file
dapars_df_bed.to_csv('DaPars2_output.bed', sep='\t', header=False, index=False)
del(dapars_df_bed)

#Isolate dapars genes that overlap proviral fragments
annotation_df_bed = annotation_df[['Chromosome','Start','End']].copy()
annotation_df_bed.to_csv('prv_fragments.bed', sep='\t', header=False, index=False)

del(annotation_df_bed)

prv_fragments_bed = BedTool('prv_fragments.bed')


dapars_bed = BedTool('DaPars2_output.bed')
intersected_bed = dapars_bed.intersect(prv_fragments_bed, wa=True, wb=True)

del(dapars_bed)
del(prv_fragments_bed)

# Load intersected into a DataFrame
intersected_df = pd.read_table(intersected_bed.fn, header=None)
del(intersected_bed)

# Convert numerical columns to integers in intersected_df
intersected_df[1] = intersected_df[1].astype(int)
intersected_df[2] = intersected_df[2].astype(int)
intersected_df[4] = intersected_df[4].astype(int)
intersected_df[5] = intersected_df[5].astype(int)


# Specify the columns you want to keep from each dataframe
dapars_columns_to_keep = ['Gene', 'fit_value', 'Predicted_Proximal_APA', 'Loci', 'Control_median_PDUI','Disease_median_PDUI', 'T-Test Raw p_value', 'T-Test Adjusted p_value', 'Kruskal-Wallis p_value', 'chr', 'coords', '3-prime UTR start', '3-prime UTR end', 'Transcript_ID', 'Gene_ID', 'Chr', 'Gene_Strand','start_coordinate','end_coordinate']

# Drop unwanted columns
dapars_df = dapars_df[dapars_columns_to_keep]
                     
# Assign column names to the dataframes

annotation_df.columns = ['Chromosome','Frag_Start','Frag_End','Provirus_ID','Provirus_Sense','Gene_Strand','gene_id','gene_type','gene_name','transcript_id','transcript_type','transcript_name','transcript_support','exon_number','exon_id','is_first_exon','is_last_exon','Terminal_exon_overlap','relative_orientation','counts','total_counts','fractional_counts']
annotation_columns_to_keep = ['Chromosome', 'Frag_Start', 'Frag_End', 'Provirus_ID', 'Provirus_Sense', 'Gene_Strand','transcript_type','gene_name','exon_number','Terminal_exon_overlap','relative_orientation']

## TEST MERGE ##

#intersected_df_test = intersected_df[1:1000]
#dapars_df_test = dapars_df[1:100]
#annotation_df_test = annotation_df[1:100]

#merged_test1 = intersected_df_test.merge(dapars_df_test[dapars_columns_to_keep], left_on=[0, 1, 2], right_on=['chr', 'start_coordinate', 'end_coordinate'], how='left')
#merged_test2 = merged_test1.merge(annotation_df_test[annotation_columns_to_keep], left_on=[3, 4, 5], right_on=['Chromosome', 'Frag_Start', 'Frag_End'], how='left')


# Perform the merge, but only keep specified columns from each dataframe
merged_all = intersected_df.merge(dapars_df[dapars_columns_to_keep], left_on=[0, 1, 2], right_on=['chr', 'start_coordinate', 'end_coordinate'], how='left')
merged_all = merged_all.merge(annotation_df[annotation_columns_to_keep], left_on=[3, 4, 5], right_on=['Chromosome', 'Frag_Start', 'Frag_End'], how='left')

del(intersected_df,dapars_df,annotation_df,dapars_columns_to_keep,annotation_columns_to_keep)

# Cleaning up the column names

merged_all.columns = ['3-prime UTR chr', '3-prime UTR start', '3-prime UTR end', 'Provirus_chr_intersect', 'Provirus_start_intersect', 'Provirus_end_intersect', 'Gene', 'fit_value', 'Predicted_Proximal_APA', 'Loci', 'Control_median_PDUI', 'Disease_median_PDUI', 'T-Test Raw p_value', 'T-Test Adjusted p_value', 'Kruskal-Wallis p_value', 'DaPars_chr', 'DaPars_coords', 'DaPars_3-prime UTR start', 'DaPars_3-prime UTR end', 'Transcript_ID', 'Gene_name', 'DaPars_Chr', 'Gene_Strand', 'Gene_Start','Gene_End', 'Provirus_chr', 'Provirus_Frag_Start', 'Provirus_Frag_End', 'Provirus_ID', 'Provirus_Sense','Gene_Strand_dup','transcript_type','gene_name','Overlapping_exon_number','Terminal_exon_overlap','Sense_Strand_Match']


# Drop the now unnecessary info columns
merged_all.drop(['DaPars_chr', 'DaPars_coords','Loci','3-prime UTR chr','3-prime UTR start','3-prime UTR end','Gene','Provirus_chr_intersect','Provirus_start_intersect','Provirus_end_intersect','Gene_Strand_dup','Gene','fit_value','Loci','DaPars_3-prime UTR start', 'DaPars_3-prime UTR end', ], axis=1, inplace=True)

#merged_all['Sense_Strand_Match'] = merged_all['Provirus_Sense'] == merged_all['Gene_Strand']

columns_order = ['Gene_name', 'Gene_Start', 'Gene_End', 'Provirus_chr', 'Provirus_Frag_Start', 'Provirus_Frag_End', 'Provirus_Sense', 'Gene_Strand', 'Provirus_ID', 'Sense_Strand_Match', 'DaPars_3-prime UTR start', 'DaPars_3-prime UTR end', 'fit_value', 'Predicted_Proximal_APA', 'Control_median_PDUI', 'Disease_median_PDUI', 'T-Test Raw p_value', 'T-Test Adjusted p_value', 'Kruskal-Wallis p_value', 'Transcript_ID', 'Overlapping_exon_number','Terminal_exon_overlap']

merged_all = merged_all.reindex(columns=columns_order)

# Save the result
merged_all.to_csv('Intersect_ALLprvs_Dapars_raw', sep='\t', index=False)


# Fill NA/NaN values with 0
merged_all['Overlapping_exon_number'] = merged_all['Overlapping_exon_number'].fillna(0)

# Convert 'Overlapping_exon_number' to int for sorting
merged_all['Overlapping_exon_number'] = merged_all['Overlapping_exon_number'].astype(int)


# Create a column that represents the difference between Control_median_PDUI and Disease_median_PDUI
merged_all['PDUI_difference'] = (merged_all['Control_median_PDUI'] - merged_all['Disease_median_PDUI'])

# Explicitly create a new DataFrame instead of a view
simplified_df = merged_all.sort_values(['Terminal_exon_overlap', 'Overlapping_exon_number', 'PDUI_difference', 'T-Test Raw p_value'], ascending=[False, False, False, True]).drop_duplicates(['Provirus_ID', 'Gene_name'], keep='first').copy()

# Replace both integer 0s and string '0's with NaN
simplified_df['Overlapping_exon_number'] = simplified_df['Overlapping_exon_number'].replace([0, '0'], np.nan)



# Save to a TSV file
simplified_df.to_csv("simplified_sorted_ALLprvs.tsv", sep="\t", index=False)

