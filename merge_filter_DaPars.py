#!/usr/bin/env python3

import sys
import pandas as pd
import numpy as np
import os
import glob
from scipy import stats
from statsmodels.stats.multitest import multipletests

def read_key(filename):
    df = pd.read_csv(filename, sep=',', index_col=0, header=None)
    df.index = df.index.str.replace('_trimmed_hisat2_hg38_align-telescope_report.tsv', '', regex=True)
    df[1] = df[1].str.replace('-telescope_report.tsv', '', regex=True)
    return df[1].to_dict()

def process_files(rename_key):
    df_list = []  # This list will store all processed dataframes
    for filename in glob.glob('DaPars2_data.txt_result_temp.*.txt'):
        print(f'Processing file: {filename}')
        df = pd.read_csv(filename, sep='\t')
        columns = list(df.columns)
        print('Renaming columns...')
        for i in range(len(columns)):
            split_column = columns[i].split('/')
            if len(split_column) > 1:
                sample_name = split_column[-1].split('_')[0]
                if sample_name in rename_key:
                    columns[i] = rename_key[sample_name] + "_PDUI"
                else:
                    columns[i] = sample_name + "_PDUI"
        df.columns = columns

        print('Calculating median PDUI values...')
        gtex_columns = [col for col in columns if "GTEX" in col]
        normal_tcga_columns = [col for col in columns if "Normal" in col]
        primary_tcga_columns = [col for col in columns if "Primary" in col]

        control_columns = gtex_columns + normal_tcga_columns
        disease_columns = primary_tcga_columns

        df['Control_median_PDUI'] = df[control_columns].median(axis=1)
        df['Disease_median_PDUI'] = df[disease_columns].median(axis=1)
        df_list.append(df)  # Add the processed dataframe to the list

    # Concatenate all processed dataframes into a single dataframe
    df_all = pd.concat(df_list, ignore_index=True, verify_integrity=True).drop_duplicates()
    
    # Check if there are any duplicate column names
    
    print('Performing statistical tests...')
    ttest_pvalues = []
    kruskal_pvalues = []
    for idx, row in df_all.iterrows():
        control_values = row[control_columns].dropna().values.astype(float)
        disease_values = row[disease_columns].dropna().values.astype(float)
              
        if (len(control_values) > 1 and np.var(control_values) > 0) and \
        (len(disease_values) > 1 and np.var(disease_values) > 0):
            try:
                _, p_ttest = stats.ttest_ind(control_values, disease_values, equal_var = False)
                _, p_kruskal = stats.kruskal(control_values, disease_values)
                kruskal_pvalues.append(p_kruskal)
                ttest_pvalues.append(p_ttest)
            except ValueError:
                kruskal_pvalues.append(np.nan)
                ttest_pvalues.append(np.nan)
        else:
            ttest_pvalues.append(np.nan)
            kruskal_pvalues.append(np.nan)

    # Apply Benjamini-Hochberg procedure to t-test p-values
    mask_invalid_ttest_pvalues = np.isnan(ttest_pvalues)
    valid_ttest_pvalues = [p for p in ttest_pvalues if p == p]  # Excludes NaNs

    _, valid_ttest_pvalues_corrected, _, _ = multipletests(valid_ttest_pvalues, method='fdr_bh')

    # Reinsert corrected p-values back into their original positions
    ttest_pvalues_corrected = np.empty_like(ttest_pvalues)
    ttest_pvalues_corrected[mask_invalid_ttest_pvalues] = np.nan
    ttest_pvalues_corrected[~mask_invalid_ttest_pvalues] = valid_ttest_pvalues_corrected

    df_all['Kruskal-Wallis p_value'] = kruskal_pvalues
    df_all['T-Test Raw p_value'] = ttest_pvalues  # Raw p-values
    df_all['T-Test Adjusted p_value'] = ttest_pvalues_corrected  # Adjusted p-values
        
    # Reorder the columns of df_all
    selected_columns = ['Gene', 'fit_value', 'Predicted_Proximal_APA', 'Loci', 'Control_median_PDUI', 'Disease_median_PDUI', 'T-Test Raw p_value', 'T-Test Adjusted p_value', 'Kruskal-Wallis p_value']
    all_columns = df_all.columns.tolist()
    extra_columns = [col for col in all_columns if col not in selected_columns]
    df_all = df_all[selected_columns + extra_columns]
    
        
    # Write the resulting dataframe to the output
    outputname = 'DaPars2_output_merged.tsv'
    print(f'Writing to output file: {outputname}\n')
    df_all.to_csv(outputname, sep='\t', index=False)


#### FILTER STEP ####
    # Filter rows based on adjusted p-value
    #filtered_df = df_all[(df_all['T-Test Adjusted p_value'] <= 0.05) & (~df_all['T-Test Adjusted p_value'].isnull())]
###########
#### UNFILTERED ####
    # Filter rows based on adjusted p-value
    filtered_df = df_all

##########
    # Write the filtered dataframe to the output file
    filtered_outputname = 'DaPars2_output_merged_filtered.tsv'
    print(f'Writing filtered output file: {filtered_outputname}\n')
    filtered_df.to_csv(filtered_outputname, sep='\t', index=False)

    
        
if __name__ == "__main__":
    rename_key = read_key('/Users/Daniel/Documents/Med_School/G1-G4/Coffin_Lab/Thesis/lncRNA_Project/Scripts/report_renaming_key.csv')
    process_files(rename_key)
    
