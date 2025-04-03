#!/usr/bin/env python3
import pandas as pd
import pickle 

#Reading in overlaps
#ALL HERVS
#merged_df = pd.read_csv('/cluster/tufts/coffinlab/dworst01/lncrna_project/rnaseq/ref/merged_df_overlaps_ALLprvs.csv',sep=',', header=0)
#DE HERVS
merged_df = pd.read_csv('/cluster/tufts/coffinlab/dworst01/lncrna_project/rnaseq/ref/overlaps_DEprvs_genes.csv',sep=',', header=0)

# Count the number of occurrences of each 'Provirus'
provirus_counts = merged_df['Provirus'].value_counts()
provirus_list = provirus_counts.index.values.astype(str)

#Create dictionary for each provirus
provirus_dict = {provirus: merged_df[merged_df['Provirus'] == provirus] for provirus in provirus_list}

#Save the dictionary
with open('PRV_Gene_OL_dict_DEprvs.pkl', 'wb') as f:
    pickle.dump(provirus_dict, f)