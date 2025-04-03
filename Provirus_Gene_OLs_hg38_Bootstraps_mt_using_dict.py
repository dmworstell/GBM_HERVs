#!/usr/bin/env python3
import pandas as pd
import numpy as np
from tqdm import tqdm
from multiprocessing import Pool
from bootstrap_fun_using_dict_GenePrv import bootstrap
from weighted_random_choice import weighted_random_choice
import pickle

#Reading in overlaps
#merged_df = pd.read_csv('/cluster/tufts/coffinlab/dworst01/lncrna_project/rnaseq/ref/merged_df_overlaps_ALLprvs.csv',sep=',', header=0)
merged_df = pd.read_csv('/cluster/tufts/coffinlab/dworst01/lncrna_project/rnaseq/ref/fixed_overlaps_DEprvs_genes.csv',sep=',', header=0)

# Count the number of occurrences of each 'Provirus'
provirus_counts = merged_df['Provirus'].value_counts()
provirus_list = provirus_counts.index.values.astype(str)

##Begin monte carlo simulation + bootstrapping
###Randomly select one provirus###

print("BEGIN MONTE CARLO SIMULATION")
# Number of bootstrap iterations
n_iterations = 1000

# Initialize lists to store the results
mRNA_results = []
lncRNA_results = []
rng = np.random.default_rng()

index = ['same', 'opposite']
columns = ['first_exon','middle_exon','last_exon','single_exon','intron']
data = {col: [0, 0] for col in columns}

#Read in previously created dictionary for each provirus
#with open('PRV_Gene_OL_dict.pkl', 'rb') as f:
with open('PRV_Gene_OL_dict_DEprvs.pkl', 'rb') as f:

    provirus_dict = pickle.load(f)

if __name__ == "__main__":
    # Number of processes to use
    n_processes = 12

    # Create a pool of worker processes
    with Pool(n_processes) as p:
        # Run the function for each iteration
        results = p.starmap(bootstrap, [(provirus_dict, provirus_list, rng, data, columns, index, weighted_random_choice) for _ in range(n_iterations)])
    
    # The results variable now contains the results from all processes
    mRNA_results, lncRNA_results = zip(*results)

# Convert the lists to numpy arrays for easier manipulation
mRNA_results = np.array(mRNA_results)
lncRNA_results = np.array(lncRNA_results)


mRNA_mc_counts = pd.DataFrame(data, columns=columns, index=index)
lncRNA_mc_counts = pd.DataFrame(data, columns=columns, index=index)

x_lncRNA = [f'{index} | {column}'.replace('first', 'First').replace('middle', 'Middle').replace('last', 'Last').replace('single', 'Single').replace('intron', 'Intron').replace('opposite', '-').replace('same', '+').replace('_',' ').replace('exon','') for index in lncRNA_mc_counts.index for column in lncRNA_mc_counts.columns]
x_mRNA = [f'{index} | {column}'.replace('first', 'First').replace('middle', 'Middle').replace('last', 'Last').replace('single', 'Single').replace('intron', 'Intron').replace('opposite', '-').replace('same', '+').replace('_',' ').replace('exon','') for index in mRNA_mc_counts.index for column in mRNA_mc_counts.columns]

##Write out the files for use later
mRNA_results_with_headers = pd.DataFrame(mRNA_results)
mRNA_results_with_headers.columns = x_mRNA

lncRNA_results_with_headers = pd.DataFrame(lncRNA_results)
lncRNA_results_with_headers.columns = x_lncRNA

mRNA_results_with_headers.to_csv("/cluster/tufts/coffinlab/dworst01/lncrna_project/rnaseq/results/bootstrapping_work/mRNA_proviruses_1000_bootstraps_7146_MT_100percentExon_DEHervs.csv", header=True)
lncRNA_results_with_headers.to_csv("/cluster/tufts/coffinlab/dworst01/lncrna_project/rnaseq/results/bootstrapping_work/lncRNA_proviruses_1000_bootstraps_7146_MT_100percentExon_DEHervs.csv", header=True)
