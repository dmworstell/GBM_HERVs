#!/usr/bin/bash
#SBATCH --job-name=gene_prv_OLs		# Job name
#SBATCH --mail-type=FAIL		# Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=Daniel.Murimi_Worstell@tufts.edu
#SBATCH --partition=batch,largemem,gpu
#SBATCH --nodes=1			# Run a single task
#SBATCH --cpus-per-task=12		# Number of CPU cores per task
#SBATCH --mem-per-cpu=4gb
#SBATCH --time=4-00:00:00
#SBATCH --output=%j.out
#SBATCH --error=%j.err

echo LOADING ANACONDA
module load anaconda/2024.06-py312

echo ACTIVATING ENVIRONMENT
source activate /cluster/tufts/coffinlab/dworst01/condaenv/newest_python

echo CHANGING DIRECTORY
cd ..

echo RUNNING SCRIPT
#python Provirus_Gene_OLs_hg38_Bootstraps_mt.py
python Provirus_Gene_OLs_hg38_Bootstraps_mt_using_dict.py
echo DONE