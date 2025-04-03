#!/usr/bin/bash
#SBATCH --job-name=gen_dict	# Job name
#SBATCH --mail-type=ALL		# Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=Daniel.Murimi_Worstell@tufts.edu
#SBATCH --partition=batch,largemem,gpu
#SBATCH --nodes=1			# Run a single task
#SBATCH --cpus-per-task=1		# Number of CPU cores per task
#SBATCH --mem-per-cpu=24gb
#SBATCH --time=3-00:00:00
#SBATCH --output=%j.out
#SBATCH --error=%j.err

module load anaconda/2023.07
source activate newest_python

cd ..
echo RUNNING SCRIPT
python Provirus_Gene_Overlaps_Provirus_Dict_Generator.sh
echo DONE