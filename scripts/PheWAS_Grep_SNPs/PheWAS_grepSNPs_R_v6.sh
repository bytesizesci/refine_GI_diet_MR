#!/bin/bash
#SBATCH --job-name=getSNPs_BothSex_v6
#SBATCH --output=/pl/active/colelab/users/kjames/refinedMR/config/job_%j.out
#SBATCH --error=/pl/active/colelab/users/kjames/refinedMR/config/job_%j.err
#SBATCH --partition=amilan
#SBATCH --array=3995-4586 #4586 files total
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=7G        # Adjust the memory limit as per your needs
#SBATCH --time=00:30:00  # Specify the time limit for your allocation

#https://researchcomputing.princeton.edu/support/knowledge-base/slurm#terminology 
#helps explain resource terminology
#this is a serial job


# Be in "/pl/active/colelab/users/kjames/refinedMR" when you run this sbatch
# To run, 
# sbatch scripts/PheWAS_Grep_SNPs/PheWAS_grepSNPs_R_v6.sh

# Just in case, set wd using full path
cd /pl/active/colelab/users/kjames/refinedMR

# Load R
module load R/4.2.2

# Set variables 

# run Rscript
Rscript /pl/active/colelab/users/kjames/refinedMR/scripts/PheWAS_Grep_SNPs/PheWAS_grepSNPs_R_v6.R $SLURM_ARRAY_TASK_ID 