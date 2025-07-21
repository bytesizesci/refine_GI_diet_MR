# Author: KJS
# Date: June 23, 2025
# Project: Refining diet MR

# Goal: After getting all of the clumped SNPs from the different traits we need
# to run clumping again to make sure they are all independent.

# To run: Be in /pl/active/colelab/users/kjames/refinedMR/scripts/MVMR
# Make sure to be on a worker node - acompile (4 hrs)
# acompile --time=04:00:00
# Might have to change permissions the first time - chmod +x <scriptname>
# ./MVMR_plink_clump_alcohol_cirrhosis_clumpedGI_250623.sh

# Get plink
module load plink2/2.00a2.3

# MANUALLY NAME THIS
diet_outcome="alcohol_cirrhosis"
file="/pl/active/colelab/users/kjames/refinedMR/interim_data/MVMR/${diet_outcome}/all_clumped_snps_${diet_outcome}_MVMR_GIclump.txt"

# Define this folder (doesn't change)
GI_clumped="GI_clumped"

# Make a folder for the clumped results if it doesn't exist
if [ ! -d "/pl/active/colelab/users/kjames/refinedMR/interim_data/MVMR/${diet_outcome}/${GI_clumped}" ]; then
  mkdir /pl/active/colelab/users/kjames/refinedMR/interim_data/MVMR/${diet_outcome}/${GI_clumped}
fi

# Use plink v1.9 for this; just type plink instead of plink2
for chr in {1..22}; do
  plink \
    --bfile /pl/active/colelab/users/kjames/refinedMR/interim_data/1KGP3_HG19/1KGP3_HG19_files_processed_EUR_Jenkai_chr${chr}_chrbpID \
    --clump ${file} \
    --clump-snp-field SNP \
    --clump-field P \
    --clump-p1 5e-8 \
    --clump-r2 0.01 \
    --clump-kb 10000 \
    --out /pl/active/colelab/users/kjames/refinedMR/interim_data/MVMR/${diet_outcome}/${GI_clumped}/GI_clumped_results_${diet_outcome}_chr${chr}.tsv
done







