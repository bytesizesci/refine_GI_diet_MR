#!/bin/bash
#SBATCH --job-name=08_plink_clump_GI_step2
#SBATCH --output=config/plink_clump_GI2_%j.out      
#SBATCH --error=config/plink_clump_GI2_%j.err       
#SBATCH --partition=your_partition_name                 # CHANGE: partition name for your system
#SBATCH --nodes=1
#SBATCH --qos=normal                                    # REMOVE or modify if your system doesn't use QoS
#SBATCH --cpus-per-task=1
#SBATCH --mail-type=ALL
#SBATCH --mail-user=your.email@domain.edu               # CHANGE: add your email address
#SBATCH --time=01:00:00                                 # Max run time. Specify the time limit for your allocation
#SBATCH --mem=4G                                        # Memory per job

# Goal: After getting all of the clumped SNPs from the different traits we need
# to run clumping again to make sure they are all independent. However, there are 
# many SNPs that are missing in the reference... we need to try switching the alleles 
# to see if that solves the problem

# ===================== USER CONFIGURATION =====================
# Define this folder (doesn't change)
GI_clumped="GI_clumped"
GI_temp_clumped="GI_temp_clumped"

# Define diet outcome variable
DIET_OUTCOME="alcohol_ALT" 

# Define other paths
DATA_DIR= working_example/data
REF_DIR="${DATA_DIR}/reference_1KGP3_HG19"
OUT_DIR= working_example/results/plink_clumped
GI_TEMP_DIR="${OUT_DIR}/${DIET_OUTCOME}/${GI_temp_clumped}" #where the log files with the missing SNP IDs are
RES_DIR="${OUT_DIR}/${DIET_OUTCOME}/${GI_clumped}" # where the new clumping results will go

# Define input file
file="working_example/results/GI_temp/${DIET_OUTCOME}/all_clumped_snps_${DIET_OUTCOME}_MVMR_GI.txt"

# Load plink
# Customized to Alpine HPC
# Modify for your system
module load plink2/2.00a2.3
# ==============================================================

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Step 1: Create flipped allele .txt and combine with original
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# This code expands to all .log files in the directory
# strips directory and log suffix to get the base name
# writes the extracted SNPs to a missing_snps.txt file in the same directory
for logfile in "$GI_TEMP_DIR"/*.log; do
  base=$(basename "$logfile" .log)
  outname="$GI_TEMP_DIR/${base}.missing_snps.txt"
  grep "Warning:" "$logfile" | awk -F"'" '{print $2}' > "$outname"
done

# Filter the missing SNPs from the original file (that has SNP and P value info)
# First, combine missing SNPs into 1 file
cat "$GI_TEMP_DIR"/*.missing_snps.txt | sort -u > $GI_TEMP_DIR/combined_missing_snps.txt

awk '
  NR==FNR { flip[$1]; next } 
  NR==1 { print; next } 
  $1 in flip {
    split($1, parts, ":")
    flipped_snp = parts[1] ":" parts[2] ":" parts[4] ":" parts[3]
    print flipped_snp, $2
    next
  }
  { print }
' $GI_TEMP_DIR/combined_missing_snps.txt $file > $GI_TEMP_DIR/flipped_summary_file.txt

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Step 2: Rerun plink and write over original plink files
# Should have dramatically fewer warnings about missing SNPs
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Use plink v1.9 for this; just type plink instead of plink2
for chr in {1..22}; do
  plink \
    --bfile "${REF_DIR}/1KGP3_HG19_files_processed_EUR_Jenkai_chr${chr}_chrbpID" \
    --clump "$GI_TEMP_DIR/flipped_summary_file.txt" \
    --clump-snp-field SNP \
    --clump-field P \
    --clump-p1 5e-8 \
    --clump-r2 0.01 \
    --clump-kb 10000 \
    --out "${RES_DIR}/GI_flipped_clumped_results_chr${chr}.tsv"
done



