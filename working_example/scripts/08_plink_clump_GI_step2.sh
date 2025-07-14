#!/bin/bash
#SBATCH --job-name=07_plink_clump_GI
#SBATCH --output=/your/path/to/logs/ldsc_rg_%j.out      # CHANGE: set to your own logs directory
#SBATCH --error=/your/path/to/logs/ldsc_rg_%j.err       # CHANGE: set to your own logs directory
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
diet_outcome="alcohol_ALT"

file="~/working_example/data/GI_temp/${diet_outcome}/all_clumped_snps_${diet_outcome}_MVMR_GI.txt"

# Define this folder (doesn't change)
GI_clumped="GI_clumped"

DATA_DIR= "~/working_example/data" 
REF_DIR="${DATA_DIR}/reference_1KGP3_HG19"
OUTDIR="${DATA_DIR}/plink_clumped"
DIET_OUTCOME="alcohol_ALT" #diet_outcome
trait_outdir="${DATA_DIR}/${DIET_OUTCOME}/${GI_clumped}"

# Define where the log files with the missing SNP IDs are
logdir="${DATA_DIR}/${DIET_OUTCOME}/${GI_clumped}"

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
for logfile in "$logdir"/*.log; do
  base=$(basename "$logfile" .log)
  outname="$logdir/${base}.missing_snps.txt"
  grep "Warning:" "$logfile" | awk -F"'" '{print $2}' > "$outname"
done

# Filter the missing SNPs from the original file (that has SNP and P value info)
# First, combine missing SNPs into 1 file
cat "$logdir"/*.missing_snps.txt | sort -u > $out_path/combined_missing_snps.txt

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
' $out_path/combined_missing_snps.txt $file > $out_path/flipped_summary_file.txt

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Step 2: Rerun plink and write over original plink files
# Should have dramatically fewer warnings about missing SNPs
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Use plink v1.9 for this; just type plink instead of plink2
for chr in {1..22}; do
  plink \
    --bfile "${REF_DIR}/1KGP3_HG19_files_processed_EUR_Jenkai_chr${chr}_chrbpID" \
    --clump "$out_path/flipped_summary_file.txt" \
    --clump-snp-field SNP \
    --clump-field P \
    --clump-p1 5e-8 \
    --clump-r2 0.01 \
    --clump-kb 10000 \
    --out "${trait_outdir}/GI_clumped_results_chr${chr}.tsv"
done



