# Author: KJS
# Date: June 23, 2025
# Project: Refining diet MR

# Goal: After getting all of the clumped SNPs from the different traits we need
# to run clumping again to make sure they are all independent. However, there are 
# many SNPs that are missing in the reference... we need to try switching the alleles 
# to see if that solves the problem

# To run: Be in /pl/active/colelab/users/kjames/refinedMR/scripts/MVMR
# Make sure to be on a worker node - acompile (4 hrs)
# acompile --time=04:00:00
# Might have to change permissions the first time - chmod +x <scriptname>
# ./MVMR_plink_clump_<name>.sh

# Get plink
module load plink2/2.00a2.3

# MANUALLY NAME THIS
diet_outcome="alcohol_cirrhosis"
file="/pl/active/colelab/users/kjames/refinedMR/interim_data/MVMR/${diet_outcome}/all_clumped_snps_${diet_outcome}_MVMR_GIclump.txt"
out_path="/pl/active/colelab/users/kjames/refinedMR/interim_data/MVMR/${diet_outcome}"

# Define this folder (doesn't change)
GI_clumped="GI_clumped"

# Define where the log files with the missing SNP IDs are
logdir="/pl/active/colelab/users/kjames/refinedMR/interim_data/MVMR/${diet_outcome}/${GI_clumped}"

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



# Second, grep them from the original file; this preserves the heading
#(head -n 1 $file && grep -Ff $out_path/combined_missing_snps.txt $file | tail -n +2) > $out_path/filtered_missing_snps.txt

# Third, flip the alleles and rerun clumping
#outfile="all_missing_snps_flipped.txt"    # final combined output file

# NEED TO SWAP THE ORIGINAL FILE WITH P VALUES ETC
# Combine all .missing_snps.txt files, remove duplicates, flip REF and ALT
#cat "$logdir"/*.missing_snps.txt | sort -u | \
#awk -F':' 'BEGIN {OFS=":"} {print $1, $2, $4, $3}' > "$out_path/${outfile}"

# Use plink v1.9 for this; just type plink instead of plink2
for chr in {1..22}; do
  plink \
    --bfile /pl/active/colelab/users/kjames/refinedMR/interim_data/1KGP3_HG19/1KGP3_HG19_files_processed_EUR_Jenkai_chr${chr}_chrbpID \
    --clump "$out_path/flipped_summary_file.txt" \
    --clump-snp-field SNP \
    --clump-field P \
    --clump-p1 5e-8 \
    --clump-r2 0.01 \
    --clump-kb 10000 \
    --out /pl/active/colelab/users/kjames/refinedMR/interim_data/MVMR/${diet_outcome}/${GI_clumped}/GI_clumped_results_${diet_outcome}_chr${chr}.tsv
done



