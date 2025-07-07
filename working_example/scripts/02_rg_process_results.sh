#!/bin/bash

# Make the script executable with: chmod +x 02_rg_process_results.sh
# Run the script: ./02_rg_process_results.sh

# ===================== USER CONFIGURATION =====================
# Define the directory and output file
output_dir="/your/path/to/ldsc_results" #"/pl/active/colelab/users/kjames/refinedMR/interim_data/ldsc/PheWAS_Ttest_ldsc_results/all_gcor_files"  # Directory with results .log files
output_file="/your/path/to/processed_results" #"/pl/active/colelab/users/kjames/refinedMR/interim_data/ldsc/processed_results"

# Set the path to be removed (user should update this to match their system)
path_to_remove="/your/full/path/to/LDSC/files/" #/pl/active/colelab/common/published_gwas/Neale_UKB_GWAS_round2/Both_sexes/LDSC/
# ==============================================================

# Ensure the output file is empty
> "$output_file"

# Loop through all files with "ldsc_results" in their names
for file in "$output_dir"/ldsc_results*; do
  # Check if the file exists (in case of no matches)
  if [ -f "$file" ]; then
    echo "Processing $file..."
    grep -A122 'p1' "$file" | grep '/pl' >> $output_file/gcor_PheWASttest_results_all.txt
  fi
done

echo "Results saved to $output_file"

# Clean up file names
# Replace the specified path with an empty string
sed "s|$path_to_remove||g" "$output_file/gcor_PheWASttest_results_all.txt" > "$output_file/temp_results.txt"

# Overwrite the original output file with the modified results
mv $output_file/temp_results.txt $output_file/gcor_PheWASttest_results_all.txt
