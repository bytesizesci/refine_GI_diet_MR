#!/bin/bash

# Define the directory and output file
results_dir="/pl/active/colelab/users/kjames/refinedMR/interim_data/ldsc/results"  # Directory with results .log files
output_file="/pl/active/colelab/users/kjames/refinedMR/interim_data/ldsc/processed_results"

# Ensure the output file is empty
> "$output_file"

# Loop through all files with "ldsc_results" in their names
for file in "$results_dir"/ldsc_results*; do
  # Check if the file exists (in case of no matches)
  if [ -f "$file" ]; then
    echo "Processing $file..."
    grep -A90 'p1' "$file" | grep '/pl' >> $output_file/all_results.txt
  fi
done

echo "Results saved to $output_file"

# Replace the specified path with an empty string in the output file
sed 's|/pl/active/colelab/common/published_gwas/Neale_UKB_GWAS_round2/Both_sexes/LDSC/||g' $output_file/all_results.txt > $output_file/temp_results.txt

# Overwrite the original output file with the modified results
mv $output_file/temp_results.txt $output_file/all_results.txt

