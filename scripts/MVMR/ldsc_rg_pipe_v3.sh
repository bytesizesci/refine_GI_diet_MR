#!/bin/bash
# Make the script executable with: chmod +x run_ldsc.sh
# Run the script: ./run_ldsc.sh

#load ldsc module
module use --append /pl/active/colelab/common/software/lmod-files
module load ldsc


# Step 1: Define the text file containing the file paths
in_dir=/pl/active/colelab/users/kjames/refinedMR/interim_data/ldsc/PheWAS_Ttest_LDSC_input

file_list="$in_dir/oilyfish_GCST90239664_TG_Graham_descending_TTest_filenames.txt" # Contains ldsc file paths, one per line

# Step 2: Define output directory and reference directory
output_dir="/pl/active/colelab/users/kjames/refinedMR/interim_data/ldsc/PheWAS_Ttest_ldsc_results/oilyfish_tg"  # Directory to save results
ref_dir="/pl/active/colelab/common/eur_w_ld_chr/" # Reference directory for LDSC

# Step 3: Ensure output directory exists
mkdir -p $output_dir  # Create the output directory if it doesn't exist

# Step 4: Read the file paths into an array
mapfile -t file_paths < $file_list
# `mapfile` reads the file paths from `file_list` into an array named `file_paths`
# Each row in the file becomes an element in the array (named file_list)
# to look at the array: echo "${file_paths[@]}"

# Step 5: Get the total number of files
total_files=${#file_paths[@]}  # ${#file_paths[@]} gives the number of elements in the array

# Step 6: Loop through the array to create subsets of files
for ((i=0; i<total_files-1; i++)); do
  # Create a comma-separated list of files starting from index `i`
  subset_files=$(IFS=,; echo "${file_paths[*]:i}") 
  # Explanation of the above line:
  # - IFS=, sets the internal field separator to a comma
  # - "${file_paths[@]:i}" extracts elements from index `i` to the end of the array
  # - `echo` joins the extracted elements with commas
  
  # Define output file name based on the starting file in the subset
  output_file="$output_dir/ldsc_results_starting_with_${i}.txt"

  # Step 7: Run the `ldsc.py` command
  ldsc.py \
    --rg $subset_files \
    --ref-ld-chr $ref_dir \
    --w-ld-chr $ref_dir \
    --out $output_file
  
  # Print status update
  echo "LDSC run complete for subset starting with file $((i+1))"
done

# Final status update
echo "All LDSC commands have been submitted."
