#!/bin/bash
#SBATCH --job-name=ldsc_rg
#SBATCH --output=config/ldsc_rg_%j.out                  
#SBATCH --error=config/ldsc_rg_%j.err                   
#SBATCH --partition=your_partition_name                 # CHANGE: partition name for your system
#SBATCH --nodes=1
#SBATCH --qos=long                                      # REMOVE or modify if your system doesn't use QoS
#SBATCH --cpus-per-task=1
#SBATCH --mail-type=ALL
#SBATCH --mail-user=your.email@domain.edu               # CHANGE: add your email address
#SBATCH --time=144:00:00                                # Max run time. Specify the time limit for your allocation
#SBATCH --mem=4G                                        # Memory per job

# Make the script executable with: chmod +x 01_ldsc_rg.sh
# Run the script: sbatch 01_ldsc_rg.sh

# ===================== USER CONFIGURATION =====================
# Step 1: Define paths 

# Input path (update as needed)
in_dir="working_example/data" #/pl/active/colelab/users/kjames/refinedMR/interim_data/ldsc/PheWAS_Ttest_LDSC_input

# Contains ldsc file paths, one per line
# Note, depending on where your ldsc files are saved, may have to modify the content in LDSC_filepaths to correctly reflect the location
file_list="$in_dir/LDSC_filepaths.txt" #$in_dir/PheWAS_Ttest_ALLFILES_for_gcor.txt" 

# Directory to save results
# This is created by this script
output_dir="working_example/ldsc_results/rg" #"/pl/active/colelab/users/kjames/refinedMR/interim_data/ldsc/PheWAS_Ttest_ldsc_results/all_gcor_files" 

# Reference directory for LDSC
ref_dir="working_example/data/reference_1KGP3_HG19" #"/pl/active/colelab/common/eur_w_ld_chr/" 
# ==============================================================

# Step 2: Load LDSC
# Customized to Alpine HPC, uncomment for Alpine
#module use --append /pl/active/colelab/common/software/lmod-files
# Modify for your system 
module load ldsc

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
