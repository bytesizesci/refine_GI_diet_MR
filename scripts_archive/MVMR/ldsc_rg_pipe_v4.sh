#!/bin/bash
#SBATCH --job-name=LDSC_bothsexes
#SBATCH --output=/pl/active/colelab/users/kjames/refinedMR/config/download_job_%j.out
#SBATCH --error=/pl/active/colelab/users/kjames/refinedMR/config/download_job_%j.err
#SBATCH --partition=amilan
#SBATCH --nodes=1
#SBATCH --qos=long
#SBATCH --cpus-per-task=1
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kristen.sutton@cuanschutz.edu
#SBATCH --time=144:00:00  # Specify the time limit for your allocation; if qos=long then can have 7 days
#SBATCH --mem=4G        # Adjust the memory limit as per your needs

# Make the script executable with: chmod +x filename.sh
# Run the script: sbatch filename.sh

#load ldsc module
module use --append /pl/active/colelab/common/software/lmod-files
module load ldsc


# Step 1: Define the text file containing the file paths
in_dir=/pl/active/colelab/users/kjames/refinedMR/interim_data/ldsc/PheWAS_Ttest_LDSC_input

# PICK THE CORRECT FILE!!!
file_list="$in_dir/PheWAS_Ttest_ALLFILES_for_gcor.txt" # Contains ldsc file paths, one per line
# Submitted tracker: 
# whiteXwholegrain_LDL
# oilyfish_CVD

#~~~~~~~~~~~~~~~~~~~~~
# "$in_dir/oilyfish_GCST90239664_TG_Graham_descending_TTest_filenames.txt"
# "$in_dir/oilyfish_GCST90132314_CVD_Aragam_descending_TTest_filenames.txt"
# "$in_dir/oilyfish_GCST006901_Height_Yengo_2018_descending_TTest_filenames.txt"

# "$in_dir/whiteXwholegrain_GCST90239658_LDL_Graham_descending_TTest_filenames.txt"
# "$in_dir/whiteXwholegrain_GCST90132314_CVD_Aragam_descending_TTest_filenames.txt"
# "$in_dir/whiteXwholegrain_GCST006901_Height_Yengo_2018_descending_TTest_filenames.txt"

# "$in_dir/alcohol_GCST90013405_ALT_Pazoki_descending_TTest_filenames.txt"
# "$in_dir/alcohol_GCST90319877_cir_Ghouse_descending_TTest_filenames.txt"
# "$in_dir/alcohol_GCST006901_Height_Yengo_2018_descending_TTest_filenames.txt"

# Contains all unique files in one!
# "$in_dir/PheWAS_Ttest_ALLFILES_for_gcor.txt"
#~~~~~~~~~~~~~~~~~~~~~

# Step 2: Define output directory and reference directory
# NAME THE CORRECt OUTPUT LOCATION!!!
# SHOULD CORRESPOND TO file_list VARIABLE!!
# Directory to save results
output_dir="/pl/active/colelab/users/kjames/refinedMR/interim_data/ldsc/PheWAS_Ttest_ldsc_results/all_gcor_files"  

#~~~~~~~~~~~~~~~~~~~~~
# "/pl/active/colelab/users/kjames/refinedMR/interim_data/ldsc/PheWAS_Ttest_ldsc_results/oilyfish_TG" 
# "/pl/active/colelab/users/kjames/refinedMR/interim_data/ldsc/PheWAS_Ttest_ldsc_results/oilyfish_CVD" 
# "/pl/active/colelab/users/kjames/refinedMR/interim_data/ldsc/PheWAS_Ttest_ldsc_results/oilyfish_Height" 

# "/pl/active/colelab/users/kjames/refinedMR/interim_data/ldsc/PheWAS_Ttest_ldsc_results/whiteXwholegrain_LDL" 
# "/pl/active/colelab/users/kjames/refinedMR/interim_data/ldsc/PheWAS_Ttest_ldsc_results/whiteXwholegrain_CVD" 
# "/pl/active/colelab/users/kjames/refinedMR/interim_data/ldsc/PheWAS_Ttest_ldsc_results/whiteXwholegrain_Height" 

# "/pl/active/colelab/users/kjames/refinedMR/interim_data/ldsc/PheWAS_Ttest_ldsc_results/alcohol_ALT" 
# "/pl/active/colelab/users/kjames/refinedMR/interim_data/ldsc/PheWAS_Ttest_ldsc_results/alcohol_cir" 
# "/pl/active/colelab/users/kjames/refinedMR/interim_data/ldsc/PheWAS_Ttest_ldsc_results/alcohol_Height" 
#~~~~~~~~~~~~~~~~~~~~~

# Reference directory for LDSC
ref_dir="/pl/active/colelab/common/eur_w_ld_chr/" 

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
