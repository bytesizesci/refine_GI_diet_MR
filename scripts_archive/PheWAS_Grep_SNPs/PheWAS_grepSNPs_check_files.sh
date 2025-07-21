#!/bin/bash

# Submit from enviroMR with the next line of code
# from refinedMR
# ./scripts/PheWAS_Grep_SNPs/PheWAS_grepSNPs_check_files.sh

# Define the directory containing the files
#directory="/pl/active/colelab/users/kjames/enviroMR/results_array_dietPA"
#directory="/pl/active/colelab/users/kjames/refinedMR/interim_data/PheWAS_Grep_SNPs/PheWAS_SNPs_Oilyfish_WhiteXWhole_Alcochol_UpdatedProxies"
directory="/pl/active/colelab/users/kjames/refinedMR/interim_data/PheWAS_Grep_SNPs/PheWAS_SNPs_Height_Proxies"

# Create an array with all expected file numbers
expected_numbers=($(seq 1 4587))

# Create an array to store missing file numbers
missing_numbers=()

# Iterate over each expected file number
for number in "${expected_numbers[@]}"; do
    # Check if the file exists
    if [ ! -f "$directory/PheWAS_SNPs_Height_Proxies_${number}.csv" ]; then
        missing_numbers+=($number)
    fi
done

# Check if any files are missing
if [ ${#missing_numbers[@]} -eq 0 ]; then
    echo "All files are present."
else
    # Write the missing file numbers to a text file
    printf "%s\n" "${missing_numbers[@]}" > temp/missing_files_PheWAS_SNPs_Height_Proxies.txt
    echo "Missing files have been written to temp/missing_files_PheWAS_SNPs_Height_Proxies.txt."
fi