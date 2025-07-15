#!/bin/bash
#SBATCH --job-name=05_plink_clump
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

# Goal: Use plink clumping to find independent loci. Results will be
# used in MVMR.

# ===================== USER CONFIGURATION =====================
# Define key directories
DATA_DIR= ~/working_example/data #"/pl/active/colelab/users/kjames/refinedMR/interim_data"
REF_DIR="${DATA_DIR}/reference_1KGP3_HG19"
GWAS_DIR="${DATA_DIR}/GWAS_NL" #this is created below
OUT_DIR= ~/working_example/results/plink_clumped #this is created below
DIET_OUTCOME="alcohol_ALT" #diet_outcome

# List of traits (one per line)
# Defined from 04_snp_overlap.R
TRAITS=(
  "30080_irnt"
  "30270_irnt"
  "23098_irnt"
  "20015_irnt"
  "23106_irnt"
)

# Load plink
# Customized to Alpine HPC
# Modify for your system
module load plink2/2.00a2.3

# ==============================================================
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Step 1: Get the GWAS NL files
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Create folders if they doesn't exist
mkdir -p "${GWAS_DIR}"
mkdir -p "${OUT_DIR}"

# See https://www.nealelab.is/uk-biobank for source files
BASE_URL="https://broad-ukb-sumstats-us-east-1.s3.amazonaws.com/round2/additive-tsvs"

for trait in "${TRAITS[@]}"; do
  url="${BASE_URL}/${trait}.gwas.imputed_v3.both_sexes.tsv.bgz"
  output_path="${GWAS_DIR}/${trait}.gwas.imputed_v3.both_sexes.tsv.bgz"
  echo "Downloading: $trait to $output_path"
  wget "$url" -O "$output_path
done

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Step 2: Unzip the reference files if needed
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Unzip the reference files
for file in ~/working_example/data/reference_1KGP3_HG19/*.gz; do
  [ -e "$file" ] || continue  # Skip loop if no .bgz files found
  gunzip "$file"
done

# Unzip them, remove the .bgz file after unzipping it
for file in ~/working_example/data/GWAS_NL/*.bgz; do
  [ -e "$file" ] || continue  # Skip loop if no .bgz files found
  zcat "$file" > "${file%.bgz}" && rm "$file"
done

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Step 3: Run plink clumping
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
for trait in "${TRAITS[@]}"; do
  gwas_file="${GWAS_DIR}/${trait}.gwas.imputed_v3.both_sexes.tsv"

  trait_outdir="${OUT_DIR}/${DIET_OUTCOME}/${trait}"
  mkdir -p "${trait_outdir}"

  for chr in {1..22}; do
    plink \
      --bfile "${REF_DIR}/1KGP3_HG19_files_processed_EUR_Jenkai_chr${chr}_chrbpID" \
      --clump "${gwas_file}" \
      --clump-snp-field variant \
      --clump-field pval \
      --clump-p1 5e-8 \
      --clump-kb 500 \
      --clump-r2 0.1 \
      --out "${trait_outdir}/clumped_results_${trait}_chr${chr}.tsv"
  done
done


