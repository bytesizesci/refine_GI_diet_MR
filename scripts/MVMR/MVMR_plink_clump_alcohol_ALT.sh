# Author: KJS
# Date: May 12, 2025
# Project: Refining diet MR

# Goal: Copy the NL GWAS results to local project from shared location.
# Unzip them. Then use plink clumping to find independent loci. Results will be
# used in MVMR.

# To run: Be in /pl/active/colelab/users/kjames/refinedMR/scripts/MVMR
# Make sure to be on a worker node - acompile (4 hrs)
# acompile --time=04:00:00
# Might have to change permissions the first time - chmod +x <scriptname>
# ./MVMR_plink_clump_alcohol_ALT.sh

# Get plink
module load plink2/2.00a2.3

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ONLY RUN THIS CHUNK ONCE - comment out after you have the files!
# Copy reference files to a temp directory
# Get chrpos format (not rsID) to match Neale Lab files
#cp /pl/active/colelab/common/1KGP3_geno/1KGP3_HG19_files_processed_EUR_Jenkai_chr*_chrbpID* /pl/active/colelab/users/kjames/refinedMR/interim_data/1KGP3_HG19/

# Unzip the reference files
#gunzip /pl/active/colelab/users/kjames/refinedMR/interim_data/1KGP3_HG19/*.gz
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Get the traits for the MVMR
# WILL BE TRAIT-OUTCOME SPECIFIC; MANUAL JOB
# These correspond to alcohol and ALT
# Coded to add them manually from dat3 object in gcor_PheWASttest_data_reduction script
cp /pl/active/colelab/common/published_gwas/Neale_UKB_GWAS_round2/Both_sexes/GWAS2/30080_irnt.gwas.imputed_v3.both_sexes.tsv.bgz /pl/active/colelab/users/kjames/refinedMR/interim_data/ldsc/NL_bothsexes_temp/
cp /pl/active/colelab/common/published_gwas/Neale_UKB_GWAS_round2/Both_sexes/GWAS2/30270_irnt.gwas.imputed_v3.both_sexes.tsv.bgz /pl/active/colelab/users/kjames/refinedMR/interim_data/ldsc/NL_bothsexes_temp/
cp /pl/active/colelab/common/published_gwas/Neale_UKB_GWAS_round2/Both_sexes/GWAS2/23098_irnt.gwas.imputed_v3.both_sexes.tsv.bgz /pl/active/colelab/users/kjames/refinedMR/interim_data/ldsc/NL_bothsexes_temp/
cp /pl/active/colelab/common/published_gwas/Neale_UKB_GWAS_round2/Both_sexes/GWAS2/20015_irnt.gwas.imputed_v3.both_sexes.tsv.bgz /pl/active/colelab/users/kjames/refinedMR/interim_data/ldsc/NL_bothsexes_temp/
cp /pl/active/colelab/common/published_gwas/Neale_UKB_GWAS_round2/Both_sexes/GWAS2/23106_irnt.gwas.imputed_v3.both_sexes.tsv.bgz /pl/active/colelab/users/kjames/refinedMR/interim_data/ldsc/NL_bothsexes_temp/

# Unzip them, remove the .bgz file after unzipping it
for file in /pl/active/colelab/users/kjames/refinedMR/interim_data/ldsc/NL_bothsexes_temp/*.bgz; do
  zcat "$file" > "${file%.bgz}" && rm "$file"
done
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


# Currently very manual - for 1 trait at a time...
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 1.
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Clump the gwas results for a given trait
# MANUALLY DEFINE TRAIT AND DIET_OUTCOME
file="/pl/active/colelab/users/kjames/refinedMR/interim_data/ldsc/NL_bothsexes_temp/30080_irnt.gwas.imputed_v3.both_sexes.tsv"
trait="30080_irnt"
diet_outcome="alcohol_ALT"

# Make a folder for the diet_outcome if it doesn't exist
if [ ! -d "/pl/active/colelab/users/kjames/refinedMR/interim_data/ldsc/plink_clumped_1KGP_HG19/${diet_outcome}" ]; then
  mkdir /pl/active/colelab/users/kjames/refinedMR/interim_data/ldsc/plink_clumped_1KGP_HG19/${diet_outcome}
fi

# Within the diet_outcome folder, make a folder for the trait if it doesn't exist
if [ ! -d "/pl/active/colelab/users/kjames/refinedMR/interim_data/ldsc/plink_clumped_1KGP_HG19/${diet_outcome}/${trait}" ]; then
  mkdir /pl/active/colelab/users/kjames/refinedMR/interim_data/ldsc/plink_clumped_1KGP_HG19/${diet_outcome}/${trait}
fi

# Use plink v1.9 for this; just type plink instead of plink2
for chr in {1..22}; do
  plink \
    --bfile /pl/active/colelab/users/kjames/refinedMR/interim_data/1KGP3_HG19/1KGP3_HG19_files_processed_EUR_Jenkai_chr${chr}_chrbpID \
    --clump ${file} \
    --clump-snp-field variant \
    --clump-field pval \
    --clump-p1 5e-8 \
    --clump-kb 500 \
    --clump-r2 0.1 \
    --out /pl/active/colelab/users/kjames/refinedMR/interim_data/ldsc/plink_clumped_1KGP_HG19/${diet_outcome}/${trait}/clumped_results_${trait}_chr${chr}.tsv
done

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 2.
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Clump the gwas results for a given trait
# Manually define trait
file="/pl/active/colelab/users/kjames/refinedMR/interim_data/ldsc/NL_bothsexes_temp/30270_irnt.gwas.imputed_v3.both_sexes.tsv"
trait="30270_irnt"
diet_outcome="alcohol_ALT"

# Make a folder for the diet_outcome if it doesn't exist
if [ ! -d "/pl/active/colelab/users/kjames/refinedMR/interim_data/ldsc/plink_clumped_1KGP_HG19/${diet_outcome}" ]; then
  mkdir /pl/active/colelab/users/kjames/refinedMR/interim_data/ldsc/plink_clumped_1KGP_HG19/${diet_outcome}
fi

# Within the diet_outcome folder, make a folder for the trait if it doesn't exist
if [ ! -d "/pl/active/colelab/users/kjames/refinedMR/interim_data/ldsc/plink_clumped_1KGP_HG19/${diet_outcome}/${trait}" ]; then
  mkdir /pl/active/colelab/users/kjames/refinedMR/interim_data/ldsc/plink_clumped_1KGP_HG19/${diet_outcome}/${trait}
fi

# Use plink v1.9 for this; just type plink instead of plink2
for chr in {1..22}; do
  plink \
    --bfile /pl/active/colelab/users/kjames/refinedMR/interim_data/1KGP3_HG19/1KGP3_HG19_files_processed_EUR_Jenkai_chr${chr}_chrbpID \
    --clump ${file} \
    --clump-snp-field variant \
    --clump-field pval \
    --clump-p1 5e-8 \
    --clump-kb 500 \
    --clump-r2 0.1 \
    --out /pl/active/colelab/users/kjames/refinedMR/interim_data/ldsc/plink_clumped_1KGP_HG19/${diet_outcome}/${trait}/clumped_results_${trait}_chr${chr}.tsv
done

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 3.
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Clump the gwas results for a given trait
# Manually define trait
file="/pl/active/colelab/users/kjames/refinedMR/interim_data/ldsc/NL_bothsexes_temp/23098_irnt.gwas.imputed_v3.both_sexes.tsv"
trait="23098_irnt"
diet_outcome="alcohol_ALT"

# Make a folder for the diet_outcome if it doesn't exist
if [ ! -d "/pl/active/colelab/users/kjames/refinedMR/interim_data/ldsc/plink_clumped_1KGP_HG19/${diet_outcome}" ]; then
  mkdir /pl/active/colelab/users/kjames/refinedMR/interim_data/ldsc/plink_clumped_1KGP_HG19/${diet_outcome}
fi

# Within the diet_outcome folder, make a folder for the trait if it doesn't exist
if [ ! -d "/pl/active/colelab/users/kjames/refinedMR/interim_data/ldsc/plink_clumped_1KGP_HG19/${diet_outcome}/${trait}" ]; then
  mkdir /pl/active/colelab/users/kjames/refinedMR/interim_data/ldsc/plink_clumped_1KGP_HG19/${diet_outcome}/${trait}
fi

# Use plink v1.9 for this; just type plink instead of plink2
for chr in {1..22}; do
  plink \
    --bfile /pl/active/colelab/users/kjames/refinedMR/interim_data/1KGP3_HG19/1KGP3_HG19_files_processed_EUR_Jenkai_chr${chr}_chrbpID \
    --clump ${file} \
    --clump-snp-field variant \
    --clump-field pval \
    --clump-p1 5e-8 \
    --clump-kb 500 \
    --clump-r2 0.1 \
    --out /pl/active/colelab/users/kjames/refinedMR/interim_data/ldsc/plink_clumped_1KGP_HG19/${diet_outcome}/${trait}/clumped_results_${trait}_chr${chr}.tsv
done

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 4.
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Clump the gwas results for a given trait
# Manually define trait
file="/pl/active/colelab/users/kjames/refinedMR/interim_data/ldsc/NL_bothsexes_temp/20015_irnt.gwas.imputed_v3.both_sexes.tsv"
trait="20015_irnt"
diet_outcome="alcohol_ALT"

# Make a folder for the diet_outcome if it doesn't exist
if [ ! -d "/pl/active/colelab/users/kjames/refinedMR/interim_data/ldsc/plink_clumped_1KGP_HG19/${diet_outcome}" ]; then
  mkdir /pl/active/colelab/users/kjames/refinedMR/interim_data/ldsc/plink_clumped_1KGP_HG19/${diet_outcome}
fi

# Within the diet_outcome folder, make a folder for the trait if it doesn't exist
if [ ! -d "/pl/active/colelab/users/kjames/refinedMR/interim_data/ldsc/plink_clumped_1KGP_HG19/${diet_outcome}/${trait}" ]; then
  mkdir /pl/active/colelab/users/kjames/refinedMR/interim_data/ldsc/plink_clumped_1KGP_HG19/${diet_outcome}/${trait}
fi

# Use plink v1.9 for this; just type plink instead of plink2
for chr in {1..22}; do
  plink \
    --bfile /pl/active/colelab/users/kjames/refinedMR/interim_data/1KGP3_HG19/1KGP3_HG19_files_processed_EUR_Jenkai_chr${chr}_chrbpID \
    --clump ${file} \
    --clump-snp-field variant \
    --clump-field pval \
    --clump-p1 5e-8 \
    --clump-kb 500 \
    --clump-r2 0.1 \
    --out /pl/active/colelab/users/kjames/refinedMR/interim_data/ldsc/plink_clumped_1KGP_HG19/${diet_outcome}/${trait}/clumped_results_${trait}_chr${chr}.tsv
done

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 5.
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Clump the gwas results for a given trait
# Manually define trait
file="/pl/active/colelab/users/kjames/refinedMR/interim_data/ldsc/NL_bothsexes_temp/23106_irnt.gwas.imputed_v3.both_sexes.tsv"
trait="23106_irnt"
diet_outcome="alcohol_ALT"

# Make a folder for the diet_outcome if it doesn't exist
if [ ! -d "/pl/active/colelab/users/kjames/refinedMR/interim_data/ldsc/plink_clumped_1KGP_HG19/${diet_outcome}" ]; then
  mkdir /pl/active/colelab/users/kjames/refinedMR/interim_data/ldsc/plink_clumped_1KGP_HG19/${diet_outcome}
fi

# Within the diet_outcome folder, make a folder for the trait if it doesn't exist
if [ ! -d "/pl/active/colelab/users/kjames/refinedMR/interim_data/ldsc/plink_clumped_1KGP_HG19/${diet_outcome}/${trait}" ]; then
  mkdir /pl/active/colelab/users/kjames/refinedMR/interim_data/ldsc/plink_clumped_1KGP_HG19/${diet_outcome}/${trait}
fi

# Use plink v1.9 for this; just type plink instead of plink2
for chr in {1..22}; do
  plink \
    --bfile /pl/active/colelab/users/kjames/refinedMR/interim_data/1KGP3_HG19/1KGP3_HG19_files_processed_EUR_Jenkai_chr${chr}_chrbpID \
    --clump ${file} \
    --clump-snp-field variant \
    --clump-field pval \
    --clump-p1 5e-8 \
    --clump-kb 500 \
    --clump-r2 0.1 \
    --out /pl/active/colelab/users/kjames/refinedMR/interim_data/ldsc/plink_clumped_1KGP_HG19/${diet_outcome}/${trait}/clumped_results_${trait}_chr${chr}.tsv
done
