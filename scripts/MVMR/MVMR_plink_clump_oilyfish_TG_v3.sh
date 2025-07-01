# KJS
# March 25, 2025

# Picking up from MVMR_plink_clump_oilyfish_TG.sh
# That script copied and unzipped 8 files to include in the MVMR
# It also ran plink clumping on the trait 4124_irnt
# I'm now working on automating the plink clumping process. I don't need
# to recopy the files over so I'm working on putting the clumping in this script.

# To run: Be in /pl/active/colelab/users/kjames/refinedMR/scripts/MVMR
# Make sure to be on a worker node - acompile
# Might have to change permissions the first time - chmod +x <scriptname>
# ./MVMR_plink_clump_oilyfish_TG_v3.sh

# Get plink
module load plink2/2.00a2.3

# Files with oilyfish - TG to clump with plink
#/pl/active/colelab/users/kjames/refinedMR/interim_data/ldsc/NL_bothsexes_temp/4124_irnt.gwas.imputed_v3.both_sexes.tsv
#/pl/active/colelab/users/kjames/refinedMR/interim_data/ldsc/NL_bothsexes_temp/30100_irnt.gwas.imputed_v3.both_sexes.tsv
#/pl/active/colelab/users/kjames/refinedMR/interim_data/ldsc/NL_bothsexes_temp/4100_irnt.gwas.imputed_v3.both_sexes.tsv
#/pl/active/colelab/users/kjames/refinedMR/interim_data/ldsc/NL_bothsexes_temp/50_irnt.gwas.imputed_v3.both_sexes.tsv
#/pl/active/colelab/users/kjames/refinedMR/interim_data/ldsc/NL_bothsexes_temp/5084_irnt.gwas.imputed_v3.both_sexes.tsv
#/pl/active/colelab/users/kjames/refinedMR/interim_data/ldsc/NL_bothsexes_temp/23104_irnt.gwas.imputed_v3.both_sexes.tsv
#/pl/active/colelab/users/kjames/refinedMR/interim_data/ldsc/NL_bothsexes_temp/23106_irnt.gwas.imputed_v3.both_sexes.tsv
#/pl/active/colelab/users/kjames/refinedMR/interim_data/ldsc/NL_bothsexes_temp/23101_irnt.gwas.imputed_v3.both_sexes.tsv

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 1.
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Clump the gwas results for a given trait
# Define trait
file="/pl/active/colelab/users/kjames/refinedMR/interim_data/ldsc/NL_bothsexes_temp/4124_irnt.gwas.imputed_v3.both_sexes.tsv"
trait="4124_irnt"

if [ ! -d "/pl/active/colelab/users/kjames/refinedMR/interim_data/ldsc/plink_clumped_1KGP_HG19/oilyfish_tg/${trait}" ]; then
  mkdir /pl/active/colelab/users/kjames/refinedMR/interim_data/ldsc/plink_clumped_1KGP_HG19/oilyfish_tg/${trait}
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
    --out /pl/active/colelab/users/kjames/refinedMR/interim_data/ldsc/plink_clumped_1KGP_HG19/oilyfish_tg/${trait}/clumped_results_${trait}_chr${chr}.tsv
done
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 2.
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Clump the gwas results for a given trait
# Define trait
file="/pl/active/colelab/users/kjames/refinedMR/interim_data/ldsc/NL_bothsexes_temp/30100_irnt.gwas.imputed_v3.both_sexes.tsv"
trait="30100_irnt"

if [ ! -d "/pl/active/colelab/users/kjames/refinedMR/interim_data/ldsc/plink_clumped_1KGP_HG19/oilyfish_tg/${trait}" ]; then
  mkdir /pl/active/colelab/users/kjames/refinedMR/interim_data/ldsc/plink_clumped_1KGP_HG19/oilyfish_tg/${trait}
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
    --out /pl/active/colelab/users/kjames/refinedMR/interim_data/ldsc/plink_clumped_1KGP_HG19/oilyfish_tg/${trait}/clumped_results_${trait}_chr${chr}.tsv
done
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 3.
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
file="/pl/active/colelab/users/kjames/refinedMR/interim_data/ldsc/NL_bothsexes_temp/4100_irnt.gwas.imputed_v3.both_sexes.tsv"
trait="4100_irnt"

if [ ! -d "/pl/active/colelab/users/kjames/refinedMR/interim_data/ldsc/plink_clumped_1KGP_HG19/oilyfish_tg/${trait}" ]; then
  mkdir /pl/active/colelab/users/kjames/refinedMR/interim_data/ldsc/plink_clumped_1KGP_HG19/oilyfish_tg/${trait}
fi

for chr in {1..22}; do
  plink \
    --bfile /pl/active/colelab/users/kjames/refinedMR/interim_data/1KGP3_HG19/1KGP3_HG19_files_processed_EUR_Jenkai_chr${chr}_chrbpID \
    --clump ${file} \
    --clump-snp-field variant \
    --clump-field pval \
    --clump-p1 5e-8 \
    --clump-kb 500 \
    --clump-r2 0.1 \
    --out /pl/active/colelab/users/kjames/refinedMR/interim_data/ldsc/plink_clumped_1KGP_HG19/oilyfish_tg/${trait}/clumped_results_${trait}_chr${chr}.tsv
done
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 4.
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
file="/pl/active/colelab/users/kjames/refinedMR/interim_data/ldsc/NL_bothsexes_temp/50_irnt.gwas.imputed_v3.both_sexes.tsv"
trait="50_irnt"

if [ ! -d "/pl/active/colelab/users/kjames/refinedMR/interim_data/ldsc/plink_clumped_1KGP_HG19/oilyfish_tg/${trait}" ]; then
  mkdir /pl/active/colelab/users/kjames/refinedMR/interim_data/ldsc/plink_clumped_1KGP_HG19/oilyfish_tg/${trait}
fi

for chr in {1..22}; do
  plink \
    --bfile /pl/active/colelab/users/kjames/refinedMR/interim_data/1KGP3_HG19/1KGP3_HG19_files_processed_EUR_Jenkai_chr${chr}_chrbpID \
    --clump ${file} \
    --clump-snp-field variant \
    --clump-field pval \
    --clump-p1 5e-8 \
    --clump-kb 500 \
    --clump-r2 0.1 \
    --out /pl/active/colelab/users/kjames/refinedMR/interim_data/ldsc/plink_clumped_1KGP_HG19/oilyfish_tg/${trait}/clumped_results_${trait}_chr${chr}.tsv
done
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 5.
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
file="/pl/active/colelab/users/kjames/refinedMR/interim_data/ldsc/NL_bothsexes_temp/5084_irnt.gwas.imputed_v3.both_sexes.tsv"
trait="5084_irnt"

if [ ! -d "/pl/active/colelab/users/kjames/refinedMR/interim_data/ldsc/plink_clumped_1KGP_HG19/oilyfish_tg/${trait}" ]; then
  mkdir /pl/active/colelab/users/kjames/refinedMR/interim_data/ldsc/plink_clumped_1KGP_HG19/oilyfish_tg/${trait}
fi

for chr in {1..22}; do
  plink \
    --bfile /pl/active/colelab/users/kjames/refinedMR/interim_data/1KGP3_HG19/1KGP3_HG19_files_processed_EUR_Jenkai_chr${chr}_chrbpID \
    --clump ${file} \
    --clump-snp-field variant \
    --clump-field pval \
    --clump-p1 5e-8 \
    --clump-kb 500 \
    --clump-r2 0.1 \
    --out /pl/active/colelab/users/kjames/refinedMR/interim_data/ldsc/plink_clumped_1KGP_HG19/oilyfish_tg/${trait}/clumped_results_${trait}_chr${chr}.tsv
done
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 6.
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
file="/pl/active/colelab/users/kjames/refinedMR/interim_data/ldsc/NL_bothsexes_temp/23104_irnt.gwas.imputed_v3.both_sexes.tsv"
trait="23104_irnt"

if [ ! -d "/pl/active/colelab/users/kjames/refinedMR/interim_data/ldsc/plink_clumped_1KGP_HG19/oilyfish_tg/${trait}" ]; then
  mkdir /pl/active/colelab/users/kjames/refinedMR/interim_data/ldsc/plink_clumped_1KGP_HG19/oilyfish_tg/${trait}
fi

for chr in {1..22}; do
  plink \
    --bfile /pl/active/colelab/users/kjames/refinedMR/interim_data/1KGP3_HG19/1KGP3_HG19_files_processed_EUR_Jenkai_chr${chr}_chrbpID \
    --clump ${file} \
    --clump-snp-field variant \
    --clump-field pval \
    --clump-p1 5e-8 \
    --clump-kb 500 \
    --clump-r2 0.1 \
    --out /pl/active/colelab/users/kjames/refinedMR/interim_data/ldsc/plink_clumped_1KGP_HG19/oilyfish_tg/${trait}/clumped_results_${trait}_chr${chr}.tsv
done
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 7.
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
file="/pl/active/colelab/users/kjames/refinedMR/interim_data/ldsc/NL_bothsexes_temp/23106_irnt.gwas.imputed_v3.both_sexes.tsv"
trait="23106_irnt"

if [ ! -d "/pl/active/colelab/users/kjames/refinedMR/interim_data/ldsc/plink_clumped_1KGP_HG19/oilyfish_tg/${trait}" ]; then
  mkdir /pl/active/colelab/users/kjames/refinedMR/interim_data/ldsc/plink_clumped_1KGP_HG19/oilyfish_tg/${trait}
fi

for chr in {1..22}; do
  plink \
    --bfile /pl/active/colelab/users/kjames/refinedMR/interim_data/1KGP3_HG19/1KGP3_HG19_files_processed_EUR_Jenkai_chr${chr}_chrbpID \
    --clump ${file} \
    --clump-snp-field variant \
    --clump-field pval \
    --clump-p1 5e-8 \
    --clump-kb 500 \
    --clump-r2 0.1 \
    --out /pl/active/colelab/users/kjames/refinedMR/interim_data/ldsc/plink_clumped_1KGP_HG19/oilyfish_tg/${trait}/clumped_results_${trait}_chr${chr}.tsv
done
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 8.
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
file="/pl/active/colelab/users/kjames/refinedMR/interim_data/ldsc/NL_bothsexes_temp/23101_irnt.gwas.imputed_v3.both_sexes.tsv"
trait="23101_irnt"

if [ ! -d "/pl/active/colelab/users/kjames/refinedMR/interim_data/ldsc/plink_clumped_1KGP_HG19/oilyfish_tg/${trait}" ]; then
  mkdir /pl/active/colelab/users/kjames/refinedMR/interim_data/ldsc/plink_clumped_1KGP_HG19/oilyfish_tg/${trait}
fi

for chr in {1..22}; do
  plink \
    --bfile /pl/active/colelab/users/kjames/refinedMR/interim_data/1KGP3_HG19/1KGP3_HG19_files_processed_EUR_Jenkai_chr${chr}_chrbpID \
    --clump ${file} \
    --clump-snp-field variant \
    --clump-field pval \
    --clump-p1 5e-8 \
    --clump-kb 500 \
    --clump-r2 0.1 \
    --out /pl/active/colelab/users/kjames/refinedMR/interim_data/ldsc/plink_clumped_1KGP_HG19/oilyfish_tg/${trait}/clumped_results_${trait}_chr${chr}.tsv
done