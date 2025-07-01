





# Get plink
module load plink2/2.00a2.3

# Copy reference files to a temp directory
# Get chrpos format (not rsID) to match Neale Lab files
cp /pl/active/colelab/common/1KGP3_geno/1KGP3_HG19_files_processed_EUR_Jenkai_chr*_chrbpID* /pl/active/colelab/users/kjames/refinedMR/interim_data/1KGP3_HG19/

# Unzip the reference files
gunzip /pl/active/colelab/users/kjames/refinedMR/interim_data/1KGP3_HG19/*.gz

# Get the traits for the MVMR
# These correspond to oily fish and TG
# Coded to add them manually from dat3 object in gcor_PheWASttest_data_reduction script
cp /pl/active/colelab/common/published_gwas/Neale_UKB_GWAS_round2/Both_sexes/GWAS2/4124_irnt.gwas.imputed_v3.both_sexes.tsv.bgz /pl/active/colelab/users/kjames/refinedMR/interim_data/ldsc/NL_bothsexes_temp/
cp /pl/active/colelab/common/published_gwas/Neale_UKB_GWAS_round2/Both_sexes/GWAS2/30100_irnt.gwas.imputed_v3.both_sexes.tsv.bgz /pl/active/colelab/users/kjames/refinedMR/interim_data/ldsc/NL_bothsexes_temp/
cp /pl/active/colelab/common/published_gwas/Neale_UKB_GWAS_round2/Both_sexes/GWAS2/4100_irnt.gwas.imputed_v3.both_sexes.tsv.bgz /pl/active/colelab/users/kjames/refinedMR/interim_data/ldsc/NL_bothsexes_temp/
cp /pl/active/colelab/common/published_gwas/Neale_UKB_GWAS_round2/Both_sexes/GWAS2/50_irnt.gwas.imputed_v3.both_sexes.tsv.bgz /pl/active/colelab/users/kjames/refinedMR/interim_data/ldsc/NL_bothsexes_temp/
cp /pl/active/colelab/common/published_gwas/Neale_UKB_GWAS_round2/Both_sexes/GWAS2/5084_irnt.gwas.imputed_v3.both_sexes.tsv.bgz /pl/active/colelab/users/kjames/refinedMR/interim_data/ldsc/NL_bothsexes_temp/
cp /pl/active/colelab/common/published_gwas/Neale_UKB_GWAS_round2/Both_sexes/GWAS2/23104_irnt.gwas.imputed_v3.both_sexes.tsv.bgz /pl/active/colelab/users/kjames/refinedMR/interim_data/ldsc/NL_bothsexes_temp/
cp /pl/active/colelab/common/published_gwas/Neale_UKB_GWAS_round2/Both_sexes/GWAS2/23106_irnt.gwas.imputed_v3.both_sexes.tsv.bgz /pl/active/colelab/users/kjames/refinedMR/interim_data/ldsc/NL_bothsexes_temp/
cp /pl/active/colelab/common/published_gwas/Neale_UKB_GWAS_round2/Both_sexes/GWAS2/23101_irnt.gwas.imputed_v3.both_sexes.tsv.bgz /pl/active/colelab/users/kjames/refinedMR/interim_data/ldsc/NL_bothsexes_temp/

# Unzip them
for file in /pl/active/colelab/users/kjames/refinedMR/interim_data/ldsc/NL_bothsexes_temp/*.bgz; do
  zcat "$file" > "${file%.bgz}"
done

# Currently very manual - for 1 trait at a time...
# Clump the gwas results for a given trait
# Define trait
file="/pl/active/colelab/users/kjames/refinedMR/interim_data/ldsc/NL_bothsexes_temp/4124_irnt.gwas.imputed_v3.both_sexes.tsv"
trait="4124_irnt"

# Use plink v1.9 for this; just type plink instead of plink2
for chr in {1..22}; do
  plink \
    --bfile /pl/active/colelab/users/kjames/refinedMR/interim_data/1KGP3_HG19/1KGP3_HG19_files_processed_EUR_Jenkai_chr${chr}_chrbpID \
    --clump ${file} \
    --clump-snp-field variant \
    --clump-field pval \
    --clump-p1 5e-8 \
    --clump-p2 1e-2 \
    --clump-kb 250 \
    --out /pl/active/colelab/users/kjames/refinedMR/interim_data/ldsc/plink_clumped_1KGP_HG19/oilyfish_tg/clumped_results_${trait}_chr${chr}.tsv
done
