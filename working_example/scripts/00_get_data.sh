# Goal: Get data for working_example. 

# To submit,
# Make the script executable with: chmod +x 00_get_data.sh
# Run the script: sbatch 01_ldsc_rg.sh
# NOTE: For working_example, submit this file from working_example/scripts

# Create a config folder in working_example
mkdir ../config 

# 30080_irnt
wget https://broad-ukb-sumstats-us-east-1.s3.amazonaws.com/round2/additive-tsvs/30080_irnt.gwas.imputed_v3.both_sexes.tsv.bgz -O 30080_irnt.gwas.imputed_v3.both_sexes.tsv.bgz

# 30270_irnt
wget https://broad-ukb-sumstats-us-east-1.s3.amazonaws.com/round2/additive-tsvs/30270_irnt.gwas.imputed_v3.both_sexes.tsv.bgz -O 30270_irnt.gwas.imputed_v3.both_sexes.tsv.bgz

# 23098_irnt
wget https://broad-ukb-sumstats-us-east-1.s3.amazonaws.com/round2/additive-tsvs/23098_irnt.gwas.imputed_v3.both_sexes.tsv.bgz -O 23098_irnt.gwas.imputed_v3.both_sexes.tsv.bgz

# 20015_irnt
wget https://broad-ukb-sumstats-us-east-1.s3.amazonaws.com/round2/additive-tsvs/20015_irnt.gwas.imputed_v3.both_sexes.tsv.bgz -O 20015_irnt.gwas.imputed_v3.both_sexes.tsv.bgz

# 23106_irnt
wget https://broad-ukb-sumstats-us-east-1.s3.amazonaws.com/round2/additive-tsvs/23106_irnt.gwas.imputed_v3.both_sexes.tsv.bgz -O 23106_irnt.gwas.imputed_v3.both_sexes.tsv.bgz

# ALT 
wget https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90013001-GCST90014000/GCST90013405/harmonised/33972514-GCST90013405-EFO_0004735-Build37.f.tsv.gz

# Alcohol, https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90132001-GCST90133000/GCST90132992/
wget https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90132001-GCST90133000/GCST90132992/BOLTlmm_UKB_genoQCEURN455146_v3_diet_anyalcohol_glassespermonth.derived.average_QT_BgenSnps_mac20_maf0.005_info0.6.gz

# Unzip the files
zcat 30080_irnt.gwas.imputed_v3.both_sexes.tsv.bgz 30080_irnt.gwas.imputed_v3.both_sexes.tsv
zcat 30270_irnt.gwas.imputed_v3.both_sexes.tsv.bgz 30270_irnt.gwas.imputed_v3.both_sexes.tsv
zcat 23098_irnt.gwas.imputed_v3.both_sexes.tsv.bgz 23098_irnt.gwas.imputed_v3.both_sexes.tsv
zcat 20015_irnt.gwas.imputed_v3.both_sexes.tsv.bgz 20015_irnt.gwas.imputed_v3.both_sexes.tsv
zcat 23106_irnt.gwas.imputed_v3.both_sexes.tsv.bgz 23106_irnt.gwas.imputed_v3.both_sexes.tsv
zcat 33972514-GCST90013405-EFO_0004735-Build37.f.tsv.gz 33972514-GCST90013405-EFO_0004735-Build37.f.tsv
zcat BOLTlmm_UKB_genoQCEURN455146_v3_diet_anyalcohol_glassespermonth.derived.average_QT_BgenSnps_mac20_maf0.005_info0.6.gz BOLTlmm_UKB_genoQCEURN455146_v3_diet_anyalcohol_glassespermonth.derived.average_QT_BgenSnps_mac20_maf0.005_info0.6


