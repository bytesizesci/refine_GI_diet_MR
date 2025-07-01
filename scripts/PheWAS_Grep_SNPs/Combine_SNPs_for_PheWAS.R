# Date: 11/12/24 
# Goal: Get new list of SNPs to grep for PheWAS (with exposure proxies)
# Author: KJS

# Set wd
wdir <- "/pl/active/colelab/users/kjames/refinedMR/"
setwd(wdir)

# Load mr data
df1 <- readRDS("interim_data/merged_GI/mrdat_exp_out_proxies_oilyfish_GCST90239664_TG_Graham_GRCh37.rds")
df2 <- readRDS("interim_data/merged_GI/mrdat_exp_out_proxies_oilyfish_GCST90132314_CVD_Aragam_GRCh37.rds")
df3 <- readRDS("interim_data/merged_GI/mrdat_exp_out_proxies_whiteXwholegrain_GCST90239658_LDL_Graham.rds")
df4 <- readRDS("interim_data/merged_GI/mrdat_exp_out_proxies_whiteXwholegrain_GCST90132314_CVD_Aragam.rds")
df5 <- readRDS("interim_data/merged_GI/mrdat_exp_out_proxies_alcohol_GCST90013405_ALT_Pazoki.rds")
df6 <- readRDS("interim_data/merged_GI/mrdat_exp_out_proxies_alcohol_GCST90319877_cir_Ghouse.rds")

# Get exposure snps
exp1 <- df1$variant
exp2 <- df2$variant
exp3 <- df3$variant
exp4 <- df4$variant
exp5 <- df5$variant
exp6 <- df6$variant

# Combine
exp_snps <- c(exp1, exp2, exp3, exp4, exp5, exp6)

# Save
write.csv(exp_snps, "interim_data/merged_GI/PheWAS_variants_oilyfish_whiteXwholegrain_alcohol_updated_exposure_proxies.csv")
