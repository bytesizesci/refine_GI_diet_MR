# Date: 11/19/24 
# Goal: Get new list of SNPs to grep for PheWAS with height (with exposure proxies)
# Author: KJS

# Set wd
wdir <- "/pl/active/colelab/users/kjames/refinedMR/"
setwd(wdir)

# Load mr data
df1 <- readRDS("interim_data/merged_GI/mrdat_exp_out_proxies_oilyfish_GCST006901_Height_Yengo_2018.rds")
df2 <- readRDS("interim_data/merged_GI/mrdat_exp_out_proxies_whiteXwholegrain_GCST006901_Height_Yengo_2018.rds")
df3 <- readRDS("interim_data/merged_GI/mrdat_exp_out_proxies_alcohol_GCST006901_Height_Yengo_2018.rds")


# Get exposure snps
exp1 <- df1$variant
exp2 <- df2$variant
exp3 <- df3$variant

# Combine
exp_snps <- c(exp1, exp2, exp3)

# Save
write.csv(exp_snps, "interim_data/merged_GI/PheWAS_variants_height_outcome_proxies.csv")
