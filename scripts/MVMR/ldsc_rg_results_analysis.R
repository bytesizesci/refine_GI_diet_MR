# MVMR Analysis: Process ldsc rg results
# Author: KJS
# Date: 11/26/24 

# Goal: Perform percent variance explained thresholding on genetic instrument 

# Set wd
wdir <- "/pl/active/colelab/users/kjames/refinedMR/"
setwd(wdir)

# Load libraries
library(data.table)
library(dplyr)
library(ggplot2)
library(ggrepel)

# Get ldsc rg results
gc <- read.table("interim_data/ldsc/processed_results/all_results.txt", sep = "", header = FALSE)

# Define column names
gc_cols <- c("p1", "p2", "rg", "se", "z", "p", "h2_obs", "h2_obs_se", "h2_int", "h2_int_se", "gcov_int", "gcov_int_se")
colnames(gc) <- gc_cols

# Get trait id from filename
gc$p1_name <- gsub("\\..*", "", gc$p1)
gc$p2_name <- gsub("\\..*", "", gc$p2)

# Get the manifest key for PheWAS
key <- read.csv("/pl/active/colelab/common/published_gwas/Neale_UKB_GWAS_round2/metadata/UKBB GWAS Imputed v3 - File Manifest Release 20180731 - Description Lookup.csv",header = FALSE)
key <- key[!duplicated(key),] #keeps the first occurrence of each duplicate

# 1st Merge
gc_named <- merge(gc, key, by.x = "p1_name", by.y = "V1",all.x = TRUE, all.y = FALSE)
# After 1st merge, rename "V2" col to "p1_description"
gc_named <- gc_named %>% rename("p1_desc" = "V2")

# 2nd Merge
gc_named <- merge(gc_named, key, by.x = "p2_name", by.y = "V1", all.x = TRUE, all.y = FALSE)
# After 2nd merge, rename "V2" col to "p2_description"
gc_named <- gc_named %>% rename("p2_desc" = "V2")

# Save table
write.csv(gc_named, "interim_data/ldsc/processed_results/LDSC_rg_whiteXwholegrain_GCST90239658_LDL_Graham.csv")
