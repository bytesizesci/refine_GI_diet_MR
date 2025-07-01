# Date: 05/05/2025 
# Goal: Find SNP overlap between the phenos in the clusters and diet traits.
# This builds from the "gcor_kmeans_alcohol_ALT_v2.R" script, which builds clusters.
#
# Author: KJS

# For each cluster, we need to find the trait that has the most overlap with diet
# To do this, 
# - get list of traits that were assigned to clusters
# - get the list of SNP-trait connections from PheWAS-TTest method
# - count how many SNPs assigned to each trait
# - For each cluster, see which trait had the most SNPs associated to it
# - Boom! Those are the traits to move to the MVMR

#Set working direcotry
setwd("/pl/active/colelab/users/kjames/refinedMR/")

# Libraries
library(dplyr)
library(ggplot2)
library(tidyr)
library(RColorBrewer)

# Define exposure and outcome of analysis
EXP_pheno <- "alcohol"
OUT_pheno <- "ALT"

# Load the cluster results 
clst <- read.csv("interim_data/MVMR/alcohol_ALT/hclust_complete_cut5.csv")
clst$trait_id <- gsub("ldsc", "gwas", clst$trait_id)
clst$trait_id <- gsub(".bgz", "", clst$trait_id)

# Load list of SNP-trait connections
sus <- read.csv("results/PheWAS_Cluster/tables/MR_PheWASCluster_SusSNPs_RelatedTraits_FilePaths_alcohol_GCST90013405_ALT_Pazoki.csv")

# combine clst and sus
df <- merge(clst, sus, by.x = "trait_id", by.y = "rmd_traitnames")

# Count how many SNPs relate to each trait
df_counts <- df %>%
  distinct(trait_id, rmd_SNPnames) %>%
  count(trait_id, name = "num_snps")

# Merge with clst
df_counts2 <- merge(df_counts, clst, by.x = "trait_id", by.y = "trait_id")

# Find trait with most SNPs per cluster
top_traits <- df_counts2 %>%
  group_by(cluster) %>%
  slice_max(num_snps, with_ties = TRUE) %>%
  ungroup()

# Merge traits with their description
key <- read.csv("/pl/active/colelab/common/published_gwas/Neale_UKB_GWAS_round2/metadata/UKBB GWAS Imputed v3 - File Manifest Release 20180731 - Description Lookup.csv",header = FALSE)
key <- unique(key)
# Isolate the ID portion from file name
top_traits$trait_id_num <- sub("\\..*","", top_traits$trait_id)
top_traits2 <- merge(top_traits, key, by.x = "trait_id_num", by.y = "V1")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# We still have ties in the top_traits2 object. Test whether heritability can distinguish which traits to pick from the ties
# File with gcors from traits identified in PheWAS-Ttest method
gcor <- read.table("/pl/active/colelab/users/kjames/refinedMR/interim_data/ldsc/processed_results/gcor_PheWASttest_results_all.txt")

# Header
gcor_header <- c("p1", "p2", "rg", "se", "z", "p", "h2_obs", "h2_obs_se", "h2_int", "h2_int_se", "gcov_int", "gcov_int_se")

# Assign header to df
colnames(gcor) <- gcor_header

# Subset
gcor <- gcor %>% select(p2, h2_obs, h2_obs_se)

# Get trait ID
gcor$prefix <- sub("\\..*", "", gcor$p2)

# Get unique
gcor <- unique(gcor)

# Manually add heritability for 1448 trait (1st trait that gets left out in heritability table)
# See file: ldsc_results_starting_with_0.txt.log
# Heritability of phenotype 1
# ---------------------------
# Total Observed scale h2: 0.0047 (0.0015)
manual_row <- data.frame("p2" = "1448_4.ldsc.imputed_v3.both_sexes.tsv.bgz",
                         "h2_obs" = 0.0047,
                         "h2_obs_se"= 0.0015,
                         "prefix"= "1448_4")
gcor <- rbind(manual_row, gcor)

# Merge with dat
top_traits3 <- merge(top_traits2, gcor, by.x = "trait_id_num", by.y="prefix")

# Find the trait with the highest h2_obs within each Cluster
top_traits4 <- top_traits3 %>%
  group_by(cluster) %>%
  slice_max(h2_obs, with_ties = TRUE)

print(top_traits4)

# Save dat3
write.csv(top_traits4, paste0("interim_data/MVMR/",EXP_pheno,"_",OUT_pheno,"/hclust_maxh2obs_cut5.csv"))
