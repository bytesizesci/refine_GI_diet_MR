# Goal: Manuscript vignettes
# Author: KJS
# Date: 12/05/24

# Libraries 
library(data.table)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(VennDiagram)

# Set wd
wdir <- "/pl/active/colelab/users/kjames/refinedMR/"
setwd(wdir)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 1. Format and save tables with instrument information
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
mrobj <- readRDS("interim_data/merged_GI/MRobject_steiger_PheWAStt_PheWASPVE_sensitivity_241203.rds")

# Add Steiger keep column
# Using `lapply()` to add a new stgfilt_keep column
mrobj <- lapply(mrobj, function(df) {
  # If the direction is as hyposthesized (TRUE) and the evidence is strong (P<0.05) keep SNP
  df$Steiger_keep <- ifelse(df$steiger_dir == "TRUE" & df$steiger_pval <= 0.05, TRUE, FALSE )
  return(df)
})

cols2save <- c("exposure","outcome", "SNP","chr.exposure","pos.exposure","effect_allele.exposure","other_allele.exposure","oaf.exposure","se.exposure",
               "pval.exposure", "chr.outcome", "pos.outcome","effect_allele.outcome","other_allele.outcome","eaf.outcome",
               "beta.outcome","se.outcome","pval.outcome", "variant", "proxy","original_exposure_SNP", "PheWAS_ttest_keep",
               "PheWAS_ttest_sensitivity_keep", "PheWAS_PVE_keep", "PheWAS_PVE_sensitivity_keep", "Steiger_keep") 

# Using `lapply()` to add a new stgfilt_keep column
mrobj_sub <- lapply(mrobj, function(df) {
  # If the direction is as hyposthesized (TRUE) and the evidence is strong (P<0.05) keep SNP
  df <- df[,colnames(df) %in% cols2save]
  return(df)
})

# Merge such that if a col is missing in a list slot it gets NA for those rows
gi_table <- bind_rows(mrobj_sub)

# Rearrange
gi_table <- gi_table %>% select(exposure, outcome, variant, SNP, proxy, original_exposure_SNP, everything()) %>%
  rename("rsID" = "SNP")
table(gi_table$proxy)

# Due to varrying formats of GWAS outcome data, the proxy/original SNP columns were not ubiquitously needed
# Thus, they are blank for some trait pairs.
# To make the table complete, if you are NA in proxy, that means you are original
gi_table$proxy [is.na(gi_table$proxy)] <- "original_SNP"
table(gi_table$proxy)

# Save
write.csv(gi_table, "results/general/tables/GI_variant_supplemental_table.csv")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 2. What SNP is removed from bread --> LDL with Steiger filtering
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
bread_ldl <- mrobj$whiteXwholegrain_LDL
weirdSNP <- bread_ldl[bread_ldl$Steiger_keep == FALSE,]
weirdSNP

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 3. Format PheWAS-PVE table for supplement
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
pve1 <- read.csv("results/PheWAS_PVE/tables/MR_PheWAS_PVE_RemovedVariants_top1_oilyfish_GCST90239664_TG_Graham.csv")
pve1$exp_out <- "Oily fish on Tryglycerides"
pve2 <- read.csv("results/PheWAS_PVE/tables/MR_PheWAS_PVE_RemovedVariants_top1_oilyfish_GCST90132314_CVD_Aragam.csv")
pve2$exp_out <- "Oily fish on Cardiovascular Disease"
pve3 <- read.csv("results/PheWAS_PVE/tables/MR_PheWAS_PVE_RemovedVariants_top1_oilyfish_GIANT_Height_Yengo.csv")
pve3$exp_out <- "Oily fish on Height"
pve4 <- read.csv("results/PheWAS_PVE/tables/MR_PheWAS_PVE_RemovedVariants_top1_whiteXwholegrain_GCST90239658_LDL_Graham.csv")
pve4$exp_out <- "White vs Whole Grain or Brown Bread on LDL Cholesterol"
pve5 <- read.csv("results/PheWAS_PVE/tables/MR_PheWAS_PVE_RemovedVariants_top1_whiteXwholegrain_GCST90132314_CVD_Aragam.csv")
pve5$exp_out <-  "White vs Whole Grain or Brown Bread on Cardiovascular Disease"
pve6 <- read.csv("results/PheWAS_PVE/tables/MR_PheWAS_PVE_RemovedVariants_top1_whiteXwholegrain_height.csv")
pve6$exp_out <-  "White vs Whole Grain or Brown Bread on Height"
pve7 <- read.csv("results/PheWAS_PVE/tables/MR_PheWAS_PVE_RemovedVariants_top1_alcohol_GCST90013405_ALT_Pazoki.csv")
pve7$exp_out <- "Alcohol on Alanine Aminotransferase"
pve8 <- read.csv("results/PheWAS_PVE/tables/MR_PheWAS_PVE_RemovedVariants_top1_alcohol_GCST90319877_cir_Ghouse.csv")
pve8$exp_out <- "Alcohol on Liver Cirrhosis"
pve9 <- read.csv("results/PheWAS_PVE/tables/MR_PheWAS_PVE_RemovedVariants_top1_alcohol_height.csv")
pve9$exp_out <- "Alcohol on Height"

# Keep select cols
pve_keep_cols <- c("exp_out", "trait", "description", "variant", "PVE")

# Subset
pve1_s <- pve1[,colnames(pve1) %in% pve_keep_cols]
pve2_s <- pve2[,colnames(pve2) %in% pve_keep_cols]
pve3_s <- pve3[,colnames(pve3) %in% pve_keep_cols]
pve4_s <- pve4[,colnames(pve4) %in% pve_keep_cols]
pve5_s <- pve5[,colnames(pve5) %in% pve_keep_cols]
pve6_s <- pve6[,colnames(pve6) %in% pve_keep_cols]
pve7_s <- pve7[,colnames(pve7) %in% pve_keep_cols]
pve8_s <- pve8[,colnames(pve8) %in% pve_keep_cols]
pve9_s <- pve9[,colnames(pve9) %in% pve_keep_cols]

# Merge 
pve_all <- rbind(pve1_s,pve2_s,pve3_s,pve4_s,pve5_s,pve6_s,pve7_s,pve8_s,pve9_s)

# Wrangle
pve_all <- pve_all %>% select(exp_out, trait, description, variant, PVE) %>% 
  rename("PheWAS_Trait_ID" = "trait",
         "Description"="description",
         "Exposure-Outcome" = "exp_out",
         "Variant" = "variant") 

# Save
write.csv(pve_all, "results/general/tables/PVE_removedtraits_all.csv")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 4. Repeat:Sensitivity - Format PheWAS-PVE sensitivity table for supplement
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
pve1 <- read.csv("results/PheWAS_PVE_Sensitivity/tables/MR_PheWAS_PVE_RemovedVariants_top1_oilyfish_GCST90239664_TG_Graham_Sensitivity.csv")
pve1$exp_out <- "Oily fish on Tryglycerides"
pve2 <- read.csv("results/PheWAS_PVE_Sensitivity/tables/MR_PheWAS_PVE_RemovedVariants_top1_oilyfish_GCST90132314_CVD_Aragam_Sensitivity.csv")
pve2$exp_out <- "Oily fish on Cardiovascular Disease"
pve3 <- read.csv("results/PheWAS_PVE_Sensitivity/tables/MR_PheWAS_PVE_RemovedVariants_top1_oilyfish_GIANT_Height_Yengo_Sensitivity.csv")
pve3$exp_out <- "Oily fish on Height"
pve4 <- read.csv("results/PheWAS_PVE_Sensitivity/tables/MR_PheWAS_PVE_RemovedVariants_top1_whiteXwholegrain_GCST90239658_LDL_Graham_Sensitivity.csv")
pve4$exp_out <- "White vs Whole Grain or Brown Bread on LDL Cholesterol"
pve5 <- read.csv("results/PheWAS_PVE_Sensitivity/tables/MR_PheWAS_PVE_RemovedVariants_top1_whiteXwholegrain_GCST90132314_CVD_Aragam_Sensitivity.csv")
pve5$exp_out <-  "White vs Whole Grain or Brown Bread on Cardiovascular Disease"
pve6 <- read.csv("results/PheWAS_PVE_Sensitivity/tables/MR_PheWAS_PVE_RemovedVariants_top1_whiteXwholegrain_height_Sensitivity.csv")
pve6$exp_out <-  "White vs Whole Grain or Brown Bread on Height"
pve7 <- read.csv("results/PheWAS_PVE_Sensitivity/tables/MR_PheWAS_PVE_RemovedVariants_top1_alcohol_GCST90013405_ALT_Pazoki_Sensitivity.csv")
pve7$exp_out <- "Alcohol on Alanine Aminotransferase"
pve8 <- read.csv("results/PheWAS_PVE_Sensitivity/tables/MR_PheWAS_PVE_RemovedVariants_top1_alcohol_GCST90319877_cir_Ghouse_Sensitivity.csv")
pve8$exp_out <- "Alcohol on Liver Cirrhosis"
pve9 <- read.csv("results/PheWAS_PVE_Sensitivity/tables/MR_PheWAS_PVE_RemovedVariants_top1_alcohol_height_Sensitivity.csv")
pve9$exp_out <- "Alcohol on Height"

# Keep select cols
# Keep select cols
pve_keep_cols <- c("exp_out", "trait", "description", "variant", "PVE")

# Subset
pve1_s <- pve1[,colnames(pve1) %in% pve_keep_cols]
pve2_s <- pve2[,colnames(pve2) %in% pve_keep_cols]
pve3_s <- pve3[,colnames(pve3) %in% pve_keep_cols]
pve4_s <- pve4[,colnames(pve4) %in% pve_keep_cols]
pve5_s <- pve5[,colnames(pve5) %in% pve_keep_cols]
pve6_s <- pve6[,colnames(pve6) %in% pve_keep_cols]
pve7_s <- pve7[,colnames(pve7) %in% pve_keep_cols]
pve8_s <- pve8[,colnames(pve8) %in% pve_keep_cols]
pve9_s <- pve9[,colnames(pve9) %in% pve_keep_cols]

# Merge 
pve_all <- rbind(pve1_s,pve2_s,pve3_s,pve4_s,pve5_s,pve6_s,pve7_s,pve8_s,pve9_s)

# Wrangle
# Merge 
pve_all <- rbind(pve1_s,pve2_s,pve3_s,pve4_s,pve5_s,pve6_s,pve7_s,pve8_s,pve9_s)

# Wrangle
pve_all <- pve_all %>% select(exp_out, trait, description, variant, PVE) %>% 
  rename("PheWAS_Trait_ID" = "trait",
         "Description"="description",
         "Exposure-Outcome" = "exp_out",
         "Variant" = "variant") 

# Save
write.csv(pve_all, "results/general/tables/PVE_removedtraits_all_sensitivity.csv")
