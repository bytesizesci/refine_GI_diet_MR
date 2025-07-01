# Date: 11/19/24 
# Goal: Develop trait selection pipeline; get lists of files in descending order to compare using ldsc
# Author: KJS

#Set working direcotry
setwd("/pl/active/colelab/users/kjames/refinedMR/")

library(dplyr)

# Load for loop
df_names <- data.frame(EXP_pheno = c("oilyfish", "oilyfish", "oilyfish", "whiteXwholegrain", "whiteXwholegrain",  "whiteXwholegrain", "alcohol", "alcohol","alcohol"),
                       OUT_pheno = c("GCST90239664_TG_Graham", "GCST90132314_CVD_Aragam", "GCST006901_Height_Yengo_2018", "GCST90239658_LDL_Graham","GCST90132314_CVD_Aragam", "GCST006901_Height_Yengo_2018", "GCST90013405_ALT_Pazoki", "GCST90319877_cir_Ghouse", "GCST006901_Height_Yengo_2018"),
                       EXP_names = c("Oily Fish", "Oily fish", "Oily fish", "White vs Whole Grain Bread", "White vs Whole Grain Bread", "White vs Whole Grain Bread","Alcohol", "Alcohol", "Alcohol"),
                       OUT_names = c("Triglycerides", "Cardiovascular Disease", "Height", "LDL Cholesterol", "Cardiovascular Disease", "Height", "Alanine Aminotransferase", "Liver Cirrhosis", "Height"))

# PVE file paths for loop
file_paths <- data.frame(pve_file_paths = c("/pl/active/colelab/users/kjames/refinedMR/results/PheWAS_PVE/tables/MR_PheWAS_PVE_RemovedVariants_top1_oilyfish_GCST90239664_TG_Graham.csv",
                                            "/pl/active/colelab/users/kjames/refinedMR/results/PheWAS_PVE/tables/MR_PheWAS_PVE_RemovedVariants_top1_oilyfish_GCST90132314_CVD_Aragam.csv",
                                            "/pl/active/colelab/users/kjames/refinedMR/results/PheWAS_PVE/tables/MR_PheWAS_PVE_RemovedVariants_top1_oilyfish_GIANT_Height_Yengo.csv",
                                            "/pl/active/colelab/users/kjames/refinedMR/results/PheWAS_PVE/tables/MR_PheWAS_PVE_RemovedVariants_top1_whiteXwholegrain_GCST90239658_LDL_Graham.csv",
                                            "/pl/active/colelab/users/kjames/refinedMR/results/PheWAS_PVE/tables/MR_PheWAS_PVE_RemovedVariants_top1_whiteXwholegrain_GCST90132314_CVD_Aragam.csv",
                                            "/pl/active/colelab/users/kjames/refinedMR/results/PheWAS_PVE/tables/MR_PheWAS_PVE_RemovedVariants_top1_whiteXwholegrain_height.csv",
                                            "/pl/active/colelab/users/kjames/refinedMR/results/PheWAS_PVE/tables/MR_PheWAS_PVE_RemovedVariants_top1_alcohol_GCST90013405_ALT_Pazoki.csv",
                                            "/pl/active/colelab/users/kjames/refinedMR/results/PheWAS_PVE/tables/MR_PheWAS_PVE_RemovedVariants_top1_alcohol_GCST90319877_cir_Ghouse.csv",
                                            "/pl/active/colelab/users/kjames/refinedMR/results/PheWAS_PVE/tables/MR_PheWAS_PVE_RemovedVariants_top1_alcohol_height.csv"))

# Need file paths to link traits to ldsc files so load PheWAS results
pwas <- read.csv("interim_data/PheWAS_Grep_SNPs/PheWAS_Oilyfish_WhiteXWhole_Alcochol_UpdatedProxies.csv")
# Make file name prettier. Get rid of extra junk.
pwas$trait <- sub("\\..*$", "", pwas$file_name)
# Subset to cols of interest
pwas <- pwas %>% select(file_name, trait)
# Remove duplicates, keeping the first occurrence
pwas_nr <- pwas[!duplicated(pwas$trait), ]

# Enter loop
for(i in seq_len(nrow(file_paths))){
  # Load a PheWAS-PVE file
  pve_temp <- read.csv(file_paths[i,])
  # Get names
  EXP_pheno <- df_names[i,1]
  OUT_pheno <- df_names[i,2]
  EXP_name <- df_names[i,3]
  OUT_name <- df_names[i,4]
  
  # merge to get file ID
  pve2 <- merge(pve_temp, pwas_nr, by = "trait", all.x = TRUE, all.y = FALSE)
  
  # Arrange rows by descending order of pve
  pve2_sorted <- pve2 %>%
    arrange(desc(PVE))
  
  # Just want file name in descending order
  filenames_descending <- pve2_sorted$file_name.x
  
  # Gsub gwas with ldsc
  filenames_descending <- gsub("gwas", "ldsc", filenames_descending)
  
  # Append "/pl/active/colelab/users/kjames/refinedMR/interim_data/ldsc/PVE_descending_input" to the front of each string
  path_to_add <- "/pl/active/colelab/common/published_gwas/Neale_UKB_GWAS_round2/Both_sexes/LDSC/"
  filenames_descending <- paste0(path_to_add, filenames_descending, ".bgz")
  
  # Save
  write.table(filenames_descending, paste0("interim_data/ldsc/PVE_descending_input/", EXP_pheno, "_", OUT_pheno, "_descendingPVE_filenames.txt"), quote = FALSE, row.names = FALSE, col.names = FALSE)
}

