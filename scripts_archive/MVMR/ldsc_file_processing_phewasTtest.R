# Date: 01/30/2025
# Goal: Develop trait selection pipeline; get lists of files in descending order to compare using ldsc
# Using PheWAS-Ttest results
# Author: KJS

#Set working direcotry
setwd("/pl/active/colelab/users/kjames/refinedMR/")

library(dplyr)

# Load for loop
df_names <- data.frame(EXP_pheno = c("oilyfish", "oilyfish", "oilyfish", "whiteXwholegrain", "whiteXwholegrain",  "whiteXwholegrain", "alcohol", "alcohol","alcohol"),
                       OUT_pheno = c("GCST90239664_TG_Graham", "GCST90132314_CVD_Aragam", "GCST006901_Height_Yengo_2018", "GCST90239658_LDL_Graham","GCST90132314_CVD_Aragam", "GCST006901_Height_Yengo_2018", "GCST90013405_ALT_Pazoki", "GCST90319877_cir_Ghouse", "GCST006901_Height_Yengo_2018"),
                       EXP_names = c("Oily Fish", "Oily fish", "Oily fish", "White vs Whole Grain Bread", "White vs Whole Grain Bread", "White vs Whole Grain Bread","Alcohol", "Alcohol", "Alcohol"),
                       OUT_names = c("Triglycerides", "Cardiovascular Disease", "Height", "LDL Cholesterol", "Cardiovascular Disease", "Height", "Alanine Aminotransferase", "Liver Cirrhosis", "Height"))


# Going to get snp and file names from PheWAS-TTest workflow
file_paths <- data.frame(pve_file_paths = c("/pl/active/colelab/users/kjames/refinedMR/results/PheWAS_Cluster/tables/MR_PheWASCluster_SusSNPs_RelatedTraits_FilePaths_oilyfish_GCST90239664_TG_Graham.csv",
                                            "/pl/active/colelab/users/kjames/refinedMR/results/PheWAS_Cluster/tables/MR_PheWASCluster_SusSNPs_RelatedTraits_FilePaths_oilyfish_GCST90132314_CVD_Aragam.csv",
                                            "/pl/active/colelab/users/kjames/refinedMR/results/PheWAS_Cluster/tables/MR_PheWASCluster_SusSNPs_RelatedTraits_FilePaths_oilyfish_GIANT_Height_Yengo.csv",
                                            "/pl/active/colelab/users/kjames/refinedMR/results/PheWAS_Cluster/tables/MR_PheWASCluster_SusSNPs_RelatedTraits_FilePaths_whiteXwholegrain_GCST90239658_LDL_Graham.csv",
                                            "/pl/active/colelab/users/kjames/refinedMR/results/PheWAS_Cluster/tables/MR_PheWASCluster_SusSNPs_RelatedTraits_FilePaths_whiteXwholegrain_GCST90132314_CVD_Aragam.csv",
                                            "/pl/active/colelab/users/kjames/refinedMR/results/PheWAS_Cluster/tables/MR_PheWASCluster_SusSNPs_RelatedTraits_FilePaths_whiteXwholegrain_GIANT_Height_Yengo.csv",
                                            "/pl/active/colelab/users/kjames/refinedMR/results/PheWAS_Cluster/tables/MR_PheWASCluster_SusSNPs_RelatedTraits_FilePaths_alcohol_GCST90013405_ALT_Pazoki.csv",
                                            "/pl/active/colelab/users/kjames/refinedMR/results/PheWAS_Cluster/tables/MR_PheWASCluster_SusSNPs_RelatedTraits_FilePaths_alcohol_GCST90319877_cir_Ghouse.csv",
                                            "/pl/active/colelab/users/kjames/refinedMR/results/PheWAS_Cluster/tables/MR_PheWASCluster_SusSNPs_RelatedTraits_FilePaths_alcohol_GIANT_Height_Yengo.csv"))


# Testing, comment out when not testing
#df <- read.csv("/pl/active/colelab/users/kjames/refinedMR/results/PheWAS_Cluster/tables/MR_PheWASCluster_SusSNPs_RelatedTraits_FilePaths_oilyfish_GCST90239664_TG_Graham.csv")
#i=1

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
  #pve2 <- merge(pve_temp, pwas_nr, by = "trait", all.x = TRUE, all.y = FALSE)
  pve_temp <- pve_temp %>% select(rmd_SNPnames, rmd_traitnames)
  
  # Arrange rows by alphabetical order 
  pve_sorted <- pve_temp %>%
    arrange(rmd_traitnames)
  
  # Just want file name in descending order
  filenames <- pve_sorted$rmd_traitnames
  
  # Gsub gwas with ldsc
  filenames <- gsub("gwas", "ldsc", filenames)
  
  # Append "/pl/active/colelab/users/kjames/refinedMR/interim_data/ldsc/PVE_descending_input" to the front of each string
  path_to_add <- "/pl/active/colelab/common/published_gwas/Neale_UKB_GWAS_round2/Both_sexes/LDSC/"
  filenames_descending <- paste0(path_to_add, filenames, ".bgz")
  
  # Save
  write.table(filenames_descending, paste0("interim_data/ldsc/PheWAS_Ttest_LDSC_input/", EXP_pheno, "_", OUT_pheno, "_descending_TTest_filenames.txt"), quote = FALSE, row.names = FALSE, col.names = FALSE)
}


# Note, you can combine all of the file paths into one file to submit one master gcor run for this project...
# Load results
df1 <- read.csv("interim_data/ldsc/PheWAS_Ttest_LDSC_input/oilyfish_GCST90239664_TG_Graham_descending_TTest_filenames.txt", header = FALSE)
df2 <- read.csv("interim_data/ldsc/PheWAS_Ttest_LDSC_input/oilyfish_GCST90132314_CVD_Aragam_descending_TTest_filenames.txt", header = FALSE)
df3 <- read.csv("interim_data/ldsc/PheWAS_Ttest_LDSC_input/oilyfish_GCST006901_Height_Yengo_2018_descending_TTest_filenames.txt", header = FALSE)
df4 <- read.csv("interim_data/ldsc/PheWAS_Ttest_LDSC_input/whiteXwholegrain_GCST90239658_LDL_Graham_descending_TTest_filenames.txt", header = FALSE)
df5 <- read.csv("interim_data/ldsc/PheWAS_Ttest_LDSC_input/whiteXwholegrain_GCST90132314_CVD_Aragam_descending_TTest_filenames.txt", header = FALSE)
df6 <- read.csv("interim_data/ldsc/PheWAS_Ttest_LDSC_input/whiteXwholegrain_GCST006901_Height_Yengo_2018_descending_TTest_filenames.txt", header = FALSE)
df7 <- read.csv("interim_data/ldsc/PheWAS_Ttest_LDSC_input/alcohol_GCST90013405_ALT_Pazoki_descending_TTest_filenames.txt", header = FALSE)
df8 <- read.csv("interim_data/ldsc/PheWAS_Ttest_LDSC_input/alcohol_GCST90319877_cir_Ghouse_descending_TTest_filenames.txt", header = FALSE)
df9 <- read.csv("interim_data/ldsc/PheWAS_Ttest_LDSC_input/alcohol_GCST006901_Height_Yengo_2018_descending_TTest_filenames.txt", header = FALSE)

#rowbind
dfall <- rbind(df1, df2, df3, df4, df5, df6, df7, df8, df9)

# Unique
dfall_uniq <- unique(dfall$V1)
length(dfall_uniq) #122

# Save
write.table(dfall_uniq, "interim_data/ldsc/PheWAS_Ttest_LDSC_input/PheWAS_Ttest_ALLFILES_for_gcor.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)
