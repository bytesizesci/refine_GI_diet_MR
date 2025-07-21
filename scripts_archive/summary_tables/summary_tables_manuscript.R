# Goal: Create summary tables from results for manuscript
# Author: KJS
# Date: November 20, 2024

# IGNORE CODE BELOW AND LOOK AT Nonfiltered_GI/MR_Nonfiltered_GI.R instrument



# Set wd
wdir <- "/pl/active/colelab/users/kjames/refinedMR/"
setwd(wdir)

# Libraries 
library(data.table)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(VennDiagram)

# NON-FILTERED GI

# Load Results UPDATE THESE PATHS
df1 <- read.csv(file = "results/PheWAS_Cluster/tables/MR_PheWASCluster_Filtered_Results_oilyfish_GCST90239664_TG_Graham.csv")
df2 <- read.csv(file = "results/PheWAS_Cluster/tables/MR_PheWASCluster_Filtered_Results_oilyfish_GCST90132314_CVD_Aragam.csv")
df3 <- read.csv(file = "results/PheWAS_Cluster/tables/MR_PheWASCluster_Filtered_Results_oilyfish_GIANT_Height_Yengo.csv")

df4 <- read.csv(file = "results/PheWAS_Cluster/tables/MR_PheWASCluster_Filtered_Results_whiteXwholegrain_GCST90239658_LDL_Graham.csv")
df5 <- read.csv(file = "results/PheWAS_Cluster/tables/MR_PheWASCluster_Filtered_Results_whiteXwholegrain_GCST90132314_CVD_Aragam.csv")
df6 <- read.csv(file = "results/PheWAS_Cluster/tables/MR_PheWASCluster_Filtered_Results_whiteXwholegrain_GIANT_Height_Yengo.csv")

df7 <- read.csv(file = "results/PheWAS_Cluster/tables/MR_PheWASCluster_Filtered_Results_alcohol_GCST90013405_ALT_Pazoki.csv")
df8 <- read.csv(file = "results/PheWAS_Cluster/tables/MR_PheWASCluster_Filtered_Results_alcohol_GCST90319877_cir_Ghouse.csv")
df9 <- read.csv(file = "results/PheWAS_Cluster/tables/MR_PheWASCluster_Filtered_Results_alcohol_GIANT_Height_Yengo.csv")

# Add trait description
df1$exposure <- "Oily Fish"
df1$outcome <- "Triglycerides"

df2$exposure <- "Oily Fish"
df2$outcome <- "Cardiovascular Disease"

df3$exposure <- "Oily Fish"
df4$outcome <- "Height"

df4$exposure <- "White vs Whole Grain Bread"
df4$outcome <- "LDL Cholesterol"

df5$exposure <- "White vs Whole Grain Bread"
df5$outcome <- "Cardiovascular Disease"

df6$exposure <- "White vs Whole Grain Bread"
df6$outcome <- "Height"

df7$exposure <- "Alcohol"
df7$outcome <- "Alanine Aminotransferase"

df8$exposure <- "Alcohol"
df8$outcome <- "Liver Cirrhosis"

df9$exposure <- "Alcohol"
df9$outcome <- "Height"


# Add MR method
df1$mr_method <- "Unfiltered GI"
df2$mr_method <- "Unfiltered GI"
df3$mr_method <- "Unfiltered GI"
df4$mr_method <- "Unfiltered GI"
df5$mr_method <- "Unfiltered GI"
df6$mr_method <- "Unfiltered GI"

# IVW method table first
df1_ivw <- df1[df1$Method == "IVW",]
df2_ivw <- df2[df2$Method == "IVW",]
df3_ivw <- df3[df3$Method == "IVW",]
df4_ivw <- df4[df4$Method == "IVW",]
df5_ivw <- df5[df5$Method == "IVW",]
df6_ivw <- df6[df6$Method == "IVW",]

# Merge
#all_ivw <- rbind(df1_ivw, df2_ivw, df3_ivw, df4_ivw, df5_ivw, df6_ivw)

# Save
#write.csv(all_ivw, "results/julie_collab/summary_tables/allpairs_nonfilteredGI_IVW_table.csv")


# Load Steiger Results
st1 <- read.csv(file = "alpine_from/JulieCollab/Tables/SteigerFiltered_GI/SteigerfilteredGI_MR_oilyfish_Triglycerides.csv")
st2 <- read.csv(file = "alpine_from/JulieCollab/Tables/SteigerFiltered_GI/SteigerfilteredGI_MR_oilyfish_Cardiovascular Disease.csv")
st3 <- read.csv(file = "alpine_from/JulieCollab/Tables/SteigerFiltered_GI/SteigerfilteredGI_MR_whiteXwholegrain_LDL Cholesterol.csv")
st4 <- read.csv(file = "alpine_from/JulieCollab/Tables/SteigerFiltered_GI/SteigerfilteredGI_MR_whiteXwholegrain_Cardiovascular Disease.csv")
st5 <- read.csv(file = "alpine_from/JulieCollab/Tables/SteigerFiltered_GI/SteigerfilteredGI_MR_alcohol_Alanine Aminotransferase.csv")
st6 <- read.csv(file = "alpine_from/JulieCollab/Tables/SteigerFiltered_GI/SteigerfilteredGI_MR_alcohol_Liver Cirrhosis.csv")

# Add trait description
st1$exposure <- "Oily Fish"
st1$outcome <- "Triglycerides"

st2$exposure <- "Oily Fish"
st2$outcome <- "Cardiovascular Disease"

st3$exposure <- "White vs Whole Grain Bread"
st3$outcome <- "LDL Cholesterol"

st4$exposure <- "White vs Whole Grain Bread"
st4$outcome <- "Cardiovascular Disease"

st5$exposure <- "Alcohol"
st5$outcome <- "Alanine Aminotransferase"

st6$exposure <- "Alcohol"
st6$outcome <- "Liver Cirrhosis"

# Add MR method
st1$mr_method <- "Steiger filtered GI"
st2$mr_method <- "Steiger filtered GI"
st3$mr_method <- "Steiger filtered GI"
st4$mr_method <- "Steiger filtered GI"
st5$mr_method <- "Steiger filtered GI"
st6$mr_method <- "Steiger filtered GI"

# IVW method table first
st1_ivw <- st1[st1$Method == "IVW",]
st2_ivw <- st2[st2$Method == "IVW",]
st3_ivw <- st3[st3$Method == "IVW",]
st4_ivw <- st4[st4$Method == "IVW",]
st5_ivw <- st5[st5$Method == "IVW",]
st6_ivw <- st6[st6$Method == "IVW",]

# Load PheWAS Cluster Results
# PheWAS T-test (aka cluster)
# Load Results
ptt1 <- read.csv(file = "results/PheWAS_Cluster/tables/MR_PheWASCluster_Filtered_Results_oilyfish_GCST90239664_TG_Graham.csv")
ptt2 <- read.csv(file = "results/PheWAS_Cluster/tables/MR_PheWASCluster_Filtered_Results_oilyfish_GCST90132314_CVD_Aragam.csv")
ptt3 <- read.csv(file = "results/PheWAS_Cluster/tables/MR_PheWASCluster_Filtered_Results_oilyfish_GIANT_Height_Yengo.csv")

ptt4 <- read.csv(file = "results/PheWAS_Cluster/tables/MR_PheWASCluster_Filtered_Results_whiteXwholegrain_GCST90239658_LDL_Graham.csv")
ptt5 <- read.csv(file = "results/PheWAS_Cluster/tables/MR_PheWASCluster_Filtered_Results_whiteXwholegrain_GCST90132314_CVD_Aragam.csv")
ptt6 <- read.csv(file = "results/PheWAS_Cluster/tables/MR_PheWASCluster_Filtered_Results_whiteXwholegrain_GIANT_Height_Yengo.csv")

ptt7 <- read.csv(file = "results/PheWAS_Cluster/tables/MR_PheWASCluster_Filtered_Results_alcohol_GCST90013405_ALT_Pazoki.csv")
ptt8 <- read.csv(file = "results/PheWAS_Cluster/tables/MR_PheWASCluster_Filtered_Results_alcohol_GCST90319877_cir_Ghouse.csv")
ptt9 <- read.csv(file = "results/PheWAS_Cluster/tables/MR_PheWASCluster_Filtered_Results_alcohol_GIANT_Height_Yengo.csv")

# Add trait description
ptt1$exposure <- "Oily Fish"
ptt1$outcome <- "Triglycerides"

ptt2$exposure <- "Oily Fish"
ptt2$outcome <- "Cardiovascular Disease"

ptt3$exposure <- "Oily Fish"
ptt4$outcome <- "Height"

ptt4$exposure <- "White vs Whole Grain Bread"
ptt4$outcome <- "LDL Cholesterol"

ptt5$exposure <- "White vs Whole Grain Bread"
ptt5$outcome <- "Cardiovascular Disease"

ptt6$exposure <- "White vs Whole Grain Bread"
ptt6$outcome <- "Height"

ptt7$exposure <- "Alcohol"
ptt7$outcome <- "Alanine Aminotransferase"

ptt8$exposure <- "Alcohol"
ptt8$outcome <- "Liver Cirrhosis"

ptt9$exposure <- "Alcohol"
ptt9$outcome <- "Height"

# Add MR method
ptt1$mr_method <- "PheWAS-Cluster filtered GI"
ptt2$mr_method <- "PheWAS-Cluster filtered GI"
ptt3$mr_method <- "PheWAS-Cluster filtered GI"
ptt4$mr_method <- "PheWAS-Cluster filtered GI"
ptt5$mr_method <- "PheWAS-Cluster filtered GI"
ptt6$mr_method <- "PheWAS-Cluster filtered GI"

# IVW method table first
ptt1_ivw <- ptt1[ptt1$Method == "IVW",]
ptt2_ivw <- ptt2[ptt2$Method == "IVW",]
ptt3_ivw <- ptt3[ptt3$Method == "IVW",]
ptt4_ivw <- ptt4[ptt4$Method == "IVW",]
ptt5_ivw <- ptt5[ptt5$Method == "IVW",]
ptt6_ivw <- ptt6[ptt6$Method == "IVW",]

# Load PheWAS Cluster Results
pve1 <- read.csv(file = "alpine_from/JulieCollab/Tables/PheWASPVE/tables/MR_PheWAS_PVE_Filtered_Results_oilyfish_GCST90239664_TG_Graham.csv")
pve2 <- read.csv(file = "alpine_from/JulieCollab/Tables/PheWASPVE/tables/MR_PheWAS_PVE_Filtered_Results_oilyfish_GCST90132314_CVD_Aragam.csv")
pve3 <- read.csv(file = "alpine_from/JulieCollab/Tables/PheWASPVE/tables/MR_PheWAS_PVE_Filtered_Results_whiteXwholegrain_GCST90239658_LDL_Graham.csv")
pve4 <- read.csv(file = "alpine_from/JulieCollab/Tables/PheWASPVE/tables/MR_PheWAS_PVE_Filtered_Results_whiteXwholegrain_GCST90132314_CVD_Aragam.csv")
pve5 <- read.csv(file = "alpine_from/JulieCollab/Tables/PheWASPVE/tables/MR_PheWAS_PVE_Filtered_Results_alcohol_GCST90013405_ALT_Pazoki.csv")
pve6 <- read.csv(file = "alpine_from/JulieCollab/Tables/PheWASPVE/tables/MR_PheWAS_PVE_Filtered_Results_alcohol_GCST90319877_cir_Ghouse.csv")

# Add trait description
pve1$exposure <- "Oily Fish"
pve1$outcome <- "Triglycerides"

pve2$exposure <- "Oily Fish"
pve2$outcome <- "Cardiovascular Disease"

pve3$exposure <- "White vs Whole Grain Bread"
pve3$outcome <- "LDL Cholesterol"

pve4$exposure <- "White vs Whole Grain Bread"
pve4$outcome <- "Cardiovascular Disease"

pve5$exposure <- "Alcohol"
pve5$outcome <- "Alanine Aminotransferase"

pve6$exposure <- "Alcohol"
pve6$outcome <- "Liver Cirrhosis"

# Add MR method
pve1$mr_method <- "PheWAS-PVE filtered GI"
pve2$mr_method <- "PheWAS-PVE filtered GI"
pve3$mr_method <- "PheWAS-PVE filtered GI"
pve4$mr_method <- "PheWAS-PVE filtered GI"
pve5$mr_method <- "PheWAS-PVE filtered GI"
pve6$mr_method <- "PheWAS-PVE filtered GI"

# IVW method table first
pve1_ivw <- pve1[pve1$Method == "IVW",]
pve2_ivw <- pve2[pve2$Method == "IVW",]
pve3_ivw <- pve3[pve3$Method == "IVW",]
pve4_ivw <- pve4[pve4$Method == "IVW",]
pve5_ivw <- pve5[pve5$Method == "IVW",]
pve6_ivw <- pve6[pve6$Method == "IVW",]

# rbind
filtered_all_ivw <- rbind(df1_ivw, st1_ivw, pc1_ivw, pve1_ivw,
                          df2_ivw, st2_ivw,  pc2_ivw, pve2_ivw,
                          df3_ivw, st3_ivw,  pc3_ivw, pve3_ivw,
                          df4_ivw, st4_ivw,  pc4_ivw, pve4_ivw,
                          df5_ivw, st5_ivw,  pc5_ivw, pve5_ivw,
                          df6_ivw, st6_ivw,  pc6_ivw, pve6_ivw)

# Save
write.csv(filtered_all_ivw, "results/tables/allpairs_IVW_table_241028.csv")

# Weigted median table
df1_wm <- df1[df1$Method == "Weighted median",]
df2_wm <- df2[df2$Method == "Weighted median",]
df3_wm <- df3[df3$Method == "Weighted median",]
df4_wm <- df4[df4$Method == "Weighted median",]
df5_wm <- df5[df5$Method == "Weighted median",]
df6_wm <- df6[df6$Method == "Weighted median",]

# Weigted median table
st1_wm <- st1[st1$Method == "Weighted median",]
st2_wm <- st2[st2$Method == "Weighted median",]
st3_wm <- st3[st3$Method == "Weighted median",]
st4_wm <- st4[st4$Method == "Weighted median",]
st5_wm <- st5[st5$Method == "Weighted median",]
st6_wm <- st6[st6$Method == "Weighted median",]

# wm method table first
pc1_wm <- pc1[pc1$Method == "Weighted median",]
pc2_wm <- pc2[pc2$Method == "Weighted median",]
pc3_wm <- pc3[pc3$Method == "Weighted median",]
pc4_wm <- pc4[pc4$Method == "Weighted median",]
pc5_wm <- pc5[pc5$Method == "Weighted median",]
pc6_wm <- pc6[pc6$Method == "Weighted median",]

# wm method table first
pve1_wm <- pve1[pve1$Method == "Weighted median",]
pve2_wm <- pve2[pve2$Method == "Weighted median",]
pve3_wm <- pve3[pve3$Method == "Weighted median",]
pve4_wm <- pve4[pve4$Method == "Weighted median",]
pve5_wm <- pve5[pve5$Method == "Weighted median",]
pve6_wm <- pve6[pve6$Method == "Weighted median",]

# rbind
filtered_all_wm <- rbind(df1_wm, st1_wm, pc1_wm, pve1_wm,
                         df2_wm, st2_wm,  pc2_wm, pve2_wm,
                         df3_wm, st3_wm,  pc3_wm, pve3_wm,
                         df4_wm, st4_wm,  pc4_wm, pve4_wm,
                         df5_wm, st5_wm,  pc5_wm, pve5_wm,
                         df6_wm, st6_wm,  pc6_wm, pve6_wm)


# Save
write.csv(filtered_all_wm, "results/tables/allpairs_WM_table_241028.csv")
#
#
#
# Now, get F values, Q-stat, and P-value of Q-stat
# Load unfiltered
df1_het <- read.csv("alpine_from/JulieCollab/Tables/Nonfiltered_GI/NonfilteredGI_MR_IVWOily fish intake_Triglycerides.csv")
df2_het <- read.csv("alpine_from/JulieCollab/Tables/Nonfiltered_GI/NonfilteredGI_MR_IVWOily fish intake_CVD.csv")
df3_het <- read.csv("alpine_from/JulieCollab/Tables/Nonfiltered_GI/NonfilteredGI_MR_IVWWhite vs Whole Grain Bread_LDL.csv")
df4_het <- read.csv("alpine_from/JulieCollab/Tables/Nonfiltered_GI/NonfilteredGI_MR_IVWWhite vs Whole Grain Bread_CVD.csv")
df5_het <- read.csv("alpine_from/JulieCollab/Tables/Nonfiltered_GI/NonfilteredGI_MR_IVWAlcohol_Alanine Aminotransferase.csv")
df6_het <- read.csv("alpine_from/JulieCollab/Tables/Nonfiltered_GI/NonfilteredGI_MR_IVWAlcohol_Cirrhosis.csv")
# Add variable names
df1_het$pair <- "Oily Fish on Tryglicerides"
df2_het$pair <- "Oily Fish on Cardiovascular Disease"
df3_het$pair <- "White vs Whole Grain Bread on LDL Cholesterol"
df4_het$pair <- "White vs Whole Grain Bread on Cardiovascular Disease"
df5_het$pair <- "Alcohol on Alanine Aminotransferase"
df6_het$pair <- "Alcohol on Liver Cirrhosis"

# Load Steiger
st1_het <- read.csv("alpine_from/JulieCollab/Tables/SteigerFiltered_GI/MR_Steiger_Filtered_IVWres_oilyfish_GCST90239664_TG_Graham.csv")
st2_het <- read.csv("alpine_from/JulieCollab/Tables/SteigerFiltered_GI/MR_Steiger_Filtered_IVWres_oilyfish_GCST90132314_CVD_Aragam.csv")
st3_het <- read.csv("alpine_from/JulieCollab/Tables/SteigerFiltered_GI/MR_Steiger_Filtered_IVWres_whiteXwholegrain_GCST90239658_LDL_Graham.csv")
st4_het <- read.csv("alpine_from/JulieCollab/Tables/SteigerFiltered_GI/MR_Steiger_Filtered_IVWres_whiteXwholegrain_GCST90132314_CVD_Aragam.csv")
st5_het <- read.csv("alpine_from/JulieCollab/Tables/SteigerFiltered_GI/MR_Steiger_Filtered_IVWres_alcohol_GCST90013405_ALT_Pazoki.csv")
st6_het <- read.csv("alpine_from/JulieCollab/Tables/SteigerFiltered_GI/MR_Steiger_Filtered_IVWres_alcohol_GCST90319877_cir_Ghouse.csv")
# Add variable names
st1_het$pair <- "Oily Fish on Tryglicerides"
st2_het$pair <- "Oily Fish on Cardiovascular Disease"
st3_het$pair <- "White vs Whole Grain Bread on LDL Cholesterol"
st4_het$pair <- "White vs Whole Grain Bread on Cardiovascular Disease"
st5_het$pair <- "Alcohol on Alanine Aminotransferase"
st6_het$pair <- "Alcohol on Liver Cirrhosis"

# Load PheWASCluster
pc1_het <- read.csv("alpine_from/JulieCollab/Tables/PheWAS_Cluster/MR_PheWASCluster_Filtered_IVWres_oilyfish_GCST90239664_TG_Graham.csv")
pc2_het <- read.csv("alpine_from/JulieCollab/Tables/PheWAS_Cluster/MR_PheWASCluster_Filtered_IVWres_oilyfish_GCST90132314_CVD_Aragam.csv")
pc3_het <- read.csv("alpine_from/JulieCollab/Tables/PheWAS_Cluster/MR_PheWASCluster_Filtered_IVWres_whiteXwholegrain_GCST90239658_LDL_Graham.csv")
pc4_het <- read.csv("alpine_from/JulieCollab/Tables/PheWAS_Cluster/MR_PheWASCluster_Filtered_IVWres_whiteXwholegrain_GCST90132314_CVD_Aragam.csv")
pc5_het <- read.csv("alpine_from/JulieCollab/Tables/PheWAS_Cluster/MR_PheWASCluster_Filtered_IVWres_alcohol_GCST90013405_ALT_Pazoki.csv")
pc6_het <- read.csv("alpine_from/JulieCollab/Tables/PheWAS_Cluster/MR_PheWASCluster_Filtered_IVWres_alcohol_GCST90319877_cir_Ghouse.csv")
# Add variable names
pc1_het$pair <- "Oily Fish on Tryglicerides"
pc2_het$pair <- "Oily Fish on Cardiovascular Disease"
pc3_het$pair <- "White vs Whole Grain Bread on LDL Cholesterol"
pc4_het$pair <- "White vs Whole Grain Bread on Cardiovascular Disease"
pc5_het$pair <- "Alcohol on Alanine Aminotransferase"
pc6_het$pair <- "Alcohol on Liver Cirrhosis"

# Load PheWAS PVE 
pve1_het <- read.csv("alpine_from/JulieCollab/Tables/PheWASPVE/tables/MR_PheWAS_PVE_Filtered_IVWres_oilyfish_GCST90239664_TG_Graham.csv")
pve2_het <- read.csv("alpine_from/JulieCollab/Tables/PheWASPVE/tables/MR_PheWAS_PVE_Filtered_IVWres_oilyfish_GCST90132314_CVD_Aragam.csv")
pve3_het <- read.csv("alpine_from/JulieCollab/Tables/PheWASPVE/tables/MR_PheWAS_PVE_Filtered_IVWres_whiteXwholegrain_GCST90239658_LDL_Graham.csv")
pve4_het <- read.csv("alpine_from/JulieCollab/Tables/PheWASPVE/tables/MR_PheWAS_PVE_Filtered_IVWres_whiteXwholegrain_GCST90132314_CVD_Aragam.csv")
pve5_het <- read.csv("alpine_from/JulieCollab/Tables/PheWASPVE/tables/MR_PheWAS_PVE_Filtered_IVWres_alcohol_GCST90013405_ALT_Pazoki.csv")
pve6_het <- read.csv("alpine_from/JulieCollab/Tables/PheWASPVE/tables/MR_PheWAS_PVE_Filtered_IVWres_alcohol_GCST90319877_cir_Ghouse.csv")
# Add variable names
pve1_het$pair <- "Oily Fish on Tryglicerides"
pve2_het$pair <- "Oily Fish on Cardiovascular Disease"
pve3_het$pair <- "White vs Whole Grain Bread on LDL Cholesterol"
pve4_het$pair <- "White vs Whole Grain Bread on Cardiovascular Disease"
pve5_het$pair <- "Alcohol on Alanine Aminotransferase"
pve6_het$pair <- "Alcohol on Liver Cirrhosis"

# rbind
hetstat_all_ivw <- rbind(df1_het, st1_het, pc1_het, pve1_het,
                         df2_het, st2_het, pc2_het, pve2_het,
                         df3_het, st3_het, pc3_het, pve3_het,
                         df4_het, st4_het, pc4_het, pve4_het,
                         df5_het, st5_het, pc5_het, pve5_het,
                         df6_het, st6_het, pc6_het, pve6_het)

# Save
write.csv(hetstat_all_ivw, "results/tables/allpairs_IVW_fstat_het_table_241029.csv")

