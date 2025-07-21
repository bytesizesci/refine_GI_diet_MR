# Goal: Combine PheWAS-Ttest removed variant-trait tables into one
# Author: KJS
# Date: 12/18/24

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
# 1. Load output tables
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
ttest1 <- read.csv("results/PheWAS_Cluster/tables/MR_PheWASCluster_SusSNPs_RelatedTraits_oilyfish_GCST90239664_TG_Graham.csv")
ttest2 <- read.csv("results/PheWAS_Cluster/tables/MR_PheWASCluster_SusSNPs_RelatedTraits_oilyfish_GCST90132314_CVD_Aragam.csv")
ttest3 <- read.csv("results/PheWAS_Cluster/tables/MR_PheWASCluster_SusSNPs_RelatedTraits_oilyfish_GIANT_Height_Yengo.csv")

ttest4 <- read.csv("results/PheWAS_Cluster/tables/MR_PheWASCluster_SusSNPs_RelatedTraits_whiteXwholegrain_GCST90239658_LDL_Graham.csv")
ttest5 <- read.csv("results/PheWAS_Cluster/tables/MR_PheWASCluster_SusSNPs_RelatedTraits_whiteXwholegrain_GCST90132314_CVD_Aragam.csv")
ttest6 <- read.csv("results/PheWAS_Cluster/tables/MR_PheWASCluster_SusSNPs_RelatedTraits_whiteXwholegrain_GIANT_Height_Yengo.csv")

ttest7 <- read.csv("results/PheWAS_Cluster/tables/MR_PheWASCluster_SusSNPs_RelatedTraits_alcohol_GCST90013405_ALT_Pazoki.csv")
ttest8 <- read.csv("results/PheWAS_Cluster/tables/MR_PheWASCluster_SusSNPs_RelatedTraits_alcohol_GCST90319877_cir_Ghouse.csv")
ttest9 <- read.csv("results/PheWAS_Cluster/tables/MR_PheWASCluster_SusSNPs_RelatedTraits_alcohol_GIANT_Height_Yengo.csv")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 2. Combine
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
ttest_all <- rbind(ttest1, ttest2, ttest3, ttest4, ttest5, ttest6, ttest7, ttest8, ttest9)
ttest_all <- ttest_all %>% select(-X)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 3. Save
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
write.csv(ttest_all, "results/general/tables/PheWAS_Ttest_SusSNPS_RelatedTraits_AllPairs.csv")
