# Goal: Update MR Object with inclusion/exclusion columns based on filtering
# Author: KJS
# Date: December 3, 2024

# Set wd
wdir <- "/pl/active/colelab/users/kjames/refinedMR/"
setwd(wdir)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Load MR Object with Steiger results
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~
mrobj <- readRDS("interim_data/merged_GI/MRobject_steiger_241113.rds")
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# PheWAS T-test results 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~
# Oilyfish - TG
#~~~~~~~~~~~
# USER INPUT - select correct data; used for obtaining filtering labels
load("interim_data/PheWAS_Cluster/PheWASCluster_objects_oilyfish_GCST90239664_TG_Graham.Rdata")

# USER INPUT - select correct slot
mrobj_slot <- mrobj$oilyfish_TG

# Make column with sus_snps
mrobj_slot$chr_pos <- regmatches(mrobj_slot$variant, regexpr("[a-z0-9]+\\:[a-z0-9]+", mrobj_slot$variant, ignore.case = TRUE))
sus_SNPs_chr_pos <- regmatches(sus_SNPs, regexpr("[a-z0-9]+\\:[a-z0-9]+", sus_SNPs, ignore.case = TRUE))
mrobj_slot$PheWAS_ttest_keep <- ifelse(mrobj_slot$chr_pos %in% sus_SNPs_chr_pos, FALSE, TRUE) # keep the "TRUE"
table(mrobj_slot$PheWAS_ttest_keep)

# USER INPUT - select correct slot
mrobj$oilyfish_TG <- mrobj_slot

# Clean up - get rid of everything other than mrobj
rm(list = c("df", "exp_b", "exp_n", "exp_pval", "exp_se", "exp_stdB", "exp_stdSE",
            "mrobj_slot", "sig_ind", "stdBeta_df", "stdBeta_df_noEXP", 
            "stdBeta_EXP", "sus_traits_df", "tstat_df", "tstat_dif",
            "tstat_pval_df", "sus_SNPs", "sus_SNPs_chr_pos", "sus_traits"))


#~~~~~~~~~~~
# Oilyfish - CVD
#~~~~~~~~~~~
# USER INPUT - select correct data; used for obtaining filtering labels
load("interim_data/PheWAS_Cluster/PheWASCluster_objects_oilyfish_GCST90132314_CVD_Aragam.Rdata")

# USER INPUT - select correct slot
mrobj_slot <- mrobj$oilyfish_CVD

# Make column with sus_snps
mrobj_slot$chr_pos <- regmatches(mrobj_slot$variant, regexpr("[a-z0-9]+\\:[a-z0-9]+", mrobj_slot$variant, ignore.case = TRUE))
sus_SNPs_chr_pos <- regmatches(sus_SNPs, regexpr("[a-z0-9]+\\:[a-z0-9]+", sus_SNPs, ignore.case = TRUE))
mrobj_slot$PheWAS_ttest_keep <- ifelse(mrobj_slot$chr_pos %in% sus_SNPs_chr_pos, FALSE, TRUE) # keep the "TRUE"
table(mrobj_slot$PheWAS_ttest_keep)

# USER INPUT - select correct slot
mrobj$oilyfish_CVD <- mrobj_slot

# Clean up - get rid of everything other than mrobj
rm(list = c("df", "exp_b", "exp_n", "exp_pval", "exp_se", "exp_stdB", "exp_stdSE",
            "mrobj_slot", "sig_ind", "stdBeta_df", "stdBeta_df_noEXP", 
            "stdBeta_EXP", "sus_traits_df", "tstat_df", "tstat_dif",
            "tstat_pval_df", "sus_SNPs", "sus_SNPs_chr_pos", "sus_traits"))


#~~~~~~~~~~~
# Oilyfish - Height
#~~~~~~~~~~~
# USER INPUT - select correct data; used for obtaining filtering labels
load("interim_data/PheWAS_Cluster/PheWASCluster_objects_oilyfish_GIANT_Height_Yengo.Rdata")

# USER INPUT - select correct slot
mrobj_slot <- mrobj$oilyfish_height

# Make column with sus_snps
mrobj_slot$chr_pos <- regmatches(mrobj_slot$variant, regexpr("[a-z0-9]+\\:[a-z0-9]+", mrobj_slot$variant, ignore.case = TRUE))
sus_SNPs_chr_pos <- regmatches(sus_SNPs, regexpr("[a-z0-9]+\\:[a-z0-9]+", sus_SNPs, ignore.case = TRUE))
mrobj_slot$PheWAS_ttest_keep <- ifelse(mrobj_slot$chr_pos %in% sus_SNPs_chr_pos, FALSE, TRUE) # keep the "TRUE"
table(mrobj_slot$PheWAS_ttest_keep)

# USER INPUT - select correct slot
mrobj$oilyfish_height <- mrobj_slot

# Clean up - get rid of everything other than mrobj
rm(list = c("df", "exp_b", "exp_n", "exp_pval", "exp_se", "exp_stdB", "exp_stdSE",
            "mrobj_slot", "sig_ind", "stdBeta_df", "stdBeta_df_noEXP", 
            "stdBeta_EXP", "sus_traits_df", "tstat_df", "tstat_dif",
            "tstat_pval_df", "sus_SNPs", "sus_SNPs_chr_pos", "sus_traits"))

#~~~~~~~~~~~
# White X whole grain/brown bread - LDL
#~~~~~~~~~~~
# USER INPUT - select correct data; used for obtaining filtering labels
load("interim_data/PheWAS_Cluster/PheWASCluster_objects_whiteXwholegrain_GCST90239658_LDL_Graham.Rdata")

# USER INPUT - select correct slot
mrobj_slot <- mrobj$whiteXwholegrain_LDL

# Make column with sus_snps
mrobj_slot$chr_pos <- regmatches(mrobj_slot$variant, regexpr("[a-z0-9]+\\:[a-z0-9]+", mrobj_slot$variant, ignore.case = TRUE))
sus_SNPs_chr_pos <- regmatches(sus_SNPs, regexpr("[a-z0-9]+\\:[a-z0-9]+", sus_SNPs, ignore.case = TRUE))
mrobj_slot$PheWAS_ttest_keep <- ifelse(mrobj_slot$chr_pos %in% sus_SNPs_chr_pos, FALSE, TRUE) # keep the "TRUE"
table(mrobj_slot$PheWAS_ttest_keep)

# USER INPUT - select correct slot
mrobj$whiteXwholegrain_LDL <- mrobj_slot

# Clean up - get rid of everything other than mrobj
rm(list = c("df", "exp_b", "exp_n", "exp_pval", "exp_se", "exp_stdB", "exp_stdSE",
            "mrobj_slot", "sig_ind", "stdBeta_df", "stdBeta_df_noEXP", 
            "stdBeta_EXP", "sus_traits_df", "tstat_df", "tstat_dif",
            "tstat_pval_df", "sus_SNPs", "sus_SNPs_chr_pos", "sus_traits"))


#~~~~~~~~~~~
# White X whole grain/brown bread - CVD
#~~~~~~~~~~~
# USER INPUT - select correct data; used for obtaining filtering labels
load("interim_data/PheWAS_Cluster/PheWASCluster_objects_whiteXwholegrain_GCST90132314_CVD_Aragam.Rdata")

# USER INPUT - select correct slot
mrobj_slot <- mrobj$whiteXwholegrain_CVD

# Make column with sus_snps
mrobj_slot$chr_pos <- regmatches(mrobj_slot$variant, regexpr("[a-z0-9]+\\:[a-z0-9]+", mrobj_slot$variant, ignore.case = TRUE))
sus_SNPs_chr_pos <- regmatches(sus_SNPs, regexpr("[a-z0-9]+\\:[a-z0-9]+", sus_SNPs, ignore.case = TRUE))
mrobj_slot$PheWAS_ttest_keep <- ifelse(mrobj_slot$chr_pos %in% sus_SNPs_chr_pos, FALSE, TRUE) # keep the "TRUE"
table(mrobj_slot$PheWAS_ttest_keep)

# USER INPUT - select correct slot
mrobj$whiteXwholegrain_CVD <- mrobj_slot

# Clean up - get rid of everything other than mrobj
rm(list = c("df", "exp_b", "exp_n", "exp_pval", "exp_se", "exp_stdB", "exp_stdSE",
            "mrobj_slot", "sig_ind", "stdBeta_df", "stdBeta_df_noEXP", 
            "stdBeta_EXP", "sus_traits_df", "tstat_df", "tstat_dif",
            "tstat_pval_df", "sus_SNPs", "sus_SNPs_chr_pos", "sus_traits"))

#~~~~~~~~~~~
# White X whole grain/brown bread - Height
#~~~~~~~~~~~
# USER INPUT - select correct data; used for obtaining filtering labels
load("interim_data/PheWAS_Cluster/PheWASCluster_objects_whiteXwholegrain_GIANT_Height_Yengo.Rdata")

# USER INPUT - select correct slot
mrobj_slot <- mrobj$whiteXwholegrain_height

# Make column with sus_snps
mrobj_slot$chr_pos <- regmatches(mrobj_slot$variant, regexpr("[a-z0-9]+\\:[a-z0-9]+", mrobj_slot$variant, ignore.case = TRUE))
sus_SNPs_chr_pos <- regmatches(sus_SNPs, regexpr("[a-z0-9]+\\:[a-z0-9]+", sus_SNPs, ignore.case = TRUE))
mrobj_slot$PheWAS_ttest_keep <- ifelse(mrobj_slot$chr_pos %in% sus_SNPs_chr_pos, FALSE, TRUE) # keep the "TRUE"
table(mrobj_slot$PheWAS_ttest_keep)

# USER INPUT - select correct slot
mrobj$whiteXwholegrain_height <- mrobj_slot

# Clean up - get rid of everything other than mrobj
rm(list = c("df", "exp_b", "exp_n", "exp_pval", "exp_se", "exp_stdB", "exp_stdSE",
            "mrobj_slot", "sig_ind", "stdBeta_df", "stdBeta_df_noEXP", 
            "stdBeta_EXP", "sus_traits_df", "tstat_df", "tstat_dif",
            "tstat_pval_df", "sus_SNPs", "sus_SNPs_chr_pos", "sus_traits"))


#~~~~~~~~~~~
# Alcohol - ALT
#~~~~~~~~~~~
# USER INPUT - select correct data; used for obtaining filtering labels
load("interim_data/PheWAS_Cluster/PheWASCluster_objects_alcohol_GCST90013405_ALT_Pazoki.Rdata")

# USER INPUT - select correct slot
mrobj_slot <- mrobj$alcohol_ALT

# Make column with sus_snps
mrobj_slot$chr_pos <- regmatches(mrobj_slot$variant, regexpr("[a-z0-9]+\\:[a-z0-9]+", mrobj_slot$variant, ignore.case = TRUE))
sus_SNPs_chr_pos <- regmatches(sus_SNPs, regexpr("[a-z0-9]+\\:[a-z0-9]+", sus_SNPs, ignore.case = TRUE))
mrobj_slot$PheWAS_ttest_keep <- ifelse(mrobj_slot$chr_pos %in% sus_SNPs_chr_pos, FALSE, TRUE) # keep the "TRUE"
table(mrobj_slot$PheWAS_ttest_keep)

# USER INPUT - select correct slot
mrobj$alcohol_ALT <- mrobj_slot

# Clean up - get rid of everything other than mrobj
rm(list = c("df", "exp_b", "exp_n", "exp_pval", "exp_se", "exp_stdB", "exp_stdSE",
            "mrobj_slot", "sig_ind", "stdBeta_df", "stdBeta_df_noEXP", 
            "stdBeta_EXP", "sus_traits_df", "tstat_df", "tstat_dif",
            "tstat_pval_df", "sus_SNPs", "sus_SNPs_chr_pos", "sus_traits"))


#~~~~~~~~~~~
# Alcohol - Liver cirrhosis
#~~~~~~~~~~~
# USER INPUT - select correct data; used for obtaining filtering labels
load("interim_data/PheWAS_Cluster/PheWASCluster_objects_alcohol_GCST90319877_cir_Ghouse.Rdata")

# USER INPUT - select correct slot
mrobj_slot <- mrobj$alcohol_cir

# Make column with sus_snps
mrobj_slot$chr_pos <- regmatches(mrobj_slot$variant, regexpr("[a-z0-9]+\\:[a-z0-9]+", mrobj_slot$variant, ignore.case = TRUE))
sus_SNPs_chr_pos <- regmatches(sus_SNPs, regexpr("[a-z0-9]+\\:[a-z0-9]+", sus_SNPs, ignore.case = TRUE))
mrobj_slot$PheWAS_ttest_keep <- ifelse(mrobj_slot$chr_pos %in% sus_SNPs_chr_pos, FALSE, TRUE) # keep the "TRUE"
table(mrobj_slot$PheWAS_ttest_keep)

# USER INPUT - select correct slot
mrobj$alcohol_cir <- mrobj_slot

# Clean up - get rid of everything other than mrobj
rm(list = c("df", "exp_b", "exp_n", "exp_pval", "exp_se", "exp_stdB", "exp_stdSE",
            "mrobj_slot", "sig_ind", "stdBeta_df", "stdBeta_df_noEXP", 
            "stdBeta_EXP", "sus_traits_df", "tstat_df", "tstat_dif",
            "tstat_pval_df", "sus_SNPs", "sus_SNPs_chr_pos", "sus_traits"))


#~~~~~~~~~~~
# Alcohol - Height
#~~~~~~~~~~~
# USER INPUT - select correct data; used for obtaining filtering labels
load("interim_data/PheWAS_Cluster/PheWASCluster_objects_alcohol_GIANT_Height_Yengo.Rdata")

# USER INPUT - select correct slot
mrobj_slot <- mrobj$alcohol_height

# Make column with sus_snps
mrobj_slot$chr_pos <- regmatches(mrobj_slot$variant, regexpr("[a-z0-9]+\\:[a-z0-9]+", mrobj_slot$variant, ignore.case = TRUE))
sus_SNPs_chr_pos <- regmatches(sus_SNPs, regexpr("[a-z0-9]+\\:[a-z0-9]+", sus_SNPs, ignore.case = TRUE))
mrobj_slot$PheWAS_ttest_keep <- ifelse(mrobj_slot$chr_pos %in% sus_SNPs_chr_pos, FALSE, TRUE) # keep the "TRUE"
table(mrobj_slot$PheWAS_ttest_keep)

# USER INPUT - select correct slot
mrobj$alcohol_height <- mrobj_slot

# Clean up - get rid of everything other than mrobj
rm(list = c("df", "exp_b", "exp_n", "exp_pval", "exp_se", "exp_stdB", "exp_stdSE",
            "mrobj_slot", "sig_ind", "stdBeta_df", "stdBeta_df_noEXP", 
            "stdBeta_EXP", "sus_traits_df", "tstat_df", "tstat_dif",
            "tstat_pval_df", "sus_SNPs", "sus_SNPs_chr_pos", "sus_traits"))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# PheWAS T-test results - SENSITIVITY ANALYSIS
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~
# Oilyfish - TG
#~~~~~~~~~~~
# USER INPUT - select correct data; used for obtaining filtering labels
load("interim_data/PheWAS_Ttest_Sensitivity/PheWASCluster_objects_oilyfish_GCST90239664_TG_Graham_Sensitivity.Rdata")

# USER INPUT - select correct slot
mrobj_slot <- mrobj$oilyfish_TG

# Make column with sus_snps
mrobj_slot$chr_pos <- regmatches(mrobj_slot$variant, regexpr("[a-z0-9]+\\:[a-z0-9]+", mrobj_slot$variant, ignore.case = TRUE))
sus_SNPs_chr_pos <- regmatches(sus_SNPs, regexpr("[a-z0-9]+\\:[a-z0-9]+", sus_SNPs, ignore.case = TRUE))
mrobj_slot$PheWAS_ttest_sensitivity_keep <- ifelse(mrobj_slot$chr_pos %in% sus_SNPs_chr_pos, FALSE, TRUE) # keep the "TRUE"
table(mrobj_slot$PheWAS_ttest_sensitivity_keep)

# USER INPUT - select correct slot
mrobj$oilyfish_TG <- mrobj_slot

# Clean up - get rid of everything other than mrobj
rm(list = c("df", "exp_b", "exp_n", "exp_pval", "exp_se", "exp_stdB", "exp_stdSE",
            "mrobj_slot", "sig_ind", "stdBeta_df", "stdBeta_df_noEXP", 
            "stdBeta_EXP", "sus_traits_df", "tstat_df", "tstat_dif",
            "tstat_pval_df", "sus_SNPs", "sus_SNPs_chr_pos", "sus_traits"))


#~~~~~~~~~~~
# Oilyfish - CVD
#~~~~~~~~~~~
# USER INPUT - select correct data; used for obtaining filtering labels
load("interim_data/PheWAS_Ttest_Sensitivity/PheWASCluster_objects_oilyfish_GCST90132314_CVD_Aragam_Sensitivity.Rdata")

# USER INPUT - select correct slot
mrobj_slot <- mrobj$oilyfish_CVD

# Make column with sus_snps
mrobj_slot$chr_pos <- regmatches(mrobj_slot$variant, regexpr("[a-z0-9]+\\:[a-z0-9]+", mrobj_slot$variant, ignore.case = TRUE))
sus_SNPs_chr_pos <- regmatches(sus_SNPs, regexpr("[a-z0-9]+\\:[a-z0-9]+", sus_SNPs, ignore.case = TRUE))
mrobj_slot$PheWAS_ttest_sensitivity_keep <- ifelse(mrobj_slot$chr_pos %in% sus_SNPs_chr_pos, FALSE, TRUE) # keep the "TRUE"
table(mrobj_slot$PheWAS_ttest_sensitivity_keep)

# USER INPUT - select correct slot
mrobj$oilyfish_CVD <- mrobj_slot

# Clean up - get rid of everything other than mrobj
rm(list = c("df", "exp_b", "exp_n", "exp_pval", "exp_se", "exp_stdB", "exp_stdSE",
            "mrobj_slot", "sig_ind", "stdBeta_df", "stdBeta_df_noEXP", 
            "stdBeta_EXP", "sus_traits_df", "tstat_df", "tstat_dif",
            "tstat_pval_df", "sus_SNPs", "sus_SNPs_chr_pos", "sus_traits"))


#~~~~~~~~~~~
# Oilyfish - Height
#~~~~~~~~~~~
# USER INPUT - select correct data; used for obtaining filtering labels
load("interim_data/PheWAS_Ttest_Sensitivity/PheWASCluster_objects_oilyfish_GIANT_Height_Yengo_Sensitivity.Rdata")

# USER INPUT - select correct slot
mrobj_slot <- mrobj$oilyfish_height

# Make column with sus_snps
mrobj_slot$chr_pos <- regmatches(mrobj_slot$variant, regexpr("[a-z0-9]+\\:[a-z0-9]+", mrobj_slot$variant, ignore.case = TRUE))
sus_SNPs_chr_pos <- regmatches(sus_SNPs, regexpr("[a-z0-9]+\\:[a-z0-9]+", sus_SNPs, ignore.case = TRUE))
mrobj_slot$PheWAS_ttest_sensitivity_keep <- ifelse(mrobj_slot$chr_pos %in% sus_SNPs_chr_pos, FALSE, TRUE) # keep the "TRUE"
table(mrobj_slot$PheWAS_ttest_sensitivity_keep)

# USER INPUT - select correct slot
mrobj$oilyfish_height <- mrobj_slot

# Clean up - get rid of everything other than mrobj
rm(list = c("df", "exp_b", "exp_n", "exp_pval", "exp_se", "exp_stdB", "exp_stdSE",
            "mrobj_slot", "sig_ind", "stdBeta_df", "stdBeta_df_noEXP", 
            "stdBeta_EXP", "sus_traits_df", "tstat_df", "tstat_dif",
            "tstat_pval_df", "sus_SNPs", "sus_SNPs_chr_pos", "sus_traits"))

#~~~~~~~~~~~
# White X whole grain/brown bread - LDL
#~~~~~~~~~~~
# USER INPUT - select correct data; used for obtaining filtering labels
load("interim_data/PheWAS_Ttest_Sensitivity/PheWASCluster_objects_whiteXwholegrain_GCST90239658_LDL_Graham_Sensitivity.Rdata")

# USER INPUT - select correct slot
mrobj_slot <- mrobj$whiteXwholegrain_LDL

# Make column with sus_snps
mrobj_slot$chr_pos <- regmatches(mrobj_slot$variant, regexpr("[a-z0-9]+\\:[a-z0-9]+", mrobj_slot$variant, ignore.case = TRUE))
sus_SNPs_chr_pos <- regmatches(sus_SNPs, regexpr("[a-z0-9]+\\:[a-z0-9]+", sus_SNPs, ignore.case = TRUE))
mrobj_slot$PheWAS_ttest_sensitivity_keep <- ifelse(mrobj_slot$chr_pos %in% sus_SNPs_chr_pos, FALSE, TRUE) # keep the "TRUE"
table(mrobj_slot$PheWAS_ttest_sensitivity_keep)

# USER INPUT - select correct slot
mrobj$whiteXwholegrain_LDL <- mrobj_slot

# Clean up - get rid of everything other than mrobj
rm(list = c("df", "exp_b", "exp_n", "exp_pval", "exp_se", "exp_stdB", "exp_stdSE",
            "mrobj_slot", "sig_ind", "stdBeta_df", "stdBeta_df_noEXP", 
            "stdBeta_EXP", "sus_traits_df", "tstat_df", "tstat_dif",
            "tstat_pval_df", "sus_SNPs", "sus_SNPs_chr_pos", "sus_traits"))


#~~~~~~~~~~~
# White X whole grain/brown bread - CVD
#~~~~~~~~~~~
# USER INPUT - select correct data; used for obtaining filtering labels
load("interim_data/PheWAS_Ttest_Sensitivity/PheWASCluster_objects_whiteXwholegrain_GCST90132314_CVD_Aragam_Sensitivity.Rdata")

# USER INPUT - select correct slot
mrobj_slot <- mrobj$whiteXwholegrain_CVD

# Make column with sus_snps
mrobj_slot$chr_pos <- regmatches(mrobj_slot$variant, regexpr("[a-z0-9]+\\:[a-z0-9]+", mrobj_slot$variant, ignore.case = TRUE))
sus_SNPs_chr_pos <- regmatches(sus_SNPs, regexpr("[a-z0-9]+\\:[a-z0-9]+", sus_SNPs, ignore.case = TRUE))
mrobj_slot$PheWAS_ttest_sensitivity_keep <- ifelse(mrobj_slot$chr_pos %in% sus_SNPs_chr_pos, FALSE, TRUE) # keep the "TRUE"
table(mrobj_slot$PheWAS_ttest_sensitivity_keep)

# USER INPUT - select correct slot
mrobj$whiteXwholegrain_CVD <- mrobj_slot

# Clean up - get rid of everything other than mrobj
rm(list = c("df", "exp_b", "exp_n", "exp_pval", "exp_se", "exp_stdB", "exp_stdSE",
            "mrobj_slot", "sig_ind", "stdBeta_df", "stdBeta_df_noEXP", 
            "stdBeta_EXP", "sus_traits_df", "tstat_df", "tstat_dif",
            "tstat_pval_df", "sus_SNPs", "sus_SNPs_chr_pos", "sus_traits"))

#~~~~~~~~~~~
# White X whole grain/brown bread - Height
#~~~~~~~~~~~
# USER INPUT - select correct data; used for obtaining filtering labels
load("interim_data/PheWAS_Ttest_Sensitivity/PheWASCluster_objects_whiteXwholegrain_GIANT_Height_Yengo_Sensitivity.Rdata")

# USER INPUT - select correct slot
mrobj_slot <- mrobj$whiteXwholegrain_height

# Make column with sus_snps
mrobj_slot$chr_pos <- regmatches(mrobj_slot$variant, regexpr("[a-z0-9]+\\:[a-z0-9]+", mrobj_slot$variant, ignore.case = TRUE))
sus_SNPs_chr_pos <- regmatches(sus_SNPs, regexpr("[a-z0-9]+\\:[a-z0-9]+", sus_SNPs, ignore.case = TRUE))
mrobj_slot$PheWAS_ttest_sensitivity_keep <- ifelse(mrobj_slot$chr_pos %in% sus_SNPs_chr_pos, FALSE, TRUE) # keep the "TRUE"
table(mrobj_slot$PheWAS_ttest_sensitivity_keep)

# USER INPUT - select correct slot
mrobj$whiteXwholegrain_height <- mrobj_slot

# Clean up - get rid of everything other than mrobj
rm(list = c("df", "exp_b", "exp_n", "exp_pval", "exp_se", "exp_stdB", "exp_stdSE",
            "mrobj_slot", "sig_ind", "stdBeta_df", "stdBeta_df_noEXP", 
            "stdBeta_EXP", "sus_traits_df", "tstat_df", "tstat_dif",
            "tstat_pval_df", "sus_SNPs", "sus_SNPs_chr_pos", "sus_traits"))


#~~~~~~~~~~~
# Alcohol - ALT
#~~~~~~~~~~~
# USER INPUT - select correct data; used for obtaining filtering labels
load("interim_data/PheWAS_Ttest_Sensitivity/PheWASCluster_objects_alcohol_GCST90013405_ALT_Pazoki_Sensitivity.Rdata")

# USER INPUT - select correct slot
mrobj_slot <- mrobj$alcohol_ALT

# Make column with sus_snps
mrobj_slot$chr_pos <- regmatches(mrobj_slot$variant, regexpr("[a-z0-9]+\\:[a-z0-9]+", mrobj_slot$variant, ignore.case = TRUE))
sus_SNPs_chr_pos <- regmatches(sus_SNPs, regexpr("[a-z0-9]+\\:[a-z0-9]+", sus_SNPs, ignore.case = TRUE))
mrobj_slot$PheWAS_ttest_sensitivity_keep <- ifelse(mrobj_slot$chr_pos %in% sus_SNPs_chr_pos, FALSE, TRUE) # keep the "TRUE"
table(mrobj_slot$PheWAS_ttest_sensitivity_keep)

# USER INPUT - select correct slot
mrobj$alcohol_ALT <- mrobj_slot

# Clean up - get rid of everything other than mrobj
rm(list = c("df", "exp_b", "exp_n", "exp_pval", "exp_se", "exp_stdB", "exp_stdSE",
            "mrobj_slot", "sig_ind", "stdBeta_df", "stdBeta_df_noEXP", 
            "stdBeta_EXP", "sus_traits_df", "tstat_df", "tstat_dif",
            "tstat_pval_df", "sus_SNPs", "sus_SNPs_chr_pos", "sus_traits"))


#~~~~~~~~~~~
# Alcohol - Liver cirrhosis
#~~~~~~~~~~~
# USER INPUT - select correct data; used for obtaining filtering labels
load("interim_data/PheWAS_Ttest_Sensitivity/PheWASCluster_objects_alcohol_GCST90319877_cir_Ghouse_Sensitivity.Rdata")

# USER INPUT - select correct slot
mrobj_slot <- mrobj$alcohol_cir

# Make column with sus_snps
mrobj_slot$chr_pos <- regmatches(mrobj_slot$variant, regexpr("[a-z0-9]+\\:[a-z0-9]+", mrobj_slot$variant, ignore.case = TRUE))
sus_SNPs_chr_pos <- regmatches(sus_SNPs, regexpr("[a-z0-9]+\\:[a-z0-9]+", sus_SNPs, ignore.case = TRUE))
mrobj_slot$PheWAS_ttest_sensitivity_keep <- ifelse(mrobj_slot$chr_pos %in% sus_SNPs_chr_pos, FALSE, TRUE) # keep the "TRUE"
table(mrobj_slot$PheWAS_ttest_sensitivity_keep)

# USER INPUT - select correct slot
mrobj$alcohol_cir <- mrobj_slot

# Clean up - get rid of everything other than mrobj
rm(list = c("df", "exp_b", "exp_n", "exp_pval", "exp_se", "exp_stdB", "exp_stdSE",
            "mrobj_slot", "sig_ind", "stdBeta_df", "stdBeta_df_noEXP", 
            "stdBeta_EXP", "sus_traits_df", "tstat_df", "tstat_dif",
            "tstat_pval_df", "sus_SNPs", "sus_SNPs_chr_pos", "sus_traits"))


#~~~~~~~~~~~
# Alcohol - Height
#~~~~~~~~~~~
# USER INPUT - select correct data; used for obtaining filtering labels
load("interim_data/PheWAS_Ttest_Sensitivity/PheWASCluster_objects_alcohol_GIANT_Height_Yengo_Sensitivity.Rdata")

# USER INPUT - select correct slot
mrobj_slot <- mrobj$alcohol_height

# Make column with sus_snps
mrobj_slot$chr_pos <- regmatches(mrobj_slot$variant, regexpr("[a-z0-9]+\\:[a-z0-9]+", mrobj_slot$variant, ignore.case = TRUE))
sus_SNPs_chr_pos <- regmatches(sus_SNPs, regexpr("[a-z0-9]+\\:[a-z0-9]+", sus_SNPs, ignore.case = TRUE))
mrobj_slot$PheWAS_ttest_sensitivity_keep <- ifelse(mrobj_slot$chr_pos %in% sus_SNPs_chr_pos, FALSE, TRUE) # keep the "TRUE"
table(mrobj_slot$PheWAS_ttest_sensitivity_keep)

# USER INPUT - select correct slot
mrobj$alcohol_height <- mrobj_slot

# Clean up - get rid of everything other than mrobj
rm(list = c("df", "exp_b", "exp_n", "exp_pval", "exp_se", "exp_stdB", "exp_stdSE",
            "mrobj_slot", "sig_ind", "stdBeta_df", "stdBeta_df_noEXP", 
            "stdBeta_EXP", "sus_traits_df", "tstat_df", "tstat_dif",
            "tstat_pval_df", "sus_SNPs", "sus_SNPs_chr_pos", "sus_traits"))






#~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# PheWAS PVE results
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~
# Oilyfish - TG
#~~~~~~~~~~~
# USER INPUT - select correct data; used for obtaining filtering labels
load("interim_data/PheWAS_PVE/PheWAS_PVE_interimObjects_oilyfish_GCST90239664_TG_Graham.RData")

# USER INPUT - select correct slot
mrobj_slot <- mrobj$oilyfish_TG

# Filter label
mrobj_slot$chr_pos <- regmatches(mrobj_slot$variant, regexpr("[a-z0-9]+\\:[a-z0-9]+", mrobj_slot$variant, ignore.case = TRUE))
mrobj_slot$PheWAS_PVE_keep <- ifelse(mrobj_slot$chr_pos %in% remove_variants_names, FALSE, TRUE)
table(mrobj_slot$PheWAS_PVE_keep)

# USER INPUT - select correct slot
mrobj$oilyfish_TG <- mrobj_slot

# Clean up - get rid of everything other than mrobj
rm(list = c("filtered_pve", "keep_variants", "mrobj_slot", "phewas_thresh_results", "remove_variants", "keep_variant_names", "remove_variants_names"))

#~~~~~~~~~~~
# Oilyfish - CVD
#~~~~~~~~~~~
# USER INPUT - select correct data; used for obtaining filtering labels
load("interim_data/PheWAS_PVE/PheWAS_PVE_interimObjects_oilyfish_GCST90132314_CVD_Aragam.RData")

# USER INPUT - select correct slot
mrobj_slot <- mrobj$oilyfish_CVD

# Filter label
mrobj_slot$chr_pos <- regmatches(mrobj_slot$variant, regexpr("[a-z0-9]+\\:[a-z0-9]+", mrobj_slot$variant, ignore.case = TRUE))
mrobj_slot$PheWAS_PVE_keep <- ifelse(mrobj_slot$chr_pos %in% remove_variants_names, FALSE, TRUE)
table(mrobj_slot$PheWAS_PVE_keep)

# USER INPUT - select correct slot
mrobj$oilyfish_CVD <- mrobj_slot

# Clean up - get rid of everything other than mrobj
rm(list = c("filtered_pve", "keep_variants", "mrobj_slot", "phewas_thresh_results", "remove_variants", "keep_variant_names", "remove_variants_names"))


#~~~~~~~~~~~
# Oilyfish - Height
#~~~~~~~~~~~
# USER INPUT - select correct data; used for obtaining filtering labels
load("interim_data/PheWAS_PVE/PheWAS_PVE_interimObjects_oilyfish_GIANT_Height_Yengo.RData")

# USER INPUT - select correct slot
mrobj_slot <- mrobj$oilyfish_height

# Filter label
mrobj_slot$chr_pos <- regmatches(mrobj_slot$variant, regexpr("[a-z0-9]+\\:[a-z0-9]+", mrobj_slot$variant, ignore.case = TRUE))
mrobj_slot$PheWAS_PVE_keep <- ifelse(mrobj_slot$chr_pos %in% remove_variants_names, FALSE, TRUE)
table(mrobj_slot$PheWAS_PVE_keep)

# USER INPUT - select correct slot
mrobj$oilyfish_height <- mrobj_slot

# Clean up - get rid of everything other than mrobj
rm(list = c("filtered_pve", "keep_variants", "mrobj_slot", "phewas_thresh_results", "remove_variants", "keep_variant_names", "remove_variants_names"))


#~~~~~~~~~~~
# White X whole grain/brown bread - LDL
#~~~~~~~~~~~
# USER INPUT - select correct data; used for obtaining filtering labels
load("interim_data/PheWAS_PVE/PheWAS_PVE_interimObjects_whiteXwholegrain_GCST90239658_LDL_Graham.RData")

# USER INPUT - select correct slot
mrobj_slot <- mrobj$whiteXwholegrain_LDL

# Filter label
mrobj_slot$chr_pos <- regmatches(mrobj_slot$variant, regexpr("[a-z0-9]+\\:[a-z0-9]+", mrobj_slot$variant, ignore.case = TRUE))
mrobj_slot$PheWAS_PVE_keep <- ifelse(mrobj_slot$chr_pos %in% remove_variants_names, FALSE, TRUE)
table(mrobj_slot$PheWAS_PVE_keep)

# USER INPUT - select correct slot
mrobj$whiteXwholegrain_LDL <- mrobj_slot

# Clean up - get rid of everything other than mrobj
rm(list = c("filtered_pve", "keep_variants", "mrobj_slot", "phewas_thresh_results", "remove_variants", "keep_variant_names", "remove_variants_names"))



#~~~~~~~~~~~
# White X whole grain/brown bread - CVD
#~~~~~~~~~~~
# USER INPUT - select correct data; used for obtaining filtering labels
load("interim_data/PheWAS_PVE/PheWAS_PVE_interimObjects_whiteXwholegrain_GCST90132314_CVD_Aragam.RData")

# USER INPUT - select correct slot
mrobj_slot <- mrobj$whiteXwholegrain_CVD

# Filter label
mrobj_slot$chr_pos <- regmatches(mrobj_slot$variant, regexpr("[a-z0-9]+\\:[a-z0-9]+", mrobj_slot$variant, ignore.case = TRUE))
mrobj_slot$PheWAS_PVE_keep <- ifelse(mrobj_slot$chr_pos %in% remove_variants_names, FALSE, TRUE)
table(mrobj_slot$PheWAS_PVE_keep)

# USER INPUT - select correct slot
mrobj$whiteXwholegrain_CVD <- mrobj_slot

# Clean up - get rid of everything other than mrobj
rm(list = c("filtered_pve", "keep_variants", "mrobj_slot", "phewas_thresh_results", "remove_variants", "keep_variant_names", "remove_variants_names"))


#~~~~~~~~~~~
# White X whole grain/brown bread - Height
#~~~~~~~~~~~
# USER INPUT - select correct data; used for obtaining filtering labels
load("interim_data/PheWAS_PVE/PheWAS_PVE_interimObjects_whiteXwholegrain_height.RData")

# USER INPUT - select correct slot
mrobj_slot <- mrobj$whiteXwholegrain_height

# Filter label
mrobj_slot$chr_pos <- regmatches(mrobj_slot$variant, regexpr("[a-z0-9]+\\:[a-z0-9]+", mrobj_slot$variant, ignore.case = TRUE))
mrobj_slot$PheWAS_PVE_keep <- ifelse(mrobj_slot$chr_pos %in% remove_variants_names, FALSE, TRUE)
table(mrobj_slot$PheWAS_PVE_keep)

# USER INPUT - select correct slot
mrobj$whiteXwholegrain_height <- mrobj_slot

# Clean up - get rid of everything other than mrobj
rm(list = c("filtered_pve", "keep_variants", "mrobj_slot", "phewas_thresh_results", "remove_variants", "keep_variant_names", "remove_variants_names"))


#~~~~~~~~~~~
# Alcohol - ALT
#~~~~~~~~~~~
# USER INPUT - select correct data; used for obtaining filtering labels
load("interim_data/PheWAS_PVE/PheWAS_PVE_interimObjects_alcohol_GCST90013405_ALT_Pazoki.RData")

# USER INPUT - select correct slot
mrobj_slot <- mrobj$alcohol_ALT

# Filter label
mrobj_slot$chr_pos <- regmatches(mrobj_slot$variant, regexpr("[a-z0-9]+\\:[a-z0-9]+", mrobj_slot$variant, ignore.case = TRUE))
mrobj_slot$PheWAS_PVE_keep <- ifelse(mrobj_slot$chr_pos %in% remove_variants_names, FALSE, TRUE)
table(mrobj_slot$PheWAS_PVE_keep)

# USER INPUT - select correct slot
mrobj$alcohol_ALT <- mrobj_slot

# Clean up - get rid of everything other than mrobj
rm(list = c("filtered_pve", "keep_variants", "mrobj_slot", "phewas_thresh_results", "remove_variants", "keep_variant_names", "remove_variants_names"))



#~~~~~~~~~~~
# Alcohol - Liver cirrhosis
#~~~~~~~~~~~
# USER INPUT - select correct data; used for obtaining filtering labels
load("interim_data/PheWAS_PVE/PheWAS_PVE_interimObjects_alcohol_GCST90319877_cir_Ghouse.RData")

# USER INPUT - select correct slot
mrobj_slot <- mrobj$alcohol_cir

# Filter label
mrobj_slot$chr_pos <- regmatches(mrobj_slot$variant, regexpr("[a-z0-9]+\\:[a-z0-9]+", mrobj_slot$variant, ignore.case = TRUE))
mrobj_slot$PheWAS_PVE_keep <- ifelse(mrobj_slot$chr_pos %in% remove_variants_names, FALSE, TRUE)
table(mrobj_slot$PheWAS_PVE_keep)

# USER INPUT - select correct slot
mrobj$alcohol_cir <- mrobj_slot

# Clean up - get rid of everything other than mrobj
rm(list = c("filtered_pve", "keep_variants", "mrobj_slot", "phewas_thresh_results", "remove_variants", "keep_variant_names", "remove_variants_names"))


#~~~~~~~~~~~
# Alcohol - Height
#~~~~~~~~~~~
# USER INPUT - select correct data; used for obtaining filtering labels
load("interim_data/PheWAS_PVE/PheWAS_PVE_interimObjects_alcohol_height.RData")

# USER INPUT - select correct slot
mrobj_slot <- mrobj$alcohol_height

# Filter label
mrobj_slot$chr_pos <- regmatches(mrobj_slot$variant, regexpr("[a-z0-9]+\\:[a-z0-9]+", mrobj_slot$variant, ignore.case = TRUE))
mrobj_slot$PheWAS_PVE_keep <- ifelse(mrobj_slot$chr_pos %in% remove_variants_names, FALSE, TRUE)
table(mrobj_slot$PheWAS_PVE_keep)

# USER INPUT - select correct slot
mrobj$alcohol_height <- mrobj_slot

# Clean up - get rid of everything other than mrobj
rm(list = c("filtered_pve", "keep_variants", "mrobj_slot", "phewas_thresh_results", "remove_variants", "keep_variant_names", "remove_variants_names"))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# PheWAS PVE results - SENSITIVITY PheWAS
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~
# Oilyfish - TG
#~~~~~~~~~~~
# USER INPUT - select correct data; used for obtaining filtering labels
load("interim_data/PheWAS_PVE_Sensitivity/PheWAS_PVE_interimObjects_oilyfish_GCST90239664_TG_Graham_Sensitivity.RData")

# USER INPUT - select correct slot
mrobj_slot <- mrobj$oilyfish_TG

# Filter label
mrobj_slot$chr_pos <- regmatches(mrobj_slot$variant, regexpr("[a-z0-9]+\\:[a-z0-9]+", mrobj_slot$variant, ignore.case = TRUE))
mrobj_slot$PheWAS_PVE_sensitivity_keep <- ifelse(mrobj_slot$chr_pos %in% remove_variants_names, FALSE, TRUE)
table(mrobj_slot$PheWAS_PVE_sensitivity_keep)

# USER INPUT - select correct slot
mrobj$oilyfish_TG <- mrobj_slot

# Clean up - get rid of everything other than mrobj
rm(list = c("filtered_pve", "keep_variants", "mrobj_slot", "phewas_thresh_results", "remove_variants", "keep_variant_names", "remove_variants_names"))

#~~~~~~~~~~~
# Oilyfish - CVD
#~~~~~~~~~~~
# USER INPUT - select correct data; used for obtaining filtering labels
load("interim_data/PheWAS_PVE_Sensitivity/PheWAS_PVE_interimObjects_oilyfish_GCST90132314_CVD_Aragam_Sensitivity.RData")

# USER INPUT - select correct slot
mrobj_slot <- mrobj$oilyfish_CVD

# Filter label
mrobj_slot$chr_pos <- regmatches(mrobj_slot$variant, regexpr("[a-z0-9]+\\:[a-z0-9]+", mrobj_slot$variant, ignore.case = TRUE))
mrobj_slot$PheWAS_PVE_sensitivity_keep <- ifelse(mrobj_slot$chr_pos %in% remove_variants_names, FALSE, TRUE)
table(mrobj_slot$PheWAS_PVE_sensitivity_keep)

# USER INPUT - select correct slot
mrobj$oilyfish_CVD <- mrobj_slot

# Clean up - get rid of everything other than mrobj
rm(list = c("filtered_pve", "keep_variants", "mrobj_slot", "phewas_thresh_results", "remove_variants", "keep_variant_names", "remove_variants_names"))


#~~~~~~~~~~~
# Oilyfish - Height
#~~~~~~~~~~~
# USER INPUT - select correct data; used for obtaining filtering labels
load("interim_data/PheWAS_PVE_Sensitivity/PheWAS_PVE_interimObjects_oilyfish_GIANT_Height_Yengo_Sensitivity.RData")

# USER INPUT - select correct slot
mrobj_slot <- mrobj$oilyfish_height

# Filter label
mrobj_slot$chr_pos <- regmatches(mrobj_slot$variant, regexpr("[a-z0-9]+\\:[a-z0-9]+", mrobj_slot$variant, ignore.case = TRUE))
mrobj_slot$PheWAS_PVE_sensitivity_keep <- ifelse(mrobj_slot$chr_pos %in% remove_variants_names, FALSE, TRUE)
table(mrobj_slot$PheWAS_PVE_sensitivity_keep)

# USER INPUT - select correct slot
mrobj$oilyfish_height <- mrobj_slot

# Clean up - get rid of everything other than mrobj
rm(list = c("filtered_pve", "keep_variants", "mrobj_slot", "phewas_thresh_results", "remove_variants", "keep_variant_names", "remove_variants_names"))


#~~~~~~~~~~~
# White X whole grain/brown bread - LDL
#~~~~~~~~~~~
# USER INPUT - select correct data; used for obtaining filtering labels
load("interim_data/PheWAS_PVE_Sensitivity/PheWAS_PVE_interimObjects_whiteXwholegrain_GCST90239658_LDL_Graham_Sensitivity.RData")

# USER INPUT - select correct slot
mrobj_slot <- mrobj$whiteXwholegrain_LDL

# Filter label
mrobj_slot$chr_pos <- regmatches(mrobj_slot$variant, regexpr("[a-z0-9]+\\:[a-z0-9]+", mrobj_slot$variant, ignore.case = TRUE))
mrobj_slot$PheWAS_PVE_sensitivity_keep <- ifelse(mrobj_slot$chr_pos %in% remove_variants_names, FALSE, TRUE)
table(mrobj_slot$PheWAS_PVE_sensitivity_keep)

# USER INPUT - select correct slot
mrobj$whiteXwholegrain_LDL <- mrobj_slot

# Clean up - get rid of everything other than mrobj
rm(list = c("filtered_pve", "keep_variants", "mrobj_slot", "phewas_thresh_results", "remove_variants", "keep_variant_names", "remove_variants_names"))



#~~~~~~~~~~~
# White X whole grain/brown bread - CVD
#~~~~~~~~~~~
# USER INPUT - select correct data; used for obtaining filtering labels
load("interim_data/PheWAS_PVE_Sensitivity/PheWAS_PVE_interimObjects_whiteXwholegrain_GCST90132314_CVD_Aragam_Sensitivity.RData")

# USER INPUT - select correct slot
mrobj_slot <- mrobj$whiteXwholegrain_CVD

# Filter label
mrobj_slot$chr_pos <- regmatches(mrobj_slot$variant, regexpr("[a-z0-9]+\\:[a-z0-9]+", mrobj_slot$variant, ignore.case = TRUE))
mrobj_slot$PheWAS_PVE_sensitivity_keep <- ifelse(mrobj_slot$chr_pos %in% remove_variants_names, FALSE, TRUE)
table(mrobj_slot$PheWAS_PVE_sensitivity_keep)

# USER INPUT - select correct slot
mrobj$whiteXwholegrain_CVD <- mrobj_slot

# Clean up - get rid of everything other than mrobj
rm(list = c("filtered_pve", "keep_variants", "mrobj_slot", "phewas_thresh_results", "remove_variants", "keep_variant_names", "remove_variants_names"))


#~~~~~~~~~~~
# White X whole grain/brown bread - Height
#~~~~~~~~~~~
# USER INPUT - select correct data; used for obtaining filtering labels
load("interim_data/PheWAS_PVE_Sensitivity/PheWAS_PVE_interimObjects_whiteXwholegrain_height_Sensitivity.RData")

# USER INPUT - select correct slot
mrobj_slot <- mrobj$whiteXwholegrain_height

# Filter label
mrobj_slot$chr_pos <- regmatches(mrobj_slot$variant, regexpr("[a-z0-9]+\\:[a-z0-9]+", mrobj_slot$variant, ignore.case = TRUE))
mrobj_slot$PheWAS_PVE_sensitivity_keep <- ifelse(mrobj_slot$chr_pos %in% remove_variants_names, FALSE, TRUE)
table(mrobj_slot$PheWAS_PVE_sensitivity_keep)

# USER INPUT - select correct slot
mrobj$whiteXwholegrain_height <- mrobj_slot

# Clean up - get rid of everything other than mrobj
rm(list = c("filtered_pve", "keep_variants", "mrobj_slot", "phewas_thresh_results", "remove_variants", "keep_variant_names", "remove_variants_names"))


#~~~~~~~~~~~
# Alcohol - ALT
#~~~~~~~~~~~
# USER INPUT - select correct data; used for obtaining filtering labels
load("interim_data/PheWAS_PVE_Sensitivity/PheWAS_PVE_interimObjects_alcohol_GCST90013405_ALT_Pazoki_Sensitivity.RData")

# USER INPUT - select correct slot
mrobj_slot <- mrobj$alcohol_ALT

# Filter label
mrobj_slot$chr_pos <- regmatches(mrobj_slot$variant, regexpr("[a-z0-9]+\\:[a-z0-9]+", mrobj_slot$variant, ignore.case = TRUE))
mrobj_slot$PheWAS_PVE_sensitivity_keep <- ifelse(mrobj_slot$chr_pos %in% remove_variants_names, FALSE, TRUE)
table(mrobj_slot$PheWAS_PVE_sensitivity_keep)

# USER INPUT - select correct slot
mrobj$alcohol_ALT <- mrobj_slot

# Clean up - get rid of everything other than mrobj
rm(list = c("filtered_pve", "keep_variants", "mrobj_slot", "phewas_thresh_results", "remove_variants", "keep_variant_names", "remove_variants_names"))



#~~~~~~~~~~~
# Alcohol - Liver cirrhosis
#~~~~~~~~~~~
# USER INPUT - select correct data; used for obtaining filtering labels
load("interim_data/PheWAS_PVE_Sensitivity/PheWAS_PVE_interimObjects_alcohol_GCST90319877_cir_Ghouse_Sensitivity.RData")

# USER INPUT - select correct slot
mrobj_slot <- mrobj$alcohol_cir

# Filter label
mrobj_slot$chr_pos <- regmatches(mrobj_slot$variant, regexpr("[a-z0-9]+\\:[a-z0-9]+", mrobj_slot$variant, ignore.case = TRUE))
mrobj_slot$PheWAS_PVE_sensitivity_keep <- ifelse(mrobj_slot$chr_pos %in% remove_variants_names, FALSE, TRUE)
table(mrobj_slot$PheWAS_PVE_sensitivity_keep)

# USER INPUT - select correct slot
mrobj$alcohol_cir <- mrobj_slot

# Clean up - get rid of everything other than mrobj
rm(list = c("filtered_pve", "keep_variants", "mrobj_slot", "phewas_thresh_results", "remove_variants", "keep_variant_names", "remove_variants_names"))


#~~~~~~~~~~~
# Alcohol - Height
#~~~~~~~~~~~
# USER INPUT - select correct data; used for obtaining filtering labels
load("interim_data/PheWAS_PVE_Sensitivity/PheWAS_PVE_interimObjects_alcohol_height_Sensitivity.RData")

# USER INPUT - select correct slot
mrobj_slot <- mrobj$alcohol_height

# Filter label
mrobj_slot$chr_pos <- regmatches(mrobj_slot$variant, regexpr("[a-z0-9]+\\:[a-z0-9]+", mrobj_slot$variant, ignore.case = TRUE))
mrobj_slot$PheWAS_PVE_sensitivity_keep <- ifelse(mrobj_slot$chr_pos %in% remove_variants_names, FALSE, TRUE)
table(mrobj_slot$PheWAS_PVE_sensitivity_keep)

# USER INPUT - select correct slot
mrobj$alcohol_height <- mrobj_slot

# Clean up - get rid of everything other than mrobj
rm(list = c("filtered_pve", "keep_variants", "mrobj_slot", "phewas_thresh_results", "remove_variants", "keep_variant_names", "remove_variants_names"))




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Save
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#saveRDS(mrobj, "interim_data/merged_GI/MRobject_steiger_PheWAStt_PheWASPVE_241120.rds")
saveRDS(mrobj, "interim_data/merged_GI/MRobject_steiger_PheWAStt_PheWASPVE_sensitivity_241203.rds")

