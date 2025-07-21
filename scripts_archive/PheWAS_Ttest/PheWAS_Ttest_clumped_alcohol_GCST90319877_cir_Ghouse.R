# MR Analysis: PheWAS-Cluster - updated with new proxies for exposure
# Author: KJS
# Date: 11/18/24 

# Exposure: alcohol (from clustering)
# Outcome: GCST90319877_cir_Ghouse

# Pre-requisite: 
# - Get PheWAS data from all exposure SNPs
# - Have GI for exposure/outcome

# Overview of steps:
# 0. Set names - USER INPUT NEEDED
# 1. Load PheWAS data
# 2. Filter PheWAS traits to those with an effective sample size of 50K - USER INPUT NEEDED
# 3. Filter PheWAS data to traits that have a genetic correlation of less than 75% with the exposure
# 4. Load exposure data 
# 5. Wrangle exposure data to get a variant x statistic matrix for n, se, pval, b, t-stat
# 6. Wrangle PheWAS data to get a variant x statistic matrix for n, se, pval, b, t-stat 
# 7. Prep for PheWAS filtering
# 8. Apply PheWAS filtering
# 9. Perform MR
# 10. Plot
# 11. SKIP - Modify MRobject in its own script 
# 12. Investigate SNPs-traits removed

# Set wd
wdir <- "/pl/active/colelab/users/kjames/refinedMR/"
setwd(wdir)

# Load libraries
library(data.table)
library(dplyr)
library(ggplot2)
library(ggrepel)

#~~~~~~~~~~~~~~~~~~~~~
# 0. Set names
#~~~~~~~~~~~~~~~~~~~~~

# Set traits for naming files 
exposure <- "alcohol"
outcome <- "GCST90319877_cir_Ghouse"
outcome_short <- "cir"
OUT_pheno <- outcome # used for some naming
# EXP_pheno defined below

# For plots
EXP_name <- "Alcohol"
OUT_name <- "Liver Cirrhosis"

#~~~~~~~~~~~~~~~~~~~~~
# 1. Load PheWAS data
#~~~~~~~~~~~~~~~~~~~~~
pwas <- read.csv("interim_data/PheWAS_Grep_SNPs/PheWAS_Oilyfish_WhiteXWhole_Alcochol_UpdatedProxies.csv")
pwas <- pwas %>% select(-c(X, V1))

# Make file name prettier. Get rid of extra junk.
pwas$trait <- sub("\\..*$", "", pwas$file_name)

# PWAS includes some duplicated traits that are in both raw and 
# inverse rank normalized transformed states. 
# In these cases, we want to keep the IRNT versions. 
# Remove the traits from pwas that end in _raw
pwas <- pwas[!grepl("_raw$", pwas$trait),]

#~~~~~~~~~~~~~~~~~~~~~
# 2. Filter PheWAS traits to those with an effective sample size of 50K
# For binary traits:
# Neff = (4 x Ncases x Ncontrols)/(Ncases + Ncontrols)
# For quantitative traits: 
# Neff = samplesize [... so this step is geared towards binary traits]
#~~~~~~~~~~~~~~~~~~~~~
# NOTE - THE quote="" is VERY IMPORTANT; without it, it will truncate your file!
md <- read.table(gzfile("/pl/active/colelab/common/published_gwas/Neale_UKB_GWAS_round2/Both_sexes/metadata/phenotypes.both_sexes.v2.tsv.bgz"),header=TRUE, sep="\t", fill=TRUE, quote = "")

# Look at numbers of case/control, binary, ordinal, continuous traits
# Compare to this numbers in "imputed-v3 Phenotypes" section here, https://github.com/Nealelab/UK_Biobank_GWAS?tab=readme-ov-file#imputed-v3-association-model
# binary traits are high at 3531 vs the expected 2346 from Phesant and 633 from ICD10 codes
table(md$variable_type)

#~~~~~~~~~~~~~~~~~~~~~
# 2.1. Filter based on general sample size and ordinal variable sample size
#~~~~~~~~~~~~~~~~~~~~~
# from, https://github.com/LizaDarrous/PheWAS-cluster/blob/main/R/getSNPs_PheWAS.R
n_thresh <- 50000
md$n_eff = md$n_non_missing
ncont_ind = which(!is.na(md$n_cases)) #which rows have values in n_cases 
md$n_eff[ncont_ind] = (4*md$n_cases[ncont_ind]*md$n_controls[ncont_ind]) / (md$n_cases[ncont_ind]+md$n_controls[ncont_ind]) # (4*(control*case))/totN
small_n = which(md$n_eff < n_thresh)
md = md[-small_n,]

# Only keep traits in PheWAS that have sample size >50K
pwas2 <- pwas %>% filter((trait %in% md$phenotype))
length(unique(pwas2$trait)) #829 unique traits

# Save list of trait IDs
pwas2_trait_IDs <- unique(pwas2$trait)
write.csv(pwas2_trait_IDs, "/pl/active/colelab/users/kjames/refinedMR/interim_data/PheWAS_Grep_SNPs/pwas2_trait_ids.csv", row.names = FALSE, quote = FALSE)
#~~~~~~~~~~~~~~~~~~~~~
# 2.2. Merge pwas data with other meta data
#~~~~~~~~~~~~~~~~~~~~~
# Get the manifest key for PheWAS
key <- read.csv("/pl/active/colelab/common/published_gwas/Neale_UKB_GWAS_round2/metadata/UKBB GWAS Imputed v3 - File Manifest Release 20180731 - Description Lookup.csv", header = FALSE)

# Merge to make a description column. Also, name the data frame traits to match LD's workflow.
trait_info <- merge(pwas2, key, by.x = "trait", by.y = "V1", all.x = TRUE, all.y = FALSE)
# new_name = old_name
trait_info = trait_info %>% dplyr::rename("description" = "V2")

# !ATTENTION!
# USER INPUT (used in section XXX)
# Find the ID of the initial exposures of interest
key_term <- "Alcohol consumed" #"Bread consumed" #"Oily fish intake"
trait_ID <- key[grep(key_term, key$V2),] 
trait_ID <- trait_ID[1,1]
#~~~~~~~~~~~~~~~~~~~~~
# 3. Filter PheWAS data to traits that have a genetic correlation of less than 75% with the exposure
#~~~~~~~~~~~~~~~~~~~~~
# Genetic correlation data
gcor <- read.csv("/pl/active/colelab/users/kjames/enviroMR/inputs/UKBB_genetic_correlation_data-2024-03-20.csv")

# How many of the pwas traits are in the gcor data?
gcor_pwas_sub <- pwas2[pwas2$trait %in% gcor$ID1,]
length(unique(gcor_pwas_sub$trait)) #404, cite this as a limitation to this qc filter

# Make phenotype names easier to read using the column with the url
gcor$phenotype.1.name <- stringr::str_extract(gcor$Phenotype.1, "(?<=\\>).*(?=\\<)")
gcor$phenotype.2.name <- stringr::str_extract(gcor$Phenotype.2, "(?<=\\>).*(?=\\<)")

# USER INPUT: Define genetic correlation threshold
gcorr_thresh <- 0.75
# USER INPUT: Define trait of interest
cat(key_term) # check 
toi_gcorr <- gcor[grep(key_term, gcor$phenotype.1.name),]
# Filter
high_gcorr = which(abs(toi_gcorr$rg)>=gcorr_thresh) #note, gcorr are + and -

# Function that checks the genetic correlations
check_gcorr <- function(high_gcorr) {
  if (length(high_gcorr) > 0) {
    cat(paste0("Remove these phenotypes:\n", high_gcorr))
  } else {
    cat(paste0("There are no phenotypes with a genetic correlation greater than:\n", gcorr_thresh ))
  }
}

# Apply function
check_gcorr(high_gcorr)

# Set up: 
toi_gcorr_rm <- toi_gcorr[high_gcorr,]
trait_rm <- toi_gcorr_rm$ID2

# Get rows of pwas2 that should be removed due to gcorr
rows <- pwas2$trait == trait_rm
# Remove from pwas2; does nothing if there are none
if (length(rows) == 0) {
  # Do nothing
} else {
  # Remove rows with indices
  pwas2=pwas2[!rows,]
}

# Get rows of trait_info that should be removed due to gcorr
rows <- trait_info$trait == trait_rm
# Remove from trait info; does nothing if there are none
if (length(rows) == 0) {
  # Do nothing
} else {
  # Remove rows with indices
  trait_info=trait_info[!rows,]
}

#~~~~~~~~~~~~~~~~~~~~~
# 4. Load exposure data
#~~~~~~~~~~~~~~~~~~~~~
# Load merged, harmonized GI data
mrobj <- readRDS("interim_data/merged_GI/MRobject_steiger_241113.rds")

# Programatically get the right component of the list
exposure_outcome <- paste0(exposure, "_", outcome_short)
df <- mrobj[[exposure_outcome]]

#~~~~~~~~~~~~~~~~~~~~~
# 5. Wrangle exposure data to get a variant x statistic matrix for n, se, pval, b, t-stat
#~~~~~~~~~~~~~~~~~~~~~
# Bread + LDL - have a duplicate SNP, keep first
#duplicates <- df %>% group_by(variant) %>% filter(n() > 1) %>% ungroup() #affecting 2 SNPs rs35811711, rs66495454
#print(duplicates)
# Filter, keep first
#df <- df %>% distinct(variant, .keep_all = TRUE)

# Set row names
rownames(df) <- df$variant
# n
exp_n <- df[,"samplesize.exposure", drop=FALSE]
# se
exp_se <- df[,"se.exposure", drop=FALSE]
# pval
exp_pval <- df[,"pval.exposure", drop=FALSE]
# beta
exp_b <- df[,"beta.exposure", drop=FALSE]

#~~~~~~~~~~~~~~~~~~~~~
# 6. Wrangle PheWAS data to get a variant x statistic matrix for n, se, pval, b, t-stat 
#~~~~~~~~~~~~~~~~~~~~~
# Filter pwas with chr_pos info (without allele info)
pwas2$chr_pos <- regmatches(pwas2$variant, regexpr("[a-z0-9]+\\:[a-z0-9]+", pwas2$variant, ignore.case = TRUE)) 
exp_n$chr_pos <- regmatches(rownames(exp_n), regexpr("[a-z0-9]+\\:[a-z0-9]+", rownames(exp_n), ignore.case = TRUE)) 
pwas2_sub <- pwas2[pwas2$chr_pos %in% exp_n$chr_pos,]

# Remove chr_pos variable from exp_n
exp_n <- exp_n %>% select(-chr_pos)

# Identify duplicates, if any. If present, add code to remove.
dups = pwas2_sub %>%
  dplyr::group_by(variant, file_name) %>%
  dplyr::summarise(n = dplyr::n(), .groups = "drop") %>%
  dplyr::filter(n > 1L)

# If dups is greater than 0, run the next chunk of code
# 100001 is a repeat
# Identify and retain only the first occurrence of each duplicate based on `variant` and `trait`
pwas2_sub <- pwas2_sub %>%
  mutate(dup_flag = duplicated(pwas2_sub %>% select(variant, trait))) %>%
  filter(!duplicated(select(., variant, trait))) 

# Get variant x n matrix
pwas_sub_n <- tidyr::pivot_wider(pwas2_sub, id_cols = c(variant), names_from = file_name, values_from = n_complete_samples)

# Get variant x se matrix
pwas_sub_se <- tidyr::pivot_wider(pwas2_sub, id_cols = c(variant), names_from = file_name, values_from = se)

# Get variant x pval matrix
pwas_sub_pval <- tidyr::pivot_wider(pwas2_sub, id_cols = c(variant), names_from = file_name, values_from = pval)

# Get variant x b matrix
pwas_sub_beta <- tidyr::pivot_wider(pwas2_sub, id_cols = c(variant), names_from = file_name, values_from = beta)

# Get variant x t-statistic matrix
pwas_sub_tstat <- tidyr::pivot_wider(pwas2_sub, id_cols = c(variant), names_from = file_name, values_from = tstat)

# Check that row names and column names are the same for each of the data frames
identical(colnames(pwas_sub_n), colnames(pwas_sub_se))
identical(rownames(pwas_sub_n), rownames(pwas_sub_se))
# TRUE

#~~~~~~~~~~~~~~~~~~~~~
# 7. Prep for PheWAS filtering
#~~~~~~~~~~~~~~~~~~~~~
# Name files to match LD objects
unstdBeta_df <- pwas_sub_beta
unstdSE_df <- pwas_sub_se
tstat_df <- pwas_sub_tstat
pval_df <- pwas_sub_pval

# filter out traits that are NA (albumin and other measurements)
na_traits = which(is.na(trait_info$description)==T)
# not all are also NA in beta file
na_effect = which(is.na(unstdBeta_df)==T, arr.ind=T)
na_both = intersect(na_traits,unique(as.numeric(na_effect[,2])));
if(length(na_both)>0){
  unstdBeta_df = unstdBeta_df[,-na_both]
  unstdSE_df = unstdSE_df[,-na_both]
  tstat_df = tstat_df[,-na_both]
  pval_df = pval_df[,-na_both]
  trait_info = trait_info[-na_both,]
}

# Make row names objects (if you set them now they will get erased)
rnames_unstdBeta_df <- unstdBeta_df$variant
rnames_unstdSE_df <- unstdSE_df$variant
rnames_tstat_df <- tstat_df$variant
rnames_pval_df <- pval_df$variant

# Get rid of variant column in first position. Will set row names after math.
unstdBeta_df <- unstdBeta_df %>% select(-variant)
unstdSE_df <- unstdSE_df %>% select(-variant)
tstat_df <- tstat_df %>% select(-variant)
pval_df <- pval_df %>% select(-variant)



#Find the ID of the initial exposures of interest
EXP_pheno <-  trait_ID

# Find indices of NA values
na_indices <- which(is.na(tstat_df), arr.ind = TRUE)
na_indices_uniq <- unique(na_indices[,2]) 
# Remove NA
# If the value is null/NA and there are indices, then filter
if (is.null(na_indices_uniq) && length(na_indices_uniq) > 0) {
  unstdBeta_df = unstdBeta_df[ , -na_indices_uniq]
}

if (is.null(na_indices_uniq) && length(na_indices_uniq) > 0) {
  unstdSE_df = unstdSE_df[ , -na_indices_uniq]
}

if (is.null(na_indices_uniq) && length(na_indices_uniq) > 0) {
  tstat_df = tstat_df[ , -na_indices_uniq]
}

if (is.null(na_indices_uniq) && length(na_indices_uniq) > 0) {
  pval_df = pval_df[ , -na_indices_uniq]
}

if (is.null(na_indices_uniq) && length(na_indices_uniq) > 0) {
  trait_info4 = trait_info4[-na_indices_uniq,]
}

print(paste0("Dimension of unstd SNP-trait matrix: ", dim(unstdBeta_df)[1],", ",dim(unstdBeta_df)[2]))

#~~~~~~~~~~~~~~~~~~~~~
# 8. Apply PheWAS filtering
#~~~~~~~~~~~~~~~~~~~~~
#11/18/24 - adjusted instances of trait_info4 to use pwas2_sub data
identical(colnames(tstat_df), colnames(pwas_sub_n[,-1]))
phewas_n <- as.data.frame(t(pwas_sub_n[1,-1]))
phewas_n <- phewas_n %>% rename("N"="V1")

# Getting the standardized effects (tstat/sqrt(N)) 
# Performs matrix multiplication
stdBeta_df = as.matrix(tstat_df) %*% diag(1/(sqrt(phewas_n$N)))   #mat %*% diag(1 / dev)

# Manual check to test if you understand matrix multiplication step
tstat_df[1,1]/sqrt(phewas_n$N[1]) # gives the same value as stdBeta_df, good check
# Manual check to test if you understand matrix multiplication step
tstat_df[1,1] * 1/sqrt(phewas_n[1,1])
stdBeta_df[1,1]
# These are the same value. 
# So, yes, the standardization for beta is tstat/sqrt(N). 

# Make a data frame and reassign names
stdBeta_df <- as.data.frame(stdBeta_df)
colnames(stdBeta_df) = colnames(unstdBeta_df)
rownames(stdBeta_df) = rnames_tstat_df

# Name and standardize non-exposure containing objects
stdBeta_df_noEXP = stdBeta_df # if the exposure is not in the PheWAS then these objects are the same
stdSE_noEXP = 1/(sqrt(phewas_n$N)) #1/sqrt(N)

# The SNPs in your exposure data and the PheWAS data are supposed to be the same
# But sometimes there are minor differences. This causes a problem in the next step
# Filter to make the pwas SNPs and the exposure SNPs overlap completely
exp_n$chr_pos <- regmatches(rownames(exp_n), regexpr("[a-z0-9]+\\:[a-z0-9]+", rownames(exp_n), ignore.case = TRUE)) 
exp_b$chr_pos <- regmatches(rownames(exp_b), regexpr("[a-z0-9]+\\:[a-z0-9]+", rownames(exp_b), ignore.case = TRUE)) 
exp_pval$chr_pos <- regmatches(rownames(exp_pval), regexpr("[a-z0-9]+\\:[a-z0-9]+", rownames(exp_pval), ignore.case = TRUE)) 
exp_se$chr_pos <- regmatches(rownames(exp_se), regexpr("[a-z0-9]+\\:[a-z0-9]+", rownames(exp_se), ignore.case = TRUE)) 

stdBeta_df_noEXP$chr_pos <- regmatches(rownames(stdBeta_df_noEXP), regexpr("[a-z0-9]+\\:[a-z0-9]+", rownames(stdBeta_df_noEXP), ignore.case = TRUE)) 

exp_n <- exp_n %>% filter(chr_pos %in% stdBeta_df_noEXP$chr_pos)
exp_b <- exp_b %>% filter(chr_pos %in% stdBeta_df_noEXP$chr_pos)
exp_pval <- exp_pval %>% filter(chr_pos %in% stdBeta_df_noEXP$chr_pos)
exp_se <- exp_se %>% filter(chr_pos %in% stdBeta_df_noEXP$chr_pos)

# Remove the 'chr_pos' column and set row names back
exp_n <- exp_n %>% select(-chr_pos)
exp_b <- exp_b %>% select(-chr_pos)
exp_pval <- exp_pval %>% select(-chr_pos)
exp_se <- exp_se %>% select(-chr_pos)

# Standardize exposure beta and se - corresponds to your data/exposure
# Beta
exp_stdB <- (exp_b/exp_se)/sqrt(exp_n)
# SE
exp_stdSE <- exp_se/sqrt(exp_n) 

# Rename standardized exposure objects
# Rows = n variants in instrument 
stdBeta_EXP <- exp_stdB
stdSE_EXP <- exp_stdSE


# Order the rows the same between your exposure SNPs and pwas SNPs
# Beta
# Convert rownames to a column for both data frames
stdBeta_EXP$chr_pos <- regmatches(rownames(stdBeta_EXP), regexpr("[a-z0-9]+\\:[a-z0-9]+", rownames(stdBeta_EXP), ignore.case = TRUE)) 
stdBeta_df_noEXP$chr_pos <- regmatches(rownames(stdBeta_df_noEXP), regexpr("[a-z0-9]+\\:[a-z0-9]+", rownames(stdBeta_df_noEXP), ignore.case = TRUE)) 
# arrange
stdBeta_EXP_ordered = stdBeta_EXP %>% arrange(match(chr_pos, stdBeta_df_noEXP$chr_pos))

# Check that it worked
identical(stdBeta_EXP_ordered$chr_pos, stdBeta_df_noEXP$chr_pos) # TRUE - stop if not true!

# Remove the 'rowname' column and set row names back
stdBeta_EXP <- stdBeta_EXP %>% select(-chr_pos)
stdBeta_df_noEXP <- stdBeta_df_noEXP %>% select(-chr_pos)


# SE
# use rownames of stdBeta (exposure) to check variants as sdSE_noEXP corresponds to each PheWAS trait
#stdSE_EXP <- stdSE_EXP[which(rownames(stdSE_EXP) %in% rownames(stdBeta_df_noEXP)),, drop = FALSE]
identical(rownames(stdSE_EXP), rownames(stdBeta_EXP)) # TRUE - stop if not true!
# PWAS data
stdSE_EXP <- 1/(sqrt(exp_n$samplesize.exposure [1])) #1/sqrt(N)

# Get rid of exposure from df
#trait_info_noEXP = trait_info4[-which(grepl(paste0("^",EXP_pheno), trait_info4$file_name)) ,] #its not in there
print(paste0("Dimension of std SNP-trait matrix: ", dim(stdBeta_df)[1],", ",dim(stdBeta_df)[2]))
print(paste0("Number of traits: ", dim(pwas_sub_n)[2]))

# filtering for SNPs more associated with other traits using t-stat
tstat_dif = matrix(NA, nrow = nrow(stdBeta_df_noEXP), ncol = ncol(stdBeta_df_noEXP)) #makes empty matrix with specified dimensions
for(i in 1:nrow(stdBeta_df_noEXP)){ #1:i SNPs
  for(j in 1:ncol(stdBeta_df_noEXP)){ #1:j traits
    tstat_dif[i,j] = (abs(stdBeta_EXP[i,]) - abs(stdBeta_df_noEXP[i,j])) / #trait of interest, snp i - trait j, snp i 
      sqrt(stdSE_EXP^2 + stdSE_noEXP[j]^2)
  }
}
# beta of exposure corresponding to SNP i - <minus> beta of trait j corresponding to SNP i 
# / <divided by> 
# square root of (the standard error of the exposure squared + <plus> the standard error of trait j squared)
# result is a matrix with the tstat_difference

# Set up t-test
tstat_pval_df = pnorm(tstat_dif)  # larger than on EXP
tstat_thresh = 0.05/dim(phewas_n)[1] # accuracy: divided by nSNPs* nIndp.traits - multiple test correction step

# Format
tstat_dif = as.data.frame(tstat_dif) # turn matrix into df
colnames(tstat_dif) = colnames(stdBeta_df_noEXP)
rownames(tstat_dif) = rownames(stdBeta_df_noEXP)

# Apply test/filter
sig_ind = which(tstat_pval_df<tstat_thresh, arr.ind = T) #return array indices (df with columns named: row, column)

# Identify SNPS and traits that are more related to other traits
# Which SNPs?
sus_SNPs=rownames(tstat_dif)[unique(sig_ind[,1])] #get the row names of SNPs from tstat_dif that are less than the threshold (indicating more related to snp than exposure)

# Which traits?
sus_traits = unique(rownames(phewas_n)[(sig_ind[,2])]) 
sus_traits_df = as.data.frame(sus_traits)
sus_traits_df$trait = gsub("\\..*$", "", sus_traits_df$sus_traits)
sus_traits_df = merge(sus_traits_df, key, by.x = "trait", by.y = "V1", all.x = TRUE, all.y = FALSE )
print(paste0("SNP filtering revealed ", length(sus_SNPs)," SNPs more strongly associated to traits other than exposure trait ", EXP_pheno,"."))

# Print results?
print(paste0("The traits more with SNPs more strongly associated with traits other than the exposure include:", sus_traits_df$V2))
print(paste0("The SNPs more strongly associated with traits other than the exposure include:", sus_SNPs))  


# Save output objects
interim_obj_path <- "/pl/active/colelab/users/kjames/refinedMR/interim_data/PheWAS_Cluster/"
save(list=c("exp_b", "exp_n", "exp_pval", "exp_se", "exp_stdB", "exp_stdSE", 
            "stdBeta_df", "stdBeta_df_noEXP", "stdBeta_EXP", #"trait_info4",
            "tstat_df", "tstat_dif", "tstat_pval_df",
            "df", "sus_SNPs", "sus_traits","sig_ind", "sus_traits_df"),
     file = paste0(interim_obj_path,"/PheWASCluster_objects_",exposure,"_", outcome,".Rdata"))

# Load if needed
load(file = paste0(interim_obj_path,"/PheWASCluster_objects_",exposure,"_", outcome,".Rdata"))
#~~~~~~~~~~~~~~~~~~~~~
# 9. Perform MR
#~~~~~~~~~~~~~~~~~~~~~
# Get variant column
head(df)

# Try filtering with just chr:pos
df$chr_pos <- regmatches(df$variant, regexpr("[a-z0-9]+\\:[a-z0-9]+", df$variant, ignore.case = TRUE)) 
sus_SNPs_chr_pos <- regmatches(sus_SNPs, regexpr("[a-z0-9]+\\:[a-z0-9]+", sus_SNPs, ignore.case = TRUE)) #14

# Filter using this variable
df_filt <- df[!df$chr_pos %in% sus_SNPs_chr_pos,]
dim(df_filt)

# Raw, unfiltered
# This step can be used to compare with the unstandardized results
MRobj_raw <- MendelianRandomization::mr_input(bx = df_filt$beta.exposure,
                                              bxse = df_filt$se.exposure,
                                              by = df_filt$beta.outcome,
                                              byse = df_filt$se.outcome)

MRres_raw <- MendelianRandomization::mr_allmethods(MRobj_raw)
MRres_raw

# Get df
MRres_df <- MRres_raw@Values

# Save tables
# Path
phewas_table_dir <- "/pl/active/colelab/users/kjames/refinedMR/results/PheWAS_Cluster/tables"
# All MR results
write.csv(MRres_df, paste0(phewas_table_dir, "/MR_PheWASCluster_Filtered_Results_", exposure,"_", outcome, ".csv"))

# Run IVW 
mr_ivw <- MendelianRandomization::mr_ivw(MRobj_raw)
# Create a table of results for saving
mr_ivw_df <- data.frame(test = "mr_ivw",
                        model = mr_ivw@Model,
                        robust = mr_ivw@Robust,
                        penalized = mr_ivw@Penalized,
                        estimate = mr_ivw@Estimate,
                        std_error = mr_ivw@StdError,
                        cilower = mr_ivw@CILower,
                        ciupper = mr_ivw@CIUpper,
                        pvalue = mr_ivw@Pvalue,
                        heter_Q = mr_ivw@Heter.Stat[1],
                        heter_pvalue = mr_ivw@Heter.Stat[2],
                        fstat = mr_ivw@Fstat)

mr_ivw_df

# Save table
write.csv(mr_ivw_df, paste0(phewas_table_dir, "/MR_PheWASCluster_Filtered_IVWres_", exposure,"_", outcome, ".csv"))

#~~~~~~~~~~~~~~~~~~~~~
# 10. Plot
#~~~~~~~~~~~~~~~~~~~~~
# Plot with all methods
MRres_p <- MendelianRandomization::mr_plot(MRres_raw, error=TRUE, orientate=FALSE, line="ivw", labels=TRUE) +
  ggplot2::theme_bw() +
  ggplot2::theme(
    #legend.position="none", 
    axis.title.y=ggplot2::element_text(size=14),
    axis.text.y=ggplot2::element_text(size=14), 
    axis.ticks.y=ggplot2::element_line(linewidth=0),
    axis.title.x=ggplot2::element_text(size=14),
    axis.text.x=ggplot2::element_text(size=14),
    axis.ticks.x=ggplot2::element_line(linewidth=0),
    panel.border = element_rect(color = "black", 
                                fill = NA, 
                                linewidth = 3)) +
  labs(x = paste0("Genetic Effect on ", EXP_name),
       y = paste0( "Genetic Effect on ", OUT_name)) +
  ggtitle(paste0("PheWAS-Cluster Filtered Variants (n=", nrow(df_filt), ")")) 
#scale_color_viridis(discrete = TRUE)

MRres_p

# USER INPUT: Check paths and variable names
plots_dir <- "/pl/active/colelab/users/kjames/refinedMR/results/PheWAS_Cluster/plots"

# Save
ggsave(paste0(plots_dir,"/PheWAS_Cluster_MR_",exposure, "_", outcome,"_scatterplot_allmethods.png"),
       plot = MRres_p,
       device = NULL,
       path = NULL,
       scale = 1,
       width = 6,
       height = 4,
       units = c("in"),
       dpi = 300,
       limitsize = TRUE)

# Get values and subset to weighted median and IVW for custom plot below
mrRes_df_subsetMethods <- MRres_df %>% filter(Method %in% c("Weighted median", "IVW"))

# Custom plot for external facing meetings
custom_MRscatter <- ggplot(df_filt, aes(x = beta.exposure, y = beta.outcome)) +
  geom_point() +
  ggplot2::geom_errorbar(ggplot2::aes(ymin=beta.outcome-se.outcome, ymax=beta.outcome+se.outcome), colour="grey", width=0) +
  ggplot2::geom_errorbarh(ggplot2::aes(xmin=beta.exposure-se.exposure, xmax=beta.exposure+se.exposure), colour="grey", height=0) +
  ggplot2::geom_abline(data = mrRes_df_subsetMethods, # Specify correct data here
                       ggplot2::aes(intercept=0, slope=Estimate, colour=Method), size = 1.75, show.legend=TRUE) + #0 if not MR-Egger
  ggplot2::theme_bw() +
  ggplot2::theme(
    #legend.position="none", 
    axis.title.y=ggplot2::element_text(size=14),
    axis.text.y=ggplot2::element_text(size=14), 
    axis.ticks.y=ggplot2::element_line(size=0),
    axis.title.x=ggplot2::element_text(size=14),
    axis.text.x=ggplot2::element_text(size=14),
    axis.ticks.x=ggplot2::element_line(size=0),
    panel.border = element_rect(color = "black", 
                                fill = NA, 
                                size = 3)) +
  labs(x = paste0("Genetic Effect on ", EXP_name),
       y = paste0( "Genetic Effect on ", OUT_name)) +
  ggtitle(paste0("PheWAS-Cluster Genetic Instrument (n=", nrow(df_filt), ")")) +
  scale_color_manual(values = c("#3E4A89FF", "#6DCD59FF"))

custom_MRscatter

ggsave(paste0(plots_dir,"/PheWAS_Cluster_MR_",exposure, "_", outcome,"_scatterplot_IVW_weightedmdn.png"),
       plot = custom_MRscatter,
       device = NULL,
       path = NULL,
       scale = 1,
       width = 6,
       height = 4,
       units = c("in"),
       dpi = 300,
       limitsize = TRUE)

#~~~~~~~~~~~~~~~~~~~~~
# 11. Load MR object, add PheWAS Cluster screen to slot
#~~~~~~~~~~~~~~~~~~~~~
# Get master object
#mrobj <- readRDS("interim_data/merged_GI/standardized_shared/std_mrobject_steiger_240918.rds")

# USER INPUT - select correct slot
# Get oily fish TG data
#mrobj_slot <- mrobj$oilyfish_TG

# Make column with sus_snps
#mrobj_slot$chr_pos <- regmatches(mrobj_slot$variant.exposure, regexpr("[a-z0-9]+\\:[a-z0-9]+", mrobj_slot$variant.exposure, ignore.case = TRUE))
#sus_SNPs_chr_pos <- regmatches(sus_SNPs, regexpr("[a-z0-9]+\\:[a-z0-9]+", sus_SNPs, ignore.case = TRUE))
#robj_slot$PheWAS_clust_keep <- ifelse(mrobj_slot$chr_pos %in% sus_SNPs_chr_pos, FALSE, TRUE) # keep the "TRUE"
#table(mrobj_slot$PheWAS_clust_keep)

# USER INPUT - select correct slot
#mrobj$oilyfish_TG <- mrobj_slot

# Save
#saveRDS(mrobj, "interim_data/merged_GI/standardized_shared/std_mrobject_steiger_PheWASClust_240924.rds")

#~~~~~~~~~~~~~~~~~~~~~
# 12. Investigate SNPs removed
#~~~~~~~~~~~~~~~~~~~~~
# Get the row names of SNPs from tstat_dif that are less than the threshold (indicating more related to snp than exposure)
rmd_SNPnames=rownames(tstat_dif)[sig_ind[,1]] 
sig_ind_readable= cbind(rmd_SNPnames, sig_ind)
# Get the colnames names of traits from tstat_dif that correspond to the SNPs with a stronger relationship than the threshold
rmd_traitnames=colnames(tstat_dif)[sig_ind[,2]]
sig_ind_readable=cbind(sig_ind_readable, rmd_traitnames)

# Save
# Path
phewas_table_dir <- "/pl/active/colelab/users/kjames/refinedMR/results/PheWAS_Cluster/tables"
write.csv(sig_ind_readable, paste0(phewas_table_dir, "/MR_PheWASCluster_SusSNPs_RelatedTraits_FilePaths_", exposure,"_", outcome, ".csv"))

#Get description
traits_unique <- trait_info %>% select(file_name, description)
traits_unique <- unique(traits_unique)
sig_ind_readable = merge(sig_ind_readable, traits_unique, by.x = "rmd_traitnames", by.y = "file_name", all.x = TRUE, all.y = FALSE)

# Format
sig_ind_readable <- sig_ind_readable %>% 
  select(rmd_SNPnames, description) %>% 
  dplyr::rename("variant" = "rmd_SNPnames")

# Add column to remind you the exposure and outcome
sig_ind_readable$exposure <- exposure
sig_ind_readable$outcome <- outcome

# Save tables
# Path
phewas_table_dir <- "/pl/active/colelab/users/kjames/refinedMR/results/PheWAS_Cluster/tables"
# MR results
write.csv(sig_ind_readable, paste0(phewas_table_dir, "/MR_PheWASCluster_SusSNPs_RelatedTraits_", exposure,"_", outcome, ".csv"))


