# Author: KJS
# Date: June 18, 2025
# Project: Refining diet MR

# Goal: Make the GI for MVMR using the clumped files from Plink.

# Load libraries
library(dplyr)

#Set working direcotry
setwd("/pl/active/colelab/users/kjames/refinedMR/")

# SET TRAITS MANUALLY
EXP_pheno <- "oilyfish"
OUT_pheno <- "TG"
diet_trait <- paste0(EXP_pheno, "_", OUT_pheno)

# Load the clumped results in 
#tst <- read.table("interim_data/ldsc/plink_clumped_1KGP_HG19/oilyfish_tg/clumped_results_4124_irnt_chr1.tsv.clumped", header = TRUE)

# Load programatically. HOWEVER, PHEWAS_TRAIT MANUALLY
phewas_trait <- "30110_irnt"
file_list <- paste0("interim_data/ldsc/plink_clumped_1KGP_HG19/", diet_trait, "/",phewas_trait,"/clumped_results_",phewas_trait,"_chr", 1:22,".tsv.clumped")
all_data1 <- lapply(file_list, function(f) read.table(f, header = TRUE)) %>% bind_rows()
all_data1$pheno <- phewas_trait

phewas_trait <- "30000_irnt"
file_list <- paste0("interim_data/ldsc/plink_clumped_1KGP_HG19/oilyfsih_TG/",phewas_trait,"/clumped_results_",phewas_trait,"_chr", 1:22,".tsv.clumped") #had a typo in last script here
all_data2 <- lapply(file_list, function(f) read.table(f, header = TRUE)) %>% bind_rows()
all_data2$pheno <- phewas_trait

phewas_trait <- "30050_irnt"
file_list <- paste0("interim_data/ldsc/plink_clumped_1KGP_HG19/", diet_trait, "/",phewas_trait,"/clumped_results_",phewas_trait,"_chr", 1:22,".tsv.clumped")
all_data3 <- lapply(file_list, function(f) read.table(f, header = TRUE)) %>% bind_rows()
all_data3$pheno <- phewas_trait

# Combine all data
conf_snp <- rbind(all_data1, all_data2, all_data3)
dim(conf_snp) 

# Save
write.table(conf_snp, paste0("interim_data/MVMR/oilyfish_tg/exp_clumped_snps_",diet_trait,"_MVMR_GIclump.txt"), row.names = FALSE, quote = FALSE)

#~~~~~~~~~~~~~~~~~~~
# Get diet SNPs and add to other exposures for clumping
# Oilyfish
# Load Exposure 2 top hits
diet_sig <- readRDS("interim_data/exposure/formatted_oilyfish_clumped_noLD_nobiallelic.rds")
# Subset to SNP/variant and P
diet_sig <- diet_sig %>% select(variant, pval.exposure)
colnames(diet_sig) <- c("SNP", "P")

# Subset conf_snp to SNP and P
conf_snp_sub <- conf_snp %>% select(SNP, P)
# Rbind
mvmr_snps <- rbind(conf_snp_sub, diet_sig)

# Save
#write.table(mvmr_snps, paste0("interim_data/MVMR/oilyfish_tg/all_clumped_snps_",diet_trait,"_MVMR_GIclump.txt"), row.names = FALSE, quote = FALSE)
write.table(mvmr_snps, paste0("interim_data/MVMR/oilyfish_tg/all_clumped_snps_oilyfish_tg_MVMR_GIclump.txt"), row.names = FALSE, quote = FALSE)
# Custom here due to name discrepancy erg

#~~~~~~~~~~~~~~~~~~~
# any repeats?
length(unique(conf_snp$SNP)) #4445
unq_conf_snp <- unique(conf_snp$SNP)
unq_conf_snp <- as.data.frame(unq_conf_snp)
names(unq_conf_snp) <- "variant"

# Save
write.table(unq_conf_snp, paste0("interim_data/MVMR/oilyfish_tg/unique_snps_",diet_trait,"_MVMR_GIclump.txt"), row.names = FALSE, quote = FALSE) #doesn't have diet


#~~~~~~~~~~~~~~~~~~~
# NEED TO CLUMP THESE
# Source: https://mrcieu.github.io/TwoSampleMR/articles/perform_mr.html#multivariable-mr
# Which states,
# Get instruments for each exposure
# Combine these into a set of all instruments
# Clump these to avoid the possibility that e.g. a variant for exposure 1 is in LD with a variant for exposure 2
# Re-extract all the final clumped SNPs from all of the exposures
# Harmonise them all to be on the same effect allele
# Use the multivariable MR method against these harmonised data

# Go run MVMR_plink_clump_oilyfish_TG_clumpGI.sh
# Come back, load results
phewas_trait <- "GI_clumped"
diet_trait <- "oilyfish_tg"
file_list <- paste0("interim_data/MVMR/", diet_trait, "/GI_clumped/",phewas_trait,"_results_",diet_trait,"_chr", 1:22,".tsv.clumped")
clumped_data <- lapply(file_list, function(f) read.table(f, header = TRUE)) %>% bind_rows()
length(unique(clumped_data$SNP)) # all unique, 1528
# This is what we will use for MVMR!
# Now we need to re-extract all the final clumped SNPs from the exposures

# Save
write.table(clumped_data, paste0("interim_data/MVMR/oilyfish_tg/clumped_",diet_trait,"_MVMR_GIsnps.txt"), row.names = FALSE, quote = FALSE)

# Load if needed
clumped_data <- read.table(paste0("interim_data/MVMR/oilyfish_tg/clumped_",diet_trait,"_MVMR_GIsnps.txt"), header = TRUE)

snp_list <- clumped_data$SNP
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Prep for MVMR - developing pipeline
# https://wspiller.github.io/MVMR/articles/MVMR.html
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Strategy:
# - Gather all significant SNPs from the clumped files
# - Subset each full GWAS summary stat file by those SNPs
# - Column bind to make the GI

#~~~~~~~~~~~~~~~~~~~~
# Subset
# - Load GWAS sum-stat file, subset, save subset, unload GWAS sum-stat file
# - Repeat with next file

# 30110_irnt
# MANUALLY DEFINE TRAIT
trait <- "30110_irnt"
# Full summary stats
df_full <- read.table(paste0("/pl/active/colelab/users/kjames/refinedMR/interim_data/ldsc/NL_bothsexes_temp/",trait,".gwas.imputed_v3.both_sexes.tsv"), header = TRUE)
# Extract SNPs from GWAS
df_subset <- df_full[df_full$variant %in% snp_list,]
# Save subsetted SNPs in case R runs out of memory and closes
saveRDS(df_subset, paste0("interim_data/MVMR/",diet_trait,"/df_",trait,"_subset.rds"))
# Remove full summary stats to clear space
rm(df_full)

# 30000_irnt
# MANUALLY DEFINE TRAIT
trait <- "30000_irnt"
# Full summary stats
df_full <- read.table(paste0("/pl/active/colelab/users/kjames/refinedMR/interim_data/ldsc/NL_bothsexes_temp/",trait,".gwas.imputed_v3.both_sexes.tsv"), header = TRUE)
# Extract SNPs from GWAS
df_subset <- df_full[df_full$variant %in% snp_list,]
# Save subsetted SNPs in case R runs out of memory and closes
saveRDS(df_subset, paste0("interim_data/MVMR/",diet_trait,"/df_",trait,"_subset.rds"))
# Remove full summary stats to clear space
rm(df_full)

# 30050_irnt
# MANUALLY DEFINE TRAIT
trait <- "30050_irnt"
# Full summary stats
df_full <- read.table(paste0("/pl/active/colelab/users/kjames/refinedMR/interim_data/ldsc/NL_bothsexes_temp/",trait,".gwas.imputed_v3.both_sexes.tsv"), header = TRUE)
# Extract SNPs from GWAS
df_subset <- df_full[df_full$variant %in% snp_list,]
# Save subsetted SNPs in case R runs out of memory and closes
saveRDS(df_subset, paste0("interim_data/MVMR/",diet_trait,"/df_",trait,"_subset.rds"))
# Remove full summary stats to clear space
rm(df_full)

# Diet full sum stats
# Load full exposure 2 data 
trait <- "oilyfish"
diet_full <- readRDS("interim_data/exposure/formatted_oilyfish_clumped_full_sumstats.rds")
# Extract SNPs from GWAS
df_diet_subset <- diet_full[diet_full$variant %in% snp_list,]
# Save subsetted SNPs in case R runs out of memory and closes
saveRDS(df_diet_subset, paste0("interim_data/MVMR/",diet_trait,"/df_",trait,"_subset.rds"))
# Remove full summary stats to clear space
rm(diet_full)

# Outcome - NEEDS SPECIAL ATTENTION
# Now add in TG (outcome)
trait <- "TG"
out <- readRDS("interim_data/outcome/formatted_GCST90239664_TG_Graham_GRCh37.rds") 
# Extract SNPs from GWAS
df_out_subset <- out[out$variant %in% snp_list,] # gave 0, investigate

# Investigating
snp_list_noA0A1 <- stringr::str_replace(snp_list, "^([^:]+:[0-9]+):.*", "\\1")
out$chr_pos <- stringr::str_replace(out$variant, "^([^:]+:[0-9]+):.*", "\\1")
df_out_subset <- out[out$chr_pos %in% snp_list_noA0A1,] #2309, pretty close

# Save subsetted SNPs in case R runs out of memory and closes
saveRDS(df_out_subset, paste0("interim_data/MVMR/",diet_trait,"/df_",trait,"_subset.rds"))
# Remove full summary stats to clear space
rm(out)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
diet_trait2 <- "oilyfish_tg"
# Read in subsets
trait <- "30110_irnt"
df_subset1 <- readRDS(paste0("interim_data/MVMR/",diet_trait2,"/df_",trait,"_subset.rds"))
colnames(df_subset1) <- paste0(colnames(df_subset1), ".",trait)
df_subset1 <- df_subset1 %>% rename(variant = variant.30110_irnt) #update
df_subset1$chr_pos <- stringr::str_replace(df_subset1$variant, "^([^:]+:[0-9]+):.*", "\\1")

trait <- "30000_irnt"
df_subset2 <- readRDS(paste0("interim_data/MVMR/",diet_trait2,"/df_",trait,"_subset.rds"))
colnames(df_subset2) <- paste0(colnames(df_subset2), ".",trait)
df_subset2 <- df_subset2 %>% rename(variant = variant.30000_irnt) #update
df_subset2$chr_pos <- stringr::str_replace(df_subset2$variant, "^([^:]+:[0-9]+):.*", "\\1")

trait <- "30050_irnt"
df_subset3 <- readRDS(paste0("interim_data/MVMR/",diet_trait2,"/df_",trait,"_subset.rds"))
colnames(df_subset3) <- paste0(colnames(df_subset3), ".",trait)
df_subset3 <- df_subset3 %>% rename(variant = variant.30050_irnt) #update
df_subset3$chr_pos <- stringr::str_replace(df_subset3$variant, "^([^:]+:[0-9]+):.*", "\\1")

trait <- "oilyfish"
df_subset6 <- readRDS(paste0("interim_data/MVMR/",diet_trait2,"/df_",trait,"_subset.rds"))
colnames(df_subset6) <- paste0(colnames(df_subset6), ".",trait)
df_subset6 <- df_subset6 %>% rename(variant = variant.oilyfish) #update
df_subset6$chr_pos <- stringr::str_replace(df_subset6$variant, "^([^:]+:[0-9]+):.*", "\\1")

trait <- "TG"
df_subset7 <- readRDS(paste0("interim_data/MVMR/",diet_trait2,"/df_",trait,"_subset.rds"))
colnames(df_subset7) <- paste0(colnames(df_subset7), ".",trait)
df_subset7 <- df_subset7 %>% rename(variant = variant.TG) #update
df_subset7$chr_pos <- stringr::str_replace(df_subset7$variant, "^([^:]+:[0-9]+):.*", "\\1")

# Combine
dfs <- list(df_subset1, df_subset2, df_subset3, df_subset6, df_subset7)
# Use left_join (prioritizes the snps in the first position)
mvmr_dat <- purrr::reduce(dfs, left_join, by = "chr_pos") #"variant")

# Complete only
mvmr_dat <- na.omit(mvmr_dat) #1039


# Try MVMR
# Note, will need to get covariance of traits from UKB RAP
# Also, need to check harmonization etc.

# quickly subset to cols with beta
beta <- mvmr_dat %>% select(starts_with("beta."))

# Load correlation matrix from the RAP
cor_mat <- read.csv("interim_data/MVMR/cor_res_mvmr_pipe_all_250618.csv", row.names = 1)
colnames(cor_mat) # order the F.data to match this 

# subset cormat to relevant traits
keep_exp <- c("p30110_i0", 
              "p30000_i0",
              "p30050_i0",
              "p1558_i0")
cor_mat <- cor_mat[keep_exp, keep_exp]

# Run MVMR
library(MVMR)
F.data <- format_mvmr(BXGs = cbind(mvmr_dat$beta.30110_irnt,
                                   mvmr_dat$beta.30000_irnt,
                                   mvmr_dat$beta.30050_irnt,
                                   mvmr_dat$beta.exposure.oilyfish), #beta of exposures
                      BYG = mvmr_dat$beta.outcome.TG, # beta of outcome
                      seBXGs = cbind(mvmr_dat$se.30110_irnt,
                                     mvmr_dat$se.30000_irnt,
                                     mvmr_dat$se.30050_irnt,
                                     mvmr_dat$se.exposure.oilyfish), # se of exposures
                      seBYG = mvmr_dat$se.outcome.TG, # se of outcome
                      RSID = mvmr_dat$variant.x)
head(F.data)
str(mvmr_dat)
seBXGs <- matrix(cbind(mvmr_dat$se.30110_irnt,
                       mvmr_dat$se.30000_irnt,
                       mvmr_dat$se.30050_irnt,
                       mvmr_dat$se.exposure.oilyfish),
                 byrow = TRUE,
                 ncol=4)
gencov <- phenocov_mvmr(cor_mat, seBXGs)
str(gencov)

# Need to turn gencov into a list with matrices. This will make the matrix multiplication work.
gencov <- lapply(gencov, as.matrix)

# Test for weak instruments
sres <- strength_mvmr(r_input = F.data, gencov = gencov)
sres
# Save
write.csv(sres, paste0("interim_data/MVMR/oilyfish_tg/",diet_trait,"_MVMR_F_GIstrength.txt"))

# Pleiotropy
pres <- pleiotropy_mvmr(r_input = F.data, gencov = gencov) 
pres
# Save
write.csv(pres, paste0("interim_data/MVMR/oilyfish_tg/",diet_trait,"_MVMR_Q_GIheterogeneity.txt"))

# Estimate causality
res <- ivw_mvmr(r_input = F.data, gencov = gencov)
res

# Get 95% CI - oily fish
low_95CI <- res[4,1] - (res[4,2]*1.96) #alcohol is the 6th exposure
high_95CI <- res[4,1] + (res[4,2]*1.96)
low_95CI
high_95CI

# weight
exp1_low_95CI <- res[1,1] - (res[1,2]*1.96) #alcohol is the 6th exposure
exp1_high_95CI <- res[1,1] + (res[1,2]*1.96)
exp1_low_95CI
exp1_high_95CI

# sitting height
exp2_low_95CI <- res[2,1] - (res[2,2]*1.96) #alcohol is the 6th exposure
exp2_high_95CI <- res[2,1] + (res[2,2]*1.96)
exp2_low_95CI
exp2_high_95CI

# Save
write.csv(res, paste0("interim_data/MVMR/oilyfish_tg/",diet_trait2,"_MVMR_IVW.txt"))




