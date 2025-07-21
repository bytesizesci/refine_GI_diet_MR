# Goal: MVMR - Harmonize exposures and outcome 
# Exposure 1: alcohol
# Exposure 2: BMI
# Outcome: alanine aminotransferase, ALT
# Author: KJS
# Date: December 30, 2024

# Load libraries
library(dplyr)
library(LDlinkR)

# Set wd
wdir <- "/pl/active/colelab/users/kjames/refinedMR/"
setwd(wdir)

# Load Exposure 1
ex1 <- readRDS("interim_data/exposure/formatted_alcohol_clumped_noLD_nobiallelic.rds")

# Load full exposure 1 data for when we need proxies
ex1_full <- readRDS("interim_data/exposure/formatted_alcohol_clumped_full_sumstats.rds")

# Load Exposure 2 top hits
ex2 <- readRDS("interim_data/exposure/formatted_top_GCST006900_BMI_Yengo_2018.rds")

# Load full exposure 2 data 
ex2_full <- readRDS("interim_data/exposure/formatted_alcohol_clumped_noLD_nobiallelic.rds")

# Load Outcome
out <- readRDS("interim_data/outcome/formatted_GCST90013405_ALT_Pazoki.rds") 

#~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Merge exposure 1 and 2
#~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Make chr_pos merge column
ex1$merge_ID <- paste0(ex1$chr.exposure,":", ex1$pos.exposure)
ex2$merge_ID <- paste0(ex2$chr.exposure,":", ex2$pos.exposure)

ex1_full$merge_ID <- paste0(ex1_full$chr.exposure,":", ex1_full$pos.exposure)
ex2_full$merge_ID <- paste0(ex2_full$chr.exposure,":", ex2_full$pos.exposure)

# Merge to make genetic instrument
#gi_full = merge(ex1_full, ex2_full, by = "merge_ID", # "SNP"
           #all.x = TRUE, all.y = TRUE)

#gi = gi_full %>% filter(merge_ID %in% ex1$merge_ID | merge_ID %in% ex2$merge_ID)

# Merge to make genetic instrument
gi_alc = merge(ex2_full, ex1, by = "merge_ID", # "SNP"
           all.x = FALSE, all.y = TRUE)

gi_bmi = merge(ex2, ex1_full, by = "merge_ID", # "SNP"
           all.x = TRUE, all.y = FALSE) 

# Check column names of gi_alc and gi_bmi. Simplify, then rbind.
colnames(gi_alc)
colnames(gi_bmi)

# Keep essential cols
keepcols <- c("merge_ID", 
              # exposure 1
              "SNP.x", "variant.x", "chr.exposure.x", "pos.exposure.x", "effect_allele.exposure.x", "other_allele.exposure.x",
              "eaf.exposure.x", "oaf.exposure.x", "beta.exposure.x", "se.exposure.x", "pval.exposure.x", "exposure.x",
              # exposure 2
              "SNP.y", "variant.y", "chr.exposure.y", "pos.exposure.y", "effect_allele.exposure.y", "other_allele.exposure.y",
              "eaf.exposure.y", "oaf.exposure.y", "beta.exposure.y", "se.exposure.y", "pval.exposure.y", "exposure.y")

gi_alc <- gi_alc[,keepcols,drop=FALSE]
gi_bmi <- gi_bmi[,keepcols,drop=FALSE]

# Rbind
gi_exp <- rbind(gi_alc, gi_bmi)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Harmonize exposure data
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Check strand ambiguity
# Use custom function; Flags A/T and C/G SNPs
source("scripts/harmonization/is_ambiguous.R")
# Identify ambiguous SNPs
gi_exp$ambiguous <- ifelse(is_ambiguous(gi_exp$effect_allele.exposure.x, gi_exp$other_allele.exposure.x) == TRUE, "ambiguous", "ok")

# Check allele frequency 
# If 0.45 < allele frequency < 0.55 then we will say we can't confidently say which strand is which and we will remove them
gi_exp$remove <- ifelse(gi_exp$ambiguous == "ambiguous" & gi_exp$eaf.exposure.x <= 0.55 & gi_exp$eaf.exposure.x >= 0.45, "remove", "ok")
table(gi_exp$remove)

# Apply filter
gi_exp <- gi_exp[gi_exp$remove == "ok",] # only keep "ok" not "remove", 21 removed leaving 1012

# Also, check that the exposure and outcome allele frequencies are within +-0.05 of each other. 
# They should be (especially if input is from the same ancestry)
# Add a new column to flag discrepancies
gi_exp$exp_out_af_flag <- ""
# Define the threshold for flaggi_expng discrepancies
threshold <- 0.05
# Function to flag discrepancies
gi_exp$exp_out_af_flag <- ifelse(abs(gi_exp$eaf.exposure.x - gi_exp$eaf.exposure.y) > threshold, "Flag", "ok")
table(gi_exp$exp_out_af_flag) 
# Flag is a sign that the data needs to be harmonized

# Get list of complete/matching SNPs
gi_exp_complete = gi_exp[complete.cases(gi_exp), ]
dim(gi_exp_complete) #1012

# Harmonize complete data
# Use custom script that accounts for MVMR column names
source("scripts/harmonization/harmonize_alleles.R")
source("scripts/harmonization/harmonize_alleles_mvmr_nameAssist.R")
gi_exp_harmon = rename_and_execute(gi_exp_complete, harmonize_alleles)

# Any changes?
identical(gi_exp_harmon, gi_exp_complete) # FALSE - that means we had some switched alleles

# Get list of missing outcome SNPs
gi_exp_needproxy <- gi_exp[!complete.cases(gi_exp), ]
dim(gi_exp_needproxy) # 0, don't need to find proxies; use gi_exp_harmon

# Rename gi_exp_harmon to make it simpler to use
# genetic instrument for exposures
giE <- gi_exp_harmon

# Save
#saveRDS(giE, "interim_data/merged_GI/MVMR_alcohol_BMI_GI_241230.rds")

# Load if needed
giE <- readRDS("interim_data/merged_GI/MVMR_alcohol_BMI_GI_241230.rds")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Merge exposures with outcome
#~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Make chr_pos merge column for outcome
out$merge_ID <- paste0(out$chr.outcome,":", out$pos.outcome)

# Merge to make genetic instrument
giEO = merge(giE, out, by = "merge_ID", # "SNP"
           all.x = TRUE, all.y = FALSE)

# Check strand ambiguity
# Use custom function; Flags A/T and C/G SNPs
source("scripts/harmonization/is_ambiguous.R")
# Identify ambiguous SNPs
giEO$ambiguous <- ifelse(is_ambiguous(giEO$effect_allele.outcome, giEO$other_allele.outcome) == TRUE, "ambiguous", "ok")

# Check allele frequency - on outcome data (the new data to the GI)
# If 0.45 < allele frequency < 0.55 then we will say we can't confidently say which strand is which and we will remove them
giEO$remove <- ifelse(giEO$ambiguous == "ambiguous" & giEO$eaf.outcome <= 0.55 & giEO$eaf.outcome >= 0.45, "remove", "ok")
table(giEO$remove)

# Apply filter
giEO <- giEO[giEO$remove == "ok",] # only keep "ok" not "remove"

# Also, check that the exposure and outcome allele frequencies are within +-0.05 of each other. 
# They should be (especially if input is from the same ancestry)
# Add a new column to flag discrepancies
giEO$exp_out_af_flag <- ""
# Define the threshold for flaggiEOng discrepancies
threshold <- 0.05
# Function to flag discrepancies
giEO$exp_out_af_flag <- ifelse(abs(giEO$eaf.exposure.x - giEO$eaf.outcome) > threshold, "Flag", "ok")
table(giEO$exp_out_af_flag) 
# Flag is a sign that the data needs to be harmonized

# Get list of complete/matching SNPs
giEO_complete = giEO[complete.cases(giEO), ]
dim(giEO_complete) #987 complete, 44 missing

# Harmonize complete data
#source("scripts/harmonization/harmonize_alleles.R")
#giEO_harmon <- harmonize_alleles(giEO_complete)
#giEO_harmon <- giEO_harmon$data

# Use custom script that accounts for MVMR column names
source("scripts/harmonization/harmonize_alleles.R")
source("scripts/harmonization/harmonize_alleles_mvmr_nameAssist_step2.R")
giEO_harmon = rename_and_execute_step2(giEO_complete, harmonize_alleles)

# Any changes?
identical(giEO_harmon, giEO_complete) # FALSE - that means we had some switched alleles

# Save
#saveRDS(giEO_harmon, "interim_data/merged_GI/MVMR_alcohol_BMI_ALT_GI_241231.rds")
# Load if needed
giEO_harmon <- readRDS("interim_data/merged_GI/MVMR_alcohol_BMI_ALT_GI_241231.rds")


# Get list of missing outcome SNPs
giEO_needproxy <- giEO[!complete.cases(giEO), ] #CHECK THIS CHECK THIS CHECK THIS
dim(giEO_needproxy) # 24x44

# Try MVMR package (installed from GitHub)
library(MVMR)
F.data <- format_mvmr(BXGs = cbind(giEO_harmon$beta.exposure.x, giEO_harmon$beta.exposure.y),
                      BYG = giEO_harmon$beta.outcome,
                      seBXGs = cbind(giEO_harmon$se.exposure.x, giEO_harmon$se.exposure.y),
                      seBYG = giEO_harmon$se.outcome,
                      RSID = giEO_harmon$merge_ID)

pcov_mat <- matrix(c(1,-72,-72,1),nrow=2,ncol=2) #-72 obtained from UKBRAP cov()
seBXGs <- matrix(cbind(giEO_harmon$se.exposure.x, giEO_harmon$se.exposure.y),ncol=2)
gencov <- phenocov_mvmr(pcov_mat, seBXGs)
# instrument strength
sres <- strength_mvmr(r_input = F.data, gencov = gencov) # really weak... probably stop...

# Test heterogeneity
pres <- pleiotropy_mvmr(r_input=F.data, gencov = gencov)
pres # heterogeneity not significant

# MVMR
res <- ivw_mvmr(r_input = F.data, gencov = gencov) # we get the same estimates as the MendelianRandomization package below

# Robust causal effect estimation - when you have weak instruments
res1 <- qhet_mvmr(F.data, pcov_mat, CI = TRUE, iterations = 100) #CI=FALSE - runs w/o CI; fails with CI
res1

#
MVMRinputobj <- MendelianRandomization::mr_mvinput(bx = cbind(giEO_harmon$beta.exposure.x, giEO_harmon$beta.exposure.y),
                                                   bxse = cbind(giEO_harmon$se.exposure.x, giEO_harmon$se.exposure.y),
                                                   by = giEO_harmon$beta.outcome,
                                                   byse = giEO_harmon$se.outcome)

MVIVWObject <- MendelianRandomization::mr_mvivw(MVMRinputobj,
                       model = "default",
                       correl = FALSE,
                       correl.x = NULL,
                       nx = NA,
                       distribution = "normal",
                       alpha = 0.05)

MVIVWObject




















##############
# In this data (alcohol and ALT) there is a variant.x and variant.y column
# Some of the SNPs are not harmonized which causes identical(variant.x, variant.y)==FALSE
# We are harmonizing again at the end so its okay... but use variant.x (exposure based) here

# Most effective way to get proxies is using chr#:### format
proxySearch <- gi_needproxy$variant.x #variant 
proxySearch <- paste0("chr", proxySearch)
# Get rid of allele parts
parts <- strsplit(proxySearch, ":")
chrpos_string <- c()
for(i in 1:length(proxySearch)){
  chrpos_string[[i]] <- paste0(parts[[i]][1], ":", parts[[i]][2])
}
chrpos_string <- do.call(rbind, chrpos_string)

# Proxy search
noLD_snps <- list()
missing_snps <- list()
no_combinedProxy <- list()
out_proxy_results <- list() # store results
for(i in 1:length(chrpos_string)){
  # Get snp ID for proxy finding
  snp_temp = chrpos_string[i]
  print(snp_temp)
  # For variant id merging later
  variantID_temp = chrpos_string[i]
  print(variantID_temp)
  # Get proxy
  out_proxies = LDproxy(snp = snp_temp,
                        pop = "GBR", #CAN SPECIFY THE REFERENCE PANEL HERE
                        r2d = "r2", 
                        genome_build = "grch37", #default
                        token = "1125e7a68581") #ADD USER'S TOKEN
  
  # If the SNP is missing, report it and skip
  if(any(grepl("error", out_proxies))){
    print(paste0("WARNING: SNP not present in 1000 genomes, must search for proxy with another method: ", snp_temp))
    missing_snps[[i]] <- snp_temp
    next
  }
  print("Working on proxy info")
  
  # From proxies, find which have data in outcome and exposure data
  proxies_in_outgw2 <- out_proxies[out_proxies$RS_Number %in% out$SNP & out_proxies$RS_Number %in% ex_full$SNP,]
  
  # If there isn't a proxy with data in exposure and outcome, skip it
  if(nrow(proxies_in_outgw2) == 0){
    no_combinedProxy[[i]] <- snp_temp
    next # Skip to the next iteration of the loop if object has 0 rows
  }
  
  # From proxies in out, take SNP (row of data) with highest r2 value
  proxy_temp <- proxies_in_outgw2[proxies_in_outgw2$R2 == max(proxies_in_outgw2$R2),]
  # Sometimes multiple SNPs have the highest r2 value
  # Also filter by closest distance to original SNP
  abs_df <- abs(proxy_temp$Distance)
  row_index <- which.min(abs_df)
  proxy_temp <- proxy_temp[row_index,]
  proxy_temp_name <- proxy_temp$RS_Number
  
  # If the r2 is <0.8, skip it
  cat("r2 value is", proxy_temp$R2)
  if(proxy_temp$R2 < 0.8){
    print(paste0("WARNING: r2 value is less than 0.8 and proxy should be skipped. r2 =", proxy_temp$R2))
    noLD_snps[[i]] <- snp_temp
    next
  }
  
  # Now that we have a proxy SNP selected, get its outcome data
  out_proxy_temp <- out %>% filter(SNP %in% proxy_temp_name)
  
  # Get the same SNP in the exposure data
  exp_proxy_temp <- ex_full %>% filter(SNP %in% proxy_temp_name)
  
  # Merge
  proxy_df <- merge(exp_proxy_temp, out_proxy_temp, by = "SNP")
  
  # Add a new column to indicate a proxy
  proxy_df$proxy <- "proxy"
  
  # Add a new column for the original SNP
  proxy_df$original_exposure_SNP <- variantID_temp
  
  # save data 
  out_proxy_results[[i]] <- proxy_df
}

# Turn list into data.frame
proxy_res <- do.call("rbind", out_proxy_results) 
noLD_snps_res <- do.call("rbind", noLD_snps) 
missing_snps_res <- do.call("rbind", missing_snps) 
no_ex_out_snps_res <- do.call("rbind", no_combinedProxy)

# Nothing we can do about the snps with no good LD proxy
cat("SNPs in the exposure that were missing in the outcome but don't have a proxy in LD:", noLD_snps_res[,1])

# Nothing we can do about the snps with no good proxy in exp AND out data sets
cat("SNPs whose outcome proxy was not in the exposure data:", no_ex_out_snps_res[,1])

# We could try to manually look up proxies for missing_snps. Its not many. Move on.
cat("SNPs in the exposure that were missing in the outcome and were missing in 1000G and/or dbSNP155:", missing_snps_res[,1])

# Save proxy result items
save(list= c("proxy_res", "noLD_snps_res", "missing_snps_res", "no_ex_out_snps_res"),
     file = "interim_data/merged_GI/proxydata_exp_out_proxies_alcohol_GCST90013405_ALT_Pazoki.Rdata") 

# ~~~~~~~~~~~~~~~
# Merge with exposure-outcome data
# ~~~~~~~~~~~~~~~
# For ALCOHOL and ALT there were no proxies - SKIP to SAVING
# Then combine proxies and originals 
#gi2 = gtools::smartbind(gi_harmon, proxy_res)

# In the proxy column, replace NA with original_SNP 
#gi2$proxy = gi2$proxy %>% tidyr::replace_na("original_SNP")

# Re-harmonize
#gi3 <- harmonize_alleles(gi2)
#gi3 <- gi3$data

# Get rid of variant.x and variant.y columns - we will rely on variant column correctly created in harmonize_alleles function
#gi3 <- gi3 %>% select(-variant.x, -variant.y, -merge_ID)
gi_harmon <- gi_harmon %>% select(-variant.x, -variant.y, -merge_ID)

# Save
saveRDS(gi_harmon, "interim_data/merged_GI/mrdat_exp_out_proxies_alcohol_GCST90013405_ALT_Pazoki.rds")
write.csv(gi_harmon, "interim_data/merged_GI/mrdat_exp_out_proxies_alcohol_GCST90013405_ALT_Pazoki.csv")

# Out of an abundance of curiosity
MRobj_std <- MendelianRandomization:: mr_input(bx = gi_harmon$beta.exposure,
                                               bxse = gi_harmon$se.exposure,
                                               by = gi_harmon$beta.outcome,
                                               byse = gi_harmon$se.outcome)
MRres_std <- MendelianRandomization::mr_allmethods(MRobj_std, method = "all")
MRres_std