# Goal: Harmonize of exposure and outcome data
# Exposure: whiteXwholegrain
# Outcome: CVD - GRCh37
# Author: KJS
# Date: Sep 3, 2024

# Load libraries
library(dplyr)
library(LDlinkR)

# Set wd
wdir <- "/pl/active/colelab/users/kjames/refinedMR/"
setwd(wdir)

# Load Exposure
ex <- readRDS("interim_data/exposure/formatted_whiteXwholegrain_clumped_noLD_nobiallelic.rds")

# Load full exposure data for when we need proxies
ex_full <- readRDS("interim_data/exposure/formatted_whiteXwholegrain_clumped_full_sumstats.rds")

# Load Outcome
out <- readRDS("interim_data/outcome/formatted_GCST90132314_CVD_Aragam_GRCh37.rds") 

# Manual processing due to uncommon structure in out data
out$effect_allele.outcome <- toupper(out$effect_allele.outcome)
out$other_allele.outcome <- toupper(out$other_allele.outcome)

# Make chr_pos merge column
ex$merge_ID <- paste0(ex$chr.exposure,":", ex$pos.exposure)
out$merge_ID <- paste0(out$chr.outcome,":", out$pos.outcome)

# Merge to make genetic instrument
gi = merge(ex, out, by = "merge_ID", # "SNP"
           all.x = TRUE, all.y = FALSE)

# Merge to make genetic instrument (chr:bp:other:effect)
#gi = merge(ex, out, by = "SNP", #"variant", 
#           all.x = TRUE, all.y = FALSE)

# Check strand ambiguity
# Use custom function; Flags A/T and C/G SNPs
source("scripts/harmonization/is_ambiguous.R")
# Identify ambiguous SNPs
gi$ambiguous <- ifelse(is_ambiguous(gi$effect_allele.exposure, gi$other_allele.exposure) == TRUE, "ambiguous", "ok")

# Check allele frequency 
# If 0.45 < allele frequency < 0.55 then we will say we can't confidently say which strand is which and we will remove them
gi$remove <- ifelse(gi$ambiguous == "ambiguous" & gi$eaf.exposure <= 0.55 & gi$eaf.exposure >= 0.45, "remove", "ok")
table(gi$remove)

# Apply filter
gi <- gi[gi$remove == "ok",] # only keep "ok" not "remove"

# Also, check that the exposure and outcome allele frequencies are within +-0.05 of each other. 
# They should be (especially if input is from the same ancestry)
# Add a new column to flag discrepancies
gi$exp_out_af_flag <- ""
# Define the threshold for flagging discrepancies
threshold <- 0.05
# Function to flag discrepancies
gi$exp_out_af_flag <- ifelse(abs(gi$eaf.exposure - gi$eaf.outcome) > threshold, "Flag", "ok")
table(gi$exp_out_af_flag) 
# Filter, if needed. All OK here.
# If needed, repeat on just complete data

# Get list of complete/matching SNPs
gi_complete = gi[complete.cases(gi), ]
dim(gi_complete) #42

# Harmonize complete data
source("scripts/harmonization/harmonize_alleles.R")
gi_harmon <- harmonize_alleles(gi_complete)
gi_harmon <- gi_harmon$data

# Any changes?
identical(gi_harmon, gi_complete) # FALSE - that means we had some switched alleles

# Get list of missing outcome SNPs
gi_needproxy <- gi[!complete.cases(gi), ]
dim(gi_needproxy) # 2

##############
# In this data (oily fish and TG) there is a variant.x and variant.y column
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

# ~~~~~~~~~~~~~~~
# Merge with exposure-outcome data
# ~~~~~~~~~~~~~~~
# For BREAD TYPE and CVD there were no proxies - SKIP to SAVING
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
saveRDS(gi_harmon, "interim_data/merged_GI/mrdat_exp_out_proxies_whiteXwholegrain_GCST90132314_CVD_Aragam.rds")
write.csv(gi_harmon, "interim_data/merged_GI/mrdat_exp_out_proxies_whiteXwholegrain_GCST90132314_CVD_Aragam.csv")

# Out of an abundance of curiosity
MRobj_std <- MendelianRandomization:: mr_input(bx = gi_harmon$beta.exposure,
                                               bxse = gi_harmon$se.exposure,
                                               by = gi_harmon$beta.outcome,
                                               byse = gi_harmon$se.outcome)
MRres_std <- MendelianRandomization::mr_allmethods(MRobj_std, method = "all")
MRres_std
