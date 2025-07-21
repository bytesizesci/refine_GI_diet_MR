# Goal: Harmonize of exposure and outcome data
# Exposure: alcohol
# Outcome: liver cirrhosis
# Author: KJS
# Date: Sep 3, 2024

# Load libraries
library(dplyr)
library(LDlinkR)

# Set wd
wdir <- "/pl/active/colelab/users/kjames/enviroMR/"
setwd(wdir)

# Load Exposure
ex <- readRDS("interim_data/exposure/formatted_alcohol_clumped_noLD_nobiallelic.rds")

# Load Outcome
out <- readRDS("interim_data/outcome/formatted_GCST90319877_cir_Ghouse.rds") 

# Check if the SNP column exists
if ("SNP" %in% colnames(out)) {
  # SNP column is present
  print("SNP column present")
  # Show the first few entries of the SNP column
  print(head(out$SNP))
} else {
  # SNP column is missing
  print("SNP column missing, using chr:pos data")
}

# Manual processing due to uncommon structur in out data
out$effect_allele.outcome <- toupper(out$effect_allele.outcome)
out$other_allele.outcome <- toupper(out$other_allele.outcome)
# Make chr_pos merge column
ex$merge_ID <- paste0(ex$chr.exposure,":", ex$pos.exposure)
out$merge_ID <- paste0(out$chr.outcome,":", out$pos.outcome)

# Merge to make genetic instrument (chr:bp:other:effect)
gi = merge(ex, out, by = "merge_ID", 
           all.x = TRUE, all.y = FALSE)

# Check strand ambiguity
# Use custom function; Flags A/T and C/G SNPs
source("scripts/is_ambiguous.R")
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
dim(gi_complete) #61

# Harmonize complete data
source("scripts/harmonize_alleles.R")
gi_harmon <- harmonize_alleles(gi_complete)
gi_harmon <- gi_harmon$data

# Any changes?
identical(gi_harmon, gi_complete) # FALSE - harmonization was needed

# Get list of missing outcome SNPs
gi_needproxy <- gi[!complete.cases(gi), ]
dim(gi_needproxy) # 29

##############
# Most effective way to get proxies is using chr#:### format
#proxySearch <- gi_needproxy$variant
#proxySearch <- paste0("chr", proxySearch)
# Get rid of allele parts
#parts <- strsplit(proxySearch, ":")
#chrpos_string <- c()
#for(i in 1:length(proxySearch)){
#  chrpos_string[[i]] <- paste0(parts[[i]][1], ":", parts[[i]][2])
#}
#chrpos_string <- do.call(rbind, chrpos_string)

# Get chr#:## from mergeID variable created above
chrpos_string <- paste0("chr", gi_needproxy$merge_ID)

# Proxy search
noLD_snps <- list()
missing_snps <- list()
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
  
  # From proxies, find which have data in out
  proxies_in_outgw2 <- out_proxies[out_proxies$RS_Number %in% out$SNP,]
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
  # Make final table
  # Want original SNP, proxy SNP, proxy SNP data
  out_proxy_temp$original_SNP <- snp_temp
  out_proxy_temp$original_variant <- variantID_temp
  out_proxy_temp <- out_proxy_temp %>% select(original_SNP, original_variant, SNP, everything()) # rearrange columns
  out_proxy_temp <- out_proxy_temp %>% rename("proxySNP" = "SNP")
  # save data 
  out_proxy_results[[i]] <- out_proxy_temp
}

# Turn list into data.frame
proxy_res <- do.call("rbind", out_proxy_results) 
noLD_snps_res <- do.call("rbind", noLD_snps) 
missing_snps_res <- do.call("rbind", missing_snps) 

# Nothing we can do about the snps with no good LD proxy
cat("SNPs in the exposure that were missing in the outcome but don't have a proxy in LD:", noLD_snps_res[,1])

# We could try to manually look up proxies for missing_snps. Its not many. Move on.
cat("SNPs in the exposure that were missing in the outcome and were missing in 1000G and/or dbSNP155:", missing_snps_res[,1])

# Merge the missing_snps_res snps using rsid to see if you get any matches
ex_test <- ex[paste0("chr", ex$merge_ID) %in% missing_snps_res,]
gi_rsid <- merge(ex_test, out, by = "SNP")
gi_rsid <- harmonize_alleles(gi_rsid)
gi_rsid <- gi_rsid$data
gi_rsid <- gi_rsid %>% rename(#"SNP" = "SNP.x",
                              "merge_ID"="merge_ID.x") %>%
          select(-c(merge_ID.y))

gi_test = gtools::smartbind(gi_harmon, gi_rsid)
gi_test <- gi_test %>% 
  dplyr::mutate(SNP.x = coalesce(SNP, SNP.x)) %>%
  select(-c(SNP))

gi_harmon <- gi_test

# ~~~~~~~~~~~~~~~
# Merge with exposure-outcome data
# ~~~~~~~~~~~~~~~
# First...
# Recall the SNPs that needed a proxy
dim(gi_needproxy) # 29x38
# Get rid of NA cols in the original "full" data
gi_needproxy <- gi_needproxy[,1:19]
# Create a merge ID col
parts <- strsplit(gi_needproxy$variant, ":")
chrpos_string <- c()
for(i in 1:nrow(gi_needproxy)){
  chrpos_string[[i]] <- paste0(parts[[i]][1], ":", parts[[i]][2])
}
chrpos_string <- do.call(rbind, chrpos_string)
chrpos_string <- paste0("chr", chrpos_string)
gi_needproxy <- cbind(gi_needproxy, chrpos_string)

# Merge proxies with the exposure half of the data
gi_proxies <- merge(gi_needproxy, proxy_res, by.x = "chrpos_string", by.y = "original_SNP", all.x = FALSE, all.y = TRUE)

# Note, harmonizing the proxy data using effect/other alleles will not work as 
# the proxies have their own genotype.

# Add clarifying names/column to proxy data
# newname=oldname
gi_proxies = gi_proxies %>% rename(
  "chrpos.original"="chrpos_string",
  "variant.original"="variant.x",
  "rsID.original" = "SNP.x",
  "variant.outcome" = "variant") # corresponds to proxy in this case
  #"proxy_rsID" = "rs_id")

# Note, the outcome data in exp_out_proxy_res corresponds to the proxy! So the positions, alleles, rsID etc. will be different
gi_proxies$proxy <- "proxy"

# Then combine proxies and originals 
gi2 = gtools::smartbind(gi_harmon, gi_proxies)
gi2$proxy = gi2$proxy %>% tidyr::replace_na("original")

# Check 
table(gi2$proxy) # as expected!

# Note:
# SNP.x = original exposure SNP; originated from exposure data
# SNP.y = outcome SNP; originated from outcome data
# rsID = rsID of the original or proxy SNP; originated from the outcome data
# original_variant = variant ID of the original exposure variant - only for variants that needed proxies; originated from exposure data
# [OLD]variant.x = variant ID of the original exposure variants - only for the variants that needed proxies; originated from exposure data
# [OLD]variant.y = variant ID for proxies; originated from outcome data
# variant = variant ID of the original exposure variants - only for the variants with matches
# [OLD]Note, variant.x and variant columns can be combined to create the column variant.exposure
# Note, original_variant and variant columns can be combined to create the column variant.exposure

# Merge variant.x and variant
# Note, for select -c(), some variables are not in a given df but are in others
colnames(gi2)
gi3 <- gi2 %>% 
  dplyr::mutate(variant.exposure = coalesce(variant.x, variant.original)) %>%
  dplyr::mutate(snp.exposure = coalesce(SNP.x, rsID.original)) %>%
  dplyr::mutate(snp.outcome = coalesce(SNP.y, proxySNP)) %>%
  dplyr::select(-c(variant.x, variant.y, variant.original, SNP.x, rsID.original, SNP.y, proxySNP,
                   merge_ID.x, merge_ID.y, merge_ID, ambiguous, remove, exp_out_af_flag, chrpos.original,
                   variant.original, rsID.original, original_variant, variant.outcome)) %>% # proxy_rsID, rs_id)) %>%
  dplyr::select(variant.exposure, snp.exposure, snp.outcome, proxy,
         everything())

# Retest 
gi4 <- harmonize_alleles(gi3)
gi5 <- gi4[1] 
gi5 <- do.call(rbind, gi5)

# Relabel proxies (some were mislabeled due to flipping that needed to be performed)
gi5$proxy <- ifelse(gi5$snp.outcome == gi5$snp.exposure, "original", "proxy")
table(gi5$proxy)

colnames(gi5)

# Save
saveRDS(gi5, "interim_data/merged_GI/mrdat_alcohol_GCST90319877_cir_Ghouse.rds")
write.csv(gi5, "interim_data/merged_GI/mrdat_alcohol_GCST90319877_cir_Ghouse.csv")

# Out of an abundance of curiosity
MRobj_std <- MendelianRandomization:: mr_input(bx = gi5$beta.exposure,
                                               bxse = gi5$se.exposure,
                                               by = gi5$beta.outcome,
                                               byse = gi5$se.outcome)
MRres_std <- MendelianRandomization::mr_allmethods(MRobj_std, method = "all")
MRres_std
