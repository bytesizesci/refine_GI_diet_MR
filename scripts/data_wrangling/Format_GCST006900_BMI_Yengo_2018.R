# Goal: 2018 BMI Yengo et al data - to be used in MVMR
# Author: KJS
# Date: December 30, 2024

# Data download: https://portals.broadinstitute.org/collaboration/giant/index.php/GIANT_consortium_data_files#BMI_and_Height_GIANT_and_UK_BioBank_Meta-analysis_Summary_Statistics
# "Download Updated Meta-analysis Locke et al + UKBiobank 2018 GZIP"

# README: https://portals.broadinstitute.org/collaboration/giant/images/0/01/README_summary_statistics_Yengo_et_al_2018.txt
## This describes the columns of the summary statistics generated in Yengo et al. (2018)
## Meta-analysis of genome-wide association studies for height and body mass index in ~700,000 individuals of European ancestry

#----------------------
#  Columns description
#----------------------
#SNP:                        RS ID
#CHR:                        Chromosome
#POS:                        Physical position (Genome build hg19)
#Tested_Allele:              Allele corresponding to the effect size (BETA/BETA_COJO)
#Other_Allele:               Other allele
#Freq_Tested_Allele_in_HRS:  Frequency of the tested allele in the Health and Retirement Study (from 8,552 unrelated participants).
#BETA:                       Marginal SNP effect size.
#SE:                         Standard error of the effect size.
#P:                          P-value measuring the significance of the marginal effect.
#N:                          Sample size.
#BETA_COJO:                  Conditional effect size from the software GCTA version 1.9 using HRS as linkage disequilibrium panel.
#SE_COJO:                    Standard error of the conditional effect size.
#P_COJO:                     P-value measuring the significance of the conditional effect.

# Load necessary library
library(dplyr)
library(R.utils)

# Set wd
wdir <- "/pl/active/colelab/users/kjames/refinedMR/"
setwd(wdir)

#~~~~~~~~~~~~~~~
# 1) 2018 DATA - full
#~~~~~~~~~~~~~~~
# Set trait 
trait <- "BMI_Yengo_2018"

file_path_unzipped <- "/pl/active/colelab/users/kjames/refinedMR/raw_data/BMI_GIANT_Yengo/Meta-analysis_Locke_et_al+UKBiobank_2018_UPDATED.txt"

#test path for developing code interactively
df <- read.table(file_path_unzipped, header = TRUE)
head(df)

# Get colnames for renaming below
colnames(df)

# Rename columns - following TwoSampleMR format where possible
# Some columns (like eaf.exposure we need to add manually)
# Height has RSID, which we will keep calling RSID
# new_name = old_name
df2 = df %>% 
  rename("SNP" = "SNP",
         "beta.exposure"="BETA",
         "pval.exposure"="P",
         "chr.exposure"="CHR",
         "se.exposure"="SE",
         "pos.exposure"="POS",
         "effect_allele.exposure" = "Tested_Allele",
         "other_allele.exposure" = "Other_Allele",
         "eaf.exposure" = "Freq_Tested_Allele_in_HRS",
         "samplesize.exposure" = "N"
  )

# Add manual columns
# Columns are needed for twosampleMR harmonization
df2$exposure <- trait
df2$mr_keep.exposure = TRUE
df2$pval_origin.exposure = "reported"
df2$data_source.exposure = "GWAS_catalog"

# Add allele frequency columns 
# Column is needed for twosampleMR harmonization
df2$oaf.exposure = (1 - df2$eaf.exposure) # effect allele frequency for exposure

# Add sample size column if it is missing
# Numbers obtained from yaml file
#df2$samplesize.exposure <- 437438

# NOTE: JBC's effect allele and beta are opposite common formatting! 
# Switch her data so you don't have to switch every exposure you work with
#df2$beta.exposure <- df2$beta.exposure * -1
#df2 <- df2 %>% rename("effect_allele.exposure" = "other_allele.exposure",
#                      "other_allele.exposure" = "effect_allele.exposure", 
#                      "oaf.exposure" = "eaf.exposure",
#                      "eaf.exposure" = "oaf.exposure",
#                      "allelecount_other.exposure" = "allelecount_effect.exposure",
#                      "allelecount_effect.exposure" = "allelecount_other.exposure")

# Make a variant column
# We will match the format of Neale's variant variable - chr:pos:ref:alt (chr:pos:other:effect)
# where ref=aligned to the forward strand
# alt=the effect allele
df2$variant <- paste0(df2$chr.exposure, ":", df2$pos.exposure, ":", df2$other_allele.exposure,":", df2$effect_allele.exposure)

# Note: identical(df2$SNP, df2$variant) is false. The reference/alternate is different than effect/other at times. 
# As long as we are consistent in how we name things and how we check which strand the beta corresponds to we are okay.

# Get rid of extra columns
colnames(df2)
#df2 <- df2 %>% select(-variant_id, -odds_ratio, -ci_lower, -ci_upper)

# View
head(df2)

# Reorganize columns
df3 <- df2 %>% select(SNP, 
                      variant, #chr:pos:other:effect
                      #RSID,
                      chr.exposure,
                      pos.exposure,
                      effect_allele.exposure,
                      other_allele.exposure,
                      eaf.exposure,
                      oaf.exposure,
                      samplesize.exposure,
                      beta.exposure,
                      se.exposure,
                      pval.exposure,
                      #id.exposure,
                      mr_keep.exposure,
                      exposure,
                      pval_origin.exposure,
                      data_source.exposure
                      #chromosome_grch38,
                      #start_grch38,
                      #rs_id 
                      #h2g.observed,
                      #h2g.P.value
)

# Save
write.csv(df3, "interim_data/exposure/formatted_GCST006900_BMI_Yengo_2018.csv")
saveRDS(df3, "interim_data/exposure/formatted_GCST006900_BMI_Yengo_2018.rds")

#~~~~~~~~~~~~~~~
# 1) 2018 DATA
#~~~~~~~~~~~~~~~
# Set trait 
trait <- "BMI_Yengo_2018_top"

file_path_unzipped <- "/pl/active/colelab/users/kjames/refinedMR/raw_data/BMI_GIANT_Yengo/Meta-analysis_Locke_et_al+UKBiobank_2018_top_941_from_COJO_analysis_UPDATED.txt"

#test path for developing code interactively
df <- read.table(file_path_unzipped, header = TRUE)
head(df)

# Get colnames for renaming below
colnames(df)

# Rename columns - following TwoSampleMR format where possible
# Some columns (like eaf.exposure we need to add manually)
# Height has RSID, which we will keep calling RSID
# new_name = old_name
df2 = df %>% 
  rename("SNP" = "SNP",
         "beta.exposure"="BETA",
         "pval.exposure"="P",
         "chr.exposure"="CHR",
         "se.exposure"="SE",
         "pos.exposure"="POS",
         "effect_allele.exposure" = "Tested_Allele",
         "other_allele.exposure" = "Other_Allele",
         "eaf.exposure" = "Freq_Tested_Allele_in_HRS",
         "samplesize.exposure" = "N"
  )

# Add manual columns
# Columns are needed for twosampleMR harmonization
df2$exposure <- trait
df2$mr_keep.exposure = TRUE
df2$pval_origin.exposure = "reported"
df2$data_source.exposure = "GWAS_catalog"

# Add allele frequency columns 
# Column is needed for twosampleMR harmonization
df2$oaf.exposure = (1 - df2$eaf.exposure) # effect allele frequency for exposure

# Add sample size column if it is missing
# Numbers obtained from yaml file
#df2$samplesize.exposure <- 437438

# NOTE: JBC's effect allele and beta are opposite common formatting! 
# Switch her data so you don't have to switch every exposure you work with
#df2$beta.exposure <- df2$beta.exposure * -1
#df2 <- df2 %>% rename("effect_allele.exposure" = "other_allele.exposure",
#                      "other_allele.exposure" = "effect_allele.exposure", 
#                      "oaf.exposure" = "eaf.exposure",
#                      "eaf.exposure" = "oaf.exposure",
#                      "allelecount_other.exposure" = "allelecount_effect.exposure",
#                      "allelecount_effect.exposure" = "allelecount_other.exposure")

# Make a variant column
# We will match the format of Neale's variant variable - chr:pos:ref:alt (chr:pos:other:effect)
# where ref=aligned to the forward strand
# alt=the effect allele
df2$variant <- paste0(df2$chr.exposure, ":", df2$pos.exposure, ":", df2$other_allele.exposure,":", df2$effect_allele.exposure)

# Note: identical(df2$SNP, df2$variant) is false. The reference/alternate is different than effect/other at times. 
# As long as we are consistent in how we name things and how we check which strand the beta corresponds to we are okay.

# Get rid of extra columns
colnames(df2)
#df2 <- df2 %>% select(-variant_id, -odds_ratio, -ci_lower, -ci_upper)

# View
head(df2)

# Reorganize columns
df3 <- df2 %>% select(SNP, 
                      variant, #chr:pos:other:effect
                      #RSID,
                      chr.exposure,
                      pos.exposure,
                      effect_allele.exposure,
                      other_allele.exposure,
                      eaf.exposure,
                      oaf.exposure,
                      samplesize.exposure,
                      beta.exposure,
                      se.exposure,
                      pval.exposure,
                      #id.exposure,
                      mr_keep.exposure,
                      exposure,
                      pval_origin.exposure,
                      data_source.exposure
                      #chromosome_grch38,
                      #start_grch38,
                      #rs_id 
                      #h2g.observed,
                      #h2g.P.value
)

# Save
write.csv(df3, "interim_data/exposure/formatted_top_GCST006900_BMI_Yengo_2018.csv")
saveRDS(df3, "interim_data/exposure/formatted_top_GCST006900_BMI_Yengo_2018.rds")
