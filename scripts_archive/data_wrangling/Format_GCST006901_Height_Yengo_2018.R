# Goal: 2018 Height Yengo et al data - to be used as a negative control
# Author: KJS
# Date: November 1, 2024

# https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST006001-GCST007000/GCST006901/README_summary_statistics_Yengo_et_al_2018.txt
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
wdir <- "/pl/active/colelab/users/kjames/enviroMR/"
setwd(wdir)

#~~~~~~~~~~~~~~~
# 1) 2018 DATA
#~~~~~~~~~~~~~~~
# Set trait 
trait <- "Height_Yengo_2018"

file_path <- "/pl/active/colelab/users/kjames/enviroMR/inputs/GCST006901_Height_Yengo_2018/Meta-analysis_Wood_et_al+UKBiobank_2018.txt.gz"
file_path_unzipped <- "/pl/active/colelab/users/kjames/enviroMR/inputs/GCST006901_Height_Yengo_2018/Meta-analysis_Wood_et_al+UKBiobank_2018.txt"

# Unzip file if needed
if (grepl("\\.gz$", file_path, ignore.case = TRUE)) {
  # File has a .gz extension, so unzip it
  gunzip(file_path, remove = FALSE)  # `remove = FALSE` keeps the original file
  message("File has been gunzipped.")
} else {
  message("File is not a gzip file; no action taken.")
}

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
         "beta.outcome"="BETA",
         "pval.outcome"="P",
         "chr.outcome"="CHR",
         "se.outcome"="SE",
         "pos.outcome"="POS",
         "effect_allele.outcome" = "Tested_Allele",
         "other_allele.outcome" = "Other_Allele",
         "eaf.outcome" = "Freq_Tested_Allele_in_HRS",
         "samplesize.outcome" = "N"
  )

# Add manual columns
# Columns are needed for twosampleMR harmonization
df2$outcome <- trait
df2$mr_keep.outcome = TRUE
df2$pval_origin.outcome = "reported"
df2$data_source.outcome = "GWAS_catalog"

# Add allele frequency columns 
# Column is needed for twosampleMR harmonization
df2$oaf.outcome = (1 - df2$eaf.outcome) # effect allele frequency for outcome

# Add sample size column if it is missing
# Numbers obtained from yaml file
#df2$samplesize.outcome <- 437438

# NOTE: JBC's effect allele and beta are opposite common formatting! 
# Switch her data so you don't have to switch every outcome you work with
#df2$beta.outcome <- df2$beta.outcome * -1
#df2 <- df2 %>% rename("effect_allele.outcome" = "other_allele.outcome",
#                      "other_allele.outcome" = "effect_allele.outcome", 
#                      "oaf.outcome" = "eaf.outcome",
#                      "eaf.outcome" = "oaf.outcome",
#                      "allelecount_other.outcome" = "allelecount_effect.outcome",
#                      "allelecount_effect.outcome" = "allelecount_other.outcome")

# Make a variant column
# We will match the format of Neale's variant variable - chr:pos:ref:alt (chr:pos:other:effect)
# where ref=aligned to the forward strand
# alt=the effect allele
df2$variant <- paste0(df2$chr.outcome, ":", df2$pos.outcome, ":", df2$other_allele.outcome,":", df2$effect_allele.outcome)

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
                      chr.outcome,
                      pos.outcome,
                      effect_allele.outcome,
                      other_allele.outcome,
                      eaf.outcome,
                      oaf.outcome,
                      samplesize.outcome,
                      beta.outcome,
                      se.outcome,
                      pval.outcome,
                      #id.outcome,
                      mr_keep.outcome,
                      outcome,
                      pval_origin.outcome,
                      data_source.outcome
                      #chromosome_grch38,
                      #start_grch38,
                      #rs_id 
                      #h2g.observed,
                      #h2g.P.value
)

# Save
write.csv(df3, "interim_data/outcome/formatted_GCST006901_Height_Yengo_2018.csv")
saveRDS(df3, "interim_data/outcome/formatted_GCST006901_Height_Yengo_2018.rds")

