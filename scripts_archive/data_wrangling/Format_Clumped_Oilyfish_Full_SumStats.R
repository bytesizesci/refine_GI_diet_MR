# Goal: Format JBC full diet data - oily fish
# Author: KJS
# Date: November 11, 2024

# Load libraries
library(LDlinkR)
library(dplyr)
library(stringr)

# Load full exposure data for when we need proxies
ex_full <- read.table("/pl/active/colelab/common/published_gwas/diet/cole_naturecomm_2020/BOLTlmm_UKB_genoQCEURN455146_v3_diet_oilyfish_overallfreq.1329.average_QT_BgenSnps_mac20_maf0.005_info0.6.gz", header = TRUE)

# Add A1 and A0 column using strsplit() command
# A1=effect
# A0=other

# Rename columns - following TwoSampleMR format where possible
# Some columns (like eaf.exposure we need to add manually)
# new_name = old_name
df2 = ex_full %>% 
  rename("SNP" = "SNP",
         "beta.exposure"="BETA",
         "pval.exposure"="P_BOLT_LMM",
         "chr.exposure"="CHR",
         "se.exposure"="SE",
         "pos.exposure"="BP",
         "eaf.exposure" = "A1FREQ",
         #"oaf.exposure" = "A0FREQ",
         #"exposure"="TRAIT",
         "other_allele.exposure" = "ALLELE0",
         "effect_allele.exposure" = "ALLELE1" # we switch them below
  )

# Add manual columns
# Columns are needed for twosampleMR harmonization
df2$mr_keep.exposure = TRUE
df2$pval_origin.exposure = "reported"
df2$data_source.exposure = "JBC_GWAS_full"
df2$exposure <- "oilyfish_fullGWAS"

# Add allele frequency columns 
# Column is needed for twosampleMR harmonization
#df2$eaf.exposure = (1 - df2$oaf.exposure) # effect allele frequency for exposure
df2$oaf.exposure = (1 - df2$eaf.exposure) # other allele frequency for exposure

# Get N
#df2$samplesize.exposure <- round((df2$allelecount_effect.exposure+df2$allelecount_other.exposure)/2)
df2$samplesize.exposure <- 446854

# NOTE: JBC's effect allele and beta are opposite common formatting! 
# Switch her data so you don't have to switch every outcome you work with
df2$beta.exposure <- df2$beta.exposure * -1
df2 <- df2 %>% rename("effect_allele.exposure" = "other_allele.exposure",
                      "other_allele.exposure" = "effect_allele.exposure", 
                      "oaf.exposure" = "eaf.exposure",
                      "eaf.exposure" = "oaf.exposure")
                      #"allelecount_other.exposure" = "allelecount_effect.exposure",
                      #"allelecount_effect.exposure" = "allelecount_other.exposure")

# Make a variant column
# We will match the format of Neale's variant variable - chr:pos:ref:alt (chr:pos:other:effect)
# where ref=aligned to the forward strand
# alt=the effect allele
df2$variant <- paste0(df2$chr.exposure, ":", df2$pos.exposure, ":", df2$other_allele.exposure,":", df2$effect_allele.exposure)

# Get rid of CHR_BP_ALLELE1_ALLELE0 col or other unwanted columns depending on input format
df2 <- df2 %>% select(-INFO) #-CHR_BP_ALLELE1_ALLELE0)

# Reorganize columns
df3 <- df2 %>% select(SNP,
                      variant,
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
                      data_source.exposure)
                      #h2g.observed,
                      #h2g.P.value,
                      #merge_ID,
                      #allelecount_effect.exposure,
                      #allelecount_other.exposure)

# Save
write.csv(df3, "/pl/active/colelab/users/kjames/refinedMR/interim_data/exposure/formatted_oilyfish_clumped_full_sumstats.csv")
saveRDS(df3, "/pl/active/colelab/users/kjames/refinedMR/interim_data/exposure/formatted_oilyfish_clumped_full_sumstats.rds")
