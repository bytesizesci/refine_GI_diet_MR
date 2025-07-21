# Goal: Format JBC diet data - oily fish
# Author: KJS
# Date: September 1st, 2024

# Set wd
wdir <- "/pl/active/colelab/users/kjames/enviroMR/"
setwd(wdir)

# Load libraries
library(LDlinkR)
library(dplyr)
library(stringr)

# Read data
df <- read.csv("/pl/active/colelab/users/kjames/enviroMR/interim_data/exposure/oilyfish_clumped_noLD_nobiallelic.csv", row.names = 1)

# Add A1 and A0 column using strsplit() command
# We are assuming A1=effect, A0=other
col_split=df$CHR_BP_ALLELE1_ALLELE0

# Wrangle the chr:pos:Allele:Allele format 
# Make a new A1 and A0 column
parts <- strsplit(col_split, ":")
# Remove the last element from the vector
A1_string <- c()
A0_string <- c()
for(i in 1:length(col_split)){
  A1_string[[i]] <- parts[[i]][3]
  A0_string[[i]] <- parts[[i]][4]
}
A1_string <- do.call(rbind, A1_string)
A0_string <- do.call(rbind, A0_string)

# Combine to df
df$effect_allele.exposure <- A1_string
df$other_allele.exposure <- A0_string

# Get colnames for renaming below
colnames(df)

# Rename columns - following TwoSampleMR format where possible
# Some columns (like eaf.exposure we need to add manually)
# new_name = old_name
df2 = df %>% 
  rename("SNP" = "SNP",
         "beta.exposure"="BETA",
         "pval.exposure"="P",
         "chr.exposure"="CHR",
         "se.exposure"="SE",
         "pos.exposure"="BP",
         "oaf.exposure" = "A0FREQ",
         "exposure"="TRAIT",
         "allelecount_other.exposure" = "AC_A1",
         "allelecount_effect.exposure" = "AC_A0"
  )

# Add manual columns
# Columns are needed for twosampleMR harmonization
df2$mr_keep.exposure = TRUE
df2$pval_origin.exposure = "reported"
df2$data_source.exposure = "JBC_GWAS_clumped"

# Add allele frequency columns 
# Column is needed for twosampleMR harmonization
df2$eaf.exposure = (1 - df2$oaf.exposure) # effect allele frequency for exposure

# Get N
df2$samplesize.exposure <- round((df2$allelecount_effect.exposure+df2$allelecount_other.exposure)/2)

# NOTE: JBC's effect allele and beta are opposite common formatting! 
# Switch her data so you don't have to switch every outcome you work with
df2$beta.exposure <- df2$beta.exposure * -1
df2 <- df2 %>% rename("effect_allele.exposure" = "other_allele.exposure",
                      "other_allele.exposure" = "effect_allele.exposure", 
                      "oaf.exposure" = "eaf.exposure",
                      "eaf.exposure" = "oaf.exposure",
                      "allelecount_other.exposure" = "allelecount_effect.exposure",
                      "allelecount_effect.exposure" = "allelecount_other.exposure")

# Make a variant column
# We will match the format of Neale's variant variable - chr:pos:ref:alt (chr:pos:other:effect)
# where ref=aligned to the forward strand
# alt=the effect allele
df2$variant <- paste0(df2$chr.exposure, ":", df2$pos.exposure, ":", df2$other_allele.exposure,":", df2$effect_allele.exposure)

# Get rid of CHR_BP_ALLELE1_ALLELE0 col
df2 <- df2 %>% select(-CHR_BP_ALLELE1_ALLELE0)

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
                      data_source.exposure,
                      #h2g.observed,
                      #h2g.P.value,
                      merge_ID,
                      allelecount_effect.exposure,
                      allelecount_other.exposure)

# Save
write.csv(df3, "/pl/active/colelab/users/kjames/enviroMR/interim_data/exposure/formatted_oilyfish_clumped_noLD_nobiallelic.csv")
saveRDS(df3, "/pl/active/colelab/users/kjames/enviroMR/interim_data/exposure/formatted_oilyfish_clumped_noLD_nobiallelic.rds")
