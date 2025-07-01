# Goal: GCST90239658_LDL_Graham data
# Author: KJS
# Date: September 3, 2024

# Load necessary library
library(dplyr)
library(R.utils)

# Set wd
wdir <- "/pl/active/colelab/users/kjames/enviroMR/"
setwd(wdir)

# Set trait 
trait <- "GCST90239658"

#file_path <- "/pl/active/colelab/users/kjames/enviroMR/inputs/GCST90239658_LDL_Graham/GCST90239658_buildGRCh37.tsv"
file_path_unzipped <- "/pl/active/colelab/users/kjames/enviroMR/inputs/GCST90239658_LDL_Graham/GCST90239658_buildGRCh37.tsv"

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
#table(df$strand) # should all be +

# Get colnames for renaming below
colnames(df)

# Rename columns - following TwoSampleMR format where possible
# Some columns (like eaf.exposure we need to add manually)
# new_name = old_name
df2 = df %>% 
  rename("SNP" = "variant_id", #"rsid",
         "beta.outcome"="beta",
         "pval.outcome"="p_value",
         "chr.outcome"="chromosome",
         "se.outcome"="standard_error",
         "pos.outcome"="base_pair_location",
         "effect_allele.outcome" = "effect_allele",
         "other_allele.outcome" = "other_allele",
         "eaf.outcome" = "effect_allele_frequency",
         "samplesize.outcome" = "n"
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

# Add sample size column
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

# Get rid of extra columns
colnames(df2)
df2 <- df2 %>% select(-N_studies)

# View
head(df2)

# Reorganize columns
df3 <- df2 %>% select(SNP,
                      variant,
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
write.csv(df3, "interim_data/outcome/formatted_GCST90239658_LDL_Graham.csv")
saveRDS(df3, "interim_data/outcome/formatted_GCST90239658_LDL_Graham.rds")