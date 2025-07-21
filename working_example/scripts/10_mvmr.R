# Step 10: 10_mvmr.R
# Goal: Run mvmr! Prior to running MVMR, need to get correlation of traits and
# need to check harmonization of the outcome data with one exposure

# About the traits in the working example:
#23098 - exposure 1 - weight - body composition by impedance
#20015 - exposure 2 - sitting height - body size measures
#23106 - exposure 3 - impedance of whole body - body composition by impedance
#30080 - exposure 4 - platelet count - blood count
#30270 - exposure 5 - mean sphered cell volume - blood count

# Load libraries
library(MVMR)
library(dplyr)
library(tibble)
# ===================== USER CONFIGURATION p1 =====================
# Set working direcotry
#setwd("/pl/active/colelab/users/kjames/refinedMR/")
setwd("~/working_example")

# Set traits
EXP_pheno <- "alcohol"
OUT_pheno <- "ALT"
diet_trait <- paste0(EXP_pheno, "_", OUT_pheno)

# Set output file directory
output_dir <- paste0("~/working_example/results/GI/", diet_trait)

# Load data
mvmr_dat <- readRDS(output_dir, "/MVMR_GI_input_", diet_trait, ".rds")

# Identified in 04_snp_overlap.R
exposure_traits <- c("23098_irnt",
                     "20015_irnt",
                     "23106_irnt",
                     "30080_irnt",
                     "30270_irnt",
                     "1558") #alcohol intake frequency, look this up from GWAS round 2 results, manifest, https://www.nealelab.is/uk-biobank
# ==============================================================

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Step 1: Generate the correlation matrix
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
df <- read.table("~/working_example/data/gcor_PheWASttest_results_all.txt")

# Header
gcor_header <- c("p1", "p2", "rg", "se", "z", "p", "h2_obs", "h2_obs_se", "h2_int", "h2_int_se", "gcov_int", "gcov_int_se")

# Assign header to df
colnames(df) <- gcor_header

# Make sure p1-p2 and p2-p1 exists (dcast has issues if not)
df_full <- df %>%
  mutate(p1_alt = p2, p2_alt = p1) %>%
  rename(rg_alt = rg) %>%
  select(p1 = p1_alt, p2 = p2_alt, rg = rg_alt) %>%
  bind_rows(df)  # combine original and flipped

# Now cast
wide_matrix <- reshape2::dcast(df_full, p1 ~ p2, value.var = "rg")

# Set rownames
rownames(wide_matrix) <- wide_matrix[,"p1"]
wide_matrix <- wide_matrix %>% select(-p1)
str(wide_matrix)

# Add diagonal
diag(wide_matrix) <- 1

# Add ".ldsc.imputed_v3.both_sexes.tsv.bgz" to exposure_traits
phewas_traits_long <- paste0(exposure_traits, ".ldsc.imputed_v3.both_sexes.tsv.bgz")

# Now only keep rows and columns that are in the trait list
# Filter columns
wide_matrix_sub <- wide_matrix %>%
  select(all_of(phewas_traits_long)) 

# Filter rows; have to turn the df into a tibble
wide_matrix_sub <- wide_matrix_sub %>%
  rownames_to_column(var = "trait")
wide_matrix_sub <- filter(wide_matrix_sub, trait %in% phewas_traits_long)

# Sort
desired_row_order <- colnames(wide_matrix_sub)[-1]
wide_matrix_sub <- wide_matrix_sub %>%
  mutate(trait = factor(trait, levels = desired_row_order)) %>%
  arrange(trait)

# Turn back into a square matrix
wide_matrix_sub <- wide_matrix_sub %>%
  column_to_rownames("trait") %>%
  as.matrix()

# Check
str(wide_matrix_sub)  

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Step 2: Format the data for MVMR 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ===================== USER CONFIGURATION p2 =====================
# Adjust the MVMR columns to match your data!

# Run MVMR
# When creating the mvmr object, make sure the order matches the order in wide_matrix_sub
F.data <- format_mvmr(BXGs = cbind(mvmr_dat$beta.23098_irnt,
                                   mvmr_dat$beta.20015_irnt,
                                   mvmr_dat$beta.23106_irnt,
                                   mvmr_dat$beta.30080_irnt,
                                   mvmr_dat$beta.30270_irnt,
                                   mvmr_dat$beta.exposure.alcohol), #beta of exposures
                      BYG = mvmr_dat$beta.outcome.ALT, # beta of outcome
                      seBXGs = cbind(mvmr_dat$se.23098_irnt,
                                     mvmr_dat$se.20015_irnt,
                                     mvmr_dat$se.23106_irnt,
                                     mvmr_dat$se.30080_irnt,
                                     mvmr_dat$se.30270_irnt,
                                     mvmr_dat$se.exposure.alcohol), # se of exposures
                      seBYG = mvmr_dat$se.outcome.ALT, # se of outcome
                      RSID = mvmr_dat$variant.x)
head(F.data)
str(mvmr_dat)

# Adjust the MVMR se columns to match your data!
# Make se matrix
seBXGs <- matrix(cbind(mvmr_dat$se.23098_irnt,
                       mvmr_dat$se.20015_irnt,
                       mvmr_dat$se.23106_irnt,
                       mvmr_dat$se.30080_irnt,
                       mvmr_dat$se.30270_irnt,
                       mvmr_dat$se.exposure.alcohol),
                 byrow = TRUE,
                 ncol=6)

# make genetic covariance object
gencov <- phenocov_mvmr(wide_matrix_sub, seBXGs)
str(gencov)

# Need to turn gencov into a list with matrices. This will make the matrix multiplication work.
gencov <- lapply(gencov, as.matrix)
# ==============================================================

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Step 1: MVMR & Sensitivity 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Test for weak instruments
sres <- strength_mvmr(r_input = F.data, gencov = gencov) 
sres
# Save res
write.csv(sres, paste0("~/working_example/results/tables/strength_mvmr_",diet_trait,"_f_results.txt"))

# Pleiotropy
pres <- pleiotropy_mvmr(r_input = F.data, gencov = gencov)
pres
# Save pres
write.csv(pres, paste0("~/working_example/results/tables/strength_mvmr_",diet_trait,"_p_results.txt"))

# Estimate causality
res <- ivw_mvmr(r_input = F.data, gencov = gencov)
res
# Save res
write.csv(res, paste0("~/working_example/results/tables/MVMR_",diet_trait,"_results.txt"))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Step 1: Record 95% CI 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Get 95% CI - alcohol
low_95CI <- res[6,1] - (res[6,2]*1.96) #alcohol is the 6th exposure
cat("Low 95% confidence interval:", low_95CI)

high_95CI <- res[6,1] + (res[6,2]*1.96)
cat("High 95% confidence interval:", high_95CI)

# weight
exp1_low_95CI <- res[1,1] - (res[1,2]*1.96) 
cat("Low 95% confidence interval:", exp1_low_95CI)

exp1_high_95CI <- res[1,1] + (res[1,2]*1.96)
cat("High 95% confidence interval:", exp1_high_95CI)

# sitting height
exp2_low_95CI <- res[2,1] - (res[2,2]*1.96) 
cat("Low 95% confidence interval:", exp2_low_95CI)

exp2_high_95CI <- res[2,1] + (res[2,2]*1.96)
cat("High 95% confidence interval:", exp2_high_95CI)

