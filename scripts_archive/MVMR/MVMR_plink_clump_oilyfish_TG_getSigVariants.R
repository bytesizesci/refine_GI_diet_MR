

library(dplyr)

#Set working direcotry
setwd("/pl/active/colelab/users/kjames/refinedMR/")


# Load the clumped results in 
#tst <- read.table("interim_data/ldsc/plink_clumped_1KGP_HG19/oilyfish_tg/clumped_results_4124_irnt_chr1.tsv.clumped", header = TRUE)

# Load programatically
file_list <- paste0("interim_data/ldsc/plink_clumped_1KGP_HG19/oilyfish_tg/4124_irnt/clumped_results_4124_irnt_chr", 1:22,".tsv.clumped")
all_data1 <- lapply(file_list, function(f) read.table(f, header = TRUE)) %>% bind_rows()

file_list <- paste0("interim_data/ldsc/plink_clumped_1KGP_HG19/oilyfish_tg/30100_irnt/clumped_results_30100_irnt_chr", 1:22,".tsv.clumped")
all_data2 <- lapply(file_list, function(f) read.table(f, header = TRUE)) %>% bind_rows()

file_list <- paste0("interim_data/ldsc/plink_clumped_1KGP_HG19/oilyfish_tg/4100_irnt/clumped_results_4100_irnt_chr", c(1:20, 22),".tsv.clumped") # selected in kmeans 2, ankle spacing width
all_data3 <- lapply(file_list, function(f) read.table(f, header = TRUE)) %>% bind_rows()

file_list <- paste0("interim_data/ldsc/plink_clumped_1KGP_HG19/oilyfish_tg/50_irnt/clumped_results_50_irnt_chr", 1:22,".tsv.clumped") # selected in kmeans 2, standing height
all_data4 <- lapply(file_list, function(f) read.table(f, header = TRUE)) %>% bind_rows()

file_list <- paste0("interim_data/ldsc/plink_clumped_1KGP_HG19/oilyfish_tg/5084_irnt/clumped_results_5084_irnt_chr", 1:22,".tsv.clumped")
all_data5 <- lapply(file_list, function(f) read.table(f, header = TRUE)) %>% bind_rows()

file_list <- paste0("interim_data/ldsc/plink_clumped_1KGP_HG19/oilyfish_tg/23106_irnt/clumped_results_23106_irnt_chr", 1:22,".tsv.clumped")
all_data6 <- lapply(file_list, function(f) read.table(f, header = TRUE)) %>% bind_rows()

file_list <- paste0("interim_data/ldsc/plink_clumped_1KGP_HG19/oilyfish_tg/23104_irnt/clumped_results_23104_irnt_chr", 1:22,".tsv.clumped")
all_data7 <- lapply(file_list, function(f) read.table(f, header = TRUE)) %>% bind_rows()

file_list <- paste0("interim_data/ldsc/plink_clumped_1KGP_HG19/oilyfish_tg/23101_irnt/clumped_results_23101_irnt_chr", 1:22,".tsv.clumped")
all_data8 <- lapply(file_list, function(f) read.table(f, header = TRUE)) %>% bind_rows()

# So, this is an active development... If we use the AIC criteria it recommends 2 clusters with 4100 and 50 rising to the top
# Use this to keep developing the method
all_data3$pheno <- "4100_irnt"
all_data4$pheno <- "50_irnt"

# combine
conf_snp <- rbind(all_data3, all_data4)
dim(conf_snp) #4333

# any repeats?
length(unique(conf_snp$SNP)) #4315 --> so some but few repeats

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Prep for MVMR - developing pipeline
# https://wspiller.github.io/MVMR/articles/MVMR.html
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 3/27/25 - 
# Strategy:
# - Gather all significant SNPs 
# - Subset each full GWAS summary stat file by those SNPs
# - Column bind to make the GI
# - Recall, focusing on traits 4100 and 50 because they were identified from the AIC k-means strategy

#~~~~~~~~~~~~~~~~~~~~
# 4100_irnt
# Clumped, significant loci
file_list <- paste0("interim_data/ldsc/plink_clumped_1KGP_HG19/oilyfish_tg/4100_irnt/clumped_results_4100_irnt_chr", c(1:20, 22),".tsv.clumped") # selected in kmeans 2, ankle spacing width
df_4100 <- lapply(file_list, function(f) read.table(f, header = TRUE)) %>% bind_rows()

#~~~~~~~~~~~~~~~~~~~~
# 50_irnt
# Clumped, significant loci
file_list <- paste0("interim_data/ldsc/plink_clumped_1KGP_HG19/oilyfish_tg/50_irnt/clumped_results_50_irnt_chr", 1:22,".tsv.clumped") # selected in kmeans 2, standing height
df_50 <- lapply(file_list, function(f) read.table(f, header = TRUE)) %>% bind_rows()

#~~~~~~~~~~~~~~~~~~~~
# Oily fish
# Load Exposure 2 top hits
diet_sig <- readRDS("interim_data/exposure/formatted_oilyfish_clumped_noLD_nobiallelic.rds")

#~~~~~~~~~~~~~~~~~~~~
# Merge to one SNP list
tmp_4100 <- df_4100$SNP
tmp_50 <- df_50$SNP
tmp_diet <- diet_sig$variant
snp_list <- c(tmp_4100, tmp_50, tmp_diet)
# Save snp_list
saveRDS(snp_list, "interim_data/MVMR/oilyfish_tg/snp_list.rds")

#~~~~~~~~~~~~~~~~~~~~
# Subset
# - Load GWAS sum-stat file, subset, save subset, unload GWAS sum-stat file
# - Repeat with next file

# 4100_irnt
# Full summary stats
df_4100_full <- read.table("/pl/active/colelab/users/kjames/refinedMR/interim_data/ldsc/NL_bothsexes_temp/4100_irnt.gwas.imputed_v3.both_sexes.tsv", header = TRUE)
# Extract SNPs from GWAS
df_4100_subset <- df_4100_full[df_4100_full$variant %in% snp_list,]
# Save subsetted SNPs in case R runs out of memory and closes
saveRDS(df_4100_subset, "interim_data/MVMR/oilyfish_tg/df_4100_subset.rds")
# Remove full summary stats to clear space
rm(df_4100_full)

# 50_irnt
# Full summary stats
df_50_full <- read.table("/pl/active/colelab/users/kjames/refinedMR/interim_data/ldsc/NL_bothsexes_temp/50_irnt.gwas.imputed_v3.both_sexes.tsv", header=TRUE)
# Extract SNPs from GWAS
df_50_subset <- df_50_full[df_50_full$variant %in% snp_list,]
# Save subsetted SNPs in case R runs out of memory and closes
saveRDS(df_50_subset, "interim_data/MVMR/oilyfish_tg/df_50_subset.rds")
# Remove full summary stats to clear space
rm(df_50_full)

# Oily fish
# Load full exposure 2 data 
diet_full <- readRDS("interim_data/exposure/formatted_oilyfish_clumped_full_sumstats.rds")
# Extract SNPs from GWAS
df_diet_subset <- diet_full[diet_full$variant %in% snp_list,]
# Save subsetted SNPs in case R runs out of memory and closes
saveRDS(df_diet_subset, "interim_data/MVMR/oilyfish_tg/df_diet_subset.rds")
# Remove full summary stats to clear space
rm(diet_full)

# Triglycerides
# Now add in TG (outcome)
out <- readRDS("interim_data/outcome/formatted_GCST90239664_TG_Graham_GRCh37.rds") 
# Extract SNPs from GWAS
df_out_subset <- out[out$variant %in% snp_list,]
# Save subsetted SNPs in case R runs out of memory and closes
saveRDS(df_out_subset, "interim_data/MVMR/oilyfish_tg/df_out_subset.rds")
# Remove full summary stats to clear space
rm(out)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Read in subsets
df_4100_subset <- readRDS("interim_data/MVMR/oilyfish_tg/df_4100_subset.rds")
df_50_subset <- readRDS("interim_data/MVMR/oilyfish_tg/df_50_subset.rds")
df_diet_subset <- readRDS("interim_data/MVMR/oilyfish_tg/df_diet_subset.rds")
df_out_subset <- readRDS("interim_data/MVMR/oilyfish_tg/df_out_subset.rds")

# Add column extensions
colnames(df_4100_subset) <- paste0(colnames(df_4100_subset), ".4100_irnt")
colnames(df_50_subset) <- paste0(colnames(df_50_subset), ".50_irnt")
colnames(df_diet_subset) <- paste0(colnames(df_diet_subset), ".oilyfish")
colnames(df_out_subset) <- paste0(colnames(df_out_subset), ".triglycerides")

# Merge
mvmr_data <- merge(df_4100_subset, df_50_subset, by.x = "variant.4100_irnt", by.y = "variant.50_irnt")
mvmr_data2 <- merge(mvmr_data, df_diet_subset, by.x  = "variant.4100_irnt", by.y = "variant.oilyfish")
mvmr_data3 <- merge(mvmr_data2, df_out_subset, by.x  = "variant.4100_irnt", by.y = "variant.triglycerides")

# Try MVMR
# Note, will need to get covariance of traits from UKB RAP
# Also, need to check harmonization etc.
library(MVMR)
F.data <- format_mvmr(BXGs = cbind(mvmr_data3$beta.4100_irnt, mvmr_data3$beta.50_irnt, mvmr_data3$beta.exposure.oilyfish), #beta of exposures
                      BYG = mvmr_data3$beta.outcome.triglycerides, # beta of outcome
                      seBXGs = cbind(mvmr_data3$se.4100_irnt, mvmr_data3$se.50_irnt, mvmr_data3$se.exposure.oilyfish), # se of exposures
                      seBYG = mvmr_data3$se.outcome.triglycerides, # se of outcome
                      RSID = mvmr_data3$variant.4100_irnt)
head(F.data)

# Test for weak instruments
sres <- strength_mvmr(r_input = F.data, gencov = 0) # CHANGE GENCOV WHEN YOU CALCULATE IT
# We have very weak instruments...

# Pleiotropy
pres <- pleiotropy_mvmr(r_input = F.data, gencov = 0) 
# VERY heterogeneous... not good!

# Estimate causality
res <- ivw_mvmr(r_input = F.data)

# robust causal effect estimation
#res_robust <- qhet_mvmr(F.data, iterations=100, CI = FALSE)
# Need a covariance matrix, which we have to calculate in UKB





