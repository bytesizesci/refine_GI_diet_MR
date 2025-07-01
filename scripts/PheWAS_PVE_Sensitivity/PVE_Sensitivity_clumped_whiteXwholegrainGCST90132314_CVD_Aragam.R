# MR Analysis: PheWAS-PVE - Sensitivity
# Author: KJS
# Date: 12/03/24 

# Goal: Perform percent variance explained thresholding on genetic instrument 

# Set wd
wdir <- "/pl/active/colelab/users/kjames/refinedMR//"
setwd(wdir)

# Set traits for naming files 
exposure <- "whiteXwholegrain"
outcome <- "GCST90132314_CVD_Aragam"
outcome_short <- "CVD"

OUT_pheno <- outcome # used for some naming
# EXP_pheno defined below

# For plots
EXP_name <- "White vs Whole Grain Bread"
OUT_name <- "Cardiovascular Disease"

# Load custom functions
source("scripts/PheWAS_PVE/phewas_thresh_test.R")

# Load libraries
library(data.table)
library(dplyr)
library(ggplot2)
library(ggrepel)

#~~~~~~~~~~~~~~~~~~~~~
# 1. Load PheWAS data
#~~~~~~~~~~~~~~~~~~~~~
pwas <- read.csv("interim_data/PheWAS_Grep_SNPs/PheWAS_Oilyfish_WhiteXWhole_Alcochol_UpdatedProxies.csv")
pwas <- pwas %>% select(-c(X, V1))

# Make file name prettier. Get rid of extra junk.
pwas$trait <- sub("\\..*$", "", pwas$file_name)

# PWAS includes some duplicated traits that are in both raw and 
# inverse rank normalized transformed states. 
# In these cases, we want to keep the IRNT versions. 
# Remove the traits from pwas that end in _raw
pwas <- pwas[!grepl("_raw$", pwas$trait),]

#~~~~~~~~~~~~~~~~~~~~~
# 2. Filter PheWAS traits to those with an effective sample size of 50K
# For binary traits:
# Neff = (4 x Ncases x Ncontrols)/(Ncases + Ncontrols)
# For quantitative traits: 
# Neff = samplesize [... so this step is geared towards binary traits]
#~~~~~~~~~~~~~~~~~~~~~
md <- read.table(gzfile("../../../common/published_gwas/Neale_UKB_GWAS_round2/Both_sexes/metadata/phenotypes.both_sexes.v2.tsv.bgz"),header=TRUE, sep="\t", fill=TRUE, quote = "")
# quote="" important to get a full file

#~~~~~~~~~~~~~~~~~~~~~
# 2.1. Filter based on general sample size and ordinal variable sample size
#~~~~~~~~~~~~~~~~~~~~~
# from, https://github.com/LizaDarrous/PheWAS-cluster/blob/main/R/getSNPs_PheWAS.R
n_thresh <- 50000
md$n_eff = md$n_non_missing
ncont_ind = which(!is.na(md$n_cases)) #which rows have values in n_cases 
md$n_eff[ncont_ind] = (4*md$n_cases[ncont_ind]*md$n_controls[ncont_ind]) / (md$n_cases[ncont_ind]+md$n_controls[ncont_ind]) # (4*(control*case))/totN
small_n = which(md$n_eff < n_thresh)
md = md[-small_n,]

# Only keep traits in PheWAS that have sample size >50K
pwas2 <- pwas %>% filter((trait %in% md$phenotype))

# Remove traits with  (bread traits - look at md)
breadtraits = c("1448_3", "1448_1", "1448_2", "1448_4", "100940", "1438_irnt")
pwas2 <- pwas2 %>%
  filter(!(trait %in% breadtraits))

#~~~~~~~~~~~~~~~~~~~~~
# 2.1.1 Load manual curation of traits without diet
#~~~~~~~~~~~~~~~~~~~~~
# This was created in PheWAS_Sensitivity_Remove_Diet.R and then manually curated
diet_df <- read.csv("interim_data/PheWAS_PVE/PheWAS_traitinfo_manualcuration_sensitivity_rmdiet_manual_curation.csv", row.names = 1)

# Get object with traits with TRUE 
diet_related <- na.omit(diet_df)
diet_related <- diet_related$trait

# Remove traits with TRUE - they have to do with diet
pwas2 <- pwas2[!pwas2$trait %in% diet_related,]

#~~~~~~~~~~~~~~~~~~~~~
# 2.2. Merge pwas data with other meta data
#~~~~~~~~~~~~~~~~~~~~~
# Get the manifest key for PheWAS
key <- read.csv("/pl/active/colelab/common/published_gwas/Neale_UKB_GWAS_round2/metadata/UKBB GWAS Imputed v3 - File Manifest Release 20180731 - Description Lookup.csv",header = FALSE)

# Merge to make a description column. Also, name the data frame traits to match LD's workflow.
trait_info <- merge(pwas2, key, by.x = "trait", by.y = "V1", all.x = TRUE, all.y = FALSE)
# new_name = old_name
trait_info = trait_info %>% dplyr::rename("description" = "V2")

# USER INPUT (used in section XXX)
#Find the ID of the initial exposures of interest
key_term <- "Bread consumed" #"Oily fish intake", note- alcohol consumed isn't the direct metric but should be a good proxy
trait_ID <- key[grep(key_term, key$V2),] 
trait_ID <- trait_ID[1,1]
#~~~~~~~~~~~~~~~~~~~~~
# 3. Filter PheWAS data to traits that have a genetic correlation of less than 75% with the exposure
#~~~~~~~~~~~~~~~~~~~~~
# Genetic correlation data
gcor <- read.csv("/pl/active/colelab/users/kjames/enviroMR/inputs/UKBB_genetic_correlation_data-2024-03-20.csv")

# Make phenotype names easier to read using the column with the url
gcor$phenotype.1.name <- stringr::str_extract(gcor$Phenotype.1, "(?<=\\>).*(?=\\<)")
gcor$phenotype.2.name <- stringr::str_extract(gcor$Phenotype.2, "(?<=\\>).*(?=\\<)")

# USER INPUT: Define genetic correlation threshold
gcorr_thresh <- 0.75
# USER INPUT: Define trait of interest
cat(key_term) # check 
toi_gcorr <- gcor[grep(key_term, gcor$phenotype.1.name),]
# Filter
high_gcorr = which(abs(toi_gcorr$rg)>=gcorr_thresh) #note, gcorr are + and -

# Function that checks the genetic correlations
check_gcorr <- function(high_gcorr) {
  if (length(high_gcorr) > 0) {
    cat(paste0("Remove these phenotypes:\n", high_gcorr))
  } else {
    cat(paste0("There are no phenotypes with a genetic correlation greater than:\n", gcorr_thresh ))
  }
}

# Apply function
check_gcorr(high_gcorr)

# If needed, filter out the traits from objects
# Remove from gcor
#toi_gcorr2 <- toi_gcorr[-high_gcorr,]

# Set up: 
toi_gcorr_rm <- toi_gcorr[high_gcorr,]
trait_rm <- toi_gcorr_rm$ID2

# Get rows of pwas2 that should be removed due to gcorr
rows <- pwas2$trait == trait_rm
# Remove from pwas2; does nothing if there are none
if (length(rows) == 0) {
  # Do nothing
} else {
  # Remove rows with indices
  pwas2=pwas2[!rows,]
}

# Get rows of trait_info that should be removed due to gcorr
rows <- trait_info$trait == trait_rm
# Remove from trait info; does nothing if there are none
if (length(rows) == 0) {
  # Do nothing
} else {
  # Remove rows with indices
  trait_info=trait_info[!rows,]
}
#~~~~~~~~~~~~~~~~~~~~~
# 3.5 Wrangle PheWAS data to get stdZ
#~~~~~~~~~~~~~~~~~~~~~
# Calculate standard beta - make the variant trait association column Z
pwas2$std_Z <- (pwas2$beta / pwas2$se)/sqrt(pwas2$n_complete_samples)
# View
hist(pwas2$std_Z, breaks = 50)

# For each SNP, find the SNP with the greatest PVE
pwas2_top5 = pwas2 %>% arrange(desc(PVE)) %>% group_by(variant) %>% slice(1:5)
pwas2_top1 <- pwas2 %>% arrange(desc(PVE)) %>% group_by(variant) %>% slice(1:1)

#~~~~~~~~~~~~~~~~~~~~~
# 4. Load MR data and calculate PVE
#~~~~~~~~~~~~~~~~~~~~~
# Load merged, harmonized GI data
mrobj <- readRDS("interim_data/merged_GI/MRobject_steiger_241113.rds")

# Programatically get the right component of the list
exposure_outcome <- paste0(exposure, "_", outcome_short)
mrdat <- mrobj[[exposure_outcome]]

# Calculate PVE for diet trait
MAF <- mrdat$eaf.exposure
N  <- mrdat$samplesize.exposure
se <- mrdat$se.exposure
beta <- mrdat$beta.exposure

#Calculate PVE variable
numerator <- 2*(beta^2)*MAF*(1-MAF)
denominator <- (2*(beta^2)*MAF*(1-MAF)+((se^2)*2*N*MAF*(1-MAF))) 
mrdat$PVE.exposure <- numerator/denominator

#~~~~~~~~~~~~~~~~~~~~~
# 4.1 Subset PheWAS data to SNPs in mrdat
#~~~~~~~~~~~~~~~~~~~~~
#Filter phewas based on the variants in the MR genetic instrument

pwas2$chr_pos <- regmatches(pwas2$variant, regexpr("[a-z0-9]+\\:[a-z0-9]+", pwas2$variant, ignore.case = TRUE)) 
mrdat$chr_pos <- regmatches(mrdat$variant, regexpr("[a-z0-9]+\\:[a-z0-9]+", mrdat$variant, ignore.case = TRUE)) 

phewas_trait <- pwas2[pwas2$chr_pos %in% mrdat$chr_pos,]

#~~~~~~~~~~~~~~~~~~~~~
# 5. Apply PVE Filter
#~~~~~~~~~~~~~~~~~~~~~
# For a given SNP, see if its PVE value is higher than the PVE of the trait of interest. 
phewas_thresh_results = phewas_thresh_test(variant_ID_eff_alt = mrdat$chr_pos, expos_df = mrdat, phewas_sub = phewas_trait)

#name list
names(phewas_thresh_results) <- mrdat$chr_pos

#get the names of the slots of the list that have 0 rows in them (these are variants we want to keep!)
# Function to check if a data frame has dimensions n>0 x 7
#is_valid_df_keep <- function(x) {
#  is.data.frame(x) && nrow(x) > 0 && ncol(x) == 7
#}
# Get the names of the slots with valid data frames
#keep_variants <- names(phewas_thresh_results)[sapply(phewas_thresh_results, is_valid_df_keep)]

keep_variants <- phewas_thresh_results[sapply(phewas_thresh_results, function(df) is.data.frame(df) || (is.data.frame(df) && nrow(df) != 0))] #historical method... yikes - got rid of a !

# Get names
keep_variant_names = names(keep_variants) # use this to filter data
keep_variant_names #42 with this method

#Which variant(s) explain less variance?
remove_variants = phewas_thresh_results[sapply(phewas_thresh_results, function(df) !is.data.frame(df) || (is.data.frame(df) && nrow(df) == 0))]
# Get names
remove_variants_names = names(remove_variants)
remove_variants_names 

#Remove them from MR data
filtered_pve <- mrdat %>% filter(!(chr_pos %in% remove_variants_names))
dim(filtered_pve)

# Save items for future use
save(list=c("phewas_thresh_results", "keep_variants", "keep_variant_names", "remove_variants", "remove_variants_names", "filtered_pve"), file = paste0("interim_data/PheWAS_PVE_Sensitivity/PheWAS_PVE_interimObjects_",exposure,"_",outcome,"_Sensitivity.RData"))

#~~~~~~~~~~~~~~~~~~~~~
# 5.1 Add keep column to mrobj based on PVE filter
#~~~~~~~~~~~~~~~~~~~~~
#mrdat$PheWAS_PVE_keep <- ifelse(mrdat$chr_pos %in% remove_variants_names, FALSE, TRUE)

# Put back in slot/component position
#mrobj[[exposure_outcome]] <- mrdat

# Save
#saveRDS(mrobj, paste0("interim_data/merged_GI/standardized_shared/std_mrobject_steiger_PheWASClust_PVE_240930e.rds"))

#~~~~~~~~~~~~~~~~~~~~~
# 5.2 Look at what the removed SNPs were related to
#~~~~~~~~~~~~~~~~~~~~~
# Get a list of unique traits/descriptions 
trait_info_unique <- trait_info %>% select(trait, description) %>% unique()

# Filter pwas_top5 for variants that were removed
instr_rmd_top5 <- pwas2_top5 %>% filter(chr_pos %in% remove_variants_names)
instr_rmd_top5_info <- merge(instr_rmd_top5, trait_info_unique, by = "trait", all.x = TRUE, all.y = FALSE)

# Filter pwas_top1 for variants that were removed
instr_rmd_top1 <- pwas2_top1 %>% filter(chr_pos %in% remove_variants_names)
instr_rmd_top1_info <- merge(instr_rmd_top1, trait_info_unique, by = "trait", all.x = TRUE, all.y = FALSE)

# Save
phewas_table_dir <- "/pl/active/colelab/users/kjames/refinedMR/results/PheWAS_PVE_Sensitivity/tables"
# List of top 5 hits PVE for removed variants
write.csv(instr_rmd_top5_info, paste0(phewas_table_dir, "/MR_PheWAS_PVE_RemovedVariants_top5_", exposure,"_", outcome, "_Sensitivity.csv"))
# List of top 1 hit PVE for removed variants
write.csv(instr_rmd_top1_info, paste0(phewas_table_dir, "/MR_PheWAS_PVE_RemovedVariants_top1_", exposure,"_", outcome, "_Sensitivity.csv"))

#~~~~~~~~~~~~~~~~~~~~~
# 6. Perform MR
#~~~~~~~~~~~~~~~~~~~~~

# Get variant column
df_filt <- filtered_pve

# PVE filtered object
MRobj_raw <- MendelianRandomization::mr_input(bx = df_filt$beta.exposure,
                                              bxse = df_filt$se.exposure,
                                              by = df_filt$beta.outcome,
                                              byse = df_filt$se.outcome)

MRres_raw <- MendelianRandomization::mr_allmethods(MRobj_raw)
MRres_raw

# Get df
MRres_df <- MRres_raw@Values

# Save tables
# All MR results
write.csv(MRres_df, paste0(phewas_table_dir, "/MR_PheWAS_PVE_Filtered_All_Results_", exposure,"_", outcome, "_Sensitivity.csv"))

# Run IVW 
mr_ivw <- MendelianRandomization::mr_ivw(MRobj_raw)
# Create a table of results for saving
mr_ivw_df <- data.frame(test = "mr_ivw",
                        model = mr_ivw@Model,
                        robust = mr_ivw@Robust,
                        penalized = mr_ivw@Penalized,
                        estimate = mr_ivw@Estimate,
                        std_error = mr_ivw@StdError,
                        cilower = mr_ivw@CILower,
                        ciupper = mr_ivw@CIUpper,
                        pvalue = mr_ivw@Pvalue,
                        heter_Q = mr_ivw@Heter.Stat[1],
                        heter_pvalue = mr_ivw@Heter.Stat[2],
                        fstat = mr_ivw@Fstat)

# Save table
write.csv(mr_ivw_df, paste0(phewas_table_dir, "/MR_PheWAS_PVE_Filtered_IVW_Results_", exposure,"_", outcome, "_Sensitivity.csv"))

#~~~~~~~~~~~~~~~~~~~~~
# 10. Plot
#~~~~~~~~~~~~~~~~~~~~~
# Plot with all methods
MRres_p <- MendelianRandomization::mr_plot(MRres_raw, error=TRUE, orientate=FALSE, line="ivw", labels=TRUE) +
  ggplot2::theme_bw() +
  ggplot2::theme(
    #legend.position="none", 
    axis.title.y=ggplot2::element_text(size=14),
    axis.text.y=ggplot2::element_text(size=14), 
    axis.ticks.y=ggplot2::element_line(linewidth=0),
    axis.title.x=ggplot2::element_text(size=14),
    axis.text.x=ggplot2::element_text(size=14),
    axis.ticks.x=ggplot2::element_line(linewidth=0),
    panel.border = element_rect(color = "black", 
                                fill = NA, 
                                linewidth = 3)) +
  labs(x = paste0("Genetic Effect on ", EXP_name),
       y = paste0( "Genetic Effect on ", OUT_name)) +
  ggtitle(paste0("PheWAS-PVE Filtered Significant Variants (n=", nrow(df_filt), ")")) 
#scale_color_viridis(discrete = TRUE)

MRres_p

# USER INPUT: Check paths and variable names
plots_dir <- "/pl/active/colelab/users/kjames/refinedMR/results/PheWAS_PVE_Sensitivity/plots"

# Save
ggsave(paste0(plots_dir,"/PheWAS_PVE_MR_",exposure, "_", outcome,"_scatterplot_allmethods_Sensitivity.png"),
       plot = MRres_p,
       device = NULL,
       path = NULL,
       scale = 1,
       width = 6,
       height = 4,
       units = c("in"),
       dpi = 300,
       limitsize = TRUE)

# Get values and subset to weighted median and IVW for custom plot below
mrRes_df_subsetMethods <- MRres_df %>% filter(Method %in% c("Weighted median", "IVW"))

# Custom plot for external facing meetings
custom_MRscatter <- ggplot(df_filt, aes(x = beta.exposure, y = beta.outcome)) +
  geom_point() +
  ggplot2::geom_errorbar(ggplot2::aes(ymin=beta.outcome-se.outcome, ymax=beta.outcome+se.outcome), colour="grey", width=0) +
  ggplot2::geom_errorbarh(ggplot2::aes(xmin=beta.exposure-se.exposure, xmax=beta.exposure+se.exposure), colour="grey", height=0) +
  ggplot2::geom_abline(data = mrRes_df_subsetMethods, # Specify correct data here
                       ggplot2::aes(intercept=0, slope=Estimate, colour=Method), size = 1.75, show.legend=TRUE) + #0 if not MR-Egger
  ggplot2::theme_bw() +
  ggplot2::theme(
    #legend.position="none", 
    axis.title.y=ggplot2::element_text(size=14),
    axis.text.y=ggplot2::element_text(size=14), 
    axis.ticks.y=ggplot2::element_line(size=0),
    axis.title.x=ggplot2::element_text(size=14),
    axis.text.x=ggplot2::element_text(size=14),
    axis.ticks.x=ggplot2::element_line(size=0),
    panel.border = element_rect(color = "black", 
                                fill = NA, 
                                size = 3)) +
  labs(x = paste0("Genetic Effect on ", EXP_name),
       y = paste0( "Genetic Effect on ", OUT_name)) +
  ggtitle(paste0("PheWAS-PVE Filtered Genetic Instrument (n=", nrow(df_filt), ")")) +
  scale_color_manual(values = c("#3E4A89FF", "#6DCD59FF"))

custom_MRscatter

ggsave(paste0(plots_dir,"/PheWAS_PVE_MR_",exposure, "_", outcome,"_scatterplot_IVW_weightedmedn_Sensitivity.png"),
       plot = custom_MRscatter,
       device = NULL,
       path = NULL,
       scale = 1,
       width = 6,
       height = 4,
       units = c("in"),
       dpi = 300,
       limitsize = TRUE)

