# Goal: Harmonize of exposure and outcome data for bidirectional MR
# Outcome: alcohol
# Exposure: alanine aminotransferase, ALT
# Author: KJS
# Date: Oct 10, 2024

# Load libraries
library(dplyr)
library(LDlinkR)
library(ggplot2)
library(ggrepel)
library(viridis)

# Set wd
wdir <- "/pl/active/colelab/users/kjames/enviroMR/"
setwd(wdir)

exposure <- "GCST90013405_ALT_Pazoki"
outcome <- "alcohol"

# Load Exposure
ex <- readRDS("interim_data/outcome/formatted_ST6validated_GCST90013405_ALT_Pazoki.rds")

# Load Outcome
#About diet data: https://kp4cd.org/sites/default/files/READMEs/Cole_UKB_Diet_GWAS_README.pdf
# * GRCh37/hg19
# * ALLELE1 = first allele in UKB bim file, effect allele
# * ALLELE0 = second allele in UKB bim file, non-effect allele
out <- read.table(gzfile("/pl/active/colelab/common/published_gwas/diet/cole_naturecomm_2020/BOLTlmm_UKB_genoQCEURN455146_v3_diet_alcohol_overallfreq.1558.average_QT_BgenSnps_mac20_maf0.005_info0.6.gz"), header = TRUE)
colnames(out)
out$variant_id <- paste0(out$CHR, ":", out$BP, ":", out$ALLELE0, ":", out$ALLELE1)

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

# Rename Outcome (Alcohol) columns
out <- out %>%
  rename(#"SNP" = "rsid",
    "beta.outcome"="BETA",
    "pval.outcome"="P_BOLT_LMM",
    "chr.outcome"="CHR",
    "se.outcome"="SE",
    "pos.outcome"="BP",
    "effect_allele.outcome" = "ALLELE1",
    "other_allele.outcome" = "ALLELE0",
    "eaf.outcome" = "A1FREQ"
    #"samplesize.outcome" = "n"
  )

# Make chr_pos merge column
ex$merge_ID <- paste0(ex$chr.exposure,":", ex$pos.exposure)
out$merge_ID <- paste0(out$chr.outcome,":", out$pos.outcome)

# Merge to make genetic instrument (chr:bp:other:effect)
gi = merge(ex, out, by = "merge_ID", 
           all.x = TRUE, all.y = FALSE)

# Make a oaf.outcome column (needed for harmonization scripts)
gi$oaf.outcome <- 1 - gi$eaf.outcome

# Check strand ambiguity
# Use custom function; Flags A/T and C/G SNPs
source("scripts/is_ambiguous.R")
# Identify ambiguous SNPs
gi$ambiguous <- ifelse(is_ambiguous(gi$effect_allele.exposure, gi$other_allele.exposure) == TRUE, "ambiguous", "ok")
table(gi$ambiguous)

# Check allele frequency 
# If 0.45 < allele frequency < 0.55 then we will say we can't confidently say which strand is which and we will remove them
# Exposure (ALT) eaf is missing so use outcome (alcohol)
gi$remove <- ifelse(gi$ambiguous == "ambiguous" & gi$eaf.outcome <= 0.55 & gi$eaf.outcome >= 0.45, "remove", "ok")
table(gi$remove)

# Apply filter
gi <- gi[gi$remove == "ok",] # only keep "ok" not "remove"

# FOR BIDIRECTIONAL MR 0F ALT-ALCOHOL WE CANNOT CHECK THIS BC ALT DATA EAF IS MISSING
# Also, check that the exposure and outcome allele frequencies are within +-0.05 of each other. 
# They should be (especially if input is from the same ancestry)
# Add a new column to flag discrepancies
#gi$exp_out_af_flag <- ""
# Define the threshold for flagging discrepancies
#threshold <- 0.05
# Function to flag discrepancies
#gi$exp_out_af_flag <- ifelse(abs(gi$eaf.exposure - gi$eaf.outcome) > threshold, "Flag", "ok")
#table(gi$exp_out_af_flag) 
# Don't filter yet, might be fixed when harmonized fully

# Get list of complete/matching SNPs
gi_complete = gi[complete.cases(gi), ]
dim(gi_complete) #145

# Harmonize complete data
# Note, this harmonization script was adapted to account for the missing eaf.exposure variable in ALT data
source("scripts/harmonize_alleles_noEAFexp.R")
gi_harmon <- harmonize_alleles(gi_complete)
gi_harmon <- gi_harmon$data

# Any changes?
identical(gi_harmon, gi_complete) # FALSE - harmonization was needed

# Get list of missing outcome SNPs
gi_needproxy <- gi[!complete.cases(gi), ]
dim(gi_needproxy) # 0

# Save
saveRDS(gi_harmon, "interim_data/merged_GI/mrdat_bidirectional_ST6_GCST90013405_ALT_Pazoki_alcohol.rds")
write.csv(gi_harmon, "interim_data/merged_GI/mrdat_bidirectional_ST6_GCST90013405_ALT_Pazoki_alcohol.csv")

# Read if needed
gi_harmon <- readRDS("interim_data/merged_GI/mrdat_bidirectional_ST6_GCST90013405_ALT_Pazoki_alcohol.rds")

# Out of an abundance of curiosity
MRobj <- MendelianRandomization:: mr_input(bx = gi_harmon$beta.exposure,
                                               bxse = gi_harmon$se.exposure,
                                               by = gi_harmon$beta.outcome,
                                               byse = gi_harmon$se.outcome)
MRres <- MendelianRandomization::mr_allmethods(MRobj, method = "all")
MRres

EXP_name <- "Alanine Aminotransferase"
OUT_name <- "Alcohol"

# Get df
MRres_df <- MRres@Values
# Get values and subset to weighted median and IVW for custom plot below
mrRes_df_subsetMethods <- MRres_df %>% filter(Method %in% c("Weighted median", "IVW"))

# Get heterogeneity statistics 
MRres_ivw <- MendelianRandomization::mr_ivw(MRobj)
MRres_ivw_het <-  MRres_ivw@Heter.Stat # cochran's Q, p-value

# Save tables
# Path
phewas_table_dir <- "/pl/active/colelab/users/kjames/enviroMR/results/julie_collab/bidirectional_MR/tables"
# MR results
write.csv(MRres_df, paste0(phewas_table_dir, "/MR_bidirectional_Filtered_Results_", exposure,"_", outcome, ".csv"))
# Heterogeneity stats from IVW method
write.csv(MRres_ivw_het, paste0(phewas_table_dir, "/MR_bidirectional_Filtered_IVW_CochransQ_Pval_", exposure,"_", outcome, ".csv"))


# Plot all methods
MRres_p <- MendelianRandomization::mr_plot(MRres, error=TRUE, orientate=FALSE, line="ivw", labels=TRUE) +
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
  ggtitle(paste0("Bi-directional MR (n=", nrow(gi_harmon), ")")) 
#scale_color_viridis(discrete = TRUE)

MRres_p

# Custom plot for external facing meetings
custom_MRscatter <- ggplot(gi_harmon, aes(x = beta.exposure, y = beta.outcome)) +
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
  ggtitle(paste0("Bi-directional MR (n=", nrow(gi_harmon), ")")) +
  scale_color_manual(values = c("#3E4A89FF", "#6DCD59FF"))

custom_MRscatter

plots_dir <- "/pl/active/colelab/users/kjames/enviroMR/results/julie_collab/bidirectional_MR/plots/"
ggsave(paste0(plots_dir,"/bidirectionalMR_",exposure, "_", outcome,"_scatterplot_IVE_weightedmedn.png"),
       plot = custom_MRscatter,
       device = NULL,
       path = NULL,
       scale = 1,
       width = 6,
       height = 4,
       units = c("in"),
       dpi = 300,
       limitsize = TRUE)

# Load MR results in the alcohol to ALT direction, extract key results, and make plot showing results for both analyses
mrres_alc_alt <- read.csv("results/julie_collab/MR_nonfilteredGI/NonfilteredGI_MR_Alcohol_Alanine Aminotransferase.csv")
mrres_alc_alt

bidir_dat <- data.frame(
  Direction = c("Alcohol on Alanine Aminotransferase", "Alanine Aminotransferase on Alcohol"),
  IVW_Estimate = c(-2.173830E-02, -0.1710353974),
  StdError = c(0.0069438635, 0.0756395077)
)

bdmr_scatter <- ggplot(bidir_dat, aes(x = Direction, y = IVW_Estimate)) +
  geom_point(size = 3, color = "grey") +
  geom_errorbar(aes(ymin = IVW_Estimate - 1.96 * StdError, ymax = IVW_Estimate + 1.96 * StdError), 
                width = 0.2, color = "#287c8EFF") +
  labs(title = "Bi-directional MR IVW Results", 
       x = "Direction", 
       y = "IVW Effect Estimate") +
  coord_flip() + 
  theme_bw()  
  
ggsave(paste0(plots_dir,"/bidirectionalMR_scatter_",exposure, "_", outcome,"_scatterplot_IVW.png"),
       plot = bdmr_scatter,
       device = NULL,
       path = NULL,
       scale = 1,
       width = 6,
       height = 4,
       units = c("in"),
       dpi = 300,
       limitsize = TRUE)
  
  
  
  
  




