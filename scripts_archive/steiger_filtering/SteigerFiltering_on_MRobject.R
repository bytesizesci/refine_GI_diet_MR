# Goal: Perform Steiger filtering MR - updated with new proxy exposure SNPs and neg controls 
# Author: KJS
# Date: November 13, 2024

#Source of steiger code from TwoSampleMR:
# * Steiger function - https://rdrr.io/github/MRCIEU/TwoSampleMR/src/R/steiger_filtering.R
# * R-squared functions - https://rdrr.io/github/MRCIEU/TwoSampleMR/src/R/add_rsq.r 

#Use: 
# * From our mrdata structure, add or rename the exposure column to "id.exposure" and the outcome column to "id.outcome"
# * Load the functions
# * Run steiger_filtering function
# * Look at direction and p-values to filter variants that may be acting from outcome to exposure

# Set wd
wdir <- "/pl/active/colelab/users/kjames/refinedMR/"
setwd(wdir)

# Load libraries
library(dplyr)
library(ggplot2)

# Source functions
source("/pl/active/colelab/users/kjames/refinedMR/scripts/steiger_filtering/steiger_function.R")
source("/pl/active/colelab/users/kjames/refinedMR/scripts/steiger_filtering/r_square_functions.R")

# Create a data frame with naming conventions for loop
df_names <- data.frame(EXP_pheno = c("oilyfish", "oilyfish", "oilyfish", "whiteXwholegrain", "whiteXwholegrain",  "whiteXwholegrain", "alcohol", "alcohol","alcohol"),
                       OUT_pheno = c("GCST90239664_TG_Graham", "GCST90132314_CVD_Aragam", "GCST006901_Height_Yengo_2018", "GCST90239658_LDL_Graham","GCST90132314_CVD_Aragam", "GCST006901_Height_Yengo_2018", "GCST90013405_ALT_Pazoki", "GCST90319877_cir_Ghouse", "GCST006901_Height_Yengo_2018"),
                       EXP_names = c("Oily Fish", "Oily fish", "Oily fish", "White vs Whole Grain Bread", "White vs Whole Grain Bread", "White vs Whole Grain Bread","Alcohol", "Alcohol", "Alcohol"),
                       OUT_names = c("Triglycerides", "Cardiovascular Disease", "Height", "LDL Cholesterol", "Cardiovascular Disease", "Height", "Alanine Aminotransferase", "Liver Cirrhosis", "Height"))

#~~~~~~~~~~~~
# 1. 
# Run this section once
# Skip to second section if just tweaking (like saving more of) the MR output
#~~~~~~~~~~~~
# Load mrobject with harmonized and standardized genetic instruments
mrobj <- readRDS("interim_data/merged_GI/MRobject_241113.rds")

# Note, new column from steiger function include
# units.exposure
# units.outcome
# rsq.exposure
# effective_n.exposure
# rsq.outcome
# effective_n.outcome
# steiger_dir
# steiger_pval

# Build steiger function into a loop
for (i in seq_along(mrobj)){
  # Get data frame
  mrdat <- mrobj[[i]]
  # Get names
  EXP_pheno <- df_names[i,1]
  OUT_pheno <- df_names[i,2]
  EXP_name <- df_names[i,3]
  OUT_name <- df_names[i,4]
  # Check for important column names for steiger functions
  mrdat <- assign_exposure_col(mrdat, EXP_name = EXP_name, EXP_pheno = EXP_pheno)
  mrdat <- assign_outcome_col(mrdat, OUT_name = OUT_name, OUT_pheno = OUT_pheno)
  # New - Update exposure column so that they are all identical
  # The reason they are different is that if the SNP was a proxy I used a different
  # blurb to describe the exposure trait "trait_fullGWAS" instead of the blurb
  # used in the original GI. This causes an error in the steiger function.
  mrdat$exposure <- EXP_pheno
  # Apply steiger test
  mrdat_stg = steiger_filtering(mrdat)
  # Reassign to mrobj in correct component slot
  mrobj[[i]] <- mrdat_stg
}

# Save mrobj
saveRDS(mrobj, "interim_data/merged_GI/MRobject_steiger_241113.rds")

#~~~~~~~~~~~~
# 2. 
# Wanting to use the same input but save more output?
# Start here
#~~~~~~~~~~~~
# Read in, if needed
mrobj <- readRDS("interim_data/merged_GI/MRobject_steiger_241113.rds")

# Perform MR 
for (i in seq_along(mrobj)){
  # Get data frame
  mrdat <- mrobj[[i]]
  # Get trait names
  # Get names
  EXP_pheno <- df_names[i,1]
  OUT_pheno <- df_names[i,2]
  EXP_name <- df_names[i,3]
  OUT_name <- df_names[i,4]
  
  # Perform Steiger filtering
  # Keep snps that indicate they go from X->Y (TRUE) or have a p-value greater than 0.05 (we don't reject the null - no difference) 
  mrdat_stgfilt <- mrdat[mrdat$steiger_dir == "TRUE" | mrdat$steiger_pval >= 0.05,]
  # Perform MR
  mrInput_stgfilt <- MendelianRandomization:: mr_input(bx = mrdat_stgfilt$beta.exposure,
                                             bxse = mrdat_stgfilt$se.exposure,
                                             by = mrdat_stgfilt$beta.outcome,
                                             byse = mrdat_stgfilt$se.outcome)
  mrRes_stgfilt <- MendelianRandomization::mr_allmethods(mrInput_stgfilt, method = "all")
  mrRes_stgfilt_df <- mrRes_stgfilt@Values
  
  # Save table
  table_dir <- "/pl/active/colelab/users/kjames/refinedMR/results/steiger/tables/"
  write.csv(mrRes_stgfilt_df, paste0(table_dir, "/SteigerfilteredGI_MR_", EXP_pheno, "_", OUT_pheno, ".csv"))
  
  # Get values and subset to weighted median and IVW for custom plot below
  mrRes_stgfilt_df_subsetMethods <- mrRes_stgfilt_df %>% filter(Method %in% c("Weighted median", "IVW"))
  
  # Custom plot for external facing meetings
  custom_MRscatter_StgFlt <- ggplot(mrdat_stgfilt, aes(x = beta.exposure, y = beta.outcome)) +
    geom_point() +
    ggplot2::geom_errorbar(ggplot2::aes(ymin=beta.outcome-se.outcome, ymax=beta.outcome+se.outcome), colour="grey", width=0) +
    ggplot2::geom_errorbarh(ggplot2::aes(xmin=beta.exposure-se.exposure, xmax=beta.exposure+se.exposure), colour="grey", height=0) +
    ggplot2::geom_abline(data = mrRes_stgfilt_df_subsetMethods, # Specify correct data here
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
    ggtitle(paste0("Steiger-filtered Genetic Instrument (n=", nrow(mrdat_stgfilt), ")")) +
    scale_color_manual(values = c("#3E4A89FF", "#6DCD59FF"))
  
  custom_MRscatter_StgFlt
  
  plots_dir <- "/pl/active/colelab/users/kjames/refinedMR/results/steiger/plots"
  ggsave(paste0(plots_dir,"/SteigerfilteredGI_MR_scatterplot_",EXP_pheno, "_", OUT_pheno,".png"),
         plot = custom_MRscatter_StgFlt,
         device = NULL,
         path = NULL,
         scale = 1,
         width = 6,
         height = 4,
         units = c("in"),
         dpi = 300,
         limitsize = TRUE)
  
  # Run IVW to get the F stat and Het stats
  mr_ivw <- MendelianRandomization::mr_ivw(mrInput_stgfilt) #object with filtered data 
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
  write.csv(mr_ivw_df, paste0("/pl/active/colelab/users/kjames/refinedMR/results/steiger/tables/MR_Steiger_Filtered_GI_IVWres_", EXP_pheno,"_", OUT_pheno, ".csv"))
  
  # Now perform non-filtered MR and combine with steiger filtering for reporting ease
  # Perform MR on non-filtered data
  mrInput <- MendelianRandomization:: mr_input(bx = mrdat$beta.exposure,
                                                       bxse = mrdat$se.exposure,
                                                       by = mrdat$beta.outcome,
                                                       byse = mrdat$se.outcome)
  mrRes <- MendelianRandomization::mr_allmethods(mrInput, method = "all")
  mrRes_df <- mrRes@Values
  
  # Make table
  # For S4 objects, access slots with @
  # Steiger filtered results
  mrRes_stgfilt_df$GI <- "filtered"
  mrRes_stgfilt_df$N_snps <- nrow(mrdat_stgfilt)
  
  # Unfiltered results
  mrRes_df$GI <- "all"
  mrRes_df$N_snps <- nrow(mrdat)
  
  # Create a data frame that shows the non-filtered and filtered data next to each other for sharing
  # Number of rows in each data frame
  n_rows_all <- nrow(mrRes_df)
  n_rows_filt <- nrow(mrRes_stgfilt_df)
  
  # Create an empty data frame C
  df_pub <- data.frame(ID = integer(), Value = character())
  
  # Loop through rows of A and B alternatively and combine them into C
  for (j in 1:max(n_rows_all, n_rows_filt)) {
    if (j <= n_rows_all) {
      df_pub <- rbind(df_pub, mrRes_df[j, ])
    }
    if (j <= n_rows_filt) {
      df_pub <- rbind(df_pub, mrRes_stgfilt_df[j, ])
    }
  }
  
  # Save the combined table
  combined_table_dir <- "/pl/active/colelab/users/kjames/refinedMR/results/steiger/tables"
  write.csv(df_pub, paste0(combined_table_dir, "/MR_NonFiltered_and_SteigerFiltered_Results_Formatted", EXP_pheno,"_", OUT_pheno, ".csv"))
}


