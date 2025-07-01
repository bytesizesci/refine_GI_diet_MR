# Goal: Perform classic MR - update: with sensitivity PheWAS
# Author: KJS
# Date: December 3, 2024

# Set wd
wdir <- "/pl/active/colelab/users/kjames/refinedMR/"
setwd(wdir)

# Load libraries
library(dplyr)
library(ggplot2)
library(MendelianRandomization)

# Load harmonized and standardized genetic instruments
#mrobj <- readRDS("interim_data/merged_GI/MRobject_steiger_PheWAStt_PheWASPVE_241120.rds")
mrobj <- readRDS("interim_data/merged_GI/MRobject_steiger_PheWAStt_PheWASPVE_sensitivity_241203.rds") # has sensitivity PheWAS columns in it

# Create a data frame with naming conventions for loop
df_names <- data.frame(EXP_pheno = c("oilyfish", "oilyfish", "oilyfish", "whiteXwholegrain", "whiteXwholegrain",  "whiteXwholegrain", "alcohol", "alcohol","alcohol"),
                       OUT_pheno = c("GCST90239664_TG_Graham", "GCST90132314_CVD_Aragam", "GCST006901_Height_Yengo_2018", "GCST90239658_LDL_Graham","GCST90132314_CVD_Aragam", "GCST006901_Height_Yengo_2018", "GCST90013405_ALT_Pazoki", "GCST90319877_cir_Ghouse", "GCST006901_Height_Yengo_2018"),
                       EXP_names = c("Oily Fish", "Oily fish", "Oily fish", "White vs Whole Grain Bread", "White vs Whole Grain Bread", "White vs Whole Grain Bread","Alcohol", "Alcohol", "Alcohol"),
                       OUT_names = c("Triglycerides", "Cardiovascular Disease", "Height", "LDL Cholesterol", "Cardiovascular Disease", "Height", "Alanine Aminotransferase", "Liver Cirrhosis", "Height"))

# Perform MR 
mr_wm_list <- list()
for (i in seq_along(mrobj)){
  # Get data frame
  mrdat <- mrobj[[i]]
  # Get trait names
  # Get names
  EXP_pheno <- df_names[i,1]
  OUT_pheno <- df_names[i,2]
  EXP_name <- df_names[i,3]
  OUT_name <- df_names[i,4]
  
  # For no-filtering instrument, use the whole data frame. Comment out filtering
  #mrdat_stgfilt <- mrdat[mrdat$steiger_dir == "TRUE" | mrdat$steiger_pval >= 0.05,]
  # Perform MR
  mrInput <- MendelianRandomization:: mr_input(bx = mrdat$beta.exposure,
                                                       bxse = mrdat$se.exposure,
                                                       by = mrdat$beta.outcome,
                                                       byse = mrdat$se.outcome)
  mrRes <- MendelianRandomization::mr_allmethods(mrInput, method = "all")
  mrRes_df <- mrRes@Values
  
  # Save table
  table_dir <- "/pl/active/colelab/users/kjames/refinedMR/results/Nonfiltered_GI/tables/"
  write.csv(mrRes_df, paste0(table_dir, "/NonfilteredGI_MR_", EXP_pheno, "_", OUT_pheno, ".csv"))
  
  # Get values and subset to weighted median and IVW for custom plot below
  mrRes_df_subsetMethods <- mrRes_df %>% filter(Method %in% c("Weighted median", "IVW"))
  
  # Custom plot for external facing meetings
  custom_MRscatter <- ggplot(mrdat, aes(x = beta.exposure, y = beta.outcome)) +
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
    ggtitle(paste0("Non-filtered Genetic Instrument (n=", nrow(mrdat), ")")) +
    scale_color_manual(values = c("#3E4A89FF", "#6DCD59FF"))
  
  custom_MRscatter
  
  plots_dir <- "/pl/active/colelab/users/kjames/refinedMR/results/Nonfiltered_GI/plots"
  ggsave(paste0(plots_dir,"/NonfilteredGI_MR_scatterplot_",EXP_pheno, "_", OUT_pheno,".png"),
         plot = custom_MRscatter,
         device = NULL,
         path = NULL,
         scale = 1,
         width = 6,
         height = 4,
         units = c("in"),
         dpi = 300,
         limitsize = TRUE)
  #~~~~~~~~~~~~~
  # Run IVW (fixed effects - default) to get the F stat and Het stats
  #~~~~~~~~~~~~~
  mr_ivw <- MendelianRandomization::mr_ivw(mrInput) #object with filtered data 
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
  write.csv(mr_ivw_df, paste0("/pl/active/colelab/users/kjames/refinedMR/results/Nonfiltered_GI/tables/MR_NonFiltered_GI_IVWres_", EXP_pheno,"_", OUT_pheno, ".csv"))
  #~~~~~~~~~~~~~
  # Non-filtered: Run IVW with random effects due to large amounts of heterogeneity
  #~~~~~~~~~~~~~
  mr_ivw_random <- MendelianRandomization::mr_ivw(mrInput,
                                           model = "random") #object with filtered data 
  # Create a table of results for saving
  mr_ivw_random_df <- data.frame(test = "mr_ivw_randomeffects",
                          model = mr_ivw_random@Model,
                          robust = mr_ivw_random@Robust,
                          penalized = mr_ivw_random@Penalized,
                          estimate = mr_ivw_random@Estimate,
                          std_error = mr_ivw_random@StdError,
                          cilower = mr_ivw_random@CILower,
                          ciupper = mr_ivw_random@CIUpper,
                          pvalue = mr_ivw_random@Pvalue,
                          heter_Q = mr_ivw_random@Heter.Stat[1],
                          heter_pvalue = mr_ivw_random@Heter.Stat[2],
                          fstat = mr_ivw_random@Fstat)
  
  # Save table
  write.csv(mr_ivw_random_df, paste0("/pl/active/colelab/users/kjames/refinedMR/results/Nonfiltered_GI/tables/MR_NonFiltered_GI_IVWres_RANDOMeffects_", EXP_pheno,"_", OUT_pheno, ".csv"))
  
  #~~~~~~~~~~~~~
  # Now perform PheWAS-PVE MR and combine with non-filtered for reporting ease
  #~~~~~~~~~~~~~
  mrdat_pve_filt <- mrdat[mrdat$PheWAS_PVE_keep == "TRUE",]
  # Perform MR on pve-filtered data
  mrInput_pve_filt <- MendelianRandomization:: mr_input(bx = mrdat_pve_filt$beta.exposure,
                                               bxse = mrdat_pve_filt$se.exposure,
                                               by = mrdat_pve_filt$beta.outcome,
                                               byse = mrdat_pve_filt$se.outcome)
  mrRes_pve_filt <- MendelianRandomization::mr_allmethods(mrInput_pve_filt, method = "all")
  mrRes_pve_filt_df <- mrRes_pve_filt@Values
  
  #~~~~~~~~~~~~~
  # Now perform PheWAS-PVE SENSITIVITY MR and combine with non-filtered for reporting ease
  #~~~~~~~~~~~~~
  mrdat_pve_sens_filt <- mrdat[mrdat$PheWAS_PVE_sensitivity_keep == "TRUE",]
  # Perform MR on pve-filtered data
  mrInput_pve_sens_filt <- MendelianRandomization:: mr_input(bx = mrdat_pve_sens_filt$beta.exposure,
                                                        bxse = mrdat_pve_sens_filt$se.exposure,
                                                        by = mrdat_pve_sens_filt$beta.outcome,
                                                        byse = mrdat_pve_sens_filt$se.outcome)
  mrRes_pve_sens_filt <- MendelianRandomization::mr_allmethods(mrInput_pve_sens_filt, method = "all")
  mrRes_pve_sens_filt_df <- mrRes_pve_sens_filt@Values
  
  #~~~~~~~~~~~~~
  # Now perform PheWAS-Ttest MR and combine for reporting ease
  #~~~~~~~~~~~~~
  mrdat_ttest_filt <- mrdat[mrdat$PheWAS_ttest_keep == "TRUE",]
  # Perform MR on ttest-filtered data
  mrInput_ttest_filt <- MendelianRandomization:: mr_input(bx = mrdat_ttest_filt$beta.exposure,
                                                        bxse = mrdat_ttest_filt$se.exposure,
                                                        by = mrdat_ttest_filt$beta.outcome,
                                                        byse = mrdat_ttest_filt$se.outcome)
  mrRes_ttest_filt <- MendelianRandomization::mr_allmethods(mrInput_ttest_filt, method = "all")
  mrRes_ttest_filt_df <- mrRes_ttest_filt@Values
  
  #~~~~~~~~~~~~~
  # Now perform PheWAS-Ttest SENSITIVITY MR and combine for reporting ease
  #~~~~~~~~~~~~~
  mrdat_ttest_sens_filt <- mrdat[mrdat$PheWAS_ttest_sensitivity_keep == "TRUE",]
  # Perform MR on ttest_sens_filtered data
  mrInput_ttest_sens_filt <- MendelianRandomization:: mr_input(bx = mrdat_ttest_sens_filt$beta.exposure,
                                                          bxse = mrdat_ttest_sens_filt$se.exposure,
                                                          by = mrdat_ttest_sens_filt$beta.outcome,
                                                          byse = mrdat_ttest_sens_filt$se.outcome)
  mrRes_ttest_sens_filt <- MendelianRandomization::mr_allmethods(mrInput_ttest_sens_filt, method = "all")
  mrRes_ttest_sens_filt_df <- mrRes_ttest_sens_filt@Values
  
  #~~~~~~~~~~~~~
  # Now perform Steiger filtering and combine for reporting ease
  #~~~~~~~~~~~~~
  mrdat_stgfilt <- mrdat[mrdat$steiger_dir == "TRUE" | mrdat$steiger_pval >= 0.05,]
  # Perform MR on steiger-filtered data
  mrInput_stg_filt <- MendelianRandomization:: mr_input(bx = mrdat_stgfilt$beta.exposure,
                                                          bxse = mrdat_stgfilt$se.exposure,
                                                          by = mrdat_stgfilt$beta.outcome,
                                                          byse = mrdat_stgfilt$se.outcome)
  mrRes_stg_filt <- MendelianRandomization::mr_allmethods(mrInput_stg_filt, method = "all")
  mrRes_stg_filt_df <- mrRes_stg_filt@Values
  
  #~~~~~~~~~~~~~
  # Now perform "Core" MR (no sensitivity here)
  #~~~~~~~~~~~~~
  mrdat_corefilt <- mrdat[(mrdat$steiger_dir == "TRUE" | mrdat$steiger_pval >= 0.05) &
                            (mrdat$PheWAS_ttest_keep == "TRUE") &
                            (mrdat$PheWAS_PVE_keep == "TRUE"),]
  # Perform MR on steiger-filtered data
  mrInput_corefilt <- MendelianRandomization:: mr_input(bx = mrdat_corefilt$beta.exposure,
                                                        bxse = mrdat_corefilt$se.exposure,
                                                        by = mrdat_corefilt$beta.outcome,
                                                        byse = mrdat_corefilt$se.outcome)
  mrRes_corefilt <- MendelianRandomization::mr_allmethods(mrInput_corefilt, method = "all")
  mrRes_corefilt_df <- mrRes_corefilt@Values
  
  #~~~~~~~~~~~~~
  # Make Table - all filtering methods grouped by MR test 
  #~~~~~~~~~~~~~
  # Make table
  # For S4 objects, access slots with @
  
  # Unfiltered results
  mrRes_df$GI <- "all"
  mrRes_df$N_snps <- nrow(mrdat)
  
  # Steiger filtered results
  mrRes_stg_filt_df$GI <- "Steiger Filtered"
  mrRes_stg_filt_df$N_snps <- nrow(mrdat_stgfilt)
  
  # PheWAS T-test results
  mrRes_ttest_filt_df$GI <- "PheWAS T-test Filtered"
  mrRes_ttest_filt_df$N_snps <- nrow(mrdat_ttest_filt)
  
  # PheWAS T-test Sensitivity results
  mrRes_ttest_sens_filt_df$GI <- "PheWAS T-test Filtered - Sensitivity"
  mrRes_ttest_sens_filt_df$N_snps <- nrow(mrdat_ttest_sens_filt)
  
  # PheWAS PVE results
  mrRes_pve_filt_df$GI <- "PheWAS PVE Filtered"
  mrRes_pve_filt_df$N_snps <- nrow(mrdat_pve_filt)
  
  # PheWAS PVE Sensitivity results
  mrRes_pve_sens_filt_df$GI <- "PheWAS PVE Filtered - Sensitivity"
  mrRes_pve_sens_filt_df$N_snps <- nrow(mrdat_pve_sens_filt)
  
  # CORE results
  mrRes_corefilt_df$GI <- "Core GI"
  mrRes_corefilt_df$N_snps <- nrow(mrdat_corefilt)
  
  # Create a data frame that shows the non-filtered and filtered data next to each other for sharing
  # Number of rows in each data frame
  n_rows_all <- nrow(mrRes_df)
  n_rows_stgfilt <- nrow(mrRes_stg_filt_df)
  n_rows_ttestfilt <- nrow(mrRes_ttest_filt_df)
  n_rows_ttest_sens_filt <- nrow(mrRes_ttest_sens_filt_df)
  n_rows_pvefilt <- nrow(mrRes_pve_filt_df)
  n_rows_pve_sens_filt <- nrow(mrRes_pve_sens_filt_df)
  n_rows_corefilt <- nrow(mrRes_corefilt_df)
  
  # Create an empty data frame C
  df_pub <- data.frame(ID = integer(), Value = character())
  
  # Loop through rows of A and B alternatively and combine them into C
  for (j in 1:max(n_rows_all, n_rows_stgfilt, n_rows_ttestfilt, n_rows_pvefilt, n_rows_ttest_sens_filt, n_rows_pve_sens_filt)) {
    if (j <= n_rows_all) {
      df_pub <- rbind(df_pub, mrRes_df[j, ])
    }
    if (j <= n_rows_stgfilt) {
      df_pub <- rbind(df_pub, mrRes_stg_filt_df[j, ])
    }
    if (j <= n_rows_ttestfilt) {
      df_pub <- rbind(df_pub, mrRes_ttest_filt_df[j,])
    }
    if (j <= n_rows_ttest_sens_filt) {
      df_pub <- rbind(df_pub, mrRes_ttest_sens_filt_df[j,])
    }
    if (j <= n_rows_pvefilt) {
      df_pub <- rbind(df_pub, mrRes_pve_filt_df[j,])
    }
    if (j <= n_rows_pve_sens_filt) {
      df_pub <- rbind(df_pub, mrRes_pve_sens_filt_df[j,])
    }
    if (j <= n_rows_corefilt) {
      df_pub <- rbind(df_pub, mrRes_corefilt_df[j,])
    }
  }
  
  # Save the combined table
  combined_table_dir <- "/pl/active/colelab/users/kjames/refinedMR/results/Nonfiltered_GI/tables"
  write.csv(df_pub, paste0(combined_table_dir, "/MR_NonFiltered_and_Filtered_Results_Formatted", EXP_pheno,"_", OUT_pheno, "_sensitivity.csv"))
  
  #~~~~~~~~~~~~~
  # Prep Table - Weighted median method with all trait-pairs and filtering methods
  #~~~~~~~~~~~~~
  mrRes_df_wm <- mrRes_df %>% filter(Method %in% c("Weighted median"))
  mrRes_stg_filt_df_wm <- mrRes_stg_filt_df %>% filter(Method %in% c("Weighted median"))
  mrRes_ttest_filt_df_wm <- mrRes_ttest_filt_df %>% filter(Method %in% c("Weighted median"))
  mrRes_ttest_sens_filt_df_wm <- mrRes_ttest_sens_filt_df %>% filter(Method %in% c("Weighted median"))
  mrRes_pve_filt_df_wm <- mrRes_pve_filt_df %>% filter(Method %in% c("Weighted median"))
  mrRes_pve_sens_filt_df_wm <- mrRes_pve_sens_filt_df %>% filter(Method %in% c("Weighted median"))
  mrRes_core_filt_df_wm <- mrRes_corefilt_df %>% filter(Method %in% c("Weighted median"))
  
  # Combine
  mr_wm_df <- rbind(mrRes_df_wm,
                    mrRes_stg_filt_df_wm,
                    mrRes_ttest_filt_df_wm,
                    mrRes_ttest_sens_filt_df_wm,
                    mrRes_pve_filt_df_wm,
                    mrRes_pve_sens_filt_df_wm,
                    mrRes_core_filt_df_wm)
  # Name exposure/outcome
  mr_wm_df$exposure <- EXP_name
  mr_wm_df$outcome <- OUT_name
  # Rearrange
  mr_wm_df <- mr_wm_df %>% select(exposure, outcome, everything())
  # Save in list slot
  mr_wm_list [[i]] <- mr_wm_df
  
  # Save the combined table
  general_table_dir <- "/pl/active/colelab/users/kjames/refinedMR/results/general/tables"
  write.csv(mr_wm_df, paste0(general_table_dir, "/MR_NonFiltered_and_Filtered_WeightedMedian_Results_", EXP_pheno,"_", OUT_pheno, "_sensitivity.csv"))
}

# Process list into all-trait table
mr_wm_full_table <- do.call("rbind", mr_wm_list)
general_table_dir <- "/pl/active/colelab/users/kjames/refinedMR/results/general/tables"
write.csv(mr_wm_full_table, paste0(general_table_dir, "/MR_NonFiltered_and_Filtered_WeightedMedian_AllTraitPairs_ProcessedManuscriptReady_withSensitivity.csv"))

#~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~
# IVW - Random Effects: 
#~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~
# Due to high heterogeneity and a project all about pleiotropy...
# Run the IVW random effects model for all subsets of GI
ivw_re_list <- list()
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
  #~~~~~~~~~~~~~
  # Non-filtered: Run IVW with random effects due to large amounts of heterogeneity
  #~~~~~~~~~~~~~
  # Define MR object
  mrInput <- MendelianRandomization:: mr_input(bx = mrdat$beta.exposure,
                                               bxse = mrdat$se.exposure,
                                               by = mrdat$beta.outcome,
                                               byse = mrdat$se.outcome)
  # Run MR-IVW-random effect model
  mr_ivw_random <- MendelianRandomization::mr_ivw(mrInput,
                                                  model = "random") #object with filtered data 
  # Create a table of results for saving
  mr_ivw_random_df <- data.frame(test = "mr_ivw_randomeffects",
                                 model = mr_ivw_random@Model,
                                 robust = mr_ivw_random@Robust,
                                 penalized = mr_ivw_random@Penalized,
                                 estimate = mr_ivw_random@Estimate,
                                 std_error = mr_ivw_random@StdError,
                                 cilower = mr_ivw_random@CILower,
                                 ciupper = mr_ivw_random@CIUpper,
                                 pvalue = mr_ivw_random@Pvalue,
                                 heter_Q = mr_ivw_random@Heter.Stat[1],
                                 heter_pvalue = mr_ivw_random@Heter.Stat[2],
                                 fstat = mr_ivw_random@Fstat)
  
  # Save table
  write.csv(mr_ivw_random_df, paste0("/pl/active/colelab/users/kjames/refinedMR/results/Nonfiltered_GI/tables/MR_NonFiltered_GI_IVWres_RANDOMeffects_", EXP_pheno,"_", OUT_pheno, ".csv"))
  
  #~~~~~~~~~~~~~
  # PheWAS-PVE:  Run IVW with random effects due to large amounts of heterogeneity
  #~~~~~~~~~~~~~
  mrdat_pve_filt <- mrdat[mrdat$PheWAS_PVE_keep == "TRUE",]
  
  # Define MR object
  mrInput_pve <- MendelianRandomization:: mr_input(bx = mrdat_pve_filt$beta.exposure,
                                               bxse = mrdat_pve_filt$se.exposure,
                                               by = mrdat_pve_filt$beta.outcome,
                                               byse = mrdat_pve_filt$se.outcome)
  # Run MR-IVW-random effect model
  mr_pve_ivw_random <- MendelianRandomization::mr_ivw(mrInput_pve,
                                                  model = "random") #object with filtered data 
  
  # Create a table of results for saving
  mr_pve_ivw_random_df <- data.frame(test = "mr_pve_ivw_randomeffects",
                                 model = mr_pve_ivw_random@Model,
                                 robust = mr_pve_ivw_random@Robust,
                                 penalized = mr_pve_ivw_random@Penalized,
                                 estimate = mr_pve_ivw_random@Estimate,
                                 std_error = mr_pve_ivw_random@StdError,
                                 cilower = mr_pve_ivw_random@CILower,
                                 ciupper = mr_pve_ivw_random@CIUpper,
                                 pvalue = mr_pve_ivw_random@Pvalue,
                                 heter_Q = mr_pve_ivw_random@Heter.Stat[1],
                                 heter_pvalue = mr_pve_ivw_random@Heter.Stat[2],
                                 fstat = mr_pve_ivw_random@Fstat)
  
  # Save table
  write.csv(mr_pve_ivw_random_df, paste0("/pl/active/colelab/users/kjames/refinedMR/results/PheWAS_PVE/tables/MR_NonFiltered_GI_IVWres_RANDOMeffects_", EXP_pheno,"_", OUT_pheno, ".csv"))
  
  #~~~~~~~~~~~~~
  # Now perform PheWAS-PVE SENSITIVITY MR and combine with non-filtered for reporting ease
  #~~~~~~~~~~~~~
  mrdat_pve_sens_filt <- mrdat[mrdat$PheWAS_PVE_sensitivity_keep == "TRUE",]
  # Perform MR on pve-filtered data
  mrInput_pve_sens_filt <- MendelianRandomization:: mr_input(bx = mrdat_pve_sens_filt$beta.exposure,
                                                             bxse = mrdat_pve_sens_filt$se.exposure,
                                                             by = mrdat_pve_sens_filt$beta.outcome,
                                                             byse = mrdat_pve_sens_filt$se.outcome)
  
  # Run MR-IVW-random effect model
  mr_pve_sens_ivw_random <- MendelianRandomization::mr_ivw(mrInput_pve_sens_filt,
                                                      model = "random") #object with filtered data 
  
  # Create a table of results for saving
  mr_pve_sens_ivw_random_df <- data.frame(test = "mr_pve_sensitivity_ivw_randomeffects",
                                     model = mr_pve_sens_ivw_random@Model,
                                     robust = mr_pve_sens_ivw_random@Robust,
                                     penalized = mr_pve_sens_ivw_random@Penalized,
                                     estimate = mr_pve_sens_ivw_random@Estimate,
                                     std_error = mr_pve_sens_ivw_random@StdError,
                                     cilower = mr_pve_sens_ivw_random@CILower,
                                     ciupper = mr_pve_sens_ivw_random@CIUpper,
                                     pvalue = mr_pve_sens_ivw_random@Pvalue,
                                     heter_Q = mr_pve_sens_ivw_random@Heter.Stat[1],
                                     heter_pvalue = mr_pve_sens_ivw_random@Heter.Stat[2],
                                     fstat = mr_pve_sens_ivw_random@Fstat)
  
  # Save table
  write.csv(mr_pve_sens_ivw_random_df, paste0("/pl/active/colelab/users/kjames/refinedMR/results/PheWAS_PVE/tables/MR_NonFiltered_GI_IVWres_RANDOMeffects_", EXP_pheno,"_", OUT_pheno, "_sensitivity.csv"))
  
  #~~~~~~~~~~~~~
  # PheWAS-Ttest: Run IVW with random effects due to large amounts of heterogeneity
  #~~~~~~~~~~~~~
  mrdat_ttest_filt <- mrdat[mrdat$PheWAS_ttest_keep == "TRUE",]
  # Define MR object
  mrInput_ttest <- MendelianRandomization:: mr_input(bx = mrdat_ttest_filt$beta.exposure,
                                               bxse = mrdat_ttest_filt$se.exposure,
                                               by = mrdat_ttest_filt$beta.outcome,
                                               byse = mrdat_ttest_filt$se.outcome)
  # Run MR-IVW-random effect model
  mr_ttest_ivw_random <- MendelianRandomization::mr_ivw(mrInput_ttest,
                                                  model = "random") #object with filtered data 
  # Create a table of results for saving
  mr_ttest_ivw_random_df <- data.frame(test = "mr_ttest_ivw_randomeffects",
                                 model = mr_ttest_ivw_random@Model,
                                 robust = mr_ttest_ivw_random@Robust,
                                 penalized = mr_ttest_ivw_random@Penalized,
                                 estimate = mr_ttest_ivw_random@Estimate,
                                 std_error = mr_ttest_ivw_random@StdError,
                                 cilower = mr_ttest_ivw_random@CILower,
                                 ciupper = mr_ttest_ivw_random@CIUpper,
                                 pvalue = mr_ttest_ivw_random@Pvalue,
                                 heter_Q = mr_ttest_ivw_random@Heter.Stat[1],
                                 heter_pvalue = mr_ttest_ivw_random@Heter.Stat[2],
                                 fstat = mr_ttest_ivw_random@Fstat)
  
  # Save table
  write.csv(mr_ttest_ivw_random_df, paste0("/pl/active/colelab/users/kjames/refinedMR/results/PheWAS_Cluster/tables/MR_NonFiltered_GI_IVWres_RANDOMeffects_", EXP_pheno,"_", OUT_pheno, ".csv"))
  
  #~~~~~~~~~~~~~
  # Now perform PheWAS-Ttest SENSITIVITY MR and combine for reporting ease
  #~~~~~~~~~~~~~
  mrdat_ttest_sens_filt <- mrdat[mrdat$PheWAS_ttest_sensitivity_keep == "TRUE",]
  # Perform MR on ttest_sens_filtered data
  mrInput_ttest_sens_filt <- MendelianRandomization:: mr_input(bx = mrdat_ttest_sens_filt$beta.exposure,
                                                               bxse = mrdat_ttest_sens_filt$se.exposure,
                                                               by = mrdat_ttest_sens_filt$beta.outcome,
                                                               byse = mrdat_ttest_sens_filt$se.outcome)
  
  # Run MR-IVW-random effect model
  mr_ttest_sens_ivw_random <- MendelianRandomization::mr_ivw(mrInput_ttest_sens_filt,
                                                        model = "random") #object with filtered data 
  # Create a table of results for saving
  mr_ttest_sens_ivw_random_df <- data.frame(test = "mr_ttest_sensitivity_ivw_randomeffects",
                                       model = mr_ttest_sens_ivw_random@Model,
                                       robust = mr_ttest_sens_ivw_random@Robust,
                                       penalized = mr_ttest_sens_ivw_random@Penalized,
                                       estimate = mr_ttest_sens_ivw_random@Estimate,
                                       std_error = mr_ttest_sens_ivw_random@StdError,
                                       cilower = mr_ttest_sens_ivw_random@CILower,
                                       ciupper = mr_ttest_sens_ivw_random@CIUpper,
                                       pvalue = mr_ttest_sens_ivw_random@Pvalue,
                                       heter_Q = mr_ttest_sens_ivw_random@Heter.Stat[1],
                                       heter_pvalue = mr_ttest_sens_ivw_random@Heter.Stat[2],
                                       fstat = mr_ttest_sens_ivw_random@Fstat)
  
  # Save table
  write.csv(mr_ttest_sens_ivw_random_df, paste0("/pl/active/colelab/users/kjames/refinedMR/results/PheWAS_Cluster/tables/MR_NonFiltered_GI_IVWres_RANDOMeffects_", EXP_pheno,"_", OUT_pheno, "_sensitivity.csv"))
  
  
  #~~~~~~~~~~~~~
  # Steiger filtering: Run IVW with random effects due to large amounts of heterogeneity
  #~~~~~~~~~~~~~
  mrdat_stgfilt <- mrdat[mrdat$steiger_dir == "TRUE" | mrdat$steiger_pval >= 0.05,]
  # Define MR object
  mrInput_stg <- MendelianRandomization:: mr_input(bx = mrdat_stgfilt$beta.exposure,
                                               bxse = mrdat_stgfilt$se.exposure,
                                               by = mrdat_stgfilt$beta.outcome,
                                               byse = mrdat_stgfilt$se.outcome)
  # Run MR-IVW-random effect model
  mr_stg_ivw_random <- MendelianRandomization::mr_ivw(mrInput_stg,
                                                  model = "random") #object with filtered data 
  # Create a table of results for saving
  mr_stg_ivw_random_df <- data.frame(test = "mr_stg_ivw_randomeffects",
                                 model = mr_stg_ivw_random@Model,
                                 robust = mr_stg_ivw_random@Robust,
                                 penalized = mr_stg_ivw_random@Penalized,
                                 estimate = mr_stg_ivw_random@Estimate,
                                 std_error = mr_stg_ivw_random@StdError,
                                 cilower = mr_stg_ivw_random@CILower,
                                 ciupper = mr_stg_ivw_random@CIUpper,
                                 pvalue = mr_stg_ivw_random@Pvalue,
                                 heter_Q = mr_stg_ivw_random@Heter.Stat[1],
                                 heter_pvalue = mr_stg_ivw_random@Heter.Stat[2],
                                 fstat = mr_stg_ivw_random@Fstat)
  
  # Save table
  write.csv(mr_stg_ivw_random_df, paste0("/pl/active/colelab/users/kjames/refinedMR/results/steiger/tables/MR_NonFiltered_GI_IVWres_RANDOMeffects_", EXP_pheno,"_", OUT_pheno, ".csv"))
  
  #~~~~~~~~~~~~~
  # "Core" MR: Run IVW with random effects due to large amounts of heterogeneity; no sensitivity here
  #~~~~~~~~~~~~~
  mrdat_corefilt <- mrdat[(mrdat$steiger_dir == "TRUE" | mrdat$steiger_pval >= 0.05) &
                            (mrdat$PheWAS_ttest_keep == "TRUE") &
                            (mrdat$PheWAS_PVE_keep == "TRUE"),]

  # Define MR object
  mrInput_core <- MendelianRandomization:: mr_input(bx = mrdat_corefilt$beta.exposure,
                                               bxse = mrdat_corefilt$se.exposure,
                                               by = mrdat_corefilt$beta.outcome,
                                               byse = mrdat_corefilt$se.outcome)
  # Run MR-IVW-random effect model
  mr_core_ivw_random <- MendelianRandomization::mr_ivw(mrInput_core,
                                                  model = "random") #object with filtered data 
  # Create a table of results for saving
  mr_core_ivw_random_df <- data.frame(test = "mr_core_ivw_randomeffects",
                                 model = mr_core_ivw_random@Model,
                                 robust = mr_core_ivw_random@Robust,
                                 penalized = mr_core_ivw_random@Penalized,
                                 estimate = mr_core_ivw_random@Estimate,
                                 std_error = mr_core_ivw_random@StdError,
                                 cilower = mr_core_ivw_random@CILower,
                                 ciupper = mr_core_ivw_random@CIUpper,
                                 pvalue = mr_core_ivw_random@Pvalue,
                                 heter_Q = mr_core_ivw_random@Heter.Stat[1],
                                 heter_pvalue = mr_core_ivw_random@Heter.Stat[2],
                                 fstat = mr_core_ivw_random@Fstat)
  
  # Save table
  write.csv(mr_core_ivw_random_df, paste0("/pl/active/colelab/users/kjames/refinedMR/results/core/tables/MR_NonFiltered_GI_IVWres_RANDOMeffects_", EXP_pheno,"_", OUT_pheno, ".csv"))
  
  #~~~~~~~~~~~~~
  # Make Table
  #~~~~~~~~~~~~~
  mr_ivw_randomeffects_df <- rbind(mr_ivw_random_df,
                                   mr_stg_ivw_random_df,
                                   mr_ttest_ivw_random_df,
                                   mr_ttest_sens_ivw_random_df,
                                   mr_pve_ivw_random_df,
                                   mr_pve_sens_ivw_random_df,
                                   mr_core_ivw_random_df)
  

  # Save the combined table
  general_table_dir <- "/pl/active/colelab/users/kjames/refinedMR/results/general/tables"
  write.csv(mr_ivw_randomeffects_df, paste0(general_table_dir, "/MR_NonFiltered_and_Filtered_IVW_RandomEffects_Results_", EXP_pheno,"_", OUT_pheno, "_withSensitivity.csv"))

  #~~~~~~~~~~~~~
  # Make combined trait pair list for table generation outside of loop
  #~~~~~~~~~~~~~
  mr_ivw_randomeffects_df$exposure <- EXP_name
  mr_ivw_randomeffects_df$outcome <- OUT_name
  ivw_re_list[[i]] <- mr_ivw_randomeffects_df
  
}

# Get table
ivw_re_df <- do.call("rbind",ivw_re_list)
# Save whole table
write.csv(ivw_re_df, paste0(general_table_dir, "/MR_NonFiltered_and_Filtered_IVW_RandomEffects_Results_AllTraitPairs_withSensitivity.csv"))

# Process table
# Edit for easy input into manuscript
ivw_re_df_processed <- ivw_re_df
# For columns with numbers, only keep 3 significant figures
ivw_re_df_processed[sapply(ivw_re_df, is.numeric)] <- lapply(ivw_re_df[sapply(ivw_re_df, is.numeric)], signif, digits = 3)

# Rename test column with manuscript labels
ivw_re_df_processed$test[ivw_re_df_processed$test == "mr_core_ivw_randomeffects"] <- "Core"
ivw_re_df_processed$test[ivw_re_df_processed$test == "mr_ivw_randomeffects"] <- "All"
ivw_re_df_processed$test[ivw_re_df_processed$test == "mr_pve_ivw_randomeffects"] <- "PheWAS-PVE"
ivw_re_df_processed$test[ivw_re_df_processed$test == "mr_pve_sensitivity_ivw_randomeffects"] <- "PheWAS-PVE Sensitivity (no diet)"
ivw_re_df_processed$test[ivw_re_df_processed$test == "mr_stg_ivw_randomeffects"] <- "Steiger"
ivw_re_df_processed$test[ivw_re_df_processed$test == "mr_ttest_ivw_randomeffects"] <- "PheWAS-TT"
ivw_re_df_processed$test[ivw_re_df_processed$test == "mr_ttest_sensitivity_ivw_randomeffects"] <- "PheWAS-TT Sensitivity (no diet)"

# Rearrange/select columns for manuscript
ivw_re_df_processed_2 <- ivw_re_df_processed %>% select(exposure, outcome, test, estimate,
                                                        std_error, cilower, ciupper, pvalue,
                                                        heter_Q, heter_pvalue, fstat)

# Save processed table
write.csv(ivw_re_df_processed_2, paste0(general_table_dir, "/MR_NonFiltered_and_Filtered_IVW_RandomEffects_Results_AllTraitPairs_ProcessedManuscriptReady_withSensitivity.csv"))

