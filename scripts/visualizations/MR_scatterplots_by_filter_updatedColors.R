# Goal: Develop visualizations for manuscript 
# Author: KJS
# Date: November 13, 2024

# Load ggplot2 library
library(ggplot2)
library(viridis)
library(dplyr)
library(RColorBrewer)

# Set wd
wdir <- "/pl/active/colelab/users/kjames/refinedMR/"
setwd(wdir)

# Load harmonized and standardized genetic instruments
mrobj <- readRDS("interim_data/merged_GI/MRobject_steiger_PheWAStt_PheWASPVE_241120.rds")

# Get colors for plots
colors <- viridis(12)
print(colors)

# Get updated colors for plots
paired_colors <- brewer.pal(12, "Paired")
print(paired_colors)
colors <- c("#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C", "#FDBF6F", "#FF7F00", "#CAB2D6", "#6A3D9A", "#FFFF99", "#B15928")

# Create a data frame with naming conventions for loop
df_names <- data.frame(EXP_pheno = c("oilyfish", "oilyfish", "oilyfish", "whiteXwholegrain", "whiteXwholegrain",  "whiteXwholegrain", "alcohol", "alcohol","alcohol"),
                       OUT_pheno = c("GCST90239664_TG_Graham", "GCST90132314_CVD_Aragam", "GCST006901_Height_Yengo_2018", "GCST90239658_LDL_Graham","GCST90132314_CVD_Aragam", "GCST006901_Height_Yengo_2018", "GCST90013405_ALT_Pazoki", "GCST90319877_cir_Ghouse", "GCST006901_Height_Yengo_2018"),
                       EXP_names = c("Oily Fish", "Oily fish", "Oily fish", "White vs Whole Grain Bread", "White vs Whole Grain Bread", "White vs Whole Grain Bread","Alcohol", "Alcohol", "Alcohol"),
                       OUT_names = c("Triglycerides", "Cardiovascular Disease", "Height", "LDL Cholesterol", "Cardiovascular Disease", "Height", "Alanine Aminotransferase", "Liver Cirrhosis", "Height"))

# Test with first slot, comment out when not testing parts of the loop
#mrdat <- mrobj[[1]]

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

  # Make mrstg column
  mrdat$stgfilt_keep <- ifelse(mrdat$steiger_dir == "TRUE" & mrdat$steiger_pval <= 0.05, TRUE, FALSE )
  
  # Make a column for core SNPs kept with all methods
  mrdat$allfilt_keep <- ifelse(mrdat$stgfilt_keep == TRUE & mrdat$PheWAS_PVE_keep == TRUE & mrdat$PheWAS_ttest_keep == TRUE, TRUE, FALSE)
  
  # Create new columns for the graph
  mrdat$stg_plot <- ifelse(mrdat$stgfilt_keep == TRUE, "Steiger-Kept", "Steiger-Filtered")
  mrdat$ttest_plot <- ifelse(mrdat$PheWAS_ttest_keep == TRUE, "Ttest-Kept", "Ttest-Filtered")
  mrdat$pve_plot <- ifelse(mrdat$PheWAS_PVE_keep == TRUE, "PVE-Kept", "PVE-Filtered")
  
  # Create a new column that combines approach_1 and approach_2
  mrdat$color_group <- with(mrdat, paste(stg_plot, ttest_plot, pve_plot, sep = ", "))
  table(mrdat$color_group)
  
  # Plot the scatterplot
  prettyplot <- ggplot(mrdat, aes(x = beta.exposure, y = beta.outcome)) +
    # Add points
    geom_point(aes(color = color_group), size = 2) +
    # Add a trend line for each GI
    geom_smooth(
      data = subset(mrdat, allfilt_keep == "TRUE"),
      method = "lm",
      se = TRUE,
      aes(#color = "Core", 
          linetype = "Core"),
      alpha = 0.2,
      color = "black"
    ) +
    geom_smooth(
      data = subset(mrdat, stgfilt_keep == "TRUE"),
      method = "lm",
      aes(#color = "Steiger",
          linetype = "Steiger"),
      se = TRUE,
      alpha = 0.2,
      color = "black"
    ) +
    geom_smooth(
      data = subset(mrdat, PheWAS_ttest_keep == "TRUE"),
      method = "lm",
      se = TRUE,
      aes(#color = "Ttest",
          linetype = "Ttest"),
      alpha = 0.2,
      color = "black"
    ) +
    geom_smooth(
      data = subset(mrdat, PheWAS_PVE_keep == "TRUE"),
      method = "lm",
      se = TRUE,
      aes(#color = "PVE",
          linetype = "PVE"),
      alpha = 0.2,
      color = "black"
    ) +
    geom_smooth(
      data = mrdat, #all data
      method = "lm",
      se = TRUE,
      aes(#color = "No filtering",
          linetype = "No filtering"),
      alpha = 0.2,
      color = "black"
    ) +
    # Add manual colors for points based on inclusion/exclusion across methods
    #"#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C", "#FDBF6F", "#FF7F00", "#CAB2D6", "#6A3D9A", "#FFFF99", "#B15928"
    scale_color_manual(values = c(
      "Steiger-Kept, Ttest-Kept, PVE-Kept" = "#A6CEE3", #"#440154FF",
      "Steiger-Kept, Ttest-Kept, PVE-Filtered" = "#1F78B4", #"#26828EFF",
      "Steiger-Kept, Ttest-Filtered, PVE-Kept" = "#B2DF8A", #"#6DCD59FF",
      "Steiger-Kept, Ttest-Filtered, PVE-Filtered" = "#33A02C", #"#3E4A89FF",
      "Steiger-Filtered, Ttest-Kept, PVE-Filtered" = "#FB9A99", #"#35B779FF",
      "Steiger-Filtered, Ttest-Kept, PVE-Kept" = "#E31A1C", #"#B4DE2CFF",
      "Steiger-Filtered, Ttest-Filtered, PVE-Kept" = "#FDBF6F" #"#FDE725FF"
      #"Steiger-Filtered, Ttest-Filtered, PVE-Filtered"= "#FF7F00" #this group (all false- won't occur because it would get removed!)
    )) +
    # Set linetype legend explicitly
    scale_linetype_manual(values = c(
      "Core" = "solid",
      "Steiger" = "dashed",
      "Ttest" = "dotted",
      "PVE" = "dotdash",
      "No filtering" = "longdash"
    )) +
    # Legend and axis info
    labs(
      #title = "MR Scatterplot",
      x = paste0("Genetic Liability of ", EXP_name),
      y = paste0("Genetic Liability of ", OUT_name),
      color = "Filtering",
      linetype = "Method"
    ) +
    theme_bw() +
    # Make the linetype sample bigger in the legend
    guides(
      linetype = guide_legend(
        keywidth = 3,  # Wider linetype sample
        keyheight = 0.5
      ))
  prettyplot
  # Save
  # USER INPUT: Check paths and variable names
  plots_dir <- "/pl/active/colelab/users/kjames/refinedMR/results/general/plots/"
  
  # Save
  ggsave(paste0(plots_dir,"/MR_",EXP_pheno, "_", OUT_pheno,"_scatterplot_allmethods_inclusion_exclusion_updatedColors.png"),
         plot = prettyplot,
         device = NULL,
         path = NULL,
         scale = 1,
         width = 6,
         height = 4,
         units = c("in"),
         dpi = 300,
         limitsize = TRUE)
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Simplify the plots
# Scatter plots with GWS IVW and PheWAS-Ttest IVW
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~




