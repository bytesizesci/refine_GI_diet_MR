# Goal: Forest plot code for MR-diet project
# Author: KJS
# Date: 11/27/2024

# Set wd
wdir <- "/pl/active/colelab/users/kjames/refinedMR/"
setwd(wdir)

# Load libraries
library(dplyr)
library(ggplot2)
library(MendelianRandomization)

# Load IVW results
data <- read.csv("results/general/tables/MR_NonFiltered_and_Filtered_IVW_RandomEffects_Results_AllTraitPairs.csv", row.names = 1)

# Make label for traits
data$exp_out <- paste0(data$exposure, " \u2192 ", data$outcome) #\u2192 prints an arrow
# Order exp_out 
str(data$exp_out)
desired_order <- c("Alcohol → Height",
                   "Alcohol → Liver Cirrhosis",
                   "Alcohol → Alanine Aminotransferase",
                   "White vs Whole Grain Bread → Height",
                   "White vs Whole Grain Bread → Cardiovascular Disease",
                   "White vs Whole Grain Bread → LDL Cholesterol",
                   "Oily fish → Height",
                   "Oily fish → Cardiovascular Disease",
                   "Oily Fish → Triglycerides")
data$exp_out <- factor(data$exp_out, levels = desired_order)

# Edit for plot
# Rename test column with manuscript labels
data$test[data$test == "mr_core_ivw_randomeffects"] <- "Core"
data$test[data$test == "mr_ivw_randomeffects"] <- "All"
data$test[data$test == "mr_pve_ivw_randomeffects"] <- "PheWAS-PVE"
data$test[data$test == "mr_stg_ivw_randomeffects"] <- "Steiger"
data$test[data$test == "mr_ttest_ivw_randomeffects"] <- "PheWAS-TT"

# Add colors
my_colors <- c("#d73027","#fc8d59", "#fee090", "#e0f3f8", "#91bfdb", "#4575b4")

# Create the forest plot
forestplot <- ggplot(data, aes(x = estimate, y = exp_out, color = test)) +
  geom_point(position = position_dodge(width = 0.5), size = 3) + # Points for estimates
  geom_errorbar(aes(xmin = cilower, xmax = ciupper),
                position = position_dodge(width = 0.5), width = 0.3) + # Error bars
  scale_color_manual(values = my_colors) +
  #scale_color_brewer(palette = "Spectral") + # Color palette for methods
  labs(
    #title = "Forest Plot of Methods per Category",
    x = "Estimate",
    y = "Category",
    color = "Method"
  ) +
  theme_bw() +
  theme(
    axis.text.y = element_text(size = 12),
    legend.position = "bottom"
  )

forestplot

# Save
plots_dir <- "/pl/active/colelab/users/kjames/refinedMR/results/general/plots/"
ggsave(paste0(plots_dir,"/MR_IVW_forestplot_alltraits_allmethods.png"),
       plot = forestplot,
       device = NULL,
       path = NULL,
       scale = 1,
       width = 10,
       height = 4,
       units = c("in"),
       dpi = 300,
       limitsize = TRUE)
