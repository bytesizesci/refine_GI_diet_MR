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
library(RColorBrewer)
library(cowplot)

# Get updated colors for plots
paired_colors <- brewer.pal(12, "Paired")
print(paired_colors)
colors <- c("#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C", "#FDBF6F", "#FF7F00", "#CAB2D6", "#6A3D9A", "#FFFF99", "#B15928")

# Load IVW results
data <- read.csv("results/general/tables/MR_NonFiltered_and_Filtered_IVW_RandomEffects_Results_AllTraitPairs.csv", row.names = 1)

# Make label for traits
data$exp_out <- paste0(data$exposure, " \u2192 ", data$outcome) #\u2192 prints an arrow

# Factor exp_out 
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

# Check
levels(data$exp_out)

# Now shorten them for plotting purposes
desired_order_short <- c("Alcohol → Height",
                   "Alcohol → Liver Cirrhosis",
                   "Alcohol → ALT",
                   "White vs Whole Grain Bread → Height",
                   "White vs Whole Grain Bread → CVD",
                   "White vs Whole Grain Bread → LDL-C",
                   "Oily fish → Height",
                   "Oily fish → CVD",
                   "Oily Fish → TG")
levels(data$exp_out) <- desired_order_short

# Edit for plot
# Rename test column with manuscript labels
data$test[data$test == "mr_core_ivw_randomeffects"] <- "Core"
data$test[data$test == "mr_ivw_randomeffects"] <- "All"
data$test[data$test == "mr_pve_ivw_randomeffects"] <- "PheWAS-PVE"
data$test[data$test == "mr_stg_ivw_randomeffects"] <- "Steiger"
data$test[data$test == "mr_ttest_ivw_randomeffects"] <- "PheWAS-TT"

# Add colors
#my_colors <- c("#d73027","#fc8d59", "#fee090", "#e0f3f8", "#91bfdb", "#4575b4")
my_colors <- c("#33A02C","#FB9A99", "#E31A1C", "#A6CEE3", "#1F78B4")

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
  geom_vline(xintercept = 0, linetype = "dashed", color = "darkgrey", linewidth = 1) +
  theme_bw() +
  theme(
    axis.text.y = element_text(size = 12),
    legend.position = "bottom"
  )

forestplot

# Save
plots_dir <- "/pl/active/colelab/users/kjames/refinedMR/results/general/plots/"
ggsave(paste0(plots_dir,"/MR_IVW_forestplot_alltraits_allmethods_updatedColors.png"),
       plot = forestplot,
       device = NULL,
       path = NULL,
       scale = 1,
       width = 10,
       height = 4,
       units = c("in"),
       dpi = 300,
       limitsize = TRUE)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Updated way to show the same data
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
pos_df <- data %>% filter(outcome != "Height")
neg_df <- data %>% filter(outcome == "Height")

# Create the forest plot for positive controls
forestplot_pos <- ggplot(pos_df, aes(x = estimate, y = exp_out, color = test)) +
  geom_point(position = position_dodge(width = 0.5), size = 3) + # Points for estimates
  geom_errorbar(aes(xmin = cilower, xmax = ciupper),
                position = position_dodge(width = 0.5), width = 0.3) + # Error bars
  scale_color_manual(values = my_colors) +
  #scale_color_brewer(palette = "Spectral") + # Color palette for methods
  labs(
    #title = "Forest Plot of Methods per Category",
    x = "Estimate",
    y = "",
    color = "Method"
  ) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "darkgrey", linewidth = 1) +
  theme_bw() +
  theme(
    axis.text.y = element_text(size = 14),
    axis.text.x = element_text(size = 14),
    legend.position = "bottom",
    legend.justification = "left"
  )

forestplot_pos

# Create the forest plot for negative controls
forestplot_neg <- ggplot(neg_df, aes(x = estimate, y = exp_out, color = test)) +
  geom_point(position = position_dodge(width = 0.5), size = 3) + # Points for estimates
  geom_errorbar(aes(xmin = cilower, xmax = ciupper),
                position = position_dodge(width = 0.5), width = 0.3) + # Error bars
  scale_color_manual(values = my_colors) +
  #scale_color_brewer(palette = "Spectral") + # Color palette for methods
  labs(
    #title = "Forest Plot of Methods per Category",
    x = "Estimate",
    y = "",
    color = "Method"
  ) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "darkgrey", linewidth = 1) +
  theme_bw() +
  theme(
    axis.text.y = element_text(size = 14),
    axis.text.x = element_text(size = 14),
    legend.position = "bottom",
    legend.justification = "left"
  )

forestplot_neg

# Arrange with cowplot
# Arrange plots
combined_plot <- cowplot::plot_grid(
  forestplot_pos, forestplot_neg,
  ncol = 1,
  align = "v", # horizontal alignment (v=stacked)
  axis = "l", # align the left (y-axis) of the plots
  labels = c("A", "B"),  # Auto-label plots
  label_size = 16
)

combined_plot

# Save
plots_dir <- "/pl/active/colelab/users/kjames/refinedMR/results/general/plots/"
ggsave(paste0(plots_dir,"/MR_IVW_forestplot_posnegFacet_allmethods_updatedColors.png"),
       plot = combined_plot,
       device = NULL,
       path = NULL,
       scale = 1,
       width = 10,
       height = 10,
       units = c("in"),
       dpi = 300,
       limitsize = TRUE)
