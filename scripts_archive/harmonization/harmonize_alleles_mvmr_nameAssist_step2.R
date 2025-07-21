# Harmonization assisting function for MVMR cols
# Author: KJS
# Date: December 30, 2024


# Needs full data with no NAs
# Requires columns named
# effect_allele.exposure, other_allele.exposure
# effect_allele.outcome, other_allele.outcome
# beta.exposure
#
# This chunk adapted from Clustering.R by LD
# Checks that the effect_allele is the same in the exposure and outcome
# If it isn't, it multiplies beta by -1 and flips the alleles
# Then it only keeps variants that were originally aligned or originally identified as 
# swapped (and were fixed). If three alleles are present the SNP is booted.

rename_and_execute_step2 <- function(data, func) {
  # Rename columns a, b, c, d to e, f, g, h
  original_names <- c("effect_allele.exposure.x", "other_allele.exposure.x", "effect_allele.outcome", "other_allele.outcome", "beta.exposure.x","oaf.exposure.x","eaf.exposure.x", "chr.exposure.x", "pos.exposure.x" )
  new_names <- c("effect_allele.exposure", "other_allele.exposure", "effect_allele.outcome", "other_allele.outcome", "beta.exposure","oaf.exposure", "eaf.exposure", "chr.exposure","pos.exposure" )
  
  # Ensure columns exist
  if (!all(original_names %in% colnames(data))) {
    stop("Not all required columns (effect_allele.exposure.x, other_allele.exposure.x, effect_allele.exposure.y, other_allele.exposure.y, beta.exposure.x) are present in the data.")
  }
  
  # Rename columns
  colnames(data)[match(original_names, colnames(data))] <- new_names
  
  # Execute the second function #harmonize_alleles
  #data <- harmonize_alleles(data) # for testing
  data <- func(data) #gi_exp_complete
  data <- data$data
  
  # Rename columns back to original names
  colnames(data)[match(new_names, colnames(data))] <- original_names
  
  return(data)
}
