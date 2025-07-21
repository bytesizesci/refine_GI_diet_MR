# Harmonization function
# Author: KJS
# Date: July 31, 2024


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

harmonize_alleles <- function(gi) {
  if (all(gi$`effect_allele.exposure` == gi$`effect_allele.outcome`)) {
    print("action=1")
    action = 1
  } else {
    print("action=2/3")
    aligned <- which(gi$`effect_allele.exposure` == gi$`effect_allele.outcome` &
                       gi$`other_allele.exposure` == gi$`other_allele.outcome`)
    swapped <- which(gi$`effect_allele.exposure` == gi$`other_allele.outcome` &
                       gi$`other_allele.exposure` == gi$`effect_allele.outcome`)
    
    # Flip the effect sizes for the swapped alleles
    gi[swapped, 'beta.exposure'] <- gi[swapped, 'beta.exposure'] * -1
    
    # Swap the allele labels for consistency
    temp <- gi[swapped, 'effect_allele.exposure']
    gi[swapped, 'effect_allele.exposure'] <- gi[swapped, 'other_allele.exposure']
    gi[swapped, 'other_allele.exposure'] <- temp
    
    # New - Swap the allele frequency labels for consistency
    temp <- gi[swapped, 'eaf.exposure']
    gi[swapped, 'eaf.exposure'] <- gi[swapped, 'oaf.exposure']
    gi[swapped, 'oaf.exposure'] <- temp
    
    # New - Update the variant name column once the changes are made
    gi$variant <- paste0(gi$chr.exposure, ":", gi$pos.exposure, ":", gi$other_allele.exposure, ":", gi$effect_allele.exposure)
    
    # Keep the aligned and swapped SNPs
    #gi <- gi[c(aligned, swapped), ]
    action = 1  # Made sure all strands are okay
  }
  return(list(data = gi, action = action))
}
