# Author: Kristen Sutton
# Version 2: Update to original is merges on chr:position data excluding
# effect/other allele in name. Check the direction using another function.

# Inputs:
# expos_df = columns in mrdat that refer to the exposure
# variant_ID_eff_alt = variant ID coded with alleles as effect_alternative
# phewas_sub = phewas data frame subset to an exposure's genetic variants
# ex: phewas_sub <- pwas_all %>% filter(variant %in% variant_ID_eff_alt)


# For a given SNP, see if its PVE value is higher than the PVE of the trait of interest. 
phewas_thresh_test <- function(variant_ID_eff_alt, expos_df, phewas_sub) {
  # Place for results
  phewas_thresh_results <- list()
  # Loop
  for(i in 1:length(variant_ID_eff_alt)) {
    # Get variant id
    variantID <- variant_ID_eff_alt[i]
    
    # Get snp threshold for that variant from the exposure data
    expos_df_variant <- expos_df[expos_df$chr_pos == variantID,]
    
    # Now, get PVE threshold 
    thresh <- expos_df_variant$PVE.exposure
    
    # Subset PheWAS to that variant
    phewas_variant <- phewas_sub[phewas_sub$chr_pos == variantID,]
    
    # Check if the subset operation resulted in an empty data frame
    if(nrow(phewas_variant) == 0) {
      # If empty, store a placeholder value in the results list
      phewas_thresh_results[[i]] <- NULL
      next  # Skip further processing for this iteration
    }
    
    # Do other traits explain more PVE?
    # 0=fruit explains more PVE, keep SNP
    # 1=other trait explains more PVE, remove/question SNP
    phewas_variant$greater <- ifelse(phewas_variant$PVE <= thresh, 0, 1)
    
    # Subset cols
    phewas_res <- phewas_variant[,c("variant","chr_pos","PVE", "greater", "pval", "minor_allele", "minor_AF", "beta", "se", "tstat", "trait", "file_name")]
    
    # Only keep row if "greater" = 1 (PVE of other trait is greater than PVE food)
    phewas_res <- phewas_res[(phewas_res$greater == 1 & phewas_res$pval < 5E-6),] # Can adjust
    
    # Store results
    phewas_thresh_results[[i]] <- phewas_res
  }
  return(phewas_thresh_results)
}