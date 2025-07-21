#' Estimate r-square of each association
#'
#' Can be applied to exposure_dat, outcome_dat or harmonised_data.
#' Note that it will be beneficial in some circumstances to add the meta data to
#' the data object using [add_metadata()] before running this function. 
#' Also adds effective sample size for case control data.
#'
#' @param dat exposure_dat, outcome_dat or harmonised_data
#'
#' @export
#' @return data frame
add_rsq <- function(dat)
{
  if("id.exposure" %in% names(dat))
  {
    dat <- 	plyr::ddply(dat, c("id.exposure"), function(x)
    {
      add_rsq_one(x, "exposure")
    })
  }
  if("id.outcome" %in% names(dat))
  {
    dat <- 	plyr::ddply(dat, c("id.outcome"), function(x)
    {
      add_rsq_one(x, "outcome")
    })
  }
  return(dat)
}

add_rsq_one <- function(dat, what="exposure")
{
  print("Learning about the add_rsq_one function")
  print(paste("Unique values in", what, "column:", length(unique(dat[[what]]))))
  print(paste("Unique values in units.", what, "column:", length(unique(dat[[paste0("units.", what)]]))))
  
  if(! paste0("units.", what) %in% names(dat))
  {
    dat[[paste0("units.", what)]] <- NA
  }
  stopifnot(length(unique(dat[[what]])) == 1)
  stopifnot(length(unique(dat[[paste0("units.", what)]])) == 1)
  
  if(! paste0("rsq.", what) %in% names(dat))
  {
    dat[[paste0("pval.", what)]][dat[[paste0("pval.", what)]] < 1e-300] <- 1e-300
    if(compareNA(dat[[paste0("units.", what)]][1], "log odds"))
    {
      # message("Estimating rsq.exposure for binary trait")
      # message("Ensure that beta.exposure, eaf.exposure, ncase.exposure, ncontrol.exposure are all specified with no missing values")
      if(! paste0("prevalence.", what) %in% names(dat))
      {
        dat[[paste0("prevalence.", what)]] <- 0.1
        warning(paste0("Assuming ", what, " prevalence of 0.1. Alternatively, add prevalence.", what, " column and re-run."))
      }
      ind1 <- !is.na(dat[[paste0("beta.", what)]]) &
        !is.na(dat[[paste0("eaf.", what)]]) &
        !is.na(dat[[paste0("ncase.", what)]]) &
        !is.na(dat[[paste0("ncontrol.", what)]]) &
        !is.na(dat[[paste0("prevalence.", what)]])
      dat[[paste0("rsq.", what)]] <- NA
      if(sum(ind1) > 0)
      {
        dat[[paste0("rsq.", what)]][ind1] <- get_r_from_lor(
          dat[[paste0("beta.", what)]][ind1],
          dat[[paste0("eaf.", what)]][ind1],
          dat[[paste0("ncase.", what)]][ind1],
          dat[[paste0("ncontrol.", what)]][ind1],
          dat[[paste0("prevalence.", what)]]
        )^2
        dat[[paste0("effective_n.", what)]][ind1] <- effective_n(dat[[paste0("ncase.", what)]][ind1], dat[[paste0("ncontrol.", what)]][ind1])
      } else {
        message("Try adding metadata with add_metadata()")
      }
    } else if(all(grepl("SD", dat[[paste0("units.", what)]])) & all(!is.na(dat[[paste0("eaf.", what)]]))) {
      dat[[paste0("rsq.", what)]] <- NA
      dat[[paste0("rsq.", what)]] <- 2 * dat[[paste0("beta.", what)]]^2 * dat[[paste0("eaf.", what)]] * (1-dat[[paste0("eaf.", what)]])
      dat[[paste0("effective_n.", what)]] <- dat[[paste0("samplesize.", what)]]
    } else {
      ind1 <- !is.na(dat[[paste0("pval.", what)]]) & !is.na(dat[[paste0("samplesize.", what)]])
      dat[[paste0("rsq.", what)]] <- NA
      if(sum(ind1) > 0)
      {		
        dat[[paste0("rsq.", what)]][ind1] <- get_r_from_bsen(
          dat[[paste0("beta.", what)]][ind1],
          dat[[paste0("se.", what)]][ind1],
          dat[[paste0("samplesize.", what)]][ind1]
        )^2
        dat[[paste0("effective_n.", what)]] <- dat[[paste0("samplesize.", what)]]
      } else {
        message("Try adding metadata with add_metadata()")
      }
    }
  }
  return(dat)
}

#' Estimate R-squared from beta, standard error and sample size
#'
#' @param b Array of effect sizes
#' @param se Array of standard errors
#' @param n Array of (effective) sample sizes
#'
#' @export
#' @return Vector of signed r values
get_r_from_bsen <- function(b, se, n)
{
  Fval <- (b/se)^2
  R2 <- Fval / (n-2+Fval)
  return(sqrt(R2) * sign(b))
}

compareNA <- function(v1,v2) {
  same <- (v1 == v2) | (is.na(v1) & is.na(v2))
  same[is.na(same)] <- FALSE
  return(same)
}

# Source: KJ
# Function to check that the naming columns are present
assign_outcome_col <- function(dat, OUT_name, OUT_pheno) {
  # Check if 'outcome' column exists
  if(!"outcome" %in% names(dat)) {
    # If 'outcome' does not exist, create it and assign OUT_name to it
    dat$outcome <- OUT_name
  }
  if(!"id.outcome" %in% names(dat)) {
    # If "id.outcome" does not exist, create it and assign OUT_pheno to it
    dat$id.outcome <- OUT_pheno
  }
  return(dat)
}

assign_exposure_col <- function(dat, EXP_name, EXP_pheno) {
  # Check if 'exposure' column exists
  if(!"exposure" %in% names(dat)) {
    # If 'exposure' does not exist, create it and assign EXP_name to it
    dat$exposure <- EXP_name
  }
  if(!"id.exposure" %in% names(dat)) {
    # If "id.exposure" does not exist, create it and assign EXP_pheno to it
    dat$id.exposure <- EXP_pheno
  }
  return(dat)
}