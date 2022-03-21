#' Ratio Derivations
#'
#' Derives a ratio variable for the adlb dataset, e.g Ratio to Baseline, Ratio
#' to Analysis Range y Lower Limit or Ratio to Analysis Range y Upper Limit.
#'
#' @param dataset Input dataset
#'
#' @param numer_var Variable containing values to be used in the numerator of
#'  the ratio calculation.
#'
#' @param denom_var Variable containing values to be used in the numerator of the of
#'  the ratio calculation.
#'
#' @param ratio_var_suffix User supplied suffix to be concatenated with `"R2"`
#'
#' @param override_prefix A user wishing to remove the default prefix and supply
#'  their own complete new variable can set this to TRUE.
#'
#'  Default is `"FALSE"`.
#'
#' @param new_var Name of the new ratio variable being created.
#'
#' @details
#'
#'
#' @author Ben Straub
#'
#' @return The input dataset with a ratio variable appended
#'
#' @keywords adam bds adlb derivation
#'
#' @export
#'
#'
#' @examples
#' library(dplyr, warn.conflicts = FALSE)
#'
#' data <-tibble::tribble(
#'  ~USUBJID, ~PARAMCD, ~SEQ,  ~AVAL, ~BASE, ~ANRLO, ~ANRHI,
#'  "P01",    "ALT",     1,     27,    27,    6,      34,
#'  "P01",    "ALT",     2,     41,    27,    6,      34,
#'  "P01",    "ALT",     3,     17,    27,    6,      34,
#'  "P02",    "ALB",     1,     38,    38,    33,     49,
#'  "P02",    "ALB",     2,     39,    38,    33,     49,
#'  "P02",    "ALB",     3,     37,    38,    33,     49
#')
#'
#' data %>%
#'  derive_var_analysis_ratio(numer_var = AVAL, denom_var = BASE, ratio_var = BASE) %>%
#'  derive_var_analysis_ratio(numer_var = AVAL, denom_var = ANRLO, ratio_var = ANR1LO) %>%
#'  derive_var_analysis_ratio(numer_var = AVAL, denom_var = ANRHI, ratio_var = ANR1HI)
#'
derive_var_analysis_ratio <- function(dataset,
                                      numer_var,
                                      denom_var,
                                      ratio_var_suffix,
                                      override_prefix = FALSE,
                                      new_var
                                      ) {
  numer_var <- assert_symbol(enquo(numer_var))
  denom_var <- assert_symbol(enquo(denom_var))
  #ratio_var_suffix <- assert_symbol(enquo(ratio_var_suffix))

  assert_data_frame(dataset, required_vars = quo_c(numer_var, denom_var))

  if(!override_prefix){
  r2 <- paste0("R2", ratio_var_suffix)
  dataset <- dataset %>%
    mutate(
      !!sym(r2) := if_else(!is.na(!!numer_var) &
                            !is.na(!!denom_var) &
                            !!denom_var != 0,
                            !!numer_var / !!denom_var,
                            NA_real_)
    )
  } else {

    dataset <- dataset %>%
      mutate(
        !!sym(new_var) := if_else(!is.na(!!numer_var) &
                               !is.na(!!denom_var) &
                               !!denom_var != 0,
                             !!numer_var / !!denom_var,
                             NA_real_),
      )
  }

  dataset
}

