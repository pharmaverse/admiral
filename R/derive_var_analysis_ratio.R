#' Derive Ratio Variable
#'
#' Derives a ratio variable for a BDS dataset, e.g Ratio to Baseline `"AVAL / BASE"`,
#'  Ratio to Analysis Range y Lower Limit `"AVAL / AyLO"`,  or Ratio to Analysis
#'  Range y Upper Limit `"AVAL / AyHI"`.
#'
#' @param dataset Input dataset
#'
#' @param numer_var Variable containing values to be used in the numerator of
#'  the ratio calculation.
#'
#' @param denom_var Variable containing values to be used in the denominator of
#'  the ratio calculation.
#'
#' @param ratio_var_suffix User supplied suffix to be concatenated with `"R2"`
#'
#' @param override A user wishing to remove the default prefix and supply
#'  their own completely new variable can set this to TRUE.
#'
#'  Default is `"FALSE"`.
#'
#' @param new_var Name of the new ratio variable being created.
#'
#' @details Reference CDISC ADaM Implementation Guide Version 1.1
#' Section 3.3.4 Analysis Parameter Variables for BDS Datasets
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
#'  derive_var_analysis_ratio(numer_var = AVAL, denom_var = BASE, ratio_var_suffix = "BASE") %>%
#'  derive_var_analysis_ratio(numer_var = AVAL, denom_var = ANRLO, ratio_var_suffix = "ANR1LO") %>%
#'  derive_var_analysis_ratio(numer_var = AVAL, denom_var = ANRHI, ratio_var_suffix = "ANR1HI")
#'
derive_var_analysis_ratio <- function(dataset,
                                      numer_var,
                                      denom_var,
                                      #ratio_var_suffix,
                                      #override = FALSE,
                                      new_var = NULL) {
  numer_var <- assert_symbol(enquo(numer_var))
  denom_var <- assert_symbol(enquo(denom_var))

  assert_data_frame(dataset, required_vars = quo_c(numer_var, denom_var))

  if (is.null(new_var)) {

    std_var_list <- list(BASE = "R2BASE", ANRLO = "R2ANRLO", ANRHI = "R2ANRHI")
    std_var <-unlist(std_var_list[denom_var])

      dataset <- dataset %>%
        mutate(
          !!sym(std_var) := if_else(!is.na(!!numer_var) &
            !is.na(!!denom_var) &
            !!denom_var != 0,
          !!numer_var / !!denom_var,
          NA_real_
          )
        )
  } else {
    dataset <- dataset %>%
      mutate(
        !!sym(new_var) := if_else(!is.na(!!numer_var) &
          !is.na(!!denom_var) &
          !!denom_var != 0,
        !!numer_var / !!denom_var,
        NA_real_
        ),
      )
  }

  dataset
}
