#' Ratio Derivations
#'
#' Derives a three ratio variables for adlb dataset.
#'
#' @param dataset Input dataset
#'
#' @param soure_vars Columns `AVAL`, `BASE`, `ANRLO` and `ANRHI` are expected.
#'
#' @param ratio_vars Users can select from three values: "all", "low" and "high".
#' Selecting all will append `R2BASE`, `R2ANRLO`, `R2ANRHI` to the dataset.
#' Selecting "low" will append only `R2BASE`, `R2ANRLO` and "high" only
#' `R2BASE`, `R2ANRHI`.
#'
#'  Default: "all"
#'
#' @details
#'
#'
#' @author Ben Straub
#'
#' @return The input dataset with the ratio variables added
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
#'   derive_var_analysis_ratio()
#'
derive_var_analysis_ratio <- function(dataset,
                                      numer_var,
                                      denom_var,
                                      ratio_var
                                      ) {
  numer_var <- assert_symbol(enquo(numer_var))
  denom_var <- assert_symbol(enquo(denom_var))
  ratio_var <- assert_symbol(enquo(ratio_var))

  assert_data_frame(dataset, required_vars = quo_c(numer_var, denom_var))

  dataset <- dataset %>%
    mutate(
      !!ratio_var := if_else(!is.na(!!numer_var) &
                            !is.na(!!denom_var) &
                            !!denom_var != 0,
                            !!numer_var / !!denom_var,
                            NA_real_),
      )

  dataset
}

