#' Ratio Derivations
#'
#' Derives a three ratio variables for adlb dataset.
#'
#' @param dataset Input dataset
#'
#' @param soure_vars
#'
#' @param ratio_vars
#'
#'  Default: "all"
#'
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
#'   derive_vars_analysis_ratio(ratio_vars = "high")
#'
derive_vars_analysis_ratio <- function(dataset,
                                       soure_vars = vars(AVAL, BASE, ANRLO, ANRHI),
                                       ratio_vars = "all") {
  soure_vars <- assert_vars(soure_vars)
  ratio_vars <- assert_character_scalar(ratio_vars)

  assert_data_frame(dataset, required_vars = quo_c(soure_vars))

  dataset <- dataset %>%
    mutate(
      R2BASE = if_else(!is.na(AVAL) & !is.na(BASE) & BASE != 0, AVAL / BASE, NA_real_),
      R2ANRLO = if_else(!is.na(AVAL) & !is.na(ANRLO) & ANRLO != 0, AVAL / ANRLO, NA_real_),
      R2ANRHI = if_else(!is.na(AVAL) & !is.na(ANRHI) & ANRHI != 0, AVAL / ANRHI, NA_real_)
    )

  # Remove High/Low Variables based on user-input
  switch(ratio_vars,
    all = {
      dataset <- dataset %>%
        select(everything(), R2BASE, R2ANRLO, R2ANRHI)
    },

    low = {
      dataset <- dataset %>%
        select(everything(), R2BASE, R2ANRLO, -R2ANRHI)
    },

    high = {
      dataset <- dataset %>%
        select(everything(), R2BASE, -R2ANRLO, R2ANRHI)
    }
  )

  dataset
}

