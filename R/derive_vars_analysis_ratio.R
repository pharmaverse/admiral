#' Ratio Derivations
#'
#' Derives a character shift variable containing concatenated shift in
#' values based on user-defined pairing, e.g., shift from baseline to
#' analysis value, shift from baseline grade to analysis grade, ...
#'
#' @param dataset Input dataset
#'
#' @param
#'
#' @param
#'
#' @param
#'
#' @param
#'
#'  Default: "NULL"
#'
#' @param
#'
#'  Default: " to "
#'
#' @param
#'
#'  Default: `NULL`
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
#' data <- tibble::tribble(
#'   ~USUBJID, ~PARAMCD, ~AVAL, ~ABLFL, ~BNRIND, ~ANRIND,
#'   "P01",    "ALB",     33,    "Y",    "LOW",   "LOW",
#'   "P01",    "ALB",     38,    NA,     "LOW",   "NORMAL",
#'   "P01",    "ALB",     NA,    NA,     "LOW",   NA,
#'   "P02",    "ALB",     37,    "Y",    "NORMAL","NORMAL",
#'   "P02",    "ALB",     49,    NA,     "NORMAL","HIGH",
#'   "P02",    "SODIUM",  147,   "Y",    "HIGH",  "HIGH"
#' )
#'
#' data %>%
#'   convert_blanks_to_na() %>%
#'   derive_var_shift(
#'     new_var = SHIFT1,
#'     from_var = BNRIND,
#'     to_var = ANRIND
#'   )
#'
#' # or only populate post-baseline records
#' data %>%
#'   convert_blanks_to_na() %>%
#'   derive_var_shift(
#'     new_var = SHIFT1,
#'     from_var = BNRIND,
#'     to_var = ANRIND,
#'     filter_post_baseline = ABLFL != "Y"
#'   )
#'
derive_vars_analysis_ratio <- function(dataset,
                                       soure_vars = vars(AVAL, BASE, ANRLO, ANRHI),
                                       ratio_vars = "all") {

  soure_vars <- assert_vars(soure_vars)
  ratio_vars <- assert_character_scalar(ratio_vars)

  assert_data_frame(dataset, required_vars = quo_c(soure_vars))

  dataset  <- dataset %>%
    mutate(
      R2BASE = if_else(!is.na(AVAL) & !is.na(BASE) & BASE != 0, AVAL / BASE, NA_real_),
      R2ANRLO = if_else(!is.na(AVAL) & !is.na(ANRLO) & ANRLO != 0, AVAL / ANRLO, NA_real_),
      R2ANRHI = if_else(!is.na(AVAL) & !is.na(ANRHI) & ANRHI != 0, AVAL / ANRHI, NA_real_)
    )

  # Remove High/Low Variables based on user-input
  switch(ratio_vars,
           all  = {dataset <- dataset %>%
             select(everything(), R2BASE, R2ANRLO, R2ANRHI)},

           low  = {dataset <- dataset %>%
             select(everything(), R2BASE, R2ANRLO, -R2ANRHI)},

           high = {dataset <- dataset %>%
             select(everything(), R2BASE, -R2ANRLO, R2ANRHI)}
  )

    dataset

}

