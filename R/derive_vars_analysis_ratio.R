#' Ratio Derivations
#'
#' Derives a character shift variable containing concatenated shift in
#' values based on user-defined pairing, e.g., shift from baseline to
#' analysis value, shift from baseline grade to analysis grade, ...
#'
#' @param dataset Input dataset
#'
#' @param new_var Name of the character shift variable to create.
#'
#' @param from_var Variable containing value to shift from.
#'
#' @param to_var Variable containing value to shift to.
#'
#' @param na_val Character string to replace missing values in `from_var` or `to_var`.
#'
#'  Default: "NULL"
#'
#' @param sep_val Character string to concatenate values of `from_var` and `to_var`.
#'
#'  Default: " to "
#'
#' @param filter_post_baseline Condition to identify post-baseline records. If not
#' specified, `new_var` is populated for all records.
#'
#'  Default: `NULL`
#'
#' @details `new_var` is derived by concatenating the values of `from_var` to values of `to_var`
#' (e.g. "NORMAL to HIGH"). When `from_var` or `to_var` has missing value, the
#' missing value is replaced by `na_val` (e.g. "NORMAL to NULL"). If `filter_post_baseline` is
#' specified, `new_var` is populated only for post-baseline records. If `filter_post_baseline`
#' is `NULL` (default), `new_var` is populated for all records.
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
derive_vars_analysis_ratio <- function(dataset) {
  AVAL <- assert_symbol(enquo(AVAL))
  BASE <- assert_symbol(enquo(BASE))
  ANRLO <- assert_symbol(enquo(ANRLO))
  ANRHI <- assert_symbol(enquo(ANRHI))

  assert_data_frame(dataset, required_vars = quo_c(AVAL, BASE, ANRLO, ANRHI))

  dataset %>%
    mutate(
      R2BASE = if_else(!is.na(AVAL) & !is.na(BASE) & BASE != 0, AVAL / BASE, NA_real_),
      R2ANRLO = if_else(!is.na(AVAL) & !is.na(ANRLO) & ANRLO != 0, AVAL / ANRLO, NA_real_),
      R2ANRHI = if_else(!is.na(AVAL) & !is.na(ANRHI) & ANRHI != 0, AVAL / ANRHI, NA_real_)
    )
}

