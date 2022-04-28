#' Derive Shift
#'
#' Derives a character shift variable containing concatenated shift in
#' values based on user-defined pairing, e.g., shift from baseline to
#' analysis value, shift from baseline grade to analysis grade, ...
#'
#' @param dataset Input dataset
#'
#'   The columns specified by `from_var` and the `to_var` parameters are expected.
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
#' @author Annie Yang
#'
#' @return The input dataset with the character shift variable added
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
derive_var_shift <- function(dataset,
                             new_var,
                             from_var,
                             to_var,
                             na_val = "NULL",
                             sep_val = " to ",
                             filter_post_baseline = NULL) {
  new_var <- assert_symbol(enquo(new_var))
  from_var <- assert_symbol(enquo(from_var))
  to_var <- assert_symbol(enquo(to_var))
  na_val <- assert_character_scalar(na_val)
  sep_val <- assert_character_scalar(sep_val)
  assert_data_frame(dataset, required_vars = vars(!!from_var, !!to_var))
  filter_post_baseline <- assert_filter_cond(enquo(filter_post_baseline), optional = TRUE)

  # Derive shift variable. If from_var or to_var has missing value then set to na_val.
  dataset <- dataset %>%
    mutate(
      temp_from_var = if_else(is.na(!!from_var), !!na_val, as.character(!!from_var)),
      temp_to_var = if_else(is.na(!!to_var), !!na_val, as.character(!!to_var))
    ) %>%
    mutate(
      !!new_var := paste(temp_from_var, temp_to_var, sep = !!sep_val)
    ) %>%
    select(-starts_with("temp_"))

  # If post-baseline condition specified, then only populate shift for post-baseline records.
  if (!quo_is_null(filter_post_baseline)) {
    dataset <- dataset  %>%
      mutate(
        !!new_var := if_else(!!filter_post_baseline, !!new_var, NA_character_,
                             missing = !!new_var)
      )
  }

  dataset
}
