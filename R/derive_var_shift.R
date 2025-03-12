#' Derive Shift
#'
#' Derives a character shift variable containing concatenated shift in
#' values based on user-defined pairing, e.g., shift from baseline to
#' analysis value, shift from baseline grade to analysis grade, ...
#'
#' @param dataset
#' `r roxygen_param_dataset(expected_vars = c("from_var", "to_var"))`
#'
#' @param new_var Name of the character shift variable to create.
#'
#' @param from_var Variable containing value to shift from.
#'
#' @param to_var Variable containing value to shift to.
#'
#' @param missing_value Character string to replace missing values in `from_var` or `to_var`.
#'
#'
#' @param sep_val Character string to concatenate values of `from_var` and `to_var`.
#'
#'
#' @details `new_var` is derived by concatenating the values of `from_var` to values of `to_var`
#' (e.g. "NORMAL to HIGH"). When `from_var` or `to_var` has missing value, the
#' missing value is replaced by `missing_value` (e.g. "NORMAL to NULL").
#'
#'
#' @return The input dataset with the character shift variable added
#'
#' @family der_bds_findings
#' @keywords der_bds_findings
#'
#' @export
#'
#' @examples
#' library(tibble)
#'
#' data <- tribble(
#'   ~USUBJID, ~PARAMCD, ~AVAL, ~ABLFL, ~BNRIND,  ~ANRIND,
#'   "P01",    "ALB",       33, "Y",    "LOW",    "LOW",
#'   "P01",    "ALB",       38, NA,     "LOW",    "NORMAL",
#'   "P01",    "ALB",       NA, NA,     "LOW",    NA,
#'   "P02",    "ALB",       37, "Y",    "NORMAL", "NORMAL",
#'   "P02",    "ALB",       49, NA,     "NORMAL", "HIGH",
#'   "P02",    "SODIUM",   147, "Y",    "HIGH",   "HIGH"
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
#'   restrict_derivation(
#'     derivation = derive_var_shift,
#'     args = params(
#'       new_var = SHIFT1,
#'       from_var = BNRIND,
#'       to_var = ANRIND
#'     ),
#'     filter = is.na(ABLFL)
#'   )
derive_var_shift <- function(dataset,
                             new_var,
                             from_var,
                             to_var,
                             missing_value = "NULL",
                             sep_val = " to ") {
  new_var <- assert_symbol(enexpr(new_var))
  from_var <- assert_symbol(enexpr(from_var))
  to_var <- assert_symbol(enexpr(to_var))
  missing_value <- assert_character_scalar(missing_value)
  sep_val <- assert_character_scalar(sep_val)
  assert_data_frame(dataset, required_vars = exprs(!!from_var, !!to_var))

  # Derive shift variable. If from_var or to_var has missing value then set to missing_value.
  dataset %>%
    mutate(
      temp_from_var = if_else(is.na(!!from_var), !!missing_value, as.character(!!from_var)),
      temp_to_var = if_else(is.na(!!to_var), !!missing_value, as.character(!!to_var))
    ) %>%
    mutate(
      !!new_var := paste(temp_from_var, temp_to_var, sep = !!sep_val)
    ) %>%
    select(-starts_with("temp_"))
}
