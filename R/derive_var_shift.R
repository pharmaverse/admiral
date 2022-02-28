#' Derive Shift
#'
#' Derives a character shift variable containing concatenated shift in
#' values based on user-defined pairing, e.g., shift from baseline to
#' analysis value, shift from baseline grade to analysis grade, ...
#'
#' @param dataset Input dataset
#'
#'   The columns specified by `from_var` and the `to_var`parameters are expected.
#'
#' @param new_var Name of the character shift variable to create.
#'
#' @param from_var Variable containing value to shift from.
#'
#' @param to_var Variable containing value to shift to.
#'
#' @param post_baseline_condition Condition to identify post-baseline records. If not
#' specified, `new_var` is populated for all records.
#'
#'  Default: `NULL`
#'
#' @details `new_var` is derived by concatenating the values of `from_var` to values of `to_var`
#' (e.g. "NORMAL to HIGH"). When `from_var` or `to_var` has missing value, the
#' missing value is replaced by `NULL` (e.g. "NORMAL to NULL"). If `post_baseline_condition` is
#' specified, `new_var` is populated only for post-baseline records. If `post_baseline_condition`
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
#' data <- tibble::tribble(
#' ~USUBJID, ~PARAMCD, ~AVAL, ~ABLFL, ~BNRIND, ~ANRIND,
#' "P01",    "ALB",     33,    "Y",    "LOW",   "LOW",
#' "P01",    "ALB",     38,    "",     "LOW",   "NORMAL",
#' "P01",    "ALB",     NA,    "",     "LOW",   "",
#' "P02",    "ALB",     37,    "Y",    "NORMAL","NORMAL",
#' "P02",    "ALB",     49,    "",     "NORMAL","HIGH",
#' "P02",    "SODIUM",  147,   "Y",    "HIGH",  "HIGH"
#' )
#'
#' data %>%
#' convert_blanks_to_na %>%
#' derive_var_shift(
#'    new_var = SHIFT1,
#'    from_var = BNRIND,
#'    to_var = ANRIND
#' )
#'
#' # or only populate post-baseline records
#' data %>%
#' convert_blanks_to_na %>%
#' derive_var_shift(
#'     new_var = SHIFT1,
#'     from_var = BNRIND,
#'     to_var = ANRIND,
#'     post_baseline_condition = (ABLFL != "Y")
#' )
#'
derive_var_shift <- function(dataset,
                             new_var,
                             from_var,
                             to_var,
                             post_baseline_condition = NULL) {

  new_var <- assert_symbol(enquo(new_var))
  from_var <- assert_symbol(enquo(from_var))
  to_var <- assert_symbol(enquo(to_var))
  assert_data_frame(dataset, required_vars = vars(!!from_var, !!to_var))
  post_baseline_condition <- assert_filter_cond(enquo(post_baseline_condition), optional = TRUE)

  # Derive shift variable. If from_var or to_var has missing value then set to "NULL".
  dataset <- dataset %>%
    mutate(
      temp_from_var = case_when(
        is.na(!!from_var) ~ "NULL",
        TRUE ~ as.character(!!from_var)
      ),
      temp_to_var = case_when(
        is.na(!!to_var) ~ "NULL",
        TRUE ~ as.character(!!to_var)
      )
    ) %>%
    mutate(
      !!new_var := paste(temp_from_var, temp_to_var, sep = " to ")
    ) %>%
    select(-starts_with("temp_"))

  # If post-baseline condition specified, then only populate shift for post-baseline records.
  if (!is.null(quo_get_expr(post_baseline_condition))) {
    dataset <- dataset  %>%
      mutate(
        !!new_var := case_when(
          !eval(quo_get_expr(post_baseline_condition)) ~ NA_character_,
          TRUE ~ !!new_var
        )
      )
  }

  dataset
}
