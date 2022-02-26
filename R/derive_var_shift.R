#' Derive Shift
#'
#' Derives a character shift variable containing concatenated shift in
#' values based on user-defined pairing, e.g., shift from baseline to
#' analysis value, shift from baseline grade to analysis grade, ...
#'
#' @param dataset Input dataset
#'
#'   The columns specified by `from_var`, and the `to_var`parameters are expected.
#'
#' @param new_var Name of the character shift variable to create
#'
#' @details `new_var` is derived by concatenating the values of `from_var` to values of `to_var`
#' (e.g. "NORMAL to HIGH"). When `from_var` or `to_var` has missing value, the
#' missing value is replaced by `NULL` (e.g. "NORMAL to NULL"), if both are missing `new_var` is
#' set to missing.
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
#' "P02",    "ALB",     37,    "Y",    "NORMAL","NORMAL",
#' "P02",    "ALB",     49,    "",     "NORMAL","HIGH",
#' "P02",    "SODIUM",  147,   "Y",    "HIGH",  "HIGH"
#' )
#'
#' derive_var_shift(
#'    data,
#'    new_var = SHIFT1,
#'    from_var = BNRIND,
#'    to_var = ANRIND
#' )
#'
derive_var_shift <- function(dataset,
                              new_var,
                              from_var,
                              to_var) {

  new_var <- assert_symbol(enquo(new_var))
  from_var <- assert_symbol(enquo(from_var))
  to_var <- assert_symbol(enquo(to_var))
  assert_data_frame(dataset, required_vars = vars(!!from_var, !!to_var))

  # If from_var or to_var has missing value then set to "NULL".

  #(Note: add codes to consider when from/to_var are numeric).
  dataset <- dataset %>%
    mutate(
      temp_from_var = case_when(
        is.na(!!from_var) & !is.na(!!to_var) ~ "NULL",
        TRUE ~ !!from_var
      ),
      temp_to_var = case_when(
        is.na(!!to_var) & !is.na(!!from_var) ~ "NULL",
        TRUE ~ !!to_var
      )
    ) %>%
    mutate(
      !!new_var := case_when (
        is.na(temp_from_var) & is.na(temp_to_var) ~ NA_character_,
        TRUE ~ paste(temp_from_var, temp_to_var, sep = " to ")
        )
    ) %>%
    select(-starts_with("temp_"))

  dataset
}
