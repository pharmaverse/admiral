#' Compute Factor for Value Imputations When Character Value Contains < or >
#'
#' Function to compute factor for value imputation when character
#' value contains < or >.  The factor is calculated using the number of decimals.
#' If there are no decimals, the factor is 1, otherwise the factor = 1/10^decimal
#' place.  For example, the factor for 100 = 1, the factor for 5.4 = 1/10^1,
#' the factor for 5.44 = 1/10^2.  This results in no additional false precision
#' added to the value.  This is an intermediate function.
#'
#' Derive an imputed value
#'
#' @param character_value_decimal Character value to determine decimal precision
#'
#' @return Decimal precision value to add or subtract
#'
#' @export
#'
#' @details This is a vector-oriented helper and is not usually called directly on a data
#' frame with `%>%`.
#'
#' @keywords com_bds_findings
#' @family com_bds_findings
#'
#' @examples
#' compute_qual_imputation_dec("<40.1")
#' compute_qual_imputation_dec(c("0.35", "1"))
compute_qual_imputation_dec <- function(character_value_decimal) {
  x <- str_trim(character_value_decimal)

  decimal <- case_when(
    str_detect(x, "\\.") ~ 1 / (10^(str_length(x) - str_locate(x, "\\.")[, 1])),
    TRUE ~ 1
  )

  unname(decimal)
}

#' Impute Numeric Values from Character Results with Qualifiers
#'
#' Converts character results to numeric values and optionally shifts results
#' with strict inequality qualifiers.
#'
#' @param character_value Character vector of numeric results, optionally prefixed
#' with a qualifier, for example, `c("<40", ">=40.1", "3")`.
#'
#' @param imputation_type Method used to derive numeric results:
#'
#' - `"strip"`: Remove qualifiers and use `offset` as the additive offset for
#'   strict inequalities.
#' - `"shift"`: Remove qualifiers and use the decimal precision of each original
#'   character result as the additive offset. See Details.
#'
#' @permitted `"strip"`, `"shift"`
#'
#' @param factor `r lifecycle::badge("deprecated")` Deprecated alias for `offset`.
#' Use `offset` instead. Supplying
#' both arguments is an error. This argument retains its position for compatibility.
#'
#' @param offset Finite, nonnegative numeric scalar used when
#' `imputation_type = "strip"`. It is added for results with `>` and subtracted
#' for results with `<`. It is ignored when `imputation_type = "shift"`.
#'
#' @details
#' Qualified results must match `^[<>]=?[0-9]+[.]?[0-9]*$`, for example,
#' `"<40"`, `">=40.1"`, or `">40."`. Unqualified decimal numbers, optionally
#' signed, are also accepted, including `"-3"`, `"+3"`, and `".5"`. Whitespace,
#' scientific notation, signed qualified values, and malformed qualifiers are
#' not accepted. Invalid strings and missing values return `NA_real_` without
#' a warning. The order and length of the input are preserved.
#'
#' For results with a strict `>` qualifier, the offset is added; for results with
#' a strict `<` qualifier, it is subtracted. Results with `>=` or `<=`, and results
#' without an inequality qualifier, are converted without a shift.
#'
#' With `imputation_type = "shift"`, the offset is \eqn{10^{-d}}, where \eqn{d} is
#' the number of decimal places in the original character result, including
#' trailing zeros. For integer results, the offset is 1. For example, `"<40"`
#' becomes 39, `">40.1"` becomes 40.2, and `">40.10"` becomes 40.11.
#'
#' Numeric methods `1` and `2` are deprecated aliases for `"strip"` and `"shift"`,
#' respectively. They continue to work with a deprecation message.
#'
#' The appropriate imputation rule depends on the analysis specification.
#' A shift can cross zero; results are not restricted to nonnegative values.
#'
#' @return A numeric vector of imputed results.
#'
#' @export
#'
#' @keywords com_bds_findings
#' @family com_bds_findings
#'
#' @examples
#' values <- c("<40", ">40.1", "<=40", ">=40.1", "3")
#' compute_qual_imputation(values)
#' compute_qual_imputation(values, offset = 0.5)
#' compute_qual_imputation(values, imputation_type = "shift")
#'
#' # Trailing zeros retain the reported precision
#' compute_qual_imputation(c(">40.1", ">40.10"), imputation_type = "shift")
#'
#' # Malformed qualifiers produce missing values without warnings
#' compute_qual_imputation(c("0>1", "-<0.9>", ">100.00", "<1"), "shift")
compute_qual_imputation <- function(character_value,
                                    imputation_type = "strip",
                                    factor = NULL,
                                    offset = 0) {
  assert_character_vector(character_value)

  if (is.numeric(imputation_type)) {
    if (length(imputation_type) != 1L || is.na(imputation_type) ||
      !imputation_type %in% c(1, 2)) {
      cli::cli_abort(
        "Numeric {.arg imputation_type} must be 1 or 2.",
        class = "compute_qual_imputation_method"
      )
    }
    deprecate_inform(
      "1.6.0", "compute_qual_imputation(imputation_type = 'must be a string')",
      details = paste(
        "Use `imputation_type = \"strip\"` instead of 1",
        "and `imputation_type = \"shift\"` instead of 2."
      )
    )
    imputation_type <- c("strip", "shift")[[imputation_type]]
  }
  assert_character_scalar(imputation_type, values = c("strip", "shift"))

  if (!missing(factor)) {
    if (!missing(offset)) {
      cli::cli_abort(
        "Supply only {.arg offset}, not both {.arg factor} and {.arg offset}.",
        class = "compute_qual_imputation_arguments"
      )
    }
    assert_numeric_vector(factor, length = 1)
    deprecate_inform(
      "1.6.0", "compute_qual_imputation(factor = )", "compute_qual_imputation(offset = )"
    )
    offset <- factor
  }
  assert_numeric_vector(offset, length = 1)
  if (!is.finite(offset) || offset < 0) {
    cli::cli_abort(
      "Argument {.arg offset} must be finite and nonnegative.",
      class = "compute_qual_imputation_offset"
    )
  }

  valid <- str_detect(
    character_value,
    "^(?:[<>]=?[0-9]+\\.?[0-9]*|[+-]?(?:[0-9]+\\.?[0-9]*|\\.[0-9]+))\\z"
  )
  clean_value <- if_else(!is.na(valid) & valid, character_value, NA_character_)
  numeric_value <- as.numeric(str_remove(clean_value, "^[<>]=?"))

  if (imputation_type == "shift") {
    offset <- compute_qual_imputation_dec(clean_value)
  }
  replace_when(
    numeric_value,
    str_detect(clean_value, "^>") & !str_detect(clean_value, "^>=") ~ numeric_value + offset,
    str_detect(clean_value, "^<") & !str_detect(clean_value, "^<=") ~ numeric_value - offset
  )
}
