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

#' Function to Impute Values When Qualifier Exists in Character Result
#'
#' Derive an imputed value
#'
#' @param character_value Character version of value to be imputed
#' @param imputation_type (default value=1)
#' Valid Values:
#' 1: Strip <, >, = and convert to numeric.
#' 2: imputation_type=1 and if the character value contains a < or >, the number of
#' of decimals associated with the character value is found and then a factor of
#' 1/10^(number of decimals + 1) will be added/subtracted from the numeric value.
#' If no decimals exists, a factor of 1/10 will be added/subtracted from the value.
#' @param factor Numeric value (default=0), when using `imputation_type` = 1, this
#' value can be added or subtracted when the qualifier is removed.
#'
#' @return The imputed value
#'
#' @importFrom dplyr case_when
#' @importFrom dplyr if_else
#' @importFrom stringr str_detect
#'
#' @export
#'
#'
#' @keywords com_bds_findings
#' @family com_bds_findings
#'
#' @examples
#' compute_qual_imputation("<40")
#' compute_qual_imputation(c("3", ">30.2"))
compute_qual_imputation <- function(character_value, imputation_type = 1, factor = 0) {
  clean_value <- character_value

  # Identify strings with letters and set to NA immediately to avoid warning with
  # as.numeric() later down
  has_letters <- str_detect(clean_value, "[A-Za-z]")
  clean_value[has_letters] <- NA_character_

  # Remove symbols from the remaining valid strings
  clean_value <- gsub("=|>|<", "", clean_value)

  numeric_value <- as.numeric(clean_value)

  if (imputation_type == 1) {
    numeric_value <-
      case_when(
        str_detect(character_value, ">") & !str_detect(character_value, "=") ~
          numeric_value + factor,
        str_detect(character_value, "<") & !str_detect(character_value, "=") ~
          numeric_value - factor,
        TRUE ~ numeric_value
      )
  }

  if (imputation_type == 2) {
    numeric_value <-
      case_when(
        str_detect(character_value, ">") & !str_detect(character_value, "=") ~
          numeric_value + compute_qual_imputation_dec(character_value),
        str_detect(character_value, "<") & !str_detect(character_value, "=") ~
          numeric_value - compute_qual_imputation_dec(character_value),
        TRUE ~ numeric_value
      )
  }

  numeric_value
}
