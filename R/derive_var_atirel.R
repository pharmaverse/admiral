#' Derive Time Relative to Reference
#'
#' @description
#' `r lifecycle::badge("deprecated")`
#'
#' This function is *deprecated*, as it is deemed as too specific for admiral.
#' Derivations like this can be implemented calling `mutate()` and
#' `case_when()`.
#'
#' Derives the variable `ATIREL` to CONCOMITANT, PRIOR, PRIOR_CONCOMITANT or NULL
#' based on the relationship of cm Analysis start/end date/times to treatment start date/time
#'
#' @param dataset Input dataset
#'   The variables `TRTSDTM`, `ASTDTM`, `AENDTM` are expected
#' @param flag_var Name of the variable with Analysis Start Date Imputation Flag
#' @param new_var Name of variable to create
#'
#' @details `ATIREL` is set to:
#'    - null, if Datetime of First Exposure to Treatment is missing,
#'    - "CONCOMITANT", if the Analysis Start Date/Time is greater than or equal to Datetime of
#'       First Exposure to Treatment,
#'    - "PRIOR", if the Analysis End Date/Time is not missing and less than
#'       the Datetime of First Exposure to Treatment,
#'    - "CONCOMITANT" if the date part of Analysis Start Date/Time is equal to
#'       the date part of Datetime of First Exposure to Treatment and
#'       the Analysis Start Time Imputation Flag is 'H' or 'M',
#'    -  otherwise it is set to "PRIOR_CONCOMITANT".
#'
#' @author Teckla Akinyi
#'
#' @return A dataset containing all observations and variables of the input
#'   dataset and additionally the variable specified by the `new_var` parameter.
#'
#' @keywords deprecated
#'
#' @export
#'
derive_var_atirel <- function(dataset,
                              flag_var,
                              new_var) {
  # checks
  warn(paste(
    "`derive_var_atirel` is deprecated as of admiral 0.7.0.",
    "Please use `mutate()` and `case_when()` instead.",
    sep = "\n"
  ))

  flag_var <- assert_symbol(enquo(flag_var))
  assert_data_frame(dataset,
    required_vars = vars(STUDYID, USUBJID, TRTSDTM, ASTDTM, AENDTM, !!flag_var)
  )
  new_var <- assert_symbol(enquo(new_var))
  warn_if_vars_exist(dataset, quo_text(new_var))

  # logic to create ATIREL
  dataset %>%
    mutate(!!new_var :=
      case_when(
        is.na(TRTSDTM) ~ NA_character_,
        ASTDTM >= TRTSDTM ~ "CONCOMITANT",
        !is.na(AENDTM) & AENDTM < TRTSDTM ~ "PRIOR",
        date(ASTDTM) == date(TRTSDTM) & toupper(!!flag_var) %in% c("H", "M") ~ "CONCOMITANT",
        TRUE ~ "PRIOR_CONCOMITANT"
      ))
}
