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
  deprecate_stop("0.7.0",
                 "derive_var_atirel()",
                 details = "Please use combination of `mutate()` and `case_when()` instead."
  )

}
