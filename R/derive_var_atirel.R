#' Derive Time Relative to Reference
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
#' @keywords ADaM Relationship Var ATIREL
#'
#' @export
#'
#' @examples
#' library(dplyr, warn.conflicts = FALSE)
#' adcm <- tibble::tribble(
#'   ~STUDYID, ~USUBJID, ~TRTSDTM, ~ASTDTM, ~AENDTM, ~ASTTMF,
#'   "TEST01", "PAT01", "2012-02-25 23:00:00", "2012-02-28 19:00:00", "2012-02-25 23:00:00", "",
#'   "TEST01", "PAT01", "", "2012-02-28 19:00:00", "", "",
#'   "TEST01", "PAT01", "2017-02-25 23:00:00", "2013-02-25 19:00:00", "2014-02-25 19:00:00", "",
#'   "TEST01", "PAT01", "2017-02-25 16:00:00", "2017-02-25 14:00:00", "2017-03-25 23:00:00", "m",
#'   "TEST01", "PAT01", "2017-02-25 16:00:00", "2017-02-25 14:00:00", "2017-04-29 14:00:00", ""
#' ) %>% dplyr::mutate(
#'   TRTSDTM = lubridate::as_datetime(TRTSDTM),
#'   ASTDTM = lubridate::as_datetime(ASTDTM),
#'   AENDTM = lubridate::as_datetime(AENDTM)
#' )
#'
#' derive_var_atirel(
#'   dataset = adcm,
#'   flag_var = ASTTMF,
#'   new_var = ATIREL
#' )
#'
derive_var_atirel <- function(dataset,
                              flag_var,
                              new_var) {
  # checks
  flag_var <- assert_symbol(enquo(flag_var))
  assert_data_frame(dataset,
    required_vars = vars(STUDYID, USUBJID, TRTSDTM, ASTDTM, AENDTM, !!flag_var)
  )
  new_var <- assert_symbol(enquo(new_var))
  warn_if_vars_exist(dataset, quo_text(new_var))

  #logic to create ATIREL
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
