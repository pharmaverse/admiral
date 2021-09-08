#' Derive date variables from datetime variables
#'
#' This function creates a date as output from a datetime variable
#'
#' @param dataset Input dataset
#'
#' @param source_vars A list of datetime variables from which dates are to be extracted
#'
#' @author Teckla Akinyi
#'
#' @return A data frame containing all observations and variables of the input
#'   dataset and adds the corresponding date variable of all specified datetime variables from
#'   the source_vars option with the correct name.
#'
#' @keywords ADaM Timing Date
#'
#' @export
#'
#' @examples
#' library(dplyr, warn.conflicts = FALSE)
#' library(lubridate, warn.conflicts = FALSE)
#'
#' adcm <- tibble::tribble(
#'   ~USUBJID, ~TRTSDTM,              ~ASTDTM,               ~AENDTM,               ~ASTTMF,
#'   "PAT01",  "2012-02-25 23:00:00", "2012-02-28 19:00:00", "2012-02-25 23:00:00", NA,
#'   "PAT01",  NA,                    "2012-02-28 19:00:00", NA,                    NA,
#'   "PAT01",  "2017-02-25 23:00:00", "2013-02-25 19:00:00", "2014-02-25 19:00:00", NA,
#'   "PAT01",  "2017-02-25 16:00:00", "2017-02-25 14:00:00", "2017-03-25 23:00:00", "M",
#'   "PAT01",  "2017-02-25 16:00:00", "2017-02-25 14:00:00", "2017-04-29 14:00:00", NA
#' ) %>%
#'   mutate(
#'     TRTSDTM = as_datetime(TRTSDTM),
#'     ASTDTM = as_datetime(ASTDTM),
#'     AENDTM = as_datetime(AENDTM)
#'   )
#'
#' derive_vars_dtm_to_dt(adcm, vars(TRTSDTM, ASTDTM, AENDTM))
derive_vars_dtm_to_dt <- function(dataset, source_vars) {
  assert_vars(source_vars)
  assert_data_frame(dataset, required_vars = source_vars)

  dataset %>%
    mutate_at(source_vars, .funs = list(new = lubridate::date)) %>%
    rename_at(vars(ends_with("new")), .funs = ~str_remove(., "M_new"))
}
