#' Derive date variables from datetime variables
#'
#' This function creates a date as output from a datetime variable
#'
#'
#' @param `dataset` Input dataset
#'
#'
#' @param `source_vars` A list of datetime variables from which dates are to be extracted
#'
#'
#' @details
#'
#'
#' @author Teckla Akinyi
#'
#' @return A data frame containing all observations and variables of the input
#'   dataset and adds the corresponding date variable of all specified datetime variables from
#'   the source_vars option with the correct name.
#'
#'
#' @keywords ADaM Timing Date
#'
#' @export
#'
#' @seealso
#'
#' @examples
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
#'derive_vars_dtm_to_dt(dataset=adcm,
#'                  source_vars=vars(TRTSDTM,ASTDTM,AENDTM))

derive_vars_dtm_to_dt <- function(dataset,
                              source_vars) {
  assert_data_frame(dataset, source_vars)

  dataset %>%
    dplyr::mutate_at(.vars = source_vars,
                    .funs = list(new = ~ lubridate::date(.))) %>%
    dplyr::rename_at(.vars = vars(ends_with("new")),
                     .funs = ~ stringr::str_remove(., "M_new"))
}
