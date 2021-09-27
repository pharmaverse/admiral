#' Derive Time Variables From Datetime Variables
#'
#' This function creates a time as output from a datetime variable
#'
#' @param dataset Input dataset
#'
#' @param source_vars A list of datetime variables from which time are to be extracted
#'
#' @author Teckla Akinyi
#'
#' @return A data frame containing all observations and variables of the input
#'   dataset and adds the corresponding time variable of all specified datetime variables from
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
#'  ~STUDYID, ~USUBJID, ~TRTSDTM,               ~ASTDTM,               ~AENDTM,
#'  "TEST01", "PAT01",  "2012-02-25 23:41:10", "2012-02-28 19:03:00", "2013-02-25 23:32:16",
#'  "TEST01", "PAT01",  "",                    "2012-02-28 19:00:00", "",
#'  "TEST01", "PAT01",  "2017-02-25 23:00:02", "2013-02-25 19:00:15", "2014-02-25 19:00:56",
#'  "TEST01", "PAT01",  "2017-02-25 16:00:00", "2017-02-25 14:25:00", "2017-03-25 23:00:00",
#'  "TEST01", "PAT01",  "2017-02-25 16:05:17", "2017-02-25 14:20:00", "2018-04-29 14:06:45",
#'  ) %>% mutate(
#'  TRTSDTM = as_datetime(TRTSDTM),
#'  ASTDTM = as_datetime(ASTDTM),
#'  AENDTM = as_datetime(AENDTM)
#'  )
#'
#' derive_vars_dtm_to_tm(adcm,
#'                       vars(TRTSDTM, ASTDTM, AENDTM))

derive_vars_dtm_to_tm <- function(dataset,
                                  source_vars) {

    assert_vars(source_vars)
    assert_data_frame(dataset, required_vars = source_vars)

    dataset %>%
      mutate_at(source_vars, .funs = list(new = hms::as_hms)) %>%
      rename_at(vars(ends_with("new")), .funs = ~str_replace(., "DTM_new", "TM"))

}
