#' Derive Time Variables from Datetime Variables
#'
#' This function creates time variable(s) as output from datetime variable(s)
#'
#' @param dataset Input dataset
#'
#' @param source_vars A list of datetime variables created using `vars()` from
#'   which time is to be extracted
#'
#' @author Teckla Akinyi
#'
#' @details
#' The names of the newly added variables are automatically set by replacing the
#' `--DTM` suffix of the `source_vars` with `--TM`. The `--TM` variables are created
#' using the {hms} package.
#'
#' @return
#' A data frame containing the input dataset with the corresponding time
#' (`--TM`) variable(s) of all datetime variables (`--DTM`) specified in
#' `source_vars` with the correct name.
#'
#' @family der_date_time
#'
#' @keywords der_gen der_date_time
#'
#' @export
#'
#' @examples
#' library(tibble)
#' library(dplyr, warn.conflicts = FALSE)
#' library(lubridate)
#'
#' adcm <- tribble(
#'   ~USUBJID, ~TRTSDTM, ~ASTDTM, ~AENDTM,
#'   "PAT01", "2012-02-25 23:41:10", "2012-02-28 19:03:00", "2013-02-25 23:32:16",
#'   "PAT01", "", "2012-02-28 19:00:00", "",
#'   "PAT01", "2017-02-25 23:00:02", "2013-02-25 19:00:15", "2014-02-25 19:00:56",
#'   "PAT01", "2017-02-25 16:00:00", "2017-02-25 14:25:00", "2017-03-25 23:00:00",
#'   "PAT01", "2017-02-25 16:05:17", "2017-02-25 14:20:00", "2018-04-29 14:06:45",
#' ) %>%
#'   mutate(
#'     TRTSDTM = as_datetime(TRTSDTM),
#'     ASTDTM = as_datetime(ASTDTM),
#'     AENDTM = as_datetime(AENDTM)
#'   )
#'
#' adcm %>%
#'   derive_vars_dtm_to_tm(vars(TRTSDTM)) %>%
#'   select(USUBJID, starts_with("TRT"), everything())
#'
#' adcm %>%
#'   derive_vars_dtm_to_tm(vars(TRTSDTM, ASTDTM, AENDTM)) %>%
#'   select(USUBJID, starts_with("TRT"), starts_with("AS"), starts_with("AE"))
derive_vars_dtm_to_tm <- function(dataset, source_vars) {
  assert_vars(source_vars)
  assert_data_frame(dataset, required_vars = source_vars)

  # Warn if `--TM` variables already exist
  dtm_vars <- quo_c(source_vars)
  dtm_vars2 <- vars2chr(dtm_vars)
  n_vars <- length(dtm_vars)

  for (i in n_vars) {
    tm_vars <- str_replace(dtm_vars2[[i]], "DTM", "TM")
    warn_if_vars_exist(dataset, tm_vars)
  }

  if (n_vars > 1L) {
    dataset %>%
      mutate_at(source_vars, .funs = list(new = as_hms)) %>%
      rename_at(vars(ends_with("new")), .funs = ~ str_replace(., "DTM_new", "TM"))
  } else {
    dataset %>%
      mutate(!!sym(tm_vars) := as_hms(!!sym(dtm_vars2)))
  }
}
