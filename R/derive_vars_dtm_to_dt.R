#' Derive Date Variables from Datetime Variables
#'
#' This function creates date(s) as output from datetime variable(s)
#'
#' @param dataset Input dataset
#'
#' @param source_vars A list of datetime variables created using `vars()` from
#'   which dates are to be extracted
#'
#' @author Teckla Akinyi
#'
#' @return
#' A data frame containing the input dataset with the corresponding date (`--DT`)
#' variable(s) of all datetime variables (`--DTM`) specified in `source_vars.`
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
#'   ~USUBJID, ~TRTSDTM,              ~ASTDTM,               ~AENDTM,
#'   "PAT01",  "2012-02-25 23:00:00", "2012-02-28 19:00:00", "2012-02-25 23:00:00",
#'   "PAT01",  NA,                    "2012-02-28 19:00:00", NA,
#'   "PAT01",  "2017-02-25 23:00:00", "2013-02-25 19:00:00", "2014-02-25 19:00:00",
#'   "PAT01",  "2017-02-25 16:00:00", "2017-02-25 14:00:00", "2017-03-25 23:00:00",
#'   "PAT01",  "2017-02-25 16:00:00", "2017-02-25 14:00:00", "2017-04-29 14:00:00",
#' ) %>%
#'   mutate(
#'     TRTSDTM = as_datetime(TRTSDTM),
#'     ASTDTM = as_datetime(ASTDTM),
#'     AENDTM = as_datetime(AENDTM)
#'   )
#'
#' adcm %>%
#'   derive_vars_dtm_to_dt(vars(TRTSDTM, ASTDTM, AENDTM)) %>%
#'   select(USUBJID, starts_with("TRT"), starts_with("AST"), starts_with("AEN"))
derive_vars_dtm_to_dt <- function(dataset, source_vars) {
  assert_vars(source_vars)
  assert_data_frame(dataset, required_vars = source_vars)

  # Warn if `--TM` variables already exist
  dtm_vars <- quo_c(source_vars)
  dtm_vars2 <- vars2chr(dtm_vars)
  n_vars <- length(dtm_vars)

  for (i in n_vars) {
    dt_vars <- str_replace(dtm_vars2[[i]], "DTM", "DT")
    warn_if_vars_exist(dataset, dt_vars)
  }

  if (n_vars > 1L) {
    dataset %>%
      mutate_at(source_vars, .funs = list(new = lubridate::date)) %>%
      rename_at(vars(ends_with("new")), .funs = ~ str_replace(., "DTM_new", "DT"))
  } else {
    dataset %>%
      mutate(!!sym(dt_vars) := lubridate::date(!!sym(dtm_vars2)))
  }
}
