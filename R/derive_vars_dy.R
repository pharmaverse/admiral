#' Derive Analysis Relative Day
#'
#' Adds the analysis relative day (`--DY`) to the dataset, i.e., study
#' day of analysis date.
#'
#' @param dataset Input dataset
#'
#'   The columns specified by the `reference_date` and the `source_vars` parameter are
#'   expected.
#'
#' @param reference_date The start date column, e.g., date of first treatment
#'
#'   A date or date-time object column is expected.
#'
#' @param source_vars A list of datetime or date variables created using `vars()` from
#'   which dates are to be extracted. This can either be a list of date(time) variables
#'   or named `--DY` variables and corresponding --DT(M) variables e.g.
#'    vars(TRTSDTM, ASTDTM, AENDT) or vars(TRTSDT, ASTDTM, AENDT,DEATHDY=DTHDT)
#'
#' @author Teckla Akinyi
#'
#' @details The study day is derived as number of days from the reference date
#'   to the end date. If it is nonnegative, one is added. I.e., the study day of the
#'   reference date is 1. The input ---DT(M) is converted to ---DY
#'
#' @return The input dataset with `--DY` corresponding to the `--DTM` or `--DT`
#'     source variable(s) added
#'
#' @keywords derivation ADaM timing
#'
#' @export
#'
#' @examples
#' library(lubridate)
#' library(dplyr)
#'
#' datain <- tibble::tribble(
#'  ~STUDYID, ~USUBJID, ~TRTSDTM, ~ASTDTM, ~AENDT,
#'  "TEST01", "PAT01", "2014-01-17T23:59:59", "2014-01-18T13:09:O9", "2014-01-20"
#' ) %>%
#'  mutate(TRTSDTM = lubridate::as_datetime(TRTSDTM),
#'         ASTDTM = lubridate::as_datetime(ASTDTM),
#'         AENDT = lubridate::ymd(AENDT))
#'
#' derive_vars_dy(datain, reference_date = TRTSDTM, source_vars = vars(TRTSDTM, ASTDTM, AENDT))
#'
#'  datain <- tibble::tribble(
#'  ~STUDYID, ~USUBJID, ~TRTSDT, ~ASTDTM, ~AENDT, ~DTHDT,
#'  "TEST01", "PAT01", "2014-01-17", "2014-01-18T13:09:O9", "2014-01-20", "2014-02-01"
#'  ) %>%
#'  mutate(TRTSDT = lubridate::ymd(TRTSDT),
#'         ASTDTM = lubridate::as_datetime(ASTDTM),
#'         AENDT = lubridate::ymd(AENDT),
#'         DTHDT = lubridate::ymd(DTHDT))
#'
#'  derive_vars_dy(datain,
#'                 reference_date = TRTSDT,
#'                 source_vars = vars(TRTSDT, ASTDTM, AENDT, DEATHDY = DTHDT))
derive_vars_dy <- function(dataset,
                           reference_date,
                           source_vars) {
  #assertions
  reference_date <- assert_symbol(enquo(reference_date))
  assert_vars(source_vars)
  assert_data_frame(dataset, required_vars = quo_c(source_vars, reference_date))

  #Warn if `--DY` variables already exist
  n_vars <- length(source_vars)
  source_names <- names(source_vars)
  dy_vars <- if_else(source_names == "",
                    stringr::str_replace_all(vars2chr(source_vars), "(DT|DTM)$", "DY"),
                    source_names)
  warn_if_vars_exist(dataset, dy_vars)

  if (n_vars > 1L) {
    dataset %>%
      mutate_at(.vars = source_vars,
                .funs = list(temp = ~ compute_duration(start_date = !!reference_date, end_date = .))
                ) %>%
      rename_at(vars(ends_with("temp")),
                .funs = ~ stringr::str_replace_all(., c(DTM_temp = "DY",
                                                        DT_temp = "DY",
                                                        DY_temp = "DY")))
  } else {
    dataset <-   dataset %>%
      mutate(!!sym(dy_vars) := compute_duration(!!reference_date, !!!source_vars)
             )
  }

}
