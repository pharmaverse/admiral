#' Derive On-Treatment Flag Variable
#'
#' Derive on-treatment flag (`ONTRTFL`) in an ADaM BDS dataset with a single
#' assessment date
#'
#' @param dataset `data.frame`.
#'
#' @param start_date The assessment date (e.g. `ADT` for the date of a VS test)
#'   Required; A date or date-time object column is expected
#'
#' @param ref_start_date The lower bound of the on-treatment period
#'   Required; A date or date-time object column is expected.
#'
#' @param ref_end_date The upper bound of the on-treatment period
#'   A date or date-time object column is expected.
#'   Optional; This can be null everything after `ref_start_date` is considered
#'   on-treatment.
#'
#' @param ref_end_window A window to add to the upper bound `ref_end_date`.
#'   (e.g. 7 if 7 days should be added to the upper bound)
#'   Optional; default is 0
#'   `ref_end_window` and `ref_end_window_units` should be used together
#'
#' @param ref_end_window_units The units of `ref_end_window`
#'   (e.g. 7 days would have `ref_end_window_units` set to days)
#'   Optional; default is "days", valid values: days, weeks, years
#'   `ref_end_window` and `ref_end_window_units` should be used together
#'
#' @param timepoint_var Timepoint variable used to exclude observations
#'   that occur on the `ref_start_date`.
#'   Optional; default is NA
#'   `timepoint_var` and `start_timepoint_pre_value` should be used together
#'
#' @param start_timepoint_pre_value If a baseline timepoint is used to
#'   differentiate observations taken on the day of treatment, those to be
#'   considered not on-treatment should be listed.
#'   Optional; default is NA
#'   Valid values can be specified as "PRE" or c("PRE", "BEFORE")
#'   `timepoint_var` and `start_timepoint_pre_value` should be used together
#'
#' @details
#' On-Treatment is calculated by determining whether the assessment date or
#' start/stop dates fall between 2 dates. The following logic is used to
#' assign on-treatment = Y:
#'   1. `start_date` is missing and `ref_start_date`is non-missing
#'   2. No timepoint is provided (`timepoint_var`) and both `start_date` and
#'      `ref_start_date` are non-missing and `start_date` = `ref_start_date`
#'   3. A timepoint is provided (`timepoint_var`) and both `start_date` and
#'      `ref_start_date` are non-missing and `start_date` = `ref_start_date` and
#'      `timepoint_var` value is not in `start_timepoint_pre_value` list
#'   4. `ref_end_date` is not provided in the function call and `ref_start_date` <
#'      `start_date`
#'   5. `ref_end_date` is provided in the function call and `ref_start_date` <
#'      `start_date` <= `ref_end_date` + `ref_end_window`.
#'
#' Any date imputations needed should be done prior to calling this function.
#'
#' @author Alice Ehmann
#'
#' @return The input dataset with an additional column named by default
#' `ONTRTFL` with a value of Y or NA

#' @export
#'
#' @examples
#' advs <- tibble::tribble(
#'   ~USUBJID, ~ADT, ~TRTSDT, ~TRTEDT, ~ONTRTFL
#'   "P01",    lubridate::ymd("2020-02-24"), lubridate::ymd("2020-01-01"),
#'             lubridate::ymd("2020-03-01"), "Y",
#'   "P02",    lubridate::ymd("2020-01-01"), lubridate::ymd("2020-01-01"),
#'             lubridate::ymd("2020-03-01"), "Y",
#'   "P03",    lubridate::ymd("2019-12-31"), lubridate::ymd("2020-01-01"),
#'             lubridate::ymd("2020-03-01"), "N")
#' )
#' #' derive_var_ontrtfl(advs, ADT, TRTSDT, TRTEDT)
#'
#' advs <- tibble::tribble(
#'   ~USUBJID, ~ADT, ~TRTSDT, ~TRTEDT, ~ONTRTFL
#'   "P01",    lubridate::ymd("2020-07-01"), lubridate::ymd("2020-01-01"),
#'             lubridate::ymd("2020-03-01"), "N",
#'   "P02",    lubridate::ymd("2020-04-30"), lubridate::ymd("2020-01-01"),
#'             lubridate::ymd("2020-03-01"), "Y",
#'   "P03",    lubridate::ymd("2020-03-15"), lubridate::ymd("2020-01-01"),
#'             lubridate::ymd("2020-03-01"), "Y")
#' )
#' derive_var_ontrtfl(advs, ADT, TRTSDT, TRTEDT, ref_end_window=60,
#'   ref_end_window_units="days")
#'
#' advs <- tibble::tribble(
#' ~USUBJID, ~ADTM, ~TRTSDTM, ~TRTEDTM, ~TPT, ~ONTRTFL,
#' "P01",    lubridate::ymd("2020-01-02T12:00"),
#'           lubridate::ymd_hm("2020-01-01T12:00"),
#'           lubridate::ymd_hm("2020-03-01T12:00"), "", "Y",
#' "P02",    lubridate::ymd("2020-01-01"),
#'           lubridate::ymd_hm("2020-01-01T12:00"),
#'           lubridate::ymd_hm("2020-03-01T12:00"), "PRE", NA,
#' "P03",    lubridate::ymd("2019-12-31"),
#'           lubridate::ymd_hm("2020-01-01T12:00"),
#'           lubridate::ymd_hm("2020-03-01T12:00"), "", NA
#' )
#' derive_var_ontrtfl(advs, ADTM, TRTSDTM, TRTEDTM, timepoint_var=TPT,
#'   start_timepoint_pre_value="PRE")
#'
derive_var_ontrtfl <- function(dataset,
                               start_date,
                               ref_start_date,
                               ref_end_date,
                               ref_end_window = 0,
                               ref_end_window_units = "days",
                               timepoint_var = NA,
                               start_timepoint_pre_value = NA) {

   derive_ontrt(
      dataset,
      new_var = ONTRTFL,
      start_date = !!enquo(start_date),
      ref_start_date = !!enquo(ref_start_date),
      ref_end_present = c(deparse(substitute(ref_end_date))),
      ref_end_date = !!enquo(ref_end_date),
      ref_end_window = ref_end_window,
      ref_end_window_units = ref_end_window_units,
      timepoint_var = !!enquo(timepoint_var),
      start_timepoint_pre_value = start_timepoint_pre_value
   )
}
