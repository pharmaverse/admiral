#' Derive On-Treatment Flag Variable
#'
#' Derive on-treatment flag (`ONTRTFL`) in an ADaM dataset with a single
#' assessment date
#'
#' @param dataset `data.frame`.
#'
#' @param date The assessment date (e.g. `ADT` for the date of a VS test)
#'   Required; A date or date-time object column is expected
#'
#' @param ref_start_date The lower bound of the on-treatment period
#'   Required; A date or date-time object column is expected.
#'
#' @param ref_end_date The upper bound of the on-treatment period
#'   A date or date-time object column is expected.
#'   Optional; This can be null and everything after `ref_start_date` will be
#'   considered on-treatment.
#'   Default is NULL.
#'
#' @param ref_end_window A window to add to the upper bound `ref_end_date`
#'   measured in days
#'   (e.g. 7 if 7 days should be added to the upper bound)
#'   Optional; default is 0
#'
#' @param filter_pre_timepoint An expression to filter observtions as not
#' on-treatment when `date` = `ref_start_date`. For example, if
#' observations where `VSTPT = PRE` should not be considered on-treatment when
#' `date = ref_start_date`, `filter_pre_timepoint` should be used
#' to denote when the on-treatment flag should be set to null.
#' Optional; default = NULL
#'
#' @details
#' On-Treatment is calculated by determining whether the assessment date or
#' start/stop dates fall between 2 dates. The following logic is used to
#' assign on-treatment = Y:
#'   1. `date` is missing and `ref_start_date`is non-missing
#'   2. No timepoint filter is provided (`filter_pre_timepoint`) and both
#'      `date` and `ref_start_date` are non-missing and `date` =
#'      `ref_start_date`
#'   3. A timepoint is provided (`filter_pre_timepoint`) and both `date`
#'      and `ref_start_date` are non-missing and `date = ref_start_date`
#'      and the filter provided in `filter_pre_timepoint` is not true.
#'   4. `ref_end_date` is not provided and `ref_start_date < date`
#'   5. `ref_end_date` is provided and `ref_start_date < date` <=
#'      `ref_end_date + ref_end_window`.
#'
#' Any date imputations needed should be done prior to calling this function.
#'
#' @author Alice Ehmann
#'
#' @keywords bds derivation
#'
#' @return The input dataset with an additional column named
#' `ONTRTFL` with a value of `"Y"` or `NA`
#'
#' @export
#'
#' @examples
#' library(lubridate, warn.conflict = FALSE)
#'
#' advs <- tibble::tribble(
#'   ~USUBJID, ~ADT,              ~TRTSDT,           ~TRTEDT,
#'   "P01",    ymd("2020-02-24"), ymd("2020-01-01"), ymd("2020-03-01"),
#'   "P02",    ymd("2020-01-01"), ymd("2020-01-01"), ymd("2020-03-01"),
#'   "P03",    ymd("2019-12-31"), ymd("2020-01-01"), ymd("2020-03-01")
#' )
#' derive_var_ontrtfl(
#'   advs,
#'   date = ADT,
#'   ref_start_date = TRTSDT,
#'   ref_end_date = TRTEDT
#' )
#'
#' advs <- tibble::tribble(
#'   ~USUBJID, ~ADT,              ~TRTSDT,           ~TRTEDT,
#'   "P01",    ymd("2020-07-01"), ymd("2020-01-01"), ymd("2020-03-01"),
#'   "P02",    ymd("2020-04-30"), ymd("2020-01-01"), ymd("2020-03-01"),
#'   "P03",    ymd("2020-03-15"), ymd("2020-01-01"), ymd("2020-03-01")
#' )
#' derive_var_ontrtfl(
#'   advs,
#'   date = ADT,
#'   ref_start_date = TRTSDT,
#'   ref_end_date = TRTEDT,
#'   ref_end_window = 60
#' )
#'
#' advs <- tibble::tribble(
#'   ~USUBJID, ~ADTM,                   ~TRTSDTM,                   ~TRTEDTM,                   ~TPT,
#'   "P01",    ymd("2020-01-02T12:00"), ymd_hm("2020-01-01T12:00"), ymd_hm("2020-03-01T12:00"), "",
#'   "P02",    ymd("2020-01-01"),       ymd_hm("2020-01-01T12:00"), ymd_hm("2020-03-01T12:00"), "PRE",
#'   "P03",    ymd("2019-12-31"),       ymd_hm("2020-01-01T12:00"), ymd_hm("2020-03-01T12:00"), ""
#' )
#' derive_var_ontrtfl(
#'   advs,
#'   date = ADTM,
#'   ref_start_date = TRTSDTM,
#'   ref_end_date = TRTEDTM,
#'   filter_pre_timepoint = expr(TPT == "PRE")
#' )
derive_var_ontrtfl <- function(dataset,
                               date,
                               ref_start_date,
                               ref_end_date = NULL,
                               ref_end_window = 0,
                               filter_pre_timepoint = NULL) {

   new_var <- expr(ONTRTFL)

   warn_if_vars_exist(
      dataset,
      c(deparse(new_var))
   )

   assert_that(
      is.data.frame(dataset),
      is.numeric(ref_end_window)
   )
   assert_has_variables(dataset, c(deparse(substitute(date)),
                                   deparse(substitute(ref_start_date))))

   date <- enquo(date)
   ref_start_date <- enquo(ref_start_date)
   ref_end_date <- enquo(ref_end_date)

   dataset <- dataset %>%
      mutate(!!new_var :=
                case_when(is.na(!!date) & !is.na(!!ref_start_date) ~ "Y",
                          !is.na(!!date) & !is.na(!!ref_start_date) &
                             !!ref_start_date == !!date ~ "Y"))

   if (!is.null(filter_pre_timepoint)) {
      dataset <- dataset  %>%
         mutate(!!new_var := if_else(!!!filter_pre_timepoint, NA_character_, !!new_var))
   }

   if (quo_is_null(ref_end_date)) {
      #Scenario 1: No treatment end date is passed
      dataset <- dataset %>%
         mutate(!!new_var := if_else(!is.na(!!ref_start_date) &
                                        !is.na(!!date) &
                                        !!ref_start_date < !!date,
                                     "Y", !!new_var))
   } else {
      #Scenario 2: Treatment end date is passed, window added above
      dataset <- dataset %>%
         mutate(!!new_var :=
                   if_else(!is.na(!!ref_start_date) &
                              !is.na(!!date) & !!ref_start_date < !!date &
                              !is.na(!!ref_end_date) &
                              !!date <= (!!ref_end_date + days(x = !!ref_end_window)), "Y", !!new_var))
   }

   dataset

}
