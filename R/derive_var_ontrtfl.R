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
#' @param end_date The end date of assessment/event
#'  (e.g. `AENDT` for the date of a CM test)
#'   A date or date-time object column is expected.
#'   Optional; Default is null. If the used and date value is missing
#'   on an observation, it is assumed the medication is ongoing and
#'   ONTRTFL is set to Y.#'
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
#' @param filter_pre_timepoint An expression to filter observations as not
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
#' @author Alice Ehmann, Teckla Akinyi
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
#'   ~USUBJID, ~ADTM,                ~TRTSDTM,                   ~TRTEDTM,                   ~TPT,
#'   "P01", ymd("2020-01-02T12:00"), ymd_hm("2020-01-01T12:00"), ymd_hm("2020-03-01T12:00"), "",
#'   "P02", ymd("2020-01-01"),       ymd_hm("2020-01-01T12:00"), ymd_hm("2020-03-01T12:00"), "PRE",
#'   "P03", ymd("2019-12-31"),       ymd_hm("2020-01-01T12:00"), ymd_hm("2020-03-01T12:00"), ""
#' )
#' derive_var_ontrtfl(
#'   advs,
#'   date = ADTM,
#'   ref_start_date = TRTSDTM,
#'   ref_end_date = TRTEDTM,
#'   filter_pre_timepoint = TPT == "PRE"
#' )
#'
#' advs <- tibble::tribble(
#'   ~USUBJID, ~ASTDT,              ~TRTSDT,           ~TRTEDT,           ~AENDT,
#'   "P01",    ymd("2020-03-15"), ymd("2020-01-01"), ymd("2020-03-01"), ymd("2020-12-01"),
#'   "P02",    ymd("2019-04-30"), ymd("2020-01-01"), ymd("2020-03-01"), ymd("2020-03-15"),
#'   "P03",    ymd("2020-07-01"), ymd("2020-01-01"), ymd("2020-03-01"), ymd("2021-01-01")
#' )
#'
#' derive_var_ontrtfl(
#'  advs,
#'  date = ASTDT,
#'  end_date = AENDT,
#'  ref_start_date = TRTSDT,
#'  ref_end_date = TRTEDT,
#'  ref_end_window = 60
#' )

derive_var_ontrtfl <- function(dataset,
                               start_date,
                               end_date = NULL,
                               ref_start_date,
                               ref_end_date = NULL,
                               ref_end_window = 0,
                               filter_pre_timepoint = NULL,
                               span_period= NULL) {

  start_date <- assert_symbol(enquo(start_date))
  end_date <-  assert_symbol(enquo(end_date), optional = TRUE)
  ref_start_date <- assert_symbol(enquo(ref_start_date))
  ref_end_date <- assert_symbol(enquo(ref_end_date), optional = TRUE)
  ref_end_window <- assert_integer_scalar(ref_end_window, "non-negative")
  filter_pre_timepoint <- assert_filter_cond(enquo(filter_pre_timepoint), optional = TRUE)
  assert_character_scalar(span_period,values=c("Y","y"), optional = TRUE)

  assert_data_frame(
    dataset,
    required_vars = quo_c(start_date, ref_start_date, ref_end_date, end_date)
  )

  dataset <- mutate(
    dataset,
    ONTRTFL = if_else(
      is.na(!!start_date) & !is.na(!!ref_start_date) | !!ref_start_date == !!start_date,
      "Y",
      NA_character_,
      missing = NA_character_)
  )

  if (!quo_is_null(filter_pre_timepoint)) {
    dataset <- mutate(
      dataset,
      ONTRTFL = if_else(!!filter_pre_timepoint, NA_character_, ONTRTFL)
    )
  }

  if (quo_is_null(ref_end_date)) {
    # Scenario 1: No treatment end date is passed
    dataset <- mutate(
      dataset,
      ONTRTFL = if_else(
        !is.na(!!ref_start_date) & !is.na(!!start_date) & !!ref_start_date < !!start_date,
        "Y",
        ONTRTFL)
    )
  } else {
    # Scenario 2: Treatment end date is passed, window added above
    dataset <- mutate(
      dataset,
      ONTRTFL = if_else(
        !is.na(!!ref_start_date) & !is.na(!!start_date) & !!ref_start_date < !!start_date &
          !is.na(!!ref_end_date) & !!start_date <= (!!ref_end_date + days(!!ref_end_window)),
        "Y",
        ONTRTFL
      )
    )
  }

    #scenario 3: end_date is parsed
  if (!quo_is_null(end_date)) {
    dataset <- mutate(
      dataset,
      ONTRTFL = if_else(!!end_date >= !!ref_start_date,
                        "Y",
                        ONTRTFL
                        )
                      )
  }

  #scenario 4: end_date and span_period are parsed
  if (span_period) {
    dataset <- mutate(
      dataset,
      ONTRTFL = if_else(
        !is.na(!!start_date)  &
          !!start_date <= (!!ref_end_date + days(!!ref_end_window)) &
          is.na(!!end_date) |
          !!end_date >= !!ref_start_date,
        "Y",
        ONTRTFL
      )
    )
  }

  dataset
}



bds4 <- tibble::tribble(
  ~USUBJID, ~ASTDT,              ~TRTSDT,           ~TRTEDT,           ~AENDT, ~TPT,
  "P01",    ymd("2020-03-01"), ymd("2020-01-01"), ymd("2020-03-01"), ymd("2020-12-01"), NA,
  "P02",    ymd("2019-04-30"), ymd("2020-01-01"), ymd("2020-03-01"), ymd("2020-03-15"), NA,
  "P03",    NA,                 ymd("2020-01-01"), ymd("2020-03-01"), ymd("2020-12-01"), NA, #endate is after trtsdt & start is NA
  "P03",    NA,                 ymd("2020-01-01"), ymd("2020-03-01"), ymd("2019-12-01"), NA, #endate is before trtsdt & start is NA
  "P03",    NA,                 ymd("2020-01-01"), ymd("2020-03-01"), ymd("2020-04-01"), NA, #endate is after trtsdt & start is NA
  "try3",     ymd("2019-01-01"),ymd("2020-01-01"), ymd("2020-03-01"), NA, NA, #endate is NA & start is before trtsdt
  "PAT04",  ymd("2020-03-01"), ymd("2020-01-01"), ymd("2021-02-01"), ymd("2021-01-01"), NA,
  "PAT01",  ymd("2020-01-01"), ymd("2020-01-01"),ymd("2020-03-01"), ymd("2020-03-01"), "PRE",
  "PAT02",  ymd("2020-01-01"), ymd("2020-01-01"), ymd("2020-03-01") ,ymd("2020-03-01"), "POST"


)

derive_var_ontrtfl(
  bds4,
  start_date =  ASTDT,
  end_date =  AENDT,
  ref_start_date = TRTSDT,
  ref_end_date = TRTEDT,
  ref_end_window = 160,
  filter_pre_timepoint = TPT == "PRE",
  span_period="y"
)

