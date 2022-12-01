#' Derive On-Treatment Flag Variable
#'
#' Derive on-treatment flag (`ONTRTFL`) in an ADaM dataset with a single
#' assessment date (e.g `ADT`) or event start and end dates (e.g.
#' `ASTDT`/`AENDT`).
#'
#' @param dataset Input dataset.
#'
#'   Required columns are `start_date`, `end_date`, `ref_start_date` and
#'   `ref_end_date`.
#'
#' @param new_var On-treatment flag variable name to be created.
#'
#'   Default is `ONTRTFL`.
#'
#' @param start_date The start date (e.g. `AESDT`) or assessment date (e.g.
#'   `ADT`) Required; A date or date-time object column is expected.
#'
#'   Refer to `derive_vars_dt()` to impute and derive a date from a date
#'   character vector to a date object.
#'
#' @param end_date The end date of assessment/event (e.g. `AENDT`) A date or
#'   date-time object column is expected.
#'
#'   Refer to `derive_vars_dt()` to impute and derive a date from a date
#'   character vector to a date object.
#'
#'   Optional; Default is null. If the used and date value is missing on an
#'   observation, it is assumed the medication is ongoing and `ONTRTFL` is set
#'   to `"Y"`.
#'
#' @param ref_start_date The lower bound of the on-treatment period Required; A
#'   date or date-time object column is expected.
#'
#'   Refer to `derive_vars_dt()` to impute and derive a date from a date
#'   character vector to a date object.
#'
#' @param ref_end_date The upper bound of the on-treatment period A date or
#'   date-time object column is expected.
#'
#'   Refer to `derive_vars_dt()` to impute and derive a date from a date
#'   character vector to a date object.
#'
#'   Optional; This can be null and everything after `ref_start_date` will be
#'   considered on-treatment. Default is `NULL`.
#'
#' @param ref_end_window A window to add to the upper bound `ref_end_date`
#'   measured in days (e.g. 7 if 7 days should be added to the upper bound)
#'   Optional; default is 0.
#'
#' @param ignore_time_for_ref_end_date
#'
#'   If the argument is set to `TRUE`, the time part is ignored for checking if
#'   the event occurred more than `ref_end_window` days after reference end
#'   date.
#'
#'   *Permitted Values:* `TRUE`, `FALSE`
#'
#' @param filter_pre_timepoint An expression to filter observations as not
#'   on-treatment when `date` = `ref_start_date`. For example, if observations
#'   where `VSTPT = PRE` should not be considered on-treatment when `date =
#'   ref_start_date`, `filter_pre_timepoint` should be used to denote when the
#'   on-treatment flag should be set to null. Optional; default is `NULL`.
#'
#' @param span_period A `"Y"` scalar character. If `"Y"`, events that started
#'   prior to the `ref_start_date`and are ongoing or end after the
#'   `ref_start_date` are flagged as `"Y"`. Optional; default is `NULL`.
#'
#' @details On-Treatment is calculated by determining whether the assessment
#'   date or start/stop dates fall between 2 dates. The following logic is used
#'   to assign on-treatment = `"Y"`:
#'   1. `start_date` is missing and `ref_start_date`is non-missing
#'   2. No timepoint filter is provided (`filter_pre_timepoint`) and both
#'   `start_date` and `ref_start_date` are non-missing and `start_date` =
#'   `ref_start_date`
#'   3. A timepoint is provided (`filter_pre_timepoint`) and both `start_date`
#'   and `ref_start_date` are non-missing and `start_date = ref_start_date` and
#'   the filter provided in `filter_pre_timepoint` is not true.
#'   4. `ref_end_date` is not provided and `ref_start_date < start_date`
#'   5. `ref_end_date` is provided and `ref_start_date < start_date` <=
#'   `ref_end_date + ref_end_window`.
#'
#'   If the `end_date` is provided and the `end_date` < ref_start_date then the
#'   `ONTRTFL` is set to `NULL`.This would be applicable to cases where the
#'   `start_date` is missing and `ONTRTFL` has been assigned as `"Y"` above.
#'
#'   If the `span_period` is specified as `"Y"`, this allows the user to assign
#'   `ONTRTFL` as `"Y"` to cases where the record started prior to the
#'   `ref_start_date` and was ongoing or ended after the `ref_start_date`.
#'
#'   Any date imputations needed should be done prior to calling this function.
#'
#' @author Alice Ehmann, Teckla Akinyi
#'
#' @family der_bds_findings
#' @keywords der_bds_findings
#'
#' @return The input dataset with an additional column named `ONTRTFL` with a
#'   value of `"Y"` or `NA`
#'
#' @export
#'
#' @examples
#' library(tibble)
#' library(dplyr, warn.conflicts = FALSE)
#' library(lubridate, warn.conflicts = FALSE)
#'
#' advs <- tribble(
#'   ~USUBJID, ~ADT,              ~TRTSDT,           ~TRTEDT,
#'   "P01",    ymd("2020-02-24"), ymd("2020-01-01"), ymd("2020-03-01"),
#'   "P02",    ymd("2020-01-01"), ymd("2020-01-01"), ymd("2020-03-01"),
#'   "P03",    ymd("2019-12-31"), ymd("2020-01-01"), ymd("2020-03-01")
#' )
#' derive_var_ontrtfl(
#'   advs,
#'   start_date = ADT,
#'   ref_start_date = TRTSDT,
#'   ref_end_date = TRTEDT
#' )
#'
#' advs <- tribble(
#'   ~USUBJID, ~ADT,              ~TRTSDT,           ~TRTEDT,
#'   "P01",    ymd("2020-07-01"), ymd("2020-01-01"), ymd("2020-03-01"),
#'   "P02",    ymd("2020-04-30"), ymd("2020-01-01"), ymd("2020-03-01"),
#'   "P03",    ymd("2020-03-15"), ymd("2020-01-01"), ymd("2020-03-01")
#' )
#' derive_var_ontrtfl(
#'   advs,
#'   start_date = ADT,
#'   ref_start_date = TRTSDT,
#'   ref_end_date = TRTEDT,
#'   ref_end_window = 60
#' )
#'
#' advs <- tribble(
#'   ~USUBJID, ~ADTM,                      ~TRTSDTM,                   ~TRTEDTM,
#'   "P01",    ymd_hm("2020-01-02T12:00"), ymd_hm("2020-01-01T12:00"), ymd_hm("2020-03-01T12:00"),
#'   "P02",    ymd("2020-01-01"),          ymd_hm("2020-01-01T12:00"), ymd_hm("2020-03-01T12:00"),
#'   "P03",    ymd("2019-12-31"),          ymd_hm("2020-01-01T12:00"), ymd_hm("2020-03-01T12:00"),
#' ) %>%
#'   mutate(TPT = c(NA, "PRE", NA))
#' derive_var_ontrtfl(
#'   advs,
#'   start_date = ADTM,
#'   ref_start_date = TRTSDTM,
#'   ref_end_date = TRTEDTM,
#'   filter_pre_timepoint = TPT == "PRE"
#' )
#'
#' advs <- tribble(
#'   ~USUBJID, ~ASTDT,            ~TRTSDT,           ~TRTEDT,           ~AENDT,
#'   "P01",    ymd("2020-03-15"), ymd("2020-01-01"), ymd("2020-03-01"), ymd("2020-12-01"),
#'   "P02",    ymd("2019-04-30"), ymd("2020-01-01"), ymd("2020-03-01"), ymd("2020-03-15"),
#'   "P03",    ymd("2019-04-30"), ymd("2020-01-01"), ymd("2020-03-01"), NA,
#' )
#' derive_var_ontrtfl(
#'   advs,
#'   start_date = ASTDT,
#'   end_date = AENDT,
#'   ref_start_date = TRTSDT,
#'   ref_end_date = TRTEDT,
#'   ref_end_window = 60,
#'   span_period = "Y"
#' )
#'
#' advs <- tribble(
#'   ~USUBJID, ~ASTDT,            ~AP01SDT,          ~AP01EDT,          ~AENDT,
#'   "P01",    ymd("2020-03-15"), ymd("2020-01-01"), ymd("2020-03-01"), ymd("2020-12-01"),
#'   "P02",    ymd("2019-04-30"), ymd("2020-01-01"), ymd("2020-03-01"), ymd("2020-03-15"),
#'   "P03",    ymd("2019-04-30"), ymd("2020-01-01"), ymd("2020-03-01"), NA,
#' )
#' derive_var_ontrtfl(
#'   advs,
#'   new_var = ONTR01FL,
#'   start_date = ASTDT,
#'   end_date = AENDT,
#'   ref_start_date = AP01SDT,
#'   ref_end_date = AP01EDT,
#'   span_period = "Y"
#' )
derive_var_ontrtfl <- function(dataset,
                               new_var = ONTRTFL,
                               start_date,
                               end_date = NULL,
                               ref_start_date,
                               ref_end_date = NULL,
                               ref_end_window = 0,
                               ignore_time_for_ref_end_date = TRUE,
                               filter_pre_timepoint = NULL,
                               span_period = NULL) {
  new_var <- assert_symbol(enquo(new_var))
  start_date <- assert_symbol(enquo(start_date))
  end_date <- assert_symbol(enquo(end_date), optional = TRUE)
  ref_start_date <- assert_symbol(enquo(ref_start_date))
  ref_end_date <- assert_symbol(enquo(ref_end_date), optional = TRUE)
  assert_data_frame(
    dataset,
    required_vars = quo_c(start_date, end_date, ref_start_date, ref_end_date)
  )
  warn_if_vars_exist(dataset, quo_text(new_var))

  ref_end_window <- assert_integer_scalar(ref_end_window, "non-negative")
  assert_logical_scalar(ignore_time_for_ref_end_date)
  filter_pre_timepoint <- assert_filter_cond(enquo(filter_pre_timepoint), optional = TRUE)
  assert_character_scalar(span_period, values = c("Y", "y"), optional = TRUE)

  dataset <- mutate(
    dataset,
    !!new_var := if_else(
      is.na(!!start_date) & !is.na(!!ref_start_date) | !!ref_start_date == !!start_date,
      "Y",
      NA_character_,
      missing = NA_character_
    )
  )

  if (!quo_is_null(filter_pre_timepoint)) {
    dataset <- mutate(
      dataset,
      !!new_var := if_else(!!filter_pre_timepoint, NA_character_, !!new_var,
        missing = !!new_var
      )
    )
  }

  if (quo_is_null(ref_end_date)) {
    # Scenario 1: No treatment end date is passed
    dataset <- mutate(
      dataset,
      !!new_var := if_else(
        !is.na(!!ref_start_date) & !is.na(!!start_date) & !!ref_start_date < !!start_date,
        "Y",
        !!new_var,
        missing = !!new_var
      )
    )
  } else {
    # Scenario 2: Treatment end date is passed, window added above
    if (ignore_time_for_ref_end_date) {
      end_cond <- expr(date(!!start_date) <= date(!!ref_end_date) + days(!!ref_end_window))
    } else {
      end_cond <- expr(!!start_date <= !!ref_end_date + days(!!ref_end_window))
    }
    dataset <- mutate(
      dataset,
      !!new_var := if_else(
        !!ref_start_date < !!start_date & !!end_cond,
        "Y",
        !!new_var,
        missing = !!new_var
      )
    )
  }

  # scenario 3: end_date is passed
  if (!quo_is_null(end_date)) {
    dataset <- mutate(
      dataset,
      !!new_var := if_else(
        !!end_date < !!ref_start_date,
        NA_character_,
        !!new_var,
        missing = !!new_var
      )
    )
  }

  # scenario 4: end_date and span_period are passed
  if (!is.null(span_period)) {
    dataset <- mutate(
      dataset,
      !!new_var := if_else(
        !!end_cond &
          (is.na(!!end_date) | !!end_date >= !!ref_start_date),
        "Y",
        !!new_var,
        missing = !!new_var
      )
    )
  }

  dataset
}
