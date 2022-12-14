#' Derive Treatment-emergent Flag
#'
#' Derive treatment emergent analysis flag (e.g., `TRTEMFL`).
#'
#' @param dataset Input dataset
#'
#'   The variables specified by `start_date`, `end_date`, `trt_start_date`,
#'   `trt_end_date`, `initial_intensity`, and `intensity` are expected.
#'
#' @param new_var New variable
#'
#' @param start_date Event start date
#'
#'   *Permitted Values:* A symbol referring to a date or datetime variable of
#'   the input dataset
#'
#' @param end_date Event end date
#'
#'   *Permitted Values:* A symbol referring to a date or datetime variable of
#'   the input dataset
#'
#' @param trt_start_date Treatment start date
#'
#'   *Permitted Values:* A symbol referring to a date or datetime variable of
#'   the input dataset
#'
#' @param trt_end_date Treatment end date
#'
#'   *Permitted Values:* A symbol referring to a date or datetime variable of
#'   the input dataset or `NULL`
#'
#' @param end_window
#'
#'    If the argument is specified, events starting more than the specified
#'    number of days after end of treatment, are not flagged.
#'
#'    *Permitted Values:* A non-negative integer or `NULL`
#'
#' @param ignore_time_for_trt_end
#'
#'   If the argument is set to `TRUE`, the time part is ignored for checking if
#'   the event occurred more than `end_window` days after end of treatment.
#'
#'   *Permitted Values:* `TRUE`, `FALSE`
#'
#' @param initial_intensity Initial severity/intensity or toxicity
#'
#'   This derivation assumes AE data collection method as single record per AE
#'   with “initial” and “most extreme” severity/intensity recorded separately.
#'
#'   If the argument is specified, events which start before treatment start and
#'   end after treatment start (or are ongoing) and worsened (i.e., the
#'   intensity is greater than the initial intensity), are flagged.
#'
#'   The values of the specified variable must be comparable with the usual
#'   comparison operators. I.e., if the intensity is greater than the initial
#'   intensity `initial_intensity < intensity` must evaluate to `TRUE`.
#'
#'   *Permitted Values:* A symbol referring to a variable of the input dataset
#'   or `NULL`
#'
#' @param intensity Severity/intensity or toxicity
#'
#'   If the argument is specified, events which start before treatment start and
#'   end after treatment start (or are ongoing) and worsened (i.e., the
#'   intensity is greater than the initial intensity), are flagged.
#'
#'   The values of the specified variable must be comparable with the usual
#'   comparison operators. I.e., if the intensity is greater than the initial
#'   intensity `initial_intensity < intensity` must evaluate to `TRUE`.
#'
#'   *Permitted Values:* A symbol referring to a variable of the input dataset
#'   or `NULL`
#'
#' @author Stefan Bundfuss
#'
#' @details For the derivation of the new variable the following cases are
#'   considered in this order. The first case which applies, defines the value
#'   of the variable.
#'
#' - *not treated*: If `trt_start_date` is `NA`, it is set to `NA_character_`.
#' - *event before treatment*: If `end_date` is before `trt_start_date` (and
#' `end_date` is not `NA`), it is set to `NA_character_`.
#' - *no event date*: If `start_date` is `NA`, it is set to `"Y"` as in such
#' cases it is usually considered more conservative to assume the event was
#' treatment-emergent.
#' - *event started during treatment*:
#'     - if `end_window` is not specified:
#'       if `start_date` is on or after `trt_start_date`, it is set to `"Y"`,
#'     - if `end_window` is specified:
#'       if `start_date` is on or after `trt_start_date` and `start_date` is on
#'       or before `trt_end_date` + `end_window` days, it is set to `"Y"`,
#'  - *event started before treatment and (possibly) worsened on treatment*:
#'    - if `initial_intensity` and `intensity` is specified: if
#'    `initial_intensity < intensity` and `start_date` is before
#'    `trt_start_date` and `end_date` is on or after `trt_start_date` or
#'    `end_date` is `NA`, it is set to `"Y"`.
#'  - Otherwise it is set to `NA_character_`.
#'
#' @return The input dataset with the variable specified by `new_var` added
#'
#' @keywords der_occds
#' @family der_occds
#'
#' @export
#'
#' @examples
#'
#' library(tibble)
#' library(dplyr, warn.conflicts = FALSE)
#' library(lubridate)
#'
#' adae <- expected <- tribble(
#'   ~USUBJID, ~ASTDTM,            ~AENDTM,            ~AEITOXGR, ~AETOXGR,
#'   # before treatment
#'   "1",      "2021-12-13T20:15", "2021-12-15T12:45", "1",       "1",
#'   "1",      "2021-12-14T20:15", "2021-12-14T22:00", "1",       "3",
#'   # starting before treatment and ending during treatment
#'   "1",      "2021-12-30T20:00", "2022-01-14T11:00", "1",       "3",
#'   "1",      "2021-12-31T20:15", "2022-01-01T01:23", "1",       "1",
#'   # starting during treatment
#'   "1",      "2022-01-01T12:00", "2022-01-02T23:25", "3",       "4",
#'   # after treatment
#'   "1",      "2022-05-10T11:00", "2022-05-10T13:05", "2",       "2",
#'   "1",      "2022-05-11T11:00", "2022-05-11T13:05", "2",       "2",
#'   # missing dates
#'   "1",      "",                 "",                 "3",       "4",
#'   "1",      "2021-12-30T09:00", "",                 "3",       "4",
#'   "1",      "2021-12-30T11:00", "",                 "3",       "3",
#'   "1",      "",                 "2022-01-04T09:00", "3",       "4",
#'   "1",      "",                 "2021-12-24T19:00", "3",       "4",
#'   "1",      "",                 "2022-06-04T09:00", "3",       "4",
#'   # without treatment
#'   "2",      "",                 "2021-12-03T12:00", "1",       "2",
#'   "2",      "2021-12-01T12:00", "2021-12-03T12:00", "1",       "2",
#'   "2",      "2021-12-06T18:00", "",                 "1",       "2"
#' ) %>%
#'   mutate(
#'     ASTDTM = ymd_hm(ASTDTM),
#'     AENDTM = ymd_hm(AENDTM),
#'     TRTSDTM = if_else(USUBJID == "1", ymd_hm("2022-01-01T01:01"), ymd_hms("")),
#'     TRTEDTM = if_else(USUBJID == "1", ymd_hm("2022-04-30T23:59"), ymd_hms(""))
#'   )
#'
#' # derive TRTEMFL without considering treatment end and worsening
#' derive_var_trtemfl(adae) %>% select(ASTDTM, AENDTM, TRTSDTM, TRTEMFL)
#'
#' # derive TRTEM2FL taking treatment end and worsening into account
#' derive_var_trtemfl(
#'   adae,
#'   new_var = TRTEM2FL,
#'   trt_end_date = TRTEDTM,
#'   end_window = 10,
#'   initial_intensity = AEITOXGR,
#'   intensity = AETOXGR
#' ) %>% select(ASTDTM, AENDTM, AEITOXGR, AETOXGR, TRTEM2FL)
derive_var_trtemfl <- function(dataset,
                               new_var = TRTEMFL,
                               start_date = ASTDTM,
                               end_date = AENDTM,
                               trt_start_date = TRTSDTM,
                               trt_end_date = NULL,
                               end_window = NULL,
                               ignore_time_for_trt_end = TRUE,
                               initial_intensity = NULL,
                               intensity = NULL) {
  new_var <- assert_symbol(enquo(new_var))
  start_date <- assert_symbol(enquo(start_date))
  end_date <- assert_symbol(enquo(end_date))
  trt_start_date <- assert_symbol(enquo(trt_start_date))
  trt_end_date <- assert_symbol(enquo(trt_end_date), optional = TRUE)
  assert_integer_scalar(end_window, subset = "non-negative", optional = TRUE)
  assert_logical_scalar(ignore_time_for_trt_end)
  initial_intensity <- assert_symbol(enquo(initial_intensity), optional = TRUE)
  intensity <- assert_symbol(enquo(intensity), optional = TRUE)
  if (quo_is_null(initial_intensity) && !quo_is_null(intensity)) {
    abort(paste(
      "`intensity` argument was specified but not `initial_intensity`",
      "Either both or none of them must be specified.",
      sep = "\n"
    ))
  }
  if (!quo_is_null(initial_intensity) && quo_is_null(intensity)) {
    abort(paste(
      "`initial_intensity` argument was specified but not `intensity`",
      "Either both or none of them must be specified.",
      sep = "\n"
    ))
  }
  assert_data_frame(
    dataset,
    required_vars = quo_c(
      start_date,
      end_date,
      trt_start_date,
      trt_end_date,
      initial_intensity,
      intensity
    )
  )
  assert_date_var(dataset, var = !!start_date)
  assert_date_var(dataset, var = !!end_date)
  assert_date_var(dataset, var = !!trt_start_date)
  if (!quo_is_null(trt_end_date)) {
    assert_date_var(dataset, var = !!trt_end_date)
  }

  if (is.null(end_window)) {
    end_cond <- expr(TRUE)
  } else {
    if (quo_is_null(trt_end_date)) {
      abort(paste(
        "`end_window` argument was specified but not `trt_end_date`",
        "Either both or none of them must be specified.",
        sep = "\n"
      ))
    }
    if (ignore_time_for_trt_end) {
      end_cond <- expr(is.na(!!trt_end_date) |
        date(!!start_date) <= date(!!trt_end_date) + days(end_window))
    } else {
      end_cond <- expr(is.na(!!trt_end_date) | !!start_date <= !!trt_end_date + days(end_window))
    }
  }

  if (quo_is_null(intensity)) {
    worsening_cond <- expr(FALSE)
  } else {
    worsening_cond <-
      expr(!!start_date < !!trt_start_date &
        (!!initial_intensity < !!intensity | is.na(!!initial_intensity) | is.na(!!intensity)))
  }

  dataset %>%
    mutate(!!new_var := case_when(
      is.na(!!trt_start_date) ~ NA_character_,
      !!end_date < !!trt_start_date ~ NA_character_,
      is.na(!!start_date) ~ "Y",
      !!start_date >= !!trt_start_date & !!end_cond ~ "Y",
      !!worsening_cond ~ "Y"
    ))
}
