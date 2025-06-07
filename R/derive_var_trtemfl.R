#' Derive Treatment-emergent Flag
#'
#' Derive treatment emergent analysis flag (e.g., `TRTEMFL`).
#'
#' @param dataset `r roxygen_param_dataset()`
#'
#'   The variables specified by `start_date`, `end_date`, `trt_start_date`,
#'   `trt_end_date`, `initial_intensity`, and `intensity` are expected.
#'
#' @permitted [dataset]
#'
#' @param new_var New variable
#'
#' @permitted [var]
#'
#' @param start_date Event start date
#'
#' @permitted [date]
#'
#' @param end_date Event end date
#'
#' @permitted [date]
#'
#' @param trt_start_date Treatment start date
#'
#' @permitted [date]
#'
#' @param trt_end_date Treatment end date
#'
#' @permitted [date]
#'
#' @param end_window
#'
#'    If the argument is specified (in 'days'), events starting more than the specified
#'    number of days after end of treatment, are not flagged.
#'
#' @permitted [pos_int]
#'
#' @param ignore_time_for_trt_end
#'
#'   If the argument is set to `TRUE`, the time part is ignored for checking if
#'   the event occurred more than `end_window` days after end of treatment.
#'
#' @permitted [boolean]
#'
#' @param initial_intensity Initial severity/intensity or toxicity
#'
#'   `initial_intensity` is ignored when `group_var` is specified.
#'
#'   If this argument is specified and `group_var` is `NULL`, events which start
#'   before treatment start and end after treatment start (or are ongoing) and
#'   worsened (i.e., the intensity is greater than the initial intensity), are
#'   flagged.
#'
#'   The values of the specified variable must be comparable with the usual
#'   comparison operators. I.e., if the intensity is greater than the initial
#'   intensity `initial_intensity < intensity` must evaluate to `TRUE`.
#'
#' @permitted [var]
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
#' @permitted [var]
#'
#' @param group_var Grouping variable
#'
#'  If the argument is specified, it assumes that AEs are recorded as one episode
#'  of AE with multiple lines using a grouping variable.
#'
#'  Events starting during treatment or before treatment and worsening afterward
#'  are flagged. Once an AE record in a group is flagged, all subsequent records
#'  in the treatment window are flagged regardless of severity.
#'
#' @permitted [var]
#'
#' @param subject_keys Variables to uniquely identify a subject.
#'
#'   This argument is only used when `group_var` is specified.
#'
#' @permitted [var_list]
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
#'    - if `initial_intensity`, `intensity` is specified and `group_var` is not specified:
#'      if `initial_intensity < intensity` and `start_date` is before `trt_start_date`
#'      and `end_date` is on or after `trt_start_date` or `end_date` is `NA`, it
#'      is set to `"Y"`;
#'    - if `group_var` is specified:
#'      if `intensity` at treatment start < `intensity` and `start_date` is after
#'      `trt_start_date` and `end_date` is on or after `trt_start_date` or
#'      `end_date` is `NA`, it is set to `"Y"`;
#'
#'  - Otherwise it is set to `NA_character_`.
#'
#'  The behavior of `derive_var_trtemfl()` is aligned with the proposed
#'  treatment-emergent AE assignment in the following
# nolint start
#'  [PHUSE White Paper](https://phuse.s3.eu-central-1.amazonaws.com/Deliverables/Safety+Analytics/WP-087+Recommended+Definition+of++Treatment-Emergent+Adverse+Events+in+Clinical+Trials+.pdf).
# nolint end
#'  See the final example in the examples section below.
#'
#' @return The input dataset with the variable specified by `new_var` added
#'
#' @keywords der_occds
#' @family der_occds
#'
#' @export
#'
#' @examplesx
#'
#' @caption Basic treatment-emergent flag
#' @info Derive `TRTEMFL` without considering treatment end and worsening
#'
#' - For this basic example, all we are using are AE start/end dates and
#'   comparing those against treatment start date.
#' - If the AE started on or after treatment then we flag as treatment-emergent
#'   (e.g. records 5-7).
#' - If missing AE start date then we flag as treatment-emergent as worst case
#'   (e.g. records 8, 11 and 13), unless we know that the AE end date was before
#'   treatment so we can rule out this being treatment-emergent (e.g. record 12).
#' - Any not treated subject would not get their AEs flagged as treatment-emergent
#'   (e.g. records 14-16).
#' @code
#' library(tibble)
#' library(dplyr, warn.conflicts = FALSE)
#' library(lubridate)
#'
#' adae <- tribble(
#'   ~USUBJID, ~ASTDT,            ~AENDT,            ~AEITOXGR, ~AETOXGR,
#'   # before treatment
#'   "1",      ymd("2021-12-13"), ymd("2021-12-15"), "1",       "1",
#'   "1",      ymd("2021-12-14"), ymd("2021-12-14"), "1",       "3",
#'   # starting before treatment and ending during treatment
#'   "1",      ymd("2021-12-30"), ymd("2022-01-14"), "1",       "3",
#'   "1",      ymd("2021-12-31"), ymd("2022-01-01"), "1",       "1",
#'   # starting during treatment
#'   "1",      ymd("2022-01-01"), ymd("2022-01-02"), "3",       "4",
#'   # after treatment
#'   "1",      ymd("2022-05-10"), ymd("2022-05-10"), "2",       "2",
#'   "1",      ymd("2022-05-11"), ymd("2022-05-11"), "2",       "2",
#'   # missing dates
#'   "1",      NA,                NA,                "3",       "4",
#'   "1",      ymd("2021-12-30"), NA,                "3",       "4",
#'   "1",      ymd("2021-12-31"), NA,                "3",       "3",
#'   "1",      NA,                ymd("2022-01-04"), "3",       "4",
#'   "1",      NA,                ymd("2021-12-24"), "3",       "4",
#'   "1",      NA,                ymd("2022-06-04"), "3",       "4",
#'   # without treatment
#'   "2",      NA,                ymd("2021-12-03"), "1",       "2",
#'   "2",      ymd("2021-12-01"), ymd("2021-12-03"), "1",       "2",
#'   "2",      ymd("2021-12-06"), NA,                "1",       "2"
#' ) %>%
#'   mutate(
#'     STUDYID = "AB42",
#'     TRTSDT = if_else(USUBJID == "1", ymd("2022-01-01"), NA),
#'     TRTEDT = if_else(USUBJID == "1", ymd("2022-04-30"), NA)
#'   )
#'
#' derive_var_trtemfl(
#'   adae,
#'   start_date = ASTDT,
#'   end_date = AENDT,
#'   trt_start_date = TRTSDT
#' ) %>% select(USUBJID, TRTSDT, ASTDT, AENDT, TRTEMFL)
#'
#' @caption Considering treatment end date (`trt_end_date` and `end_window`)
#' @info Derive `TRTEMFL` taking a treatment end window into account
#'
#' - In addition to the treatment-emergent checks explained in the above
#'   example, we now supply a treatment end date, `trt_end_date = TRTEDT` and
#'   an end window, `end_window = 10`. With these, any AE which started on or
#'   before treatment end date + 10 days is considered as treatment-emergent.
#'   Otherwise, those starting after the treatment end window are no longer
#'   flagged as treatment-emergent (e.g. record 7).
#' @code
#' derive_var_trtemfl(
#'   adae,
#'   start_date = ASTDT,
#'   end_date = AENDT,
#'   trt_start_date = TRTSDT,
#'   trt_end_date = TRTEDT,
#'   end_window = 10
#' ) %>% select(USUBJID, TRTSDT, TRTEDT, ASTDT, AENDT, TRTEMFL)
#'
#' @caption Considering treatment worsening (`initial_intensity` and `intensity`)
#' @info Derive a new variable named `TRTEM2FL` taking worsening after treatment
#'   start into account
#'
#' - We also now start look at changes in intensity following treatment start using
#'   the `initial_intensity` and `intensity` arguments. This only impacts AEs
#'   starting before treatment, and ending on or after treatment (or with missing
#'   AE end date). We can additionally consider treatment-emergence for an
#'   AE that was ongoing at the start of treatment which may have worsened
#'   as a result of treatment,  i.e. the most extreme intensity is greater than
#'   the initial intensity (e.g. records 3 and 9).
#' @code
#' derive_var_trtemfl(
#'   adae,
#'   new_var = TRTEM2FL,
#'   start_date = ASTDT,
#'   end_date = AENDT,
#'   trt_start_date = TRTSDT,
#'   trt_end_date = TRTEDT,
#'   end_window = 10,
#'   initial_intensity = AEITOXGR,
#'   intensity = AETOXGR
#' ) %>% select(USUBJID, TRTSDT, ASTDT, AENDT, AEITOXGR, AETOXGR, TRTEM2FL)
#'
#' @caption Worsening when the same AE is collected over multiple records
#'   (`intensity` and `group_var`)
#' @info Derive `TRTEMFL` taking worsening after treatment into account within a
#'   grouping variable
#'
#' - Firstly, to understand which records correspond to the same AE, we need
#'   to supply a grouping variable (`group_var`). Then this example works in a
#'   similar way to the above one, but here we don't have an initial intensity
#'   so we have to use the intensity of the AE at the time of treatment start.
#'   If an ongoing AE increases intensity after treatment start (i.e. worsens),
#'   then from that point on the records are considered treatment-emergent,
#'   unless after the treatment end window (e.g. records 4, 6 and 7).
#' @code
#' adae2 <- tribble(
#'   ~USUBJID, ~ASTDT,            ~AENDT,            ~AETOXGR, ~AEGRPID,
#'   # ongoing AE where intensity drops after treatment start
#'   "1",      ymd("2021-12-31"), ymd("2022-01-01"), "3",      "1",
#'   "1",      ymd("2022-01-02"), ymd("2022-01-11"), "2",      "1",
#'   # ongoing AE where intensity increases after treatment start
#'   "1",      ymd("2021-12-31"), ymd("2022-01-01"), "1",      "2",
#'   "1",      ymd("2022-01-02"), ymd("2022-01-11"), "2",      "2",
#'   # ongoing AE where intensity increases after treatment start and then drops
#'   "1",      ymd("2021-12-31"), ymd("2022-01-01"), "1",      "3",
#'   "1",      ymd("2022-01-02"), ymd("2022-01-11"), "2",      "3",
#'   "1",      ymd("2022-01-12"), ymd("2022-01-15"), "1",      "3"
#' ) %>%
#'   mutate(
#'     STUDYID = "AB42",
#'     TRTSDT = if_else(USUBJID == "1", ymd("2022-01-01"), NA),
#'     TRTEDT = if_else(USUBJID == "1", ymd("2022-04-30"), NA)
#'   )
#'
#' derive_var_trtemfl(
#'   adae2,
#'   start_date = ASTDT,
#'   end_date = AENDT,
#'   trt_start_date = TRTSDT,
#'   trt_end_date = TRTEDT,
#'   end_window = 10,
#'   intensity = AETOXGR,
#'   group_var = AEGRPID
#' ) %>% select(USUBJID, TRTSDT, ASTDT, AENDT, AETOXGR, AEGRPID, TRTEMFL)
#'
#' @caption Further Examples from PHUSE White Paper
#' @info Here we present more cases (some new, some similar to the examples above)
#'   which are aligned one-to-one with the scenarios in the
# nolint start
#'   [PHUSE White Paper](https://phuse.s3.eu-central-1.amazonaws.com/Deliverables/Safety+Analytics/WP-087+Recommended+Definition+of++Treatment-Emergent+Adverse+Events+in+Clinical+Trials+.pdf)
# nolint end
#' @code
#' adae3 <- tribble(
#'   ~USUBJID, ~TRTSDTM, ~TRTEDTM, ~ASTDTM, ~AENDTM, ~AEITOXGR, ~AETOXGR,
#'   # Patient 1: Pre-treatment AE
#'   "1", "2021-01-01", "2021-12-31", "2020-12-20", "2020-12-21", "2", "2",
#'   # Patient 2: On-treatment AE
#'   "2", "2021-01-01", "2021-12-31", "2021-12-20", "2021-12-21", "2", "2",
#'   # Patient 3: Pre-treatment AE, then on-treatment AE at same intensity
#'   "3", "2021-01-01", "2021-12-31", "2020-12-20", "2020-12-21", "2", "2",
#'   "3", "2021-01-01", "2021-12-31", "2021-12-20", "2021-12-21", "2", "2",
#'   # Patient 4: Pre-treatment AE, then on-treatment AE at wors. intensity
#'   "4", "2021-01-01", "2021-12-31", "2020-12-20", "2020-12-21", "2", "2",
#'   "4", "2021-01-01", "2021-12-31", "2021-12-20", "2021-12-21", "2", "3",
#'   # Patient 5: Pre-treatment AE, then on-treatment AE at impr. intensity
#'   "5", "2021-01-01", "2021-12-31", "2020-12-20", "2020-12-21", "2", "2",
#'   "5", "2021-01-01", "2021-12-31", "2021-12-20", "2021-12-21", "2", "1",
#'   # Patient 6: AE starting pre-treatment, continuing on-treatment, then 2nd AE at same intensity
#'   "6", "2021-01-01", "2021-12-31", "2020-12-23", "2021-01-21", "2", "2",
#'   "6", "2021-01-01", "2021-12-31", "2021-12-20", "2021-12-21", "2", "2",
#'   # Patient 7: AE starting pre-treatment, continuing on-treatment, then 2nd AE at wors. intensity
#'   "7", "2021-01-01", "2021-12-31", "2020-12-23", "2021-01-21", "2", "2",
#'   "7", "2021-01-01", "2021-12-31", "2021-12-20", "2021-12-21", "2", "3",
#'   # Patient 8: AE starting pre-treatment, continuing on-treatment, then 2nd AE at impr. intensity
#'   "8", "2021-01-01", "2021-12-31", "2020-12-23", "2021-01-21", "2", "2",
#'   "8", "2021-01-01", "2021-12-31", "2021-12-20", "2021-12-21", "2", "1",
#'   # Patient 9: AE starting pre-treatment, continuing on-treatment, and no change in intensity
#'   "9", "2021-01-01", "2021-12-31", "2020-12-23", "2021-01-21", "2", "2",
#'   # Patient 10: AE starting pre-treatment, continuing on-treatment, and wors. intensity
#'   "10", "2021-01-01", "2021-12-31", "2020-12-23", "2021-01-21", "2", "4",
#'   # Patient 11: AE starting pre-treatment, continuing on-treatment, and impr. intensity
#'   "11", "2021-01-01", "2021-12-31", "2020-12-23", "2021-01-21", "2", "1",
#'   # Patient 12: AE starting pre-treatment, worsening, then improving
#'   "12", "2021-01-01", "2021-12-31", "2020-12-23", "2021-01-21", "3", "2",
#'   # Patient 13: AE starting pre-treatment, improving, then worsening
#'   "13", "2021-01-01", "2021-12-31", "2020-12-23", "2021-01-21", "1", "2",
#' ) %>%
#'   mutate(
#'     ASTDTM = ymd(ASTDTM),
#'     AENDTM = ymd(AENDTM),
#'     TRTSDTM = ymd(TRTSDTM),
#'     TRTEDTM = ymd(TRTEDTM),
#'   )
#'
#' derive_var_trtemfl(
#'   adae3,
#'   new_var = TRTEMFL,
#'   trt_end_date = TRTEDTM,
#'   end_window = 0,
#'   initial_intensity = AEITOXGR,
#'   intensity = AETOXGR,
#'   subject_keys = exprs(USUBJID)
#' ) %>%
#'   select(USUBJID, TRTSDTM, TRTEDTM, ASTDTM, AENDTM, AEITOXGR, AETOXGR, TRTEMFL)
#'
derive_var_trtemfl <- function(dataset,
                               new_var = TRTEMFL,
                               start_date = ASTDTM,
                               end_date = AENDTM,
                               trt_start_date = TRTSDTM,
                               trt_end_date = NULL,
                               end_window = NULL,
                               ignore_time_for_trt_end = TRUE,
                               initial_intensity = NULL,
                               intensity = NULL,
                               group_var = NULL,
                               subject_keys = get_admiral_option("subject_keys")) {
  # Convert inputs to symbols
  new_var <- assert_symbol(enexpr(new_var))
  start_date <- assert_symbol(enexpr(start_date))
  end_date <- assert_symbol(enexpr(end_date))
  trt_start_date <- assert_symbol(enexpr(trt_start_date))
  trt_end_date <-
    assert_symbol(enexpr(trt_end_date), optional = TRUE)
  assert_integer_scalar(end_window, subset = "non-negative", optional = TRUE)
  assert_logical_scalar(ignore_time_for_trt_end)
  initial_intensity <-
    assert_symbol(enexpr(initial_intensity), optional = TRUE)
  intensity <- assert_symbol(enexpr(intensity), optional = TRUE)
  group_var <- assert_symbol(enexpr(group_var), optional = TRUE)

  # group_var is not specified
  # Check if both initial_intensity and intensity are provided
  if (is.null(group_var)) {
    if (is.null(initial_intensity) && !is.null(intensity)) {
      cli_abort(c(
        "{.arg intensity} argument was specified but not {.arg initial_intensity}",
        "Either both or none of them must be specified."
      ))
    }
    if (!is.null(initial_intensity) && is.null(intensity)) {
      cli_abort(c(
        "{.arg initial_intensity} argument was specified but not {.arg intensity}",
        "Either both or none of them must be specified."
      ))
    }
    # group_var is specified
  } else {
    if (!is.null(initial_intensity)) {
      cli_warn(c(
        "{.arg initial_intensity} argument is ignored when {.arg group_var} is specified",
        "Please only specify one of them."
      ))
    }
    if (is.null(subject_keys)) {
      cli_abort(c(
        "{.arg group_var} argument was specified but not {.arg subject_keys}",
        "{.arg subject_keys} argument must be provided when {.arg group_var} is specified."
      ))
    }
    assert_vars(subject_keys)
  }

  # Assert required variables
  required_vars <-
    expr_c(
      start_date,
      end_date,
      trt_start_date,
      trt_end_date,
      initial_intensity,
      intensity
    )
  if (!is.null(group_var)) {
    required_vars <- c(required_vars, group_var)
  }
  assert_data_frame(dataset, required_vars = required_vars)

  # Assert date variables
  assert_date_var(dataset, var = !!start_date)
  assert_date_var(dataset, var = !!end_date)
  assert_date_var(dataset, var = !!trt_start_date)
  if (!is.null(trt_end_date)) {
    assert_date_var(dataset, var = !!trt_end_date)
  }

  # end window condition
  if (is.null(end_window)) {
    end_cond <- expr(TRUE)
  } else {
    if (is.null(trt_end_date)) {
      cli_abort(c(
        "{.arg end_window} argument was specified but not {.arg trt_end_date}",
        "Either both or none of them must be specified."
      ))
    }
    if (ignore_time_for_trt_end) {
      end_cond <- expr(
        (is.na(!!trt_end_date) |
          date(!!start_date) <= date(!!trt_end_date) + days(end_window))
      )
    } else {
      end_cond <-
        expr(
          (is.na(!!trt_end_date) |
            !!start_date <= !!trt_end_date + days(end_window))
        )
    }
  }


  # new_ae_cond: Y - new AE, N - AE exists before trt_start_date
  new_ae_cond <- get_new_tmp_var(dataset)
  if (is.null(group_var)) {
    dataset <- dataset %>%
      mutate(
        !!new_ae_cond := if_else(!!start_date >= !!trt_start_date, "Y", "N")
      )
  } else {
    dataset <- dataset %>%
      derive_vars_merged(
        dataset_add = dataset,
        by_vars = expr_c(subject_keys, group_var),
        order = exprs(!!start_date),
        new_vars = exprs(!!new_ae_cond := "N"),
        filter_add = !!start_date < !!trt_start_date,
        mode = "last"
      ) %>%
      mutate(
        !!new_ae_cond := if_else(is.na(!!new_ae_cond), "Y", "N")
      )
  }

  if (is.null(intensity)) {
    worsening_cond <- expr(FALSE)
  } else {
    if (is.null(group_var)) {
      worsening_cond <-
        expr(
          !!start_date < !!trt_start_date &
            (!!initial_intensity < !!intensity | is.na(!!initial_intensity) |
              is.na(!!intensity))
        )
    } else {
      prev_intensity <- get_new_tmp_var(dataset)
      worsen_date <- get_new_tmp_var(dataset)

      dataset <- dataset %>%
        arrange(USUBJID, !!group_var, !!start_date) %>%
        group_by(USUBJID, !!group_var) %>%
        mutate(
          !!prev_intensity := lag(!!intensity),
          !!worsen_date :=
            case_when(
              !is.na(!!start_date) & !is.na(!!trt_start_date) & !is.na(!!prev_intensity) &
                !!start_date >= !!trt_start_date &
                (!!intensity > !!prev_intensity) ~ !!start_date,
              TRUE ~ NA
            )
        ) %>%
        fill(!!worsen_date, .direction = "down") %>%
        ungroup()

      worsening_cond <- expr(!is.na(!!worsen_date))
    }
  }


  # Derive TRTEMFL based on conditions

  dataset <- dataset %>%
    mutate(
      !!new_var := case_when(
        is.na(!!trt_start_date) ~ NA_character_,
        !!end_date < !!trt_start_date ~ NA_character_,
        is.na(!!start_date) ~ "Y",
        !!new_ae_cond == "Y" & !!end_cond ~ "Y", # new AE
        !!worsening_cond ~ "Y" # worsened AE
      )
    ) %>%
    # Remove temporary variable
    remove_tmp_vars()
}
