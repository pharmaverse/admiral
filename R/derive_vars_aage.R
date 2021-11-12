#' Derive Analysis Age
#'
#' Derives analysis age (`AAGE`) and analysis age unit (`AAGEU`)
#'
#' @param dataset Input dataset
#'
#'   The columns specified by the `start_date` and the `end_date` parameter are
#'   expected.
#'
#' @param start_date The start date
#'
#'   A date or date-time object is expected.
#'
#'   Default: `BRTHDT`
#'
#' @param end_date The end date
#'
#'   A date or date-time object is expected.
#'
#'   Default: `RANDDT`
#'
#' @param unit Unit
#'
#'   The age is derived in the specified unit
#'
#'   Default: 'years'
#'
#'   Permitted Values: 'years', 'months', 'days', 'hours', 'minutes', 'seconds'
#'
#' @details The age is derived as the integer part of the duration from start to
#'   end date in the specified unit.
#'
#' @author Stefan Bundfuss
#'
#' @return The input dataset with ``AAGE`` and ``AAGEU`` added
#'
#' @export
#'
#' @seealso [derive_vars_duration()]
#'
#' @examples
#' data <- tibble::tribble(
#'   ~BRTHDT, ~RANDDT,
#'   lubridate::ymd("1984-09-06"), lubridate::ymd("2020-02-24")
#' )
#'
#' derive_vars_aage(data)
derive_vars_aage <- function(dataset,
                             start_date = BRTHDT,
                             end_date = RANDDT,
                             unit = "years") {

  start_date <- assert_symbol(enquo(start_date))
  end_date <- assert_symbol(enquo(end_date))
  assert_data_frame(dataset, required_vars = quo_c(start_date, end_date))
  assert_character_scalar(
    unit,
    values = c("years", "months", "days", "hours", "minutes", "seconds")
  )

  derive_vars_duration(
    dataset,
    new_var = AAGE,
    new_var_unit = AAGEU,
    start_date = !!start_date,
    end_date = !!end_date,
    out_unit = unit,
    add_one = FALSE,
    trunc_out = TRUE
  )
}

#' Derive Analysis Age
#'
#' `derive_aage()` was renamed to `derive_vars_aage()` to create a
#' more consistent API.
#'
#' @keywords internal
#'
#' @export
derive_aage <- function(dataset,
                        start_date = BRTHDT,
                        end_date = RANDDT,
                        unit = "years") {
  deprecate_warn("0.3.0", "derive_aage()", "derive_vars_aage()")
  derive_vars_aage(
    dataset,
    start_date = !!enquo(start_date),
    end_date = !!enquo(end_date),
    unit = unit
  )
}

#' Derive Age Groups
#'
#' Functions for deriving standardized age groups.
#'
#' @param dataset Input dataset.
#' @param age_var AGE variable.
#' @param new_var New variable to be created.
#'
#' @return `dataset` with new column `new_var` of class factor.
#'
#' @author Ondrej Slama
#'
#' @name derive_agegr_fda
NULL

#' @rdname derive_agegr_fda
#' @export
#' @details `derive_agegr_fda` Derive age groups according to FDA
#' (\url{https://prsinfo.clinicaltrials.gov/results_definitions.html} ->
#' Baseline Measure Information).
#' @examples
#' library(dplyr, warn.conflicts = FALSE)
#' library(admiral.test)
#' data(dm)
#'
#' dm %>%
#'   derive_agegr_fda(AGE, AGEGR1) %>%
#'   select(SUBJID, AGE, AGEGR1)
#'
#' data.frame(AGE = 1:100) %>%
#'   derive_agegr_fda(age_var = AGE, new_var = AGEGR1)
#'
derive_agegr_fda <- function(dataset, age_var, new_var) {

  age_var <- assert_symbol(enquo(age_var))
  new_var <- assert_symbol(enquo(new_var))
  assert_data_frame(dataset, required_vars = quo_c(age_var))

  out <- mutate(
    dataset,
    !!new_var := cut(
      x = !!age_var,
      breaks = c(0, 19, 65, Inf),
      labels = c("<=18", "19-64", ">=65"),
      include.lowest = TRUE,
      right = FALSE
    )
  )
  if (anyNA(dplyr::pull(out, !!new_var))) {
    out <- mutate(out, !!new_var := addNA(!!new_var))
  }
  out
}

#' @rdname derive_agegr_fda
#' @export
#' @param adults Logical, should the age group be valid for adults studies or pediatric studies?
#' Defaults to `TRUE` (adult studies).
#' @details `derive_agegr_ema` Derive age groups according to EMA
#' (\url{https://eudract.ema.europa.eu/result.html} -> Results - Data Dictionary -> Age range).
#' @examples
#' library(dplyr, warn.conflicts = FALSE)
#' library(admiral.test)
#' data(dm)
#'
#' dm %>%
#'   derive_agegr_ema(AGE, AGEGR1) %>%
#'   select(SUBJID, AGE, AGEGR1)
#'
#' data.frame(AGE = 1:100) %>%
#'   derive_agegr_ema(age_var = AGE, new_var = AGEGR1)
#'
#' data.frame(AGE = 1:20) %>%
#'   derive_agegr_ema(age_var = AGE, new_var = AGEGR1, adults = FALSE)
derive_agegr_ema <- function(dataset, age_var, new_var, adults = TRUE) {

  age_var <- assert_symbol(enquo(age_var))
  new_var <- assert_symbol(enquo(new_var))
  assert_data_frame(dataset, required_vars = quo_c(age_var))
  assert_logical_scalar(adults)

  if (adults) {
    out <- mutate(
      dataset,
      !!new_var := cut(
        x = !!age_var,
        breaks = c(18, 65, 85, Inf),
        labels = c("18-64", "65-84", ">=85"),
        include.lowest = TRUE,
        right = FALSE
      )
    )
  } else {
    out <- mutate(
      dataset,
      !!new_var := cut(
        x = !!age_var,
        breaks = c(-Inf, 2, 12, 18),
        labels = c("0-1 (Newborns / Infants / Toddlers)",
                   "2-11 (Children)", "12-17 (Adolescents)"),
        include.lowest = FALSE,
        right = FALSE
      )
    )
  }
  if (anyNA(dplyr::pull(out, !!new_var))) {
    out <- mutate(out, !!new_var := addNA(!!new_var))
  }
  out
}
