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
#' @seealso [derive_duration()]
#'
#' @examples
#' data <- tibble::tribble(
#'   ~BRTHDT, ~RANDDT,
#'   lubridate::ymd("1984-09-06"), lubridate::ymd("2020-02-24")
#' )
#'
#' derive_aage(data)
derive_aage <- function(dataset,
                        start_date = BRTHDT,
                        end_date = RANDDT,
                        unit = "years") {
  derive_duration(
    dataset,
    new_var = AAGE,
    new_var_unit = AAGEU,
    start_date = !!enquo(start_date),
    end_date = !!enquo(end_date),
    out_unit = unit,
    add_one = FALSE,
    trunc_out = TRUE
  )
}


#' Derive age groups
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
#' (\url{https://prsinfo.clinicaltrials.gov/results_definitions.html} -> Baseline Measure Information).
#' @examples
#' derive_agegr_fda(data.frame(age = 1:100), age_var = age, new_var = agegr1)
derive_agegr_fda <- function(dataset, age_var, new_var) {
  assert_that(is.data.frame(dataset))

  age_var <- enquo(age_var)
  new_var <- enquo(new_var)

  assert_has_variables(dataset, as_string(rlang::quo_get_expr(age_var)))

  mutate(
    dataset,
    !!new_var := cut(
      x = !!age_var,
      breaks = c(-Inf, 19, 65, Inf),
      labels = c(NA_character_, "19-64", ">=65"),
      include.lowest = T,
      right = F
    )
  )
}

#' @rdname derive_agegr_fda
#' @export
#' @param adults Logical, should the age group be valid for adults studies or pediatric studies?
#' Defaults to `TRUE` (adult studies).
#' @details `derive_agegr_ema` Derive age groups according to EMA
#' (\url{https://eudract.ema.europa.eu/result.html} -> Results - Data Dictionary -> Age range).
#' @examples
#' derive_agegr_ema(data.frame(age = 1:100), age_var = age, new_var = agegr1)
#' derive_agegr_ema(data.frame(age = 1:20), age_var = age, new_var = agegr1, adults = F)
derive_agegr_ema <- function(dataset, age_var, new_var, adults = TRUE) {
  assert_that(
    is.data.frame(dataset),
    rlang::is_scalar_logical(adults)
  )

  age_var <- enquo(age_var)
  new_var <- enquo(new_var)

  assert_has_variables(dataset, as_string(rlang::quo_get_expr(age_var)))

  if (adults) {
    mutate(
      dataset,
      !!new_var := cut(
        x = !!age_var,
        breaks = c(-Inf, 18, 65, 85, Inf),
        labels = c(NA_character_, "18-64", "65-84", ">=85"),
        include.lowest = T,
        right = F
      )
    )
  } else {
    mutate(
      dataset,
      !!new_var := cut(
        x = !!age_var,
        breaks = c(-Inf, 2, 12, 18, Inf),
        labels = c("0-1 (Newborns / Infants / Toddlers)",
                   "2-11 (Children)", "12-17 (Adolescents)", ">=18"),
        include.lowest = T,
        right = F
      )
    )
  }
}
