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
#'   Permitted Values: 'years', 'months', 'weeks', 'days', 'hours', 'minutes', 'seconds'
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
    values = c("years", "months", "weeks", "days", "hours", "minutes", "seconds")
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

#' Derive Age in Years
#'
#' @details This function is used to convert age variables into years.
#' These can then be used to create age groups.
#'
#' @param dataset Input dataset.
#' @param age_var AGE variable.
#' @param age_unit AGE unit variable.
#'
#'   The AGE unit variable is used to convert AGE to 'years' so that grouping can occur.
#'   This is only used when the age_var variable does not have a corresponding unit in the dataset.
#'
#'   Default: NULL
#'
#'   Permitted Values: 'years', 'months', 'weeks', 'days', 'hours', 'minutes', 'seconds'
#'
#' @param new_var New AGE variable to be created in years.
#'
#' @author Michael Thorpe
#'
#' @return The input dataset with new_var paramater added in years.
#'
#' @export
#'
#' @examples
#'
#' library(dplyr, warn.conflicts = FALSE)
#'
#' data <- data.frame(AGE = c(27, 24, 3, 4, 1),
#'                    AGEU = c("days", "months", "years", "weeks", "years"))
#'
#' data %>%
#'      derive_var_age_years(., AGE, new_var = AAGE)
#'
#' data.frame(AGE = c(12, 24, 36, 48)) %>%
#'  derive_var_age_years(., AGE, age_unit = "months", new_var = AAGE)
#'
derive_var_age_years <- function(dataset, age_var, age_unit = NULL, new_var) {

  age_variable <- assert_symbol(enquo(age_var))
  assert_data_frame(dataset, required_vars = quo_c(age_variable))

  age_var <- pull(dataset, !!age_variable)
  assert_numeric_vector(age_var)

  age_var <- age_variable
  unit_var <- paste0(quo_get_expr(age_var), "U")

  new_var <- assert_symbol(enquo(new_var))
  warn_if_vars_exist(dataset, quo_text(new_var))

  if (!unit_var %in% colnames(dataset)) {

    if (is.null(age_unit)) {

      err_msg <- paste(
        "There is no variable unit:", unit_var, "associated with", quo_get_expr(age_var),
        "and the argument `age_unit` is missing. Please specify a value for `age_unit`"
      )
      abort(err_msg)
    } else{
      assert_character_scalar(tolower(age_unit), values = c("years", "months", "weeks", "days",
                                                            "hours", "minutes", "seconds"))
      ds <- dataset %>%
        mutate(!!new_var := time_length(duration(!!age_var, units = tolower(age_unit)),
                                        unit = "years"))
    }
  } else {

    unit <- tolower(unique(pull(dataset, !!sym(unit_var))))
    assert_character_vector(unit,
                            values = c(NA, "years", "months", "weeks", "days",
                                       "hours", "minutes", "seconds"))

    if (!is.null(age_unit)) {

      if (length(unit) > 1) {

        msg <- paste(
          "The variable unit", unit_var, "is associated with", quo_get_expr(age_var),
          "and contatins multiple values but the argument `age_unit`
          has been specified with a single different value.",
          "The `age_unit` argument is ignored and the grouping will based on",
          unit_var)
        warn(msg)

      }

      else if (unit != tolower(age_unit)) {

        msg <- paste(
          "The variable unit", unit_var, "is associated with", quo_get_expr(age_var),
          "but the argument `age_unit` has been specified with a different value.",
          "The `age_unit` argument is ignored and the grouping will based on",
          unit_var)
        warn(msg)
      }

    }

    average_durations <- c(seconds = 365.25 * 24 * 60 * 60,
                           minutes = 365.25 * 24 * 60,
                           hours = 365.25 * 24,
                           days = 365.25,
                           weeks = 365.25 / 7,
                           months = 12,
                           years = 1)

    ds <- dataset %>%
      mutate(!!new_var := !!age_var / unname(average_durations[tolower(!!sym(unit_var))])
      )
  }

}

#' Derive Age Groups
#'
#' Functions for deriving standardized age groups.
#'
#' @param dataset Input dataset.
#' @param age_var AGE variable.
#' @param age_unit AGE unit variable.
#'
#'   The AGE unit variable is used to convert AGE to 'years' so that grouping can occur.
#'   This is only used when the age_var variable does not have a corresponding unit in the dataset.
#'
#'   Default: NULL
#'
#'   Permitted Values: 'years', 'months', 'weeks', 'days', 'hours', 'minutes', 'seconds'
#'
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
#' @details `derive_agegr_fda` Derive age groups according to FDA. `age_var` will
#'  be split in categories: <18, 18-64, >=65.
#' @examples
#' library(dplyr, warn.conflicts = FALSE)
#' library(admiral.test)
#' data(dm)
#'
#' dm %>%
#'   derive_agegr_fda(age_var = AGE, new_var = AGEGR1) %>%
#'   select(SUBJID, AGE, AGEGR1)
#'
#' data <- tibble::tribble(
#'   ~BRTHDT, ~RANDDT,
#'   lubridate::ymd("1984-09-06"), lubridate::ymd("2020-02-24")
#'   )
#'
#' data %>%
#'   derive_vars_aage(unit = "months") %>%
#'   derive_agegr_fda(AAGE, age_unit = NULL, AGEGR1)
#'
#' data.frame(AGE = 1:100) %>%
#'   derive_agegr_fda(age_var = AGE, age_unit = "years", new_var = AGEGR1)
#'
derive_agegr_fda <- function(dataset, age_var, age_unit = NULL, new_var) {
  deprecate_warn("0.6.0", "derive_agegr_fda()", "derive_var_agegr_fda()")
  derive_var_agegr_fda(dataset = dataset,
                       age_var = !!enquo(age_var),
                       age_unit = age_unit,
                       new_var = !!enquo(new_var))
}

#' @rdname derive_agegr_fda
#' @export
#' @details `derive_agegr_ema` Derive age groups according to EMA
#' (\url{https://eudract.ema.europa.eu/result.html} -> Results - Data Dictionary -> Age range).
#' `age_var` will be split into categories: 0-27 days (Newborns), 28 days to
#' 23 months (Infants and Toddlers), 2-11 (Children), 12-17 (Adolescents), 18-64,
#'  65-84, >=85.
#' @examples
#' library(dplyr, warn.conflicts = FALSE)
#' library(admiral.test)
#' data(dm)
#'
#' dm %>%
#'   derive_agegr_ema(age_var = AGE, new_var = AGEGR1) %>%
#'   select(SUBJID, AGE, AGEGR1)
#'
#' data.frame(AGE = 1:100) %>%
#'   derive_agegr_ema(age_var = AGE, age_unit = "years", new_var = AGEGR1)
#'
#' data.frame(AGE = 1:20) %>%
#'   derive_agegr_ema(age_var = AGE, age_unit = "years", new_var = AGEGR1)
derive_agegr_ema <- function(dataset, age_var, age_unit = NULL, new_var) {
  deprecate_warn("0.6.0", "derive_agegr_ema()", "derive_var_agegr_ema()")
  derive_var_agegr_ema(dataset = dataset,
                       age_var = !!enquo(age_var),
                       age_unit = age_unit,
                       new_var = !!enquo(new_var))
}

#' Derive Age Groups
#'
#' Functions for deriving standardized age groups.
#'
#' @param dataset Input dataset.
#' @param age_var AGE variable.
#' @param age_unit AGE unit variable.
#'
#'   The AGE unit variable is used to convert AGE to 'years' so that grouping can occur.
#'   This is only used when the age_var variable does not have a corresponding unit in the dataset.
#'
#'   Default: NULL
#'
#'   Permitted Values: 'years', 'months', 'weeks', 'days', 'hours', 'minutes', 'seconds'
#'
#' @param new_var New variable to be created.
#'
#' @return `dataset` with new column `new_var` of class factor.
#'
#' @author Ondrej Slama
#'
#' @name derive_var_agegr_fda
NULL

#' @rdname derive_var_agegr_fda
#' @export
#' @details `derive_var_agegr_fda` Derive age groups according to FDA. `age_var`
#' will be split in categories: <18, 18-64, >=65.
#' @examples
#' library(dplyr, warn.conflicts = FALSE)
#' library(admiral.test)
#' data(dm)
#'
#' dm %>%
#'   derive_var_agegr_fda(age_var = AGE, new_var = AGEGR1) %>%
#'   select(SUBJID, AGE, AGEGR1)
#'
#' data <- tibble::tribble(
#'   ~BRTHDT, ~RANDDT,
#'   lubridate::ymd("1984-09-06"), lubridate::ymd("2020-02-24")
#'   )
#'
#' data %>%
#'   derive_vars_aage(unit = "months") %>%
#'   derive_var_agegr_fda(AAGE, age_unit = NULL, AGEGR1)
#'
#' data.frame(AGE = 1:100) %>%
#'   derive_var_agegr_fda(age_var = AGE, age_unit = "years", new_var = AGEGR1)
#'
derive_var_agegr_fda <- function(dataset, age_var, age_unit = NULL, new_var) {

  age_var <- assert_symbol(enquo(age_var))
  new_var <- assert_symbol(enquo(new_var))
  warn_if_vars_exist(dataset, quo_text(new_var))

  ds <- derive_var_age_years(dataset, !!age_var, age_unit, new_var = temp_age)

  out <- ds %>%
    mutate(
      !!new_var := cut(
        x = temp_age,
        breaks = c(0, 18, 65, Inf),
        labels = c("<18", "18-64", ">=65"),
        include.lowest = TRUE,
        right = FALSE
      )
    ) %>%
    select(-temp_age)

  if (anyNA(dplyr::pull(out, !!new_var))) {
    out <- mutate(out, !!new_var := addNA(!!new_var))
  }
  out
}

#' @rdname derive_var_agegr_fda
#' @export
#' @details `derive_var_agegr_ema` Derive age groups according to EMA
#' (\url{https://eudract.ema.europa.eu/result.html} -> Results - Data Dictionary -> Age range).
#' `age_var` will be split into categories: 0-27 days (Newborns), 28 days to
#' 23 months (Infants and Toddlers), 2-11 (Children), 12-17 (Adolescents), 18-64,
#'  65-84, >=85.
#' @examples
#' library(dplyr, warn.conflicts = FALSE)
#' library(admiral.test)
#' data(dm)
#'
#' dm %>%
#'   derive_var_agegr_ema(age_var = AGE, new_var = AGEGR1) %>%
#'   select(SUBJID, AGE, AGEGR1)
#'
#' data.frame(AGE = 1:100) %>%
#'   derive_var_agegr_ema(age_var = AGE, age_unit = "years", new_var = AGEGR1)
#'
#' data.frame(AGE = 1:20) %>%
#'   derive_var_agegr_ema(age_var = AGE, age_unit = "years", new_var = AGEGR1)
derive_var_agegr_ema <- function(dataset, age_var, age_unit = NULL, new_var) {

  age_var <- assert_symbol(enquo(age_var))
  new_var <- assert_symbol(enquo(new_var))
  warn_if_vars_exist(dataset, quo_text(new_var))

  ds <-   derive_var_age_years(dataset, !!age_var, age_unit, new_var = temp_age)

  out <- mutate(
    ds,
    !!new_var := cut(
      x = temp_age,
      breaks = c(-Inf, (28 / 365.25), 2, 12, 18, 65, 85, Inf),
      labels = c("0-27 days (Newborns)", "28 days to 23 months (Infants and Toddlers)",
                 "2-11 (Children)", "12-17 (Adolescents)", "18-64", "65-84", ">=85"),
      include.lowest = FALSE,
      right = FALSE
    )
  ) %>%
    select(-temp_age)

  if (anyNA(dplyr::pull(out, !!new_var))) {
    out <- mutate(out, !!new_var := addNA(!!new_var))
  }
  out
}
