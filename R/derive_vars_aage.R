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
#'   Refer to `derive_vars_dt()` to impute and derive a date from a date character
#'   vector to a date object.
#'
#'   Default: `BRTHDT`
#'
#' @param end_date The end date
#'
#'   A date or date-time object is expected.
#'
#'   Refer to `derive_vars_dt()` to impute and derive a date from a date character
#'   vector to a date object.
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
#'   end date in the specified unit. When 'years' or 'months' are specified in the `out_unit`
#'   parameter, because of the underlying `lubridate::time_length()` function that is used
#'   here, results are calculated based on the actual calendar length of months or years
#'   rather than assuming equal days every month (30.4375 days) or every year (365.25 days).
#'
#' @author Stefan Bundfuss
#'
#' @return The input dataset with ``AAGE`` and ``AAGEU`` added
#'
#' @family der_adsl
#' @keywords der_adsl
#'
#' @export
#'
#' @seealso [derive_vars_duration()]
#'
#' @examples
#' library(tibble)
#' library(lubridate)
#'
#' data <- tribble(
#'   ~BRTHDT, ~RANDDT,
#'   ymd("1984-09-06"), ymd("2020-02-24")
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
#' @family der_adsl
#' @keywords der_adsl
#'
#' @author Michael Thorpe
#'
#' @return The input dataset with new_var parameter added in years.
#'
#' @export
#'
#' @examples
#'
#' data <- data.frame(
#'   AGE = c(27, 24, 3, 4, 1),
#'   AGEU = c("days", "months", "years", "weeks", "years")
#' )
#'
#' data %>%
#'   derive_var_age_years(., AGE, new_var = AAGE)
#'
#' data.frame(AGE = c(12, 24, 36, 48)) %>%
#'   derive_var_age_years(., AGE, age_unit = "months", new_var = AAGE)
derive_var_age_years <- function(dataset, age_var, age_unit = NULL, new_var) {
  age_variable <- assert_symbol(enquo(age_var))
  assert_data_frame(dataset, required_vars = quo_c(age_variable))

  age_var <- pull(dataset, !!age_variable)
  assert_numeric_vector(age_var)

  age_var <- age_variable
  unit_var <- paste0(quo_get_expr(age_var), "U")

  age_unit <- assert_character_scalar(
    age_unit,
    values = c(
      "years", "months", "weeks", "days",
      "hours", "minutes", "seconds"
    ),
    case_sensitive = FALSE,
    optional = TRUE
  )

  new_var <- assert_symbol(enquo(new_var))
  warn_if_vars_exist(dataset, quo_text(new_var))

  if (!unit_var %in% colnames(dataset)) {
    if (is.null(age_unit)) {
      err_msg <- paste(
        "There is no variable unit:", unit_var, "associated with", quo_get_expr(age_var),
        "and the argument `age_unit` is missing. Please specify a value for `age_unit`"
      )
      abort(err_msg)
    } else {
      ds <- dataset %>%
        mutate(
          !!new_var := time_length(duration(!!age_var, units = age_unit),
            unit = "years"
          )
        )
    }
  } else {
    unit <- unique(tolower(pull(dataset, !!sym(unit_var))))
    assert_character_vector(
      unit,
      values = c(
        NA, "years", "months", "weeks", "days",
        "hours", "minutes", "seconds"
      )
    )

    if (!is.null(age_unit)) {
      if (length(unit) > 1) {
        msg <- paste(
          "The variable unit", unit_var, "is associated with", quo_get_expr(age_var),
          "and contatins multiple values but the argument `age_unit`
          has been specified with a single different value.",
          "The `age_unit` argument is ignored and the grouping will based on",
          unit_var
        )
        warn(msg)
      } else if (unit != age_unit) {
        msg <- paste(
          "The variable unit", unit_var, "is associated with", quo_get_expr(age_var),
          "but the argument `age_unit` has been specified with a different value.",
          "The `age_unit` argument is ignored and the grouping will based on",
          unit_var
        )
        warn(msg)
      }
    }

    average_durations <- c(
      seconds = 365.25 * 24 * 60 * 60,
      minutes = 365.25 * 24 * 60,
      hours = 365.25 * 24,
      days = 365.25,
      weeks = 365.25 / 7,
      months = 12,
      years = 1
    )

    ds <- dataset %>%
      mutate(!!new_var := !!age_var / unname(average_durations[tolower(!!sym(unit_var))]))
  }
}


#' Derive Age Groups
#'
#' @description
#' `r lifecycle::badge("deprecated")`
#'
#' These functions are *deprecated*.
#'
#' @param dataset Input dataset
#'
#' @param age_var AGE variable
#'
#' @param age_unit AGE unit variable
#'
#' @param new_var New variable to create inside `dataset`
#'
#' @keywords deprecated
#'
#' @author Ondrej Slama
#'
#' @name derive_var_agegr_fda
NULL

#' @rdname derive_var_agegr_fda
#'
#' @keywords deprecated
#'
#' @export
derive_var_agegr_fda <- function(dataset, age_var, age_unit = NULL, new_var) {
  deprecate_warn("0.8.0", "derive_var_agegr_fda()", details = "Please create a user defined function instead.")

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
#'
#' @keywords deprecated
#'
#' @export
derive_var_agegr_ema <- function(dataset, age_var, age_unit = NULL, new_var) {
  deprecate_warn("0.8.0", "derive_var_agegr_ema()", details = "Please create a user defined function instead.")

  age_var <- assert_symbol(enquo(age_var))
  new_var <- assert_symbol(enquo(new_var))
  warn_if_vars_exist(dataset, quo_text(new_var))

  ds <- derive_var_age_years(dataset, !!age_var, age_unit, new_var = temp_age)

  out <- mutate(
    ds,
    !!new_var := cut(
      x = temp_age,
      breaks = c(-Inf, (28 / 365.25), 2, 12, 18, 65, 85, Inf),
      labels = c(
        "0-27 days (Newborns)", "28 days to 23 months (Infants and Toddlers)",
        "2-11 (Children)", "12-17 (Adolescents)", "18-64", "65-84", ">=85"
      ),
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
