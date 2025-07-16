#' Derive Relative Day Variables
#'
#' Adds relative day variables (`--DY`) to the dataset, e.g., `ASTDY` and
#' `AENDY`.
#'
#' @param dataset
#'   `r roxygen_param_dataset(expected_vars = c("reference_date", "source_vars"))`
#'
#' @param reference_date A date or date-time column, e.g., date of first treatment
#'   or date-time of last exposure to treatment.
#'
#'   Refer to `derive_vars_dt()` to impute and derive a date from a date
#'   character vector to a date object.
#'
#' @param source_vars A list of datetime or date variables created using
#'   `exprs()` from which dates are to be extracted. This can either be a list of
#'   date(time) variables or named `--DY` variables and corresponding --DT(M)
#'   variables e.g. `exprs(TRTSDTM, ASTDTM, AENDT)` or `exprs(TRTSDT, ASTDTM,
#'   AENDT, DEATHDY = DTHDT)`. If the source variable does not end in --DT(M), a
#'   name for the resulting `--DY` variable must be provided.
#'
#'
#' @details The relative day is derived as number of days from the reference
#'   date to the end date. If it is nonnegative, one is added. I.e., the
#'   relative day of the reference date is 1. Unless a name is explicitly
#'   specified, the name of the resulting relative day variable is generated
#'   from the source variable name by replacing DT (or DTM as appropriate) with
#'   DY.
#'
#' @return The input dataset with `--DY` corresponding to the `--DTM` or `--DT`
#'   source variable(s) added
#'
#' @family der_date_time
#' @keywords der_gen der_date_time
#'
#' @export
#'
#' @examples
#' library(tibble)
#' library(lubridate)
#' library(dplyr, warn.conflicts = FALSE)
#'
#' datain <- tribble(
#'   ~TRTSDTM, ~ASTDTM, ~AENDT,
#'   "2014-01-17T23:59:59", "2014-01-18T13:09:O9", "2014-01-20"
#' ) %>%
#'   mutate(
#'     TRTSDTM = as_datetime(TRTSDTM),
#'     ASTDTM = as_datetime(ASTDTM),
#'     AENDT = ymd(AENDT)
#'   )
#'
#' derive_vars_dy(
#'   datain,
#'   reference_date = TRTSDTM,
#'   source_vars = exprs(TRTSDTM, ASTDTM, AENDT)
#' )
#'
#' # specifying name of new variables
#' datain <- tribble(
#'   ~TRTSDT, ~DTHDT,
#'   "2014-01-17", "2014-02-01"
#' ) %>%
#'   mutate(
#'     TRTSDT = ymd(TRTSDT),
#'     DTHDT = ymd(DTHDT)
#'   )
#'
#' derive_vars_dy(
#'   datain,
#'   reference_date = TRTSDT,
#'   source_vars = exprs(TRTSDT, DEATHDY = DTHDT)
#' )
derive_vars_dy <- function(dataset,
                           reference_date,
                           source_vars) {
  # assertions
  reference_date <- assert_symbol(enexpr(reference_date))
  assert_vars(source_vars)
  assert_data_frame(dataset, required_vars = expr_c(source_vars, reference_date))

  # Warn if `--DY` variables already exist
  n_vars <- length(source_vars)
  source_names <- names(source_vars)

  bad_vars <- vars2chr(source_vars)[(
    (source_names == "" | source_names == vars2chr(source_vars)) &
      !str_detect(vars2chr(source_vars), "(DT|DTM)$"))]

  if (length(bad_vars > 0)) {
    cli_abort(c(
      "{.arg source_vars} must end in DT or DTM or be explicitly and uniquely named.",
      i = "Please name or rename the following source_vars: {.var {bad_vars}}"
    ))
  }

  # named vector passed to `.names` in `across()` to derive name of dy_vars
  dy_vars <- set_names(if_else(
    source_names == "",
    str_replace_all(vars2chr(source_vars), "(DT|DTM)$", "DY"),
    source_names
  ), vars2chr(source_vars))

  warn_if_vars_exist(dataset, dy_vars)

  dataset %>%
    mutate(
      across(
        .cols = vars2chr(unname(source_vars)),
        .fns = ~ compute_duration(start_date = !!reference_date, end_date = .x),
        .names = "{dy_vars}"
      )
    )
}
