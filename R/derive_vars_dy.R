#' Derive Relative Day Variables
#'
#' Adds relative day variables (`--DY`) to the dataset, e.g., `ASTDY` and
#' `AENDY`.
#'
#' @param dataset Input dataset
#'
#'   The columns specified by the `reference_date` and the `source_vars`
#'   parameter are expected.
#'
#' @param reference_date The start date column, e.g., date of first treatment
#'
#'   A date or date-time object column is expected.
#'
#'   Refer to `derive_vars_dt()` to impute and derive a date from a date
#'   character vector to a date object.
#'
#' @param source_vars A list of datetime or date variables created using
#'   `vars()` from which dates are to be extracted. This can either be a list of
#'   date(time) variables or named `--DY` variables and corresponding --DT(M)
#'   variables e.g. `vars(TRTSDTM, ASTDTM, AENDT)` or `vars(TRTSDT, ASTDTM,
#'   AENDT, DEATHDY = DTHDT)`. If the source variable does not end in --DT(M), a
#'   name for the resulting `--DY` variable must be provided.
#'
#' @author Teckla Akinyi
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
#'   source_vars = vars(TRTSDTM, ASTDTM, AENDT)
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
#'   source_vars = vars(TRTSDT, DEATHDY = DTHDT)
#' )
derive_vars_dy <- function(dataset,
                           reference_date,
                           source_vars) {
  # assertions
  reference_date <- assert_symbol(enquo(reference_date))
  assert_vars(source_vars)
  assert_data_frame(dataset, required_vars = quo_c(source_vars, reference_date))

  # Warn if `--DY` variables already exist
  n_vars <- length(source_vars)
  source_names <- names(source_vars)

  bad_vars <- vars2chr(source_vars)[(
    (source_names == "" | source_names == vars2chr(source_vars)) &
      !str_detect(vars2chr(source_vars), "(DT|DTM)$"))]

  if (length(bad_vars > 0)) {
    err_msg <-
      sprintf(
        paste0(
          "source_vars must end in DT or DTM or be explicitly and uniquely named.\n",
          "Please name or rename the following source_vars:\n", "%s"
        ),
        paste0(bad_vars, collapse = ", ")
      )
    abort(err_msg)
  }

  dy_vars <- if_else(
    source_names == "",
    stringr::str_replace_all(vars2chr(source_vars), "(DT|DTM)$", "DY"),
    source_names
  )
  warn_if_vars_exist(dataset, dy_vars)

  if (n_vars > 1L) {
    dataset %>%
      mutate_at(
        .vars = source_vars,
        .funs = list(temp = ~
          compute_duration(start_date = eval(reference_date), end_date = .))
      ) %>%
      rename_at(
        vars(ends_with("temp")),
        ~dy_vars
      )
  } else {
    dataset %>%
      mutate(
        !!sym(dy_vars) :=
          compute_duration(start_date = !!reference_date, end_date = !!source_vars[[1]])
      )
  }
}
