#' Call a Single Derivation Multiple Times
#'
#' Call a single derivation multiple times with some parameters/arguments being fixed across
#' iterations and others varying.
#'
#' @param dataset  `r roxygen_param_dataset()`
#' @param derivation The derivation function to call
#'
#'   A function that performs a specific derivation is expected. A derivation
#'   adds variables or observations to a dataset. The first argument of a
#'   derivation must expect a dataset and the derivation must return a dataset.
#'   All expected arguments for the derivation function must be provided through
#'   the `params()` objects passed to the `variable_params` and `...` arguments.
#'
#' @param variable_params A `list` of function arguments that are different across iterations.
#'   Each set of function arguments must be created using [`params()`].
#' @param ... Any number of *named* function arguments that stay the same across iterations.
#'   If a function argument is specified both inside `variable_params` and `...` then
#'   the value in `variable_params` overwrites the one in `...`.
#'
#'  @details
#'
#'   It is also possible to pass functions from outside the `{admiral}` package
#'   to `call_derivation()`, e.g. an extension package function, or
#'   `dplyr::mutate()`. The only requirement for a function being passed to `derivation` is that
#'   it must take a dataset as its first argument and return a dataset.
#'
#' @return
#' The input dataset with additional records/variables added depending on
#' which `derivation` has been used.
#'
#' @family high_order_function
#' @keywords high_order_function
#'
#' @export
#'
#' @seealso [params()] [restrict_derivation()] [call_derivation()]
#'
#' @examples
#' library(dplyr, warn.conflicts = FALSE)
#' adsl <- tribble(
#'   ~STUDYID,   ~USUBJID,      ~TRTSDT,      ~TRTEDT,
#'   "PILOT01", "01-1307",           NA,           NA,
#'   "PILOT01", "05-1377", "2014-01-04", "2014-01-25",
#'   "PILOT01", "06-1384", "2012-09-15", "2012-09-24",
#'   "PILOT01", "15-1085", "2013-02-16", "2013-08-18",
#'   "PILOT01", "16-1298", "2013-04-08", "2013-06-28"
#' ) %>%
#'   mutate(
#'     across(TRTSDT:TRTEDT, as.Date)
#'   )
#'
#' ae <- tribble(
#'   ~STUDYID,  ~DOMAIN,  ~USUBJID,     ~AESTDTC,     ~AEENDTC,
#'   "PILOT01",    "AE", "06-1384", "2012-09-15", "2012-09-29",
#'   "PILOT01",    "AE", "06-1384", "2012-09-15", "2012-09-29",
#'   "PILOT01",    "AE", "06-1384", "2012-09-23", "2012-09-29",
#'   "PILOT01",    "AE", "06-1384", "2012-09-23", "2012-09-29",
#'   "PILOT01",    "AE", "06-1384", "2012-09-15", "2012-09-29",
#'   "PILOT01",    "AE", "06-1384", "2012-09-15", "2012-09-29",
#'   "PILOT01",    "AE", "06-1384", "2012-09-15", "2012-09-29",
#'   "PILOT01",    "AE", "06-1384", "2012-09-15", "2012-09-29",
#'   "PILOT01",    "AE", "06-1384", "2012-09-23", "2012-09-29",
#'   "PILOT01",    "AE", "06-1384", "2012-09-23", "2012-09-29",
#'   "PILOT01",    "AE", "16-1298", "2013-06-08", "2013-07-06",
#'   "PILOT01",    "AE", "16-1298", "2013-06-08", "2013-07-06",
#'   "PILOT01",    "AE", "16-1298", "2013-04-22", "2013-07-06",
#'   "PILOT01",    "AE", "16-1298", "2013-04-22", "2013-07-06",
#'   "PILOT01",    "AE", "16-1298", "2013-04-22", "2013-07-06",
#'   "PILOT01",    "AE", "16-1298", "2013-04-22", "2013-07-06"
#' )
#'
#' adae <- ae %>%
#'   derive_vars_merged(
#'     dataset_add = adsl,
#'     new_vars = exprs(TRTSDT, TRTEDT),
#'     by_vars = exprs(USUBJID)
#'   )
#'
#' ## While `derive_vars_dt()` can only add one variable at a time, using `call_derivation()`
#' ## one can add multiple variables in one go
#' call_derivation(
#'   dataset = adae,
#'   derivation = derive_vars_dt,
#'   variable_params = list(
#'     params(dtc = AESTDTC, date_imputation = "first", new_vars_prefix = "AST"),
#'     params(dtc = AEENDTC, date_imputation = "last", new_vars_prefix = "AEN")
#'   ),
#'   min_dates = exprs(TRTSDT),
#'   max_dates = exprs(TRTEDT)
#' )
#'
#' ## The above call using `call_derivation()` is equivalent to the following
#' adae %>%
#'   derive_vars_dt(
#'     new_vars_prefix = "AST",
#'     dtc = AESTDTC,
#'     date_imputation = "first",
#'     min_dates = exprs(TRTSDT),
#'     max_dates = exprs(TRTEDT)
#'   ) %>%
#'   derive_vars_dt(
#'     new_vars_prefix = "AEN",
#'     dtc = AEENDTC,
#'     date_imputation = "last",
#'     min_dates = exprs(TRTSDT),
#'     max_dates = exprs(TRTEDT)
#'   )
call_derivation <- function(dataset = NULL, derivation, variable_params, ...) {
  assert_data_frame(dataset, optional = TRUE)
  assert_s3_class(derivation, "function")
  assert_list_of(variable_params, "params")

  fixed_params <- enexprs(...)
  if (length(fixed_params) == 0L) {
    cli_abort("At least one argument must be set inside {.arg ...}.")
  }
  if (!is_named(fixed_params)) {
    cli_abort("All arguments inside {.arg ...} must be named.")
  }

  all_params <- union(unlist(map(variable_params, names)), names(fixed_params))
  assert_function(derivation, all_params)

  for (i in seq_along(variable_params)) {
    fixed_params_ <- fixed_params[names(fixed_params) %notin% names(variable_params[[i]])]
    call <- call2(derivation, expr(dataset), !!!variable_params[[i]], !!!fixed_params_)
    fixed_env <- caller_env()
    variable_env <- attr(variable_params[[i]], "env")
    if (!identical(fixed_env, variable_env)) {
      # prefer objects in the variable environment to objects in fixed environment
      # Note: objects in any of the parent environments of the variable environment are ignored.
      eval_env <- new_environment(
        data = c(list(dataset = dataset), as.list(variable_env)),
        parent = fixed_env
      )
    } else {
      eval_env <- new_environment(
        data = list(dataset = dataset),
        parent = fixed_env
      )
    }

    dataset <- eval_tidy(call, env = eval_env)
  }

  dataset
}

#' Create a Set of Parameters
#'
#' Create a set of variable parameters/function arguments to be used in [`call_derivation()`].
#'
#' @param ... One or more named arguments
#'
#'
#' @return An object of class `params`
#'
#' @family other_advanced
#' @keywords other_advanced
#'
#' @export
#'
#' @seealso [call_derivation()]
#'
#' @examples
#' library(dplyr, warn.conflicts = FALSE)
#'
#' adsl <- tribble(
#'   ~STUDYID,   ~USUBJID,      ~TRTSDT,      ~TRTEDT,
#'   "PILOT01", "01-1307",           NA,           NA,
#'   "PILOT01", "05-1377", "2014-01-04", "2014-01-25",
#'   "PILOT01", "06-1384", "2012-09-15", "2012-09-24",
#'   "PILOT01", "15-1085", "2013-02-16", "2013-08-18",
#'   "PILOT01", "16-1298", "2013-04-08", "2013-06-28"
#' ) %>%
#'   mutate(
#'     across(TRTSDT:TRTEDT, as.Date)
#'   )
#'
#' ae <- tribble(
#'   ~STUDYID,  ~DOMAIN,  ~USUBJID,     ~AESTDTC,     ~AEENDTC,
#'   "PILOT01",    "AE", "06-1384", "2012-09-15", "2012-09-29",
#'   "PILOT01",    "AE", "06-1384", "2012-09-15", "2012-09-29",
#'   "PILOT01",    "AE", "06-1384", "2012-09-23", "2012-09-29",
#'   "PILOT01",    "AE", "06-1384", "2012-09-23", "2012-09-29",
#'   "PILOT01",    "AE", "06-1384", "2012-09-15", "2012-09-29",
#'   "PILOT01",    "AE", "06-1384", "2012-09-15", "2012-09-29",
#'   "PILOT01",    "AE", "06-1384", "2012-09-15", "2012-09-29",
#'   "PILOT01",    "AE", "06-1384", "2012-09-15", "2012-09-29",
#'   "PILOT01",    "AE", "06-1384", "2012-09-23", "2012-09-29",
#'   "PILOT01",    "AE", "06-1384", "2012-09-23", "2012-09-29",
#'   "PILOT01",    "AE", "16-1298", "2013-06-08", "2013-07-06",
#'   "PILOT01",    "AE", "16-1298", "2013-06-08", "2013-07-06",
#'   "PILOT01",    "AE", "16-1298", "2013-04-22", "2013-07-06",
#'   "PILOT01",    "AE", "16-1298", "2013-04-22", "2013-07-06",
#'   "PILOT01",    "AE", "16-1298", "2013-04-22", "2013-07-06",
#'   "PILOT01",    "AE", "16-1298", "2013-04-22", "2013-07-06"
#' )
#'
#' adae <- ae %>%
#'   select(USUBJID, AESTDTC, AEENDTC) %>%
#'   derive_vars_merged(
#'     dataset_add = adsl,
#'     new_vars = exprs(TRTSDT, TRTEDT),
#'     by_vars = exprs(USUBJID)
#'   )
#'
#' ## In order to derive both `ASTDT` and `AENDT` in `ADAE`, one can use `derive_vars_dt()`
#' adae %>%
#'   derive_vars_dt(
#'     new_vars_prefix = "AST",
#'     dtc = AESTDTC,
#'     date_imputation = "first",
#'     min_dates = exprs(TRTSDT),
#'     max_dates = exprs(TRTEDT)
#'   ) %>%
#'   derive_vars_dt(
#'     new_vars_prefix = "AEN",
#'     dtc = AEENDTC,
#'     date_imputation = "last",
#'     min_dates = exprs(TRTSDT),
#'     max_dates = exprs(TRTEDT)
#'   )
#'
#'
#' ## While `derive_vars_dt()` can only add one variable at a time, using `call_derivation()`
#' ## one can add multiple variables in one go.
#' ## The function arguments which are different from a variable to another (e.g. `new_vars_prefix`,
#' ## `dtc`, and `date_imputation`) are specified as a list of `params()` in the `variable_params`
#' ## argument of `call_derivation()`. All other arguments which are common to all variables
#' ## (e.g. `min_dates` and `max_dates`) are specified outside of `variable_params` (i.e. in `...`).
#' call_derivation(
#'   dataset = adae,
#'   derivation = derive_vars_dt,
#'   variable_params = list(
#'     params(dtc = AESTDTC, date_imputation = "first", new_vars_prefix = "AST"),
#'     params(dtc = AEENDTC, date_imputation = "last", new_vars_prefix = "AEN")
#'   ),
#'   min_dates = exprs(TRTSDT),
#'   max_dates = exprs(TRTEDT)
#' )
#'
#' ## The above call using `call_derivation()` is equivalent to the call using `derive_vars_dt()`
#' ## to derive variables `ASTDT` and `AENDT` separately at the beginning.
params <- function(...) {
  args <- enexprs(...)
  attr(args, "env") <- caller_env()
  if (length(args) == 0L) {
    cli_abort("At least one argument must be provided.")
  }
  if (!is_named(args)) {
    cli_abort("All arguments passed to {.fun params} must be named.")
  }
  duplicate_params <- get_duplicates(names(args))
  if (length(duplicate_params) >= 1L) {
    cli_abort(paste(
      "The following argument{?s} {?has/have} been specified more than once:",
      "{.val {duplicate_params}}."
    ))
  }
  structure(args, class = c("params", "source", "list")) # nolint: undesirable_function_linter
}
