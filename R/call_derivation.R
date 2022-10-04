#' Call a Single Derivation Multiple Times
#'
#' Call a single derivation multiple times with some parameters/arguments being fixed across
#' iterations and others varying.
#'
#' @param dataset The input dataset
#' @param derivation The derivation function to call
#' @param variable_params A `list` of function arguments that are different across iterations.
#'   Each set of function arguments must be created using [`params()`].
#' @param ... Any number of *named* function arguments that stay the same across iterations.
#'   If a function argument is specified both inside `variable_params` and `...` then
#'   the value in `variable_params` overwrites the one in `...`
#'
#' @author Thomas Neitmann, Stefan Bundfuss, Tracey Wang
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
#' @seealso [params()]
#'
#' @examples
#' library(dplyr, warn.conflicts = FALSE)
#' library(admiral.test)
#' data(admiral_ae)
#' data(admiral_adsl)
#'
#' adae <-
#'   select(admiral_ae[sample(1:nrow(admiral_ae), 1000), ], USUBJID, AESTDTC, AEENDTC) %>%
#'   derive_vars_merged(
#'     dataset_add = admiral_adsl,
#'     new_vars = vars(TRTSDT, TRTEDT),
#'     by_vars = vars(USUBJID)
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
#'   min_dates = vars(TRTSDT),
#'   max_dates = vars(TRTEDT)
#' )
#'
#' ## The above call using `call_derivation()` is equivalent to the following
#' adae %>%
#'   derive_vars_dt(
#'     new_vars_prefix = "AST",
#'     dtc = AESTDTC,
#'     date_imputation = "first",
#'     min_dates = vars(TRTSDT),
#'     max_dates = vars(TRTEDT)
#'   ) %>%
#'   derive_vars_dt(
#'     new_vars_prefix = "AEN",
#'     dtc = AEENDTC,
#'     date_imputation = "last",
#'     min_dates = vars(TRTSDT),
#'     max_dates = vars(TRTEDT)
#'   )
call_derivation <- function(dataset = NULL, derivation, variable_params, ...) {
  assert_data_frame(dataset, optional = TRUE)
  assert_s3_class(derivation, "function")
  assert_list_of(variable_params, "params")

  fixed_params <- eval(substitute(alist(...)))
  if (length(fixed_params) == 0L) {
    abort("At least one argument must be set inside `...`")
  }
  if (!is_named(fixed_params)) {
    abort("All arguments inside `...` must be named")
  }

  all_params <- base::union(unlist(map(variable_params, names)), names(fixed_params))
  assert_function_param(deparse(substitute(derivation)), all_params)

  for (i in seq_along(variable_params)) {
    fixed_params_ <- fixed_params[names(fixed_params) %notin% names(variable_params[[i]])]
    args <- c(quote(dataset), variable_params[[i]], fixed_params_)
    call <- as.call(c(substitute(derivation), args))
    dataset <- eval(call, envir = list(dataset = dataset), enclos = parent.frame())
  }

  dataset
}

#' Create a Set of Parameters
#'
#' Create a set of variable parameters/function arguments to be used in [`call_derivation()`].
#'
#' @param ... One or more named arguments
#'
#' @author Thomas Neitmann, Tracey Wang
#'
#' @return An object of class `params`
#'
#' @family source_specifications
#' @keywords source_specifications
#'
#' @export
#'
#' @seealso [call_derivation()]
#'
#' @examples
#' library(dplyr, warn.conflicts = FALSE)
#' library(admiral.test)
#' data(admiral_ae)
#' data(admiral_adsl)
#'
#' adae <- admiral_ae[sample(1:nrow(admiral_ae), 1000), ] %>%
#'   select(USUBJID, AESTDTC, AEENDTC) %>%
#'   derive_vars_merged(
#'     dataset_add = admiral_adsl,
#'     new_vars = vars(TRTSDT, TRTEDT),
#'     by_vars = vars(USUBJID)
#'   )
#'
#' ## In order to derive both `ASTDT` and `AENDT` in `ADAE`, one can use `derive_vars_dt()`
#' adae %>%
#'   derive_vars_dt(
#'     new_vars_prefix = "AST",
#'     dtc = AESTDTC,
#'     date_imputation = "first",
#'     min_dates = vars(TRTSDT),
#'     max_dates = vars(TRTEDT)
#'   ) %>%
#'   derive_vars_dt(
#'     new_vars_prefix = "AEN",
#'     dtc = AEENDTC,
#'     date_imputation = "last",
#'     min_dates = vars(TRTSDT),
#'     max_dates = vars(TRTEDT)
#'   )
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
#'   min_dates = vars(TRTSDT),
#'   max_dates = vars(TRTEDT)
#' )
#'
#' ## The above call using `call_derivation()` is equivalent to the call using `derive_vars_dt()`
#' ## to derive variables `ASTDT` and `AENDT` separately at the beginning.
params <- function(...) {
  args <- eval(substitute(alist(...)))
  if (length(args) == 0L) {
    abort("At least one argument must be provided")
  }
  if (!is_named(args)) {
    abort("All arguments passed to `params()` must be named")
  }
  duplicate_params <- get_duplicates(names(args))
  if (length(duplicate_params) >= 1L) {
    err_msg <- sprintf(
      "The following parameters have been specified more than once: %s",
      enumerate(duplicate_params)
    )
    abort(err_msg)
  }
  structure(args, class = c("params", "source", "list"))
}
