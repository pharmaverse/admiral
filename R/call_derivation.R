#' Call a Single Derivation Multiple Times
#'
#' Call a single derivation multiple times with some parameters being fixed across
#' iterations and others varying.
#'
#' @param dataset The input dataset
#' @param derivation The derivation function to call
#' @param variable_params A `list` of arguments that are different across iterations.
#'   Each set of arguments must be created using [`params()`].
#' @param ... Any number of *named* arguments that are fixed across iterations.
#'   If a parameter is specified both inside `variable_params` and `...` then
#'   the value in `variable_params` overwrites the one in `...`
#'
#' @author Thomas Neitmann, Stefan Bundfuss
#'
#' @return
#' The input dataset with additional records/variables added depending on
#' which `derivation` has been used.
#'
#' @keywords user_utility
#'
#' @export
#'
#' @seealso params
#'
#' @examples
#' library(dplyr, warn.conflicts = FALSE)
#' data(ae)
#' data(adsl)
#'
#' adae <- ae[sample(1:nrow(ae), 1000), ] %>%
#'   left_join(adsl, by = "USUBJID") %>%
#'   select(USUBJID, AESTDTC, AEENDTC, TRTSDT, TRTEDT)
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
#'   min_dates = list(TRTSDT),
#'   max_dates = list(TRTEDT)
#' )
#'
#' ## The above call using `call_derivation()` is equivalent to the following
#' adae %>%
#'   derive_vars_dt(
#'     new_vars_prefix = "AST",
#'     dtc = AESTDTC,
#'     date_imputation = "first",
#'     min_dates = list(TRTSDT),
#'     max_dates = list(TRTEDT)
#'   ) %>%
#'   derive_vars_dt(
#'     new_vars_prefix = "AEN",
#'     dtc = AEENDTC,
#'     date_imputation = "last",
#'     min_dates = list(TRTSDT),
#'     max_dates = list(TRTEDT)
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
#' Create a set of variable parameters to be used in an iteration of [`call_derivation()`]
#'
#' @param ... One or more named arguments
#'
#' @author Thomas Neitmann
#'
#' @return An object of class `params`
#'
#' @keywords source_specifications
#'
#' @export
#'
#' @examples
#' params(
#'   fns = list(VSSTRESN ~ mean(., na.rm = TRUE)),
#'   set_values_to = vars(DTYPE = "AVERAGE")
#' )
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
  structure(args, class = c("params", "list"))
}
