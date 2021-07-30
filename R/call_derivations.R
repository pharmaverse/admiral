#' Call a Single Derivation Multiple Times
#'
#' Call a single derivation multiple times with some parameters being fixed across
#' iterations and others varying.
#'
#' @param dataset The input dataset
#' @param derivation The derivation function to call
#' @param variable_params A `list` of arguments that are different across iterations.
#'   Each set of arguments must be created using [`params()`].
#' @param ... Any number of arguments that are fixed across iterations
#'
#' @author Thomas Neitmann, Stefan Bundfuss
#'
#' @return
#' The input dataset with additional records/variables added depending on
#' which `derivation` has been used.
#'
#' @export
#' @seealso params
#'
#' @examples
#' library(dplyr, warn.conflicts = FLASE)
#' data(ae)
#' data(adsl)
#'
#' input <- ae[sample(1:nrow(ae), 1000), ] %>%
#'   left_join(adsl, by = "USUBJID")
#'
#' ## Call the same derivation twice in a row
#' expected_output <- input %>%
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
#'
#' ## Call the same derivation in one go
#' actual_output <- call_derivation(
#'   dataset = input,
#'   derivation = derive_vars_dt,
#'   variable_params = list(
#'     params(dtc = AESTDTC, date_imputation = "first", new_vars_prefix = "AST"),
#'     params(dtc = AEENDTC, date_imputation = "last", new_vars_prefix = "AEN")
#'   ),
#'   min_dates = list(TRTSDT),
#'   max_dates = list(TRTEDT)
#' )
call_derivation <- function(dataset, derivation, variable_params, ...) {
  assert_data_frame(dataset)
  assert_s3_class(derivation, "function")
  assert_list_of(variable_params, "params")

  fixed_params <- eval(substitute(alist(...)))
  if (!is_named(fixed_params)) {
    abort("All arguments inside `...` must be named")
  }

  all_params <- base::union(unlist(map(variable_params, names)), names(fixed_params))
  assert_function_param(deparse(substitute(derivation)), all_params)

  for (i in seq_along(variable_params)) {
    args <- c(quote(dataset), variable_params[[i]], fixed_params)
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
#' @return An object of class `params`
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
  structure(args, class = c("params", "list"))
}
