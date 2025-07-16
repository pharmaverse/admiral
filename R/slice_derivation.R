#' Execute a Derivation with Different Arguments for Subsets of the Input Dataset
#'
#' The input dataset is split into slices (subsets) and for each slice the
#' derivation is called separately. Some or all arguments of the derivation
#' may vary depending on the slice.
#'
#' @param dataset `r roxygen_param_dataset()`
#'
#' @param derivation Derivation
#'
#'   A function that performs a specific derivation is expected. A derivation
#'   adds variables or observations to a dataset. The first argument of a
#'   derivation must expect a dataset and the derivation must return a dataset.
#'   All expected arguments for the derivation function must be provided through
#'   the `params()` object passed to the `args` argument or be provided in _every_
#'   `derivation_slice()`.
#'
#' @param args Arguments of the derivation
#'
#'   A `params()` object is expected.
#'
#' @param ... A `derivation_slice()` object is expected
#'
#'   Each slice defines a subset of the input dataset and some of the parameters
#'   for the derivation. The derivation is called on the subset with the
#'   parameters specified by the `args` parameter and the `args` field of the
#'   `derivation_slice()` object. If a parameter is specified for both, the
#'   value in `derivation_slice()` overwrites the one in `args`.
#'
#' @details
#'
#'   For each slice the derivation is called on the subset defined by the
#'   `filter` field of the `derivation_slice()` object and with the parameters
#'   specified by the `args` parameter and the `args` field of the
#'   `derivation_slice()` object. If a parameter is specified for both, the
#'   value in `derivation_slice()` overwrites the one in `args`.
#'
#'   - Observations that match with more than one slice are only considered for
#'   the first matching slice.
#'
#'   - The derivation is called for slices with no observations.
#'
#'   - Observations with no match to any of the slices are included in the
#'   output dataset but the derivation is not called for them.
#'
#'   It is also possible to pass functions from outside the `{admiral}` package
#'   to `slice_derivation()`, e.g. an extension package function, or
#'   `dplyr::mutate()`. The only requirement for a function being passed to `derivation` is that
#'   it must take a dataset as its first argument and return a dataset.
#'
#' @return The input dataset with the variables derived by the derivation added
#'
#' @family high_order_function
#' @keywords high_order_function
#'
#'
#' @seealso [params()] [restrict_derivation()] [call_derivation()]
#'
#' @export
#'
#' @examples
#' library(tibble)
#' library(stringr)
#' advs <- tribble(
#'   ~USUBJID, ~VSDTC,       ~VSTPT,
#'   "1",      "2020-04-16", NA_character_,
#'   "1",      "2020-04-16", "BEFORE TREATMENT"
#' )
#'
#' # For the second slice filter is set to TRUE. Thus derive_vars_dtm is called
#' # with time_imputation = "last" for all observations which do not match for the
#' # first slice.
#' slice_derivation(
#'   advs,
#'   derivation = derive_vars_dtm,
#'   args = params(
#'     dtc = VSDTC,
#'     new_vars_prefix = "A"
#'   ),
#'   derivation_slice(
#'     filter = str_detect(VSTPT, "PRE|BEFORE"),
#'     args = params(time_imputation = "first")
#'   ),
#'   derivation_slice(
#'     filter = TRUE,
#'     args = params(time_imputation = "last")
#'   )
#' )
#'
slice_derivation <- function(dataset,
                             derivation,
                             args = NULL,
                             ...) {
  # check input
  assert_data_frame(dataset)
  assert_function(derivation)
  assert_s3_class(args, "params", optional = TRUE)
  if (!is.null(args)) {
    assert_function(derivation, names(args))
  }
  slices <- list2(...)
  assert_list_of(slices, "derivation_slice")

  # Check that every mandatory argument to the derivation function is either passed
  # inside args or present in every slice
  mandatory_args <- formals(derivation) %>%
    keep(is_missing) %>%
    names()

  if (length(mandatory_args > 1)) {
    # ignore first item, as that is the dataset
    mandatory_args <- mandatory_args[-1]
    # ignore the ... argument too, if present
    mandatory_args <- mandatory_args[mandatory_args != "..."]

    if (length(mandatory_args > 1)) {
      for (mandatory_arg in mandatory_args) {
        if (!mandatory_arg %in% names(args)) {
          check <- vapply(
            slices,
            function(x) mandatory_arg %in% names(x$args),
            FUN.VALUE = length(slices)
          )

          if (!all(check)) {
            cli_abort(
              "Issue with the mandatory argument `{mandatory_arg}` of derivation
              function `{deparse(substitute(derivation))}`. It must: (1) be passed
              to the `args` argument, or (2) be passed to all derivation slices."
            )
          }
        }
      }
    }
  }

  # the variable temp_slicenr is added to the dataset which indicates to which
  # slice the observation belongs. Observations which match to more than one
  # slice are assigned to the first matching slice.
  # Observations which does not match to any slice are assigned NA. For these
  # the derivation is not called.
  cases <- vector("list", length(slices))
  for (i in seq_along(slices)) {
    cases[[i]] <- expr(!!slices[[i]]$filter ~ !!i)
  }
  slice_call <- call2("case_when", !!!cases)
  dataset <- mutate(
    dataset,
    temp_slicenr = !!slice_call
  )
  ret <- list()
  # call derivation for each slice
  for (i in seq_along(slices)) {
    # call derivation on subset
    # remove global arguments which were specified by the slice
    act_args <- args[names(args) %notin% names(slices[[i]]$args)]

    call <- call2(derivation, expr(data), !!!act_args, !!!slices[[i]]$args)
    obsnr <- which(dataset$temp_slicenr == i)
    # call the derivation for non-empty slices only
    # create environment in which the call to the derivation is evaluated
    act_env <- attr(args, "env")
    slice_env <- attr(slices[[i]]$args, "env")
    if (!identical(act_env, slice_env)) {
      # prefer objects in the slice environment to object in args environment
      # Note: objects in any of the parent environments of the slice environment are ignored.
      eval_env <- new_environment(
        data = c(list(data = dataset[obsnr, , drop = FALSE]), as.list(slice_env)),
        parent = act_env
      )
    } else {
      eval_env <- new_environment(
        data = list(data = dataset[obsnr, , drop = FALSE]),
        parent = act_env
      )
    }

    ret[[i]] <- eval_tidy(call, env = eval_env)
  }
  # put datasets together again
  bind_rows(ret, dataset[is.na(dataset$temp_slicenr), , drop = FALSE]) %>%
    select(-temp_slicenr)
}

#' Create a `derivation_slice` Object
#'
#' Create a `derivation_slice` object as input for `slice_derivation()`.
#'
#' @param filter An unquoted condition for defining the observations of the
#'   slice
#'
#' @param args Arguments of the derivation to be used for the slice
#'
#'   A `params()` object is expected.
#'
#' @return An object of class `derivation_slice`

#'
#' @seealso [slice_derivation()], [params()]
#'
#' @family high_order_function
#' @keywords high_order_function
#'
#' @export
derivation_slice <- function(filter,
                             args = NULL) {
  out <- list(
    filter = assert_filter_cond(enexpr(filter)),
    args = assert_s3_class(args, "params", optional = TRUE)
  )
  class(out) <- c("derivation_slice", "source", "list")
  out
}
