#' Execute a Derivation with Different Arguments for Subsets of the Input Dataset
#'
#' The input dataset is split into slices (subsets) and for each slice the
#' derivation is called separately. Some or all arguments of the derivation
#' may vary depending on the slice.
#'
#' @param dataset Input dataset
#'
#' @param derivation Derivation
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
#'   - Observations with no match to any of the slices are included in the
#'   output dataset but the derivation is not called for them.
#'
#' @return The input dataset with the variables derived by the derivation added
#'
#' @family high_order_function
#' @keywords high_order_function
#'
#' @author Stefan Bundfuss
#'
#' @seealso [params()] [restrict_derivation()]
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
slice_derivation <- function(dataset,
                             derivation,
                             args = NULL,
                             ...) {
  # check input
  assert_data_frame(dataset)
  assert_function(derivation, params = c("dataset"))
  assert_s3_class(args, "params", optional = TRUE)
  if (!is.null(args)) {
    assert_function_param(deparse(substitute(derivation)), names(args))
  }
  slices <- rlang::list2(...)
  assert_list_of(slices, "derivation_slice")

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

  # split dataset into slices
  dataset_split <- dataset %>%
    group_by(temp_slicenr) %>%
    nest()

  # call derivation for each slice
  for (i in seq_along(slices)) {
    # call derivation on subset
    # remove global arguments which were specified by the slice
    act_args <- args[names(args) %notin% names(slices[[i]]$args)]

    call <- as.call(c(substitute(derivation), c(quote(data), act_args, slices[[i]]$args)))
    obsnr <- which(dataset_split$temp_slicenr == i)
    if (length(obsnr) > 0) {
      # call the derivation for non-empty slices only
      dataset_split$data[[obsnr]] <-
        eval(call, envir = list(
          data = dataset_split$data[[obsnr]],
          enclos = parent.frame()
        ))
    }
  }

  # put datasets together again
  dataset_split %>%
    unnest(cols = c(data)) %>%
    ungroup() %>%
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

#' @author Stefan Bundfuss
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
    filter = assert_filter_cond(enquo(filter)),
    args = assert_s3_class(args, "params", optional = TRUE)
  )
  class(out) <- c("derivation_slice", "source", "list")
  out
}
