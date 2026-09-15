#' Adds a Parameter Based on First or Last Record from Multiple Sources
#'
#' @description
#' `r lifecycle::badge("deprecated")` The `derive_param_extreme_record()`
#' function has been deprecated in favor of `derive_extreme_event()`.
#'
#' Generates parameter based on the first or last observation from multiple
#' source datasets, based on user-defined filter, order and by group criteria.
#' All variables of the selected observation are kept.
#'
#' @param dataset `r roxygen_param_dataset()`
#'
#' @param sources Sources
#'
#'    A list of `records_source()` objects is expected.
#'
#' @param source_datasets Source datasets
#'
#'   A named list of datasets is expected. The `dataset_name` field of
#'   `records_source()` refers to the dataset provided in the list. The variables
#'   specified by the `order` and the `by_vars` arguments are expected after applying `new_vars`.
#'
#' @param by_vars Grouping variables
#'
#'   If the argument is specified, for each by group the observations are
#'   selected separately.
#'
#'   `r roxygen_param_by_vars()`
#'
#' @param order Sort order
#'
#'   If the argument is set to a non-null value, for each by group the first or
#'   last observation from the source datasets is selected with respect to
#'   the specified order. Variables created via `new_vars` e.g., imputed date variables,
#'   can be specified as well (see examples below).
#'
#'   Please note that `NA` is considered as the last value. I.e., if a order
#'   variable is `NA` and `mode = "last"`, this observation is chosen while for
#'   `mode = "first"` the observation is chosen only if there are no
#'   observations where the variable is not `NA`.
#'
#' @permitted list of expressions created by `exprs()`, e.g.,
#'   `exprs(ADT, desc(AVAL))`
#'
#' @param mode Selection mode (first or last)
#'
#'   If `"first"` is specified, for each by group the first observation with
#'   respect to `order` is included in the output dataset. If `"last"` is
#'   specified, the last observation is included in the output dataset.
#'
#' @permitted  `"first"`, `"last"`
#'
#' @param set_values_to Variables to be set
#'
#'   The specified variables are set to the specified values for the new
#'   observations.
#'
#'   A list of variable name-value pairs is expected.
#'   + LHS refers to a variable.
#'   + RHS refers to the values to set to the variable. This can be a string, a
#'   symbol, a numeric value or `NA`, e.g., `exprs(PARAMCD = "PD", PARAM =
#'   "First Progressive Disease")`.
#'
#' @details The following steps are performed to create the output dataset:
#'
#'   \enumerate{
#'   \item For each source dataset the observations as specified by
#'   the `filter` element are selected.
#'
#'   \item Variables specified by `new_vars` are created for each source dataset.
#'
#'   \item The first or last observation (with respect to the
#'   `order` variable) for each by group (specified by `by_vars`) from multiple sources
#'   is selected and added to the input dataset. }
#'
#' @return
#' The input dataset with the first or last observation of each by group
#' added as new observations.
#'
#' @family deprecated
#' @keywords deprecated
#'
#' @export
derive_param_extreme_record <- function(dataset = NULL,
                                        sources,
                                        source_datasets,
                                        by_vars = NULL,
                                        order,
                                        mode,
                                        set_values_to) {
  deprecate_stop(
    when = "1.2.0",
    what = "derive_param_extreme_record()",
    with = "derive_extreme_event()"
  )
}
