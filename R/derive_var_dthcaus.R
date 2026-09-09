#' Derive Death Cause
#'
#' @description
#' `r lifecycle::badge("deprecated")` The `derive_var_dthcaus()`
#' function has been deprecated in favor of `derive_vars_extreme_event()`.
#'
#' Derive death cause (`DTHCAUS`) and add traceability variables if required.
#'
#' @param dataset
#'   `r roxygen_param_dataset(expected_vars = c("subject_keys"))`
#'
#' @param source_datasets A named `list` containing datasets in which to search for the
#'   death cause
#'
#' @param ... Objects of class "dthcaus_source" created by [`dthcaus_source()`].
#'
#' @param subject_keys Variables to uniquely identify a subject
#'
#' A list of expressions where the expressions are symbols as returned by
#' `exprs()` is expected.
#'
#' @details
#' This function derives `DTHCAUS` along with the user-defined traceability
#' variables, if required. If a subject has death info from multiple sources,
#' the one from the source with the earliest death date will be used. If dates are
#' equivalent, the first source will be kept, so the user should provide the inputs in
#' the preferred order.
#'
#' @family deprecated
#' @keywords deprecated
#'
#' @return The input dataset with `DTHCAUS` variable added.
#'
#' @export
#'
#' @seealso [dthcaus_source()]
#'
derive_var_dthcaus <- function(dataset,
                               ...,
                               source_datasets,
                               subject_keys = get_admiral_option("subject_keys")) {
  deprecate_stop(
    when = "1.2.0",
    what = "derive_var_dthcaus()",
    with = "derive_vars_extreme_event()"
  )
}

#' Create a `dthcaus_source` Object
#'
#' @description
#' `r lifecycle::badge("deprecated")` The `dthcaus_source()`
#' function and `dthcaus_source()` have been deprecated in favor of
#' `event()`.
#'
#' @param dataset_name The name of the dataset, i.e. a string, used to search for
#'   the death cause.
#'
#' @param filter An expression used for filtering `dataset`.
#'
#' @param date A date or datetime variable or an expression to be used for
#'   sorting `dataset`.
#'
#' @param order Sort order
#'
#'   Additional variables/expressions to be used for sorting the `dataset`. The
#'   dataset is ordered by `date` and `order`. Can be used to avoid duplicate
#'   record warning.
#'
#' @permitted list of expressions created by `exprs()`, e.g.,
#'   `exprs(ADT, desc(AVAL))` or `NULL`
#'
#' @param mode One of `"first"` or `"last"`.
#' Either the `"first"` or `"last"` observation is preserved from the `dataset`
#' which is ordered by `date`.
#'
#' @param dthcaus A variable name, an expression, or a string literal
#'
#'   If a variable name is specified, e.g., `AEDECOD`, it is the variable in the
#'   source dataset to be used to assign values to `DTHCAUS`; if an expression,
#'   e.g., `str_to_upper(AEDECOD)`, it is evaluated in the source dataset and
#'   the results is assigned to `DTHCAUS`; if a string literal, e.g. `"Adverse
#'   Event"`, it is the fixed value to be assigned to `DTHCAUS`.
#'
#' @param set_values_to Variables to be set to trace the source dataset
#'
#' @family deprecated
#' @keywords deprecated
#'
#'
#' @export
#'
#' @seealso [derive_var_dthcaus()]
#'
#' @return An object of class "dthcaus_source".
#'
dthcaus_source <- function(dataset_name,
                           filter,
                           date,
                           order = NULL,
                           mode = "first",
                           dthcaus,
                           set_values_to = NULL) {
  deprecate_stop(
    when = "1.2.0",
    what = "dthcaus_source()",
    with = "event()"
  )
}
