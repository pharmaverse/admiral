#' Derive First or Last Datetime from Multiple Sources
#'
#' @description
#' `r lifecycle::badge("deprecated")`
#'
#' The `derive_var_extreme_dtm()`
#' function has been deprecated in favor of `derive_vars_extreme_event()`.
#'
#' Add the first or last datetime from multiple sources to the dataset, e.g.,
#' the last known alive datetime (`LSTALVDTM`).
#'
#' @param dataset
#'   `r roxygen_param_dataset(expected_vars = c("subject_keys"))`
#'
#' @param new_var Name of variable to create
#'
#' @param source_datasets A named `list` containing datasets in which to search
#'   for the first or last date
#'
#' @param ... Source(s) of dates. One or more `date_source()` objects are
#'   expected.
#'
#' @param mode Selection mode (first or last)
#'
#'   If `"first"` is specified, the first date for each subject is selected. If
#'   `"last"` is specified, the last date for each subject is selected.
#'
#' @permitted  `"first"`, `"last"`
#'
#' @param subject_keys Variables to uniquely identify a subject
#'
#'   A list of expressions where the expressions are symbols as returned by
#'   `exprs()` is expected.
#'
#' @details The following steps are performed to create the output dataset:
#'
#'   1. For each source dataset the observations as specified by the `filter`
#'   element are selected and observations where `date` is `NA` are removed.
#'   Then for each patient the first or last observation (with respect to `date`
#'   and `mode`) is selected.
#'
#'   1. The new variable is set to the variable or expression specified by the
#'   `date` element. If this is a date variable (rather than datetime), then the
#'   time is imputed as `"00:00:00"`.
#'
#'   1. The variables specified by the `set_values_to` element are added.
#'
#'   1. The selected observations of all source datasets are combined into a
#'   single dataset.
#'
#'   1. For each patient the first or last observation (with respect to the new
#'   variable and `mode`) from the single dataset is selected and the new
#'   variable is merged to the input dataset.
#'
#' @return The input dataset with the new variable added.
#'
#'
#' @family deprecated
#' @keywords deprecated
#'
#' @seealso [date_source()], [derive_var_extreme_dt()],
#'   [derive_vars_merged()]
#'
#' @export
#'
derive_var_extreme_dtm <- function(dataset,
                                   new_var,
                                   ...,
                                   source_datasets,
                                   mode,
                                   subject_keys = get_admiral_option("subject_keys")) {
  deprecate_stop(
    when = "1.2.0",
    what = "derive_var_extreme_dtm()",
    with = "derive_vars_extreme_event()"
  )
}

#' Derive First or Last Date from Multiple Sources
#'
#' @description
#' `r lifecycle::badge("deprecated")` The `derive_var_extreme_dt()`
#' function has been deprecated in favor of `derive_vars_extreme_event()`.
#'
#' Add the first or last date from multiple sources to the
#' dataset, e.g., the last known alive date (`LSTALVDT`).
#'
#' **Note:** This is a wrapper function for the function `derive_var_extreme_dtm()`.
#'
#' @inheritParams derive_var_extreme_dtm
#'
#' @details The following steps are performed to create the output dataset:
#'
#'   1. For each source dataset the observations as specified by the `filter`
#'   element are selected and observations where `date` is `NA` are removed.
#'   Then for each patient the first or last observation (with respect to `date`
#'   and `mode`) is selected.
#'
#'   1. The new variable is set to the variable or expression specified by the
#'   `date` element.
#'
#'   1. The variables specified by the `set_values_to` element are added.
#'
#'   1. The selected observations of all source datasets are combined into a
#'   single dataset.
#'
#'   1. For each patient the first or last observation (with respect to the new
#'   variable and `mode`) from the single dataset is selected and the new
#'   variable is merged to the input dataset.
#'
#'   1. The time part is removed from the new variable.
#'
#' @return The input dataset with the new variable added.
#'
#'
#' @family deprecated
#' @keywords deprecated
#'
#' @seealso [date_source()], [derive_var_extreme_dtm()], [derive_vars_merged()]
#'
#' @export
#'
derive_var_extreme_dt <- function(dataset,
                                  new_var,
                                  ...,
                                  source_datasets,
                                  mode,
                                  subject_keys = get_admiral_option("subject_keys")) {
  deprecate_stop(
    when = "1.2.0",
    what = "derive_var_extreme_dt()",
    with = "derive_vars_extreme_event()"
  )
}

#' Create a `date_source` object
#'
#' @description
#' `r lifecycle::badge("deprecated")` The `date_source()`
#' function has been deprecated in favor of `event()`.
#'
#' Create a `date_source` object as input for `derive_var_extreme_dt()` and
#' `derive_var_extreme_dtm()`.
#'
#' @param dataset_name The name of the dataset, i.e. a string, used to search for
#'   the date.
#'
#' @param filter An unquoted condition for filtering `dataset`.
#'
#' @param date A variable or an expression providing a date. A date or a
#'   datetime can be specified. An unquoted symbol or expression is expected.
#'
#' @param set_values_to Variables to be set
#'
#' @seealso [derive_var_extreme_dtm()], [derive_var_extreme_dt()]
#'
#' @family deprecated
#' @keywords deprecated
#'
#' @export
#'
#' @return An object of class `date_source`.
#'
date_source <- function(dataset_name,
                        filter = NULL,
                        date,
                        set_values_to = NULL) {
  deprecate_stop(
    when = "1.2.0",
    what = "date_source()",
    with = "event()"
  )
}
