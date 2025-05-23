#' Get Duplicate Records that Led to a Prior Error
#'
#' @export
#'
#'
#' @details
#' Many `{admiral}` function check that the input dataset contains only one record
#' per `by_vars` group and throw an error otherwise. The `get_duplicates_dataset()`
#' function allows one to retrieve the duplicate records that lead to an error.
#'
#' Note that the function always returns the dataset of duplicates from the last
#' error that has been thrown in the current R session. Thus, after restarting the
#' R sessions `get_duplicates_dataset()` will return `NULL` and after a second error
#' has been thrown, the dataset of the first error can no longer be accessed (unless
#' it has been saved in a variable).
#'
#' @return A `data.frame` or `NULL`
#'
#' @family utils_ds_chk
#'
#' @keywords utils_ds_chk
#'
#' @examples
#' data(admiral_adsl)
#'
#' # Duplicate the first record
#' adsl <- rbind(admiral_adsl[1L, ], admiral_adsl)
#'
#' signal_duplicate_records(adsl, exprs(USUBJID), cnd_type = "warning")
#'
#' get_duplicates_dataset()
get_duplicates_dataset <- function() {
  admiral_environment$duplicates
}

#' Extract Duplicate Records
#'
#' @param dataset `r roxygen_param_dataset(expected_vars = c("by_vars"))`
#' @param by_vars Grouping variables
#'
#'  Defines groups of records in which to look for duplicates.
#'  If omitted, all variables in the input dataset are used in the by group.
#'
#'  **Note:**  Omitting `by_vars` will increase the function's run-time, so it is
#'  recommended to specify the necessary grouping variables for large datasets
#'  whenever possible.
#'
#' `r roxygen_param_by_vars()`
#'
#' @permitted [var_list]
#'
#' @return A `data.frame` of duplicate records within `dataset`
#'
#' @export
#' @family internal
#' @keywords internal
#'
#' @examples
#' data(admiral_adsl)
#'
#' # Duplicate the first record
#' adsl <- rbind(admiral_adsl[1L, ], admiral_adsl)
#'
#' extract_duplicate_records(adsl, exprs(USUBJID))
extract_duplicate_records <- function(dataset, by_vars = NULL) {
  if (is.null(by_vars)) {
    assert_data_frame(dataset, check_is_grouped = FALSE)
    by_vars <- exprs(!!!parse_exprs(names(dataset)))
  } else {
    assert_expr_list(by_vars)
    assert_data_frame(dataset, required_vars = extract_vars(by_vars), check_is_grouped = FALSE)
  }

  data_by <- dataset %>%
    ungroup() %>%
    # evaluate expressions in by_vars
    mutate(!!!by_vars,
      .keep = "none"
    )

  is_duplicate <- duplicated(data_by) | duplicated(data_by, fromLast = TRUE)

  dataset %>%
    ungroup() %>%
    # evaluate expressions in by_vars
    mutate(!!!by_vars) %>%
    # move by variables to the beginning
    # if by_vars includes unnamed expressions, the unevaluated expression is
    # used as variable name
    select(!!!syms(map(replace_values_by_names(by_vars), as_label)), everything()) %>%
    filter(is_duplicate) %>%
    arrange(!!!by_vars)
}

#' Signal Duplicate Records
#'
#' @param dataset `r roxygen_param_dataset(expected_vars = c("by_vars"))`
#' @param by_vars Grouping variables
#'
#'  Defines groups of records in which to look for duplicates.
#'
#' `r roxygen_param_by_vars()`
#'
#' @permitted [var_list]
#'
#' @param msg The condition message
#' @param cnd_type Type of condition to signal when detecting duplicate records.
#'
#' @permitted `"message"`, `"warning"`, or `"error"`
#'
#' @param class Class of the condition
#'
#'   The specified classes are added to the classes of the condition.
#'   `c("duplicate_records", "assert-admiral")` is always added.
#'
#' @return No return value, called for side effects
#'
#' @export
#' @family internal
#' @keywords internal
#'
#' @examples
#' data(admiral_adsl)
#'
#' # Duplicate the first record
#' adsl <- rbind(admiral_adsl[1L, ], admiral_adsl)
#'
#' signal_duplicate_records(adsl, exprs(USUBJID), cnd_type = "message")
signal_duplicate_records <- function(dataset,
                                     by_vars,
                                     msg = paste(
                                       "Dataset contains duplicate records",
                                       "with respect to",
                                       "{.var {replace_values_by_names(by_vars)}}"
                                     ),
                                     cnd_type = "error",
                                     class = NULL) {
  assert_expr_list(by_vars)
  assert_data_frame(dataset, required_vars = extract_vars(by_vars), check_is_grouped = FALSE)
  assert_character_vector(msg)
  assert_character_scalar(cnd_type, values = c("message", "warning", "error"))
  assert_character_vector(class, optional = TRUE)

  cnd_funs <- list(message = cli_inform, warning = cli_warn, error = cli_abort)

  duplicate_records <- extract_duplicate_records(dataset, by_vars)
  if (nrow(duplicate_records) >= 1L) {
    # nolint start: undesirable_function_linter
    admiral_environment$duplicates <- structure(
      duplicate_records,
      class = union("duplicates", class(duplicate_records)),
      by_vars = replace_values_by_names(by_vars)
    )
    # nolint end
    full_msg <- c(
      msg,
      i = "Run {.run admiral::get_duplicates_dataset()} to access the duplicate records"
    )
    cnd_funs[[cnd_type]](
      full_msg,
      class = c(class, "duplicate_records", "assert-admiral"),
      by_vars = by_vars
    )
  }
}

#' Print `duplicates` Objects
#'
#' @param x A `duplicates` object
#' @param ... Not used
#'
#' @return No return value, called for side effects
#'
#'
#' @keywords internal
#' @family utils_print
#'
#' @export
print.duplicates <- function(x, ...) {
  cli_text(
    "Duplicate records with respect to {.var {attr(x, \"by_vars\")}}."
  )
  NextMethod()
}
