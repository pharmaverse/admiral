#' Turn a List of Quosures into a Character Vector
#'
#' @param quosures A `list` of `quosures` created using [`vars()`]
#'
#' @return A character vector
#'
#' @author Thomas Neitmann
#'
#' @export
#'
#' @keywords user_utility
#'
#' @examples
#' vars2chr(vars(USUBJID, AVAL))
vars2chr <- function(quosures) {
  rlang::set_names(
    map_chr(quosures, ~ as_string(quo_get_expr(.x))),
    names(quosures)
  )
}

#' Negate List of Variables
#'
#' The function adds a minus sign as prefix to each variable.
#'
#' This is useful if a list of variables should be removed from a dataset,
#' e.g., `select(!!!negate_vars(by_vars))` removes all by variables.
#'
#' @param vars List of variables created by `vars()`
#'
#' @return A list of `quosures`
#'
#' @author Stefan Bundfuss
#'
#' @export
#'
#' @keywords user_utility
#'
#' @examples
#' negate_vars(vars(USUBJID, STUDYID))
negate_vars <- function(vars = NULL) {
  assert_vars(vars, optional = TRUE)
  if (is.null(vars)) {
    NULL
  } else {
    lapply(vars, function(var) expr(-!!quo_get_expr(var)))
  }
}

#' Optional Filter
#'
#' Filters the input dataset if the provided expression is not `NULL`
#'
#' @param dataset Input dataset
#' @param filter A filter condition. Must be a quosure.
#'
#' @return A `data.frame` containing all rows in `dataset` matching `filter` or
#' just `dataset` if `filter` is `NULL`
#'
#' @author Thomas Neitmann
#'
#' @export
#'
#' @keywords user_utility
#'
#' @examples
#' library(admiraltest)
#' library(admiral)
#' data(admiral_vs)
#'
#' admiral::filter_if(admiral_vs, rlang::quo(NULL))
#' admiral::filter_if(admiral_vs, rlang::quo(VSTESTCD == "WEIGHT"))
filter_if <- function(dataset, filter) {
  assert_data_frame(dataset)
  assert_filter_cond(filter, optional = TRUE)

  if (quo_is_null(filter)) {
    dataset
  } else {
    filter(dataset, !!filter)
  }
}

#' Extract Unit From Parameter Description
#'
#' Extract the unit of a parameter from a description like "Param (unit)".
#'
#' @param x A parameter description
#'
#' @return A string
#'
#' @export
#'
#' @keywords user_utility
#'
#' @examples
#' extract_unit("Height (cm)")
#'
#' extract_unit("Diastolic Blood Pressure (mmHg)")
extract_unit <- function(x) {
  assert_character_vector(x)

  x %>%
    str_extract("\\(.+\\)") %>%
    str_remove_all("\\(|\\)")
}

#' Convert Blank Strings Into NAs
#'
#' Turn SAS blank strings into proper R `NA`s.
#'
#' @param x Any R object
#'
#' @details
#' The default methods simply returns its input unchanged. The `character` method
#' turns every instance of `""` into `NA_character_` while preserving *all* attributes.
#' When given a data frame as input the function keeps all non-character columns
#' as is and applies the just described logic to `character` columns. Once again
#' all attributes such as labels are preserved.
#'
#' @return An object of the same class as the input
#'
#' @author Thomas Neitmann
#'
#' @keywords user_utility
#'
#' @export
#'
#' @examples
#' convert_blanks_to_na(c("a", "b", "", "d", ""))
#'
#' df <- tibble::tibble(
#'   a = structure(c("a", "b", "", "c"), label = "A"),
#'   b = structure(c(1, NA, 21, 9), label = "B"),
#'   c = structure(c(TRUE, FALSE, TRUE, TRUE), label = "C"),
#'   d = structure(c("", "", "s", "q"), label = "D")
#' )
#' print(df)
#' convert_blanks_to_na(df)
convert_blanks_to_na <- function(x) {
  UseMethod("convert_blanks_to_na")
}

#' @export
#' @rdname convert_blanks_to_na
convert_blanks_to_na.default <- function(x) {
  x
}

#' @export
#' @rdname convert_blanks_to_na
convert_blanks_to_na.character <- function(x) {
  do.call(structure, c(list(if_else(x == "", NA_character_, x)), attributes(x)))
}

#' @export
#' @rdname convert_blanks_to_na
convert_blanks_to_na.list <- function(x) {
  lapply(x, convert_blanks_to_na)
}

#' @export
#' @rdname convert_blanks_to_na
convert_blanks_to_na.data.frame <- function(x) { # nolint
  x[] <- lapply(x, convert_blanks_to_na)
  x
}

#' Get One to Many Values that Led to a Prior Error
#'
#' @export
#'
#' @author Stefan Bundfuss
#'
#' @details
#' If `assert_one_to_one()` detects an issue, the one to many values are stored
#' in a dataset. This dataset can be retrieved by `get_one_to_many_dataset()`.
#'
#' Note that the function always returns the one to many values from the last
#' error that has been thrown in the current R session. Thus, after restarting
#' the R sessions `get_one_to_many_dataset()` will return `NULL` and after a
#' second error has been thrown, the dataset of the first error can no longer be
#' accessed (unless it has been saved in a variable).
#'
#' @return A `data.frame` or `NULL`
#'
#' @keywords user_utility
#'
#' @examples
#' data(admiral_adsl)
#'
#' try(
#'   assert_one_to_one(admiral_adsl, vars(STUDYID), vars(SITEID))
#' )
#'
#' get_one_to_many_dataset()
get_one_to_many_dataset <- function() {
  .datasets$one_to_many
}

#' Get Many to One Values that Led to a Prior Error
#'
#' @export
#'
#' @author Stefan Bundfuss
#'
#' @details
#' If `assert_one_to_one()` detects an issue, the many to one values are stored
#' in a dataset. This dataset can be retrieved by `get_many_to_one_dataset()`.
#'
#' Note that the function always returns the many to one values from the last
#' error that has been thrown in the current R session. Thus, after restarting
#' the R sessions `get_many_to_one_dataset()` will return `NULL` and after a
#' second error has been thrown, the dataset of the first error can no longer be
#' accessed (unless it has been saved in a variable).
#'
#' @return A `data.frame` or `NULL`
#'
#' @keywords user_utility
#'
#' @examples
#' data(admiral_adsl)
#'
#' try(
#'   assert_one_to_one(admiral_adsl, vars(SITEID), vars(STUDYID))
#' )
#'
#' get_many_to_one_dataset()
get_many_to_one_dataset <- function() {
  .datasets$many_to_one
}
