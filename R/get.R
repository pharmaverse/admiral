#' Get Constant Variables
#'
#' @param dataset A data frame.
#' @param by_vars By variables
#'   The groups defined by the by variables are considered separately. I.e., if
#'   a variable is constant within each by group, it is returned.
#'
#' @param ignore_vars Variables to ignore
#'   The specified variables are not considered, i.e., they are not returned
#'   even if they are constant (unless they are included in the by variables).
#'
#'   *Permitted Values:* A list of variable names or selector function calls
#'   like `starts_with("EX")`
#'
#' @keywords get
#' @family get
#'
#'
#' @return Variable vector.
#' @export
get_constant_vars <- function(dataset, by_vars, ignore_vars = NULL) {
  assert_data_frame(dataset, optional = FALSE)
  assert_vars(by_vars, optional = FALSE)
  assert_vars(ignore_vars, optional = TRUE)

  non_by_vars <- setdiff(names(dataset), vars2chr(by_vars))

  if (!is.null(ignore_vars)) {
    non_by_vars <- setdiff(
      non_by_vars,
      vars_select(non_by_vars, !!!ignore_vars)
    )
  }

  # get unique values within each group by variables
  unique_count <- dataset %>%
    group_by(!!!by_vars) %>%
    summarise(across(!!non_by_vars, n_distinct)) %>%
    ungroup() %>%
    select(!!!syms(non_by_vars))

  # determine variables which are constant within each by group
  constant_vars <- unique_count %>%
    map_lgl(~ all(.x == 1)) %>%
    which() %>%
    names() %>%
    syms()

  exprs(!!!by_vars, !!!constant_vars)
}


#' Get Duplicates From a Vector
#'
#' @param x An atomic vector
#'
#' @return A vector of the same type as `x` contain duplicate values
#' @export
#'
#' @keywords get
#' @family get
#'
#' @export
#'
#' @examples
#' get_duplicates(1:10)
#'
#' get_duplicates(c("a", "a", "b", "c", "d", "d"))
get_duplicates <- function(x) {
  assert_atomic_vector(x)

  unique(x[duplicated(x)])
}

#' Get Source Variables from a List of Expressions
#'
#' @param expressions A list of expressions
#'
#' @param quosures *Deprecated*, please use `expressions` instead.
#'
#'
#' @keywords get
#' @family get
#'
#' @return A list of expressions
#' @export
get_source_vars <- function(expressions, quosures) {
  if (!missing(quosures)) {
    deprecate_stop(
      "0.10.0",
      "get_source_vars(quosures = )",
      "get_source_vars(expressions = )"
    )
  }
  assert_varval_list(expressions, optional = TRUE)

  source_vars <- expr_c(expressions)[lapply(expr_c(expressions), is.symbol) == TRUE]

  if (length(source_vars) == 0) {
    NULL
  } else {
    source_vars
  }
}
#' Retrieve a Dataset from the `admiraldev_environment` environment
#'
#' @details
#'
#' Sometimes, developers may want to provide information to users which does not fit into a
#' warning or error message. For example, if the input dataset of a function contains unexpected
#' records, these can be stored in a separate dataset, which users can access to investigate
#' the issue.
#'
#' To achieve this, R has a data structure known as an 'environment'. These environment objects
#' are created at build time, but can be populated with values after the package has been loaded
#' and update those values over the course of an R session.
#'
#' As so, the establishment of `admiraldev_environment` allows us to create dynamic data/objects
#' based on user-inputs that need modification. The purpose of `get_dataset` is to
#' retrieve the datasets contained inside `admiraldev_environment`.
#'
#' Currently we only support two datasets inside our `admiraldev_environment` object:
#'  - `one_to_many`
#'  - `many_to_one`
#'
#'
#' @param name The name of the dataset to retrieve
#'
#' @return A `data.frame`
#'
#'
#' @keywords get
#' @family get
#'
#' @export
get_dataset <- function(name) {
  assert_character_scalar(name, values = c("one_to_many", "many_to_one"))

  admiraldev_environment[[name]]
}
