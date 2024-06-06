#' Process `set_values_to` Argument
#'
#' The function creates the variables specified by the `set_values_to` argument,
#' catches errors, provides user friendly error messages, and optionally checks
#' the type of the created variables.
#'
#' @param dataset Input dataset
#'
#' @param set_values_to Variables to set
#'
#'   A named list returned by `exprs()` defining the variables to be set, e.g.
#'   `exprs(PARAMCD = "OS", PARAM = "Overall Survival")` is expected. The values
#'   must be symbols, character strings, numeric values, expressions, or `NA`.
#'
#' @param expected_types
#'
#'   If the argument is specified, the specified variables are checked whether
#'   the specified type matches the type of the variables created by
#'   `set_values_to`.
#'
#'   *Permitted Values*: A character vector with values `"numeric"` or
#'   `"character"`
#'
#'
#' @return The input dataset with the variables specified by `set_values_to`
#'   added/updated
#'
#' @family utils_help
#' @keywords utils_help
#'
#' @export
#'
#' @examples
#' library(dplyr)
#' data <- tribble(
#'   ~AVAL,
#'   20
#' )
#'
#' try(
#'   process_set_values_to(
#'     data,
#'     set_values_to = exprs(
#'       PARAMCD = BMI
#'     )
#'   )
#' )
#'
#' try(
#'   process_set_values_to(
#'     data,
#'     set_values_to = exprs(
#'       PARAMCD = 42
#'     ),
#'     expected_types = c(PARAMCD = "character")
#'   )
#' )
process_set_values_to <- function(dataset,
                                  set_values_to = NULL,
                                  expected_types = NULL) {
  assert_data_frame(dataset)
  assert_varval_list(set_values_to, optional = TRUE)
  assert_character_vector(
    expected_types,
    values = c("numeric", "character"),
    named = TRUE,
    optional = TRUE
  )

  tryCatch(
    result <- mutate(dataset, !!!set_values_to),
    error = function(cnd) {
      cli_abort(
        message =
          c("Assigning variables failed!",
            "*" = "{.code set_values_to = exprs({paste(names(set_values_to), '=',
                    set_values_to, collapse = ', ')})}",
            "See error message below:",
            conditionMessage(cnd)
          ),
        call = parent.frame(n = 4)
      )
    }
  )
  if (!is.null(expected_types)) {
    types <- map_chr(result, typeof) %>%
      map_chr(function(x) if_else(x %in% c("integer", "double"), "numeric", x))
    vars_to_check <- intersect(names(set_values_to), names(expected_types))
    if (length(vars_to_check) > 0) {
      actual <- types[vars_to_check]
      expected <- expected_types[vars_to_check]
      unexpected <- actual != expected
      if (any(unexpected)) {
        cli_abort(
          message =
            "The following variables have an unexpected type:" |>
              c(stats::setNames(
                paste0(
                  names(actual[unexpected]),
                  ": expected is {.cls ",
                  expected[unexpected],
                  "}, but it is {.cls ",
                  actual[unexpected],
                  "}."
                ),
                nm = rlang::rep_along(actual[unexpected], "*")
              ))
        )
      }
    }
  }
  result
}
