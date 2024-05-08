#' Compute Scale Parameters
#'
#' Computes the average of a set of source values and transforms the result
#' from the source range to the target range. For example, for calculating the
#' average of a set of questionnaire response scores and re-coding the average
#' response to obtain a subscale score.
#'
#' @param source A vector of values to be scaled
#'
#'   A numeric vector is expected.
#'
#' @param source_range The permitted source range
#'
#'   A numeric vector containing two elements is expected, representing the
#'   lower and upper bounds of the permitted source range. Alternatively, if no
#'   argument is specified for `source_range` and `target_range`, no
#'   transformation will be performed.
#'
#' @param target_range The target range
#'
#'   A numeric vector containing two elements is expected, representing the
#'   lower and upper bounds of the target range. Alternatively, if no
#'   argument is specified for `source_range` and `target_range`, no
#'   transformation will be performed.
#'
#' @param flip_direction Flip direction of the scale?
#'
#'   The transformed values will be reversed within the target range, e.g.
#'   within the range 0 to 100, 25 would be reversed to 75.
#'
#'   This argument will be ignored if `source_range` and `target_range` aren't
#'   specified.
#'
#'   Default: `FALSE`
#'
#'   Permitted Values: `TRUE`, `FALSE`
#'
#' @param min_n Minimum number of values for computation
#'
#'   The minimum number of non-missing values in source for the computation to
#'   be carried out. If the number of non-missing values is below `min_n`,
#'   the result will be set to missing, i.e. `NA`.
#'
#'   A positive integer is expected.
#'
#'   Default: 1
#'
#' @details Returns a numeric value. If source contains less than `min_n` values,
#'   the result is set to `NA`. If `source_range` and `target_range` aren't
#'   specified, the mean will be computed without any transformation being
#'   performed.
#'
#' @return The average of source transformed to the target range or `NA` if
#'   source doesn't contain `min_n` values.
#'
#' @keywords com_bds_findings
#'
#' @family com_bds_findings
#'
#' @export
#'
#' @examples
#' compute_scale(
#'   source = c(1, 4, 3, 5),
#'   source_range = c(1, 5),
#'   target_range = c(0, 100),
#'   flip_direction = TRUE,
#'   min_n = 3
#' )
#'
compute_scale <- function(source,
                          source_range = NULL,
                          target_range = NULL,
                          flip_direction = FALSE,
                          min_n = 1) {
  # Function argument checks
  assert_numeric_vector(source) # nolint: undesirable_function_linter
  assert_numeric_vector(source_range, optional = TRUE)
  if (!is.null(target_range) && is.null(source_range)) {
    cli_abort(
      c("Argument {.arg source_range} is missing with no default
         and {.arg target_range} is not missing.",
        "i" = "Either both or neither arguments should be specified."
      )
    )
  }
  assert_numeric_vector(target_range, optional = TRUE)
  if (!is.null(source_range) && is.null(target_range)) {
    cli_abort(
      c("Argument {.arg target_range} is missing with no default
         and {.arg source_range} is not missing.",
        "i" = "Either both or neither arguments should be specified."
      )
    )
  }
  assert_logical_scalar(flip_direction)
  assert_integer_scalar(min_n, subset = "positive")

  # Computation
  if (sum(!is.na(source)) >= min_n) { # nolint: undesirable_function_linter
    target <- mean(source, na.rm = TRUE) # nolint: undesirable_function_linter

    if (!is.null(source_range) && !is.null(target_range)) {
      scale_constant <- min(target_range) - min(source_range)
      scale_coefficient <- (max(target_range) - min(target_range)) /
        (max(source_range) - min(source_range))

      target <- (target + scale_constant) * scale_coefficient

      if (flip_direction == TRUE) {
        target <- max(target_range) - target
      }
    }
  } else {
    target <- NA
  }

  target
}
