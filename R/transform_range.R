#' Transform Range
#'
#' Transforms results from the source range to the target range. For example,
#' for transforming source values 1, 2, 3, 4, 5 to 0, 25, 50, 75, 100.
#'
#' @param source A vector of values to be transformed
#'
#'   A numeric vector is expected.
#'
#' @param source_range The permitted source range
#'
#'   A numeric vector containing two elements is expected, representing the
#'   lower and upper bounds of the permitted source range.
#'
#' @param target_range The target range
#'
#'   A numeric vector containing two elements is expected, representing the
#'   lower and upper bounds of the target range.
#'
#' @param flip_direction Flip direction of the range?
#'
#'   The transformed values will be reversed within the target range, e.g.
#'   within the range 0 to 100, 25 would be reversed to 75.
#'
#'   *Permitted Values*: `TRUE`, `FALSE`
#'
#' @param outside_range Handling of values outside the source range
#'
#'   Values outside the source range (`source_range`) are transformed to `NA`.
#'
#'   If `"warning"` or `"error"` is specified, a warning or error is issued if
#'   `source` includes any values outside the source range.
#'
#'   *Permitted Values*: `"NA"`, `"warning"`, `"error"`
#'
#' @details Returns the values of `source` linearly transformed from the source
#'   range (`source_range`) to the target range (`target_range`). Values outside
#'   the source range are set to `NA`.
#'
#' @return The source linearly transformed to the target range
#'
#' @keywords com_bds_findings
#'
#' @family com_bds_findings
#'
#' @export
#'
#' @examples
#' transform_range(
#'   source = c(1, 4, 3, 6, 5),
#'   source_range = c(1, 5),
#'   target_range = c(0, 100)
#' )
#'
#' transform_range(
#'   source = c(1, 4, 3, 6, 5),
#'   source_range = c(1, 5),
#'   target_range = c(0, 100),
#'   flip_direction = TRUE
#' )
transform_range <- function(source,
                            source_range,
                            target_range,
                            flip_direction = FALSE,
                            outside_range = "NA") {
  # Function argument checks
  assert_numeric_vector(source)
  assert_numeric_vector(source_range, length = 2)
  assert_numeric_vector(target_range, length = 2)
  assert_logical_scalar(flip_direction)
  assert_character_scalar(outside_range, values = c("NA", "error", "warning"))

  outsider <- !(between(source, source_range[[1]], source_range[[2]]) | is.na(source))
  if (any(outsider)) {
    outside_index <- which(outsider)
    outside_value <- source[outsider]
    source <- if_else(outsider, NA, source)
    msg <- c(
      paste(
        "{.arg source} contains values outside the range of {.val {source_range[[1]]}}",
        "to {.val {source_range[[2]]}}:"
      ),
      paste0("source[[", outside_index, "]] = {.val {", outside_value, "}}")
    )
    if (outside_range == "warning") {
      cli_warn(
        msg,
        class = c("outside_source_range", "assert-admiral"),
        outside_index = outside_index,
        outside_value = outside_value
      )
    } else if (outside_range == "error") {
      cli_abort(
        msg,
        class = c("outside_source_range", "assert-admiral"),
        outside_index = outside_index,
        outside_value = outside_value
      )
    }
  }

  # Computation
  range_constant <- min(target_range) - min(source_range)
  range_coefficient <- (max(target_range) - min(target_range)) /
    (max(source_range) - min(source_range))

  target <- (source + range_constant) * range_coefficient

  if (flip_direction == TRUE) {
    target <- max(target_range) - target
  }

  target
}
