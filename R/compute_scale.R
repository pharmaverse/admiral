#' Compute scale parameters
#'
#' Computes the average of source and transforms the result from the source
#' range to the target range.
#'
#' @param source A vector of values to be scaled
#'
#'   A numeric vector is expected.
#'
#' @param source_range The permitted source range
#'
#'   A numeric vector containing two elements is expected, representing the
#'   lower and upper bounds of the source range.
#'
#' @param target_range The permitted target range
#'
#'   A numeric vector containing two elements is expected, representing the
#'   lower and upper bounds of the target range.
#'
#' @param flip_direction Flip direction of the scale?
#'
#'   The transformed values will be flipped within the target range.
#'   REMOVE: Is there a better word for flip in this context?
#'
#'   Default: `FALSE`
#'
#'   Permitted Values: `TRUE` `FALSE`
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
#' @details
#' Returns a numeric value. If source contains less than min_n values, the
#' result is set to NA.
#'
#' @return A value representing the average of source transformed to the target
#'   range or `NA` if source doesn't contain `min_n` values.
#'
#' @keywords com_bds_findings
#'
#' @family com_bds_findings
#'
#' @export
#'
#' @examples
#'
compute_scale <- function(source,
                          source_range,
                          target_range,
                          flip_direction = FALSE,
                          min_n = 1){

  # Function argument checks
  assert_numeric_vector(source)
  assert_numeric_vector(source_range)
  assert_numeric_vector(target_range)
  assert_logical_scalar(flip_direction)
  assert_integer_scalar(min_n, subset = "positive")


  # Computation
  if(sum(!is.na(source)) >= min_n){
    source_mean <- mean(source, na.rm = TRUE)

    scale_constant <- min(target_range) - min(source_range)
    scale_coefficient <- (max(target_range) - min(target_range)) / (max(source_range) - min(source_range))

    target = (source_mean + scale_constant) * scale_coefficient

    if(flip_direction == TRUE){
      target = max(target_range) - target
    }
  }
  else{
    target = NA
  }

  target

}
