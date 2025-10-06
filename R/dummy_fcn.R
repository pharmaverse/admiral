#' Dummy Function Example
#'
#' @description
#' This dummy function demonstrates how to document an R function using roxygen2.
#' It takes two numeric inputs and returns their element-wise sum. Make it pretty!!
#' Making it more pretty!! TEST!! More TESTs!!
#'
#' @param a A numeric vector.
#' @param b A numeric vector.
#'
#' @return A numeric vector containing the element-wise sum of \code{a} and \code{b}.
#'
#' @examples
#' # Basic addition with scalars
#' dummy_function(1, 2)
#'
#' # Element-wise addition with vectors
#' dummy_function(c(1, 2, 3), c(4, 5, 6))
#'
#' @export
dummy_function <- function(a, b) {
  # Validate input types
  if (!is.numeric(a) || !is.numeric(b)) {
    stop("Both 'a' and 'b' must be numeric things")
  }

  # Return the element-wise sum
  a + b
}
