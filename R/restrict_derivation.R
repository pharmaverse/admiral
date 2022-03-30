#' Execute a Derivation on a Subset of the Input Dataset
#'
#' Execute a derivation on a subset of the input dataset.
#'
#' @param dataset Input dataset
#'
#' @param derivation Derivation
#'
#' @param args Arguments of the derivation
#'
#'   A `param()` object is expected.
#'
#' @param filter Filter condition
#'
#' @author Stefan Bundfuss
#'
#' @seealso [param()]
#'
#' @export
restrict_derivation <- function(dataset,
                                derivation,
                                args = NULL,
                                filter) {
  assert_data_frame(dataset)
  assert_function(derivation, params = c("dataset"))
  assert_s3_class(args, "params", optional = TRUE)
  if (!is.null(args)) {
    assert_function_param(deparse(substitute(derivation)), names(args))
  }
  filter <- assert_filter_cond(enquo(filter))

  data_ignore <- dataset %>%
    filter(!(!!filter) | is.na(!!filter))
  data <- dataset %>%
    filter(!!filter)

  call <- as.call(c(substitute(derivation), c(quote(data), args)))
  data <- eval(call, envir = list(data = data), enclos = parent.frame())

  bind_rows(data, data_ignore)
}
