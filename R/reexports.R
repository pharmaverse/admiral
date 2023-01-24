# Sometimes functions from other libraries become critically important while using admiral
# and thus should be included as an export. This applies especially to functions which are
# frequently called within `admiral` function arguments. The goal of these exports is such that
# admiral comes ready "out of the box", similar to how one might think the pipe operator, `%>%`,
# comes from `dplyr` but is actually native to `magrittr`.

#' @export
rlang::exprs

#' Create List of Quosures
#'
#' @description
#' `r lifecycle::badge("deprecated")`
#'
#' This function is *deprecated*, please use `exprs()` instead.
#'
#' @param ... List of variables
#'
#' @return List of expressions
#'
#' @keywords deprecated
#' @family deprecated
#'
#' @export
vars <- function(...) {
  deprecate_warn(
    "0.10.0",
    "vars()",
    "exprs()",
    details = paste(
      "The admiral functions no longer expect list of quosures created by `vars()`",
      "but list of expressions created by `exprs()`"
    )
  )
  exprs(...)
}

#' @export
dplyr::desc

#' @export
magrittr::`%>%`
