# Sometimes functions from other libraries become critically important while using admiral
# and thus should be included as an export. This applies especially to functions which are
# frequently called within `admiral` function arguments. The goal of these exports is such that
# admiral comes ready "out of the box", similar to how one might think the pipe operator, `%>%`,
# comes from `dplyr` but is actually native to `magrittr`.

#' rlang exprs
#'
#' See \code{rlang::\link[rlang:exprs]{exprs}} for details.
#'
#' @name exprs
#' @rdname reexport-exprs
#' @keywords reexport
#' @importFrom rlang exprs
#' @export
NULL

#' dplyr desc
#'
#' See \code{dplyr::\link[dplyr:desc]{desc}} for details.
#'
#' @name desc
#' @rdname reexport-desc
#' @keywords reexport
#' @importFrom dplyr desc
#' @export
NULL

#' Pipe operator
#'
#' See \code{magrittr::\link[magrittr:%>%]{%>%}} for more details.
#'
#' @name %>%
#' @rdname pipe
#' @param lhs A value or the magrittr placeholder.
#' @param rhs A function call using the magrittr semantics.
#' @keywords reexport
#' @export
#' @importFrom magrittr %>%
#' @usage lhs \%>\% rhs
NULL
