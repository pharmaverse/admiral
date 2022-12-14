# Sometimes functions from other libraries become critically important while using admiral
# and thus should be included as an export. This applies especially to functions which are
# frequently called within `admiral` function arguments. The goal of these exports is such that
# admiral comes ready "out of the box", similar to how one might think the pipe operator, `%>%`,
# comes from `dplyr` but is actually native to `magrittr`.

#' @export
dplyr::vars

#' @export
dplyr::desc

#' @export
magrittr::`%>%`
