#' @export
rlang::expr

#'@export
rlang::exprs

enumerate <- function(x) {
  paste(
    paste0(backquote(x[-length(x)]), collapse = ", "),
    "and",
    backquote(x[length(x)])
  )
}

backquote <- function(x) {
  paste0("`", x, "`")
}

`%!in%` <- Negate(`%in%`)
