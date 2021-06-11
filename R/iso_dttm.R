as_iso_dttm <- function(x, ...) {
  structure(x, class = union("iso_dttm", class(x)))
}

as.POSIXct.iso_dttm <- function(x, ...) {
  structure(x, class = setdiff(class(x), "iso_dttm"))
}

print.iso_dttm <- function(x, ...) {
  print(format(x, "%Y-%m-%dT%H:%M:%S"), quote = FALSE, print.gap = 2L)
}

setOldClass("iso_dttm")

#' Add Methods
#'
#' @param e1 First object
#' @param e2 Second object
#'
#' @rdname add
setMethod(
  "+",
  signature(e1 = "iso_dttm", e2 = "Period"),
  function(e1, e2) as_iso_dttm(as.POSIXct(e1) + e2)
)

#' @rdname add
setMethod(
  "+",
  signature(e1 = "Period", e2 = "iso_dttm"),
  function(e1, e2) as_iso_dttm(e1 + as.POSIXct(e2))
)

#' @rdname add
setMethod(
  "+",
  signature(e1 = "iso_dttm", e2 = "Duration"),
  function(e1, e2) as_iso_dttm(as.POSIXct(e1) + e2)
)

#' @rdname add
setMethod(
  "+",
  signature(e1 = "Duration", e2 = "iso_dttm"),
  function(e1, e2) as_iso_dttm(e1 + as.POSIXct(e2))
)
