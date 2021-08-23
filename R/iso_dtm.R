as_iso_dtm <- function(x, time_zone = Sys.timezone()) {
  dtm <- ymd_hms(x, tz = time_zone)
  class(dtm) <- union("iso_dtm", class(dtm))
  dtm
}

as.POSIXct.iso_dtm <- function(x, ...) {
  structure(x, class = setdiff(class(x), "iso_dtm"))
}

print.iso_dtm <- function(x, ...) {
  print(format(x, "%Y-%m-%d %H:%M:%S"), quote = FALSE, print.gap = 2L)
}

setOldClass("iso_dtm")

#' Add Methods
#'
#' @param e1 First object
#' @param e2 Second object
#'
#' @rdname add
setMethod(
  "+",
  signature(e1 = "iso_dtm", e2 = "Period"),
  function(e1, e2) as_iso_dtm(as.POSIXct(e1) + e2)
)

#' @rdname add
setMethod(
  "+",
  signature(e1 = "Period", e2 = "iso_dtm"),
  function(e1, e2) as_iso_dtm(e1 + as.POSIXct(e2))
)

#' @rdname add
setMethod(
  "+",
  signature(e1 = "iso_dtm", e2 = "Duration"),
  function(e1, e2) as_iso_dtm(as.POSIXct(e1) + e2)
)

#' @rdname add
setMethod(
  "+",
  signature(e1 = "Duration", e2 = "iso_dtm"),
  function(e1, e2) as_iso_dtm(e1 + as.POSIXct(e2))
)
