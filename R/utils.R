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

toString.tbl_df <- function(x, ...) {
  paste(capture.output(print(x))[-c(1L, 3L)], collapse = "\n")
}
