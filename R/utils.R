enumerate <- function(x, quote_fun = backquote) {
  paste(
    paste0(quote_fun(x[-length(x)]), collapse = ", "),
    "and",
    quote_fun(x[length(x)])
  )
}

squote <- function(x) {
  paste0("'", x, "'")
}

backquote <- function(x) {
  paste0("`", x, "`")
}

toString.tbl_df <- function(x, ...) {
  paste(capture.output(print(x))[-c(1L, 3L)], collapse = "\n")
}
