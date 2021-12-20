auto_name_list <- function(list) {
  if (is_named(list)) {
    return(list)
  }
  elements <- as.list(substitute(list)[-1L])
  names(list) <- vapply(elements, deparse, "")
  list
}
