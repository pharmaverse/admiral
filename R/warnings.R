#' Warn If a Variable Already Exists
#'
#' Warn if a variable already exists inside a dataset
#'
#' @param dataset A `data.frame`
#' @param vars `character` vector of columns to check for in `dataset`
#'
#' @author Thomas Neitmann
#'
#' @export
#'
#' @examples
#' data(dm)
#'
#' ## No warning as `AAGE` doesn't exist in `dm`
#' warn_if_vars_exist(dm, "AAGE")
#'
#' ## Issues a warning
#' warn_if_vars_exist(dm, "ARM")
#'
warn_if_vars_exist <- function(dataset, vars) {
  existing_vars <- vars[vars %in% colnames(dataset)]
  if (length(existing_vars) == 1L) {
    msg <- paste("Variable", backquote(existing_vars), "already exists in the dataset")
    warn(msg)
  } else if (length(existing_vars) > 1L) {
    msg <- paste("Variables", enumerate(existing_vars), "already exist in the dataset")
    warn(msg)
  } else {
    invisible(NULL)
  }
}
