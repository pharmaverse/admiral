#' Get a New Temporary Variable Name for a Dataset
#'
#' @param dataset The input dataset
#' @param prefix The prefix of the new temporary variable name to create
#'
#' @details
#' The function returns a new unique temporary variable name to be used inside
#' `dataset`. The temporary variable names have the structure `prefix_n` where
#' `n` is an integer, e.g. `tmp_var_1`. If there is already a variable inside
#' `datset` with a given `prefix` then the suffix is increased by 1, e.g. if
#' `tmp_var_1` already exists then `get_new_tmp_var()` will return `tmp_var_2`.
#'
#' @return The name of a new temporary variable as a symbol
#'
#' @keywords tmp_vars
#' @seealso [remove_tmp_vars()]
#'
#' @export
#'
#' @examples
#' library(dplyr, warn.conflicts = FALSE)
#' dm <- tribble(
#'   ~DOMAIN,  ~STUDYID,      ~USUBJID,
#'   "DM",    "STUDY X", "01-701-1015",
#'   "DM",    "STUDY X", "01-701-1016",
#' )
#'
#' tmp_var <- get_new_tmp_var(dm)
#' mutate(dm, !!tmp_var := NA)
get_new_tmp_var <- function(dataset, prefix = "tmp_var") {
  assert_data_frame(dataset, optional = TRUE)
  assert_character_scalar(prefix)
  if (!str_detect(prefix, "^tmp_")) {
    cli_abort("{.arg prefix} must start with {.val tmp_}.")
  }

  caller_env <- parent.frame()

  if (is.null(dataset)) {
    colnames_with_prefix <- vector("character")
  } else {
    regexp <- str_c("^", prefix, "_[0-9]{1,}$")
    colnames_with_prefix <- str_subset(colnames(dataset), regexp)
  }
  if (!is.null(caller_env$.tmp_vars)) {
    suffices <- str_extract(caller_env$.tmp_vars, "[0-9]{1,}$")
    counter <- max(as.integer(suffices)) + 1L
  } else if (length(colnames_with_prefix) > 0L) {
    suffices <- str_extract(colnames_with_prefix, "[0-9]{1,}$")
    counter <- max(as.integer(suffices)) + 1L
  } else {
    counter <- 1L
  }

  new_tmp_var <- str_c(prefix, counter, sep = "_")
  caller_env$.tmp_vars <- c(caller_env$.tmp_vars, new_tmp_var)

  sym(new_tmp_var)
}

#' Remove All Temporary Variables Created Within the Current Function Environment
#'
#' @param dataset The input dataset
#'
#' @export
#'
#' @keywords tmp_vars
#' @seealso [get_new_tmp_var()]
#'
#' @return
#' The input dataset with temporary variables removed
#'
#' @examples
#' library(dplyr, warn.conflicts = FALSE)
#' dm <- tribble(
#'   ~DOMAIN,  ~STUDYID,      ~USUBJID,
#'   "DM",    "STUDY X", "01-701-1015",
#'   "DM",    "STUDY X", "01-701-1016",
#' )
#' dm <- select(dm, USUBJID)
#' tmp_var <- get_new_tmp_var(dm)
#' dm <- mutate(dm, !!tmp_var := NA)
#'
#' ## This function creates two new temporary variables which are removed when calling
#' ## `remove_tmp_vars()`. Note that any temporary variable created outside this
#' ## function is **not** removed
#' do_something <- function(dataset) {
#'   tmp_var_1 <- get_new_tmp_var(dm)
#'   tmp_var_2 <- get_new_tmp_var(dm)
#'   dm %>%
#'     mutate(!!tmp_var_1 := NA, !!tmp_var_2 := NA) %>%
#'     print() %>%
#'     remove_tmp_vars()
#' }
#'
#' do_something(dm)
remove_tmp_vars <- function(dataset) {
  # In order to find the "correct" calling environment we have to make sure to
  # exclude all calls relating to the use of `%>%` from the call stack
  calls <- lapply(sys.calls(), function(x) paste(deparse(x), collapse = " "))
  contains_pipe <- map_lgl(calls, str_detect, "%>%")

  if (any(contains_pipe)) {
    last_pipe <- max(which(contains_pipe))
    tmp_vars <- sys.frame(last_pipe - 1L)$.tmp_vars
  } else {
    tmp_vars <- parent.frame()$.tmp_vars
  }

  if (is.null(tmp_vars)) {
    dataset
  } else {
    dataset[, colnames(dataset) %notin% tmp_vars]
  }
}
