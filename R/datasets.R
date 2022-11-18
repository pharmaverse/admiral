#' Set a Dataset in the `admiral_environment` environment
#'
#' @param dataset A `data.frame`
#' @param name A name for `dataset`
#'
#' @return
#' No return value, called for side effects
#'
#' @author Thomas Neitmann
#'
#' @keywords datasets
#' @family datasets
#'
#' @details
#' The object passed to the `dataset` argument will be assigned to `name` in
#' the `admiral_environment` environment. It can be retrieved later on using [get_dataset()]
#'
#' @export
set_dataset <- function(dataset, name) {
  assert_data_frame(dataset, check_is_grouped = FALSE)
  assert_character_scalar(name)

  admiral_environment[[name]] <- dataset
}

#' Retrieve a Dataset from the `admiral_environment` environment
#'
#' @param name The name of the dataset to retrieve
#'
#' @return A `data.frame`
#'
#' @author Thomas Neitmann
#'
#' @keywords datasets
#' @family datasets
#'
#' @export
get_dataset <- function(name) {
  assert_character_scalar(name)

  admiral_environment[[name]]
}
