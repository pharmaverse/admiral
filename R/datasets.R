.datasets <- new.env(parent = emptyenv())

#' Set a Dataset in the `.datasets` environment
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
#' the `.datasets` environment. It can be retrieved later on using [get_dataset()]
#'
#' @export
set_dataset <- function(dataset, name) {
  assert_data_frame(dataset, check_is_grouped = FALSE)
  assert_character_scalar(name)

  .datasets[[name]] <- dataset
}

#' Retrieve a Dataset from the `.datasets` environment
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

  .datasets[[name]]
}
