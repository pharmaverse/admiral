.datasets <- new.env(parent = emptyenv())

#' Set a Dataset in the `.datasets` environment
#'
#' @param dataset A `data.frame`
#' @param name A name for `dataset`
#'
#' @author Thomas Neitmann
#'
#' @details
#' The object passed to the `dataset` argument will be assigned to `name` in
#' the `.datasets` environment. It can be retrived later on using [get_dataset()]
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
#' @author Thomas Neitmann
#'
#' @export
get_dataset <- function(name) {
  assert_character_scalar(name)

  .datasets[[name]]
}
