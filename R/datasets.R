.datasets <- new.env(parent = emptyenv())

#' @export
set_dataset <- function(dataset, name) {
  assert_data_frame(dataset, check_is_grouped = FALSE)
  assert_character_scalar(name)

  .datasets[[name]] <- dataset
}

#' @export
get_dataset <- function(name) {
  assert_character_scalar(name)

  .datasets[[name]]
}
