#' Retrieve a Dataset from the `admiraldev_environment` environment
#'
#' @param name The name of the dataset to retrieve
#'
#' @return A `data.frame`
#'
#'
#' @keywords datasets
#' @family datasets
#'
#' @export
get_dataset <- function(name) {
  assert_character_scalar(name, values = c("one_to_many", "many_to_one"))

  admiraldev_environment[[name]]
}
