filter_if <- function (dataset, filter) {
  assert_data_frame(dataset, check_is_grouped = FALSE)
  assert_filter_cond(filter, optional = TRUE)
  if (quo_is_null(filter)) {
    dataset
  } else {
    filter(dataset, !!filter)
  }
}
