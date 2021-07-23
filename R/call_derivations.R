call_derivation <- function(dataset, derivation, variable_params, ...) {
  assert_data_frame(dataset)
  # assert_function(derivation)
  # assert_named_list(variable_params)
  fixed_params <- list(...)

  for (i in seq_along(variable_params)) {
    args <- c(quote(dataset), variable_params[[i]], fixed_params)
    dataset <- do.call(derivation, args)
  }
  dataset
}
