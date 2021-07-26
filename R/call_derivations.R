call_derivation <- function(dataset, derivation, variable_params, ...) {
  assert_data_frame(dataset)
  # assert_function(derivation)
  # assert_named_list(variable_params)
  fixed_params <- eval(substitute(alist(...)))

  # assert_consistent_lists(variable_params)
  # assert_has_params(derivation, c(names(variable_params[[1L]], names(fixed_params))))

  for (i in seq_along(variable_params)) {
    args <- c(quote(dataset), variable_params[[i]], fixed_params)
    call <- as.call(c(substitute(derivation), args))
    dataset <- eval(call)
  }
  dataset
}

params <- function(...) {
  args <- eval(substitute(alist(...)))
  if (any(names(args) == "")) {
    abort("All arguments must be named")
  }
  args
}
