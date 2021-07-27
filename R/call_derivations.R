call_derivation <- function(dataset, derivation, variable_params, ...) {
  assert_data_frame(dataset)
  assert_s3_class(derivation, "function")
  assert_list_of(variable_params, "params")

  fixed_params <- eval(substitute(alist(...)))
  if (is.null(names(fixed_params)) || any(names(fixed_params) == "")) {
    abort("All arguments inside `...` must be named")
  }

  all_params <- base::union(unlist(map(variable_params, names)), names(fixed_params))
  assert_function_param(deparse(substitute(derivation)), all_params)

  for (i in seq_along(variable_params)) {
    args <- c(quote(dataset), variable_params[[i]], fixed_params)
    call <- as.call(c(substitute(derivation), args))
    dataset <- eval(call)
  }

  dataset
}

params <- function(...) {
  args <- eval(substitute(alist(...)))
  params <- names(args)
  if (is.null(params) || any(params == "")) {
    abort("All arguments passed to `params()` must be named")
  }
  structure(args, class = c("params", "list"))
}
