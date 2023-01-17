#' @export
process_set_values_to <- function(dataset,
                                  set_values_to,
                                  expected_types = NULL) {
  assert_data_frame(dataset)
  assert_varval_list(set_values_to, accept_expr = TRUE)
  assert_list_of(expected_types, class = "character", optional = TRUE)
  if (!is.null(expected_types)) {
    invalids <- expected_types[!expected_types %in% c("characeter", "numeric")]
    if (length(invalids) > 0) {
      abort(paste(
        "The right-hand side values of `expected_types` must be either \"character\" or \"numeric\".\n",
        "The following elements are invalid:\n",
        paste(names(invalids), invalids, sep = ": ", collapse = "\n")
      ))
    }
  }
  tryCatch(
    mutate(dataset, !!!set_values_to),
    error = function(cnd) {
      abort(
        paste0(
          "Assigning variables failed!\n",
          "set_values_to = (\n",
          paste(
            " ",
            names(set_values_to),
            "=",
            lapply(set_values_to, quo_get_expr),
            collapse = "\n"
          ),
          "\n)\nError message:\n  ",
          cnd
        )
      )
    }
  )
}
