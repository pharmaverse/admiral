.duplicates <- new.env(parent = emptyenv())

get_duplicates_dataset <- function() {
  .duplicates$ds
}

get_duplicate_records <- function(dataset,
                                  by_vars) {
  assert_that(
    is.data.frame(dataset),
    is_vars(by_vars)
  )

  data_by <- dataset %>%
    ungroup() %>%
    select(!!!by_vars)

  is_duplicate <- duplicated(data_by) | duplicated(data_by, fromLast = TRUE)

  dataset %>%
    select(!!!by_vars, dplyr::everything()) %>%
    filter(is_duplicate)
}

assert_has_unique_records <- function(dataset, by_vars, msg) {
  assert_that(
    is.data.frame(dataset),
    is_vars(by_vars),
    is.character(msg)
  )

  duplicate_records <- get_duplicate_records(dataset, by_vars)
  if (nrow(duplicate_records) >= 1L) {
    .duplicates$ds <- structure(
      duplicate_records,
      class = union("duplicates", class(duplicate_records)),
      by_vars = vars2chr(by_vars)
    )
    full_msg <- paste0(msg, "\nRun `get_duplicates_dataset()` to access the duplicate records.")
    abort(full_msg)
  }
}

print.duplicates <- function(x, ...) {
  cat(
    "Dataset contains duplicate records with respect to ",
    enumerate(attr(x, "by_vars")),
    ".\n",
    sep = ""
  )
  NextMethod()
}
