initialize <- function(dataset, metadata, source_datasets) {
  assert_that(
    is.character(dataset),
    inherits(metadata, "DAP_M3"),
    is.list(source_datasets)
  )

  if (dataset %!in% names(metadata)) {
    msg <- paste0("Metadata does not contain a dataset named ", squote(dataset))
    abort(msg)
  }

  predecessors <- metadata[[dataset]] %>%
    rename(Derivation = `Derivation / Comment`) %>%
    mutate(Derivation = str_trim(Derivation)) %>%
    filter(
      Dataset == dataset,
      Source == "Predecessor",
      grepl("^[A-Z]{2,8}\\.[A-Z]{1,8}$", Derivation)
    ) %>%
    pull(Derivation) %>%
    str_remove("^[A-Z]{2,8}\\.")

  available_vars <- reduce(map(source_datasets, names), union)
  existing_predecessors <- predecessors[predecessors %in% available_vars]
  missing_predecessors <- setdiff(predecessors, existing_predecessors)

  if (length(missing_predecessors) >= 1L) {
    msg <- paste(
      "The following predecessors are missing from the source datasets:",
      enumerate(missing_predecessors)
    )
    warn(msg)
  }

  source_datasets %>%
    reduce(full_join, by = c("STUDYID", "USUBJID"), suffix = c("", "_")) %>%
    select(!!!syms(existing_predecessors))
}
