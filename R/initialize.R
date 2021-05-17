#' Derive Predecessors Based Upon Specification
#'
#' Derive all predecessor variables of an ADaM dataset based upon metadata
#'
#' @param dataset_name The name of the ADaM dataset to create, e.g. `"ADVS"`
#' @param metadata The result of [read_dap_m3()]
#' @param source_datasets A `list` of source datasets
#'
#' @author Thomas Neitmann
#'
#' @return
#' A dataset containing all predecessor variables as listed in `metadata`
#'
#' @export
#'
#' @examples
#' \dontrun{
#' data("adsl")
#' data("vs")
#' meta <- read_dap_m3("DAP_M3.xlsx")
#' initialize("ADVS", meta, list(adsl, vs))
#' }
initialize <- function(dataset_name, metadata, source_datasets) {
  assert_that(
    is.character(dataset_name),
    inherits(metadata, "DAP_M3"),
    is.list(source_datasets)
  )

  if (dataset_name %!in% names(metadata)) {
    msg <- paste0("Metadata does not contain a dataset named ", squote(dataset_name))
    abort(msg)
  }

  predecessors <- metadata[[dataset_name]] %>%
    rename(Derivation = `Derivation / Comment`) %>%
    mutate(Derivation = str_trim(Derivation)) %>%
    filter(
      Dataset == dataset_name,
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
