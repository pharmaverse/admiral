#' Derive Death Cause and Source Domain
#'
#' Derives death cause (`DTHCAUS`) and source domain (`DTHDOM`).
#'
#' @param dataset Input dataset. `USUBJID` is an expected column.
#'
#' @param sources A list of sources.
#'
#' @keywords adsl
#'
#' @author Shimeng Huang
#'
#' @return The input dataset with `DTHCAUS` and `DTHDOM` added.
#'
#' @export
#'
#' @examples
#' adsl <- tibble::tribble(
#'  ~STUDYID, ~USUBJID,
#'  "STUDY01", "PAT01",
#'  "STUDY01", "PAT02"
#' )
#' ae <- tibble::tribble(
#'   ~STUDYID, ~USUBJID, ~AEDECOD, ~AEOUT, ~AEDTHDTC,
#'   "STUDY01", "PAT01", "SUDDEN DEATH" , "FATAL", "2021-04-04"
#' )
#' ds <- tibble::tribble(
#'   ~STUDYID, ~USUBJID, ~DSDECOD, ~DSTERM, ~DSSTDTC,
#'   "STUDY01", "PAT02", "INFORMED CONSENT OBTAINED", "INFORMED CONSENT OBTAINED", "2021-04-03",
#'   "STUDY01", "PAT02", "RANDOMIZATION", "RANDOMIZATION", "2021-04-11",
#'   "STUDY01", "PAT02", "DEATH", "DEATH DUE TO PROGRESSION OF DISEASE", "2022-02-01"
#' )
#'
#' src_ae <- list(
#'   dataset = ae,
#'   filter = exprs(AEOUT == "FATAL"),
#'   order = exprs(AEDTHDTC),
#'   mode = "first",
#'   dthdom = "AE",
#'   dthcaus = expr(AEDECOD)
#' )
#' src_ds <- list(
#'   dataset = ds,
#'   filter = exprs(DSDECOD == "DEATH", grepl("DEATH DUE TO", DSTERM)),
#'   order = exprs(DSSTDTC),
#'   mode = "first",
#'   dthdom = "DS",
#'   dthcaus = expr(DSTERM)
#' )
#'
#' derive_dthcaus(adsl, list(src_ae, src_ds))
derive_dthcaus <- function(dataset, sources) {

  # check and clean up input
  for (ii in seq_along(sources)) {

    mandatory_fields <- c("dataset", "dthcaus", "dthdom")
    if (!all(mandatory_fields %in% names(sources[[ii]]))) {
      stop(paste0("Each list in `sources` must have certain fields: ",
                  paste0(mandatory_fields, collapse = ", ")))
    }

    # default to "first" if mode is not given
    if ("mode" %!in% names(sources[[ii]])) {
      sources[[ii]]$mode <- "first"
    }

    if (typeof(sources[[ii]]$dthcaus) != "symbol") {
      stop("`dthcaus` must be a single expression.")
    }
  }

  # process each source
  add_data <- vector("list", length(sources))
  for (ii in seq_along(sources)) {

    if (!is.null(sources[[ii]]$filter)) {
      add_data[[ii]] <- sources[[ii]]$dataset %>%
        filter(!!!(sources[[ii]]$filter))
    }
    else {
      add_data[[ii]] <- sources[[ii]]$dataset
    }

    if (!is.null(sources[[ii]]$order)){
      add_data[[ii]] <- add_data[[ii]] %>%
        filter_extreme(order = sources[[ii]]$order,
                       by_vars = exprs(USUBJID),
                       mode = sources[[ii]]$mode)
    }
    add_data[[ii]] <- transmute(add_data[[ii]],
                                temp_source_nr = ii,
                                DTHDOM = sources[[ii]]$dthdom,
                                DTHCAUS := !!sources[[ii]]$dthcaus)
  }

  dataset_add <- bind_rows(add_data) %>%
    filter_extreme(by_vars = exprs(USUBJID),
                   order = exprs(temp_source_nr),
                   mode = "first") %>%
    select(-starts_with("temp_"))

  left_join(dataset, dataset_add, by = "USUBJID")
}
