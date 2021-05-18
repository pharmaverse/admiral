#' Derive Death Cause and Source Domain
#'
#' Derives death cause (`DTHCAUS`) and source domain (`DTHDOM`).
#'
#' @param dataset Input dataset. `USUBJID` is an expected column.
#'
#' @param sources A list of sources.
#'
#' @param by_vars Grouping variables.
#'
#'   Default: `exprs(USUBJID)`
#'
#' @keywords Roche adsl
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
derive_dthcaus <- function(dataset, sources, by_vars = exprs(USUBJID)) {

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
      add_data[[ii]] <- filter_extreme(add_data[[ii]],
                                       order = sources[[ii]]$order,
                                       by_vars = by_vars,
                                       mode = sources[[ii]]$mode)
    }
    add_data[[ii]] <- transmute(add_data[[ii]],
                                temp_source_nr = ii,
                                DTHDOM = sources[[ii]]$dthdom,
                                !!!by_vars,
                                DTHCAUS := !!sources[[ii]]$dthcaus)
  }

  dataset_add <- bind_rows(add_data) %>%
    filter_extreme(by_vars = by_vars,
                   order = exprs(temp_source_nr),
                   mode = "first") %>%
    select(-starts_with("temp_"))

  left_join(dataset, dataset_add, by = map_chr(by_vars, as_string))
}

# derive_dthcaus <- function(dataset, dataset_ae, dataset_ds,
#                            filter_ae = exprs(AEOUT == "FATAL"),
#                            filter_ds = exprs(DSDECOD == "DEATH", grepl("DEATH DUE TO", DSTERM)),
#                            dthcaus_ae = expr(AEDECOD),
#                            dthcaus_ds = expr(DSTERM)) {
#
#   # extract from AE
#   dataset_ae <- filter_extreme(dataset_ae %>%
#                                  filter(!!!filter_ae),
#                                order = exprs(AEDTHDTC),
#                                by_vars = exprs(USUBJID),
#                                mode = "first")
#
#   if (typeof(dthcaus_ae) == "symbol") {
#     dataset_ae <- dataset_ae %>%
#       transmute(USUBJID, DTHDOM = "AE", DTHCAUS = !!dthcaus_ae, temp_source_nr = 1)
#   } else if (typeof(dthcaus_ae) == "character") {
#     dataset_ae <- dataset_ae %>%
#       transmute(USUBJID, DTHDOM = "AE", DTHCAUS = dthcaus_ae, temp_source_nr = 1)
#   } else {
#     stop("Type of `dthcaus_ae` must be either symbol or character.")
#   }
#
#   # extract from DS
#   dataset_ds <- filter_extreme(dataset_ds %>%
#                                  filter(!!!filter_ds),
#                                order = exprs(DSSTDTC),
#                                by_vars = exprs(USUBJID),
#                                mode = "first")
#
#   if (typeof(dthcaus_ds) == "symbol") {
#     dataset_ds <- dataset_ds %>%
#       transmute(USUBJID, DTHDOM = "DS", DTHCAUS = !!dthcaus_ds, temp_source_nr = 2)
#   } else if (typeof(dthcaus_ds) == "character") {
#     dataset_ds <- dataset_ds %>%
#       transmute(USUBJID, DTHDOM = "DS", DTHCAUS = dthcaus_ds, temp_source_nr = 2)
#   } else {
#     stop("Type of `dthcaus_ds` must be either symbol or character.")
#   }
#
#   # combine together and take AE one if exists
#   dataset_add <- bind_rows(dataset_ae, dataset_ds) %>%
#     filter_extreme(by_vars = exprs(USUBJID),
#                    order = exprs(temp_source_nr),
#                    mode = "first") %>%
#     select(-starts_with("temp_"))
#
#   left_join(dataset, dataset_add, by = "USUBJID")
# }
