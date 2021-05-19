#' Derive Death Cause and Source Domain
#'
#' Derives death cause (`DTHCAUS`) and source domain (`DTHDOM`).
#'
#' @param dataset Input dataset. `USUBJID` is an expected column.
#'
#' @param sources A list of sources.
#'
#' Each element is a list containing at least
#'   the following fields:
#'   * dataset: the source dataset.
#'   * dthdom: name of the source domain as a string, e.g. "AE".
#'   * dthcaus: a symbol or a string --- if a symbol, e.g., `expr(AEDECOD)`,
#'   it is the variable in the source dataset to be used to assign values to
#'   `DTHCAUS`; if a string, it is the fixed value to be assigned to `DTHCAUS`.
#'
#' Additional fields can also be supplied:
#'   * order: variables to be used to sort the source dataset,
#'   e.g., `exprs(AEDTHDTC)`.
#'   * mode: filtering "first" or "last" when multiple records exist.
#'
#' Note: If a subject has death info from multiple sources, the one from the
#'   first source will be kept.
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
#'   filter = expr(AEOUT == "FATAL"),
#'   order = exprs(AEDTHDTC),
#'   mode = "first",
#'   dthdom = "AE",
#'   dthcaus = expr(AEDECOD)
#' )
#' src_ds <- list(
#'   dataset = ds,
#'   filter = expr(DSDECOD == "DEATH" & grepl("DEATH DUE TO", DSTERM)),
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

    if (!is.character(sources[[ii]]$dthdom)) {
      stop("`dthdom` must be a string.")
    }

    if (typeof(sources[[ii]]$dthcaus) != "symbol" &
        typeof(sources[[ii]]$dthcaus) != "character") {
      stop("`dthcaus` must be a single expression or a string.")
    }

    if ("mode" %in% names(sources[[ii]])) {
      if (sources[[ii]]$mode %!in% c("first", "last")) {
        stop("`mode` must be one of 'first' or 'last'.")
      }
    }
    else if ("mode" %!in% names(sources[[ii]])) {
      # default to "first" if mode is not given
      sources[[ii]]$mode <- "first"
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

    if (typeof(sources[[ii]]$dthcaus) == "symbol") {
      add_data[[ii]] <- add_data[[ii]] %>%
        transmute(USUBJID,
                  temp_source_nr = ii,
                  DTHDOM = sources[[ii]]$dthdom,
                  DTHCAUS := !!sources[[ii]]$dthcaus)
    }
    else {
      add_data[[ii]] <- add_data[[ii]] %>%
        transmute(USUBJID,
                  temp_source_nr = ii,
                  DTHDOM = sources[[ii]]$dthdom,
                  DTHCAUS = sources[[ii]]$dthcaus)
    }
  }

  # if a subject has multiple death info, keep the one from the first source
  dataset_add <- bind_rows(add_data) %>%
    filter_extreme(by_vars = exprs(USUBJID),
                   order = exprs(temp_source_nr),
                   mode = "first") %>%
    select(-starts_with("temp_"))

  left_join(dataset, dataset_add, by = "USUBJID")
}

#' Create an `dthcaus_source` object
#'
#' @param dataset A data.frame containint a source dataset.
#' @param filter A symbol returned by `expr` to be used for filtering `dataset`.
#' @param order Alist returned by `exprs` to be used for sorting `dataset`.
#' @param mode One of "first" or "last".
dthcaus_source <- function(dataset, filter, order, mode = "first",
                           dthdom, dthcaus) {

  out <- list(
    dataset = dataset,
    filter = filter,
    order = order,
    mode = mode,
    dthdom = dthdom,
    dthcaus = dthcaus
  )
  class(out) <- "dthcaus_source"
  validate_dthcaus_source(out)
}

validate_dthcaus_source <- function(x) {
  values <- unclass(x)
  stopifnot("data.frame" %in% class(values$dataset))
  stopifnot(is_expr(values$filter))
  stopifnot(is_unnamed_exprs(values$order))
  match.arg(values$mode, c("first", "last"))
  stopifnot(is.character(values$dthdom))
  stopifnot(!is.symbol(values$dthcaus) | !is.character(values$dthcaus))
  x
}
