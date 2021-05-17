#' Derive Death Cause and Domain.
#'
#' Derives death cause (`DTHCAUS`) and domain (`DTHDOM`).
#'
#' @param dataset Input dataset
#'
#' @param dataset_ae Source AE dataset
#'
#' @param dataset_ds Source DS dataset
#'
#' @param filter_ae An expression to filter `dataset_ae`.
#'
#' @param filter_ds An expression to filter `dataset_ds`.
#'
#' @param dthcaus_ae Variable name (as symbol) or value (as character) to
#'   be set as `DTHCAUS` from AE.
#'
#' @param dthcaus_ds Variable name (as symbol) in DS or value (as character) to
#'   be set as `DTHCAUS` from DS.
#'
#' @keywords Roche
#'
#' @author Shimeng Huang
#'
#' @return The input dataset with `DTHCAT` added.
#'
#' @export
#'
#' @examples ...
#'
derive_dthcaus <- function(dataset, dataset_ae, dataset_ds,
                           filter_ae = exprs(AEOUT == "FATAL"),
                           filter_ds = exprs(DSDECOD == "DEATH", grepl("DEATH DUE TO", DSTERM)),
                           dthcaus_ae = expr(AEDECOD),
                           dthcaus_ds = expr(DSTERM)) {

  # extract from AE
  dataset_ae <- filter_extreme(dataset_ae %>%
                                 filter(!!filter_ae),
                               order = exprs(AEDTHDTC),
                               by_vars = exprs(USUBJID),
                               mode = "first") %>%
    mutate(temp_source_nr = 1)

  if (typeof(dthcaus_ae) == "symbol") {
    dataset_ae <- dataset_ae %>%
      transmute(DTHDOM = "AE", DTHCAUS = !!dthcaus_ae)
  } else if (typeof(dthcaus_ae) == "character") {
    dataset_ae <- dataset_ae %>%
      transmute(DTHDOM = "AE", DTHCAUS = dthcaus_ae)
  } else {
    stop("Type of `dthcaus_ae` must be either symbol or character.")
  }

  # extract from DS
  dataset_ds <- filter_extreme(dataset_ae %>%
                                 filter(!!filter_ds),
                               order = exprs(DSSTDTC),
                               by_vars = exprs(USUBJID),
                               mode = "first") %>%
    mutate(temp_source_nr = 2)

  if (typeof(dthcaus_ds) == "symbol") {
    dataset_ds <- dataset_ds %>%
      transmute(DTHDOM = "DS", DTHCAUS = !!dthcaus_ds)
  } else if (typeof(dthcaus_ds) == "character") {
    dataset_ds <- dataset_ds %>%
      transmute(DTHDOM = "DS", DTHCAUS = dthcaus_ds)
  } else {
    stop("Type of `dthcaus_ds` must be either symbol or character.")
  }

  # combine together and take AE one if exists
  dataset_add <- bind_rows(dataset_ae, dataset_ds) %>%
    filter_extreme(by_vars = exprs(USUBJID),
                   order = temp_source_nr,
                   mode = "first") %>%
    select(-starts_with("temp_"))

  left_join(dataset, dataset_add, by = "USUBJID")
}
