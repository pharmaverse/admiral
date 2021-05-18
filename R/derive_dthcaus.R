#' Derive Death Cause and Source Domain
#'
#' Derives death cause (`DTHCAUS`) and source domain (`DTHDOM`).
#'
#' @param dataset Input dataset. `USUBJID` is an expected column.
#'
#' @param dataset_ae Source AE dataset. `AEDTHDTC` is an expected column.
#'
#' @param dataset_ds Source DS dataset. `DSSTDTC` is an expected column.
#'
#' @param filter_ae An expression to filter `dataset_ae`.
#'
#' @param filter_ds An expression to filter `dataset_ds`.
#'
#' @param dthcaus_ae Column name (as symbol) in AE or value (as character) to
#'   be set as `DTHCAUS` from AE.
#'
#' @param dthcaus_ds Column name (as symbol) in DS or value (as character) to
#'   be set as `DTHCAUS` from DS.
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
#' data("adsl")
#' data("ae")
#' data("ds")
#'
#' derive_dthcaus(data, ae, ds)
#'
derive_dthcaus <- function(dataset, dataset_ae, dataset_ds,
                           filter_ae = exprs(AEOUT == "FATAL"),
                           filter_ds = exprs(DSDECOD == "DEATH", grepl("DEATH DUE TO", DSTERM)),
                           dthcaus_ae = expr(AEDECOD),
                           dthcaus_ds = expr(DSTERM)) {

  # extract from AE
  dataset_ae <- filter_extreme(dataset_ae %>%
                                 filter(!!!filter_ae),
                               order = exprs(AEDTHDTC),
                               by_vars = exprs(USUBJID),
                               mode = "first")

  if (typeof(dthcaus_ae) == "symbol") {
    dataset_ae <- dataset_ae %>%
      transmute(USUBJID, DTHDOM = "AE", DTHCAUS = !!dthcaus_ae, temp_source_nr = 1)
  } else if (typeof(dthcaus_ae) == "character") {
    dataset_ae <- dataset_ae %>%
      transmute(USUBJID, DTHDOM = "AE", DTHCAUS = dthcaus_ae, temp_source_nr = 1)
  } else {
    stop("Type of `dthcaus_ae` must be either symbol or character.")
  }

  # extract from DS
  dataset_ds <- filter_extreme(dataset_ds %>%
                                 filter(!!!filter_ds),
                               order = exprs(DSSTDTC),
                               by_vars = exprs(USUBJID),
                               mode = "first")

  if (typeof(dthcaus_ds) == "symbol") {
    dataset_ds <- dataset_ds %>%
      transmute(USUBJID, DTHDOM = "DS", DTHCAUS = !!dthcaus_ds, temp_source_nr = 2)
  } else if (typeof(dthcaus_ds) == "character") {
    dataset_ds <- dataset_ds %>%
      transmute(USUBJID, DTHDOM = "DS", DTHCAUS = dthcaus_ds, temp_source_nr = 2)
  } else {
    stop("Type of `dthcaus_ds` must be either symbol or character.")
  }

  # combine together and take AE one if exists
  dataset_add <- bind_rows(dataset_ae, dataset_ds) %>%
    filter_extreme(by_vars = exprs(USUBJID),
                   order = exprs(temp_source_nr),
                   mode = "first") %>%
    select(-starts_with("temp_"))

  left_join(dataset, dataset_add, by = "USUBJID")
}
