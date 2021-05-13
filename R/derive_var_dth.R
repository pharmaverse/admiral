#' Derive Death Cause, Domain, and Cause of Death Category
#'
#' Derives death cause (`DTHCAUS`), domain (`DTHDOM`), and cause of death category (`DTHCAT`)
#'
#' @param dataset Input dataset
#'
#'   The columns "DTHDOM" and "DTHCAUS" are expected.
#'
#' @details The cause of death category is derived based on DTHDOM and DTHCAUS
#'   values.
#'
#' @keywords Roche
#'
#' @author Shimeng Huang
#'
#' @return The input dataset with `DTHCAT` added
#'
#' @export
#'
#' @examples ...
#'
derive_dth <- function(dataset, dataset_ae, dataset_ds,
                       filter_ae = expr(AEOUT == "FATAL"),
                       filter_ds = expr(DSDECOD == "DEATH")) {
  dthae <- list(dataset = ae,
                vars = exprs(DTHDOM = "AE", DTHCAUS = AEDECOD),
                filter = filter_ae)
  dthds <- list(dataset = ds,
                vars = exprs(DTHDOM = "DS", DTHCAUS = DSTERM),
                filter = filter_ds)
  dataset <- derive_merged_multisource(dataset = dataset,
                                       sources = list(dthae, dthds),
                                       by_vars = exprs(USUBJID),
                                       order = exprs(temp_source_nr),
                                       mode = "first") %>%
    filter(!is.null(DTHDOM)) %>%
    select(USUBJID, DTHCAUS, DTHDOM)

  dataset <- dataset %>%
    mutate(DTHCAT = case_when(
      DTHDOM == "" ~ "",
      DTHDOM == "ADVERSE EVENT" ~ "AE",
      stringr::str_detect(DTHCAUS, "PROGRESSIVE DISEASE") |
        stringr::str_detect(DTHCAUS, "DISEASE RELAPSE") ~ "PROGRESSIVE DISEASE",
      TRUE ~ "OTHER"
    ))

  dataset
}
