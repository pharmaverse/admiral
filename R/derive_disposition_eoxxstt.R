#' Derive a disposition status at a specific timepoint
#'
#' Derive a disposition status from the the relevant records in the disposition domain.
#'
#' @param dataset Input dataset.
#'
#' @param dataset_ds Dataset containing the disposition information (e.g.: ds).
#'
#' The variable specified in the status_var parameter must be in dataset_ds.
#'
#' @param new_var Name of the disposition date variable.
#'
#' a variable name is expected (e.g. EOSSTT).
#'
#' @param status_var The variable used to derive the disposition status.
#'
#'   A character vector is expected (e.g. DSDECOD).
#'
#' @param format_new_var The format used to derive the status.
#'
#' Default: format_eoxxstt_default defined as:
#' format_eoxxstt_default<-function(x){
#' case_when (
#'   x %in% c("COMPLETED")~"COMPLETED",
#'   !(x %in% c("COMPLETED")) & ! is.na(x)~"DISCONTINUED",
#'   TRUE ~ "ONGOING"
#' )
#' }
#' where x is the status_var.
#'
#' @param filter_ds Filter condition for the disposition data.
#'
#'   Filter used to select the relevant disposition data.
#'
#'   Permitted Values: logical expression.
#'
#' @return the input dataset with the disposition status (new_var) added.
#'
#' @keywords adsl
#'
#' @author Samia Kabi
#'
#' @export
#'
#' @examples
#' library(dplyr)
#' data("dm")
#' data("ds")
#'
#' # Default derivation: EOSSTT =
#' #- COMPLETED when status_var = COMPLETED
#' #- DISCONTINUED when status_var != COMPLETED
#' #- ONGOING otherwise
#'
#' derive_disposition_eoxxstt(
#'   dataset = dm,
#'   dataset_ds = ds,
#'   new_var = EOSSTT,
#'   status_var = DSDECOD,
#'   filter_ds = expr(DSCAT == "DISPOSITION EVENT")
#' )
#'
#' # Specific derivation: EOSSTT =
#' #- COMPLETED when status_var = COMPLETED
#' #- DISCONTINUED DUE TO AE when status_var = ADVERSE EVENT
#' #- DISCONTINUED NOT DUE TO AE when status_var != ADVERSE EVENT nor COMPLETED nor missing
#' #- ONGOING otherwise
#'
#' format_eoxxstt1 <- function(x) {
#'   case_when(
#'     x %in% c("COMPLETED") ~ "COMPLETED",
#'     x %in% c("ADVERSE EVENT") ~ "DISCONTINUED DUE TO AE",
#'     !(x %in% c("ADVERSE EVENT", "COMPLETED")) & !is.na(x) ~ "DISCONTINUED NOT DUE TO AE",
#'     TRUE ~ "ONGOING"
#'   )
#' }
#'
#' adsl2 <- derive_disposition_eoxxstt(
#'   dataset = dm,
#'   dataset_ds = ds,
#'   new_var = EOSSTT,
#'   status_var = DSDECOD,
#'   format_new_var = format_eoxxstt1,
#'   filter_ds = expr(DSCAT == "DISPOSITION EVENT")
#' )
#'
#'
#' # Default format
format_eoxxstt_default <- function(x) {
  case_when(
    x %in% c("COMPLETED") ~ "COMPLETED",
    !(x %in% c("COMPLETED")) & !is.na(x) ~ "DISCONTINUED",
    TRUE ~ "ONGOING"
  )
}

derive_disposition_eoxxstt <- function(dataset,
                                       dataset_ds,
                                       new_var,
                                       status_var,
                                       format_new_var = format_eoxxstt_default,
                                       filter_ds) {
  # Checks
  warn_if_vars_exist(dataset, deparse(substitute(new_var)))
  assert_that(is.data.frame(dataset_ds))
  assert_has_variables(dataset_ds, deparse(substitute(status_var)))

  # Process the disposition data
  ds_subset <- dataset_ds %>%
    filter(!!filter_ds) %>%
    mutate(status___ = !!enquo(status_var)) %>%
    select(STUDYID, USUBJID, status___)

  # Expect 1 record per subject in the subsetted DS - issue a warning otherwise
  has_unique_records(
    dataset = ds_subset,
    by_vars = "USUBJID",
    message_type = "warning",
    message = "The filter used for DS results in several records per patient - please check"
  )
  # Add the status var and Derive the new dispo status in the input dataset
  dataset <- dataset %>%
    left_join(ds_subset, by = c("STUDYID", "USUBJID")) %>%
    mutate(
      !!enquo(new_var) := format_new_var(status___)
    ) %>%
    select(-ends_with("___"))
}
