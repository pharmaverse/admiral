#' Derive a disposition status at a specific timepoint
#'
#' Derive a disposition status from the the relevant records in the disposition domain.
#'
#' @param dataset Input dataset.
#'
#' @param dataset_ds Dataset containing the disposition information (e.g.: ds).
#'
#' It must contains
#' - `STUDYID`, `USUBJID`,
#' - The variable(s) specified in the `status_var`
#' - The variables used in `filter_ds`.
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
#' where `x` is the `status_var.`
#'
#' @param filter_ds Filter condition for the disposition data.
#'
# 'Filter used to select the relevant disposition data.
# 'It is expected that the filter restricts `dataset_ds` such that there is at most
#' one observation per patient. An error is issued otherwise.
#'
#' Permitted Values: logical expression.
#'
#' @return the input dataset with the disposition status (`new_var`) added.
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
#' derive_disposition_status(
#'   dataset = dm,
#'   dataset_ds = ds,
#'   new_var = EOSSTT,
#'   status_var = DSDECOD,
#'   filter_ds = DSCAT == "DISPOSITION EVENT"
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
#'     x == "COMPLETED" ~ "COMPLETED",
#'     x == "ADVERSE EVENT" ~ "DISCONTINUED DUE TO AE",
#'     x %!in% c("ADVERSE EVENT", "COMPLETED") & !is.na(x) ~ "DISCONTINUED NOT DUE TO AE",
#'     TRUE ~ "ONGOING"
#'   )
#' }
#'
#' derive_disposition_status(
#'   dataset = dm,
#'   dataset_ds = ds,
#'   new_var = EOSSTT,
#'   status_var = DSDECOD,
#'   format_new_var = format_eoxxstt1,
#'   filter_ds = DSCAT == "DISPOSITION EVENT"
#' )
derive_disposition_status <- function(dataset,
                                      dataset_ds,
                                      new_var,
                                      status_var,
                                      format_new_var = format_eoxxstt_default,
                                      filter_ds) {

  new_var <- assert_symbol(enquo(new_var))
  status_var <- assert_symbol(enquo(status_var))
  filter_ds <- assert_filter_cond(enquo(filter_ds))
  assert_that(is.function(format_new_var))
  assert_data_frame(dataset)
  assert_data_frame(dataset_ds, quo_c(status_var))
  warn_if_vars_exist(dataset, quo_text(new_var))

  # Process the disposition data
  ds_subset <- dataset_ds %>%
    filter(!!filter_ds) %>%
    select(STUDYID, USUBJID, !!enquo(status_var))

  # Expect 1 record per subject in the subsetted DS - issue a warning otherwise
  signal_duplicate_records(
    ds_subset,
    by_vars = vars(STUDYID, USUBJID),
    msg = "The filter used for DS results in multiple records per patient"
  )

  # Add the status variable and derive the new dispo status in the input dataset
  dataset %>%
    left_join(ds_subset, by = c("STUDYID", "USUBJID")) %>%
    mutate(!!enquo(new_var) := format_new_var(!!enquo(status_var))) %>%
    select(-!!enquo(status_var))
}

#' Default format for the disposition status
#'
#' Define a function to map the disposition status.
#'
#' @param x the disposition variable used for the mapping (e.g. `DSDECOD`).
format_eoxxstt_default <- function(x) {
  case_when(
    x == "COMPLETED" ~ "COMPLETED",
    x != "COMPLETED" & !is.na(x) ~ "DISCONTINUED",
    TRUE ~ "ONGOING"
  )
}
