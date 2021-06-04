#' Derive a disposition reason at a specific timepoint
#'
#' Derive a disposition reason from the the relevant records in the disposition domain.
#'
#' @param dataset Input dataset.
#'
#' @param dataset_ds Dataset containing the disposition information (e.g.: `ds`).
#'
#' The variable specified in the `status_var` parameter must be in `dataset_ds`.
#'
#' @param new_var Name of the disposition reason variable.
#'
#' a variable name is expected (e.g. `DCSREAS`).
#'
#' @param reason_var The variable used to derive the disposition reason
#'
#'   A character vector is expected (e.g. `DSDECOD`).
#'
#' @param format_new_vars The format used to derive the reason(s)
#'
#' Default: format_reason_default defined as:
#' format_reason_default<-function(x, y=NULL){
#' case_when (
#'   x == "COMPLETED" ~ x,
#'   x %!in% c("COMPLETED") & ! is.na(x)~"DISCONTINUED",
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
#' @return the input dataset with the disposition reason (`new_var`) added.
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
derive_disposition_reason <- function(dataset,
                                      dataset_ds,
                                      new_var,
                                      status_var,
                                      new_var_spe=NULL,
                                      status_var_spe=NULL,
                                      format_new_vars = format_reason_default,
                                      filter_ds) {


  # Checks
  warn_if_vars_exist(dataset, deparse(substitute(new_var)))
  assert_that(is.data.frame(dataset_ds))
  if (!quo_is_null(enquo(new_var_spe))){
    statusvar<-c(deparse(substitute(status_var)),deparse(substitute(status_var_spe)))
  }
  else {
    statusvar<-deparse(substitute(statusvar))
  }
  assert_has_variables(dataset_ds, statusvar)
  if (!quo_is_null(enquo(new_var_spe))){
    print("Do i enter")
    print(quo_is_null(enquo(new_var_spe)))
    if (!quo_is_null(enquo(status_var_spe))){
      warn(paste("`new_var_spe` is specified but `status_var_spe` is NULL."))
    }
  }
  filter_ds <- enquo(filter_ds)

  # Process the disposition data
  ds_subset <- dataset_ds %>%
    filter(!!filter_ds) %>%
    select(STUDYID, USUBJID, !!enquo(status_var), !!enquo(status_var_spe))

  # Expect 1 record per subject in the subsetted DS - issue a warning otherwise
  has_unique_records(
    dataset = ds_subset,
    by_vars = "USUBJID",
    message_type = "warning",
    message = "The filter used for DS results in several records per patient - please check"
  )
  # Add the status variable and derive the new dispo status in the input dataset
  if (!is.null(new_var_spe)){
    dataset %>%
      left_join(ds_subset, by = c("STUDYID", "USUBJID")) %>%
      mutate(!!enquo(new_var) := format_new_vars(!!enquo(status_var))) #%>%
    mutate(!!enquo(new_var_spe) := format_new_vars(!!enquo(status_var),!!enquo(status_var_spe))) %>%
      select(-!!enquo(status_var),-!!enquo(status_var_spe) )
  }
  else{
    dataset %>%
      left_join(ds_subset, by = c("STUDYID", "USUBJID")) %>%
      mutate(!!enquo(new_var) := format_new_vars(!!enquo(status_var))) %>%
    select(-!!enquo(status_var) )
  }

}

#' Default format for the disposition reason
#'
#' Define a function to map the disposition reason
#'
#' @param x the disposition variable used for the mapping (e.g. `DSDECOD`).
#' @param x the disposition variable used for the mapping of the deatils is required (e.g. `DSTERM`).
format_reason_default <- function(x, y=NULL) {
  out<-if (is.null(y)) x else y
  case_when(
    x != "COMPLETED" & !is.na(x) ~ out,
    TRUE ~ NA_character_
  )
}


