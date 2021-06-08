#' Derive a disposition reason at a specific timepoint
#'
#' Derive a disposition reason from the the relevant records in the disposition domain.
#'
#' @param dataset Input dataset.
#'
#' @param dataset_ds Dataset containing the disposition information (e.g.: `ds`).
#'
#' The variable specified in the `reason_var` parameter must be in `dataset_ds`.
#'
#' @param new_var Name of the disposition reason variable.
#'
#' a variable name is expected (e.g. `DCSREAS`).
#'
#' @param reason_var The variable used to derive the disposition reason
#'
#'   A character vector is expected (e.g. `DSDECOD`).
#'
#' @param new_var_spe Name of the disposition reason detail variable.
#'
#' a variable name is expected (e.g. `DCSREASP`).
#'
#' @param reason_var_spe The variable used to derive the disposition reason detail
#'
#'   A character vector is expected (e.g. `DSTERM`).
#'
#' @param format_new_var The format used to derive the reason(s)
#'
#' Default: format_reason_default defined as:
#' format_reason_default<-function(x, y=NULL){
#' case_when (
#'   x == "COMPLETED" ~ x,
#'   TRUE ~ NA_character
#' )
#' }
#' where x is the reason_var.
#'
#' @param filter_ds Filter condition for the disposition data.
#'
#'   Filter used to select the relevant disposition data.
#'
#'   Permitted Values: logical expression.
#'
#' @return the input dataset with the disposition reason (`new_var`) added.
#'
#' @details
#' This functions returns the main reason for discontinuation (e.g. `DCSREAS` or `DCTREAS`).
#' The reason for discontinuation is derived based on `reason_var` (e.g. `DSDECOD`) and
#' `format_new_var`.
#' If `new_var_spe` is not NULL, then the function will also returns the details associated
#' with the reason for discontinuation (e.g. `DCSREASP`).
#' The details associated with the reason for discontinuation are derived based on
#' `reason_var_spe` (e.g. `DSTERM`), `reason_var` and `format_new_var`. see
#'
#' @seealso [`format_reason_default()`]
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
#' # Derive DCSREAS using the default format
#' derive_disposition_reason(
#'   dataset = dm,
#'   dataset_ds = ds,
#'   new_var = DCSREAS,
#'   reason_var = DSDECOD,
#'   filter_ds = DSCAT == "DISPOSITION EVENT"
#' )
#' # Derive DCSREAS and DCSREASP using a study-specific format
#' format_dcsreas <- function(x, y = NULL) {
#'   out <- if (is.null(y)) x else y
#'   case_when(
#'     x %!in% c("COMPLETED", "SCREEN FAILURE") & !is.na(x) ~ out,
#'     TRUE ~ NA_character_
#'   )
#' }
#' derive_disposition_reason(
#'   dataset = dm,
#'   dataset_ds = ds,
#'   new_var = DCSREAS,
#'   reason_var = DSDECOD,
#'   new_var_spe = DCSREASP,
#'   reason_var_spe = DSTERM,
#'   format_new_var = format_dcsreas,
#'   filter_ds = DSCAT == "DISPOSITION EVENT"
#' )
derive_disposition_reason <- function(dataset,
                                      dataset_ds,
                                      new_var,
                                      reason_var,
                                      new_var_spe = NULL,
                                      reason_var_spe = NULL,
                                      format_new_var = format_reason_default,
                                      filter_ds) {

  # Checks
  warn_if_vars_exist(dataset, deparse(substitute(new_var)))
  assert_that(is.data.frame(dataset_ds))
  if (!quo_is_null(enquo(new_var_spe))) {
    statusvar <- c(deparse(substitute(reason_var)), deparse(substitute(reason_var_spe)))
  }
  else {
    statusvar <- deparse(substitute(reason_var))
  }
  assert_has_variables(dataset_ds, statusvar)
  if (!quo_is_null(enquo(new_var_spe))) {
    warn_if_vars_exist(dataset, deparse(substitute(new_var_spe)))
    if (quo_is_null(enquo(reason_var_spe))) {
      warn(paste(
        "`new_var_spe` is specified",
        deparse(substitute(new_var_spe)),
        "but `reason_var_spe` is NULL."
      ))
    }
  }

  new_var <- enquo(new_var)
  reason_var <- enquo(reason_var)
  new_var_spe <- enquo(new_var_spe)
  reason_var_spe <- enquo(reason_var_spe)
  filter_ds <- enquo(filter_ds)

  # Process the disposition data
  ds_subset <- dataset_ds %>%
    filter(!!filter_ds) %>%
    select(STUDYID, USUBJID, !!reason_var, !!reason_var_spe)

  # Expect 1 record per subject in the subsetted DS - issue a warning otherwise
  has_unique_records(
    dataset = ds_subset,
    by_vars = "USUBJID",
    message_type = "warning",
    message = "The filter used for DS results in several records per patient - please check"
  )
  # Add the status variable and derive the new dispo reason(s) in the input dataset
  if (!quo_is_null(new_var_spe)) {
    dataset %>%
      left_join(ds_subset, by = c("STUDYID", "USUBJID")) %>%
      mutate(!!new_var := format_new_var(!!reason_var)) %>%
      mutate(!!new_var_spe := format_new_var(!!reason_var, !!reason_var_spe)) %>%
      select(-!!reason_var, -!!reason_var_spe)
  }
  else {
    dataset %>%
      left_join(ds_subset, by = c("STUDYID", "USUBJID")) %>%
      mutate(!!new_var := format_new_var(!!reason_var)) %>%
      select(-!!reason_var)
  }
}

#' Default format for the disposition reason
#'
#' Define a function to map the disposition reason
#'
#' @param x the disposition variable used for the mapping (e.g. `DSDECOD`).
#' @param y the disposition variable used for the mapping of the details if required (e.g. `DSTERM`).
format_reason_default <- function(x, y = NULL) {
  out <- if (is.null(y)) x else y
  case_when(
    x != "COMPLETED" & !is.na(x) ~ out,
    TRUE ~ NA_character_
  )
}
