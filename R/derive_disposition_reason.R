#' Derive a disposition reason at a specific timepoint
#'
#' Derive a disposition reason from the the relevant records in the disposition domain.
#'
#' @param dataset Input dataset.
#'
#' @param dataset_ds Dataset containing the disposition information (e.g.: `ds`).
#'
#' It must contains
#' - `STUDYID`, `USUBJID`,
#' - The variable specified in the `reason_var`, (and `reason_var_spe` if required)
#' - The variables used in `filter_ds`.
#'
#' @param new_var Name of the disposition reason variable.
#'
#' A variable name is expected (e.g. `DCSREAS`).
#'
#' @param reason_var The variable used to derive the disposition reason
#'
#' A variable name is expected (e.g. `DSDECOD`).
#'
#' @param new_var_spe Name of the disposition reason detail variable.
#'
#' A variable name is expected (e.g. `DCSREASP`).
#' If `new_var_spe` is specified, it is expected that `reason_var_spe` is also specified,
#' otherwise a warning is issued and
#'
#' @param reason_var_spe The variable used to derive the disposition reason detail
#'
#' A variable name is expected (e.g. `DSTERM`).
#'
#' @param format_new_vars The function used to derive the reason(s)
#'
#' Default: format_reason_default defined as:
#' format_reason_default <- function(reason, reason_spe = NULL) {
#'   out <- if ( is.null(reason_spe) ) reason else reason_spe
#'   if_else ( reason != "COMPLETED" & !is.na(reason), out, NA_character_)
#' }
#' format_reason_default(DSDECOD) returns `DSDECOD` when `DSDECOD` is not 'COMPLETED' nor NA.
#' format_reason_default(DSDECOD, DSTERM) returns `DSTERM` when `DSDECOD` is not 'COMPLETED' nor NA.
#'
#' @param filter_ds Filter condition for the disposition data.
#'
#'   Filter used to select the relevant disposition data.
#'   It is expected that the filter restricts `dataset_ds` such that there is at most
#'   one observation per patient. A warning is issued otherwise.
#'
#'   Permitted Values: logical expression.
#'
#' @return the input dataset with the disposition reason(s) (`new_var and
#' if required `new_var_spe`)  added.
#'
#' @details
#' This functions returns the main reason for discontinuation (e.g. `DCSREAS` or `DCTREAS`).
#' The reason for discontinuation is derived based on `reason_var` (e.g. `DSDECOD`) and
#' `format_new_vars`.
#' If `new_var_spe` is not NULL, then the function will also returns the details associated
#' with the reason for discontinuation (e.g. `DCSREASP`).
#' The details associated with the reason for discontinuation are derived based on
#' `reason_var_spe` (e.g. `DSTERM`), `reason_var` and `format_new_vars`.
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
#'   format_new_vars = format_dcsreas,
#'   filter_ds = DSCAT == "DISPOSITION EVENT"
#' )
derive_disposition_reason <- function(dataset,
                                      dataset_ds,
                                      new_var,
                                      reason_var,
                                      new_var_spe = NULL,
                                      reason_var_spe = NULL,
                                      format_new_vars = format_reason_default,
                                      filter_ds) {
  new_var <- enquo(new_var)
  reason_var <- enquo(reason_var)
  new_var_spe <- enquo(new_var_spe)
  reason_var_spe <- enquo(reason_var_spe)
  filter_ds <- enquo(filter_ds)

  # Checks
  warn_if_vars_exist(dataset, quo_text(new_var))
  assert_that(is.data.frame(dataset_ds))
  if (!quo_is_null(new_var_spe)) {
    if (!quo_is_null(reason_var_spe)) {
      statusvar <- c(quo_text(reason_var), quo_text(reason_var_spe))
    }
    else {
      err_msg <- paste(
        "`new_var_spe` is specified as ", quo_text(new_var_spe),
        "but `reason_var_spe` is NULL.",
        "Please specifiy `reason_var_spe` together with `new_var_spe`."
      )
      abort(err_msg)
    }
  }
  else {
    statusvar <- quo_text(reason_var)
  }
  assert_has_variables(dataset_ds, statusvar)

  # Process the disposition data
  ds_subset <- dataset_ds %>%
    filter(!!filter_ds) %>%
    select(STUDYID, USUBJID, !!reason_var, !!reason_var_spe)

  # Expect 1 record per subject in the subsetted DS - issue a warning otherwise
  has_unique_records(
    dataset = ds_subset,
    by_vars = vars(STUDYID, USUBJID),
    message_type = "warning",
    message = "The filter used for DS results in several records per patient - please check"
  )
  # Add the status variable and derive the new dispo reason(s) in the input dataset
  if (!quo_is_null(new_var_spe)) {
    dataset %>%
      left_join(ds_subset, by = c("STUDYID", "USUBJID")) %>%
      mutate(!!new_var := format_new_vars(!!reason_var)) %>%
      mutate(!!new_var_spe := format_new_vars(!!reason_var, !!reason_var_spe)) %>%
      select(-statusvar)
  }
  else {
    dataset %>%
      left_join(ds_subset, by = c("STUDYID", "USUBJID")) %>%
      mutate(!!new_var := format_new_vars(!!reason_var)) %>%
      select(-statusvar)
  }
}

#' Default format for the disposition reason
#'
#' Define a function to map the disposition reason
#'
#' @param reason the disposition variable used for the mapping (e.g. `DSDECOD`).
#' @param reason_spe the disposition variable used for the mapping of the details
#' if required (e.g. `DSTERM`).
#' @details
#' format_reason_default(DSDECOD) returns `DSDECOD` when `DSDECOD` != 'COMPLETED' nor NA.
#' format_reason_default(DSDECOD, DSTERM) returns `DSTERM` when `DSDECOD` != 'COMPLETED' nor NA.
#' e.g. DCSREAS =  format_reason_default(DSDECOD)
#' e.g. DCSREASP =  format_reason_default(DSDECOD, DSTERM)

format_reason_default <- function(reason, reason_spe = NULL) {
  out <- if (is.null(reason_spe)) reason else reason_spe
  if_else(reason != "COMPLETED" & !is.na(reason), out, NA_character_)
}
