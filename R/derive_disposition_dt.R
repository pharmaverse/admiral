#' Derive a disposition date.
#'
#' Derive a disposition status date from the the relevant records in the disposition domain.
#'
#' @param dataset Input dataset
#'
#' @param dataset_ds Datasets containing the disposition information (e.g.: ds)
#'
#' The variable specified in dtc parameter must be in dataset_ds
#'
#' @param new_var Name of the disposition date variable
#'
#' a variable name is expected
#'
#' @param dtc The character date used to derive/impute the disposition date
#'
#'   A character date is expected in a format like yyyy-mm-dd or yyyy-mm-ddThh:mm:ss.
#'   If the year part is not recorded (missing date), no imputation is performed.
#'
#' @param date_imputation The value to impute the day/month when a datepart is missing.
#'
#'   If NULL: no date imputation is performed and partial dates are returned as missing.
#'
#'   Otherwise, a character value is expected, either as a
#'   - format with day and month specified as 'dd-mm': e.g. '15-06' for the 15th of June
#'   - or as a keyword: 'FIRST', 'MID', 'LAST' to impute to the first/mid/last day/month
#'
#'   Default is NULL
#'
#' @return the input dataset with the disposition date (new_var) added
#'
#' @keywords adsl timing
#'
#' @author Samia Kabi
#'
#' @export
#'
#' @examples
#' derive_disposition_dt(
#'   dataset = dm,
#'   dataset_ds = ds,
#'   new_var = RFICDT,
#'   dtc = DSSTDTC,
#'   filter = expr(DSCAT == "PROTOCOL MILESTONE" & DSDECOD == "INFORMED CONSENT OBTAINED")
#' )
derive_disposition_dt <- function(dataset,
                                  dataset_ds,
                                  new_var,
                                  dtc,
                                  filter,
                                  date_imputation = NULL) {

  # Checks
  warn_if_vars_exist(dataset, deparse(substitute(new_var)))
  assert_that(is.data.frame(dataset_ds))
  assert_has_variables(dataset_ds, deparse(substitute(dtc)))
  # Expect 1 record per subject - issue a warning otherwise
  if (check_unique_record) {
    assert_has_unique_records(
      dataset = ds_subset,
      by_vars = "USUBJID",
      message_type = "warning",
      message = "the filter used for DS results in several records per patient - please check"
    )
  }

  # Process the disposition data
  prefix <- sub("\\DT.*", "", deparse(substitute(new_var)))
  newvar <- paste0(prefix, "DT")
  ds_subset <- dataset %>%
    filter(!!filter) %>%
    mutate(datedtc = !!enquo(dtc)) %>%
    derive_vars_dt(
      new_vars_prefix = prefix,
      dtc = datedtc,
      date_imputation = date_imputation,
      flag_imputation = FALSE
    ) %>%
    select(STUDYID, USUBJID, !!enquo(new_var) := !!sym(newvar))

  # add the new dispo date to the input dataset
  dataset %>%
    left_join(ds_subset, by = c("STUDYID", "USUBJID"))
}
