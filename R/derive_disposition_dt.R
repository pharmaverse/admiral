#' Derive --DTM (and --DTF/--TMF)
#'
#' Derive --DTM based on --DTC. --DTM is imputed based on user input
#' Derive --DTF/--TMF if needed based on --DTC and --DT
#'
#' @param dataset Input dataset
#'
#' @param dataset_ds Datasets containiung the disposition information
#' (usually: ds)
#'
#' The variable specified in the dtc parameter must be in dataset_ds
#'
#'
#' @param new_var Name of the disposition date variable
#'
#' a variable name is expected
#'
#' @param dtc The --DTC date used to derive/impute --DT
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
                                  filter = NULL,
                                  date_imputation = NULL) {
  # check if dataset_ds exists
  assert_dataset_exist(deparse(substitute(dataset_ds)))
  # Check DTC is present in dataset_ds
  assert_has_variables(dataset_ds, deparse(substitute(dtc)))

  # Warn if the variable to derive already exists in the input dataset
  warn_if_vars_exist(dataset, deparse(substitute(new_var)))

  # if DS needs to be filtered, filter
  if (!is.null(filter)) {
    ds_subset <- dataset_ds %>%
      filter(!!!filter)
  }
  else {
    ds_subset <- dataset_ds
  }

  # only 1 record per subject is expected - issue a warning otherwise
  assert_has_unique_records(
    dataset = ds_subset,
    by_vars = "USUBJID",
    message_type = "warning",
    message = "the filter used for DS results in several records per patient - please check"
  )

  # set DTC in datedtc (resolves in mutate)
  ds_subset <- ds_subset %>%
    mutate(datedtc = !!enquo(dtc))

  # Prefix to use in derive_vars_dt ("RFIC--", "ENRL--",...)
  prefix <- sub("\\DT.*", "", deparse(substitute(new_var)))
  newvar <- paste0(prefix, "DT")
  # Create the new dispo date
  ds__ <- derive_vars_dt(
    ds_subset,
    new_vars_prefix = prefix,
    dtc = datedtc,
    date_imputation = date_imputation,
    flag_imputation = FALSE
  ) %>%
    select(STUDYID, USUBJID, !!enquo(new_var) := !!sym(newvar))

  # add the new dispo date to the input dataset
  dataset %>%
    left_join(ds__, by = c("STUDYID", "USUBJID"))
}
