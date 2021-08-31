#' Derive a Disposition Date
#'
#' Derive a disposition status date from the the relevant records in the disposition domain.
#'
#' @param dataset Input dataset
#'
#' @param dataset_ds Datasets containing the disposition information (e.g.: ds)
#'
#' It must contains
#' - `STUDYID`, `USUBJID`,
#' - The variable(s) specified in the `dtc`
#' - The variables used in `filter_ds`.
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
#' @param filter_ds Filter condition for the disposition data.
#'
# 'Filter used to select the relevant disposition data.
# 'It is expected that the filter restricts `dataset_ds` such that there is at most
#' one observation per patient. An error is issued otherwise.
#'
#' Permitted Values: logical expression.
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
#' @return the input dataset with the disposition date (`new_var`) added
#'
#' @keywords adsl timing
#'
#' @author Samia Kabi
#'
#' @export
#'
#' @examples
#' library(dplyr, warn.conflicts = FALSE)
#' data("dm")
#' data("ds")
#' dm %>%
#'   derive_disposition_dt(
#'     dataset_ds = ds,
#'     new_var = FRVDT,
#'     dtc = DSSTDTC,
#'     filter_ds = DSCAT == "OTHER EVENT" & DSDECOD == "FINAL RETRIEVAL VISIT"
#'   ) %>%
#'   select(STUDYID, USUBJID, FRVDT)
derive_disposition_dt <- function(dataset,
                                  dataset_ds,
                                  new_var,
                                  dtc,
                                  filter_ds,
                                  date_imputation = NULL) {

  new_var <- assert_symbol(enquo(new_var))
  dtc <- assert_symbol(enquo(dtc))
  filter_ds <- assert_filter_cond(enquo(filter_ds))
  assert_character_scalar(date_imputation, optional = TRUE)
  assert_data_frame(dataset)
  assert_data_frame(dataset_ds, quo_c(dtc))
  warn_if_vars_exist(dataset, quo_text(new_var))

  # Process the disposition data
  prefix <- sub("\\DT.*", "", deparse(substitute(new_var)))
  newvar <- paste0(prefix, "DT")
  ds_subset <- dataset_ds %>%
    filter(!!filter_ds) %>%
    mutate(datedtc___ = !!enquo(dtc)) %>%
    derive_vars_dt(
      new_vars_prefix = prefix,
      dtc = datedtc___,
      date_imputation = date_imputation,
      flag_imputation = FALSE
    ) %>%
    select(STUDYID, USUBJID, !!enquo(new_var) := !!sym(newvar))

  # Expect 1 record per subject - issue a warning otherwise
  signal_duplicate_records(
    ds_subset,
    by_vars = vars(STUDYID, USUBJID),
    msg = "The filter used for DS results in multiple records per patient"
  )

  # add the new dispo date to the input dataset
  dataset %>%
    left_join(ds_subset, by = c("STUDYID", "USUBJID"))
}
