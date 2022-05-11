#' Derive a Disposition Date
#'
#' @description
#' `r lifecycle::badge("deprecated")`
#'
#' *Deprecated*, please use `derive_var_disposition_dt()` instead.
#'
#' Derive a disposition status date from the the relevant records in the disposition domain.
#'
#' @param dataset Input dataset
#'
#' @param dataset_ds Datasets containing the disposition information (e.g.: `ds`)
#'
#' It must contain:
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
#' Filter used to select the relevant disposition data.
#' It is expected that the filter restricts `dataset_ds` such that there is at most
#' one observation per patient. An error is issued otherwise.
#'
#' Permitted Values: logical expression.
#'
#' @param date_imputation The value to impute the day/month when a datepart is missing.
#'
#'   If `NULL`: no date imputation is performed and partial dates are returned as missing.
#'
#'   Otherwise, a character value is expected, either as a
#'   - format with day and month specified as 'mm-dd': e.g. '06-15' for the 15th
#'   of June
#'   - or as a keyword: 'FIRST', 'MID', 'LAST' to impute to the first/mid/last day/month.
#'
#'   Default is `NULL`
#'
#' @param subject_keys Variables to uniquely identify a subject
#'
#' A list of quosures where the expressions are symbols as returned by
#' `vars()` is expected.
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
#' library(admiraltest)
#' data("dm")
#' data("ds")
#'
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
                                  date_imputation = NULL,
                                  subject_keys = vars(STUDYID, USUBJID)) {
  deprecate_warn("0.6.0", "derive_disposition_dt()", "derive_var_disposition_dt()")
  derive_var_disposition_dt(dataset = dataset,
                            dataset_ds = dataset_ds,
                            new_var = !!enquo(new_var),
                            dtc = !!enquo(dtc),
                            filter_ds = !!enquo(filter_ds),
                            date_imputation = date_imputation,
                            subject_keys = subject_keys)
}

#' Derive a Disposition Date
#'
#' @description
#' `r lifecycle::badge("deprecated")`
#'
#' This function is *deprecated*, please use `derive_vars_merged_dt()` instead.
#'
#' Derive a disposition status date from the the relevant records in the disposition domain.
#'
#' @param dataset Input dataset
#'
#' @param dataset_ds Datasets containing the disposition information (e.g.: ds)
#'
#' It must contain:
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
#' Filter used to select the relevant disposition data.
#' It is expected that the filter restricts `dataset_ds` such that there is at most
#' one observation per patient. An error is issued otherwise.
#'
#' Permitted Values: logical expression.
#'
#' @param date_imputation The value to impute the day/month when a datepart is missing.
#'
#'   If `NULL`: no date imputation is performed and partial dates are returned
#'   as missing.
#'
#'   Otherwise, a character value is expected, either as a
#'   - format with day and month specified as 'mm-dd': e.g. '06-15' for the 15th
#'   of June
#'   - or as a keyword: 'FIRST', 'MID', 'LAST' to impute to the first/mid/last day/month.
#'
#'   Default is `NULL`
#'
#' @param subject_keys Variables to uniquely identify a subject
#'
#' A list of quosures where the expressions are symbols as returned by
#' `vars()` is expected.
#'
#' @inheritParams impute_dtc
#'
#' @return the input dataset with the disposition date (`new_var`) added
#'
#' @keywords adsl timing
#'
#' @author Samia Kabi
#'
#' @export
#'
derive_var_disposition_dt <- function(dataset,
                                  dataset_ds,
                                  new_var,
                                  dtc,
                                  filter_ds,
                                  date_imputation = NULL,
                                  preserve = FALSE,
                                  subject_keys = vars(STUDYID, USUBJID)) {

  new_var <- assert_symbol(enquo(new_var))
  dtc <- assert_symbol(enquo(dtc))
  filter_ds <- assert_filter_cond(enquo(filter_ds))
  assert_character_scalar(date_imputation, optional = TRUE)
  assert_data_frame(dataset)
  assert_data_frame(dataset_ds, quo_c(dtc))
  warn_if_vars_exist(dataset, quo_text(new_var))
  assert_vars(subject_keys)
  assert_logical_scalar(preserve)
  deprecate_warn("0.7.0", "derive_var_disposition_dt()", "derive_vars_merged_dt()")

  derive_vars_merged_dt(
    dataset,
    dataset_add = dataset_ds,
    filter_add = !!filter_ds,
    new_vars_prefix = "temp_",
    by_vars = subject_keys,
    dtc = !!dtc,
    date_imputation = date_imputation,
    flag_imputation = FALSE,
    preserve = preserve,
    duplicate_msg = "The filter used for DS results in multiple records per patient."
  ) %>%
    rename(!!new_var := temp_DT)
}
