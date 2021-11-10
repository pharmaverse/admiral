#' Derive Datetime of Last Exposure to Treatment
#'
#' Derives datetime of last exposure to treatment (`TRTEDTM`)
#'
#' @param dataset Input dataset
#'
#'   The variables specified by the `by_vars` parameter are expected.
#'
#' @param dataset_ex `ex` dataset
#'
#'   The variables `EXENDTC`, `EXSEQ`, and those specified by the `filter_ex`,
#'   `dtc` and `subject_keys`parameters are expected.
#'
#' @param filter_ex Filter condition for the ex dataset
#'
#'   Only observations of the ex dataset which fulfill the specified condition
#'   are considered for the treatment start date.
#'
#'   Default: NULL
#'
#'   Permitted Values: logical expression
#'
#' @param subject_keys Variables to uniquely identify a subject
#'
#' A list of quosures where the expressions are symbols as returned by
#' `vars()` is expected.
#'
#' @param new_var name of treatment end date e.g. `TRTEDTM`
#'
#' @param dtc name of variable from which end date is derived e.g. `EXENDTC`
#'
#' @param date_imputation The value to impute the day/month when a datepart is
#'   missing.
#'
#'   If `NULL`: no date imputation is performed and partial dates are returned as
#'   missing.
#'
#'   Otherwise, a character value is expected, either as a
#'   - format with month and day specified as 'mm-dd': e.g. `"06-15"` for the 15th
#'   of June,
#'   - or as a keyword: `"FIRST"`, `"MID"`, `"LAST"` to impute to the first/mid/last
#'   day/month.
#'
#'   Default is `NULL`.
#'
#' @param time_imputation The value to impute the time when a timepart is
#'   missing.
#'
#'   A character value is expected, either as a
#'   - format with hour, min and sec specified as 'hh:mm:ss': e.g. `"00:00:00"`
#'   for the start of the day,
#'   - or as a keyword: `"FIRST"`,`"LAST"` to impute to the start/end of a day.
#'
#'   Default is `"00:00:00"`.
#'
#' @param min_dates Minimum dates
#'  A list of dates is expected. It is ensured that the imputed date is not
#'  before any of the specified dates, e.g., that the imputed adverse event start
#'  date is not before the first treatment date. Only dates which are in the
#'  range of possible dates of the dtc value are considered. The possible dates
#'  are defined by the missing parts of the dtc date (see example below). This
#'  ensures that the non-missing parts of the dtc date are not changed.
#'
#' @param max_dates Maximum dates
#'
#' A list of dates is expected. It is ensured that the imputed date is not after
#' any of the specified dates, e.g., that the imputed date is not after the data
#' cut off date. Only dates which are in the range of possible dates are
#' considered.
#'
#' @param ord_vars
#'  Order of variables for which either the last or first observation
#'  is picked as treatment start date
#'  e.g. vars(EXENDTC,EXSEQ)
#'
#' @param ord_filter
#'  The mode by which the filter_extreme is chosen.
#'  If "first", the first observation of each by group is included in the output dataset.
#'  If "last", the last observation of each by group is included in the output dataset.
#'
#'  Permitted Values: "first", "last"
#'
#' @details For each group (with respect to the variables specified for the
#'   `by_vars` parameter) the first observation (with respect to the order
#'   specified for the `order` parameter) is included in the output dataset.
#'
#' @author Stefan Bundfuss, Teckla Akinyi
#'
#' @return The input dataset with start date variable added e.g. `TRTEDTM`
#'
#' @export
#'
#' @keywords adsl timing derivation
#'
#' @examples
#' library(dplyr, warn.conflicts = FALSE)
#' data("ex")
#' data("dm")
#'
#' dm %>%
#'   derive_vars_trt_end(dataset_ex = ex
#'                         filter_ex =  (EXDOSE > 0 | (EXDOSE == 0 & str_detect(EXTRT, "PLACEBO"))) & nchar(EXSTDTC) >= 10, # nolint
#'                         subject_keys = vars(STUDYID, USUBJID),
#'                         new_var=TRTEDTM,
#'                         dtc=EXENDTC,
#'                         date_imputation = "last",
#'                         time_imputation="last",
#'                         min_dates = NULL,
#'                         max_dates = NULL,
#'                         ord_vars = vars(EXENDTC,EXSEQ),
#'                         ord_filter="last") %>%
#'   select(USUBJID, TRTEDTM)

derive_vars_trt_end <- function(dataset,
                                  dataset_ex,
                                  filter_ex = NULL,
                                  subject_keys,
                                  new_var,
                                  dtc,
                                  date_imputation,
                                  time_imputation,
                                  min_dates = NULL,
                                  max_dates = NULL,
                                  ord_vars,
                                  ord_filter) {
  new_var <- assert_symbol(enquo(new_var))
  dtc <- assert_symbol(enquo(dtc))
  ord <- assert_symbol(enquo(ord))
  assert_data_frame(dataset,subject_keys)
  assert_data_frame(dataset_ex, required_vars = quo_c(subject_keys, ord_vars))
  assert_character_scalar(ord_filter, values = c("first","last"))

  filter_ex <- assert_filter_cond(enquo(filter_ex), optional = TRUE)

  add <- dataset_ex %>%
    mutate(imputed_dtc = convert_dtc_to_dtm(!!dtc,
                                            date_imputation = date_imputation,
                                            time_imputation = time_imputation,
                                            min_dates = min_dates,
                                            max_dates = max_dates)
    )
  if (!quo_is_null(filter_pre_timepoint)){
    add %>% filter_if(filter_ex) %>%
      filter_extreme(
        order = vars(!!dtc, !!ord_vars),
        by_vars = !!!subject_keys,
        mode = ord_filter) %>%
      mutate(!!new_var := imputed_dtc)
  } else {
    filter_extreme(
      order = ord,
      by_vars = !!!subject_keys,
      mode = ord_filter) %>%
      mutate(!!new_var := imputed_dtc)
  }

  left_join(dataset, add, by = vars2chr(subject_keys)) %>%
    select(-imputed_dtc)
}

