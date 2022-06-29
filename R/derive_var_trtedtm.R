#' Derive Datetime of Last Exposure to Treatment
#'
#' @description
#' `r lifecycle::badge("deprecated")`
#'
#' This function is *deprecated*, please use `derive_vars_merged_dtm()` instead.
#'
#' Derives datetime of last exposure to treatment (`TRTEDTM`)
#'
#' @param dataset Input dataset
#'
#'   The variables specified by the `by_vars` parameter are expected.
#'
#' @param dataset_ex `ex` dataset
#'
#'   The variables `EXENDTC`, `EXSEQ`, and those specified by the `filter_ex`
#'   parameter are expected.
#'
#' @param filter_ex Filter condition for the ex dataset
#'
#'   Only observations of the ex dataset which fulfill the specified condition
#'   are considered for the treatment start date.
#'
#'   Default: `EXDOSE > 0 | (EXDOSE == 0 & str_detect(EXTRT, 'PLACEBO')) & nchar(EXENDTC) >= 10`
#'
#'   Permitted Values: logical expression
#'
#' @param subject_keys Variables to uniquely identify a subject
#'
#' A list of quosures where the expressions are symbols as returned by
#' `vars()` is expected.
#'
#' @details For each group (with respect to the variables specified for the
#'   `by_vars` parameter) the first observation (with respect to the order
#'   specified for the `order` parameter) is included in the output dataset.
#'
#' @return The input dataset with `TRTEDTM` variable added
#'
#' @author Stefan Bundfuss
#'
#' @keywords adsl timing derivation
#'
#' @export
#'
derive_var_trtedtm <- function(dataset,
                               dataset_ex,
                               filter_ex = (EXDOSE > 0 | (EXDOSE == 0 & str_detect(EXTRT, "PLACEBO"))) & nchar(EXENDTC) >= 10, # nolint
                               subject_keys = vars(STUDYID, USUBJID)) {
  assert_data_frame(dataset, subject_keys)
  assert_data_frame(dataset_ex, required_vars = quo_c(subject_keys, vars(EXENDTC, EXSEQ)))
  filter_ex <- assert_filter_cond(enquo(filter_ex), optional = TRUE)
  deprecate_warn("0.7.0", "derive_var_trtedtm()", "derive_vars_merged_dtm()")

  derive_vars_merged_dtm(
    dataset,
    dataset_add = dataset_ex,
    filter_add = !!filter_ex,
    new_vars_prefix = "TRTE",
    dtc = EXENDTC,
    date_imputation = "last",
    time_imputation = "last",
    flag_imputation = "none",
    order = vars(TRTEDTM, EXSEQ),
    mode = "last",
    by_vars = subject_keys
  )
}
