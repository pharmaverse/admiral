#' Derive Datetime of First Exposure to Treatment
#'
#' @description
#' `r lifecycle::badge("questioning")`
#'
#' Derives datetime of first exposure to treatment (`TRTSDTM`)
#'
#' @param dataset Input dataset
#'
#'   The variables specified by the `by_vars` parameter are expected.
#'
#' @param dataset_ex `ex` dataset
#'
#'   The variables `EXSTDTC`, `EXSEQ`, and those specified by the `filter_ex`
#'   parameter are expected.
#'
#' @param filter_ex Filter condition for the ex dataset
#'
#'   Only observations of the ex dataset which fulfill the specified condition
#'   are considered for the treatment start date.
#'
#'   Default: `EXDOSE > 0 | (EXDOSE == 0 & str_detect(EXTRT, 'PLACEBO')) & nchar(EXSTDTC) >= 10`
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
#' @author Stefan Bundfuss
#'
#' @return The input dataset with `TRTSDTM` variable added
#'
#' @export
#'
#' @keywords adsl timing derivation
#'
#' @examples
#' library(dplyr, warn.conflicts = FALSE)
#' library(admiraltest)
#' data("admiral_ex")
#' data("admiral_dm")
#'
#' admiral_dm %>%
#'   derive_var_trtsdtm(dataset_ex = admiral_ex) %>%
#'   select(USUBJID, TRTSDTM)
derive_var_trtsdtm <- function(dataset,
                               dataset_ex,
                               filter_ex = (EXDOSE > 0 | (EXDOSE == 0 & str_detect(EXTRT, "PLACEBO"))) & nchar(EXSTDTC) >= 10, # nolint
                               subject_keys = vars(STUDYID, USUBJID)) {
  assert_data_frame(dataset, subject_keys)
  assert_data_frame(dataset_ex, required_vars = quo_c(subject_keys, vars(EXSTDTC, EXSEQ)))
  filter_ex <- assert_filter_cond(enquo(filter_ex), optional = TRUE)

  derive_vars_merged_dtm(
    dataset,
    dataset_add = dataset_ex,
    filter_add = !!filter_ex,
    new_vars_prefix = "TRTS",
    dtc = EXSTDTC,
    date_imputation = "first",
    time_imputation = "first",
    flag_imputation = "none",
    order = vars(TRTSDTM, EXSEQ),
    mode = "first",
    by_vars = subject_keys
  )
}
