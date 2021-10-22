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
#'   The variables `EXENDTC`, `EXSEQ`, and those specified by the `filter_ex`
#'   parameter are expected.
#'
#' @param filter_ex Fiter condition for the ex dataset
#'
#'   Only observations of the ex dataset which fulfill the specified condition
#'   are considered for the treatment start date.
#'
#'   Default: `EXDOSE > 0 | (EXDOSE == 0 & str_detect(EXTRT, 'PLACEBO')`
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
#' @return The input dataset with `TRTEDTM` variable added
#'
#' @keywords adsl timing derivation
#'
#' @export
#'
#' @examples
#' library(dplyr, warn.conflicts = FALSE)
#' data("ex")
#' data("dm")
#'
#' dm %>%
#'   derive_var_trtedtm(dataset_ex = ex) %>%
#'   select(USUBJID, TRTEDTM)
derive_var_trtedtm <- function(dataset,
                               dataset_ex,
                               filter_ex = (EXDOSE > 0 | (EXDOSE == 0 & str_detect(EXTRT, "PLACEBO"))) & nchar(EXENDTC) >= 10 , # nolint
                               subject_keys = vars(STUDYID, USUBJID)) {

  assert_data_frame(dataset, subject_keys)
  assert_data_frame(dataset_ex, required_vars = quo_c(subject_keys, vars(EXENDTC, EXSEQ)))
  filter_ex <- assert_filter_cond(enquo(filter_ex), optional = TRUE)

  add <- dataset_ex %>%
    filter_if(filter_ex) %>%
    filter_extreme(
      order = vars(EXENDTC, EXSEQ),
      by_vars = subject_keys,
      mode = "last"
    ) %>%
    transmute(
      !!!subject_keys,
      TRTEDTM = convert_dtc_to_dtm(EXENDTC, date_imputation = "last", time_imputation = "last")
    )

  left_join(dataset, add, by = vars2chr(subject_keys))
}
