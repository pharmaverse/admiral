#' Derive datetime of Last Exposure to Treatment
#'
#' Derives datetime of Last Exposure to Treatment (`TRTEDTM`)
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
                               filter_ex = (EXDOSE > 0 | (EXDOSE == 0 & str_detect(EXTRT, "PLACEBO"))) & nchar(EXENDTC) >= 10) { # nolint

  assert_data_frame(dataset, vars(USUBJID))
  assert_data_frame(dataset_ex, vars(USUBJID, EXENDTC, EXSEQ))
  filter_ex <- assert_filter_cond(enquo(filter_ex), optional = TRUE)

  if (!quo_is_null(filter_ex)) {
    add <- filter(dataset_ex, !!filter_ex)
  } else {
    add <- dataset_ex
  }
  add <- add %>%
    filter_extreme(
      order = vars(EXENDTC, EXSEQ),
      by_vars = vars(USUBJID),
      mode = "last"
    ) %>%
    transmute(
      USUBJID,
      TRTEDTM = convert_dtc_to_dtm(impute_dtc(EXENDTC, time_imputation = "LAST"))
    )

  left_join(dataset, add, by = c("USUBJID"))
}
