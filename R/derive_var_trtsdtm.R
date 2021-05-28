#' Derive datetime of First Exposure to Treatment
#'
#' Derives datetime of First Exposure to Treatment (`TRTSDTM`)
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
#' @param filter_ex Fiter condition for the ex dataset
#'
#'   Only observations of the ex dataset which fulfill the specified condition
#'   are considered for the treatment start date.
#'
#'   Default: `exprs(EXDOSE > 0 | (EXDOSE == 0 & str_detect(EXTRT, 'PLACEBO'))`
#'
#'   Permitted Values: logical expression
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
#' library(dplyr)
#' data("ex")
#' data("dm")
#'
#' derive_var_trtsdtm(dm, dataset_ex = ex) %>%
#'   select(USUBJID, TRTSDTM)
derive_var_trtsdtm <- function(dataset,
                               dataset_ex,
                               filter_ex = exprs((EXDOSE > 0 | (EXDOSE == 0 & str_detect(EXTRT, "PLACEBO"))) & nchar(EXSTDTC) >= 10)) { # nolint

  assert_has_variables(dataset, c("USUBJID"))
  assert_has_variables(dataset_ex, c("USUBJID", "EXSTDTC", "EXSEQ"))

  if (!is.null(filter_ex)) {
    add <- dataset_ex %>%
      filter(!!!filter_ex)
  } else {
    add <- dataset_ex
  }
  add <- add %>%
      filter_extreme(order = exprs(EXSTDTC, EXSEQ),
                     by_vars = exprs(USUBJID),
                     mode = "first") %>%
      transmute(USUBJID, TRTSDTM = convert_dtc_to_dtm(impute_dtc(EXSTDTC)))

  left_join(dataset, add, by = c("USUBJID"))
}
