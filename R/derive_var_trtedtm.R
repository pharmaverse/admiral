#' Derive datetime of Last Exposure to Treatment (TRTEDTM)
#'
#' Derives datetime of Last Exposure to Treatment (TRTEDTM)
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
#' @return The input dataset with `TRTEDTM` variable added
#'
#' @keywords adsl time derivation
#'
#' @export
#'
#' @examples
#' data("ex")
#' data("dm")
#'
#' derive_var_trtedtm(dm,
#'                    dataset_ex = ex)
#'

derive_var_trtedtm <- function(
  dataset,
  dataset_ex,
  filter_ex = exprs((EXDOSE > 0 | (EXDOSE == 0 & str_detect(EXTRT, "PLACEBO"))) & nchar(EXENDTC) >= 10)){

  derive_merged_vars(
    dataset,
    dataset_add = dataset_ex,
    filter_add = filter_ex,
    new_vars = exprs(TRTEDTM := dtc_dtm(impute_dtc(EXENDTC, time_imputation = 'LAST'))),
    filter_first_order = exprs(desc(EXENDTC), desc(EXSEQ)))
}
