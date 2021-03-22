#' Derive datetime of First Exposure to Treatment (TRTSDTM)
#'
#' Derives datetime of First Exposure to Treatment (TRTSDTM)
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
#' @examples
#' data("ex")
#' data("dm")
#'
#' derive_var_trtsdtm(dm)
#'

derive_var_trtsdtm <- function(
  dataset,
  dataset_ex,
  filter_ex = exprs(EXDOSE > 0 | (EXDOSE == 0 & str_detect(EXTRT, "PLACEBO")))){

  derive_merged_vars(dataset,
                     dataset_add = dataset_ex,
                     filter_add = filter_ex,
                     new_vars = exprs(TRTSDTM := ymd(EXSTDTC)),
                     filter_first_order = exprs(EXSTDTC, EXSEQ))
}
