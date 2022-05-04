#' Join Supplementary Qualifier Variables into the Parent SDTM Domain
#'
#'#' @description
#' `r lifecycle::badge("deprecated")`
#'
#' *Deprecated*, please use `metatools::combine_supp()` instead.
#'
#' The SDTM does not allow any new variables beside ones assigned to each SDTM
#' domain. So, Supplemental Qualifier is introduced to supplement each SDTM
#' domain to contain non standard variables. `dataset_suppqual` can be either
#' a single SUPPQUAL dataset or separate supplementary data sets (SUPP) such
#' as SUPPDM, SUPPAE, and SUPPEX. When a `dataset_suppqual` is a single SUPPQUAL
#' dataset, specify two character`domain` value.
#'
#' `derive_vars_suppqual()` expects `USUBJID`, `RDOMAIN`, `IDVAR`, `IDVARVAL`,
#'  `QNAM`, `QLABEL`, and `QVAL` variables to exist in `dataset_suppqual`.
#'
#' @param dataset A SDTM domain data set.
#'
#' @param dataset_suppqual A Supplemental Qualifier (SUPPQUAL) data set.
#'
#'
#' @param domain Two letter domain value. Used when supplemental data set is
#'   common across multiple SDTM domain.
#'
#' @return A data frame with SUPPQUAL variables appended to parent data set.
#'
#' @author Vignesh Thanikachalam
#'
#' @export
#'

derive_vars_suppqual <- function(dataset, dataset_suppqual, domain = NULL) {
  deprecate_stop("0.7.0", "derive_vars_suppqual()", "metatools::combine_supp()")
}
