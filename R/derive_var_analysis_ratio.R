#' Derive Ratio Variable
#'
#' Derives a ratio variable for a BDS dataset based on user specified variables.
#'
#' @param dataset Input dataset
#'
#' @param numer_var Variable containing numeric values to be used in the numerator of
#'  the ratio calculation.
#'
#' @param denom_var Variable containing numeric values to be used in the denominator of
#'  the ratio calculation.
#'
#' @param new_var A user-defined variable that will be appended to the dataset.
#' The default behavior will take the denominator variable and prefix it with `R2`
#' and append to the dataset. Using this argument will override this default behavior.
#'
#' Default is `NULL`.
#'
#' @details A user wishing to calculate a Ratio to Baseline, `AVAL / BASE` will
#' have returned a new variable `R2BASE` that will be appended to the input dataset.
#' Ratio to Analysis Range Lower Limit `AVAL / ANRLO` will return a new variable
#' `R2ANRLO`, and Ratio to Analysis Range  Upper Limit `AVAL / ANRHI` will return
#' a new variable `R2ANRLO`. Please note how the denominator variable has the prefix
#' `R2----`. A user can override the default returned variables by using the
#' `new_var` argument. Also, values of 0 in the denominator will return `NA` in
#' the derivation.
#'
#' Reference CDISC ADaM Implementation Guide
#' Version 1.1 Section 3.3.4 Analysis Parameter Variables for BDS Datasets
#'
#' @author Ben Straub
#'
#' @return The input dataset with a ratio variable appended
#'
#' @family der_bds_findings
#' @keywords der_bds_findings
#'
#' @export
#'
#' @examples
#' library(tibble)
#'
#' data <- tribble(
#'   ~USUBJID, ~PARAMCD, ~SEQ, ~AVAL, ~BASE, ~ANRLO, ~ANRHI,
#'   "P01", "ALT", 1, 27, 27, 6, 34,
#'   "P01", "ALT", 2, 41, 27, 6, 34,
#'   "P01", "ALT", 3, 17, 27, 6, 34,
#'   "P02", "ALB", 1, 38, 38, 33, 49,
#'   "P02", "ALB", 2, 39, 38, 33, 49,
#'   "P02", "ALB", 3, 37, 38, 33, 49
#' )
#'
#' # Returns "R2" prefixed variables
#' data %>%
#'   derive_var_analysis_ratio(numer_var = AVAL, denom_var = BASE) %>%
#'   derive_var_analysis_ratio(numer_var = AVAL, denom_var = ANRLO) %>%
#'   derive_var_analysis_ratio(numer_var = AVAL, denom_var = ANRHI)
#'
#' # Returns user-defined variables
#' data %>%
#'   derive_var_analysis_ratio(numer_var = AVAL, denom_var = BASE, new_var = R01BASE) %>%
#'   derive_var_analysis_ratio(numer_var = AVAL, denom_var = ANRLO, new_var = R01ANRLO) %>%
#'   derive_var_analysis_ratio(numer_var = AVAL, denom_var = ANRHI, new_var = R01ANRHI)
derive_var_analysis_ratio <- function(dataset,
                                      numer_var,
                                      denom_var,
                                      new_var = NULL) {
  numer_var <- assert_symbol(enquo(numer_var))
  denom_var <- assert_symbol(enquo(denom_var))

  assert_data_frame(dataset, required_vars = quo_c(numer_var, denom_var))
  new_var <- assert_symbol(enquo(new_var), optional = TRUE)

  if (quo_is_null(new_var)) {
    new_var <- sym(paste0("R2", rlang::as_name(denom_var)))
  }
  dataset <- dataset %>%
    mutate(
      !!new_var := if_else(!!denom_var == 0,
        NA_real_,
        !!numer_var / !!denom_var
      )
    )
  dataset
}
