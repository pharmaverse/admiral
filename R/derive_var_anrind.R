#' Derive Reference Range Indicator
#'
#' @param dataset
#'   `r roxygen_param_dataset()`
#'   `ANRLO`, `ANRHI`, and `AVAL` are expected and if `use_a1hia1lo` is set to `TRUE`,
#'   `A1LO` and `A1H1` are expected as well.
#'
#' @param signif_dig Number of significant digits to use when comparing values.
#'
#'   Significant digits used to avoid floating point discrepancies when comparing numeric values.
#'   See blog: [How admiral handles floating
#'   points](https://pharmaverse.github.io/blog/posts/2023-10-30_floating_point/floating_point.html)
#'
#' @param use_a1hia1lo A logical value indicating whether to use `A1H1` and `A1LO` in
#' the derivation of `ANRIND`.
#'
#' @details
#' In the case that `A1H1` and `A1LO` are to be used, `ANRIND` is set to:
#' - `"NORMAL"` if `AVAL` is greater or equal `ANRLO` and less than
#'   or equal `ANRHI`; or if `AVAL` is greater than or equal `ANRLO` and `ANRHI`
#'   is missing; or if `AVAL` is less than or equal `ANRHI` and `ANRLO` is
#'   missing
#' - `"LOW"` if `AVAL` is less than `ANRLO` and either `A1LO` is missing or `AVAL`
#'   is greater than or equal `A1LO`
#' - `"HIGH"` if `AVAL` is greater than `ANRHI` and either `A1HI` is missing or `AVAL`
#'   is less than or equal `A1HI`
#' - `"LOW LOW"` if `AVAL` is less than `A1LO`
#' - `"HIGH HIGH"` if `AVAL` is greater than `A1HI`
#'
#' In the case that `A1H1` and `A1LO` are not to be used, `ANRIND` is set to:
#' - `"NORMAL"` if `AVAL` is greater or equal `ANRLO` and less than
#'   or equal `ANRHI`; or if `AVAL` is greater than or equal `ANRLO` and `ANRHI`
#'   is missing; or if `AVAL` is less than or equal `ANRHI` and `ANRLO` is
#'   missing
#' - `"LOW"` if `AVAL` is less than `ANRLO`
#' - `"HIGH"` if `AVAL` is greater than `ANRHI`
#'
#' @return The input dataset with additional column `ANRIND`
#'
#' @family der_bds_findings
#' @keywords der_bds_findings
#'
#' @export
#'
#' @examples
#' library(tibble)
#' library(dplyr, warn.conflicts = FALSE)
#'
#' vs <- tibble::tribble(
#'   ~USUBJID, ~PARAMCD, ~AVAL, ~ANRLO, ~ANRHI, ~A1LO, ~A1HI,
#'   "P01",       "PUL",    70,     60,    100,    40,   110,
#'   "P01",       "PUL",    57,     60,    100,    40,   110,
#'   "P01",       "PUL",    60,     60,    100,    40,   110,
#'   "P01",     "DIABP",   102,     60,     80,    40,    90,
#'   "P02",       "PUL",   109,     60,    100,    40,   110,
#'   "P02",       "PUL",   100,     60,    100,    40,   110,
#'   "P02",     "DIABP",    80,     60,     80,    40,    90,
#'   "P03",       "PUL",    39,     60,    100,    40,   110,
#'   "P03",       "PUL",    40,     60,    100,    40,   110
#' )
#'
#' vs %>% derive_var_anrind(use_a1hia1lo = TRUE)
#' vs %>% derive_var_anrind(use_a1hia1lo = FALSE)
#'
derive_var_anrind <- function(dataset,
                              signif_dig = get_admiral_option("signif_digits"),
                              use_a1hia1lo = FALSE) {
  # check input parameter holding significant digits has correct value
  assert_integer_scalar(signif_dig, subset = "positive")

  aval <- "signif(AVAL, signif_dig)"
  anrlo <- "signif(ANRLO, signif_dig)"
  anrhi <- "signif(ANRHI, signif_dig)"
  a1lo <- "signif(A1LO, signif_dig)"
  a1hi <- "signif(A1HI, signif_dig)"

  if (use_a1hia1lo) {
    assert_data_frame(dataset, required_vars = exprs(ANRLO, ANRHI, A1HI, A1LO, AVAL))

    low_cond <- paste(aval, "<", anrlo, "& (is.na(A1LO) |", aval, ">=", a1lo, ")")
    high_cond <- paste(aval, ">", anrhi, "& (is.na(A1HI) |", aval, "<=", a1hi, ")")
    lowlow_cond <- paste(aval, "<", a1lo)
    highhigh_cond <- paste(aval, ">", a1hi)
  } else {
    assert_data_frame(dataset, required_vars = exprs(ANRLO, ANRHI, AVAL))

    low_cond <- paste(aval, "<", anrlo)
    high_cond <- paste(aval, ">", anrhi)
    lowlow_cond <- "FALSE"
    highhigh_cond <- "FALSE"
  }

  result <- dataset %>%
    mutate(
      ANRIND = case_when(
        eval(parse(text = paste(aval, ">=", anrlo, "& is.na(ANRHI)"))) ~ "NORMAL",
        eval(parse(text = paste(aval, "<=", anrhi, "& is.na(ANRLO)"))) ~ "NORMAL",
        eval(parse(text = paste(aval, ">=", anrlo, "&", aval, "<=", anrhi))) ~ "NORMAL",
        eval(parse(text = low_cond)) ~ "LOW",
        eval(parse(text = high_cond)) ~ "HIGH",
        eval(parse(text = lowlow_cond)) ~ "LOW LOW",
        eval(parse(text = highhigh_cond)) ~ "HIGH HIGH",
        TRUE ~ NA_character_
      )
    )

  result
}
