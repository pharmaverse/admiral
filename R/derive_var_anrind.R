#' Derive Reference Range Indicator
#'
#' @param dataset The input dataset
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
                              use_a1hia1lo = FALSE) {
  if (use_a1hia1lo) {
    assert_data_frame(dataset, required_vars = exprs(ANRLO, ANRHI, A1HI, A1LO, AVAL))

    low_cond <- "AVAL < ANRLO & (is.na(A1LO) | AVAL >= A1LO)"
    high_cond <- "AVAL > ANRHI & (is.na(A1HI) | AVAL <= A1HI)"
    lowlow_cond <- "AVAL < A1LO"
    highhigh_cond <- "AVAL > A1HI"
  } else {
    assert_data_frame(dataset, required_vars = exprs(ANRLO, ANRHI, AVAL))

    low_cond <- "AVAL < ANRLO"
    high_cond <- "AVAL > ANRHI"
    lowlow_cond <- "FALSE"
    highhigh_cond <- "FALSE"
  }

  result <- dataset %>%
    mutate(
      ANRIND = case_when(
        AVAL >= ANRLO & is.na(ANRHI) ~ "NORMAL",
        AVAL <= ANRHI & is.na(ANRLO) ~ "NORMAL",
        AVAL >= ANRLO & AVAL <= ANRHI ~ "NORMAL",
        eval(parse(text = low_cond)) ~ "LOW",
        eval(parse(text = high_cond)) ~ "HIGH",
        eval(parse(text = lowlow_cond)) ~ "LOW LOW",
        eval(parse(text = highhigh_cond)) ~ "HIGH HIGH",
        TRUE ~ NA_character_
      )
    )

  result
}
