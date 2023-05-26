#' Derive Reference Range Indicator
#'
#' @param dataset The input dataset
#' @param use_a1hia1lo Boolean value indicating whether to use `A1H1` and `A1LO` in
#' the derivation of `ANRIND`
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
#' library(admiral.test)
#' data(admiral_vs)
#'
#' ref_ranges1 <- tribble(
#'   ~PARAMCD, ~ANRLO, ~ANRHI, ~A1LO, ~A1HI,
#'   "DIABP",      60,     80,    40,    90,
#'   "PULSE",      60,    100,    40,   110
#' )
#'
#' ref_ranges2 <- tribble(
#'   ~PARAMCD, ~ANRLO, ~ANRHI,
#'   "DIABP",      60,     80,
#'   "PULSE",      60,    100
#' )
#'
#' vs <- admiral_vs %>%
#'   mutate(
#'     PARAMCD = VSTESTCD,
#'     AVAL = VSSTRESN
#'   ) %>%
#'   filter(PARAMCD %in% c("PULSE", "DIABP"))
#'
#' vs1 <- vs %>%
#'   derive_vars_merged(ref_ranges1, by_vars = exprs(PARAMCD))
#'
#' vs2 <- vs %>%
#'   derive_vars_merged(ref_ranges2, by_vars = exprs(PARAMCD))
#'
#' vs_anrind1 <- vs1
#' derive_var_anrind(use_a1hia1lo = TRUE) %>%
#'   select(USUBJID, PARAMCD, AVAL, ANRLO:ANRIND)
#'
#' vs_anrind2 <- vs2
#' derive_var_anrind(use_a1hia1lo = FALSE) %>%
#'   select(USUBJID, PARAMCD, AVAL, ANRLO:ANRIND)
#'
#' vs_anrind1
#' vs_anrind2
derive_var_anrind <- function(dataset,
                              use_a1hia1lo = TRUE) {
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
