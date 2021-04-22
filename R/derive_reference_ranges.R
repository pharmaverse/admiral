#' Derive Reference Ranges and Range Indicator
#'
#' @param dataset The input dataset
#' @param meta_ref_ranges A dataset containing the references ranges. Required
#'        variables are `by_var`, `ANRLO`, `ANRHI`. `A1LO` and `A1HI` are optional.
#' @param by_var A `character` vector of variables by which to join `dataset` and
#'        `meta_ref_ranges`. By default, `"PARAMCD"`.
#'
#' @details
#' `ANRIND` is set to
#' - `"NORMAL"` if `AVAL` is greater or equal `ANRLO` and less than
#'   or equal `ANRHI`
#' - `"LOW"` if `AVAL` is less than `ANRLO` and either `A1LO` is missing or `AVAL`
#'   is greater than or equal `A1LO`
#' - `"HIGH"` if `AVAL` is greater than `ANRHI` and either `A1HI` is missing or `AVAL`
#'   is less than or equal `A1HI`
#' - `"LOW LOW"` if `AVAL` is less than `A1LO`
#' - `"HIGH HIGH"` if `AVAL` is greater than `A1HI`
#'
#' @return The input dataset with additional columns `ANRLO`, `ANRHI`, `ANRIND`
#'         and---if present in `meta_ref_ranges`---`A1LO` as well as `A1HI`
#' @export
#'
#' @examples
#' library(dplyr)
#' data(vs)
#' ref_ranges <- tibble::tribble(
#'   ~PARAMCD, ~ANRLO, ~ANRHI, ~A1LO, ~A1HI,
#'   "DIABP",  60,      80,    40,     90,
#'   "PULSE",  60,     100,    40,    110
#' )
#'
#' vs %>%
#'   mutate(
#'     PARAMCD = VSTESTCD,
#'     AVAL = VSSTRESN
#'   ) %>%
#'   filter(PARAMCD %in% c("PULSE", "DIABP")) %>%
#'   derive_reference_ranges(ref_ranges) %>%
#'   select(USUBJID, PARAMCD, AVAL, ANRLO:ANRIND)
derive_reference_ranges <- function(dataset,
                                    meta_ref_ranges,
                                    by_var = "PARAMCD") {
  assert_that(
    is.data.frame(dataset),
    is.data.frame(meta_ref_ranges),
    is.character(by_var)
  )
  assert_has_variables(dataset, by_var)
  assert_has_variables(meta_ref_ranges, c(by_var, "ANRLO", "ANRHI"))

  warn_if_ref_ranges_missing(dataset, meta_ref_ranges, by_var)

  has_a1lo <- "A1LO" %in% colnames(meta_ref_ranges)
  has_a1hi <- "A1HI" %in% colnames(meta_ref_ranges)
  if (!has_a1lo) meta_ref_ranges$A1LO <- NA
  if (!has_a1hi) meta_ref_ranges$A1HI <- NA

  result <- dataset %>%
    left_join(meta_ref_ranges, by = by_var) %>%
    mutate(
      ANRIND = case_when(
        AVAL >= ANRLO & AVAL <= ANRHI ~ "NORMAL",
        AVAL < ANRLO & (is.na(A1LO) | AVAL >= A1LO) ~ "LOW",
        AVAL > ANRHI & (is.na(A1HI) | AVAL <= A1HI) ~ "HIGH",
        AVAL < A1LO ~ "LOW LOW",
        AVAL > A1HI ~ "HIGH HIGH",
        TRUE ~ ""
      )
    )

  if (!has_a1lo) result$A1LO <- NULL
  if (!has_a1hi) result$A1HI <- NULL

  result
}
