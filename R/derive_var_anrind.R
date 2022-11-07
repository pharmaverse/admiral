#' Derive Reference Range Indicator
#'
#' @param dataset The input dataset
#'
#' @details
#' `ANRIND` is set to
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
#' @return The input dataset with additional column `ANRIND`
#'
#' @author Thomas Neitmann
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
#' ref_ranges <- tribble(
#'   ~PARAMCD, ~ANRLO, ~ANRHI, ~A1LO, ~A1HI,
#'   "DIABP",      60,     80,    40,    90,
#'   "PULSE",      60,    100,    40,   110
#' )
#'
#' admiral_vs %>%
#'   mutate(
#'     PARAMCD = VSTESTCD,
#'     AVAL = VSSTRESN
#'   ) %>%
#'   filter(PARAMCD %in% c("PULSE", "DIABP")) %>%
#'   derive_vars_merged(ref_ranges, by_vars = vars(PARAMCD)) %>%
#'   derive_var_anrind() %>%
#'   select(USUBJID, PARAMCD, AVAL, ANRLO:ANRIND)
derive_var_anrind <- function(dataset) {
  assert_data_frame(dataset, required_vars = vars(ANRLO, ANRHI, AVAL))

  # Temporarily add these variables to the dataset if they are not included
  has_a1lo <- "A1LO" %in% colnames(dataset)
  has_a1hi <- "A1HI" %in% colnames(dataset)
  if (!has_a1lo) dataset[["A1LO"]] <- NA_character_
  if (!has_a1hi) dataset[["A1HI"]] <- NA_character_

  result <- dataset %>%
    mutate(
      ANRIND = case_when(
        AVAL >= ANRLO & is.na(ANRHI) ~ "NORMAL",
        AVAL <= ANRHI & is.na(ANRLO) ~ "NORMAL",
        AVAL >= ANRLO & AVAL <= ANRHI ~ "NORMAL",
        AVAL < ANRLO & (is.na(A1LO) | AVAL >= A1LO) ~ "LOW",
        AVAL > ANRHI & (is.na(A1HI) | AVAL <= A1HI) ~ "HIGH",
        AVAL < A1LO ~ "LOW LOW",
        AVAL > A1HI ~ "HIGH HIGH",
        TRUE ~ NA_character_
      )
    )

  # Remove the variables if they have been added above
  if (!has_a1lo) result[["A1LO"]] <- NULL
  if (!has_a1hi) result[["A1HI"]] <- NULL

  result
}
