#' Derive `AVALC`, `AVAL` and `AVALU`
#'
#' @param dataset The input dataset
#'
#' @author Thomas Neitmann
#'
#' @details
#' The SDTM variables `--STRES`, `--STRESN` and `--STRESU` are mapped to the
#' ADaM variables `AVALC`, `AVAL` and `AVALU`, respectively.
#'
#' @export
#'
#' @examples
#' library(dplyr)
#' data(vs)
#' vs %>%
#'   derive_aval() %>%
#'   select(USUBJID, PARAMCD, VSSTRES, AVALC, VSSTRESN, AVAL, VSSTRESU, AVALU)
derive_aval <- function(dataset) {
  select_col <- function(pattern) grep(pattern, colnames(dataset), value = TRUE)

  stres <- select_col("^[A-Z]{2}STRES$")
  stresn <- select_col("^[A-Z]{2}STRESN$")
  stresu <- select_col("^[A-Z]{2}STRESU$")

  mutate(
    dataset,
    AVALC = !!sym(stres),
    AVAL = !!sym(stresn),
    AVALU = !!sym(stresu)
  )
}
