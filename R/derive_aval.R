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
#'   select(USUBJID, VSTESTCD, VSSTRESC, AVALC, VSSTRESN, AVAL, VSSTRESU, AVALU)
derive_aval <- function(dataset) {
  select_col <- function(pattern) {
    full_pattern <- paste0("^[A-Z]{2}", pattern, "$")
    col <- grep(full_pattern, colnames(dataset), value = TRUE)
    if (length(col) == 0L) {
      msg <- paste0("The input dataset doesn't contain a `--", pattern, "` variable")
      abort(msg)
    } else {
      col
    }
  }

  stresc <- select_col("STRESC")
  stresn <- select_col("STRESN")
  stresu <- select_col("STRESU")

  mutate(
    dataset,
    AVALC = !!sym(stresc),
    AVAL = !!sym(stresn),
    AVALU = !!sym(stresu)
  )
}
