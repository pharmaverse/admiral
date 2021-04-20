#' Derive Cause of Death Category
#'
#' Derives cause of death category (`DTHCAT`)
#'
#' @param dataset Input dataset
#'
#'   The columns specified by the `start_date` and the `end_date` parameter are
#'   expected.
#'
#' @details The total treatment duration is derived as the number of days from start to
#'   end date plus one.
#'
#' @keywords Roche
#'
#' @author Shimeng Huang
#'
#' @return The input dataset with `DTHCAT` added
#'
#' @export
#'
#' @examples
#' data <- tibble::tribble(
#'   ~DTHDOM, ~DTHCAUS,
#'   "ADVERSE EVENT", ""
#' )
#'
#' derive_var_dthcat(data)
derive_var_dthcat <- function(dataset) {
  assert_has_variables(dataset, c("DTHDOM", "DTHCAUS"))
  dataset <- dataset %>%
    mutate(DTHCAT = case_when(
      DTHDOM == "" ~ "",
      DTHDOM == "ADVERSE EVENT" ~ "AE",
      stringr::str_detect(DTHCAUS, "PROGRESSIVE DISEASE") |
        stringr::str_detect(DTHCAUS, "DISEASE RELAPSE") ~ "PROGRESSIVE DISEASE",
      TRUE ~ "OTHER"
    ))
  dataset
}
