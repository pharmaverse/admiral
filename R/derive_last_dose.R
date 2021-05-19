
#' Title
#'
#' @param exendtc
#' @param aestdtc
#' @param exdose
#'
#' @return
calc_ldosedtm <- function(exendtc, aestdtc, exdose) {
    ifelse(
      !is.na(exendtc) & !is.na(aestdtc) & exdose >= 0,
      `if`(any(exendtc <= aestdtc),
           as.character(max(exendtc[exendtc <= aestdtc])),
           NA_character_),
      NA_character_
  )
}

#' Title
#'
#' @param dataset
#' @param dataset_ex
#' @param filter_ex
#' @param by_vars
#'
#' @return
#' @export
#'
#' @examples
#' data(ae); data(ex)
#' derive_last_dose(ae, ex)
derive_last_dose <- function(dataset,
                             dataset_ex,
                             filter_ex = exprs((EXDOSE > 0 | (EXDOSE == 0 & str_detect(EXTRT, "PLACEBO"))) & nchar(EXENDTC) >= 10),
                             by_vars = exprs(STUDYID, USUBJID)) {

  if (!is.null(filter_ex)) {
    dataset_ex <- filter(dataset_ex, !!!filter_ex)
  }

  dataset_ex <- select(dataset_ex, !!!by_vars, .data$EXENDTC, .data$EXDOSE)

  out <- dataset %>%
    mutate(DOMAIN = NULL) %>%
    inner_join(dataset_ex, by = map_chr(by_vars, as_string)) %>%
    group_by(!!!by_vars) %>%
    mutate(
      LDOSEDTM = calc_ldosedtm(
        exendtc = convert_dtc_to_dtm(
          impute_dtc(
            dtc = .data$EXENDTC,
            date_imputation = NULL,
            time_imputation = "00:00:00"
          )
        ),
        aestdtc = convert_dtc_to_dtm(
          impute_dtc(
            dtc = .data$AESTDTC,
            date_imputation = NULL,
            time_imputation = "23:59:59"
          )
        ),
        exdose = .data$EXDOSE),
      LDOSEDT = as.character(as.Date(LDOSEDTM)),
      LDOSEDTM = as.character(LDOSEDTM)) %>%
    ungroup() %>%
    mutate(EXENDTC = NULL, EXDOSE = NULL)

  attr(out$LDOSEDTM, "label") <- "End Date/Time of Last Dose"
  attr(out$LDOSEDT, "label") <- "End Date of Last Dose"

  return(out)
}

