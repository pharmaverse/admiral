
#' Helper function to calculate last dose
#'
#' @param exendtc EX.EXENDTC
#' @param aestdtc ADAE.AESTDTC
#' @param exdose EX.EXDOSE
#'
#' @return date-time vector
calc_ldosedtm <- function(exendtc, aestdtc, exdose) {
    if (any(!is.na(exendtc) & !is.na(aestdtc) & exdose >= 0) && any(exendtc <= aestdtc)) {
      max(exendtc[exendtc <= aestdtc])
    } else {
      as.POSIXct(NA)
    }
}

#' Derive last dose date(-time)
#'
#' @param dataset Input AE dataset.
#' @param dataset_ex Input EX dataset.
#' @param filter_ex Filtering condition applied to EX dataset.
#' For example, it can be used to filter for valid dose.
#' @param by_vars Variables to join by.
#'
#' @return AE dataset with additional columns `LDOSEDTM` and `LDOSEDT`.
#'
#' @author Ondrej Slama
#'
#' @export
#'
#' @examples
#' data(ae); data(ex)
#' derive_last_dose(ae, ex)
derive_last_dose <- function(dataset,
                             dataset_ex,
                             filter_ex = exprs(
                               (EXDOSE > 0 | (EXDOSE == 0 & str_detect(EXTRT, "PLACEBO"))) &  #nolint
                                 nchar(EXENDTC) >= 10),
                             by_vars = exprs(STUDYID, USUBJID)) {

  if (!is.null(filter_ex)) {
    dataset_ex <- filter(dataset_ex, !!!filter_ex)
  }

  dataset_ex <- select(dataset_ex, !!!by_vars, .data$EXENDTC, .data$EXDOSE)

  res <- dataset %>%
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
        exdose = .data$EXDOSE)) %>%
    ungroup() %>%
    mutate(
      LDOSEDT = format(.data$LDOSEDTM, "%Y-%m-%d"),
      LDOSEDTM = format(.data$LDOSEDTM, "%Y-%m-%dT%H:%M:%S"),
      EXENDTC = NULL,
      EXDOSE = NULL
    )

  attr(res$LDOSEDTM, "label") <- "End Date/Time of Last Dose"
  attr(res$LDOSEDT, "label") <- "End Date of Last Dose"

  out <- left_join(dataset,
                   dplyr::distinct(res, !!!by_vars, .data$LDOSEDTM, .data$LDOSEDT),
                   by = map_chr(by_vars, as_string))

  return(out)
}
