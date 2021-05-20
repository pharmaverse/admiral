
#' Derive last dose date(-time)
#'
#' @param dataset Input AE dataset.
#' @param dataset_ex Input EX dataset.
#' @param filter_ex Filtering condition applied to EX dataset.
#' For example, it can be used to filter for valid dose.
#' @param by_vars Variables to join by.
#' @param check_dates_only Logical. An assumption that start and end dates of treatment match is checked.
#' By default (`TRUE`), the date as well as the time component is checked.
#' If set to `FALSE`, then only the date component of those variables is checked.
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
                               (EXDOSE > 0 | (EXDOSE == 0 & str_detect(EXTRT, "PLACEBO"))) & nchar(EXENDTC) >= 10), #nolint
                             by_vars = exprs(STUDYID, USUBJID),
                             check_dates_only = FALSE) {

  stopifnot(rlang::is_scalar_logical(check_dates_only))

  # apply filtering condition
  if (!is.null(filter_ex)) {
    dataset_ex <- filter(dataset_ex, !!!filter_ex)
  }

  # assumption for last dose derivation: start and end dates (datetimes) need to match
  if (check_dates_only) {
    dataset_ex <- filter(dataset_ex, as.Date(.data$EXSTDTC) == as.Date(.data$EXENDTC))
  } else {
    dataset_ex <- filter(dataset_ex, .data$EXSTDTC == .data$EXENDTC)
  }

  # select only a subset of columns
  dataset_ex <- select(dataset_ex, !!!by_vars, .data$EXENDTC, .data$EXDOSE)

  # calculate last dose date
  res <- dataset %>%
    mutate(DOMAIN = NULL) %>%
    inner_join(dataset_ex, by = map_chr(by_vars, as_string)) %>%
    group_by(!!!by_vars) %>%
    mutate(
      tmp_exendtc = impute_dtc(dtc = .data$EXENDTC,
                               date_imputation = NULL,
                               time_imputation = "00:00:00") %>%
        convert_dtc_to_dtm(),
      tmp_aestdtc = impute_dtc(dtc = .data$AESTDTC,
                               date_imputation = NULL,
                               time_imputation = "23:59:59") %>%
        convert_dtc_to_dtm(),
      LDOSEDTM = compute_ldosedtm(exendtc = .data$tmp_exendtc,
                                  aestdtc = .data$tmp_aestdtc,
                                  exdose = .data$EXDOSE)) %>%
    ungroup() %>%
    mutate(LDOSEDT = as.Date(.data$LDOSEDTM))

  left_join(dataset,
            dplyr::distinct(res, !!!by_vars, .data$LDOSEDTM, .data$LDOSEDT),
            by = map_chr(by_vars, as_string))
}

#' Helper function to calculate last dose
#'
#' @param exendtc EX.EXENDTC
#' @param aestdtc ADAE.AESTDTC
#' @param exdose EX.EXDOSE
#'
#' @return date-time vector
compute_ldosedtm <- function(exendtc, aestdtc, exdose) {
  if (any(!is.na(exendtc) & !is.na(aestdtc) & exdose >= 0) && any(exendtc <= aestdtc)) {
    max(exendtc[exendtc <= aestdtc])
  } else {
    as.POSIXct(NA)
  }
}
