
#' Derive last dose date(-time)
#'
#' @param dataset Input AE dataset.
#' @param dataset_ex Input EX dataset.
#' @param filter_ex Filtering condition applied to EX dataset.
#' For example, it can be used to filter for valid dose.
#' @param by_vars Variables to join by.
#' @param dose_start The dose start date variable.
#' @param dose_end The dose end date variable.
#' @param analysis_date The analysis date variable.
#' @param output_var The output variable.
#' @param output_datetime Logical. Should only date or date-time variable be returned?
#' Defaults to `TRUE` (i.e. date-time variable).
#' @param check_dates_only Logical.
#' An assumption that start and end dates of treatment match is checked.
#' By default (`FALSE`), the date as well as the time component is checked.
#' If set to `TRUE`, then only the date component of those variables is checked.
#'
#' @details All date (date-time) variables need to be characters in standard ISO format.
#' See [`impute_dtc`], parameter `dtc` for further details.
#'
#' The last dose date is derived as follows:
#' Firstly, the `dataset_ex` is filtered using `filter_ex`, if provided.
#' This is useful for, for example, filtering for valid dose only.
#' Secondly, the datasets `dataset` and `dataset_ex` are joined using `by_vars`.
#' Thirdly, the last dose date is derived using provided variables.
#' The last dose date is the maximum date where `dose_end` is lower to or equal to
#' `analysis_date`, subject to both date values are non-NA.
#' Lastly, the last dose date is appended to the `dataset` and returned to the user.
#'
#' @return AE dataset with additional column `output_var`.
#'
#' @author Ondrej Slama
#'
#' @export
#'
#' @examples
#' data(ae); data(ex)
#' derive_last_dose(
#'   ae,
#'   ex,
#'   filter_ex = exprs(
#'     (EXDOSE > 0 | (EXDOSE == 0 & str_detect(EXTRT, "PLACEBO"))) & nchar(EXENDTC) >= 10
#'   ),
#'   dose_start = EXSTDTC,
#'   dose_end = EXENDTC,
#'   analysis_date = AESTDTC,
#'   output_var = LDOSEDTM,
#'   output_datetime = TRUE,
#'   check_dates_only = FALSE
#' )
derive_last_dose <- function(dataset,
                             dataset_ex,
                             filter_ex,
                             by_vars = exprs(STUDYID, USUBJID),
                             dose_start,
                             dose_end,
                             analysis_date,
                             output_var,
                             output_datetime = TRUE,
                             check_dates_only = FALSE) {

  stopifnot(rlang::is_scalar_logical(check_dates_only))

  # apply filtering condition
  if (!is.null(filter_ex)) {
    dataset_ex <- filter(dataset_ex, !!!filter_ex)
  }

  dose_start <- enquo(dose_start)
  dose_end <- enquo(dose_end)
  analysis_date <- enquo(analysis_date)
  output_var <- enquo(output_var)

  # assumption for last dose derivation: start and end dates (datetimes) need to match
  if (check_dates_only) {
    dataset_ex <- filter(dataset_ex, as.Date(!!dose_start) == as.Date(!!dose_end))
  } else {
    dataset_ex <- filter(dataset_ex, !!dose_start == !!dose_end)
  }

  # select only a subset of columns
  dataset_ex <- select(dataset_ex, !!!by_vars, !!dose_end)

  # calculate last dose date
  res <- dataset %>%
    mutate(DOMAIN = NULL) %>%
    inner_join(dataset_ex, by = map_chr(by_vars, as_string)) %>%
    group_by(!!!by_vars) %>%
    mutate(
      tmp_exendtc = impute_dtc(dtc = !!dose_end,
                               date_imputation = NULL,
                               time_imputation = "00:00:00") %>%
        convert_dtc_to_dtm(),
      tmp_aestdtc = impute_dtc(dtc = !!analysis_date,
                               date_imputation = NULL,
                               time_imputation = "23:59:59") %>%
        convert_dtc_to_dtm(),
      !!output_var := compute_ldosedtm(exendtc = .data$tmp_exendtc,
                                       aestdtc = .data$tmp_aestdtc)) %>%
    ungroup()

  # return either date or date-time variable
  if (!output_datetime) {
    res <- mutate(res, !!output_var := as.Date(!!output_var))
  }

  # return dataset with additional column
  left_join(dataset,
            dplyr::distinct(res, !!!by_vars, !!output_var),
            by = map_chr(by_vars, as_string))
}

#' Helper function to calculate last dose
#'
#' @param exendtc dose end date
#' @param aestdtc analysis date
#'
#' @return date-time vector
compute_ldosedtm <- function(exendtc, aestdtc) {
  if (any(!is.na(exendtc) & !is.na(aestdtc)) && any(exendtc <= aestdtc)) {
    max(exendtc[exendtc <= aestdtc])
  } else {
    as.POSIXct(NA)
  }
}
