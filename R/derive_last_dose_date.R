#' Derive Last Dose Date(-time)
#'
#' @param dataset Input dataset.
#' @param dataset_ex Input EX dataset.
#' @param filter_ex Filtering condition applied to EX dataset.
#' For example, it can be used to filter for valid dose.
#' Defaults to NULL.
#' @param by_vars Variables to join by (created by `dplyr::vars`).
#' @param dose_start The dose start date variable.
#' @param dose_end The dose end date variable.
#' @param analysis_date The analysis date variable.
#' @param dataset_seq_var The sequence variable
#' (this together with `by_vars` creates the keys of `dataset`).
#' @param new_var The output variable.
#' @param output_datetime Logical. Should only date or date-time variable be returned?
#' Defaults to `TRUE` (i.e. date-time variable).
#' @param check_dates_only Logical.
#' An assumption that start and end dates of treatment match is checked.
#' By default (`FALSE`), the date as well as the time component is checked.
#' If set to `TRUE`, then only the date component of those variables is checked.
#' @param traceability_vars A named list returned by [`vars()`] listing the traceability variables,
#' e.g. `vars(LDOSEDOM = "EX", LDOSESEQ = EXSEQ)`.
#' The left-hand side (names of the list elements) gives the names of the traceability variables
#' in the returned dataset.
#' The right-hand side (values of the list elements) gives the values of the traceability variables
#' in the returned dataset.
#' These can be either strings or symbols referring to existing variables.
#'
#' @details All date (date-time) variables can be characters in standard ISO format or
#' of date / date-time class.
#' For ISO format, see [`impute_dtc`] - parameter `dtc` for further details.
#' Date-time imputations are done as follows:
#' * `dose_end`: no date imputation, time imputation to `00:00:00` if time is missing.
#' * `analysis_date`: no date imputation, time imputation to `23:59:59` if time is missing.
#'
#' The last dose date is derived as follows:
#' 1. The `dataset_ex` is filtered using `filter_ex`, if provided.
#' This is useful for, for example, filtering for valid dose only.
#' 2. The datasets `dataset` and `dataset_ex` are joined using `by_vars`.
#' 3. The last dose date is derived:
#' the last dose date is the maximum date where `dose_end` is lower to or equal to
#' `analysis_date`, subject to both date values are non-NA.
#' The last dose date is derived per `by_vars` and `dataset_seq_var`.
#' 4. The last dose date is appended to the `dataset` and returned to the user.
#'
#' Furthermore, the following assumption is checked: start and end dates (datetimes) need to match.
#' Use `check_dates_only` to control whether only dates or whole date-times need to be equal.
#'
#' @return Input dataset with additional column `new_var`.
#'
#' @author Ondrej Slama
#'
#' @keywords adae derivation
#'
#' @export
#'
#' @examples
#' library(dplyr, warn.conflicts = FALSE)
#' data(ae)
#' data(ex_single)
#'
#' ae %>%
#'   head(100) %>%
#'   derive_last_dose(
#'     head(ex_single, 100),
#'     filter_ex = (EXDOSE > 0 | (EXDOSE == 0 & grepl("PLACEBO", EXTRT))) &
#'       nchar(EXENDTC) >= 10,
#'     dose_start = EXSTDTC,
#'     dose_end = EXENDTC,
#'     analysis_date = AESTDTC,
#'     dataset_seq_var = AESEQ,
#'     new_var = LDOSEDTM,
#'     output_datetime = TRUE,
#'     check_dates_only = FALSE
#'   ) %>%
#'   select(STUDYID, USUBJID, AESEQ, AESTDTC, LDOSEDTM)

derive_last_dose_date <- function(dataset,
                                  dataset_ex,
                                  filter_ex,
                                  by_vars ,
                                  dose_start,
                                  dose_end,
                                  new_var,
                                  analysis_date,
                                  dataset_seq_var,
                                  check_dates_only,
                                  traceability_vars ){

  filter_ex <- assert_filter_cond(enquo(filter_ex), optional = TRUE)
  by_vars <- assert_vars(by_vars)
  dose_start <- assert_symbol(enquo(dose_start))
  dose_end <- assert_symbol(enquo(dose_end))
  analysis_date <- assert_symbol(enquo(analysis_date))
  dataset_seq_var <- assert_symbol(enquo(dataset_seq_var))
  new_var <- assert_symbol(enquo(new_var))
  assert_logical_scalar(check_dates_only)
  stopifnot(is_quosures(traceability_vars) | is.null(traceability_vars))
  assert_data_frame(dataset, quo_c(by_vars, analysis_date, dataset_seq_var))
  assert_data_frame(dataset_ex, quo_c(by_vars, dose_start, dose_end))

  res <- get_last_dose(dataset = dataset,
                dataset_ex = dataset_ex,
                filter_ex = !!filter_ex,
                by_vars = by_vars ,
                dose_start = !!dose_start,
                dose_end  = !!dose_end,
                analysis_date = !!analysis_date,
                dataset_seq_var = !!dataset_seq_var,
                check_dates_only = !!check_dates_only,
                traceability_vars = traceability_vars) %>%
  select(colnames(dataset), !!dose_start)

  if (!is.null(new_var)){res %>% rename(!!new_var := !!dose_start)}

  else {res}

}

