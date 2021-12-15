#' Derive Last Dose Date(-time)
#'
#' *Deprecated*, please use `derive_var_last_dose()` instead.
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
#' library(admiral.test)
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
#'
#' # or with traceability variables
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
#'     check_dates_only = FALSE,
#'     traceability_vars = dplyr::vars(LDOSEDOM = "EX", LDOSESEQ = EXSEQ, LDOSEVAR = "EXSTDTC")
#'   ) %>%
#'   select(STUDYID, USUBJID, AESEQ, AESTDTC, LDOSEDTM, LDOSEDOM, LDOSESEQ, LDOSEVAR)
#'
derive_last_dose <- function(dataset,
                             dataset_ex,
                             filter_ex = NULL,
                             by_vars = vars(STUDYID, USUBJID),
                             dose_start,
                             dose_end,
                             analysis_date,
                             dataset_seq_var,
                             new_var,
                             output_datetime = TRUE,
                             check_dates_only = FALSE,
                             traceability_vars = NULL) {
  deprecate_warn("0.6.0", "derive_last_dose()", "derive_var_last_dose()")
  filter_ex <- assert_filter_cond(enquo(filter_ex), optional = TRUE)
  by_vars <- assert_vars(by_vars)
  dose_start <- assert_symbol(enquo(dose_start))
  dose_end <- assert_symbol(enquo(dose_end))
  analysis_date <- assert_symbol(enquo(analysis_date))
  dataset_seq_var <- assert_symbol(enquo(dataset_seq_var))
  new_var <- assert_symbol(enquo(new_var))
  assert_logical_scalar(output_datetime)
  assert_logical_scalar(check_dates_only)
  stopifnot(is_quosures(traceability_vars) | is.null(traceability_vars))
  assert_data_frame(dataset, quo_c(by_vars, analysis_date, dataset_seq_var))
  assert_data_frame(dataset_ex, quo_c(by_vars, dose_start, dose_end))

  # by_vars converted to string
  by_vars_str <- vars2chr(by_vars)

  # assumption for last dose derivation: start and end dates (datetimes) need to match
  if (check_dates_only) {
    check_cond <- summarise(dataset_ex,
                            all_equal = all(as.Date(!!dose_start) == as.Date(!!dose_end)))
  } else {
    check_cond <- summarise(dataset_ex,
                            all_equal = all(!!dose_start == !!dose_end))
  }
  if (!check_cond$all_equal) {
    stop(paste(
      "Not all values of", as_string(quo_get_expr(dose_start)),
      "are equal to", as_string(quo_get_expr(dose_end))
    ))
  }

  # run traceability if requested
  if (!is.null(traceability_vars)) {
    trace_vars_str <- names(traceability_vars)
    dataset_ex <- mutate(dataset_ex, !!!traceability_vars)
  } else {
    trace_vars_str <- character(0)
  }

  # select only a subset of columns
  dataset_ex <- dataset_ex %>%
    filter_if(filter_ex) %>%
    select(!!!by_vars, !!dose_end, trace_vars_str)

  # calculate the last dose date
  res <- dataset %>%
    mutate(DOMAIN = NULL) %>%
    inner_join(dataset_ex, by = by_vars_str) %>%
    mutate_at(vars(!!dose_end, !!analysis_date),
              ~ `if`(is_date(.), convert_dtm_to_dtc(.), .)) %>%
    mutate(
      tmp_dose_end_date = convert_dtc_to_dtm(
        dtc = !!dose_end,
        date_imputation = NULL,
        time_imputation = "00:00:00"
      ),
      tmp_analysis_date = convert_dtc_to_dtm(
        dtc = !!analysis_date,
        date_imputation = NULL,
        time_imputation = "23:59:59"
      )
    ) %>%
    group_by(!!!by_vars, !!dataset_seq_var)

  # if no traceability variables are required, simply calculate the last dose date
  if (is.null(traceability_vars)) {
    res <- res %>%
      summarise(
        ldose_idx = compute_ldose_idx(dose_end = .data$tmp_dose_end_date,
                                      analysis_date = .data$tmp_analysis_date),
        !!new_var := as.POSIXct(as.character(.data$tmp_dose_end_date[.data$ldose_idx]),
                                tz = lubridate::tz(.data$tmp_dose_end_date))
      ) %>%
      ungroup()
  } else {
    # calculate the last dose date and get the appropriate traceability variables
    res <- res %>%
      mutate(
        ldose_idx = compute_ldose_idx(dose_end = .data$tmp_dose_end_date,
                                      analysis_date = .data$tmp_analysis_date),
        !!new_var := as.POSIXct(as.character(.data$tmp_dose_end_date[.data$ldose_idx]),
                                tz = lubridate::tz(.data$tmp_dose_end_date))
      ) %>%
      mutate_at(trace_vars_str, list(~ .[.data$ldose_idx])) %>%
      distinct(!!new_var, !!!syms(trace_vars_str)) %>%
      ungroup()
  }

  # return either date or date-time variable
  if (!output_datetime) {
    res <- mutate(res, !!new_var := as.Date(!!new_var))
  }

  # return dataset with additional column
  left_join(dataset,
            distinct(res,
                     !!!by_vars, !!dataset_seq_var, !!new_var,
                     !!!syms(trace_vars_str)),
            by = c(by_vars_str, as_string(quo_get_expr(dataset_seq_var))))
}


#' Helper function to get the index of last dose date
#'
#' @param dose_end dose end date
#' @param analysis_date analysis date
#'
#' @noRd
#'
#' @return index. The last dose date is then `dose_end[return_value]`
compute_ldose_idx <- function(dose_end, analysis_date) {
  if (any(!is.na(dose_end) & !is.na(analysis_date)) && any(dose_end <= analysis_date)) {
    which.max(dose_end[dose_end <= analysis_date])
  } else {
    NA_integer_
  }
}

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
#' library(admiral.test)
#' data(ae)
#' data(ex_single)
#'
#' ae %>%
#'   head(100) %>%
#'   derive_var_last_dose(
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
#'
#' # or with traceability variables
#' ae %>%
#'   head(100) %>%
#'   derive_var_last_dose(
#'     head(ex_single, 100),
#'     filter_ex = (EXDOSE > 0 | (EXDOSE == 0 & grepl("PLACEBO", EXTRT))) &
#'       nchar(EXENDTC) >= 10,
#'     dose_start = EXSTDTC,
#'     dose_end = EXENDTC,
#'     analysis_date = AESTDTC,
#'     dataset_seq_var = AESEQ,
#'     new_var = LDOSEDTM,
#'     output_datetime = TRUE,
#'     check_dates_only = FALSE,
#'     traceability_vars = dplyr::vars(LDOSEDOM = "EX", LDOSESEQ = EXSEQ, LDOSEVAR = "EXSTDTC")
#'   ) %>%
#'   select(STUDYID, USUBJID, AESEQ, AESTDTC, LDOSEDTM, LDOSEDOM, LDOSESEQ, LDOSEVAR)
#'
derive_var_last_dose <- function(dataset,
                             dataset_ex,
                             filter_ex = NULL,
                             by_vars = vars(STUDYID, USUBJID),
                             dose_start,
                             dose_end,
                             analysis_date,
                             dataset_seq_var,
                             new_var,
                             output_datetime = TRUE,
                             check_dates_only = FALSE,
                             traceability_vars = NULL) {
  filter_ex <- assert_filter_cond(enquo(filter_ex), optional = TRUE)
  by_vars <- assert_vars(by_vars)
  dose_start <- assert_symbol(enquo(dose_start))
  dose_end <- assert_symbol(enquo(dose_end))
  analysis_date <- assert_symbol(enquo(analysis_date))
  dataset_seq_var <- assert_symbol(enquo(dataset_seq_var))
  new_var <- assert_symbol(enquo(new_var))
  assert_logical_scalar(output_datetime)
  assert_logical_scalar(check_dates_only)
  stopifnot(is_quosures(traceability_vars) | is.null(traceability_vars))
  assert_data_frame(dataset, quo_c(by_vars, analysis_date, dataset_seq_var))
  assert_data_frame(dataset_ex, quo_c(by_vars, dose_start, dose_end))

  # by_vars converted to string
  by_vars_str <- vars2chr(by_vars)

  # assumption for last dose derivation: start and end dates (datetimes) need to match
  if (check_dates_only) {
    check_cond <- summarise(dataset_ex,
                            all_equal = all(as.Date(!!dose_start) == as.Date(!!dose_end)))
  } else {
    check_cond <- summarise(dataset_ex,
                            all_equal = all(!!dose_start == !!dose_end))
  }
  if (!check_cond$all_equal) {
    stop(paste(
      "Not all values of", as_string(quo_get_expr(dose_start)),
      "are equal to", as_string(quo_get_expr(dose_end))
    ))
  }

  # run traceability if requested
  if (!is.null(traceability_vars)) {
    trace_vars_str <- names(traceability_vars)
    dataset_ex <- mutate(dataset_ex, !!!traceability_vars)
  } else {
    trace_vars_str <- character(0)
  }

  # select only a subset of columns
  dataset_ex <- dataset_ex %>%
    filter_if(filter_ex) %>%
    select(!!!by_vars, !!dose_end, trace_vars_str)

  # calculate the last dose date
  res <- dataset %>%
    mutate(DOMAIN = NULL) %>%
    inner_join(dataset_ex, by = by_vars_str) %>%
    mutate_at(vars(!!dose_end, !!analysis_date),
              ~ `if`(is_date(.), convert_dtm_to_dtc(.), .)) %>%
    mutate(
      tmp_dose_end_date = convert_dtc_to_dtm(
        dtc = !!dose_end,
        date_imputation = NULL,
        time_imputation = "00:00:00"
      ),
      tmp_analysis_date = convert_dtc_to_dtm(
        dtc = !!analysis_date,
        date_imputation = NULL,
        time_imputation = "23:59:59"
      )
    ) %>%
    group_by(!!!by_vars, !!dataset_seq_var)

  # if no traceability variables are required, simply calculate the last dose date
  if (is.null(traceability_vars)) {
    res <- res %>%
      summarise(
        ldose_idx = compute_ldose_idx(dose_end = .data$tmp_dose_end_date,
                                      analysis_date = .data$tmp_analysis_date),
        !!new_var := as.POSIXct(as.character(.data$tmp_dose_end_date[.data$ldose_idx]),
                                tz = lubridate::tz(.data$tmp_dose_end_date))
      ) %>%
      ungroup()
  } else {
    # calculate the last dose date and get the appropriate traceability variables
    res <- res %>%
      mutate(
        ldose_idx = compute_ldose_idx(dose_end = .data$tmp_dose_end_date,
                                      analysis_date = .data$tmp_analysis_date),
        !!new_var := as.POSIXct(as.character(.data$tmp_dose_end_date[.data$ldose_idx]),
                                tz = lubridate::tz(.data$tmp_dose_end_date))
      ) %>%
      mutate_at(trace_vars_str, list(~ .[.data$ldose_idx])) %>%
      distinct(!!new_var, !!!syms(trace_vars_str)) %>%
      ungroup()
  }

  # return either date or date-time variable
  if (!output_datetime) {
    res <- mutate(res, !!new_var := as.Date(!!new_var))
  }

  # return dataset with additional column
  left_join(dataset,
            distinct(res,
                     !!!by_vars, !!dataset_seq_var, !!new_var,
                     !!!syms(trace_vars_str)),
            by = c(by_vars_str, as_string(quo_get_expr(dataset_seq_var))))
}
