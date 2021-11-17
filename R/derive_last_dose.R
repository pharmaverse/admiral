#' Derive Last Dose
#'
#' @param dataset Input dataset.
#' @param dataset_ex Input EX dataset.
#' @param filter_ex Filtering condition applied to EX dataset.
#' For example, it can be used to filter for valid dose.
#' Defaults to NULL.
#' @param by_vars Variables to join by (created by `dplyr::vars`).
#' @param dose_id Variables to identify unique dose (created by `dplyr::vars`).
#' Defaults to empty `vars()`.
#' @param ex_keep_vars Variables to keep from `ex_dataset`, with the option to rename. Can either
#' be variables created by `dplyr::vars` (e.g. `vars(VISIT)`), or named list returned by [`vars()`]
#' (e.g. `vars(LSTEXVIS = VISIT)`). If set to `NULL`, then all variables from `ex_dataset` are
#' kept without renaming. Defaults to `NULL`.
#' @param dose_start The dose start date variable.
#' @param dose_end The dose end date variable.
#' @param analysis_date The analysis date variable.
#' @param dataset_seq_var The sequence variable
#' (this together with `by_vars` creates the keys of `dataset`).
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
#' When doing date comparison to identify last dose, date-time imputations are done as follows:
#' * `dose_end`: no date imputation, time imputation to `00:00:00` if time is missing.
#' * `analysis_date`: no date imputation, time imputation to `23:59:59` if time is missing.
#'
#' The last dose records are identified as follows:
#' 1. The `dataset_ex` is filtered using `filter_ex`, if provided.
#' This is useful for, for example, filtering for valid dose only.
#' 2. The datasets `dataset` and `dataset_ex` are joined using `by_vars`.
#' 3. The last dose is identified:
#' the last dose is the EX record with maximum date where `dose_end` is lower to or equal to
#' `analysis_date`, subject to both date values are non-NA.
#' The last dose is identified per `by_vars` and `dataset_seq_var`.
#' If multiple such EX records exist for the same `dose_end` date, then either `dose_id`
#' needs to be supplied (e.g. `dose_id = vars(EXSEQ)`) to identify unique records,
#' or an error is issued.
#' 4. The EX source variables from last dose are appended to the `dataset` and returned to the user.
#'
#' Furthermore, the following assumption is checked: start and end dates (datetimes) need to match.
#' Use `check_dates_only` to control whether only dates or whole date-times need to be equal. This
#' is required because if start and end dates (datetimes) don't match, the `analysis_date` can occur
#' between `dose_start` and `dose_end`. When this happens, the function will choose the dose with a
#' latest `dose_end` date prior to `analysis_date`, as opposed to the actual last dose.
#'
#' If variables (other than those specified in `by_vars`) exist in both `dataset` and `ex_dataset`,
#' then join cannot be performed properly and an error is issued. To resolve the error, use
#' `ex_keep_vars` to either keep variables unique to `ex_dataset`, or use this option to rename
#' variables from `ex_dataset` (e.g. `ex_keep_vars = vars(LSTEXVIS = VISIT)`).
#'
#' @return Input dataset with EX source variables from last dose added.
#'
#' @author Ondrej Slama, Annie Yang
#'
#' @keywords adam derivation user_utility
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
#'     ex_keep_vars = vars(EXSTDTC, EXENDTC, EXDOSE, EXTRT, EXSEQ, VISIT),
#'     dose_start = EXSTDTC,
#'     dose_end = EXENDTC,
#'     analysis_date = AESTDTC,
#'     dataset_seq_var = AESEQ,
#'     check_dates_only = FALSE
#'   ) %>%
#'   select(STUDYID, USUBJID, AESEQ, AESTDTC, EXDOSE, EXTRT, EXENDTC)
#'
#' # or with traceability variables
#' ae %>%
#'   head(100) %>%
#'   derive_last_dose(
#'     head(ex_single, 100),
#'     filter_ex = (EXDOSE > 0 | (EXDOSE == 0 & grepl("PLACEBO", EXTRT))) &
#'       nchar(EXENDTC) >= 10,
#'     ex_keep_vars = vars(EXSTDTC, EXENDTC, EXDOSE, EXTRT, EXSEQ, VISIT),
#'     dose_start = EXSTDTC,
#'     dose_end = EXENDTC,
#'     analysis_date = AESTDTC,
#'     dataset_seq_var = AESEQ,
#'     check_dates_only = FALSE,
#'     traceability_vars = dplyr::vars(LDOSEDOM = "EX", LDOSESEQ = EXSEQ, LDOSEVAR = "EXENDTC")
#'   ) %>%
#'   select(STUDYID, USUBJID, AESEQ, AESTDTC, EXDOSE, EXTRT, EXENDTC, LDOSEDOM, LDOSESEQ, LDOSEVAR)
#'
derive_last_dose <- function(dataset,
                             dataset_ex,
                             filter_ex = NULL,
                             by_vars = vars(STUDYID, USUBJID),
                             dose_id = vars(),
                             ex_keep_vars = NULL,
                             dose_start,
                             dose_end,
                             analysis_date,
                             dataset_seq_var,
                             check_dates_only = FALSE,
                             traceability_vars = NULL) {

  filter_ex <- assert_filter_cond(enquo(filter_ex), optional = TRUE)
  by_vars <- assert_vars(by_vars)
  dose_id <- assert_vars(dose_id)
  dose_start <- assert_symbol(enquo(dose_start))
  dose_end <- assert_symbol(enquo(dose_end))
  analysis_date <- assert_symbol(enquo(analysis_date))
  dataset_seq_var <- assert_symbol(enquo(dataset_seq_var))
  assert_logical_scalar(check_dates_only)
  stopifnot(is_quosures(traceability_vars) | is.null(traceability_vars))
  stopifnot(is_quosures(ex_keep_vars) | is.null(ex_keep_vars))
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

  # filter based on user-specified condition
  if (!is.null(quo_get_expr(filter_ex))) {
    dataset_ex <- dataset_ex %>%
      filter_if(filter_ex)
  }

  # keep user-specified variables from EX. If no variables specified, all EX variables are kept.
  if (is.null(ex_keep_vars)) {
    ex_keep_vars_str <- colnames(dataset_ex)
  } else {
    ex_keep_vars_str <- ifelse(names(ex_keep_vars) == "",
                               vars2chr(ex_keep_vars),
                               names(ex_keep_vars))
  }

  dataset_ex <- mutate(dataset_ex, !!!ex_keep_vars) %>%
    select(!!!by_vars, !!!syms(ex_keep_vars_str), !!!syms(trace_vars_str))

  # issue an error when variables exist in both dataset and dataset_ex, except for by_vars.
  dataset_var <- colnames(dataset)[!colnames(dataset) %in% by_vars_str]
  dataset_ex_var <- colnames(dataset_ex)[!colnames(dataset_ex) %in% by_vars_str]
  dup_var <- paste(intersect(dataset_var, dataset_ex_var), collapse = " ")

  if (dup_var != ""){
    stop("Variable(s) ", paste(dup_var, "found in both datasets, cannot perform join"))
  }


  # join dataset with ex and create tmp numeric date to enable comparison
  res <- dataset %>%
    inner_join(dataset_ex, by = by_vars_str) %>%
    mutate_at(vars(!!dose_end, !!analysis_date),
              list(tmp = ~ `if`(is_date(.), convert_dtm_to_dtc(.), .))) %>%
    rename(tmp_dose_end_date = paste(quo_get_expr(dose_end), "tmp", sep = "_"),
           tmp_analysis_date = paste(quo_get_expr(analysis_date), "tmp", sep = "_")) %>%
    mutate(
      tmp_dose_end_date = convert_dtc_to_dtm(
        dtc = tmp_dose_end_date,
        date_imputation = NULL,
        time_imputation = "00:00:00"
      ),
      tmp_analysis_date = convert_dtc_to_dtm(
        dtc = tmp_analysis_date,
        date_imputation = NULL,
        time_imputation = "23:59:59"
      )
    ) %>%
    group_by(!!!by_vars, !!dataset_seq_var)

  # filter last dose records with dose_end before or on analysis_date
  res <- res %>%
    mutate(
      tmp_ldose_idx = compute_ldose_idx(dose_end = .data$tmp_dose_end_date,
                                        analysis_date = .data$tmp_analysis_date),
      tmp_ldose_dt = as.POSIXct(as.character(.data$tmp_dose_end_date[.data$tmp_ldose_idx]),
                                tz = lubridate::tz(.data$tmp_dose_end_date))
    ) %>%
    filter(tmp_dose_end_date == tmp_ldose_dt)

  # issue an error if multiple last dose on the same date exist
  signal_duplicate_records(res, c(by_vars, dataset_seq_var, dose_end, dose_id),
                           "Last dose is not unique.")

  # get unique last dose records
  res <- res %>%
    arrange(!!!by_vars, !!dataset_seq_var, tmp_dose_end_date, !!!dose_id) %>%
    slice(n()) %>%
    select(!!!syms(trace_vars_str), !!!by_vars, !!dataset_seq_var, !!!syms(ex_keep_vars_str)) %>%
    ungroup()


  # return dataset with added EX source variables from last dose record
  left_join(dataset,
            distinct(res,
                     !!!by_vars, !!dataset_seq_var, .keep_all = TRUE),
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
