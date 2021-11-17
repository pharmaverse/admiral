#' Derive Last Dose
#'
#' Add EX source variables from last dose to the input dataset.
#' @param dataset Input dataset.
#' @param dataset_ex Input EX dataset.
#' @param filter_ex Filtering condition applied to EX dataset.
#' For example, it can be used to filter for valid dose.
#' Defaults to NULL.
#' @param by_vars Variables to join by (created by `dplyr::vars`).
#' @param dose_id Variables to identify unique dose (created by `dplyr::vars`).
#' Defaults to empty `vars()`.
#' @param ex_keep_vars Variables to keep from `dataset_ex`, with the option to rename. Can either
#' be variables created by `dplyr::vars` (e.g. `vars(VISIT)`), or named list returned by [`vars()`]
#' (e.g. `vars(LSTEXVIS = VISIT)`). If set to `NULL`, then all variables from `dataset_ex` are
#' kept without renaming.
#' Defaults to `NULL`.
#' @param dose_date The EX dose date variable.
#' @param analysis_date The analysis date variable.
#' @param dataset_seq_var The sequence variable from `dataset`.
#' (this together with `by_vars` creates the keys of `dataset`).
#' @param single_dose_condition The condition for checking if `dataset_ex` is single dose. An error
#' is issued if the condition is not true. Defaults to `(EXDOSFRQ == "ONCE")`.
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
#' * `dose_date`: no date imputation, time imputation to `00:00:00` if time is missing.
#' * `analysis_date`: no date imputation, time imputation to `23:59:59` if time is missing.
#'
#' The last dose records are identified as follows:
#' 1. The `dataset_ex` is filtered using `filter_ex`, if provided.
#' This is useful for, for example, filtering for valid dose only.
#' 2. The datasets `dataset` and `dataset_ex` are joined using `by_vars`.
#' 3. The last dose is identified:
#' the last dose is the EX record with maximum date where `dose_date` is lower to or equal to
#' `analysis_date`, subject to both date values are non-NA.
#' The last dose is identified per `by_vars` and `dataset_seq_var`.
#' If multiple EX records exist for the same `dose_date`, then either `dose_id`
#' needs to be supplied (e.g. `dose_id = vars(EXSEQ)`) to identify unique records,
#' or an error is issued. When `dose_id` is supplied, the last EX record from the same `dose_date`
#' sorted by `dose_id` will be used to identify last dose.
#' 4. The EX source variables (as specified in `ex_keep_vars`) from last dose are appended to the
#' `dataset` and returned to the user.
#'
#' This function only works correctly for EX dataset with a structure of single dose per row.
#' If your study EX dataset has multiple doses per row, use `expansion_function_name??` to
#' transform the EX dataset into single dose per row structure before calling `derive_last_dose`.
#'
#' If variables (other than those specified in `by_vars`) exist in both `dataset` and `dataset_ex`,
#' then join cannot be performed properly and an error is issued. To resolve the error, use
#' `ex_keep_vars` to either keep variables unique to `dataset_ex`, or use this option to rename
#' variables from `dataset_ex` (e.g. `ex_keep_vars = vars(LSTEXVIS = VISIT)`).
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
#'     ex_keep_vars = vars(EXDOSE, EXTRT, EXSEQ, EXENDTC, VISIT),
#'     dose_date = EXENDTC,
#'     analysis_date = AESTDTC,
#'     single_dose_condition = (EXSTDTC == EXENDTC),
#'     dataset_seq_var = AESEQ
#'   ) %>%
#'   select(STUDYID, USUBJID, AESEQ, AESTDTC, EXDOSE, EXTRT, EXENDTC, EXSEQ, VISIT)
#'
#' # or with traceability variables
#' ae %>%
#'   head(100) %>%
#'   derive_last_dose(
#'     head(ex_single, 100),
#'     filter_ex = (EXDOSE > 0 | (EXDOSE == 0 & grepl("PLACEBO", EXTRT))) &
#'       nchar(EXENDTC) >= 10,
#'     ex_keep_vars = vars(EXDOSE, EXTRT, EXSEQ, EXENDTC, VISIT),
#'     dose_date = EXENDTC,
#'     analysis_date = AESTDTC,
#'     single_dose_condition = (EXSTDTC == EXENDTC),
#'     dataset_seq_var = AESEQ,
#'     traceability_vars = dplyr::vars(LDOSEDOM = "EX", LDOSESEQ = EXSEQ, LDOSEVAR = "EXENDTC")
#'   ) %>%
#'   select(STUDYID, USUBJID, AESEQ, AESTDTC, EXDOSE, EXTRT, EXENDTC, LDOSEDOM, LDOSESEQ, LDOSEVAR)
#'
derive_last_dose <- function(dataset,
                             dataset_ex,
                             filter_ex = NULL,
                             by_vars = vars(STUDYID, USUBJID),
                             dose_id = vars(),
                             dose_date,
                             analysis_date,
                             dataset_seq_var,
                             single_dose_condition = (EXDOSFRQ == "ONCE"),
                             ex_keep_vars = NULL,
                             traceability_vars = NULL) {

  filter_ex <- assert_filter_cond(enquo(filter_ex), optional = TRUE)
  by_vars <- assert_vars(by_vars)
  dose_id <- assert_vars(dose_id)
  dose_date <- assert_symbol(enquo(dose_date))
  analysis_date <- assert_symbol(enquo(analysis_date))
  dataset_seq_var <- assert_symbol(enquo(dataset_seq_var))
  single_dose_condition <- assert_filter_cond(enquo(single_dose_condition))
  assert_varval_list(ex_keep_vars, optional = TRUE, accept_var = TRUE)
  assert_varval_list(traceability_vars, optional = TRUE)
  assert_data_frame(dataset, quo_c(by_vars, analysis_date, dataset_seq_var))
  assert_data_frame(dataset_ex, quo_c(by_vars, dose_date))

  # vars converted to string
  by_vars_str <- vars2chr(by_vars)

  # check if single_dose_condition is true for all EX records
  single_dose_eval <- dataset_ex %>%
    with(eval_tidy(quo_get_expr(single_dose_condition))) %>%
    all()

  if (!single_dose_eval){
    stop("Specified single_dose_condition is not satisfied.")
  }

  # check if doses are unique based on dose_date and dose_id
  signal_duplicate_records(dataset_ex, c(by_vars, dose_date, dose_id),
  "Multiple doses exists for the same dose_date. Update dose_id to identify unique doses.")

  # filter EX based on user-specified condition
  if (!is.null(quo_get_expr(filter_ex))) {
    dataset_ex <- dataset_ex %>%
      filter_if(filter_ex)
  }

  # create traceability vars if requested
  if (!is.null(traceability_vars)) {
    trace_vars_str <- names(traceability_vars)
    dataset_ex <- mutate(dataset_ex, !!!traceability_vars)
  } else {
    trace_vars_str <- character(0)
  }

  # keep user-specified variables from EX, if no variables specified all EX variables are kept
  if (is.null(ex_keep_vars)) {
    ex_keep_vars_str <- colnames(dataset_ex)
  } else {
    ex_keep_vars_str <- ifelse(names(ex_keep_vars) == "",
                               vars2chr(ex_keep_vars),
                               names(ex_keep_vars))
  }

  dataset_ex <- dataset_ex %>% mutate(!!!ex_keep_vars) %>%
    select(!!!by_vars, !!!dose_id, !!!dose_date, !!!syms(ex_keep_vars_str), !!!syms(trace_vars_str))

  # check if any variable exist in both dataset and dataset_ex (except for by_vars) before join
  dataset_var <- colnames(dataset)[!colnames(dataset) %in% by_vars_str]
  dataset_ex_var <- colnames(dataset_ex)[!colnames(dataset_ex) %in% by_vars_str]
  dup_var <- enumerate(intersect(dataset_var, dataset_ex_var))

  if (length(intersect(dataset_var, dataset_ex_var)) != 0){
    stop("Variable(s) ", paste(dup_var, "found in both datasets, cannot perform join"))
  }

  # create tmp numeric date to enable date comparison for identifying last dose
  dataset_ex <- dataset_ex %>%
    mutate(
      tmp_dose_date = convert_date_to_dtm(
        dt = !!dose_date,
        date_imputation = NULL,
        time_imputation = "00:00:00"
      )
    )

   dataset_tmp <- dataset %>%
    mutate(
      tmp_analysis_date = convert_date_to_dtm(
        dt = !!analysis_date,
        date_imputation = NULL,
        time_imputation = "23:59:59"
      )
    )

  # join datasets and keep unique last dose records (where dose_date is before or on analysis_date)
  # for each by_var and dataset_seq_var
  res <- dataset_tmp %>%
    inner_join(dataset_ex, by = by_vars_str) %>%
    group_by(!!!by_vars, !!dataset_seq_var) %>%
    mutate(
      tmp_ldose_idx = compute_ldose_idx(dose_date = .data$tmp_dose_date,
                                        analysis_date = .data$tmp_analysis_date),
      tmp_ldose_dt = as.POSIXct(as.character(.data$tmp_dose_date[.data$tmp_ldose_idx]),
                                tz = lubridate::tz(.data$tmp_dose_date))
    ) %>%
    ungroup() %>%
    filter(tmp_dose_date == tmp_ldose_dt) %>%
    filter_extreme(by_vars = c(by_vars, vars(!!dataset_seq_var, tmp_dose_date)),
                   order = c(by_vars, vars(!!dataset_seq_var, tmp_dose_date), dose_id),
                   mode = "last") %>%
    select(!!!by_vars, !!dataset_seq_var, !!!syms(ex_keep_vars_str), !!!syms(trace_vars_str))

  # return original dataset and add EX source variables from last dose record
  left_join(dataset, res,
            by = c(by_vars_str, as_string(quo_get_expr(dataset_seq_var))))

}


#' Helper function to get the index of last dose date
#'
#' @param dose_date dose end date
#' @param analysis_date analysis date
#'
#' @noRd
#'
#' @return index. The last dose date is then `dose_date[return_value]`
compute_ldose_idx <- function(dose_date, analysis_date) {
  if (any(!is.na(dose_date) & !is.na(analysis_date)) && any(dose_date <= analysis_date)) {
    which.max(dose_date[dose_date <= analysis_date])
  } else {
    NA_integer_
  }
}
