
#' Derive last dose date(-time)
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
#' @param traceability_vars A named list returned by `vars` listing the traceability variables,
#' e.g. `vars(LDOSEDOM = "EX", LDOSESEQ = EXSEQ)`.
#' The left-hand side (names of the list elements) gives the names of the traceability variables
#' in the returned dataset.
#' The right-hand side (values of the list elements) gives the values of the traceability variables
#' in the returned dataset.
#' Those can be either strings or symbols referring to existing variables.
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
#' @export
#'
#' @examples
#' data(ae); data(ex_single)
#' derive_last_dose(
#'   head(ae, 100),
#'   head(ex_single, 100),
#'   filter_ex = (EXDOSE > 0 | (EXDOSE == 0 & stringr::str_detect(EXTRT, "PLACEBO"))) &
#'     nchar(as.character(EXENDTC)) >= 10,
#'   dose_start = EXSTDTC,
#'   dose_end = EXENDTC,
#'   analysis_date = AESTDTC,
#'   dataset_seq_var = AESEQ,
#'   new_var = LDOSEDTM,
#'   output_datetime = TRUE,
#'   check_dates_only = FALSE
#' )
#'
#' # or with traceability variables
#' derive_last_dose(
#'   head(ae, 100),
#'   head(ex_single, 100),
#'   filter_ex = (EXDOSE > 0 | (EXDOSE == 0 & stringr::str_detect(EXTRT, "PLACEBO"))) &
#'     nchar(as.character(EXENDTC)) >= 10,
#'   dose_start = EXSTDTC,
#'   dose_end = EXENDTC,
#'   analysis_date = AESTDTC,
#'   dataset_seq_var = AESEQ,
#'   new_var = LDOSEDTM,
#'   output_datetime = TRUE,
#'   check_dates_only = FALSE,
#'   traceability_vars = dplyr::vars(LDOSEDOM = "EX", LDOSESEQ = EXSEQ, LDOSEVAR = "EXSTDTC")
#' )
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
                             traceability_vars = vars()) {

  assert_that(
    !dplyr::is_grouped_df(dataset),
    !dplyr::is_grouped_df(dataset_ex),
    is_quosures(by_vars),
    rlang::is_scalar_logical(output_datetime),
    rlang::is_scalar_logical(check_dates_only),
    is_quosures(traceability_vars)
  )

  # apply filtering condition
  filter_ex <- enquo(filter_ex)
  if (!is.null(quo_get_expr(filter_ex))) {
    dataset_ex <- filter(dataset_ex, !!filter_ex)
  }

  # substitute dataset variables
  dose_start <- enquo(dose_start)
  dose_end <- enquo(dose_end)
  analysis_date <- enquo(analysis_date)
  dataset_seq_var <- enquo(dataset_seq_var)
  new_var <- enquo(new_var)

  # check if variables are symbols
  assert_that(
    assert_is_symbol(dose_start),
    assert_is_symbol(dose_end),
    assert_is_symbol(analysis_date),
    assert_is_symbol(dataset_seq_var)
  )

  # by_vars converted to string
  by_vars_str <- vars2chr(by_vars)

  # check variables existence - dataset
  assert_has_variables(
    dataset,
    c(by_vars_str,
      as_string(quo_get_expr(analysis_date)),
      as_string(quo_get_expr(dataset_seq_var))
    )
  )

  # check variables existence - dataset_ex
  assert_has_variables(
    dataset_ex,
    c(by_vars_str,
      as_string(quo_get_expr(dose_start)),
      as_string(quo_get_expr(dose_end)))
  )

  # check variables existence - new_var
  if (as_string(quo_get_expr(new_var)) == "") {
    stop("Argument 'new_var' must be specified.")
  }

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
  if (length(traceability_vars) != 0) {
    trace_vars_str <- names(traceability_vars)
    dataset_ex <- mutate(dataset_ex, !!!traceability_vars)
  } else {
    trace_vars_str <- character(0)
  }

  # select only a subset of columns
  dataset_ex <- select(dataset_ex, !!!by_vars, !!dose_end, trace_vars_str)

  # calculate the last dose date
  res <- dataset %>%
    mutate(DOMAIN = NULL) %>%
    inner_join(dataset_ex, by = by_vars_str) %>%
    mutate_at(vars(!!dose_end, !!analysis_date),
              ~ `if`(is_date(.), convert_dtm_to_dtc(.), .)) %>%
    group_by(!!!by_vars, !!dataset_seq_var) %>%
    mutate(
      tmp_dose_end_date = impute_dtc(dtc = !!dose_end,
                                     date_imputation = NULL,
                                     time_imputation = "00:00:00") %>%
        convert_dtc_to_dtm(),
      tmp_analysis_date = impute_dtc(dtc = !!analysis_date,
                                     date_imputation = NULL,
                                     time_imputation = "23:59:59") %>%
        convert_dtc_to_dtm())

  # if no traceability variables are required, simply calculate the last dose date
  if (length(traceability_vars) == 0) {
    res <- res %>%
      summarise(ldose_idx = compute_ldose_idx(dose_end = .data$tmp_dose_end_date,
                                              analysis_date = .data$tmp_analysis_date),
                !!new_var := .data$tmp_dose_end_date[.data$ldose_idx]) %>%
      ungroup()
  } else {
    # calculate the last dose date and get the appropriate traceability variables
    res <- res %>%
      mutate(ldose_idx = compute_ldose_idx(dose_end = .data$tmp_dose_end_date,
                                           analysis_date = .data$tmp_analysis_date),
             !!new_var := .data$tmp_dose_end_date[.data$ldose_idx]) %>%
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
#' @return index. The last dose date is then `dose_end[return_value]`
compute_ldose_idx <- function(dose_end, analysis_date) {
  if (any(!is.na(dose_end) & !is.na(analysis_date)) && any(dose_end <= analysis_date)) {
    which.max(dose_end[dose_end <= analysis_date])
  } else {
    NA_integer_
  }
}
