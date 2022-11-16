#' Derive Last Dose
#'
#' Add EX source variables from last dose to the input dataset.
#'
#' @param dataset Input dataset.
#' The variables specified by the `by_vars` and `analysis_date` parameters are expected.
#'
#' @param dataset_ex Input EX dataset.
#' The variables specified by the `by_vars`, `dose_date`, `new_vars` parameters,
#' and source variables from `traceability_vars` parameter are expected.
#'
#' @param filter_ex Filtering condition applied to EX dataset.
#' For example, it can be used to filter for valid dose.
#' Defaults to NULL.
#'
#' @param by_vars Variables to join by (created by `dplyr::vars`).
#'
#' @param dose_id Variables to identify unique dose (created by `dplyr::vars`).
#' Defaults to empty `vars()`.
#'
#' @param new_vars Variables to keep from `dataset_ex`, with the option to rename. Can either
#' be variables created by `dplyr::vars` (e.g. `vars(VISIT)`), or named list returned by [`vars()`]
#' (e.g. `vars(LSTEXVIS = VISIT)`). If set to `NULL`, then all variables from `dataset_ex` are
#' kept without renaming.
#' Defaults to `NULL`.
#'
#' @param dose_date The EX dose date variable. A date or date-time object is expected.
#'
#' @param analysis_date The analysis date variable. A date or date-time object is expected.
#'
#' @param single_dose_condition The condition for checking if `dataset_ex` is single dose. An error
#' is issued if the condition is not true. Defaults to `(EXDOSFRQ == "ONCE")`.
#'
#' @param traceability_vars A named list returned by [`vars()`] listing the traceability variables,
#' e.g. `vars(LDOSEDOM = "EX", LDOSESEQ = EXSEQ)`.
#' The left-hand side (names of the list elements) gives the names of the traceability variables
#' in the returned dataset.
#' The right-hand side (values of the list elements) gives the values of the traceability variables
#' in the returned dataset.
#' These can be either strings or symbols referring to existing variables.
#'
#' @details
#' When doing date comparison to identify last dose, date-time imputations are done as follows:
#' * `dose_date`: time is imputed to `00:00:00` if the variable is a date variable
#' * `analysis_date`: time is imputed to `23:59:59` if the variable is a date variable
#'
#' The last dose records are identified as follows:
#'
#' 1. The `dataset_ex` is filtered using `filter_ex`, if provided.
#' This is useful for, for example, filtering for valid dose only.
#' 2. The datasets `dataset` and `dataset_ex` are joined using `by_vars`.
#' 3. The last dose is identified:
#' the last dose is the EX record with maximum date where `dose_date` is lower to or equal to
#' `analysis_date`, subject to both date values are non-NA.
#' The last dose is identified per `by_vars`.
#' If multiple EX records exist for the same `dose_date`, then either `dose_id`
#' needs to be supplied (e.g. `dose_id = vars(EXSEQ)`) to identify unique records,
#' or an error is issued. When `dose_id` is supplied, the last EX record from the same `dose_date`
#' sorted by `dose_id` will be used to identify last dose.
#' 4. The EX source variables (as specified in `new_vars`) from last dose are appended to the
#' `dataset` and returned to the user.
#'
#' This function only works correctly for EX dataset with a structure of single dose per row.
#' If your study EX dataset has multiple doses per row, use [`create_single_dose_dataset()`] to
#' transform the EX dataset into single dose per row structure before calling
#' `derive_vars_last_dose()`.
#'
#' If variables (other than those specified in `by_vars`) exist in both `dataset` and `dataset_ex`,
#' then join cannot be performed properly and an error is issued. To resolve the error, use
#' `new_vars` to either keep variables unique to `dataset_ex`, or use this option to rename
#' variables from `dataset_ex` (e.g. `new_vars = vars(LSTEXVIS = VISIT)`).
#'
#' @return Input dataset with EX source variables from last dose added.
#'
#' @author Ondrej Slama, Annie Yang
#'
#' @family der_gen
#' @keywords der_gen
#'
#' @seealso [derive_var_last_dose_amt()], [derive_var_last_dose_date()],
#'   [derive_var_last_dose_grp()], [create_single_dose_dataset()]
#'
#' @export
#'
#' @examples
#' library(dplyr, warn.conflicts = FALSE)
#' library(admiral.test)
#' data(admiral_ae)
#' data(ex_single)
#'
#' # create datetime variables in input datasets
#' ex_single <- derive_vars_dtm(
#'   head(ex_single, 100),
#'   dtc = EXENDTC,
#'   new_vars_prefix = "EXEN",
#'   flag_imputation = "none"
#' )
#'
#' adae <- admiral_ae %>%
#'   head(100) %>%
#'   derive_vars_dtm(
#'     dtc = AESTDTC,
#'     new_vars_prefix = "AST",
#'     highest_imputation = "M"
#'   )
#'
#' # add last dose vars
#' adae %>%
#'   derive_vars_last_dose(
#'     dataset_ex = ex_single,
#'     filter_ex = (EXDOSE > 0 | (EXDOSE == 0 & grepl("PLACEBO", EXTRT))) &
#'       !is.na(EXENDTM),
#'     new_vars = vars(EXDOSE, EXTRT, EXSEQ, EXENDTC, VISIT),
#'     dose_date = EXENDTM,
#'     analysis_date = ASTDTM
#'   ) %>%
#'   select(STUDYID, USUBJID, AESEQ, AESTDTC, EXDOSE, EXTRT, EXENDTC, EXSEQ, VISIT)
#'
#' # or with traceability variables
#' adae %>%
#'   derive_vars_last_dose(
#'     dataset_ex = ex_single,
#'     filter_ex = (EXDOSE > 0 | (EXDOSE == 0 & grepl("PLACEBO", EXTRT))) &
#'       !is.na(EXENDTM),
#'     new_vars = vars(EXDOSE, EXTRT, EXSEQ, EXENDTC, VISIT),
#'     dose_date = EXENDTM,
#'     analysis_date = ASTDTM,
#'     traceability_vars = vars(LDOSEDOM = "EX", LDOSESEQ = EXSEQ, LDOSEVAR = "EXENDTC")
#'   ) %>%
#'   select(STUDYID, USUBJID, AESEQ, AESTDTC, EXDOSE, EXTRT, EXENDTC, LDOSEDOM, LDOSESEQ, LDOSEVAR)
derive_vars_last_dose <- function(dataset,
                                  dataset_ex,
                                  filter_ex = NULL,
                                  by_vars = vars(STUDYID, USUBJID),
                                  dose_id = vars(),
                                  dose_date,
                                  analysis_date,
                                  single_dose_condition = EXDOSFRQ == "ONCE",
                                  new_vars = NULL,
                                  traceability_vars = NULL) {
  filter_ex <- assert_filter_cond(enquo(filter_ex), optional = TRUE)
  by_vars <- assert_vars(by_vars)
  dose_id <- assert_vars(dose_id)
  dose_date <- assert_symbol(enquo(dose_date))
  analysis_date <- assert_symbol(enquo(analysis_date))
  single_dose_condition <- assert_filter_cond(enquo(single_dose_condition))
  assert_varval_list(new_vars, optional = TRUE, accept_var = TRUE)
  assert_varval_list(traceability_vars, optional = TRUE)
  assert_data_frame(dataset, quo_c(by_vars, analysis_date))
  if (as_name(dose_date) %in% names(new_vars)) {
    required_vars <- quo_c(by_vars, new_vars, get_source_vars(traceability_vars))
    dose_date_res <- new_vars[[as_name(dose_date)]]
  } else {
    required_vars <- quo_c(by_vars, dose_date, new_vars, get_source_vars(traceability_vars))
    dose_date_res <- dose_date
  }
  assert_data_frame(dataset_ex, required_vars)
  assert_date_var(
    dataset = dataset,
    var = !!analysis_date
  )
  assert_date_var(
    dataset = dataset_ex,
    var = !!dose_date_res
  )

  # vars converted to string
  by_vars_str <- vars2chr(by_vars)
  dose_id_str <- vars2chr(dose_id)

  # check if single_dose_condition is true for all EX records
  single_dose_eval <- dataset_ex %>%
    with(eval_tidy(quo_get_expr(single_dose_condition))) %>%
    all()

  if (!single_dose_eval) {
    stop("Specified `single_dose_condition` is not satisfied.")
  }

  # check if doses are unique based on `dose_date` and `dose_id`
  if (as_name(dose_date) %in% names(new_vars)) {
    unique_by <- c(by_vars, new_vars[[as_name(dose_date)]], dose_id)
  } else {
    unique_by <- c(by_vars, dose_date, dose_id)
  }
  signal_duplicate_records(
    dataset_ex,
    unique_by,
    "Multiple doses exist for the same `dose_date`. Update `dose_id` to identify unique doses."
  )

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
  if (!is.null(new_vars)) {
    new_vars_name <- replace_values_by_names(new_vars)
    dataset_ex <- dataset_ex %>%
      mutate(!!!new_vars) %>%
      select(
        !!!by_vars, !!!syms(dose_id_str), !!dose_date,
        !!!new_vars_name, !!!syms(trace_vars_str)
      )
  } else {
    new_vars_name <- syms(colnames(dataset_ex))
  }

  # check if any variable exist in both dataset and dataset_ex (except for by_vars) before join
  dataset_var <- colnames(dataset)[!colnames(dataset) %in% by_vars_str]
  dataset_ex_var <- colnames(dataset_ex)[!colnames(dataset_ex) %in% by_vars_str]
  dup_var <- enumerate(intersect(dataset_var, dataset_ex_var))

  if (length(intersect(dataset_var, dataset_ex_var)) != 0) {
    stop("Variable(s) ", paste(dup_var, "found in both datasets, cannot perform join"))
  }

  # create temporary observation number and temporary numeric date to identify last dose
  dataset <- dataset %>%
    derive_var_obs_number(
      order = vars(USUBJID),
      new_var = tmp_seq_var
    ) %>%
    mutate(
      tmp_analysis_date = convert_date_to_dtm(
        dt = !!analysis_date,
        time_imputation = "last"
      )
    )

  # create tmp numeric date to enable date comparison for identifying last dose
  dataset_ex <- dataset_ex %>%
    mutate(
      tmp_dose_date = convert_date_to_dtm(
        dt = !!dose_date
      )
    )


  # join datasets and keep unique last dose records (where dose_date is before or on analysis_date)
  res <- dataset %>%
    select(!!!by_vars, tmp_seq_var, tmp_analysis_date) %>%
    inner_join(dataset_ex, by = by_vars_str) %>%
    filter(!is.na(tmp_dose_date) & !is.na(tmp_analysis_date) &
      tmp_dose_date <= tmp_analysis_date) %>%
    filter_extreme(
      by_vars = vars(tmp_seq_var),
      order = c(vars(tmp_dose_date), dose_id),
      mode = "last"
    ) %>%
    select(tmp_seq_var, !!!new_vars_name, !!!syms(trace_vars_str), -by_vars_str)

  # return observations from original dataset with last dose variables added
  derive_vars_merged(
    dataset,
    dataset_add = res,
    by_vars = vars(tmp_seq_var)
  ) %>% select(-starts_with("tmp_"))
}
