#' Derive Last Dose
#'
#' Add EX source variables from last dose to the input dataset.
#'
#' @description
#' `r lifecycle::badge("deprecated")`
#'
#' This function is *deprecated*, please use `derive_vars_joined()` instead.
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
#' @param by_vars Variables to join by (created by `rlang::exprs`).
#'
#' @param dose_id Variables to identify unique dose (created by `rlang::exprs`).
#' Defaults to empty `exprs()`.
#'
#' @param new_vars Variables to keep from `dataset_ex`, with the option to
#'   rename. Can either be variables created by `rlang::exprs` (e.g.
#'   `exprs(VISIT)`), or named list returned by [`exprs()`] (e.g.
#'   `exprs(LSTEXVIS = VISIT)`). If set to `NULL`, then all variables from
#'   `dataset_ex` are kept without renaming. Defaults to `NULL`.
#'
#' @param dose_date The EX dose date variable. A date or date-time object is expected.
#'
#' @param analysis_date The analysis date variable. A date or date-time object is expected.
#'
#' @param single_dose_condition The condition for checking if `dataset_ex` is single dose. An error
#' is issued if the condition is not true. Defaults to `(EXDOSFRQ == "ONCE")`.
#'
#' @param traceability_vars A named list returned by [`exprs()`] listing the traceability variables,
#' e.g. `exprs(LDOSEDOM = "EX", LDOSESEQ = EXSEQ)`.
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
#' needs to be supplied (e.g. `dose_id = exprs(EXSEQ)`) to identify unique records,
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
#' variables from `dataset_ex` (e.g. `new_vars = exprs(LSTEXVIS = VISIT)`).
#'
#' @return Input dataset with EX source variables from last dose added.
#'
#'
#' @family deprecated
#' @keywords deprecated
#'
#' @export
derive_vars_last_dose <- function(dataset,
                                  dataset_ex,
                                  filter_ex = NULL,
                                  by_vars = exprs(STUDYID, USUBJID),
                                  dose_id = exprs(),
                                  dose_date,
                                  analysis_date,
                                  single_dose_condition = EXDOSFRQ == "ONCE",
                                  new_vars = NULL,
                                  traceability_vars = NULL) {
  deprecate_stop("0.11.0", "derive_vars_last_dose()", "derive_vars_joined()")
}
