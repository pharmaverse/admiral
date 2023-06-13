#' Derive Last Dose Date-Time
#'
#' @description Add a variable for the dose date or datetime of the last dose to
#' the input dataset.
#'
#' **Note:** This is a wrapper function for the function `derive_vars_last_dose()`.
#'
#' @description
#' `r lifecycle::badge("deprecated")`
#'
#' This function is *deprecated*, please use `derive_vars_joined()` instead.
#'
#' @inheritParams derive_vars_last_dose
#' @param new_var The new date or datetime variable added to `dataset`.
#' @param output_datetime  Display `new_var` as datetime or as date only. Defaults to `TRUE`.
#'
#' @details The last dose date is derived as the maximum dose date where the
#'   `dose_date` is lower to or equal to the `analysis_date` per `by_vars` for
#'   each observation in `dataset`. When `output_datetime` is `TRUE` and time is
#'   missing, then the last dose date time is imputed to `00:00:00`. However, if
#'   date is missing, then no imputation is done.
#'
#'  If dose information is aggregated (i.e. is a dosing frequency other than `"ONCE"`
#'  over a period defined by a start and end date) the function
#'  `create_single_dose_dataset()` can be used to generate single doses from
#'  aggregate dose information and satisfy `single_dose_condition`.
#'
#' @return Input dataset with additional column `new_var`.
#'
#'
#' @family deprecated
#' @keywords deprecated
#'
#' @export
#'
derive_var_last_dose_date <- function(dataset,
                                      dataset_ex,
                                      filter_ex = NULL,
                                      by_vars = exprs(STUDYID, USUBJID),
                                      dose_id = exprs(),
                                      dose_date,
                                      analysis_date,
                                      single_dose_condition = (EXDOSFRQ == "ONCE"),
                                      new_var,
                                      output_datetime = TRUE,
                                      traceability_vars = NULL) {
  deprecate_stop("0.11.0", "derive_var_last_dose_date()", "derive_vars_joined()")
}
