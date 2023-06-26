#' Derive Last Dose Amount
#'
#' @description Add a variable for dose amount from the last dose to the input dataset.
#'
#' **Note:** This is a wrapper function for the function `derive_vars_last_dose()`.
#'
#' @description
#' `r lifecycle::badge("deprecated")`
#'
#' This function is *deprecated*, please use `derive_vars_joined()` instead.
#'
#' @inheritParams derive_vars_last_dose
#' @param new_var The new variable added to `dataset`.
#' @param dose_var The EX source dose amount variable. Defaults to `EXDOSE`.
#'
#' @details The last dose amount is derived as the dose amount where the maximum `dose_date` is
#' lower to or equal to the `analysis_date` per `by_vars` for each observation in `dataset`.
#'
#' If dose information is aggregated (i.e. is a dosing frequency other than `"ONCE"`
#' over a period defined by a start and end date) the function
#' `create_single_dose_dataset()` can be used to generate single doses from
#' aggregate dose information and satisfy `single_dose_condition`.
#'
#' @return Input dataset with additional column `new_var`.
#'
#'
#' @family deprecated
#' @keywords deprecated
#'
#' @export
#'
derive_var_last_dose_amt <- function(dataset,
                                     dataset_ex,
                                     filter_ex = NULL,
                                     by_vars = exprs(STUDYID, USUBJID),
                                     dose_id = exprs(),
                                     dose_date,
                                     analysis_date,
                                     single_dose_condition = (EXDOSFRQ == "ONCE"),
                                     new_var,
                                     dose_var = EXDOSE,
                                     traceability_vars = NULL) {
  deprecate_stop("0.11.0", "derive_var_last_dose_amt()", "derive_vars_joined()")
}
