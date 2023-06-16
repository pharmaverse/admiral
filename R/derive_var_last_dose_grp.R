#' Derive Last Dose with User-Defined Groupings
#'
#' @description Add a variable for user-defined dose grouping of the last dose
#' to the input dataset.
#'
#' **Note:** This is a wrapper function for the function `derive_vars_last_dose()`.
#'
#' @description
#' `r lifecycle::badge("deprecated")`
#'
#' This function is *deprecated*, please use `derive_vars_joined()` instead.
#'
#' @inheritParams derive_vars_last_dose
#' @param new_var The output variable defined by the user.
#' @param grp_brks User supplied breaks to apply to groups.
#' Refer to `breaks` parameter in `cut()` for details.
#' @param grp_lbls User supplied labels to apply to groups.
#' Refer to `labels` parameter in `cut()` for details.
#' @param dose_var The source dose amount variable. Defaults to `EXDOSE`.
#' @param include_lowest logical, indicating if a value equal to the lowest
#' (or highest, for right = FALSE) ‘breaks’ value should be included.
#' Refer to `include.lowest` parameter in `cut()` for details.
#' @param right Logical, indicating if the intervals should be closed on the right
#' (and open on the left) or vice versa.
#' Refer to `right` parameter in `cut()` for details.
#'
#' @details Last dose is the dose with maximum `dose_date` that is lower to or equal to the
#' `analysis_date` per `by_vars` for each observation in `dataset`.
#' The last dose group is then derived by user-defined grouping, which groups
#' `dose_var` as specified in `grp_brks`, and returns `grp_lbls` as the values for `new_var`.
#'
#' If dose information is aggregated (i.e. is a dosing frequency other than `"ONCE"`
#' over a period defined by a start and end date) the function
#' `create_single_dose_dataset()` can be used to generate single doses from
#' aggregate dose information and satisfy `single_dose_condition`.
#' @return Input dataset with additional column `new_var`.
#'
#'
#' @family deprecated
#' @keywords deprecated
#'
#' @export
#'
derive_var_last_dose_grp <- function(dataset,
                                     dataset_ex,
                                     filter_ex = NULL,
                                     by_vars = exprs(STUDYID, USUBJID),
                                     dose_id = exprs(),
                                     dose_date,
                                     analysis_date,
                                     single_dose_condition = (EXDOSFRQ == "ONCE"),
                                     new_var,
                                     grp_brks,
                                     grp_lbls,
                                     include_lowest = TRUE,
                                     right = TRUE,
                                     dose_var = EXDOSE,
                                     traceability_vars = NULL) {
  deprecate_stop("0.11.0", "derive_var_last_dose_grp()", "derive_vars_joined()")
}
