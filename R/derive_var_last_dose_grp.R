#' Derive Last Dose with User-Defined Groupings
#'
#' Add a variable for user-defined dose grouping of the last dose to the input dataset.
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
#' @author Ben Straub
#'
#' @family der_gen
#' @keywords der_gen
#'
#' @export
#'
#' @seealso [derive_vars_last_dose()], [cut()], [create_single_dose_dataset()]
#'
#' @examples
#' library(dplyr, warn.conflicts = FALSE)
#' library(admiral.test)
#' data(admiral_ae)
#' data(ex_single)
#'
#' admiral_ae %>%
#'   head(100) %>%
#'   derive_var_last_dose_grp(
#'     head(ex_single, 100),
#'     filter_ex = (EXDOSE > 0 | (EXDOSE == 0 & grepl("PLACEBO", EXTRT))) &
#'       nchar(EXENDTC) >= 10,
#'     by_vars = vars(STUDYID, USUBJID),
#'     dose_date = EXSTDTC,
#'     new_var = LDGRP,
#'     grp_brks = c(0, 20, 40, 60),
#'     grp_lbls = c("Low", "Medium", "High"),
#'     include_lowest = TRUE,
#'     right = TRUE,
#'     dose_var = EXDOSE,
#'     analysis_date = AESTDTC,
#'     single_dose_condition = (EXSTDTC == EXENDTC),
#'     traceability_vars = vars(LDOSEDOM = "EX", LDOSESEQ = EXSEQ, LDOSEVAR = "EXENDTC")
#'   ) %>%
#'   select(USUBJID, LDGRP, LDOSEDOM, LDOSESEQ, LDOSEVAR)
derive_var_last_dose_grp <- function(dataset,
                                     dataset_ex,
                                     filter_ex = NULL,
                                     by_vars = vars(STUDYID, USUBJID),
                                     dose_id = vars(),
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
  filter_ex <- assert_filter_cond(enquo(filter_ex), optional = TRUE)
  by_vars <- assert_vars(by_vars)
  dose_date <- assert_symbol(enquo(dose_date))
  analysis_date <- assert_symbol(enquo(analysis_date))
  single_dose_condition <- assert_filter_cond(enquo(single_dose_condition))
  new_var <- assert_symbol(enquo(new_var))
  dose_var <- assert_symbol(enquo(dose_var))

  derive_vars_last_dose(
    dataset = dataset,
    dataset_ex = dataset_ex,
    filter_ex = !!filter_ex,
    by_vars = by_vars,
    dose_id = dose_id,
    dose_date = !!dose_date,
    analysis_date = !!analysis_date,
    single_dose_condition = !!single_dose_condition,
    new_vars = vars(!!dose_var),
    traceability_vars = traceability_vars
  ) %>%
    mutate(
      !!new_var :=
        as.character(
          cut(
            !!dose_var,
            breaks = !!grp_brks,
            include.lowest = include_lowest,
            right = right,
            labels = !!grp_lbls
          )
        )
    ) %>%
    select(-!!dose_var, !!new_var)
}
