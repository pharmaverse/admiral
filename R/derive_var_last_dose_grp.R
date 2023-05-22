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
#' ex_single <- tribble(
#'   ~STUDYID, ~DOMAIN,  ~USUBJID, ~EXSEQ, ~EXDOSE, ~EXTRT,     ~EXENDTC,     ~EXSTDTC, ~EXDOSFRQ,
#'   "PILOT1",    "EX", "18-1066",     1L,      54, "XANO", "2013-07-07", "2013-07-07",    "ONCE",
#'   "PILOT1",    "EX", "18-1066",     2L,      54, "XANO", "2013-07-08", "2013-07-08",    "ONCE",
#'   "PILOT1",    "EX", "18-1066",     3L,      54, "XANO", "2013-07-09", "2013-07-09",    "ONCE",
#'   "PILOT1",    "EX", "18-1066",     4L,      54, "XANO", "2013-07-10", "2013-07-10",    "ONCE",
#'   "PILOT1",    "EX", "18-1066",     5L,      54, "XANO", "2013-07-11", "2013-07-11",    "ONCE",
#'   "PILOT1",    "EX", "18-1066",     6L,      54, "XANO", "2013-07-12", "2013-07-12",    "ONCE",
#'   "PILOT1",    "EX", "18-1066",     7L,      54, "XANO", "2013-07-13", "2013-07-13",    "ONCE",
#'   "PILOT1",    "EX", "18-1066",     8L,      54, "XANO", "2013-07-14", "2013-07-14",    "ONCE",
#'   "PILOT1",    "EX", "18-1066",     9L,      54, "XANO", "2013-07-15", "2013-07-15",    "ONCE",
#'   "PILOT1",    "EX", "18-1066",    10L,      54, "XANO", "2013-07-16", "2013-07-16",    "ONCE",
#'   "PILOT1",    "EX", "10-1083",     1L,       0, "PLAC", "2013-07-22", "2013-07-22",    "ONCE",
#'   "PILOT1",    "EX", "10-1083",     2L,       0, "PLAC", "2013-07-23", "2013-07-23",    "ONCE",
#'   "PILOT1",    "EX", "10-1083",     3L,       0, "PLAC", "2013-07-24", "2013-07-24",    "ONCE",
#'   "PILOT1",    "EX", "10-1083",     4L,       0, "PLAC", "2013-07-25", "2013-07-25",    "ONCE",
#'   "PILOT1",    "EX", "10-1083",     5L,       0, "PLAC", "2013-07-26", "2013-07-26",    "ONCE",
#'   "PILOT1",    "EX", "10-1083",     6L,       0, "PLAC", "2013-07-27", "2013-07-27",    "ONCE",
#'   "PILOT1",    "EX", "10-1083",     7L,       0, "PLAC", "2013-07-28", "2013-07-28",    "ONCE",
#'   "PILOT1",    "EX", "10-1083",     8L,       0, "PLAC", "2013-07-29", "2013-07-29",    "ONCE",
#'   "PILOT1",    "EX", "10-1083",     9L,       0, "PLAC", "2013-07-30", "2013-07-30",    "ONCE",
#'   "PILOT1",    "EX", "10-1083",    10L,       0, "PLAC", "2013-07-31", "2013-07-31",    "ONCE",
#'   "PILOT1",    "EX", "10-1083",    11L,       0, "PLAC", "2013-08-01", "2013-08-01",    "ONCE"
#' )
#'
#' ae <- tribble(
#'   ~STUDYID,  ~DOMAIN,  ~USUBJID, ~AESEQ,     ~AESTDTC,                 ~AETERM,
#'   "PILOT01",    "AE", "10-1083",      1, "2013-08-02", "MYOCARDIAL INFARCTION",
#'   "PILOT01",    "AE", "18-1066",      2, "2013-07-18",             "AGITATION",
#'   "PILOT01",    "AE", "18-1066",      4, "2013-07-30",          "INFLAMMATION",
#'   "PILOT01",    "AE", "18-1066",      1, "2013-07-16",               "SYNCOPE",
#'   "PILOT01",    "AE", "18-1066",      3, "2013-07-30",               "SYNCOPE"
#' )
#'
#'
#' ex_single <- derive_vars_dtm(
#'   ex_single,
#'   dtc = EXSTDTC,
#'   new_vars_prefix = "EXST",
#'   flag_imputation = "none"
#' )
#'
#' adae <- ae %>%
#'   derive_vars_dtm(
#'     dtc = AESTDTC,
#'     new_vars_prefix = "AST",
#'     highest_imputation = "M"
#'   )
#'
#' adae %>%
#'   derive_var_last_dose_grp(
#'     dataset_ex = ex_single,
#'     filter_ex = (EXDOSE > 0 | (EXDOSE == 0 & grepl("PLACEBO", EXTRT))) &
#'       !is.na(EXSTDTM),
#'     by_vars = exprs(STUDYID, USUBJID),
#'     dose_date = EXSTDTM,
#'     new_var = LDGRP,
#'     grp_brks = c(0, 20, 40, 60),
#'     grp_lbls = c("Low", "Medium", "High"),
#'     include_lowest = TRUE,
#'     right = TRUE,
#'     dose_var = EXDOSE,
#'     analysis_date = ASTDTM,
#'     traceability_vars = exprs(LDOSEDOM = "EX", LDOSESEQ = EXSEQ, LDOSEVAR = "EXENDTC")
#'   ) %>%
#'   select(USUBJID, LDGRP, LDOSEDOM, LDOSESEQ, LDOSEVAR)
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
  filter_ex <- assert_filter_cond(enexpr(filter_ex), optional = TRUE)
  by_vars <- assert_vars(by_vars)
  dose_date <- assert_symbol(enexpr(dose_date))
  analysis_date <- assert_symbol(enexpr(analysis_date))
  single_dose_condition <- assert_filter_cond(enexpr(single_dose_condition))
  new_var <- assert_symbol(enexpr(new_var))
  dose_var <- assert_symbol(enexpr(dose_var))

  derive_vars_last_dose(
    dataset = dataset,
    dataset_ex = dataset_ex,
    filter_ex = !!filter_ex,
    by_vars = by_vars,
    dose_id = dose_id,
    dose_date = !!dose_date,
    analysis_date = !!analysis_date,
    single_dose_condition = !!single_dose_condition,
    new_vars = exprs(!!dose_var),
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
