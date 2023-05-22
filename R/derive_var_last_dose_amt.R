#' Derive Last Dose Amount
#'
#' Add a variable for dose amount from the last dose to the input dataset.
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
#' @family der_gen
#' @keywords der_gen
#'
#' @export
#'
#' @seealso [derive_vars_last_dose()], [create_single_dose_dataset()]
#'
#' @examples
#' library(dplyr, warn.conflicts = FALSE)
#'
#' ex_single <- tribble(
#'   ~STUDYID,  ~DOMAIN,  ~USUBJID, ~EXSEQ,     ~EXENDTC,    ~EXTRT, ~EXDOSE, ~EXDOSFRQ,
#'   "PILOT01",    "EX", "18-1066",     1L, "2013-07-07",    "XANO",      54,    "ONCE",
#'   "PILOT01",    "EX", "18-1066",     2L, "2013-07-08",    "XANO",      54,    "ONCE",
#'   "PILOT01",    "EX", "18-1066",     3L, "2013-07-09",    "XANO",      54,    "ONCE",
#'   "PILOT01",    "EX", "18-1066",     4L, "2013-07-10",    "XANO",      54,    "ONCE",
#'   "PILOT01",    "EX", "18-1066",     5L, "2013-07-11",    "XANO",      54,    "ONCE",
#'   "PILOT01",    "EX", "18-1066",     6L, "2013-07-12",    "XANO",      54,    "ONCE",
#'   "PILOT01",    "EX", "18-1066",     7L, "2013-07-13",    "XANO",      54,    "ONCE",
#'   "PILOT01",    "EX", "18-1066",     8L, "2013-07-14",    "XANO",      54,    "ONCE",
#'   "PILOT01",    "EX", "18-1066",     9L, "2013-07-15",    "XANO",      54,    "ONCE",
#'   "PILOT01",    "EX", "18-1066",    10L, "2013-07-16",    "XANO",      54,    "ONCE",
#'   "PILOT01",    "EX", "10-1083",     1L, "2013-07-22", "PLACEBO",       0,    "ONCE",
#'   "PILOT01",    "EX", "10-1083",     2L, "2013-07-23", "PLACEBO",       0,    "ONCE",
#'   "PILOT01",    "EX", "10-1083",     3L, "2013-07-24", "PLACEBO",       0,    "ONCE",
#'   "PILOT01",    "EX", "10-1083",     4L, "2013-07-25", "PLACEBO",       0,    "ONCE",
#'   "PILOT01",    "EX", "10-1083",     5L, "2013-07-26", "PLACEBO",       0,    "ONCE",
#'   "PILOT01",    "EX", "10-1083",     6L, "2013-07-27", "PLACEBO",       0,    "ONCE",
#'   "PILOT01",    "EX", "10-1083",     7L, "2013-07-28", "PLACEBO",       0,    "ONCE",
#'   "PILOT01",    "EX", "10-1083",     8L, "2013-07-29", "PLACEBO",       0,    "ONCE",
#'   "PILOT01",    "EX", "10-1083",     9L, "2013-07-30", "PLACEBO",       0,    "ONCE",
#'   "PILOT01",    "EX", "10-1083",    10L, "2013-07-31", "PLACEBO",       0,    "ONCE",
#'   "PILOT01",    "EX", "10-1083",    11L, "2013-08-01", "PLACEBO",       0,    "ONCE"
#' )
#'
#'
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
#'   head(ex_single, 100),
#'   dtc = EXENDTC,
#'   new_vars_prefix = "EXEN",
#'   flag_imputation = "none"
#' )
#' adae <- ae %>%
#'   derive_vars_dtm(
#'     dtc = AESTDTC,
#'     new_vars_prefix = "AST",
#'     highest_imputation = "M"
#'   )
#' adae %>%
#'   derive_var_last_dose_amt(
#'     dataset_ex = ex_single,
#'     filter_ex = (EXDOSE > 0 | (EXDOSE == 0 & grepl("PLACEBO", EXTRT))) &
#'       !is.na(EXENDTM),
#'     dose_date = EXENDTM,
#'     analysis_date = ASTDTM,
#'     new_var = LDOSE,
#'     dose_var = EXDOSE
#'   ) %>%
#'   select(STUDYID, USUBJID, AESEQ, AESTDTC, LDOSE)
#' # or with traceability variables
#' adae %>%
#'   derive_var_last_dose_amt(
#'     dataset_ex = ex_single,
#'     filter_ex = (EXDOSE > 0 | (EXDOSE == 0 & grepl("PLACEBO", EXTRT))) &
#'       !is.na(EXENDTM),
#'     dose_date = EXENDTM,
#'     analysis_date = ASTDTM,
#'     new_var = LDOSE,
#'     dose_var = EXDOSE,
#'     traceability_vars = exprs(
#'       LDOSEDOM = "EX",
#'       LDOSESEQ = EXSEQ,
#'       LDOSEVAR = "EXDOSE"
#'     )
#'   ) %>%
#'   select(STUDYID, USUBJID, AESEQ, AESTDTC, LDOSEDOM, LDOSESEQ, LDOSEVAR, LDOSE)
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
  filter_ex <- assert_filter_cond(enexpr(filter_ex), optional = TRUE)
  by_vars <- assert_vars(by_vars)
  dose_id <- assert_vars(dose_id)
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
    new_vars = exprs(!!new_var := !!dose_var),
    traceability_vars = traceability_vars
  )
}
