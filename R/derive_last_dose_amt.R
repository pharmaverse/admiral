#' Derive Last Dose Amount
#'
#' @inheritParams derive_last_dose
#' @param new_var The output variable.
#' @param dose_var The source dose amount variable. Defaults to `EXDOSE`.
#'
#' @details The last dose amount is derived as the dose amount where the `end_date` is lower to or
#'    equal to the `analysis_date` per `by_vars` and `dataset_seq_var`.
#'
#' @return Input dataset with additional column `new_var`.
#'
#' @author Annie Yang
#'
#' @keywords adam derivation
#'
#' @export
#'
#' @seealso [derive_last_dose()]
#'
#' @examples
#' library(dplyr, warn.conflicts = FALSE)
#' library(cdiscpilot)
#' data(ae)
#' data(ex_single)
#'
#' ae %>%
#'   head(100) %>%
#'   derive_last_dose_amt(
#'     head(ex_single, 100),
#'     filter_ex = (EXDOSE > 0 | (EXDOSE == 0 & grepl("PLACEBO", EXTRT))) &
#'       nchar(EXENDTC) >= 10,
#'     ex_keep_vars = vars(EXSTDTC, EXENDTC, EXDOSE, EXTRT, EXSEQ, VISIT),
#'     dose_start = EXSTDTC,
#'     dose_end = EXENDTC,
#'     analysis_date = AESTDTC,
#'     dataset_seq_var = AESEQ,
#'     new_var = LDOSE,
#'     dose_var = EXDOSE,
#'     check_dates_only = FALSE
#'   ) %>%
#'   select(STUDYID, USUBJID, AESEQ, AESTDTC, LDOSE)
#'
#' # or with traceability variables
#' ae %>%
#'   head(100) %>%
#'   derive_last_dose_amt(
#'     head(ex_single, 100),
#'     filter_ex = (EXDOSE > 0 | (EXDOSE == 0 & grepl("PLACEBO", EXTRT))) &
#'       nchar(EXENDTC) >= 10,
#'     ex_keep_vars = vars(EXSTDTC, EXENDTC, EXDOSE, EXTRT, EXSEQ, VISIT),
#'     dose_start = EXSTDTC,
#'     dose_end = EXENDTC,
#'     analysis_date = AESTDTC,
#'     dataset_seq_var = AESEQ,
#'     new_var = LDOSE,
#'     dose_var = EXDOSE,
#'     check_dates_only = FALSE,
#'     traceability_vars = dplyr::vars(LDOSEDOM = "EX", LDOSESEQ = EXSEQ, LDOSEVAR = "EXDOSE")
#'   ) %>%
#'   select(STUDYID, USUBJID, AESEQ, AESTDTC, LDOSEDOM, LDOSESEQ, LDOSEVAR, LDOSE)


derive_last_dose_amt <- function(dataset,
                                 dataset_ex,
                                 filter_ex = NULL,
                                 by_vars = vars(STUDYID, USUBJID),
                                 dose_id = vars(),
                                 ex_keep_vars = NULL,
                                 dose_start,
                                 dose_end,
                                 analysis_date,
                                 dataset_seq_var,
                                 new_var,
                                 dose_var = EXDOSE,
                                 check_dates_only = FALSE,
                                 traceability_vars = NULL) {

  filter_ex <- assert_filter_cond(enquo(filter_ex), optional = TRUE)
  by_vars <- assert_vars(by_vars)
  dose_start <- assert_symbol(enquo(dose_start))
  dose_end <- assert_symbol(enquo(dose_end))
  analysis_date <- assert_symbol(enquo(analysis_date))
  dataset_seq_var <- assert_symbol(enquo(dataset_seq_var))
  new_var <- assert_symbol(enquo(new_var))
  dose_var <- assert_symbol(enquo(dose_var))
  trace_vars_str <- names(traceability_vars)

  derive_last_dose(dataset = dataset,
                dataset_ex = dataset_ex,
                filter_ex = !!filter_ex,
                by_vars = by_vars,
                dose_id = dose_id,
                ex_keep_vars = ex_keep_vars,
                dose_start = !!dose_start,
                dose_end = !!dose_end,
                analysis_date = !!analysis_date,
                dataset_seq_var = !!dataset_seq_var,
                check_dates_only = check_dates_only,
                traceability_vars = traceability_vars) %>%
    mutate(!!new_var := !!dose_var) %>%
    select(colnames(dataset), !!!syms(trace_vars_str), !!new_var)
}
