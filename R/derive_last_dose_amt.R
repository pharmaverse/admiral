#' Derive Last Dose Amount
#'
#' Add a variable for dose amount from the last dose to the input dataset.
#'
#' @inheritParams derive_last_dose
#' @param new_var The new variable added to `dataset`.
#' @param dose_var The EX source dose amount variable. Defaults to `EXDOSE`.
#'
#' @details The last dose amount is derived as the dose amount where the maximum `dose_date` is
#' lower to or equal to the `analysis_date` per `by_vars` and `dataset_seq_var`.
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
#' library(admiral.test)
#' data(ae)
#' data(ex_single)
#'
#' ae %>%
#'   head(100) %>%
#'   derive_last_dose_amt(
#'     head(ex_single, 100),
#'     filter_ex = (EXDOSE > 0 | (EXDOSE == 0 & grepl("PLACEBO", EXTRT))) &
#'       nchar(EXENDTC) >= 10,
#'     dose_date = EXENDTC,
#'     analysis_date = AESTDTC,
#'     dataset_seq_var = AESEQ,
#'     single_dose_condition = (EXSTDTC == EXENDTC),
#'     new_var = LDOSE,
#'     dose_var = EXDOSE
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
#'     dose_date = EXENDTC,
#'     analysis_date = AESTDTC,
#'     dataset_seq_var = AESEQ,
#'     single_dose_condition = (EXSTDTC == EXENDTC),
#'     new_var = LDOSE,
#'     dose_var = EXDOSE,
#'     traceability_vars = dplyr::vars(LDOSEDOM = "EX", LDOSESEQ = EXSEQ, LDOSEVAR = "EXDOSE")
#'   ) %>%
#'   select(STUDYID, USUBJID, AESEQ, AESTDTC, LDOSEDOM, LDOSESEQ, LDOSEVAR, LDOSE)


derive_last_dose_amt <- function(dataset,
                                 dataset_ex,
                                 filter_ex = NULL,
                                 by_vars = vars(STUDYID, USUBJID),
                                 dose_id = vars(),
                                 dose_date,
                                 analysis_date,
                                 dataset_seq_var,
                                 single_dose_condition = (EXDOSFRQ =="ONCE"),
                                 new_var,
                                 dose_var = EXDOSE,
                                 traceability_vars = NULL) {

  filter_ex <- assert_filter_cond(enquo(filter_ex), optional = TRUE)
  by_vars <- assert_vars(by_vars)
  dose_id <- assert_vars(dose_id)
  dose_date <- assert_symbol(enquo(dose_date))
  analysis_date <- assert_symbol(enquo(analysis_date))
  dataset_seq_var <- assert_symbol(enquo(dataset_seq_var))
  single_dose_condition <- assert_filter_cond(enquo(single_dose_condition))
  new_var <- assert_symbol(enquo(new_var))
  dose_var <- assert_symbol(enquo(dose_var))

  trace_vars_str <- names(traceability_vars)

  derive_last_dose(dataset = dataset,
                dataset_ex = dataset_ex,
                filter_ex = !!filter_ex,
                by_vars = by_vars,
                dose_id = dose_id,
                dose_date = !!dose_date,
                analysis_date = !!analysis_date,
                dataset_seq_var = !!dataset_seq_var,
                single_dose_condition = !!single_dose_condition,
                ex_keep_vars = vars(!!dose_var),
                traceability_vars = traceability_vars) %>%
    select(colnames(dataset), !!!syms(trace_vars_str), !!new_var := !!dose_var)
}
