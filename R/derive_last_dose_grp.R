#' Derive Last Dose with user-defined groupings
#'
#' @inheritParams get_last_dose
#' @param grp_var The output variable.
#' @param grp_brk Breaks to apply to groups
#' @param grp_lbl Labels to apply to groups
#' @param dose_var The source dose amount variable. Defaults to `EXDOSE`.
#'
#' @details Holding this space
#'
#' @return Input dataset with additional column `new_var`.
#'
#' @author Ben Straub
#'
#' @keywords adae derivation
#'
#' @export
#'
#' @seealso [get_last_dose()]
#'


derive_last_dose_grp <- function(dataset,
                                 dataset_ex,
                                 filter_ex = NULL,
                                 by_vars = vars(STUDYID, USUBJID),
                                 dose_start,
                                 dose_end,
                                 analysis_date,
                                 dataset_seq_var,
                                 grp_var,
                                 grp_brks,
                                 grp_lbls,
                                 dose_var = EXDOSE,
                                 check_dates_only = FALSE,
                                 traceability_vars = NULL) {

  filter_ex <- assert_filter_cond(enquo(filter_ex), optional = TRUE)
  by_vars <- assert_vars(by_vars)
  dose_start <- assert_symbol(enquo(dose_start))
  dose_end <- assert_symbol(enquo(dose_end))
  analysis_date <- assert_symbol(enquo(analysis_date))
  dataset_seq_var <- assert_symbol(enquo(dataset_seq_var))
  grp_var <- assert_symbol(enquo(grp_var))
  #grp_brks <- assert_symbol(enquo(grp_brks))
  #grp_lbls <- assert_symbol(enquo(grp_lbls))
  dose_var <- assert_symbol(enquo(dose_var))
  trace_vars_str <- names(traceability_vars)

  get_last_dose(dataset = dataset,
                dataset_ex = dataset_ex,
                filter_ex = !!filter_ex,
                by_vars = by_vars,
                dose_start = !!dose_start,
                dose_end = !!dose_end,
                analysis_date = !!analysis_date,
                dataset_seq_var = !!dataset_seq_var,
                check_dates_only = !!check_dates_only,
                traceability_vars = traceability_vars) %>%
    mutate(!!grp_var := cut(
                   EXDOSE,
                   breaks = !!grp_brks,
                   right = T,
                   include.lowest = T,
                   labels = !!grp_lbls)) %>%
    select(colnames(dataset), !!!syms(trace_vars_str), !!grp_var)
}
