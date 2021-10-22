#' Derive Last Dose with user-defined groupings
#'
#' @inheritParams get_last_dose
#' @param grp_var The output variable defined by the user.
#' @param grp_brks User supplied breaks to apply to groups
#' @param grp_lbls User supplied labels to apply to groups
#' @param dose_var The source dose amount variable. Defaults to `EXDOSE`.
#' @param include.lowest logical, indicating if a value equal to the lowest
#' (or highest, for right = FALSE) ‘breaks’ value should be included.
#' @param right logical, indicating if the intervals should be closed on the right
#' (and open on the left) or vice versa
#'
#' @details This function brings in two datasets (e.g. adex and adae), finds the most recent
#' adverse event and then finds the most recent last dose.  Users can supply custom grouping breaks
#' and grouping labels for EXDOSE exploration.
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
#' @examples
#' library(dplyr, warn.conflicts = FALSE)
#' data(ae)
#' data(ex_single)
#'
#' ex_single_new <- ex_single %>%
#'  mutate(EXDOSE_ex = sample(0:84, 22439, replace=TRUE)) %>%
#'  select(-EXDOSE, "EXDOSE" = EXDOSE_ex)
#'
#' ae %>%
#'    head(100) %>%
#'    derive_last_dose_grp(
#'    head(ex_single_new, 100),
#'    filter_ex = (EXDOSE > 0 | (EXDOSE == 0 & grepl("PLACEBO", EXTRT))) &
#'    nchar(EXENDTC) >= 10,
#'    by_vars = vars(STUDYID, USUBJID),
#'    dose_start = EXSTDTC,
#'    dose_end = EXENDTC,
#'    grp_var = LDGRP,
#'    grp_brks = c(1, 40, 60, 100),
#'    grp_lbls = c("Low", "Medium", "High"),
#'    include.lowest = TRUE,
#'    right = TRUE,
#'    dose_var = EXDOSE,
#'    analysis_date = AESTDTC,
#'    dataset_seq_var = AESEQ,
#'    check_dates_only = FALSE,
#'    traceability_vars = NULL
#'    ) %>%
#'    select(USUBJID, LDGRP)
#'
#' # or with traceability variables
#' ae %>%
#'    head(100) %>%
#'    derive_last_dose_grp(
#'    head(ex_single_new, 100),
#'    filter_ex = (EXDOSE > 0 | (EXDOSE == 0 & grepl("PLACEBO", EXTRT))) &
#'    nchar(EXENDTC) >= 10,
#'    by_vars = vars(STUDYID, USUBJID),
#'    dose_start = EXSTDTC,
#'    dose_end = EXENDTC,
#'    grp_var = LDGRP,
#'    grp_brks = c(1, 40, 60, 100),
#'    grp_lbls = c("Low", "Medium", "High"),
#'    include.lowest = TRUE,
#'    right = TRUE,
#'    dose_var = EXDOSE,
#'    analysis_date = AESTDTC,
#'    dataset_seq_var = AESEQ,
#'    check_dates_only = FALSE,
#'    traceability_vars = dplyr::vars(LDOSEDOM = "EX", LDOSESEQ = EXSEQ, LDOSEVAR = "EXENDTC")
#'    )%>%
#'    select(USUBJID, LDGRP, LDOSEDOM, LDOSESEQ, LDOSEVAR)


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
                                 include.lowest = TRUE,
                                 right = TRUE,
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
    mutate(!!grp_var :=
             as.character(
               cut(
                   EXDOSE,
                   breaks = !!grp_brks,
                   include.lowest = !!include.lowest,
                   right = !!right,
                   labels = !!grp_lbls)))  %>%
    select(colnames(dataset), !!!syms(trace_vars_str), !!grp_var)
}
