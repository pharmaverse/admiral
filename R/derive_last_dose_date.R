#' Derive Last Dose Date-Time
#'
#' Displays the start date or start datetime of the last dose with respect to the most recent adverse event
#'
#' @inheritParams derive_last_dose
#' @param new_var The output variable defined by the user.
#' @param output_datetime  Display `new_var` as datetime or as date only.
#'
#' @details The datasets `dataset` and `dataset_ex` are joined using `by_vars`. The last dose date
#' is the maximum date where `dose_end` is lower to or equal to `analysis_date`, subject to both
#' date values are non-NA. The last dose date is derived per `by_vars` and `dataset_seq_var`, and
#' is appended to the `dataset` and returned to the user as the `new_var`.
#'
#' @return Input dataset with additional column `new_var`.
#'
#' @author Ben Straub
#'
#' @keywords adam derivation
#'
#' @export
#'
#' @seealso [derive_last_dose()]
#'
#' @examples
#' library(dplyr, warn.conflicts = FALSE)
#' data(ae)
#' data(ex_single)
#'
#' ae %>%
#'   head(100) %>%
#'   derive_last_dose_date(
#'     head(ex_single, 100),
#'     filter_ex = (EXDOSE > 0 | (EXDOSE == 0 & grepl("PLACEBO", EXTRT))) &
#'       nchar(EXENDTC) >= 10,
#'     ex_keep_vars = vars(EXSTDTC, EXENDTC, EXDOSE, EXTRT, EXSEQ, VISIT),
#'     dose_start = EXSTDTC,
#'     dose_end = EXENDTC,
#'     analysis_date = AESTDTC,
#'     dataset_seq_var = AESEQ,
#'     new_var = LDOSEDTM,
#'     output_datetime = TRUE,
#'     check_dates_only = FALSE,
#'     traceability_vars = dplyr::vars(LDOSEDOM = "EX", LDOSESEQ = EXSEQ, LDOSEVAR = "EXDOSE")
#'   ) %>%
#'   select(STUDYID, USUBJID, AESEQ, AESTDTC, LDOSEDOM, LDOSESEQ, LDOSEVAR, LDOSEDTM)

derive_last_dose_date <- function(dataset,
                                  dataset_ex,
                                  filter_ex = NULL,
                                  by_vars = vars(STUDYID, USUBJID),
                                  dose_id = vars(),
                                  ex_keep_vars = NULL,
                                  dose_start,
                                  dose_end,
                                  new_var,
                                  analysis_date,
                                  dataset_seq_var,
                                  output_datetime = TRUE,
                                  check_dates_only = FALSE,
                                  traceability_vars = NULL){

  # assert functions found in assertions.R
  filter_ex <- assert_filter_cond(enquo(filter_ex), optional = TRUE)
  by_vars <- assert_vars(by_vars)
  dose_start <- assert_symbol(enquo(dose_start))
  dose_end <- assert_symbol(enquo(dose_end))
  analysis_date <- assert_symbol(enquo(analysis_date))
  dataset_seq_var <- assert_symbol(enquo(dataset_seq_var))
  new_var <- assert_symbol(enquo(new_var))
  assert_logical_scalar(check_dates_only)
  assert_logical_scalar(output_datetime)
  trace_vars_str <- names(traceability_vars)
  assert_data_frame(dataset, quo_c(by_vars, analysis_date, dataset_seq_var))
  assert_data_frame(dataset_ex, quo_c(by_vars, dose_start, dose_end))

  res <- derive_last_dose(dataset = dataset,
                dataset_ex = dataset_ex,
                filter_ex = !!filter_ex,
                by_vars = by_vars ,
                dose_id = dose_id,
                ex_keep_vars = ex_keep_vars,
                dose_start = !!dose_start,
                dose_end  = !!dose_end,
                analysis_date = !!analysis_date,
                dataset_seq_var = !!dataset_seq_var,
                check_dates_only = check_dates_only,
                traceability_vars = traceability_vars) %>%
  select(colnames(dataset), !!new_var := !!dose_start, !!!syms(trace_vars_str))

  # return either date or date-time variable
  if (!output_datetime) {
    res <- res %>%  mutate(!!new_var := as.Date(!!new_var))}
  else {
    res <- res %>% mutate(!!new_var := as.POSIXct(as.character(!!new_var)))}

  return(res)
}

