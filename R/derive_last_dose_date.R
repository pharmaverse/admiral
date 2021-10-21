#' Derive Last Dose Date-Time
#'
#' Displays the start date or start datetime from the exposure dataset relative to an adverse event
#'
#' @inheritParams get_last_dose
#' @param new_var The output variable.
#' @param output_datetime  Display `new_var` as datetime or as date
#'
#' @details The last dose date
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
#' ae %>%
#'   head(100) %>%
#'   derive_last_dose_date(
#'     head(ex_single, 100),
#'     filter_ex = (EXDOSE > 0 | (EXDOSE == 0 & grepl("PLACEBO", EXTRT))) &
#'       nchar(EXENDTC) >= 10,
#'     dose_start = EXSTDTC,
#'     dose_end = EXENDTC,
#'     analysis_date = AESTDTC,
#'     dataset_seq_var = AESEQ,
#'     new_var = LDOSEDTM,
#'     output_datetime = TRUE,
#'     check_dates_only = FALSE,
#'     traceability_vars = dplyr::vars(LDOSEDOM = "EX", LDOSESEQ = EXSEQ, LDOSEVAR = "EXDOSE")
#'   ) %>%
#'   select(STUDYID, USUBJID, AESEQ, AESTDTC, LDOSEDTM)

derive_last_dose_date <- function(dataset,
                                  dataset_ex,
                                  filter_ex,
                                  by_vars = vars(STUDYID, USUBJID),
                                  dose_start,
                                  dose_end,
                                  new_var,
                                  analysis_date,
                                  dataset_seq_var,
                                  output_datetime = TRUE,
                                  check_dates_only,
                                  traceability_vars ){

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

  res <- get_last_dose(dataset = dataset,
                dataset_ex = dataset_ex,
                filter_ex = !!filter_ex,
                by_vars = by_vars ,
                dose_start = !!dose_start,
                dose_end  = !!dose_end,
                analysis_date = !!analysis_date,
                dataset_seq_var = !!dataset_seq_var,
                check_dates_only = !!check_dates_only,
                traceability_vars = traceability_vars) %>%
  select(colnames(dataset), !!new_var := !!dose_start, !!!syms(trace_vars_str))

  # return either date or date-time variable
  if (!output_datetime) {
    res <- res %>%  mutate(!!new_var := as.Date(!!new_var))}
  else {
    res <- res %>% mutate(!!new_var := as.POSIXct(as.character(!!new_var)))}

  return(res)
}

