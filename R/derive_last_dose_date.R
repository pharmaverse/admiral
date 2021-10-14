

derive_last_dose_date <- function(dataset,
                                  dataset_ex,
                                  filter_ex,
                                  by_vars ,
                                  dose_start,
                                  dose_end,
                                  analysis_date,
                                  dataset_seq_var,
                                  check_dates_only ,
                                  traceability_vars ){

  filter_ex <- assert_filter_cond(enquo(filter_ex), optional = TRUE)
  by_vars <- assert_vars(by_vars)
  dose_start <- assert_symbol(enquo(dose_start))
  dose_end <- assert_symbol(enquo(dose_end))
  analysis_date <- assert_symbol(enquo(analysis_date))
  dataset_seq_var <- assert_symbol(enquo(dataset_seq_var))
  new_var <- assert_symbol(enquo(new_var))
  assert_logical_scalar(output_datetime)
  assert_logical_scalar(check_dates_only)
  stopifnot(is_quosures(traceability_vars) | is.null(traceability_vars))
  assert_data_frame(dataset, quo_c(by_vars, analysis_date, dataset_seq_var))
  assert_data_frame(dataset_ex, quo_c(by_vars, dose_start, dose_end))

  get_last_dose(dataset = dataset,
                dataset_ex = dataset_ex,
                filter_ex = !!filter_ex,
                by_vars = by_vars,
                dose_start = !!dose_start,
                dose_end = !!dose_end,
                analysis_date = !!analysis_date,
                dataset_seq_var = !!dataset_seq_var,
                check_dates_only = !!check_dates_only,
                traceability_vars = traceability_vars) %>% select(EXSTDTC)

}

