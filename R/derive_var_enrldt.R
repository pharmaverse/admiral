
derive_var_enrldt <- function(dataset,
                              dataset_ds,
                              filter_ds = expr(DSDECOD == "ENROLLED"),
                              date_imputation = NULL) {


  derive_merged_vars(
    dataset,
    dataset_add = dataset_ds,
    filter_add = filter_ds,
    new_vars = exprs(ENRLDT = convert_dtc_to_dt(impute_dtc(DSSTDTC, date_imputation = !!enquo(date_imputation))))
  )

}

