
derive_var_rficdt <- function(dataset,
                              dataset_ds,
                              filter_ds = expr(DSCAT == "PROTOCOL MILESTONE" & startsWith(DSSCAT, "PROTOCOL") & DSDECOD == "INFORMED CONSENT OBTAINED"),
                              date_imputation = NULL) {
  derive_merged_vars(
    dataset,
    dataset_add = dataset_ds,
    filter_add = filter_ds,
    #expr() does not work here... exprs() required?
    new_vars = exprs(RFICDT = convert_dtc_to_dt(impute_dtc(DSSTDTC, date_imputation = !!enquo(date_imputation))))
  )
}
