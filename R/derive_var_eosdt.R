
derive_var_eosdt <- function(dataset,
                              dataset_ds = ds,
                              filter_ds = exprs(DSCAT == "DISPOSITION EVENT" & DSSCAT =="STUDY COMPLETION/EARLY DISCONTINUATION"),
                              date_imputation = NULL
){

  derive_merged_vars(dataset,
                     dataset_add = dataset_ds,
                     filter_add = filter_ds,
                     by_vars = exprs(USUBJID),
                     new_vars = exprs(EOSDT = convert_dtc_to_dt(impute_dtc(DSSTDTC, date_imputation = !!enquo(date_imputation))))
  )
}




