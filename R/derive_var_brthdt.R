
derive_var_brthdt <- function(dataset) {
  derive_vars_dt(dataset,
    new_vars_prefix = "BRTH",
    dtc = BRTHDTC,
    date_imputation = "MID"
  )
}
