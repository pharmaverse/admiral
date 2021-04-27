
derive_var_brthdt <- function(dataset, date_imputation="MID") {
  derive_vars_dt(dataset,
    new_vars_prefix = "BRTH",
    dtc = BRTHDTC,
    date_imputation = date_imputation
  )
}
