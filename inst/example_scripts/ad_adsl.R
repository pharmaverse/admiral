# Name: ADSL
#
# Label: Subject Level Analysis Dataset
#
# Input: dm, ex, ds
#

library(dplyr)
library(lubridate)
library(admiral)

# read in predecessor variables from dm
data("dm")
data("ex")
data("ds")

# derive treatment variables (TRT01P, TRT01A)
adsl <- dm %>%
  mutate(
    TRT01P = ARMCD, TRT01A = ARMCD
  ) %>%

  # derive treatment start date (TRTSDTM, TRTSDT)
  derive_var_trtsdtm(dataset_ex = ex) %>%
  mutate(TRTSDT = date(TRTSDTM)) %>%

  # derive treatment end date (TRTEDTM, TRTEDT)
  derive_var_trtedtm(dataset_ex = ex) %>%
  mutate(TRTEDT = date(TRTEDTM)) %>%

  # derive treatment duration (TRTDURD)
  derive_var_trtdurd() %>%

  # Disposition date - option 1
  derive_disposition_dt(
    dataset_ds = ds,
    new_var = RFICDT,
    dtc = DSSTDTC,
    filter = expr(DSCAT == "PROTOCOL MILESTONE" & DSDECOD == "INFORMED CONSENT OBTAINED")
  ) %>%
  derive_disposition_dt(
    dataset_ds = ds,
    new_var = FRVDT,
    dtc = DSSTDTC,
    filter = expr(DSCAT == "OTHER EVENT" & DSDECOD == "FINAL RETRIEVAL VISIT")
  ) %>%

  # Option 2
  # derive study completion/discontinuation variables
  derive_merged_vars(
    dataset_add = ds,
    filter_add = exprs(DSCAT == "DISPOSITION EVENT"),
    new_vars = exprs(
      EOSDT = convert_dtc_to_dt(impute_dtc(DSSTDTC, date_imputation = "FIRST")),
      EOSSTT = if_else(DSDECOD == "COMPLETED", "COMPLETED", "DISCONTINUED"),
      DCSREAS = if_else(DSDECOD != "COMPLETED", DSDECOD, "")
    )
  )

save(adsl, file = "data/adsl.rda", compress = TRUE)
