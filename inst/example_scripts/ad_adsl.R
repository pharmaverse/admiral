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

  # Disposition date
  derive_disposition_dt(
    dataset_ds = ds,
    new_var = EOSDT,
    dtc = DSSTDTC,
    filter = expr(DSCAT == "DISPOSITION EVENT")
  ) %>%
  derive_disposition_dt(
    dataset_ds = ds,
    new_var = FRVDT,
    dtc = DSSTDTC,
    filter = expr(DSCAT == "OTHER EVENT" & DSDECOD == "FINAL RETRIEVAL VISIT")
  ) %>%

save(adsl, file = "data/adsl.rda", compress = TRUE)
