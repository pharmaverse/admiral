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
    TRT01P = ARMCD, TRT01A = ARMCD,
    # simulate date of birth
    BRTHDTC = if_else(SEX == "F", "1970", "1970-01-01")
  ) %>%

  # derive treatment start date (TRTSDTM, TRTSDT)
  derive_var_trtsdtm(dataset_ex = ex) %>%
  mutate(TRTSDT = date(TRTSDTM)) %>%

  # derive treatment end date (TRTEDTM, TRTEDT)
  derive_var_trtedtm(dataset_ex = ex) %>%
  mutate(TRTEDT = date(TRTEDTM)) %>%

  # derive treatment duration (TRTDURD)
  derive_var_trtdurd() %>%

<<<<<<< HEAD
  # Impute date of birth as needed
  derive_var_brthdt() %>%
  # Derive ade based on TRT starts
  derive_aage(start_date = BRTHDT, end_date = TRTSDT) %>%

  # Derive date of IC, rando, enrollement
  # Derived by default from dataset_ds.DSSTDTC
  # RFICDT
  # where DSCAT == "PROTOCOL MILESTONE" & startsWith(DSSCAT, "PROTOCOL") & DSDECOD == "INFORMED CONSENT OBTAINED"
  derive_var_rficdt(
    dataset_ds = ds,
    filter_ds = expr(DSCAT == "PROTOCOL MILESTONE" & DSDECOD == "INFORMED CONSENT OBTAINED")
  ) %>%

  # RANDDT
  # where DSCAT == "PROTOCOL MILESTONE" & DSDECOD == "RANDOMIZATION"
  derive_var_randdt(dataset_ds = ds) %>%

  # ENRLDT
  # where DSDECOD == "ENROLLED"
  derive_var_enrldt(dataset_ds = ds) %>%

  ## EOS date,day, status, reason
  # Derived by default from dataset_ds.DSSTDTC
  # where DSCAT == "DISPOSITION EVENT" & DSSCAT =="STUDY COMPLETION/EARLY DISCONTINUATION"
  derive_var_eosdt(
    dataset_ds = ds,
    filter_ds = expr(DSCAT == "DISPOSITION EVENT" & !(DSDECOD %in% c("COMPLETED", "SCREEN FAILURE")))
  ) %>%

  # EOS day
  derive_var_eosdy() %>%

  # Derived by default from dataset_ds.DSDECOD/DSTERM
  # where DSCAT == "DISPOSITION EVENT" & DSSCAT =="STUDY COMPLETION/EARLY DISCONTINUATION"
  derive_var_eosstt(
    dataset_ds = ds,
    filter_ds = expr(DSCAT == "DISPOSITION EVENT" & !(DSDECOD %in% c("SCREEN FAILURE")))
  ) %>%
  derive_var_dcsreas(
    dataset_ds = ds,
    filter_ds = expr(DSCAT == "DISPOSITION EVENT" & !(DSDECOD %in% c("SCREEN FAILURE")))
  ) %>%
  derive_var_dcsreasp(
    dataset_ds = ds,
    filter_ds = expr(DSCAT == "DISPOSITION EVENT" & !(DSDECOD %in% c("SCREEN FAILURE")))
  ) %>%

  # Option to discuss: 1 function can derive multiple variables vs 1 function per variable
  derive_disposition_dt(
    dataset_ds = ds,
    new_var = RFICDT2,
    dtc = DSSTDTC,
    filter = expr(DSCAT == "PROTOCOL MILESTONE" & DSDECOD == "INFORMED CONSENT OBTAINED")
  ) %>%
  derive_disposition_dt(
    dataset_ds = ds,
    new_var = FRVDT,
    dtc = DSSTDTC,
    filter = expr(DSCAT == "OTHER EVENT" & DSDECOD == "FINAL RETRIEVAL VISIT")
  )%>%
  derive_disposition_eoxxstt(
    dataset_ds = ds,
    new_var = EOS01STT,
    status_var = DSDECOD,
    filter = expr(DSCAT == "DISPOSITION EVENT" & !(DSDECOD %in% c("SCREEN FAILURE")))
  )

=======
  # derive study completion/discontinuation variables
  derive_merged_vars(
    dataset_add = ds,
    filter_add = exprs(DSCAT == "DISPOSITION EVENT"),
    new_vars = exprs(
      EOSDT = convert_dtc_to_dt(impute_dtc(DSSTDTC, date_imputation = "FIRST")),
      EOSSTT = if_else(DSDECOD == "COMPLETED", "COMPLETED", "DISCONTINUED"),
      DCSREAS = if_else(DSDECOD != "COMPLETED", DSDECOD, "")
    ))
>>>>>>> master

save(adsl, file = "data/adsl.rda", compress = TRUE)
