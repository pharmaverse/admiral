# Name: ADSL
#
# Label: Subject Level Analysis Dataset
#
# Input: dm, ex, ds
#

library(dplyr)
library(lubridate)
library(admiral)

data("dm")
data("ex")
data("ds")


# User defined functions
# Grouping
format_agegr1 <- function(x) {
  case_when(
    !is.na(x) & x < 18 ~ "< 18",
    x >= 18 & x < 65 ~ "18 - 65",
    x >= 65 ~ ">= 65"
  )
}
format_racegr1 <- function(x) {
  case_when(
    !is.na(x) & x == "WHITE" ~ "White",
    !is.na(x) & x != "WHITE" ~ "Non-white",
    TRUE ~ "Missing"
  )
}
format_region1 <- function(x) {
  case_when(
    x %in% c("CAN", "USA") ~ "NA",
    !is.na(x) ~ "RoW",
    TRUE ~ "Missing"
  )
}
format_lddthgr1 <- function(x) {
  case_when(
    !is.na(x) & x <= 30 ~ "<= 30",
    x > 30 ~ "> 30",
    TRUE ~ NA_character_
  )
}

# EOSSTT mapping
format_eoxxstt <- function(x) {
  case_when(
    x %in% c("COMPLETED") ~ "COMPLETED",
    !(x %in% c("COMPLETED", "SCREEN FAILURE")) & !is.na(x) ~ "DISCONTINUED",
    x %in% c("SCREEN FAILURE") ~ NA_character_,
    TRUE ~ "ONGOING"
  )
}
# end of user defined functions


adsl <- dm %>%
  # derive treatment variables (TRT01P, TRT01A)
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

  # Disposition dates, status
  # Screen fail date
  derive_disposition_dt(
    dataset_ds = ds,
    new_var = SCRFDT,
    dtc = DSSTDTC,
    filter = expr(DSCAT == "DISPOSITION EVENT" & DSDECOD == "SCREEN FAILURE")
  ) %>%

  # EOS date
  derive_disposition_dt(
    dataset_ds = ds,
    new_var = EOSDT,
    dtc = DSSTDTC,
    filter = expr(DSCAT == "DISPOSITION EVENT" & DSDECOD != "SCREEN FAILURE")
  ) %>%

  # EOS status
  derive_disposition_status(
    dataset_ds = ds,
    new_var = EOSSTT,
    status_var = DSDECOD,
    format_new_var = format_eoxxstt,
    filter = expr(DSCAT == "DISPOSITION EVENT")
  ) %>%

  # study specific variables

  # Last retrieval date
  derive_disposition_dt(
    dataset_ds = ds,
    new_var = FRVDT,
    dtc = DSSTDTC,
    filter = expr(DSCAT == "OTHER EVENT" & DSDECOD == "FINAL RETRIEVAL VISIT")
  ) %>%
  # Death date - impute partial date to first day/month
  derive_vars_dt(
    new_vars_prefix = "DTH",
    dtc = DTHDTC,
    date_imputation = "FIRST"
  ) %>%

  # Relative Day of Death
  derive_duration(
    new_var = DTHADY,
    start_date = TRTSDT,
    end_date = DTHDT
  ) %>%
  # Elapsed Days from Last Dose to Death
  derive_duration(
    new_var = LDDTHELD,
    start_date = TRTEDT,
    end_date = DTHDT,
    add_one = FALSE
  ) %>%

  # Groupings, populations and others variables
  mutate(
    AGEGR1 = format_agegr1(AGE),
    RACEGR1 = format_racegr1(RACE),
    REGION1 = format_region1(COUNTRY),
    LDDTHGR1 = format_lddthgr1(LDDTHELD),
    FASFL = if_else(is.na(SCRFDT), "Y", "N"),
    SAFFL = if_else(!is.na(TRTSDTM), "Y", "N"),
    DTH30FL = if_else(LDDTHGR1 == "<= 30", "Y", NA_character_),
    DTHA30FL = if_else(LDDTHGR1 == "> 30", "Y", NA_character_),
    DTHB30FL = if_else(DTHDT <= TRTSDT + 30, "Y", NA_character_)
  )

save(adsl, file = "data/adsl.rda", compress = TRUE)
