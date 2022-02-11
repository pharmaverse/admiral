# Name: ADSL
#
# Label: Subject Level Analysis Dataset
#
# Input: dm, ex, ds
library(admiral)
library(admiral.test) # Contains example datasets from the CDISC pilot project
library(dplyr)
library(lubridate)

# ---- Load source datasets ----

# Use e.g. haven::read_sas to read in .sas7bdat, or other suitable functions
# as needed and assign to the variables below.
# For illustration purposes read in admiral test data

data("dm")
data("ds")
data("ex")
data("ae")
data("lb")

dm <- convert_blanks_to_na(dm)
ds <- convert_blanks_to_na(ds)
ex <- convert_blanks_to_na(ex)
ae <- convert_blanks_to_na(ae)
lb <- convert_blanks_to_na(lb)

# ---- User defined functions ----

# Here are some examples of how you can create your own functions that
#  operates on vectors, which can be used in `mutate`.

# Grouping
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

# ---- Derivations ----

adsl <- dm %>%
  # derive treatment variables (TRT01P, TRT01A)
  mutate(TRT01P = ARMCD, TRT01A = ARMCD) %>%

  # derive treatment start date (TRTSDTM)
  derive_var_trtsdtm(dataset_ex = ex) %>%

  # derive treatment end date (TRTEDTM)
  derive_var_trtedtm(dataset_ex = ex) %>%

  # Derive treatment end/start date TRTSDT/TRTEDT
  derive_vars_dtm_to_dt(vars(TRTSDTM, TRTEDTM)) %>%

  # derive treatment duration (TRTDURD)
  derive_var_trtdurd() %>%

  # Disposition dates, status
  # Screen fail date
  derive_var_disposition_dt(
    dataset_ds = ds,
    new_var = SCRFDT,
    dtc = DSSTDTC,
    filter = DSCAT == "DISPOSITION EVENT" & DSDECOD == "SCREEN FAILURE"
  ) %>%

  derive_var_disposition_dt(
    dataset_ds = ds,
    new_var = EOSDT,
    dtc = DSSTDTC,
    filter = DSCAT == "DISPOSITION EVENT" & DSDECOD != "SCREEN FAILURE"
  ) %>%

  # EOS status
  derive_var_disposition_status(
    dataset_ds = ds,
    new_var = EOSSTT,
    status_var = DSDECOD,
    format_new_var = format_eoxxstt,
    filter = DSCAT == "DISPOSITION EVENT"
  ) %>%

  # Last retrieval date
  derive_var_disposition_dt(
    dataset_ds = ds,
    new_var = FRVDT,
    dtc = DSSTDTC,
    filter = DSCAT == "OTHER EVENT" & DSDECOD == "FINAL RETRIEVAL VISIT"
  ) %>%

  # Death date - impute partial date to first day/month
  derive_vars_dt(
    new_vars_prefix = "DTH",
    dtc = DTHDTC,
    date_imputation = "FIRST"
  ) %>%

  # Relative Day of Death
  derive_vars_duration(
    new_var = DTHADY,
    start_date = TRTSDT,
    end_date = DTHDT
  ) %>%

  # Elapsed Days from Last Dose to Death
  derive_vars_duration(
    new_var = LDDTHELD,
    start_date = TRTEDT,
    end_date = DTHDT,
    add_one = FALSE
  )

# Last known alive date
ae_start <- lstalvdt_source(
  dataset_name = "ae",
  date = AESTDTC,
  date_imputation = "first"
)
ae_end <- lstalvdt_source(
  dataset_name = "ae",
  date = AEENDTC,
  date_imputation = "first"
)
lb_date <- lstalvdt_source(
  dataset_name = "lb",
  date = LBDTC,
  filter = nchar(LBDTC) >= 10
)
adsl_date <- lstalvdt_source(dataset_name = "adsl", date = TRTEDT)

adsl <- adsl %>%

  derive_var_lstalvdt(
    ae_start, ae_end, lb_date, adsl_date,
    source_datasets = list(ae = ae, lb = lb, adsl = adsl)
  ) %>%

  # Age group
  derive_agegr_fda(
    age_var = AGE,
    new_var = AGEGR1
  ) %>%

  # Groupings, populations and others variables
  mutate(
    RACEGR1 = format_racegr1(RACE),
    REGION1 = format_region1(COUNTRY),
    LDDTHGR1 = format_lddthgr1(LDDTHELD),
    FASFL = if_else(is.na(SCRFDT), "Y", "N"),
    SAFFL = if_else(!is.na(TRTSDTM), "Y", "N"),
    DTH30FL = if_else(LDDTHGR1 == "<= 30", "Y", NA_character_),
    DTHA30FL = if_else(LDDTHGR1 == "> 30", "Y", NA_character_),
    DTHB30FL = if_else(DTHDT <= TRTSDT + 30, "Y", NA_character_),
    DOMAIN = NULL
  )

# ---- Save output ----

save(adsl, file = "data/adsl.rda", compress = "bzip2")
