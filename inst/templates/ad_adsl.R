# Name: ADSL
#
# Label: Subject Level Analysis Dataset
#
# Input: dm, ex, ds
library(admiral)
library(admiraltest) # Contains example datasets from the CDISC pilot project
library(dplyr)
library(lubridate)
library(stringr)

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
  mutate(TRT01P = ARM, TRT01A = ACTARM) %>%

  # derive treatment start date (TRTSDTM)
  derive_vars_merged_dtm(
    dataset_add = ex,
    filter_add = (EXDOSE > 0 |
                    (EXDOSE == 0 &
                       str_detect(EXTRT, "PLACEBO"))) & nchar(EXSTDTC) >= 10,
    new_vars_prefix = "TRTS",
    dtc = EXSTDTC,
    order = vars(TRTSDTM, EXSEQ),
    mode = "first",
    by_vars = vars(STUDYID, USUBJID)
  ) %>%

  # derive treatment end date (TRTEDTM)
  derive_vars_merged_dtm(
    dataset_add = ex,
    filter_add = (EXDOSE > 0 |
                    (EXDOSE == 0 &
                       str_detect(EXTRT, "PLACEBO"))) & nchar(EXENDTC) >= 10,
    new_vars_prefix = "TRTE",
    dtc = EXENDTC,
    time_imputation = "last",
    order = vars(TRTEDTM, EXSEQ),
    mode = "last",
    by_vars = vars(STUDYID, USUBJID)
  ) %>%

  # Derive treatment end/start date TRTSDT/TRTEDT
  derive_vars_dtm_to_dt(source_vars = vars(TRTSDTM, TRTEDTM)) %>%

  # derive treatment duration (TRTDURD)
  derive_var_trtdurd() %>%

  # Disposition dates, status
  # Screen fail date
  derive_vars_merged_dt(
    dataset_add = ds,
    by_vars = vars(STUDYID, USUBJID),
    new_vars_prefix = "SCRF",
    dtc = DSSTDTC,
    filter_add = DSCAT == "DISPOSITION EVENT" & DSDECOD == "SCREEN FAILURE"
  ) %>%

  derive_vars_merged_dt(
    dataset_add = ds,
    by_vars = vars(STUDYID, USUBJID),
    new_vars_prefix = "EOS",
    dtc = DSSTDTC,
    filter_add = DSCAT == "DISPOSITION EVENT" & DSDECOD != "SCREEN FAILURE"
  ) %>%

  # EOS status
  derive_var_disposition_status(
    dataset_ds = ds,
    new_var = EOSSTT,
    status_var = DSDECOD,
    format_new_var = format_eoxxstt,
    filter_ds = DSCAT == "DISPOSITION EVENT"
  ) %>%

  # Last retrieval date
  derive_vars_merged_dt(
    dataset_add = ds,
    by_vars = vars(STUDYID, USUBJID),
    new_vars_prefix = "FRV",
    dtc = DSSTDTC,
    filter_add = DSCAT == "OTHER EVENT" & DSDECOD == "FINAL RETRIEVAL VISIT"
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
ae_start <- date_source(
  dataset_name = "ae",
  date = AESTDTC,
  date_imputation = "first"
)
ae_end <- date_source(
  dataset_name = "ae",
  date = AEENDTC,
  date_imputation = "first"
)
lb_date <- date_source(
  dataset_name = "lb",
  date = LBDTC,
  filter = nchar(LBDTC) >= 10
)
adsl_date <- date_source(
  dataset_name = "adsl",
  date = TRTEDT
)

adsl <- adsl %>%

  derive_var_extreme_dt(
    new_var = LSTALVDT,
    ae_start, ae_end, lb_date, adsl_date,
    source_datasets = list(ae = ae, lb = lb, adsl = adsl),
    mode = "last"
  ) %>%

  # Age group
  derive_var_agegr_fda(
    age_var = AGE,
    new_var = AGEGR1
  ) %>%

  # Safety population
  derive_var_merged_exist_flag(
    dataset_add = ex,
    by_vars = vars(STUDYID, USUBJID),
    new_var = SAFFL,
    condition = (EXDOSE > 0 | (EXDOSE == 0 & str_detect(EXTRT, "PLACEBO")))
  ) %>%

  # Groupings and others variables
  mutate(
    RACEGR1 = format_racegr1(RACE),
    REGION1 = format_region1(COUNTRY),
    LDDTHGR1 = format_lddthgr1(LDDTHELD),
    DTH30FL = if_else(LDDTHGR1 == "<= 30", "Y", NA_character_),
    DTHA30FL = if_else(LDDTHGR1 == "> 30", "Y", NA_character_),
    DTHB30FL = if_else(DTHDT <= TRTSDT + 30, "Y", NA_character_),
    DOMAIN = NULL
  )

# ---- Save output ----

dir <- tempdir() # Change to whichever directory you want to save the dataset in
save(adsl, file = file.path(dir, "adsl.rda"), compress = "bzip2")
