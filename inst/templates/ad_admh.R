# Name: ADMH
#
# Label: Medical History Analysis Dataset
#
# Input: mh, adsl
library(admiral)
library(admiral.test) # Contains example datasets from the CDISC pilot project
library(dplyr)
library(lubridate)

# ---- Load source datasets ----

# Use e.g. haven::read_sas to read in .sas7bdat, or other suitable functions
# as needed and assign to the variables below.
# For illustration purposes read in admiral test data
data("admiral_mh")
data("admiral_adsl")

adsl <- admiral_adsl
mh <- admiral_mh

mh <- convert_blanks_to_na(mh)

# ---- Derivations ----

# Get list of ADSL vars required for derivations
adsl_vars <- vars(TRTSDT, TRTEDT)

admh <- mh %>%
  # join adsl to mh
  derive_vars_merged(
    dataset_add = adsl,
    new_vars = adsl_vars,
    by = vars(STUDYID, USUBJID)
  ) %>%
  # derive analysis start time and flag
  derive_vars_dtm(
    dtc = MHSTDTC,
    new_vars_prefix = "AST",
    time_imputation = "first"
  ) %>%

  # derive analysis end time and flag
  # derive_vars_dtm(
  #   dtc = MHENDTC,
  #   new_vars_prefix = "AEN",
  #   date_imputation = "last",
  #   time_imputation = "last",
  #   max_dates = vars(DTHDT, EOSDT, TRTSDT)
  # ) %>%

  # derive analysis end/start date
  # derive_vars_dtm_to_dt(vars(ASTDTM, AENDTM)) %>%
  derive_vars_dtm_to_dt(vars(ASTDTM)) %>%

  # derive analysis start relative day and analysis end relative day
  derive_vars_dy(
    reference_date = TRTSDT,
    source_vars = vars(ASTDT)
  ) #%>%

  # derive_vars_dy(
  #   reference_date = TRTSDT,
  #   source_vars = vars(AENDT)
  # ) %>%

  # derive analysis flag
  # mutate(ANL01FL = ifelse(MHOCCUR != "N", "Y", NA_character_))
