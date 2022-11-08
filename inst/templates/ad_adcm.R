# Name: ADCM
#
# Label: Concomitant Medications Analysis Dataset
#
# Input: cm, adsl
library(admiral)
library(admiral.test) # Contains example datasets from the CDISC pilot project
library(dplyr)
library(lubridate)

# Load source datasets ----

# Use e.g. haven::read_sas to read in .sas7bdat, or other suitable functions
# as needed and assign to the variables below.
# For illustration purposes read in admiral test data

data("admiral_cm")
data("admiral_adsl")

adsl <- admiral_adsl
cm <- admiral_cm

# When SAS datasets are imported into R using haven::read_sas(), missing
# character values from SAS appear as "" characters in R, instead of appearing
# as NA values. Further details can be obtained via the following link:
# https://pharmaverse.github.io/admiral/articles/admiral.html#handling-of-missing-values

cm <- convert_blanks_to_na(cm)

# Derivations ----

# Get list of ADSL vars required for derivations
adsl_vars <- vars(TRTSDT, TRTEDT, DTHDT, EOSDT, TRT01P, TRT01A)

adcm <- cm %>%
  # Join ADSL with CM (only ADSL vars required for derivations)
  derive_vars_merged(
    dataset_add = adsl,
    new_vars = adsl_vars,
    by = vars(STUDYID, USUBJID)
  ) %>%
  ## Derive analysis start time ----
  derive_vars_dtm(
    dtc = CMSTDTC,
    new_vars_prefix = "AST",
    highest_imputation = "M",
    min_dates = vars(TRTSDT)
  ) %>%
  ## Derive analysis end time ----
  derive_vars_dtm(
    dtc = CMENDTC,
    new_vars_prefix = "AEN",
    highest_imputation = "M",
    date_imputation = "last",
    time_imputation = "last",
    max_dates = vars(DTHDT, EOSDT)
  ) %>%
  ## Derive analysis end/start date -----
  derive_vars_dtm_to_dt(vars(ASTDTM, AENDTM)) %>%
  ## Derive analysis start relative day and analysis end relative day ----
  derive_vars_dy(
    reference_date = TRTSDT,
    source_vars = vars(ASTDT, AENDT)
  ) %>%
  ## Derive analysis duration (value and unit) ----
  derive_vars_duration(
    new_var = ADURN,
    new_var_unit = ADURU,
    start_date = ASTDT,
    end_date = AENDT,
    in_unit = "days",
    out_unit = "days",
    add_one = TRUE,
    trunc_out = FALSE
  )

## Derive flags ----
adcm <- adcm %>%
  # Derive On-Treatment flag
  # Set `span_period = "Y"` if you want occurrences that started prior to drug
  # intake and ongoing or ended after this time to be considered as on-treatment.
  derive_var_ontrtfl(
    start_date = ASTDT,
    end_date = AENDT,
    ref_start_date = TRTSDT,
    ref_end_date = TRTEDT
  ) %>%
  # Derive Pre-Treatment flag
  mutate(PREFL = if_else(ASTDT < TRTSDT, "Y", NA_character_)) %>%
  # Derive Follow-Up flag
  mutate(FUPFL = if_else(ASTDT > TRTEDT, "Y", NA_character_)) %>%
  # Derive ANL01FL
  # This variable is sponsor specific and may be used to indicate particular
  # records to be used in subsequent derivations or analysis.
  mutate(ANL01FL = if_else(ONTRTFL == "Y", "Y", NA_character_)) %>%
  # Derive 1st Occurrence of Preferred Term Flag
  restrict_derivation(
    derivation = derive_var_extreme_flag,
    args = params(
      new_var = AOCCPFL,
      by_vars = vars(USUBJID, CMDECOD),
      order = vars(ASTDTM, CMSEQ),
      mode = "first"
    ),
    filter = ANL01FL == "Y"
  )


## Derive APHASE and APHASEN Variable ----
# Other timing variable can be derived similarly.
# See also the "Visit and Period Variables" vignette
# (https://pharmaverse.github.io/admiral/articles/visits_periods.html)
adcm <- adcm %>%
  mutate(
    APHASE = case_when(
      PREFL == "Y" ~ "Pre-Treatment",
      ONTRTFL == "Y" ~ "On-Treatment",
      FUPFL == "Y" ~ "Follow-Up"
    ),
    APHASEN = case_when(
      PREFL == "Y" ~ 1,
      ONTRTFL == "Y" ~ 2,
      FUPFL == "Y" ~ 3
    )
  ) %>%
  # Assign TRTP/TRTA
  mutate(
    TRTP = TRT01P,
    TRTA = TRT01A
  )

# Join all ADSL with CM
adcm <- adcm %>%
  derive_vars_merged(
    dataset_add = select(adsl, !!!negate_vars(adsl_vars)),
    by_vars = vars(STUDYID, USUBJID)
  )


# Save output ----

dir <- tempdir() # Change to whichever directory you want to save the dataset in
saveRDS(adcm, file = file.path(dir, "adcm.rds"), compress = "bzip2")
