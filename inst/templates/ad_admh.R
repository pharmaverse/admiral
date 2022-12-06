# Name: ADMH
#
# Label: Medical History Analysis Dataset
#
# Input: mh, adsl
library(admiral)
library(admiral.test) # Contains example datasets from the CDISC pilot project
library(dplyr)
library(lubridate)

# Load source datasets ----

# Use e.g. haven::read_sas to read in .sas7bdat, or other suitable functions
# as needed and assign to the variables below.
# For illustration purposes read in admiral test data
data("admiral_mh")
data("admiral_adsl")
data("queries_mh")

adsl <- admiral_adsl
mh <- admiral_mh

mh <- convert_blanks_to_na(mh)


# Look-up tables ----

# Creating a look-up table for assigning MHTERMN (for derivation of company specific variable)
# (this is set to align with the order of pre-printed terms on the CRF)
mhtermn_lookup <- tibble::tribble(
  ~MHTERM, ~MHTERMN,
  "ALZHEIMER'S DISEASE", 1
)

# Derivations ----

# Get list of ADSL vars required for derivations
adsl_vars <- vars(TRTSDT, TRTEDT, TRT01A, TRT01P, DTHDT, EOSDT)

admh <- mh %>%
  # join ADSL with MH
  derive_vars_merged(
    dataset_add = adsl,
    new_vars = adsl_vars,
    by_vars = vars(STUDYID, USUBJID)
  ) %>%
  ## Derive dates (ASTDT, AEDT, ...) ----
  # Derive analysis start date and flag
  derive_vars_dt(
    dtc = MHSTDTC,
    new_vars_prefix = "AST",
    date_imputation = "first"
  ) %>%
  # Derive analysis end date and flag
  derive_vars_dt(
    dtc = MHENDTC,
    new_vars_prefix = "AEN",
    date_imputation = "last",
    max_dates = vars(DTHDT, EOSDT)
  ) %>%
  # Derive analysis start relative day and analysis end relative day
  derive_vars_dy(
    reference_date = TRTSDT,
    source_vars = vars(ASTDT, AENDT)
  ) %>%
  # Derive analysis date of medical history collection - ADT (company specific variable derivation)
  derive_vars_dt(
    dtc = MHDTC,
    new_vars_prefix = "A"
  ) %>%
  # Derive analysis relative day - ADY (company specific variable derivation)
  derive_vars_dy(
    reference_date = TRTSDT,
    source_vars = vars(ADT)
  ) %>%
  ## Derive query variables ----
  derive_vars_query(queries_mh) %>%
  # Assign the AHIST (company specific variable derivation)
  mutate(AHIST = case_when(
    MHENRF == "BEFORE" ~ "Past",
    MHENRF %in% c("DURING", "AFTER") ~ "Current",
    MHSTAT == "Not Done" ~ "Not Assessed"
  )) %>%
  ## Derive occurrence flags ----
  derive_var_extreme_flag(
    by_vars = vars(USUBJID),
    order = vars(ASTDT, MHSEQ),
    new_var = AOCCFL,
    mode = "first"
  ) %>%
  derive_var_extreme_flag(
    by_vars = vars(USUBJID, MHBODSYS),
    order = vars(USUBJID, MHBODSYS, MHCAT, MHDECOD, MHTERM, ASTDT, MHSEQ),
    new_var = AOCCSFL,
    mode = "first"
  ) %>%
  derive_var_extreme_flag(
    by_vars = vars(USUBJID, MHDECOD),
    order = vars(USUBJID, MHBODSYS, MHCAT, MHDECOD, MHTERM, ASTDT, MHSEQ),
    new_var = AOCCPFL,
    mode = "first"
  ) %>%
  # (company specific occurrence flag variables derivation)
  derive_var_extreme_flag(
    by_vars = vars(USUBJID),
    order = vars(USUBJID, AHIST, ASTDT, MHSEQ),
    new_var = AOCPFL,
    mode = "first"
  ) %>%
  derive_var_extreme_flag(
    by_vars = vars(USUBJID, MHBODSYS),
    order = vars(USUBJID, AHIST, MHBODSYS, MHCAT, ASTDT, MHSEQ),
    new_var = AOCPSFL,
    mode = "first"
  ) %>%
  derive_var_extreme_flag(
    by_vars = vars(USUBJID, MHDECOD),
    order = vars(USUBJID, AHIST, MHBODSYS, MHCAT, MHDECOD, MHTERM, ASTDT, MHSEQ),
    new_var = AOCPPFL,
    mode = "first"
  ) %>%
  ## Derive analysis flag (company specific variable derivation) ----
  mutate(ANL01FL = ifelse(MHOCCUR != "N", "Y", NA_character_)) %>%
  ## Assign TRTA, TRTP (company specific variables derivation) ----
  # See also the "Visit and Period Variables" vignette
  # (https://pharmaverse.github.io/admiral/articles/visits_periods.html#treatment_bds)
  mutate(
    TRTP = TRT01P,
    TRTA = TRT01A
  ) %>%
  ## Assign APHASE and APHASEN Variable (company specific variable derivation) ----
  # See also the "Visit and Period Variables" vignette
  # (https://pharmaverse.github.io/admiral/articles/visits_periods.html#periods_bds)
  mutate(
    APHASE = case_when(
      ADT < TRTSDT ~ "Screening",
      ADT > TRTEDT ~ "Post-Treatment",
      ADT <= TRTSDT & ADT >= TRTEDT ~ "On-Treatment"
    ),
    APHASEN = case_when(
      ADT < TRTSDT ~ 1,
      ADT > TRTEDT ~ 2,
      ADT <= TRTSDT & ADT >= TRTEDT ~ 3
    )
  )

# Derive MHTERMN (company specific variable derivation)
admh <- restrict_derivation(
  admh,
  derivation = derive_vars_merged,
  args = params(
    dataset_add = mhtermn_lookup,
    by_vars = vars(MHTERM),
    new_vars = vars(MHTERMN)
  ),
  filter = (MHPRESP == "Y")
)

## Add all ADSL variables ----
admh <- admh %>%
  derive_vars_merged(
    dataset_add = select(adsl, !!!negate_vars(adsl_vars)),
    by_vars = vars(STUDYID, USUBJID)
  )


# Final Steps, Select final variables and Add labels
# This process will be based on your metadata, no example given for this reason
# ...

# Save output ----

dir <- tempdir() # Change to whichever directory you want to save the dataset in
saveRDS(admh, file = file.path(dir, "admh.rds"), compress = "bzip2")
