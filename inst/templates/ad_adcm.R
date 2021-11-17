# Name: ADCM
#
# Label: Concomitant Medications Analysis Dataset
#
# Input: cm, adsl
library(admiral)
library(admiral.test) # Contains example datasets from the CDISC pilot project
library(dplyr)
library(lubridate)

# ---- Load source datasets ----

# Use e.g. haven::read_sas to read in .sas7bdat, or other suitable functions
# as needed and assign to the variables below.
# For illustration purposes read in admiral test data

data("cm")
data("adsl")

cm <- convert_blanks_to_na(cm)

# ---- Derivations ----

# Get list of ADSL vars required for derivations
adsl_vars <- vars(TRTSDT, TRTEDT, DTHDT, EOSDT, TRT01P, TRT01A)

# Derive flags
adcm <- cm %>%

  # Join supplementary qualifier variables
  # derive_vars_suppqual(suppcm) %>%

  # Join ADSL with CM (only ADSL vars required for derivations)
  left_join(
    select(adsl, STUDYID, USUBJID, !!!adsl_vars),
    by = c("STUDYID", "USUBJID")
  ) %>%

  # Derive analysis start time
  derive_vars_dtm(
    dtc = CMSTDTC,
    new_vars_prefix = "AST",
    date_imputation = "first",
    time_imputation = "first",
    min_dates = vars(TRTSDT)
  ) %>%

  # Derive analysis end time
  derive_vars_dtm(
    dtc = CMENDTC,
    new_vars_prefix = "AEN",
    date_imputation = "last",
    time_imputation = "last",
    max_dates = vars(DTHDT, EOSDT)
  ) %>%

  # Derive analysis end/start date
  derive_vars_dtm_to_dt(vars(ASTDTM, AENDTM)) %>%

  # Derive analysis start relative day
  derive_var_astdy(
    reference_date = TRTSDT,
    date = ASTDT
  ) %>%

  # Derive analysis end relative day
  derive_var_aendy(
    reference_date = TRTSDT,
    date = AENDT
  ) %>%

  # Derive analysis duration (value and unit)
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

# Derive flags
adcm <- adcm %>%

  # Derive On-Treatment flag
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
  derive_extreme_flag(
    new_var = AOCCPFL,
    by_vars = vars(USUBJID, CMDECOD),
    order = vars(ASTDTM, CMSEQ),
    filter = ANL01FL == "Y",
    mode = "first"
  )


# Derive Aphase and Aphasen Variable
# Other timing variable can be derived similarly.
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

  left_join(select(adsl, !!!admiral:::negate_vars(adsl_vars)),
            by = c("STUDYID", "USUBJID")
  )



# ---- Save output ----

save(adcm, file = "data/adcm.rda", compress = "bzip2")
