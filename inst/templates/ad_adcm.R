# Name: ADCM
#
# Label: Concomitant Medications Analysis Dataset
#
# Input: cm, adsl, suppcm, suppdm, ex
library(admiral)
library(dplyr)
library(lubridate)

# ---- Load source datasets ----

# Use e.g. haven::read_sas to read in .sas7bdat, or other suitable functions
# as needed and assign to the variables below.
# For illustration purposes read in admiral test data

data("cm")
data("adsl")
data("ex")
data("adsl")

# ---- Derivations ----
adcm <- cm %>%

  # Join supplementary qualifier variables
  # derive_vars_suppqual(suppcm)

  # Join adsl to cm
  left_join(adsl, by = c("STUDYID", "USUBJID")) %>%

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
  mutate(
    ASTDT = date(ASTDTM),
    AENDT = date(AENDTM)
  ) %>%

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
  ) %>%

  # Derive Time Relative to Reference
  derive_var_atirel(
    flag_var = ASTTMF,
    new_var = ATIREL
  ) %>%

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

  # Derive Aphase and Aphasen Variable
  # Other timing variable can be derived similarly.
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
  ) %>%

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

# ---- Save output ----

saveRDS(adcm, file = "./ADCM.rds", compress = TRUE)
