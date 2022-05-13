# Name: ADAE
#
# Label: Adverse Event Analysis Dataset
#
# Input: ae, adsl, suppae, ex_single
library(admiral)
library(admiraltest) # Contains example datasets from the CDISC pilot project
library(dplyr)
library(lubridate)

# ---- Load source datasets ----

# Use e.g. haven::read_sas to read in .sas7bdat, or other suitable functions
# as needed and assign to the variables below.
# For illustration purposes read in admiral test data

data("admiral_ae")
data("admiral_suppae")
data("admiral_adsl")
data("ex_single")

adsl <- admiral_adsl
ae <- admiral_ae
suppae <- admiral_suppae

ae <- convert_blanks_to_na(ae)
suppae <- convert_blanks_to_na(suppae)
ex <- convert_blanks_to_na(ex_single)


# ---- Derivations ----

# Get list of ADSL vars required for derivations
adsl_vars <- vars(TRTSDT, TRTEDT, DTHDT, EOSDT)

adae <- ae %>%
  # join supplementary qualifier variables
  derive_vars_suppqual(suppae) %>%
  # join adsl to ae
  derive_vars_merged(
    dataset_add = adsl,
    new_vars = adsl_vars,
    by = vars(STUDYID, USUBJID)
  ) %>%
  # derive analysis start time
  derive_vars_dtm(
    dtc = AESTDTC,
    new_vars_prefix = "AST",
    date_imputation = "first",
    time_imputation = "first",
    min_dates = vars(TRTSDT)
  ) %>%
  # derive analysis end time
  derive_vars_dtm(
    dtc = AEENDTC,
    new_vars_prefix = "AEN",
    date_imputation = "last",
    time_imputation = "last",
    max_dates = vars(DTHDT, EOSDT)
  ) %>%
  # derive analysis end/start date
  derive_vars_dtm_to_dt(vars(ASTDTM, AENDTM)) %>%
  # derive analysis start relative day and  analysis end relative day
  derive_vars_dy(
    reference_date = TRTSDT,
    source_vars = vars(ASTDT, AENDT)
  ) %>%
  # derive analysis duration (value and unit)
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


adae <- adae %>%
  # derive last dose date/time
  derive_var_last_dose_date(
    ex,
    filter_ex = (EXDOSE > 0 | (EXDOSE == 0 & grepl("PLACEBO", EXTRT))) &
      nchar(EXENDTC) >= 10,
    dose_date = EXSTDTC,
    analysis_date = ASTDT,
    new_var = LDOSEDTM,
    single_dose_condition = (EXSTDTC == EXENDTC),
    output_datetime = TRUE
  ) %>%
  # derive severity / causality / ...
  mutate(
    ASEV = AESEV,
    AREL = AEREL
  ) %>%
  # derive treatment emergent flag
  mutate(
    TRTEMFL = ifelse(ASTDT >= TRTSDT & ASTDT <= TRTEDT + days(30), "Y", NA_character_)
  ) %>%
  # derive occurrence flags: first occurence of most severe AE
  # create numeric value ASEVN for severity
  mutate(
    ASEVN = as.integer(factor(ASEV, levels = c("MILD", "MODERATE", "SEVERE", "DEATH THREATENING")))
  ) %>%
  restrict_derivation(
    derivation = derive_var_extreme_flag,
    args = params(
      by_vars = vars(USUBJID),
      order = vars(ASTDTM, AESEQ),
      new_var = AOCCIFL,
      mode = "last"
    ),
    filter = TRTEMFL == "Y"
  )

# Join all ADSL with AE
adae <- adae %>%
  derive_vars_merged(
    dataset_add = select(adsl, !!!negate_vars(adsl_vars)),
    by_vars = vars(STUDYID, USUBJID)
  )


# ---- Save output ----

dir <- tempdir() # Change to whichever directory you want to save the dataset in
save(adae, file = file.path(dir, "adae.rda"), compress = "bzip2")
