# Name: ADAE
#
# Label: Adverse Event Analysis Dataset
#
# Input: ae, adsl, ex
library(admiral)
library(pharmaversesdtm) # Contains example datasets from the CDISC pilot project
library(dplyr)
library(lubridate)

# Load source datasets ----

# Use e.g. haven::read_sas to read in .sas7bdat, or other suitable functions
# as needed and assign to the variables below.
# For illustration purposes read in admiral test data

ae <- pharmaversesdtm::ae
adsl <- admiral::admiral_adsl
ex <- pharmaversesdtm::ex

# When SAS datasets are imported into R using haven::read_sas(), missing
# character values from SAS appear as "" characters in R, instead of appearing
# as NA values. Further details can be obtained via the following link:
# https://pharmaverse.github.io/admiral/cran-release/articles/admiral.html#handling-of-missing-values # nolint
ae <- convert_blanks_to_na(ae)
ex <- convert_blanks_to_na(ex)

# Create single dose dataset for use in date of last dose derivations.
# If the exposure dataset contains multi-day dosing records (e.g., one record
# per treatment period rather than one record per dose), use
# create_single_dose_dataset() to expand them into one record per dose.
# Whether this step is necessary depends on how dosing data were collected.
ex_single <- ex %>%
  derive_vars_dtm(
    dtc = EXSTDTC,
    new_vars_prefix = "EXST",
    time_imputation = "first",
    flag_imputation = "none"
  ) %>%
  derive_vars_dtm(
    dtc = EXENDTC,
    new_vars_prefix = "EXEN",
    time_imputation = "last",
    flag_imputation = "none"
  ) %>%
  derive_vars_dtm_to_dt(exprs(EXSTDTM, EXENDTM)) %>%
  filter(!is.na(EXSTDT), !is.na(EXENDT)) %>%
  create_single_dose_dataset(
    dose_freq = EXDOSFRQ,
    start_date = EXSTDT,
    start_datetime = EXSTDTM,
    end_date = EXENDT,
    end_datetime = EXENDTM,
    keep_source_vars = exprs(
      STUDYID, USUBJID, EXTRT, EXDOSE, EXDOSU, EXDOSFRQ, EXSTDT, EXENDT, EXSTDTM, EXENDTM
    )
  )

# Derivations ----

# Get list of ADSL vars required for derivations
adsl_vars <- exprs(TRTSDT, TRTEDT, TRTEDTM, DTHDT, EOSDT)

adae <- ae %>%
  # join adsl to ae
  derive_vars_merged(
    dataset_add = adsl,
    new_vars = adsl_vars,
    by = exprs(STUDYID, USUBJID)
  ) %>%
  ## Derive analysis start time ----
  derive_vars_dtm(
    dtc = AESTDTC,
    new_vars_prefix = "AST",
    highest_imputation = "M",
    min_dates = exprs(TRTSDT)
  ) %>%
  ## Derive analysis end time ----
  derive_vars_dtm(
    dtc = AEENDTC,
    new_vars_prefix = "AEN",
    highest_imputation = "M",
    date_imputation = "last",
    time_imputation = "last",
    max_dates = exprs(DTHDT, EOSDT)
  ) %>%
  ## Derive analysis end/start date ----
  derive_vars_dtm_to_dt(exprs(ASTDTM, AENDTM)) %>%
  ## Derive analysis start relative day and analysis end relative day ----
  derive_vars_dy(
    reference_date = TRTSDT,
    source_vars = exprs(ASTDT, AENDT)
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

adae <- adae %>%
  ## Derive last dose date/time ----
  derive_vars_joined(
    dataset_add = ex_single,
    by_vars = exprs(STUDYID, USUBJID),
    new_vars = exprs(LDOSEDTM = EXSTDTM),
    join_vars = exprs(EXSTDTM),
    join_type = "all",
    order = exprs(EXSTDTM),
    filter_add = (EXDOSE > 0 | (EXDOSE == 0 & grepl("PLACEBO", EXTRT))) & !is.na(EXSTDTM),
    filter_join = EXSTDTM <= ASTDTM,
    mode = "last"
  ) %>%
  ## Derive treatment dose and unit ----
  # ex_single has been expanded to one record per day so EXSTDTM and EXENDTM
  # represent day-level treatment windows (00:00:00 to 23:59:59).
  # It is assumed that dosing intervals do not overlap; if they do,
  # derive_vars_joined() will throw an error as handling this is study-specific.
  # Drug clearance duration should be considered when matching exposure records
  # with adverse events. EXSTDTC and EXENDTC represent the administration period
  # only, not the time the drug remains in the body. The clearance buffer is
  # applied to TRTEDTM (last exposure end datetime from ADSL) rather than to
  # every exposure record, to avoid duplicate matches across daily rows.
  # Replace days(1) with the study-specific drug clearance period. If no
  # clearance buffer is required, simplify filter_join to:
  # EXSTDTM <= ASTDTM & ASTDTM <= EXENDTM
  derive_vars_joined(
    dataset_add = ex_single,
    by_vars = exprs(STUDYID, USUBJID),
    new_vars = exprs(DOSEON = EXDOSE, DOSEU = EXDOSU),
    join_vars = exprs(EXSTDTM, EXENDTM),
    join_type = "all",
    filter_add = (EXDOSE > 0 | (EXDOSE == 0 & grepl("PLACEBO", EXTRT))) & !is.na(EXSTDTM),
    filter_join = EXSTDTM <= ASTDTM & (
      ASTDTM <= EXENDTM |
        # Apply study-specific drug clearance window to last dose only.
        # Replace days(1) with the appropriate clearance period for your study.
        (EXENDTM == TRTEDTM & ASTDTM <= TRTEDTM + days(1))
    )
  ) %>%
  ## Derive severity / causality / ... ----
  mutate(
    ASEV = AESEV,
    AREL = AEREL
  ) %>%
  ## Derive treatment emergent flag ----
  derive_var_trtemfl(
    trt_start_date = TRTSDT,
    trt_end_date = TRTEDT,
    end_window = 30
  ) %>%
  ## Derive occurrence flags: first occurrence of most severe AE ----
  # create numeric value ASEVN for severity
  mutate(
    ASEVN = as.integer(factor(ASEV, levels = c("MILD", "MODERATE", "SEVERE", "DEATH THREATENING")))
  ) %>%
  restrict_derivation(
    derivation = derive_var_extreme_flag,
    args = params(
      by_vars = exprs(USUBJID),
      order = exprs(desc(ASEVN), ASTDTM, AESEQ),
      new_var = AOCCIFL,
      mode = "first"
    ),
    filter = TRTEMFL == "Y"
  )

# Join all ADSL with AE
adae <- adae %>%
  derive_vars_merged(
    dataset_add = select(adsl, !!!negate_vars(adsl_vars)),
    by_vars = exprs(STUDYID, USUBJID)
  )


# Save output ----

# Change to whichever directory you want to save the dataset in
dir <- tools::R_user_dir("admiral_templates_data", which = "cache")
if (!file.exists(dir)) {
  # Create the folder
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
}
save(adae, file = file.path(dir, "adae.rda"), compress = "bzip2")
