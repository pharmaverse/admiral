# Name: ADAE
#
# Label: Adverse Event Analysis Dataset
#
# Input: ae, adsl, ex_single
library(admiral)
library(pharmaversesdtm) # Contains example datasets from the CDISC pilot project
library(dplyr)
library(lubridate)

# Load source datasets ----

# Use e.g. haven::read_sas to read in .sas7bdat, or other suitable functions
# as needed and assign to the variables below.
# For illustration purposes read in admiral test data

data("ae")
data("admiral_adsl")
data("ex_single")
data("suppae")

adsl <- admiral_adsl

# When SAS datasets are imported into R using haven::read_sas(), missing
# character values from SAS appear as "" characters in R, instead of appearing
# as NA values. Further details can be obtained via the following link:
# https://pharmaverse.github.io/admiral/articles/admiral.html#handling-of-missing-values # nolint

ae <- convert_blanks_to_na(ae)
ex <- convert_blanks_to_na(ex_single)


# Derivations ----

# Get list of ADSL vars required for derivations
adsl_vars <- exprs(TRTSDT, TRTEDT, DTHDT, EOSDT)

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
  ## Derive analysis start relative day and  analysis end relative day ----
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

ex_ext <- derive_vars_dtm(
  ex,
  dtc = EXSTDTC,
  new_vars_prefix = "EXST",
  flag_imputation = "none"
)

adae <- adae %>%
  ## Derive last dose date/time ----
  derive_vars_joined(
    dataset_add = ex_ext,
    by_vars = exprs(STUDYID, USUBJID),
    new_vars = exprs(LDOSEDTM = EXSTDTM),
    join_vars = exprs(EXSTDTM),
    join_type = "all",
    order = exprs(EXSTDTM),
    filter_add = (EXDOSE > 0 | (EXDOSE == 0 & grepl("PLACEBO", EXTRT))) & !is.na(EXSTDTM),
    filter_join = EXSTDTM <= ASTDTM,
    mode = "last"
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
