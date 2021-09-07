# Name: ADAE
#
# Label: Adverse Event Analysis Dataset
#
# Input: ae, adsl, suppae, suppdm, ex
library(admiral)
library(dplyr)
library(lubridate)

# ---- Load source datasets ----

# Use e.g. haven::read_sas to read in .sas7bdat, or other suitable functions
# as needed and assign to the variables below.
# For illustration purposes read in admiral test data

data("ae")
data("suppae")
data("adsl")
data("ex_single")
ex <- ex_single

# ---- Derivations ----

adae <- ae %>%
  # join supplementary qualifier variables
  derive_vars_suppqual(suppae) %>%

  # join adsl to ae
  left_join(adsl, by = c("STUDYID", "USUBJID")) %>%

  # derive analysis start time
  derive_vars_dtm(
    dtc = AESTDTC,
    new_vars_prefix = "AST",
    date_imputation = "first",
    time_imputation = "first",
    min_dates = list(TRTSDT)
  ) %>%

  # derive analysis end time
  derive_vars_dtm(
    dtc = AEENDTC,
    new_vars_prefix = "AEN",
    date_imputation = "last",
    time_imputation = "last",
    max_dates = list(DTHDT, EOSDT)
  ) %>%

  # derive analysis end/start date
  mutate(
    ASTDT = date(ASTDTM),
    AENDT = date(AENDTM)
  ) %>%

  # derive analysis start relative day
  derive_var_astdy(
    reference_date = TRTSDT,
    date = ASTDT
  ) %>%

  # derive analysis end relative day
  derive_var_aendy(
    reference_date = TRTSDT,
    date = AENDT
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
  ) %>%

  # derive last dose date/time
  derive_last_dose(
    ex,
    filter_ex = (EXDOSE > 0 | (EXDOSE == 0 & grepl("PLACEBO", EXTRT))) &
      nchar(EXENDTC) >= 10,
    dose_start = EXSTDTC,
    dose_end = EXENDTC,
    analysis_date = ASTDT,
    dataset_seq_var = AESEQ,
    new_var = LDOSEDTM,
    output_datetime = TRUE,
    check_dates_only = FALSE
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

  # derive occurrence flags
  derive_extreme_flag(
    by_vars = vars(USUBJID),
    order = vars(ASTDTM, AESEQ),
    new_var = AOCCIFL,
    filter = TRTEMFL == "Y",
    mode = "last"
  )

# ---- Save output ----

saveRDS(adae, file = "./ADAE.rds", compress = TRUE)
