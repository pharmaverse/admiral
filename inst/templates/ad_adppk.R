# Name: ADPPK
#
# Label: Population PK Analysis Dataset
#
# Description: Based on simulated data, create ADPPK analysis dataset
#
# Input: pc, ex, vs, adsl
library(admiral)
library(dplyr)
library(lubridate)
library(stringr)

library(admiral.test) # Contains example datasets from the CDISC pilot project or simulated

# ---- Load source datasets ----

# Use e.g. haven::read_sas to read in .sas7bdat, or other suitable functions
# as needed and assign to the variables below.
# For illustration purposes read in admiral test data


# Load PC, EX, VS and ADSL
data("admiral_pc")
data("admiral_ex")
data("admiral_vs")

data("admiral_adsl")

adsl <- admiral_adsl

# When SAS datasets are imported into R using haven::read_sas(), missing
# character values from SAS appear as "" characters in R, instead of appearing
# as NA values. Further details can be obtained via the following link:
# https://pharmaverse.github.io/admiral/cran-release/articles/admiral.html#handling-of-missing-values # nolint

# Load EX

ex <- convert_blanks_to_na(admiral_ex)

# Load PC

pc <- convert_blanks_to_na(admiral_pc)

# Load VS for baseline height and weight

vs <- convert_blanks_to_na(admiral_vs)

# ---- Lookup tables ----
param_lookup <- tibble::tribble(
  ~PCTESTCD, ~PARAMCD, ~PARAM, ~PARAMN,
  "XAN", "XAN", "Pharmacokinetic concentration of Xanomeline", 1,
  "DOSE", "DOSE", "Xanomeline Patch Dose", 2,
)

# ---- Derivations ----

# Get list of ADSL vars required for derivations
adsl_vars <- exprs(TRTSDT, TRTSDTM, TRT01P, TRT01A)

adppk <- pc %>%
  # Join ADSL with PC (need TRTSDT for ADY derivation)
  derive_vars_merged(
    dataset_add = adsl,
    new_vars = adsl_vars,
    by_vars = exprs(STUDYID, USUBJID)
  ) %>%
  # Derive analysis date/time
  # Impute missing time to 00:00:00
  derive_vars_dtm(
    new_vars_prefix = "A",
    dtc = PCDTC,
    time_imputation = "00:00:00"
  ) %>%
  # Derive dates and times from date/times
  derive_vars_dtm_to_dt(exprs(ADTM)) %>%
  derive_vars_dtm_to_tm(exprs(ADTM)) %>%
  derive_vars_dy(reference_date = TRTSDT, source_vars = exprs(ADT)) %>%
  # Derive event ID and nominal relative time from first dose (NFRLT)
  mutate(
    EVID = 0,
    DRUG = PCTEST,
    NFRLT = if_else(PCTPTNUM < 0, 0, PCTPTNUM), .after = USUBJID
  )


# ---- Get dosing information ----

ex <- ex %>%
  derive_vars_merged(
    dataset_add = adsl,
    new_vars = adsl_vars,
    by_vars = exprs(STUDYID, USUBJID)
  ) %>%
  # Keep records with nonzero dose
  filter(EXDOSE > 0) %>%
  # Add time and set missing end date to start date
  # Impute missing time to 00:00:00
  # Note all times are missing for dosing records in this example data
  # Derive Analysis Start and End Dates
  derive_vars_dtm(
    new_vars_prefix = "AST",
    dtc = EXSTDTC,
    time_imputation = "00:00:00"
  ) %>%
  derive_vars_dtm(
    new_vars_prefix = "AEN",
    dtc = EXENDTC,
    time_imputation = "00:00:00"
  ) %>%
  # Derive event ID and nominal relative time from first dose (NFRLT)
  mutate(
    EVID = 1,
    NFRLT = 24 * (VISITDY - 1), .after = USUBJID
  ) %>%
  # Set missing end dates to start date
  mutate(AENDTM = case_when(
    is.na(AENDTM) ~ ASTDTM,
    TRUE ~ AENDTM
  )) %>%
  # Derive dates from date/times
  derive_vars_dtm_to_dt(exprs(ASTDTM)) %>%
  derive_vars_dtm_to_dt(exprs(AENDTM))


# ---- Expand dosing records between start and end dates ----
# Updated function includes nominal_time parameter

ex_exp <- ex %>%
  create_single_dose_dataset(
    dose_freq = EXDOSFRQ,
    start_date = ASTDT,
    start_datetime = ASTDTM,
    end_date = AENDT,
    end_datetime = AENDTM,
    nominal_time = NFRLT,
    lookup_table = dose_freq_lookup,
    lookup_column = CDISC_VALUE,
    keep_source_vars = exprs(
      STUDYID, USUBJID, EVID, EXDOSFRQ, EXDOSFRM,
      NFRLT, EXDOSE, EXDOSU, EXTRT, ASTDT, ASTDTM, AENDT, AENDTM,
      VISIT, VISITNUM, VISITDY,
      TRT01A, TRT01P, DOMAIN, EXSEQ, !!!adsl_vars
    )
  ) %>%
  # Derive AVISIT based on nominal relative time
  # Derive AVISITN to nominal time in whole days using integer division
  # Define AVISIT based on nominal day
  mutate(
    AVISITN = NFRLT %/% 24 + 1,
    AVISIT = paste("Day", AVISITN),
    ADTM = ASTDTM,
    DRUG = EXTRT
  ) %>%
  # Derive dates and times from datetimes
  derive_vars_dtm_to_dt(exprs(ADTM)) %>%
  derive_vars_dtm_to_tm(exprs(ADTM)) %>%
  derive_vars_dtm_to_tm(exprs(ASTDTM)) %>%
  derive_vars_dtm_to_tm(exprs(AENDTM)) %>%
  derive_vars_dy(reference_date = TRTSDT, source_vars = exprs(ADT))


# ---- Find first dose per treatment per subject ----
# ---- Join with ADPPK data and keep only subjects with dosing ----

adppk <- adppk %>%
  derive_vars_merged(
    dataset_add = ex_exp,
    filter_add = (EXDOSE > 0 & !is.na(ADTM)),
    new_vars = exprs(FANLDTM = ADTM),
    order = exprs(ADTM, EXSEQ),
    mode = "first",
    by_vars = exprs(STUDYID, USUBJID, DRUG)
  ) %>%
  filter(!is.na(FANLDTM)) %>%
  # Derive AVISIT based on nominal relative time
  # Derive AVISITN to nominal time in whole days using integer division
  # Define AVISIT based on nominal day
  mutate(
    AVISITN = NFRLT %/% 24 + 1,
    AVISIT = paste("Day", AVISITN),
  )


# ---- Find previous dose  ----
# Use derive_vars_joined for consistency with other variables
# This is equivalent to derive_vars_last_dose in this case

adppk <- adppk %>%
  derive_vars_joined(
    dataset_add = ex_exp,
    by_vars = exprs(USUBJID),
    order = exprs(ADTM),
    new_vars = exprs(
      ADTM_prev = ADTM, EXDOSE_prev = EXDOSE, AVISIT_prev = AVISIT,
      AENDTM_prev = AENDTM
    ),
    join_vars = exprs(ADTM),
    filter_add = NULL,
    filter_join = ADTM > ADTM.join,
    mode = "last",
    check_type = "none"
  )

# ---- Find next dose  ----

adppk <- adppk %>%
  derive_vars_joined(
    dataset_add = ex_exp,
    by_vars = exprs(USUBJID),
    order = exprs(ADTM),
    new_vars = exprs(
      ADTM_next = ADTM, EXDOSE_next = EXDOSE, AVISIT_next = AVISIT,
      AENDTM_next = AENDTM
    ),
    join_vars = exprs(ADTM),
    filter_add = NULL,
    filter_join = ADTM <= ADTM.join,
    mode = "first",
    check_type = "none"
  )


# ---- Find previous nominal time ----

adppk <- adppk %>%
  derive_vars_joined(
    dataset_add = ex_exp,
    by_vars = exprs(USUBJID),
    order = exprs(NFRLT),
    new_vars = exprs(NFRLT_prev = NFRLT),
    join_vars = exprs(NFRLT),
    filter_add = NULL,
    filter_join = NFRLT > NFRLT.join,
    mode = "last",
    check_type = "none"
  )

# ---- Combine ADPPK and EX data ----
# Derive Relative Time Variables

adppk <- bind_rows(adppk, ex_exp) %>%
  group_by(USUBJID, DRUG) %>%
  mutate(
    FANLDTM = min(FANLDTM, na.rm = TRUE),
    min_NFRLT = min(NFRLT_prev, na.rm = TRUE),
    maxdate = max(ADT[EVID == 0], na.rm = TRUE), .after = USUBJID
  ) %>%
  arrange(USUBJID, ADTM) %>%
  ungroup() %>%
  filter(ADT <= maxdate) %>%
  # Derive Actual Relative Time from First Dose (AFRLT)
  derive_vars_duration(
    new_var = AFRLT,
    start_date = FANLDTM,
    end_date = ADTM,
    out_unit = "hours",
    floor_in = FALSE,
    add_one = FALSE
  ) %>%
  # Derive Actual Relative Time from Reference Dose (ARRLT)
  derive_vars_duration(
    new_var = ARRLT,
    start_date = ADTM_prev,
    end_date = ADTM,
    out_unit = "hours",
    floor_in = FALSE,
    add_one = FALSE
  ) %>%
  # Derive Actual Relative Time from Next Dose
  derive_vars_duration(
    new_var = AXRLT,
    start_date = ADTM_next,
    end_date = ADTM,
    out_unit = "hours",
    floor_in = FALSE,
    add_one = FALSE
  ) %>%
  mutate(
    ARRLT = case_when(
      EVID == 1 ~ 0,
      is.na(ARRLT) ~ AXRLT,
      TRUE ~ ARRLT
    ),
    # Derive Reference Dose Date
    PCRFTDTM = case_when(
      EVID == 1 ~ ADTM,
      is.na(ADTM_prev) ~ ADTM_next,
      TRUE ~ ADTM_prev
    )
  )

# Derive Nominal Relative Time from Reference Dose (NRRLT)

adppk <- adppk %>%
  mutate(
    NRRLT = case_when(
      EVID == 1 ~ 0,
      is.na(NFRLT_prev) ~ NFRLT - min_NFRLT,
      TRUE ~ NFRLT - NFRLT_prev
    )
  )

# ---- Derive Analysis Variables ----
# Derive ATPTN, ATPT, ATPTREF, ABLFL and BASETYPE
# Derive planned dose DOSEP, actual dose DOSEA and units
# Derive PARAMCD and relative time units
# Derive AVAL, AVALU and AVALCAT1

adppk <- adppk %>%
  mutate(
    STUDYIDN = as.numeric(word(USUBJID, 1, sep = fixed("-"))),
    USUBJIDN = as.numeric(word(USUBJID, 3, sep = fixed("-"))),
    ATPTN = case_when(
      EVID == 1 ~ 0,
      TRUE ~ PCTPTNUM
    ),
    ATPT = case_when(
      EVID == 1 ~ "Dose",
      TRUE ~ PCTPT
    ),
    ATPTREF = case_when(
      EVID == 1 ~ AVISIT,
      is.na(AVISIT_prev) ~ AVISIT_next,
      TRUE ~ AVISIT_prev
    ),
    # Derive baseline flag for pre-dose records
    ABLFL = case_when(
      ATPT == "Pre-dose" ~ "Y",
      TRUE ~ NA_character_
    ),
    # Derive Actual Dose
    DOSEA = case_when(
      EVID == 1 ~ EXDOSE,
      is.na(EXDOSE_prev) ~ EXDOSE_next,
      TRUE ~ EXDOSE_prev
    ),
    # Derive Planned Dose
    DOSEP = case_when(
      TRT01P == "Xanomeline High Dose" ~ 81,
      TRT01P == "Xanomeline Low Dose" ~ 54
    ),
    # Derive ARMN
    ARMN = case_when(
      TRT01P == "Xanomeline High Dose" ~ 2,
      TRT01P == "Xanomeline Low Dose" ~ 1
    ),
    # Derive ACTARMN
    ACTARMN = case_when(
      TRT01A == "Xanomeline High Dose" ~ 2,
      TRT01A == "Xanomeline Low Dose" ~ 1
    ),
    DOSEU = "mg",
  ) %>%
  # Derive relative time units
  mutate(
    # Derive PARAMCD
    PARAMCD = coalesce(PCTESTCD, "DOSE"),
    ALLOQ = PCLLOQ,
    # Derive CMT
    CMT = case_when(
      EVID == 1 ~ 1,
      TRUE ~ 2
    ),
    # Derive AVAL
    AVAL = case_when(
      EVID == 1 ~ EXDOSE,
      PCSTRESC == "<BLQ" ~ 0,
      TRUE ~ PCSTRESN
    ),
    AMT = case_when(
      EVID == 1 ~ EXDOSE,
      TRUE ~ NA_real_
    ),
    DV = AVAL,
    DVL = case_when(
      DV != 0 ~ log(DV),
      TRUE ~ NA_real_
    ),
    MDV = case_when(
      EVID == 1 ~ 1,
      is.na(DV) ~ 1,
      TRUE ~ NA
    ),
    AVALU = case_when(
      EVID == 1 ~ EXDOSU,
      TRUE ~ PCSTRESU
    )
  )


# ---- Add ASEQ ----

adppk <- adppk %>%
  # Calculate ASEQ
  derive_var_obs_number(
    new_var = ASEQ,
    by_vars = exprs(STUDYID, USUBJID),
    order = exprs(AFRLT, EVID),
    check_type = "error"
  ) %>%
  # Derive PARAM and PARAMN
  derive_vars_merged(dataset_add = select(param_lookup, -PCTESTCD), by_vars = exprs(PARAMCD)) %>%
  # Remove temporary variables
  select(
    -DOMAIN, -starts_with("min"), -starts_with("max"), -starts_with("EX"),
    -starts_with("PC"), -ends_with("next"), -starts_with("VIS"), -AXRLT,
    -ends_with("prev"), -ends_with("DTM"), -ends_with("DT"), -ends_with("TM"),
    -ends_with("TMF"), -starts_with("TRT"), -starts_with("ATPT"), -DRUG
  )

#---- Derive Covariates ----

# ---- Add all Covariate variables ----
# Include numeric values for STUDYIDN, USUBJIDN, SEXN, RACEN etc.


# Final Steps, Select final variables and Add labels
# This process will be based on your metadata, no example given for this reason
# ...
# ---- Save output ----

dir <- tempdir() # Change to whichever directory you want to save the dataset in
saveRDS(adppk, file = file.path(dir, "adppk.rds"), compress = "bzip2")
