# Name: ADPPK
#
# Label: Population PK Analysis Data
#
# Description: Based on simulated data, create ADPPK analysis dataset
#
# Input: pc, ex, vs, lb, adsl
library(admiral)
library(dplyr)
library(lubridate)
library(stringr)

library(admiral.test) # Contains example datasets from the CDISC pilot project or simulated

# ---- Load source datasets ----

# Use e.g. haven::read_sas to read in .sas7bdat, or other suitable functions
# as needed and assign to the variables below.
# For illustration purposes read in admiral test data

# Load PC, EX, VS, LB and ADSL
data("admiral_pc")
data("admiral_ex")
data("admiral_vs")
data("admiral_lb")

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

# Load LB for baseline lab values

lb <- convert_blanks_to_na(admiral_lb)

# ---- Lookup tables ----
param_lookup <- tibble::tribble(
  ~PCTESTCD, ~PARAMCD, ~PARAM, ~PARAMN,
  "XAN", "XAN", "Pharmacokinetic concentration of Xanomeline", 1,
  "DOSE", "DOSE", "Xanomeline Patch Dose", 2,
)

# ---- Derivations ----

# Get list of ADSL vars required for derivations
adsl_vars <- exprs(TRTSDT, TRTSDTM, TRT01P, TRT01A)

pc_dates <- pc %>%
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
  # Derive event ID and nominal relative time from first dose (NFRLT)
  mutate(
    EVID = 0,
    DRUG = PCTEST,
    NFRLT = if_else(PCTPTNUM < 0, 0, PCTPTNUM), .after = USUBJID
  )

# ---- Get dosing information ----

ex_dates <- ex %>%
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

ex_exp <- ex_dates %>%
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
  derive_vars_dtm_to_tm(exprs(AENDTM))


# ---- Find first dose per treatment per subject ----
# ---- Join with ADPPK data and keep only subjects with dosing ----

adppk_first_dose <- pc_dates %>%
  derive_vars_merged(
    dataset_add = ex_exp,
    filter_add = (!is.na(ADTM)),
    new_vars = exprs(FANLDTM = ADTM, EXDOSE_first = EXDOSE),
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

adppk_prev <- adppk_first_dose %>%
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

# ---- Find previous nominal dose ----

adppk_nom_prev <- adppk_prev %>%
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

adppk_aprlt <- bind_rows(adppk_nom_prev, ex_exp) %>%
  group_by(USUBJID, DRUG) %>%
  mutate(
    FANLDTM = min(FANLDTM, na.rm = TRUE),
    min_NFRLT = min(NFRLT, na.rm = TRUE),
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
  # Derive Actual Relative Time from Reference Dose (APRLT)
  derive_vars_duration(
    new_var = APRLT,
    start_date = ADTM_prev,
    end_date = ADTM,
    out_unit = "hours",
    floor_in = FALSE,
    add_one = FALSE
  ) %>%
  # Derive APRLT
  mutate(
    APRLT = case_when(
      EVID == 1 ~ 0,
      is.na(APRLT) ~ AFRLT,
      TRUE ~ APRLT
    ),
    NPRLT = case_when(
      EVID == 1 ~ 0,
      is.na(NFRLT_prev) ~ NFRLT - min_NFRLT,
      TRUE ~ NFRLT - NFRLT_prev
    )
  )

# ---- Derive Analysis Variables ----
# Derive actual dose DOSEA and planned dose DOSEP,
# Derive AVAL and DV

adppk_aval <- adppk_aprlt %>%
  mutate(
    # Derive Actual Dose
    DOSEA = case_when(
      EVID == 1 ~ EXDOSE,
      is.na(EXDOSE_prev) ~ EXDOSE_first,
      TRUE ~ EXDOSE_prev
    ),
    # Derive Planned Dose
    DOSEP = case_when(
      TRT01P == "Xanomeline High Dose" ~ 81,
      TRT01P == "Xanomeline Low Dose" ~ 54,
      TRT01P == "Placebo" ~ 0
    ),
    # Derive PARAMCD
    PARAMCD = case_when(
      EVID == 1 ~ "DOSE",
      TRUE ~ PCTESTCD
    ),
    ALLOQ = PCLLOQ,
    # Derive CMT
    CMT = case_when(
      EVID == 1 ~ 1,
      TRUE ~ 2
    ),
    # Derive BLQFL/BLQFN
    BLQFL = case_when(
      PCSTRESC == "<BLQ" ~ "Y",
      TRUE ~ "N"
    ),
    BLQFN = case_when(
      PCSTRESC == "<BLQ" ~ 1,
      TRUE ~ 0
    ),
    AMT = case_when(
      EVID == 1 ~ EXDOSE,
      TRUE ~ NA_real_
    ),
    # Derive DV and AVAL
    DV = PCSTRESN,
    AVAL = DV,
    DVL = case_when(
      DV != 0 ~ log(DV),
      TRUE ~ NA_real_
    ),
    # Derive MDV
    MDV = case_when(
      EVID == 1 ~ 1,
      is.na(DV) ~ 1,
      TRUE ~ 0
    ),
    AVALU = case_when(
      EVID == 1 ~ NA_character_,
      TRUE ~ PCSTRESU
    ),
    UDTC = format_ISO8601(ADTM),
    II = if_else(EVID == 1, 1, 0),
    SS = if_else(EVID == 1, 1, 0)
  )

# ---- Add ASEQ ----

adppk_aseq <- adppk_aval %>%
  # Calculate ASEQ
  derive_var_obs_number(
    new_var = ASEQ,
    by_vars = exprs(STUDYID, USUBJID),
    order = exprs(AFRLT, EVID),
    check_type = "error"
  ) %>%
  # Derive PARAM and PARAMN
  derive_vars_merged(dataset_add = select(param_lookup, -PCTESTCD), by_vars = exprs(PARAMCD)) %>%
  mutate(
    PROJID = DRUG,
    PROJIDN = 1
  ) %>%
  # Remove temporary variables
  select(
    -DOMAIN, -starts_with("min"), -starts_with("max"), -starts_with("EX"),
    -starts_with("PC"), -ends_with("first"), -ends_with("prev"),
    -ends_with("DTM"), -ends_with("DT"), -ends_with("TM"), -starts_with("VISIT"),
    -starts_with("AVISIT"), -starts_with("PARAM"),
    -ends_with("TMF"), -starts_with("TRT"), -starts_with("ATPT"), -DRUG
  )

#---- Derive Covariates ----
# Include numeric values for STUDYIDN, USUBJIDN, SEXN, RACEN etc.

covar <- adsl %>%
  mutate(
    STUDYIDN = as.numeric(word(USUBJID, 1, sep = fixed("-"))),
    SITEIDN = as.numeric(word(USUBJID, 2, sep = fixed("-"))),
    USUBJIDN = as.numeric(word(USUBJID, 3, sep = fixed("-"))),
    SUBJIDN = as.numeric(SUBJID),
    SEXN = case_when(
      SEX == "M" ~ 1,
      SEX == "F" ~ 2,
      TRUE ~ 3
    ),
    RACEN = case_when(
      RACE == "AMERICAN INDIAN OR ALASKA NATIVE" ~ 1,
      RACE == "ASIAN" ~ 2,
      RACE == "BLACK OR AFRICAN AMERICAN" ~ 3,
      RACE == "NATIVE HAWAIIAN OR OTHER PACIFIC ISLANDER" ~ 4,
      RACE == "WHITE" ~ 5,
      TRUE ~ 6
    ),
    ETHNICN = case_when(
      ETHNIC == "HISPANIC OR LATINO" ~ 1,
      ETHNIC == "NOT HISPANIC OR LATINO" ~ 2,
      TRUE ~ 3
    ),
    ARMN = case_when(
      ARM == "Placebo" ~ 0,
      ARM == "Xanomeline Low Dose" ~ 1,
      ARM == "Xanomeline High Dose" ~ 2,
      TRUE ~ 3
    ),
    ACTARMN = case_when(
      ACTARM == "Placebo" ~ 0,
      ACTARM == "Xanomeline Low Dose" ~ 1,
      ACTARM == "Xanomeline High Dose" ~ 2,
      TRUE ~ 3
    ),
    COHORT = ARMN,
    COHORTC = ARM,
    ROUTE = unique(ex$EXROUTE),
    ROUTEN = case_when(
      ROUTE == "TRANSDERMAL" ~ 3,
      TRUE ~ NA_real_
    ),
    FORM = unique(ex$EXDOSFRM),
    FORMN = case_when(
      FORM == "PATCH" ~ 3,
      TRUE ~ 4
    ),
    COUNTRYN = case_when(
      COUNTRY == "USA" ~ 1,
      COUNTRY == "CAN" ~ 2,
      COUNTRY == "GBR" ~ 3
    ),
    REGION1N = COUNTRYN,
  ) %>%
  select(
    STUDYID, STUDYIDN, SITEID, SITEIDN, USUBJID, USUBJIDN,
    SUBJID, SUBJIDN, AGE, SEX, SEXN, COHORT, COHORTC, ROUTE, ROUTEN,
    RACE, RACEN, ETHNIC, ETHNICN, FORM, FORMN, COUNTRY, COUNTRYN,
    REGION1, REGION1N
  )

#---- Derive additional baselines from VS and LB ----

labsbl <- lb %>%
  filter(LBBLFL == "Y" & LBTESTCD %in% c("CREAT", "ALT", "AST", "BILI")) %>%
  mutate(LBTESTCDB = paste0(LBTESTCD, "BL")) %>%
  select(STUDYID, USUBJID, LBTESTCDB, LBSTRESN)

covar_vslb <- covar %>%
  derive_vars_merged(
    dataset_add = vs,
    filter_add = VSTESTCD == "HEIGHT",
    by_vars = exprs(STUDYID, USUBJID),
    new_vars = exprs(HTBL = VSSTRESN)
  ) %>%
  derive_vars_merged(
    dataset_add = vs,
    filter_add = VSTESTCD == "WEIGHT" & VSBLFL == "Y",
    by_vars = exprs(STUDYID, USUBJID),
    new_vars = exprs(WTBL = VSSTRESN)
  ) %>%
  derive_vars_transposed(
    dataset_merge = labsbl,
    by_vars = exprs(STUDYID, USUBJID),
    key_var = LBTESTCDB,
    value_var = LBSTRESN
  ) %>%
  mutate(
    BMIBL = compute_bmi(height = HTBL, weight = WTBL),
    BSABL = compute_bsa(
      height = HTBL,
      weight = HTBL,
      method = "Mosteller"
    ),
    CRCLBL = compute_egfr(
      creat = CREATBL, creatu = "SI", age = AGE, wt = WTBL, sex = SEX,
      method = "CRCL"
    ),
    EGFRBL = compute_egfr(
      creat = CREATBL, creatu = "SI", age = AGE, wt = WTBL, sex = SEX,
      method = "CKD-EPI"
    )
  ) %>%
  rename(TBILBL = BILIBL)

# Combine covariates with APPPK data

adppk <- adppk_aseq %>%
  derive_vars_merged(
    dataset_add = covar_vslb,
    by_vars = exprs(STUDYID, USUBJID)
  ) %>%
  arrange(STUDYIDN, USUBJIDN, AFRLT, EVID) %>%
  mutate(RECSEQ = row_number())

# Final Steps, Select final variables and Add labels
# This process will be based on your metadata, no example given for this reason
# ...
# ---- Save output ----

dir <- tempdir() # Change to whichever directory you want to save the dataset in
saveRDS(adppk, file = file.path(dir, "adppk.rds"), compress = "bzip2")
