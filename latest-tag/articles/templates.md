# Explore the admiral ADaM Templates

## Introduction

This page is a catalog of the ADaM templates available through
[admiral](https://pharmaverse.github.io/admiral/). These are lifted from
our GitHub repository. The intention is for this to be a reference for
users, should they wish to peruse the code of a specific template
without visiting the [admiral](https://pharmaverse.github.io/admiral/)
repository.

As a reminder, you can create a starter script from a template using
[`admiral::use_ad_template()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/use_ad_template.md),
e.g.,

``` r
use_ad_template(
  adam_name = "adsl",
  save_path = "./ad_adsl.R",
  package = "admiral"
)
```

## ADaM templates

ad_adab.R

```
# Name: ADAB
#
# Label: Anti-Drug Antibody Analysis Dataset
#
# Description: Experimental. Based on simulated data, create ADAB analysis dataset
#
# Input: is_ada, ex, adsl

library(admiral)
library(dplyr)
library(lubridate)
library(stringr)
library(pharmaversesdtm) # Contains example datasets from the CDISC pilot project or simulated

# ---- Load source datasets ----

# Use e.g. haven::read_sas to read in .sas7bdat, or other suitable functions
# as needed and assign to the variables below.
# For illustration purposes read in admiral test data

# Load IS, EX and ADSL from pharmaversesdtm and admiral

is <- pharmaversesdtm::is_ada
ex <- pharmaversesdtm::ex
adsl <- admiral::admiral_adsl

# When SAS datasets are imported into R using haven::read_sas(), missing
# character values from SAS appear as "" characters in R, instead of appearing
# as NA values. Further details can be obtained via the following link:
# https://pharmaverse.github.io/admiral/cran-release/articles/admiral.html#handling-of-missing-values # nolint

ex <- convert_blanks_to_na(ex)
is <- convert_blanks_to_na(is)
adsl <- convert_blanks_to_na(adsl)

# Define values for records with overall values
# Suggested are AVISIT=Overall, AVISITN=11111
overall_avisit <- "OVERALL"
overall_avisitn <- 11111

# Derivations ----

is_dates <- is %>%
  # Filter as needed (i.e. exclude ISSTAT has "NOT DONE")
  filter(!(ISSTAT %in% c("NOT DONE")) & toupper(ISBDAGNT) == "XANOMELINE") %>%
  # Initial preparation and core variables
  mutate(
    # ADATYPE and ADAPARM are assigned values to serve BY analyte processing

    # Uncomment if Setting ADATYPE based on SDTM V1.x ISTESTCD
    # ADATYPE = case_when(
    #   toupper(ISTESTCD) == "ADATEST1" ~ "ADA_BAB",
    #   toupper(ISTESTCD) == "NABTEST1" ~ "ADA_NAB",
    #   TRUE ~ NA_character_
    # ),

    # Uncomment if Setting ADAPARM based on SDTM V1.x ISTESTCD
    # ADAPARM = ISTESTCD,

    # Setting ADATYPE based on SDTM V2.x ISTESTCD (assumed to have ADA_BAB, ADA_NAB)
    # Remove or comment out if not >= SDTM V2.x
    ADATYPE = ISTESTCD,

    # Setting ADAPARM based on SDTM V2.x ISBDAGNT
    # Remove or comment out if not >= SDTM V2.x
    ADAPARM = ISBDAGNT,

    # When SDTM V1.x, Setting ISBDAGNT from ISTESTCD to work with template
    # ISBDAGNT = ISTESTCD,

    # Map the analyte test to corresponding DRUG in on EX.EXTRT
    # This is especially critical when multiple analytes and EX.EXTRT instances
    DRUG = case_when(
      toupper(ADAPARM) == "XANOMELINE" ~ "XANOMELINE",
      toupper(ADAPARM) == "OTHER_DRUG" ~ "OTHER_DRUG",
      TRUE ~ NA_character_
    ),
    # Set AVISIT and AVISITN based on VISIT and VISITNUM
    AVISIT = VISIT,
    AVISITN = VISITNUM
  ) %>%
  # Assign nominal time to NFRLT in DAYS
  # Special visits can be set to NA then can add more code to assign custom values
  #   (i.e. UNSCHEDULED to 99999, etc.)
  derive_var_nfrlt(
    new_var = NFRLT,
    new_var_unit = FRLTU,
    out_unit = "DAYS",
    tpt_var = ISTPT,
    visit_day = VISITDY,
    treatment_duration = 0,
    set_values_to_na = str_detect(toupper(VISIT), "UNSCHED") | str_detect(toupper(VISIT), "TREATMENT DISC")
  ) %>%
  mutate(
    NFRLT = case_when(
      str_detect(toupper(VISIT), "TREATMENT DISC") ~ 99997,
      str_detect(toupper(VISIT), "UNSCHED") ~ 99999,
      TRUE ~ NFRLT
    )
  ) %>%
  # Join ADSL with is (need TRTSDT for ADY derivation)
  derive_vars_merged(
    dataset_add = adsl,
    new_vars = exprs(TRTSDT),
    by_vars = exprs(STUDYID, USUBJID)
  ) %>%
  # Derive analysis date/time then compute ADY
  # Impute missing time to 00:00:00 or as desired.
  # Could replace this code with custom imputation code or function
  derive_vars_dtm(
    new_vars_prefix = "A",
    highest_imputation = "s",
    dtc = ISDTC,
    ignore_seconds_flag = FALSE,
    time_imputation = "00:00:00"
  ) %>%
  # Derive dates and times from date/times
  derive_vars_dtm_to_dt(exprs(ADTM)) %>%
  derive_vars_dtm_to_tm(exprs(ADTM)) %>%
  derive_vars_dy(reference_date = TRTSDT, source_vars = exprs(ADT))

# ---- Get dosing information ----

ex_dates <- ex %>%
  # Keep applicable desired records based on EXTRT and/or dose values (>=0, >0, etc.)
  filter(
    str_detect(toupper(EXTRT), "XANOMELINE") | str_detect(toupper(EXTRT), "PLACEBO"),
    EXDOSE >= 0
  ) %>%
  mutate(
    # DRUG is a merge variable to map and merge ADA data with EX.EXTRT
    # This will be used to merge first dose into IS working data.
    # PLACEBO example is for if/when ADA was also collected on Placebo subjects
    # or treatments are scrambled prior to database lock.
    DRUG = case_when(
      str_detect(toupper(EXTRT), "XANOMELINE") ~ "XANOMELINE",
      str_detect(toupper(EXTRT), "PLACEBO") ~ "XANOMELINE",
      str_detect(toupper(EXTRT), "OTHER_DRUG") ~ "OTHER_DRUG",
      TRUE ~ NA_character_
    )
  ) %>%
  # Assign nominal time to NFRLT in DAYS
  derive_var_nfrlt(
    new_var = NFRLT,
    new_var_unit = FRLTU,
    out_unit = "DAYS",
    visit_day = VISITDY,
    treatment_duration = 0
  ) %>%
  # Add analysis datetime variables and set missing end date to start date
  # Impute missing time to 00:00:00 or as desired.
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
  # Set missing end dates to start date or as desired
  mutate(
    AENDTM = if_else(is.na(AENDTM), ASTDTM, AENDTM)
  ) %>%
  # Derive dates from date/times
  derive_vars_dtm_to_dt(exprs(ASTDTM)) %>%
  derive_vars_dtm_to_dt(exprs(AENDTM))

# Note: This template computes only AFRLT, see ADPC template if need EX dose expansion example.

# Derive AFRLT in IS data
is_afrlt <- is_dates %>%
  derive_vars_merged(
    dataset_add = ex_dates,
    filter_add = (EXDOSE >= 0 & !is.na(ASTDTM)),
    new_vars = exprs(FANLDTM = ASTDTM, FANLTMF = ASTTMF),
    order = exprs(ASTDTM, EXSEQ),
    mode = "first",
    by_vars = exprs(STUDYID, USUBJID, DRUG)
  ) %>%
  derive_vars_dtm_to_dt(exprs(FANLDTM)) %>%
  derive_vars_dtm_to_tm(exprs(FANLDTM)) %>%
  derive_vars_duration(
    new_var = AFRLT,
    start_date = FANLDTM,
    end_date = ADTM,
    out_unit = "DAYS",
    floor_in = FALSE,
    add_one = FALSE
  )

# Compute or assign BASETYPE, APERIOD and APHASE ----------------------------------
# Add study specific code as applicable using ADEX or ADSL APxx / PHw variables

is_basetype <- is_afrlt %>%
  mutate(
    APERIOD = 1,
    APERIODC = "Period 01",
    APHASE = NA_character_,
    APHASEN = NA_integer_,
    BASETYPE = "DOUBLE_BLINDED"
  )

# Assign AVAL, AVALC, AVALU and DTYPE for each ISTESTCD and ISBDAGNT

is_aval <- is_basetype %>%
  mutate(
    MRT = case_when(
      ADATYPE == "ADA_BAB" & ADAPARM == "XANOMELINE" ~ 1.4,
      ADATYPE == "ADA_BAB" & ADAPARM == "OTHER_DRUG" ~ 9.9,
      TRUE ~ NA_real_
    ),
    DTL = case_when(
      ADATYPE == "ADA_BAB" & ADAPARM == "XANOMELINE" ~ 999,
      ADATYPE == "ADA_BAB" & ADAPARM == "OTHER_DRUG" ~ 888,
      TRUE ~ NA_real_
    ),
    RESULTC = case_when(
      toupper(ISSTRESC) %in% c(
        "NEGATIVE", "NEGATIVE SCREEN", "NEGATIVE IMMUNODEPLETION",
        "NEGATIVE CONFIRMATION"
      ) ~ "NEGATIVE",
      toupper(ISSTRESC) %in% c(
        "NEGATIVE TITER", "<1.70", "< 1.70", "<1.30", "< 1.30", "<1.40",
        "POSITIVE IMMUNODEPLETION", "POSITIVE CONFIRMATION", "POSITIVE"
      ) ~ "POSITIVE",
      ISSTRESN > 0 ~ "POSITIVE",
      TRUE ~ NA_character_
    ),
    RESULTN = case_when(
      toupper(RESULTC) == "POSITIVE" ~ 1,
      toupper(RESULTC) == "NEGATIVE" ~ 0,
      TRUE ~ NA_integer_
    ),
    AVAL = case_when(
      ADATYPE == "ADA_BAB" & toupper(RESULTC) == "POSITIVE" & !is.na(ISSTRESN) ~ ISSTRESN,
      ADATYPE == "ADA_BAB" & toupper(RESULTC) == "POSITIVE" & is.na(ISSTRESN) & !is.na(MRT) ~ MRT,
      TRUE ~ NA_real_
    ),
    AVALC = case_when(
      # NABSTAT gets ISSTRESC, Standard ADA is set to NA as AVAL are numeric original results
      ADATYPE == "ADA_NAB" ~ ISSTRESC,
      TRUE ~ NA_character_
    ),
    AVALU = case_when(
      ADATYPE == "ADA_NAB" ~ ISSTRESU,
      ADATYPE == "ADA_BAB" & toupper(RESULTC) == "POSITIVE" & !is.na(ISSTRESN) ~ ISSTRESU,
      ADATYPE == "ADA_BAB" & toupper(RESULTC) == "POSITIVE" & is.na(ISSTRESN) & !is.na(MRT)
      ~ "titer",
      TRUE ~ NA_character_
    ),
    DTYPE = case_when(
      ADATYPE == "ADA_BAB" & toupper(RESULTC) == "POSITIVE" & is.na(ISSTRESN) & !is.na(MRT) ~ "MRT",
      TRUE ~ NA_character_
    )
  )

# Begin computation of parameters -----------------------------------------

# Identify Best Baseline for each analyte and parameter type
# Baseline is NFRLT <= 0 or Unscheduled AND the ADA Date is on or before the date of first dose.

is_baseline <- is_aval %>%
  # Calculate ABLFL. If more than one record for the 'order' and 'filter' will throw a duplicate
  # record warning, user can decide how to adjust the mode, order, filter or adjust in prior step.
  restrict_derivation(
    derivation = derive_var_extreme_flag,
    args = params(
      by_vars = exprs(STUDYID, USUBJID, BASETYPE, ADATYPE, ADAPARM),
      order = exprs(ADT, NFRLT),
      new_var = ABLFL,
      mode = "last"
    ),
    # flag baseline based on time not values due to ADA parameters can in some cases permit missing.
    # If special visits (i.e. > 60000) are eligible as baselines, include those
    filter = ((NFRLT <= 0 | NFRLT > 60000) & (ADT <= FANLDT) & !is.na(BASETYPE))
  ) %>%
  mutate(
    # VALID flags for use later as applicable:
    # VALIDBASE flags non-missing values on baseline (by each ADATYPE and ADAPARM)
    #   Note: VALIDBASE is not used as this template allows a baseline to be valid as
    #        as long as its present (can be missing), adapt as needed.
    # VALIDPOST flags non-missing values on post-baseline (by each ADATYPE and ADAPARM)
    VALIDBASE = case_when(
      ABLFL == "Y" & (!is.na(AVALC) | !is.na(RESULTC) | !is.na(AVAL)) ~ "Y",
      ABLFL == "Y" & (is.na(AVALC) & is.na(RESULTC) & is.na(AVAL)) ~ "N",
      TRUE ~ NA_character_
    ),
    VALIDPOST = case_when(
      ADTM > FANLDTM & is.na(ABLFL) & (!is.na(RESULTC) | !is.na(AVAL)) ~ "Y",
      ADTM > FANLDTM & is.na(ABLFL) & (is.na(RESULTC) & is.na(AVAL)) ~ "N",
      TRUE ~ NA_character_
    )
  )

# Compute BASE and CHG
is_aval_change <- is_baseline %>%
  derive_var_base(
    by_vars = exprs(STUDYID, USUBJID, BASETYPE, ADATYPE, ADAPARM),
    source_var = AVAL,
    new_var = BASE,
    filter = ABLFL == "Y"
  ) %>%
  restrict_derivation(
    derivation = derive_var_chg,
    filter = is.na(ABLFL)
  )

# Interpreted Result Baseline
is_result_change <- is_aval_change %>%
  derive_var_base(
    by_vars = exprs(STUDYID, USUBJID, BASETYPE, ADATYPE, ADAPARM),
    source_var = RESULTN,
    new_var = BASE_RESULT,
    filter = ABLFL == "Y"
  )

# Get base only data for use later
base_data <- is_result_change %>%
  filter(ABLFL == "Y") %>%
  select(
    STUDYID, USUBJID, DRUG, BASETYPE, ADATYPE, ADAPARM, BASE_RESULT, BASE,
    ABLFL
  )

# Assign and save ADABLPFL for later use
adablpfl <- is_result_change %>%
  filter(ABLFL == "Y") %>%
  distinct(STUDYID, USUBJID, DRUG, BASETYPE,
    ADATYPE, ADAPARM,
    .keep_all = TRUE
  ) %>%
  select(STUDYID, USUBJID, DRUG, BASETYPE, ADATYPE, ADAPARM, ABLFL) %>%
  rename(ADABLPFL = ABLFL)

# Calculate the By Visit parameters
is_visit_flags <- is_result_change %>%
  mutate(
    TFLAGV = case_when(
      VALIDPOST == "Y" & is.na(ABLFL) & ADATYPE == "ADA_BAB" & (BASE_RESULT == 1 & (CHG >= 0.6))
      ~ 2,
      VALIDPOST == "Y" & is.na(ABLFL) & ADATYPE == "ADA_BAB" & BASE_RESULT == 1 &
        (((CHG < 0.6) & !is.na(AVAL)) | (RESULTN == 0)) ~ 3,
      VALIDPOST == "Y" & is.na(ABLFL) & ADATYPE == "ADA_BAB" & ((BASE_RESULT == 0 |
        is.na(BASE_RESULT)) & RESULTN == 0) ~ 0,
      VALIDPOST == "Y" & is.na(ABLFL) & ADATYPE == "ADA_BAB" & ((BASE_RESULT == 0 |
        is.na(BASE_RESULT)) & RESULTN == 1) ~ 1,
      TRUE ~ NA_integer_
    ),
    PBFLAGV = case_when(
      !is.na(TFLAGV) & TFLAGV %in% c(1, 2) ~ 1,
      !is.na(TFLAGV) & TFLAGV %in% c(0, 3) ~ 0,
      TRUE ~ NA_integer_
    ),
    ADASTATV = case_when(
      !is.na(PBFLAGV) & PBFLAGV == 1 ~ "ADA+",
      !is.na(PBFLAGV) & PBFLAGV == 0 ~ "ADA-",
      VALIDPOST == "Y" & is.na(ABLFL) & ADATYPE == "ADA_BAB" & is.na(PBFLAGV) ~ "MISSING",
      TRUE ~ "MISSING"
    ),
  )

# These next code segments create utility datasets that will
# then get merged back into the main dataset (is_visit_flags)

# Post baseline must be valid post data (result not missing)
post_data <- is_visit_flags %>%
  filter(VALIDPOST == "Y") %>%
  select(
    STUDYID, USUBJID, DRUG, BASETYPE, ADATYPE, ADAPARM, RESULTN,
    AVAL, ADTM, CHG
  ) %>%
  rename(AVAL_P = AVAL) %>%
  rename(RESULT_P = RESULTN)

# Use "post_data" to make a ADPBLPFL flag data set to merge back later in the program
#  Note:  "post_data" is if VALIDPOST="Y" (has a non missing post baseline result)
adpblpfl <- post_data %>%
  distinct(STUDYID, USUBJID, DRUG, BASETYPE, ADATYPE,
    ADAPARM,
    .keep_all = TRUE
  ) %>%
  select(STUDYID, USUBJID, DRUG, BASETYPE, ADATYPE, ADAPARM) %>%
  mutate(
    ADPBLPFL = "Y"
  )

# Compute BFLAG, TFLAG, PBFLAG ----

most_post_result <- post_data %>%
  group_by(STUDYID, USUBJID, DRUG, BASETYPE, ADATYPE, ADAPARM) %>%
  summarize(RESULT_P = max(RESULT_P)) %>%
  ungroup()

most_post_aval <- post_data %>%
  filter(!is.na(AVAL_P)) %>%
  group_by(STUDYID, USUBJID, BASETYPE, ADATYPE, ADAPARM) %>%
  summarize(AVAL_P = max(AVAL_P)) %>%
  ungroup()

most_post_chg <- post_data %>%
  filter(!is.na(AVAL_P)) %>%
  group_by(STUDYID, USUBJID, BASETYPE, ADATYPE, ADAPARM) %>%
  summarize(MAXCHG = max(CHG)) %>%
  ungroup()

# Merge the most_post_result, most_post_aval and most_post_chg into one utility data set
most_post <- most_post_result %>%
  derive_vars_merged(
    dataset_add = most_post_aval,
    by_vars = exprs(STUDYID, USUBJID, BASETYPE, ADATYPE, ADAPARM)
  ) %>%
  derive_vars_merged(
    dataset_add = most_post_chg,
    by_vars = exprs(STUDYID, USUBJID, BASETYPE, ADATYPE, ADAPARM)
  )

# Use an outer Join to combine baseline with most post results create utility dataset with flag data
flagdata_init <- full_join(base_data, most_post, by = c(
  "DRUG", "STUDYID", "USUBJID", "BASETYPE",
  "ADATYPE", "ADAPARM"
)) %>%
  mutate(
    BFLAG =
      case_when(
        BASE_RESULT == 0 ~ 0,
        BASE_RESULT == 1 ~ 1,
        TRUE ~ NA_integer_
      ),
    TFLAG =
      case_when(
        (BFLAG == 0 | is.na(BFLAG)) & RESULT_P == 0 ~ 0,
        (BFLAG == 0 | is.na(BFLAG)) & RESULT_P == 1 ~ 1,
        BFLAG == 1 & (MAXCHG >= 0.6) ~ 2,
        BFLAG == 1 & !is.na(MAXCHG) & (MAXCHG < 0.6) ~ 3,
        BFLAG == 1 & is.na(MAXCHG) & RESULT_P == 0 ~ 3,
        TRUE ~ NA_integer_
      ),
    PBFLAG =
      case_when(
        ADATYPE == "ADA_BAB" & (TFLAG == 1 | TFLAG == 2) ~ 1,
        ADATYPE == "ADA_BAB" & (TFLAG == 0 | TFLAG == 3) ~ 0,
        ADATYPE == "ADA_NAB" & RESULT_P == 0 ~ 0,
        ADATYPE == "ADA_NAB" & RESULT_P == 1 ~ 1,
        TRUE ~ NA_integer_
      ),
    ADASTAT = case_when(
      ADATYPE == "ADA_BAB" & PBFLAG == 1 ~ 1,
      ADATYPE == "ADA_BAB" & PBFLAG == 0 ~ 0,
      TRUE ~ NA_integer_
    )
  )

# Get overall ADASTAT status for computing NABSTAT
flagdata_adastat <- flagdata_init %>%
  derive_vars_merged(
    dataset_add = flagdata_init,
    filter_add = ADATYPE == "ADA_BAB",
    new_vars = exprs(ADASTAT_MAIN = ADASTAT),
    by_vars = exprs(STUDYID, USUBJID, DRUG)
  )

# Compute NABPOSTMISS onto flag_data for NABSTAT
flagdata_nab <- flagdata_adastat %>%
  derive_var_merged_exist_flag(
    dataset_add = is_visit_flags,
    new_var = NABPOSTMISS,
    condition = is.na(ABLFL) & VALIDPOST == "N" & ADATYPE == "ADA_NAB",
    by_vars = exprs(STUDYID, USUBJID, BASETYPE, ADATYPE, ADAPARM)
  )

# Compute NAB Stat using both methods
# Note: Option 2 added as a placekeeper starter method, adjust as needed for
#   study specific specs.
flagdata_final <- flagdata_nab %>%
  mutate(
    # For Option 1, if any post baseline NAB were blank without a Positive (NABPOSTMISS = Y),
    # NABSTAT = missing
    nabstat_opt1 = case_when(
      # Based on ADASTAT (EMERNEG/EMERPOS) and NAB results
      ADATYPE == "ADA_NAB" & ADASTAT_MAIN == 1 & RESULT_P == 1 ~ 1,
      ADATYPE == "ADA_NAB" & ADASTAT_MAIN == 1 & RESULT_P == 0 &
        (is.na(NABPOSTMISS) | NABPOSTMISS == "N") ~ 0,
      TRUE ~ NA_integer_
    ),
    nabstat_opt2 = case_when(
      # Based on ADASTAT (EMERNEG/EMERPOS) Only.
      ADATYPE == "ADA_NAB" & ADASTAT_MAIN == 1 ~ 1,
      ADATYPE == "ADA_NAB" & ADASTAT_MAIN == 0 ~ 0,
      TRUE ~ NA_integer_
    ),
  ) %>%
  # Drop variables no longer needed from flag_data before merging with main ADAB
  select(-BASE, -BASE_RESULT, -AVAL_P, -RESULT_P, -DRUG, -ABLFL)

# Put TFLAG, BFLAG and PBFLAG and the nabstat_ variables onto the main
# dataset (is_flagdata from is_visit_flags)
is_flagdata <- is_visit_flags %>%
  # main_aab_flagdata
  derive_vars_merged(
    dataset_add = flagdata_final,
    by_vars = exprs(STUDYID, USUBJID, BASETYPE, ADATYPE, ADAPARM)
  )

# Create a utility data set to compute PERSADA, TRANADA, INDUCED, ENHANCED related parameters

per_tran_pre <- is_flagdata %>%
  filter(VALIDPOST == "Y") %>%
  select(
    STUDYID, USUBJID, BASETYPE, ADATYPE, ADAPARM, ADTM, FANLDTM, FANLDT,
    TFLAGV, ISSEQ
  ) %>%
  group_by(STUDYID, USUBJID, BASETYPE, ADATYPE, ADAPARM) %>%
  mutate(
    MaxADTM = max(ADTM)
  ) %>%
  ungroup() %>%
  mutate(
    LFLAGPOS = case_when(
      ADTM == MaxADTM & (TFLAGV == 1) ~ 1,
      ADTM == MaxADTM & (TFLAGV == 2) ~ 1,
      TRUE ~ NA_integer_
    )
  )

# Keep TFLAGV = 1 or TFLAGV = 2 (Any Treatment Emergent)
# Regular Induced PERSADA and TRANADA will later be based on TFLAGV = 1
# TFLAGV = 2 can use for separate PERSADA/TRANADA based on Enhanced
per_tran_all <- per_tran_pre %>%
  filter(TFLAGV == 1 | TFLAGV == 2) %>%
  select(-MaxADTM, -ISSEQ) %>%
  group_by(STUDYID, USUBJID, BASETYPE, ADATYPE, ADAPARM) %>%
  mutate(
    FPPDTM = min(ADTM),
    LPPDTM = max(ADTM)
  ) %>%
  ungroup()

per_tran_inc_last <- per_tran_all %>%
  group_by(STUDYID, USUBJID, BASETYPE, ADATYPE, ADAPARM) %>%
  summarize(COUNT_INC_LAST = n()) %>%
  ungroup()

per_tran_exc_last <- per_tran_all %>%
  filter(is.na(LFLAGPOS)) %>%
  group_by(STUDYID, USUBJID, BASETYPE, ADATYPE, ADAPARM) %>%
  summarize(COUNT_EXC_LAST = n()) %>%
  ungroup()

# Reduce "per_tran_all" to one record per by_vars using ADTM = LPPDTM to get best last records
#   Then Merge in the Include Last and Exclude Last Flags
per_tran_last <- per_tran_all %>%
  filter(ADTM == LPPDTM) %>%
  derive_vars_merged(
    dataset_add = per_tran_inc_last,
    by_vars = exprs(STUDYID, USUBJID, BASETYPE, ADATYPE, ADAPARM)
  ) %>%
  derive_vars_merged(
    dataset_add = per_tran_exc_last,
    by_vars = exprs(STUDYID, USUBJID, BASETYPE, ADATYPE, ADAPARM)
  )

# Compute final parameters
per_tran_final <- per_tran_last %>%
  mutate(
    FPPDT = as_date(FPPDTM),
    LPPDT = as_date(LPPDTM),
    ADADUR = case_when(
      LPPDTM - FPPDTM == 0 & (!is.na(LPPDTM) & !is.na(FPPDTM)) ~ 1 / 7,
      !is.na(LPPDTM) & !is.na(FPPDTM)
      ~ ((as.numeric(difftime(LPPDTM, FPPDTM, units = "secs")) / (60 * 60 * 24)) + 1) / 7,
      !is.na(LPPDT) & !is.na(FPPDT)
      ~ (as.numeric(LPPDT - FPPDT) + 1) / 7,
      TRUE ~ NA_real_
    ),
    TIMADA = case_when(
      !is.na(FPPDTM) & !is.na(FANLDTM) ~ as.numeric(difftime(FPPDTM, FANLDTM, units = "weeks")),
      !is.na(FPPDT) & !is.na(FANLDT) ~ (as.numeric(FPPDT - FANLDT) + 1) / 7,
      TRUE ~ NA_real_
    ),
    tdur = case_when(
      !is.na(LPPDTM) & !is.na(FPPDTM) ~
        as.numeric(difftime(LPPDTM, FPPDTM, units = "secs")) / (7 * 3600 * 24),
      !is.na(LPPDT) & !is.na(FPPDT) ~ as.numeric(LPPDT - FPPDT + 1) / 7,
      TRUE ~ NA_real_
    ),
    # Standard TRANADA and PERSADA based on TFLAGV = 1 (Induced)
    TRANADA = case_when(
      TFLAGV == 1 & ((COUNT_EXC_LAST == 1 | (COUNT_INC_LAST >= 2 & tdur < 16)) & is.na(LFLAGPOS)) ~ 1,
      TRUE ~ NA_integer_
    ),
    PERSADA = case_when(
      TFLAGV == 1 & (COUNT_INC_LAST == 1 & (COUNT_EXC_LAST <= 0 | is.na(COUNT_EXC_LAST)) |
        (COUNT_INC_LAST >= 2 & tdur >= 16) | LFLAGPOS == 1) ~ 1,
      TRUE ~ NA_integer_
    ),
    # These TRANADAE and PERSADAE based on TFLAGV = 1 or TFLAGV=2 (Induced or Enhanced)
    TRANADAE = case_when(
      TFLAGV >= 1 & ((COUNT_EXC_LAST == 1 | (COUNT_INC_LAST >= 2 & tdur < 16)) & is.na(LFLAGPOS)) ~ 1,
      TRUE ~ NA_integer_
    ),
    PERSADAE = case_when(
      TFLAGV >= 1 & (COUNT_INC_LAST == 1 & (COUNT_EXC_LAST <= 0 | is.na(COUNT_EXC_LAST)) |
        (COUNT_INC_LAST >= 2 & tdur >= 16) | LFLAGPOS == 1) ~ 1,
      TRUE ~ NA_integer_
    ),
    INDUCED = case_when(
      TFLAGV == 1 ~ "Y",
      TRUE ~ "N"
    ),
    ENHANCED = case_when(
      TFLAGV == 2 ~ "Y",
      TRUE ~ "N"
    )
  ) %>%
  # Drop temporary variables that do not need to be merged into main ADAB
  select(-ADTM, -TFLAGV, -FANLDTM, -FANLDT, -LFLAGPOS, -FPPDT, -LPPDT)

# Put PERSADA, TRANADA, INDUCED, ENHANCED, TDUR, ADADUR onto "is_flagdata" as "main_aab_pertran"
# Note: signal_duplicate_records() error usually occurs when a subject has duplicate
# records for a given BASETYPE, ADATYPE, ADAPARM and ISDTC. Investigate then add code
# to filter it down to best one record per USUBJID, BASETYPE, ADATYPE and ADAPARM
main_aab_pertran <- is_flagdata %>%
  derive_vars_merged(
    dataset_add = per_tran_final,
    by_vars = exprs(STUDYID, USUBJID, BASETYPE, ADATYPE, ADAPARM)
  )

main_aab_rtimes <- main_aab_pertran %>%
  mutate(
    ATPT = ISTPT,
    ADADUR = round(ADADUR, digits = 4),
    TIMADA = round(TIMADA, digits = 4),
    PERSADA = case_when(
      is.na(PERSADA) ~ 0,
      TRUE ~ PERSADA
    ),
    PERSADAE = case_when(
      is.na(PERSADAE) ~ 0,
      TRUE ~ PERSADAE
    ),
    TRANADA = case_when(
      is.na(TRANADA) ~ 0,
      TRUE ~ TRANADA
    ),
    TRANADAE = case_when(
      is.na(TRANADAE) ~ 0,
      TRUE ~ TRANADAE
    ),
    NABSTAT = nabstat_opt1
  )

# Merge ADABLPFL and ADPBFPFL onto main dataset.
main_aab <- main_aab_rtimes %>%
  derive_vars_merged(
    dataset_add = adablpfl,
    new_vars = exprs(ADABLPFL),
    by_vars = exprs(STUDYID, USUBJID, BASETYPE, ADATYPE, ADAPARM)
  ) %>%
  derive_vars_merged(
    dataset_add = adpblpfl,
    new_vars = exprs(ADPBLPFL),
    by_vars = exprs(STUDYID, USUBJID, BASETYPE, ADATYPE, ADAPARM)
  )

# Begin Creation of each PARAM for the final ADAB format using main_aab --------

#  First create "core_aab" with PARCAT1 to be the input for all the parameter sub-assemblies
core_aab <- main_aab %>%
  mutate(
    # SDTM V1.x: Assign PARAM from ISTEST or customize
    PARCAT1 = ISTEST,
    # SDTM V2.x: Assign PARAM using text values plus ADAPARM
    PARCAT1 = case_when(
      ADATYPE == "ADA_BAB" ~ paste("Anti-", ADAPARM, " Antibody", sep = ""),
      ADATYPE == "ADA_NAB" ~ paste("Anti-", ADAPARM, " Neutralizing Antibody", sep = ""),
      TRUE ~ NA_character_
    ),
    # Initialize PARAMCD and PARAM
    PARAMCD = ADAPARM,
    PARAM = PARCAT1
  )

# By Visit Primary ADA ISTESTCD Titer Results
adab_titer <- core_aab %>%
  filter(ADATYPE == "ADA_BAB") %>%
  mutate(
    # For ADASTAT, append "Titer Units" to PARAM
    PARAM = paste(PARCAT1, ", Titer Units", sep = ""),
  )

# By Visit NAB ISTESTCD Results
adab_nabvis <- core_aab %>%
  filter(ADATYPE == "ADA_NAB") %>%
  # These two flags, BASE and CHG are only kept on the primary ADA ISTESTCD test
  select(-BASE, -CHG, -ADABLPFL, -ADPBLPFL)

# By Visit Main ADA Titer Interpreted RESULT data
adab_result <- core_aab %>%
  filter(ADATYPE == "ADA_BAB") %>%
  mutate(
    PARAMCD = "RESULTy",
    PARAM = paste("ADA interpreted per sample result,", PARCAT1, sep = " "),
    AVALC = toupper(RESULTC),
    AVAL = case_when(
      AVALC == "NEGATIVE" ~ 0,
      AVALC == "POSITIVE" ~ 1,
      TRUE ~ NA_integer_
    ),
    AVALU = NA_character_
  ) %>%
  # DTYPE is only kept on the primary ISTESTCD parameter, drop then recompute final BASE and CHG
  select(-DTYPE, -BASE, -CHG)

# Derive BASE and Calculate Change from Baseline on BAB RESULT records
adab_result <- adab_result %>%
  derive_var_base(
    by_vars = exprs(STUDYID, USUBJID, BASETYPE, ADATYPE, ADAPARM),
    source_var = AVAL,
    new_var = BASE,
    filter = ABLFL == "Y"
  ) %>%
  restrict_derivation(
    derivation = derive_var_chg,
    filter = is.na(ABLFL)
  )

# By Visit NAB Interpreted RESULT data
adab_nabres <- core_aab %>%
  filter(ADATYPE == "ADA_NAB") %>%
  mutate(
    PARAMCD = "RESULTy",
    PARAM = paste("Nab interpreted per sample result,", PARCAT1, sep = " "),
    AVALC = toupper(RESULTC),
    AVAL = case_when(
      AVALC == "NEGATIVE" ~ 0,
      AVALC == "POSITIVE" ~ 1,
      TRUE ~ NA_integer_
    ),
    AVALU = NA_character_
  ) %>%
  # These two flags are only kept on the primary ADA test, drop then recompute final BASE and CHG
  select(-ADABLPFL, -ADPBLPFL, -BASE, -CHG)

# Derive BASE and Calculate Change from Baseline on NAB RESULT records
adab_nabres <- adab_nabres %>%
  derive_var_base(
    by_vars = exprs(STUDYID, USUBJID, BASETYPE, ADATYPE, ADAPARM),
    source_var = AVAL,
    new_var = BASE,
    filter = ABLFL == "Y"
  ) %>%
  restrict_derivation(
    derivation = derive_var_chg,
    filter = is.na(ABLFL)
  )

# By Visit Titer TFLAGV data -----------------------

# Note: By Visit parameters do not keep CHG, MRT and DTL variables
adab_tflagv <- core_aab %>%
  filter(ADATYPE == "ADA_BAB") %>%
  mutate(
    PARAMCD = "TFLAGV",
    PARAM = paste("Treatment related ADA by Visit,", PARCAT1, sep = " "),
    AVAL = TFLAGV,
    AVALC = as.character(AVAL),
    AVALU = NA_character_
  ) %>%
  # Drop BASE, CHG, MRT and DTL and the two --FL flags (are only kept on the primary ADA test)
  select(-BASE, -CHG, -MRT, -DTL, -ADABLPFL, -ADPBLPFL)

# By Visit Titer PBFLAGV data ---------------------
# Note: By Visit parameters do not keep CHG, MRT and DTL variables
adab_pbflagv <- core_aab %>%
  filter(ADATYPE == "ADA_BAB") %>%
  mutate(
    PARAMCD = "PBFLAGV",
    PARAM = paste("Post Baseline Pos/Neg by Visit,", PARCAT1, sep = " "),
    AVALC = case_when(
      PBFLAGV == 1 | PBFLAGV == 2 ~ "POSITIVE",
      PBFLAGV == 0 | PBFLAGV == 3 ~ "NEGATIVE",
      TRUE ~ "MISSING"
    ),
    AVAL = case_when(
      AVALC == "POSITIVE" ~ 1,
      AVALC == "NEGATIVE" ~ 0,
      TRUE ~ NA_integer_
    ),
    AVALU = NA_character_
  ) %>%
  # Drop BASE, CHG, MRT and DTL and the two --FL flags (are only kept on the primary ADA test)
  select(-BASE, -CHG, -MRT, -DTL, -ADABLPFL, -ADPBLPFL)

# By Visit Titer ADASTATV data
# Note: By Visit parameters do not keep CHG, MRT and DTL variables
adab_adastatv <- core_aab %>%
  filter(ADATYPE == "ADA_BAB") %>%
  mutate(
    PARAMCD = "ADASTATV",
    PARAM = paste("ADA Status of a patient by Visit,", PARCAT1, sep = " "),
    AVALC = ADASTATV,
    AVAL = case_when(
      AVALC == "ADA+" ~ 1,
      AVALC == "ADA-" ~ 0,
      TRUE ~ NA_integer_
    ),
    AVALU = NA_character_
  ) %>%
  # Drop BASE, CHG, MRT and DTL and the two --FL flags (are only kept on the primary ADA test)
  select(-BASE, -CHG, -MRT, -DTL, -ADABLPFL, -ADPBLPFL)

# Next below are the individual params
# assign the AVISIT and AVISITN for individual params

core_aab <- core_aab %>%
  mutate(
    AVISIT = overall_avisit,
    AVISITN = overall_avisitn
  )

# Get Patient flag BFLAG
adab_bflag <- core_aab %>%
  filter(ADATYPE == "ADA_BAB") %>%
  distinct(
    STUDYID, USUBJID, PARCAT1, BASETYPE, PARAM, ADATYPE, ADAPARM,
    ISTESTCD, ISTEST, ISCAT, ISBDAGNT, AVISITN, AVISIT, ISSPEC, BFLAG
  ) %>%
  mutate(
    PARAMCD = "BFLAGy",
    PARAM = paste("Baseline Pos/Neg,", PARCAT1, sep = " "),
    AVAL = BFLAG,
    AVALC = case_when(
      AVAL == 1 ~ "Y",
      AVAL == 0 ~ "N",
      TRUE ~ NA_character_
    ),
    AVALU = NA_character_
  )

# Get Patient flag INDUCD
adab_incucd <- core_aab %>%
  filter(ADATYPE == "ADA_BAB") %>%
  distinct(
    STUDYID, USUBJID, PARCAT1, BASETYPE, PARAM, ADATYPE, ADAPARM,
    ISTESTCD, ISTEST, ISCAT, ISBDAGNT, AVISITN, AVISIT, ISSPEC, TFLAG
  ) %>%
  mutate(
    PARAMCD = "INDUCDy",
    PARAM = paste("Treatment induced ADA,", PARCAT1, sep = " "),
    AVAL = case_when(
      TFLAG == 1 ~ 1,
      TRUE ~ 0
    ),
    AVALC = case_when(
      AVAL == 1 ~ "Y",
      AVAL == 0 ~ "N",
      TRUE ~ NA_character_
    ),
    AVALU = NA_character_
  )

# Get Patient flag ENHANC
adab_enhanc <- core_aab %>%
  filter(ADATYPE == "ADA_BAB") %>%
  distinct(
    STUDYID, USUBJID, PARCAT1, BASETYPE, PARAM, ADATYPE, ADAPARM,
    ISTESTCD, ISTEST, ISCAT, ISBDAGNT, AVISITN, AVISIT, ISSPEC, TFLAG
  ) %>%
  mutate(
    PARAMCD = "ENHANCy",
    PARAM = paste("Treatment enhanced ADA,", PARCAT1, sep = " "),
    AVAL = case_when(
      TFLAG == 2 ~ 1,
      TRUE ~ 0
    ),
    AVALC = case_when(
      AVAL == 1 ~ "Y",
      AVAL == 0 ~ "N",
      TRUE ~ NA_character_
    ),
    AVALU = NA_character_
  )

# Get Patient flag EMERPOS
adab_emerpos <- core_aab %>%
  filter(ADATYPE == "ADA_BAB") %>%
  distinct(
    STUDYID, USUBJID, PARCAT1, BASETYPE, PARAM, ADATYPE, ADAPARM,
    ISTESTCD, ISTEST, ISCAT, ISBDAGNT, AVISITN, AVISIT, ISSPEC, PBFLAG
  ) %>%
  mutate(
    PARAMCD = "EMERPOSy",
    PARAM = paste("Treatment Emergent - Positive,", PARCAT1, sep = " "),
    AVAL = case_when(
      PBFLAG == 1 ~ 1,
      TRUE ~ 0
    ),
    AVALC = case_when(
      AVAL == 1 ~ "Y",
      AVAL == 0 ~ "N",
      TRUE ~ NA_character_
    ),
    AVALU = NA_character_
  )

# Get Patient flag TRUNAFF
adab_trunaff <- core_aab %>%
  filter(ADATYPE == "ADA_BAB") %>%
  distinct(
    STUDYID, USUBJID, PARCAT1, BASETYPE, PARAM, ADATYPE, ADAPARM,
    ISTESTCD, ISTEST, ISCAT, ISBDAGNT, AVISITN, AVISIT, ISSPEC, TFLAG
  ) %>%
  mutate(
    PARAMCD = "TRUNAFFy",
    PARAM = paste("Treatment unaffected,", PARCAT1, sep = " "),
    AVAL = case_when(
      TFLAG == 3 ~ 1,
      TRUE ~ 0
    ),
    AVALC = case_when(
      AVAL == 1 ~ "Y",
      AVAL == 0 ~ "N",
      TRUE ~ NA_character_
    ),
    AVALU = NA_character_
  )

# Get Patient flag EMERNEG
adab_emerneg <- core_aab %>%
  filter(ADATYPE == "ADA_BAB") %>%
  distinct(
    STUDYID, USUBJID, PARCAT1, BASETYPE, PARAM, ADATYPE, ADAPARM,
    ISTESTCD, ISTEST, ISCAT, ISBDAGNT, AVISITN, AVISIT, ISSPEC, PBFLAG
  ) %>%
  mutate(
    PARAMCD = "EMERNEGy",
    PARAM = paste("Treatment Emergent - Negative,", PARCAT1, sep = " "),
    AVAL = case_when(
      PBFLAG == 0 ~ 1,
      TRUE ~ 0
    ),
    AVALC = case_when(
      AVAL == 1 ~ "Y",
      AVAL == 0 ~ "N",
      TRUE ~ NA_character_
    ),
    AVALU = NA_character_
  )

# Get Patient flag NOTRREL
adab_notrrel <- core_aab %>%
  filter(ADATYPE == "ADA_BAB") %>%
  distinct(
    STUDYID, USUBJID, PARCAT1, BASETYPE, PARAM, ADATYPE, ADAPARM,
    ISTESTCD, ISTEST, ISCAT, ISBDAGNT, AVISITN, AVISIT, ISSPEC, PBFLAG
  ) %>%
  mutate(
    PARAMCD = "NOTRRELy",
    PARAM = paste("No treatment related ADA,", PARCAT1, sep = " "),
    AVAL = case_when(
      is.na(PBFLAG) ~ 1,
      TRUE ~ 0
    ),
    AVALC = case_when(
      AVAL == 1 ~ "Y",
      AVAL == 0 ~ "N",
      TRUE ~ NA_character_
    ),
    AVALU = NA_character_
  )

# Get Patient flag ADASTAT
adab_adastat <- core_aab %>%
  filter(ADATYPE == "ADA_BAB") %>%
  distinct(
    STUDYID, USUBJID, PARCAT1, BASETYPE, PARAM, ADATYPE, ADAPARM,
    ISTESTCD, ISTEST, ISCAT, ISBDAGNT, AVISITN, AVISIT, ISSPEC, ADASTAT
  ) %>%
  mutate(
    PARAMCD = "ADASTATy",
    PARAM = paste("ADA Status of a patient,", PARCAT1, sep = " "),
    AVAL = ADASTAT,
    AVALC = case_when(
      AVAL == 1 ~ "POSITIVE",
      AVAL == 0 ~ "NEGATIVE",
      TRUE ~ "MISSING"
    ),
    AVALU = NA_character_
  )

# Get Patient flag TIMADA
adab_timada <- core_aab %>%
  filter(ADATYPE == "ADA_BAB") %>%
  distinct(
    STUDYID, USUBJID, PARCAT1, BASETYPE, PARAM, ADATYPE, ADAPARM,
    ISTESTCD, ISTEST, ISCAT, ISBDAGNT, AVISITN, AVISIT, ISSPEC, TIMADA
  ) %>%
  mutate(
    PARAMCD = "TIMADAy",
    PARAM = paste("Time to onset of ADA (Weeks),", PARCAT1, sep = " "),
    AVAL = TIMADA,
    AVALC = NA_character_,
    AVALU = if_else(!is.na(AVAL), "WEEKS", NA_character_)
  )

# Get Patient flag PERSADA
adab_persada <- core_aab %>%
  filter(ADATYPE == "ADA_BAB") %>%
  distinct(
    STUDYID, USUBJID, PARCAT1, BASETYPE, PARAM, ADATYPE, ADAPARM,
    ISTESTCD, ISTEST, ISCAT, ISBDAGNT, AVISITN, AVISIT, ISSPEC, PERSADA
  ) %>%
  mutate(
    PARAMCD = "PERSADAy",
    PARAM = paste("Persistent ADA,", PARCAT1, sep = " "),
    AVAL = PERSADA,
    AVALC = case_when(
      AVAL == 1 ~ "Y",
      AVAL == 0 ~ "N",
      TRUE ~ NA_character_
    ),
    AVALU = NA_character_
  )

# Get Patient flag TRANADA
adab_tranada <- core_aab %>%
  filter(ADATYPE == "ADA_BAB") %>%
  distinct(
    STUDYID, USUBJID, PARCAT1, BASETYPE, PARAM, ADATYPE, ADAPARM,
    ISTESTCD, ISTEST, ISCAT, ISBDAGNT, AVISITN, AVISIT, ISSPEC, TRANADA
  ) %>%
  mutate(
    PARAMCD = "TRANADAy",
    PARAM = paste("Transient ADA,", PARCAT1, sep = " "),
    AVAL = TRANADA,
    AVALC = case_when(
      AVAL == 1 ~ "Y",
      AVAL == 0 ~ "N",
      TRUE ~ NA_character_
    ),
    AVALU = NA_character_
  )

# Get Patient flag NABSTAT
adab_nabstat <- core_aab %>%
  filter(ADATYPE == "ADA_NAB") %>%
  distinct(
    STUDYID, USUBJID, PARCAT1, BASETYPE, PARAM, ADATYPE, ADAPARM,
    ISTESTCD, ISTEST, ISCAT, ISBDAGNT, AVISITN, AVISIT, ISSPEC, NABSTAT
  ) %>%
  mutate(
    PARAMCD = "NABSTATy",
    PARAM = paste("Nab Status,", PARCAT1, sep = " "),
    AVAL = NABSTAT,
    AVALC = case_when(
      AVAL == 1 ~ "POSITIVE",
      AVAL == 0 ~ "NEGATIVE",
      TRUE ~ "MISSING"
    ),
    AVALU = NA_character_
  )

# Get Patient flag ADADUR
adab_adadur <- core_aab %>%
  filter(ADATYPE == "ADA_BAB") %>%
  distinct(
    STUDYID, USUBJID, PARCAT1, BASETYPE, PARAM, ADATYPE, ADAPARM,
    ISTESTCD, ISTEST, ISCAT, ISBDAGNT, AVISITN, AVISIT, ISSPEC, ADADUR
  ) %>%
  mutate(
    PARAMCD = "ADADURy",
    PARAM = paste("ADA Duration (Weeks),", PARCAT1, sep = " "),
    AVAL = ADADUR,
    AVALC = NA_character_,
    AVALU = if_else(!is.na(AVAL), "WEEKS", NA_character_)
  )

# Get Patient flag FPPDTM
adab_fppdtm <- core_aab %>%
  filter(ADATYPE == "ADA_BAB") %>%
  distinct(
    STUDYID, USUBJID, PARCAT1, BASETYPE, PARAM, ADATYPE, ADAPARM,
    ISTESTCD, ISTEST, ISCAT, ISBDAGNT, AVISITN, AVISIT, ISSPEC, FPPDTM
  ) %>%
  mutate(
    PARAMCD = "FPPDTMy",
    PARAM = paste("First Post Dose Positive Datetime,", PARCAT1, sep = " "),
    AVAL = (as.numeric(FPPDTM)),
    AVALC = if_else(is.na(FPPDTM), NA_character_,
      toupper(as.character(format(FPPDTM, "%d%b%Y:%H:%M:%S")))
    ),
    AVALU = NA_character_
  )

# Get Patient flag LPPDTM
adab_lppdtm <- core_aab %>%
  filter(ADATYPE == "ADA_BAB") %>%
  distinct(
    STUDYID, USUBJID, PARCAT1, BASETYPE, PARAM, ADATYPE, ADAPARM,
    ISTESTCD, ISTEST, ISCAT, ISBDAGNT, AVISITN, AVISIT, ISSPEC, LPPDTM
  ) %>%
  mutate(
    PARAMCD = "LPPDTMy",
    PARAM = paste("Last Post Dose Positive Datetime,", PARCAT1, sep = " "),
    AVAL = as.numeric(LPPDTM),
    AVALC = if_else(is.na(LPPDTM), NA_character_,
      toupper(as.character(format(LPPDTM, "%d%b%Y:%H:%M:%S")))
    ),
    AVALU = NA_character_
  )

# Set all the standard PARAM components together -----------------------------------
adab_paramcds <- bind_rows(
  adab_titer, adab_nabvis, adab_result, adab_nabres, adab_bflag,
  adab_incucd, adab_enhanc, adab_emerpos, adab_trunaff, adab_emerneg, adab_notrrel,
  adab_adastat, adab_timada, adab_persada, adab_tranada, adab_nabstat
)

# In this sample, also have BY VISIT parameters,
adab_visits <- bind_rows(adab_tflagv, adab_pbflagv, adab_adastatv)

# Create the final parameterized dataset --------------------------------------------

# To keep only standard parameters:
# adab_study <- adab_paramcds
# To include BY VISIT parameters and/or others:
adab_study <- bind_rows(adab_paramcds, adab_visits, adab_adadur, adab_fppdtm, adab_lppdtm)

# Drop temporary variables and ADA Flag variables that are now parameterized, further below is a
# final 'select' statement to to customize the final select.
adab_study <- adab_study %>%
  select(
    -TIMADA, -ADADUR, -TRANADA, -PERSADA, -TRANADAE, -PERSADAE, -INDUCED, -ENHANCED,
    -RESULTC, -RESULTN, -ADASTAT, -BFLAG, -TFLAG, -PBFLAG, -FPPDTM, -LPPDTM, -TFLAGV, -PBFLAGV,
    -ADASTATV, -nabstat_opt1, -nabstat_opt2, -NABSTAT, -MAXCHG, -VALIDBASE, -VALIDPOST,
    -tdur, -ADASTAT_MAIN, -NABPOSTMISS, -TRTSDT
  )

# Merge in ADSL Values --------------------------------
adab_adsl <- adab_study %>%
  derive_vars_merged(
    dataset_add = adsl,
    by_vars = exprs(STUDYID, USUBJID)
  )

# Compute COHORT From ADSL, other source could be ADSUB
adab_cohort <- adab_adsl %>%
  derive_vars_merged(
    dataset_add = adsl,
    new_vars = exprs(COHORT = ARMCD),
    by_vars = exprs(STUDYID, USUBJID)
  )

# Compute optional ADAFL
# Method could be ADSL.SAFFL, ADAB.ADPBLPFL by ISTESTCD, etc.
adab_adafl <- adab_cohort %>%
  derive_vars_merged(
    dataset_add = adsl,
    new_vars = exprs(ADAFL = SAFFL),
    by_vars = exprs(STUDYID, USUBJID)
  )

# Study Specific Specs Post-Processing ------------------------------------

# Create a Tibble to map above PARAMCD and PARAM to final study spec values
# Multiple NAB and ADA analytes example usage:
# If ISTESTCD is 'ADA_BAB' and ISBDAGNT is 'XANOMELINE', RESULT1, ADASTAT1, etc. PARAM_SUFFIX = '(1)'
# If ISTESTCD is 'ADA_BAB' and ISBDAGNT is 'Y012345678', RESULT2, ADASTAT2, etc. PARAM_SUFFIX = '(2)'
# If ISTESTCD is 'ADA_NAB' and ISBDAGNT is 'XANOMELINE', RESULT3, etc. PARAM_SUFFIX = '(3)'
# If ISTESTCD is 'ADA_NAB' and ISBDAGNT is 'Y012345678', RESULT4, etc. PARAM_SUFFIX = '(4)'

adab_param_data <- tribble(
  ~PARAMCD, ~ADATYPE, ~ADAPARM, ~PARAMCD_NEW, ~PARAM_SUFFIX,
  "XANOMELINE", "ADA_BAB", "XANOMELINE", "XANOMELINE", "(1)",
  "XANOMELINE", "ADA_NAB", "XANOMELINE", "XANOMELINE", "(1)",
  "RESULTy", "ADA_BAB", "XANOMELINE", "RESULT1", "(1)",
  "RESULTy", "ADA_NAB", "XANOMELINE", "RESULT2", "(2)",
  "BFLAGy", "ADA_BAB", "XANOMELINE", "BFLAG1", "(1)",
  "INDUCDy", "ADA_BAB", "XANOMELINE", "INDUCD1", "(1)",
  "ENHANCy", "ADA_BAB", "XANOMELINE", "ENHANC1", "(1)",
  "EMERPOSy", "ADA_BAB", "XANOMELINE", "EMERPOS1", "(1)",
  "TRUNAFFy", "ADA_BAB", "XANOMELINE", "TRUNAFF1", "(1)",
  "EMERNEGy", "ADA_BAB", "XANOMELINE", "EMERNEG1", "(1)",
  "NOTRRELy", "ADA_BAB", "XANOMELINE", "NOTRREL1", "(1)",
  "ADASTATy", "ADA_BAB", "XANOMELINE", "ADASTAT1", "(1)",
  "TIMADAy", "ADA_BAB", "XANOMELINE", "TIMADA1", "(1)",
  "PERSADAy", "ADA_BAB", "XANOMELINE", "PERSADA1", "(1)",
  "TRANADAy", "ADA_BAB", "XANOMELINE", "TRANADA1", "(1)",
  "NABSTATy", "ADA_NAB", "XANOMELINE", "NABSTAT1", "(1)",
  "ADASTATV", "ADA_BAB", "XANOMELINE", "ADASTTV1", "(1)",
  "TFLAGV", "ADA_BAB", "XANOMELINE", "TFLAGV1", "(1)",
  "PBFLAGV", "ADA_BAB", "XANOMELINE", "PBFLAGV1", "(1)",
  "LPPDTMy", "ADA_BAB", "XANOMELINE", "LPPDTM1", "(1)",
  "FPPDTMy", "ADA_BAB", "XANOMELINE", "FPPDTM1", "(1)",
  "ADADURy", "ADA_BAB", "XANOMELINE", "ADADUR1", "(1)",
)

# Merge the Parameter dataset into the main data
adab_params <- adab_adafl %>%
  derive_vars_merged(
    dataset_add = adab_param_data,
    by_vars = exprs(PARAMCD, ADAPARM, ADATYPE)
  ) %>%
  # for the original data, assign PARAM
  mutate(
    PARAM = case_when(
      !is.na(PARAM_SUFFIX) ~ paste(PARAM, PARAM_SUFFIX, sep = " "),
      TRUE ~ PARAM
    ),
    # Example to assign PARAMCD based on ADAPARM assigned at program start (SDTM.IS < v2.0)
    PARAMCD = case_when(
      PARAMCD == ADAPARM ~ PARAMCD,
      TRUE ~ PARAMCD_NEW
    ),
    # Example to assign PARAMCD based ISTESTCD type and last 5 chars of ISBDAGNT (SDTM.IS >= v2.0)
    PARAMCD = case_when(
      ISTESTCD == "ADA_BAB" & PARAMCD == "XANOMELINE" ~
        paste("BAB", substr(ISBDAGNT, 1, 5), sep = ""),
      ISTESTCD == "ADA_NAB" & PARAMCD == "XANOMELINE" ~
        paste("NAB", substr(ISBDAGNT, 1, 5), sep = ""),
      TRUE ~ PARAMCD
    )
  )

# Sort by the key variables then compute ASEQ
adab_prefinal <- adab_params %>%
  # Calculate ASEQ
  derive_var_obs_number(
    new_var = ASEQ,
    by_vars = exprs(STUDYID, USUBJID),
    order = exprs(PARCAT1, PARAMCD, BASETYPE, NFRLT, AFRLT, ISSEQ),
    check_type = "error"
  )

# Choose final variables to keep
# When SDTM V1.x, suggest remove ISBDAGNT
adab <- adab_prefinal %>%
  select(
    STUDYID, USUBJID, SUBJID, SITEID, ASEQ,
    REGION1, COUNTRY, ETHNIC,
    AGE, AGEU, SEX, RACE,
    SAFFL, TRT01P, TRT01A,
    TRTSDTM, TRTSDT, TRTEDTM, TRTEDT,
    ISSEQ, ISTESTCD, ISTEST, ISCAT, ISBDAGNT,
    ISSTRESC, ISSTRESN, ISSTRESU,
    ISSTAT, ISREASND, ISSPEC,
    DTL, MRT,
    VISITNUM, VISIT, VISITDY,
    EPOCH, ISDTC, ISDY, ISTPT, ISTPTNUM,
    PARAM, PARAMCD, PARCAT1,
    AVAL, AVALC, AVALU,
    BASETYPE, BASE, CHG, DTYPE,
    ADTM, ADT, ADY, ATMF,
    AVISIT, AVISITN, ATPT,
    APHASE, APHASEN, APERIOD, APERIODC,
    FANLDTM, FANLDT, FANLTM, FANLTMF,
    NFRLT, AFRLT, FRLTU,
    ABLFL, ADABLPFL, ADPBLPFL, ADAFL
  )

# Final Steps, Select final variables and Add labels
# This process will be based on your metadata, no example given for this reason
# ...
# ---- Save output ----

# Change to whichever directory you want to save the dataset in
dir <- tools::R_user_dir("admiral_templates_data", which = "cache")
if (!file.exists(dir)) {
  # Create the folder
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
}
save(adab, file = file.path(dir, "adab.rda"), compress = "bzip2")
```

ad_adae.R

```
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

ae <- pharmaversesdtm::ae
suppae <- pharmaversesdtm::suppae
adsl <- admiral::admiral_adsl
ex_single <- admiral::ex_single

# When SAS datasets are imported into R using haven::read_sas(), missing
# character values from SAS appear as "" characters in R, instead of appearing
# as NA values. Further details can be obtained via the following link:
# https://pharmaverse.github.io/admiral/cran-release/articles/admiral.html#handling-of-missing-values # nolint

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
) %>%
  derive_vars_dtm(
    dtc = EXENDTC,
    new_vars_prefix = "EXEN",
    time_imputation = "last",
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
  ## Derive treatment dose and unit ----
  derive_vars_joined(
    dataset_add = ex_ext,
    by_vars = exprs(STUDYID, USUBJID),
    new_vars = exprs(DOSEON = EXDOSE, DOSEU = EXDOSU),
    join_vars = exprs(EXSTDTM, EXENDTM),
    join_type = "all",
    filter_add = (EXDOSE > 0 | (EXDOSE == 0 & grepl("PLACEBO", EXTRT))) & !is.na(EXSTDTM),
    filter_join = EXSTDTM <= ASTDTM & (ASTDTM <= EXENDTM | is.na(EXENDTM))
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
```

ad_adcm.R

```
# Name: ADCM
#
# Label: Concomitant Medications Analysis Dataset
#
# Input: cm, adsl
library(admiral)
library(pharmaversesdtm) # Contains example datasets from the CDISC pilot project
library(dplyr)
library(lubridate)

# Load source datasets ----

# Use e.g. haven::read_sas to read in .sas7bdat, or other suitable functions
# as needed and assign to the variables below.
# For illustration purposes read in admiral test data

cm <- pharmaversesdtm::cm
adsl <- admiral::admiral_adsl

# When SAS datasets are imported into R using haven::read_sas(), missing
# character values from SAS appear as "" characters in R, instead of appearing
# as NA values. Further details can be obtained via the following link:
# https://pharmaverse.github.io/admiral/cran-release/articles/admiral.html#handling-of-missing-values # nolint

cm <- convert_blanks_to_na(cm)

# Derivations ----

# Get list of ADSL vars required for derivations
adsl_vars <- exprs(TRTSDT, TRTEDT, DTHDT, EOSDT, TRT01P, TRT01A)

adcm <- cm %>%
  # Join ADSL with CM (only ADSL vars required for derivations)
  derive_vars_merged(
    dataset_add = adsl,
    new_vars = adsl_vars,
    by = exprs(STUDYID, USUBJID)
  ) %>%
  ## Derive analysis start time ----
  derive_vars_dtm(
    dtc = CMSTDTC,
    new_vars_prefix = "AST",
    highest_imputation = "M",
    min_dates = exprs(TRTSDT)
  ) %>%
  ## Derive analysis end time ----
  derive_vars_dtm(
    dtc = CMENDTC,
    new_vars_prefix = "AEN",
    highest_imputation = "M",
    date_imputation = "last",
    time_imputation = "last",
    max_dates = exprs(DTHDT, EOSDT)
  ) %>%
  ## Derive analysis end/start date -----
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

## Derive flags ----
adcm <- adcm %>%
  # Derive On-Treatment flag
  # Set `span_period = TRUE` if you want occurrences that started prior to drug
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
  # This variable is producer specific and may be used to indicate particular
  # records to be used in subsequent derivations or analysis.
  mutate(ANL01FL = if_else(ONTRTFL == "Y", "Y", NA_character_)) %>%
  # Derive 1st Occurrence of Preferred Term Flag
  restrict_derivation(
    derivation = derive_var_extreme_flag,
    args = params(
      new_var = AOCCPFL,
      by_vars = exprs(USUBJID, CMDECOD),
      order = exprs(ASTDTM, CMSEQ),
      mode = "first"
    ),
    filter = ANL01FL == "Y"
  )


## Derive APHASE and APHASEN Variable ----
# Other timing variable can be derived similarly.
# See also the "Visit and Period Variables" vignette
# (https://pharmaverse.github.io/admiral/cran-release/articles/visits_periods.html)
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
    by_vars = exprs(STUDYID, USUBJID)
  )


# Save output ----

# Change to whichever directory you want to save the dataset in
dir <- tools::R_user_dir("admiral_templates_data", which = "cache")
if (!file.exists(dir)) {
  # Create the folder
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
}
save(adcm, file = file.path(dir, "adcm.rda"), compress = "bzip2")
```

ad_adeg.R

```
# Name: ADEG
#
# Label: Electrocardiogram Analysis Dataset
#
# Description: Based on CDISC Pilot data, create ADEG analysis dataset
#
# Input: adsl, eg
library(admiral)
library(pharmaversesdtm) # Contains example datasets from the CDISC pilot project
library(dplyr)
library(lubridate)
library(stringr)

# Load source datasets ----

# Use e.g. `haven::read_sas()` to read in .sas7bdat, or other suitable functions
# as needed and assign to the variables below.
# For illustration purposes read in admiral test data

eg <- pharmaversesdtm::eg
adsl <- admiral::admiral_adsl

# When SAS datasets are imported into R using haven::read_sas(), missing
# character values from SAS appear as "" characters in R, instead of appearing
# as NA values. Further details can be obtained via the following link:
# https://pharmaverse.github.io/admiral/cran-release/articles/admiral.html#handling-of-missing-values # nolint

eg <- convert_blanks_to_na(eg)

# Lookup tables ----

# Assign PARAMCD, PARAM, and PARAMN
param_lookup <- tibble::tribble(
  ~EGTESTCD, ~PARAMCD,                                                  ~PARAM, ~PARAMN,
  "ECGINT",  "EGINTP",                                    "ECG Interpretation",       1,
  "HR",          "HR",                                "Heart Rate (beats/min)",       2,
  "RR",          "RR",                                      "RR Duration (ms)",       3,
  "RRR",        "RRR",                            "RR Duration Rederived (ms)",       4,
  "QT",          "QT",                                      "QT Duration (ms)",      10,
  "QTCBR",    "QTCBR",     "QTcB - Bazett's Correction Formula Rederived (ms)",      11,
  "QTCFR",    "QTCFR", "QTcF - Fridericia's Correction Formula Rederived (ms)",      12,
  "QTLCR",    "QTLCR",      "QTlc - Sagie's Correction Formula Rederived (ms)",      13,
)

range_lookup <- tibble::tribble(
  ~PARAMCD, ~ANRLO, ~ANRHI,
  "EGINTP",     NA,     NA,
  "HR",         40,    100,
  "RR",        600,   1500,
  "QT",        350,    450,
  "RRR",       600,   1500,
  "QTCBR",     350,    450,
  "QTCFR",     350,    450,
  "QTLCR",     350,    450
)

# Assign AVALCAx
avalcax_lookup <- exprs(
  ~condition,                                                ~AVALCAT1, ~AVALCA1N,
  startsWith(PARAMCD, "QT") & AVAL <= 450,                 "<= 450 ms",         1,
  startsWith(PARAMCD, "QT") & AVAL > 450 & AVAL <= 480, ">450<=480 ms",         2,
  startsWith(PARAMCD, "QT") & AVAL > 480 & AVAL <= 500, ">480<=500 ms",         3,
  startsWith(PARAMCD, "QT") & AVAL > 500,                    ">500 ms",         4
)
# Assign CHGCAx
chgcax_lookup <- exprs(
  ~condition,                                           ~CHGCAT1, ~CHGCAT1N,
  startsWith(PARAMCD, "QT") & CHG <= 30,              "<= 30 ms",         1,
  startsWith(PARAMCD, "QT") & CHG > 30 & CHG <= 60, ">30<=60 ms",         2,
  startsWith(PARAMCD, "QT") & CHG > 60,                 ">60 ms",         3
)

# Derivations ----

# Get list of ADSL vars required for derivations
adsl_vars <- exprs(TRTSDT, TRTEDT, TRT01A, TRT01P)

adeg <- eg %>%
  # Join ADSL & EG (need TRTSDT for ADY derivation)
  derive_vars_merged(
    dataset_add = adsl,
    new_vars = adsl_vars,
    by_vars = exprs(STUDYID, USUBJID)
  ) %>%
  ## Calculate ADTM, ADY ----
  derive_vars_dtm(
    new_vars_prefix = "A",
    dtc = EGDTC,
  ) %>%
  derive_vars_dy(reference_date = TRTSDT, source_vars = exprs(ADTM))

adeg <- adeg %>%
  ## Add PARAMCD only (add PARAM, etc later) ----
  derive_vars_merged_lookup(
    dataset_add = param_lookup,
    new_vars = exprs(PARAMCD),
    by_vars = exprs(EGTESTCD)
  ) %>%
  ## Calculate AVAL and AVALC ----
  mutate(
    AVAL = EGSTRESN,
    AVALC = ifelse(
      is.na(EGSTRESN) | as.character(EGSTRESN) != EGSTRESC,
      EGSTRESC,
      NA
    )
  ) %>%
  ## Derive new parameters based on existing records ----
  # Note that, for the following four `derive_param_*()` functions, only the
  # variables specified in `by_vars` will be populated in the newly created
  # records.

  # Derive RRR
  derive_param_rr(
    by_vars = exprs(STUDYID, USUBJID, !!!adsl_vars, VISIT, VISITNUM, EGTPT, EGTPTNUM, ADTM, ADY),
    set_values_to = exprs(PARAMCD = "RRR"),
    hr_code = "HR",
    get_unit_expr = tolower(EGSTRESU),
    filter = EGSTAT != "NOT DONE" | is.na(EGSTAT)
  ) %>%
  # Derive QTCBR
  derive_param_qtc(
    by_vars = exprs(STUDYID, USUBJID, !!!adsl_vars, VISIT, VISITNUM, EGTPT, EGTPTNUM, ADTM, ADY),
    method = "Bazett",
    set_values_to = exprs(PARAMCD = "QTCBR"),
    qt_code = "QT",
    rr_code = "RR",
    get_unit_expr = EGSTRESU,
    filter = EGSTAT != "NOT DONE" | is.na(EGSTAT)
  ) %>%
  # Derive QTCFR
  derive_param_qtc(
    by_vars = exprs(STUDYID, USUBJID, !!!adsl_vars, VISIT, VISITNUM, EGTPT, EGTPTNUM, ADTM, ADY),
    method = "Fridericia",
    set_values_to = exprs(PARAMCD = "QTCFR"),
    qt_code = "QT",
    rr_code = "RR",
    get_unit_expr = EGSTRESU,
    filter = EGSTAT != "NOT DONE" | is.na(EGSTAT)
  ) %>%
  # Derive QTLCR
  derive_param_qtc(
    by_vars = exprs(STUDYID, USUBJID, !!!adsl_vars, VISIT, VISITNUM, EGTPT, EGTPTNUM, ADTM, ADY),
    method = "Sagie",
    set_values_to = exprs(PARAMCD = "QTLCR"),
    qt_code = "QT",
    rr_code = "RR",
    get_unit_expr = EGSTRESU,
    filter = EGSTAT != "NOT DONE" | is.na(EGSTAT)
  )

## Get visit info ----
# See also the "Visit and Period Variables" vignette
# (https://pharmaverse.github.io/admiral/cran-release/articles/visits_periods.html#visits)
adeg <- adeg %>%
  # Derive Timing
  mutate(
    ADT = date(ADTM),
    ATPTN = EGTPTNUM,
    ATPT = EGTPT,
    AVISIT = case_when(
      str_detect(VISIT, "SCREEN|UNSCHED|RETRIEVAL|AMBUL") ~ NA_character_,
      !is.na(VISIT) ~ str_to_title(VISIT),
      TRUE ~ NA_character_
    ),
    AVISITN = as.numeric(
      case_when(
        AVISIT == "Baseline" ~ "0",
        str_detect(VISIT, "WEEK") ~ str_trim(str_replace(VISIT, "WEEK", "")),
        TRUE ~ NA_character_
      )
    ),
  )

## Derive a summary records representing the mean of the triplicates at each visit ----
# (if least 2 records available) for all parameter except EGINTP
adeg <- adeg %>%
  derive_summary_records(
    dataset_add = adeg,
    by_vars = exprs(STUDYID, USUBJID, !!!adsl_vars, PARAMCD, AVISITN, AVISIT, ADT, ADY),
    filter_add = dplyr::n() >= 2 & PARAMCD != "EGINTP",
    set_values_to = exprs(
      AVAL = mean(AVAL, na.rm = TRUE),
      DTYPE = "AVERAGE"
    )
  )

adeg <- adeg %>%
  ## Calculate ONTRTFL: from trt start up to 30 days after trt ends ----
  derive_var_ontrtfl(
    start_date = ADT,
    ref_start_date = TRTSDT,
    ref_end_date = TRTEDT,
    ref_end_window = 30,
    filter_pre_timepoint = AVISIT == "Baseline"
  )

## Calculate ANRIND: requires the reference ranges ANRLO, ANRHI ----
# Also accommodates the ranges A1LO, A1HI
adeg <- adeg %>%
  derive_vars_merged(
    dataset_add = range_lookup,
    by_vars = exprs(PARAMCD)
  ) %>%
  # Calculate ANRIND
  derive_var_anrind()

## Derive baseline flags ----
adeg <- adeg %>%
  # Calculate BASETYPE
  derive_basetype_records(
    basetypes = exprs(
      "BASELINE DAY 1" = TRUE
    )
  ) %>%
  # Calculate ABLFL
  restrict_derivation(
    derivation = derive_var_extreme_flag,
    args = params(
      by_vars = exprs(STUDYID, USUBJID, BASETYPE, PARAMCD),
      order = exprs(ADT, VISITNUM, EGSEQ),
      new_var = ABLFL,
      mode = "last"
    ),
    filter = ((!is.na(AVAL) | !is.na(AVALC)) &
      ADT <= TRTSDT & !is.na(BASETYPE) & DTYPE == "AVERAGE" &
      PARAMCD != "EGINTP"
    )
  )

## Derive baseline information ----
adeg <- adeg %>%
  # Calculate BASE
  derive_var_base(
    by_vars = exprs(STUDYID, USUBJID, PARAMCD, BASETYPE),
    source_var = AVAL,
    new_var = BASE
  ) %>%
  # Calculate BASEC
  derive_var_base(
    by_vars = exprs(STUDYID, USUBJID, PARAMCD, BASETYPE),
    source_var = AVALC,
    new_var = BASEC
  ) %>%
  # Calculate BNRIND
  derive_var_base(
    by_vars = exprs(STUDYID, USUBJID, PARAMCD, BASETYPE),
    source_var = ANRIND,
    new_var = BNRIND
  ) %>%
  # Calculate CHG for post-baseline records
  # The decision on how to populate pre-baseline and baseline values of CHG is left to producer choice
  restrict_derivation(
    derivation = derive_var_chg,
    filter = AVISITN > 0
  ) %>%
  # Calculate PCHG for post-baseline records
  # The decision on how to populate pre-baseline and baseline values of PCHG is left to producer choice
  restrict_derivation(
    derivation = derive_var_pchg,
    filter = AVISITN > 0
  )

## ANL01FL: Flag last result within an AVISIT and ATPT for post-baseline records ----
adeg <- adeg %>%
  restrict_derivation(
    derivation = derive_var_extreme_flag,
    args = params(
      by_vars = exprs(USUBJID, PARAMCD, AVISIT, ATPT, DTYPE),
      order = exprs(ADT, AVAL),
      new_var = ANL01FL,
      mode = "last"
    ),
    filter = !is.na(AVISITN) & (ONTRTFL == "Y" | ABLFL == "Y") & DTYPE == "AVERAGE"
  )

## Get treatment information ----
# See also the "Visit and Period Variables" vignette
# (https://pharmaverse.github.io/admiral/cran-release/articles/visits_periods.html#treatment_bds)
adeg <- adeg %>%
  # Assign TRTA, TRTP
  mutate(TRTP = TRT01P, TRTA = TRT01A)

## Get ASEQ and AVALCAT1/CHGCAT1 and add PARAM/PARAMN ----
adeg <- adeg %>%
  # Calculate ASEQ
  derive_var_obs_number(
    new_var = ASEQ,
    by_vars = exprs(STUDYID, USUBJID),
    order = exprs(PARAMCD, ADT, AVISITN, VISITNUM, ATPTN, DTYPE),
    check_type = "error"
  ) %>%
  # Derive AVALCA1N and AVALCAT1
  derive_vars_cat(
    definition = avalcax_lookup
  ) %>%
  # Derive CHGCAT1N and CHGCAT1
  derive_vars_cat(
    definition = chgcax_lookup
  ) %>%
  # Derive PARAM and PARAMN
  derive_vars_merged(
    dataset_add = select(param_lookup, -EGTESTCD),
    by_vars = exprs(PARAMCD)
  )

# Add all ADSL variables
adeg <- adeg %>%
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
save(adeg, file = file.path(dir, "adeg.rda"), compress = "bzip2")
```

ad_adex.R

```
# Name: ADEX
#
# Label: Exposure Analysis Dataset
#
# Input: adsl, ex
#

library(admiral)
library(pharmaversesdtm) # Contains example datasets from the CDISC pilot project
library(dplyr)
library(lubridate)
library(stringr)

# Load source datasets ----
# Use e.g. haven::read_sas to read in .sas7bdat, or other suitable functions
#  as needed and assign to the variables below.
# The CDISC pilot datasets are used for demonstration purpose.

ex <- pharmaversesdtm::ex
adsl <- admiral::admiral_adsl

# When SAS datasets are imported into R using haven::read_sas(), missing
# character values from SAS appear as "" characters in R, instead of appearing
# as NA values. Further details can be obtained via the following link:
# https://pharmaverse.github.io/admiral/cran-release/articles/admiral.html#handling-of-missing-values # nolint

ex <- convert_blanks_to_na(ex)

# The CDISC pilot data does not contain EXADJ,nor a SUPPEX dataset
# add a fake EXADJ to demonstrate the derivation for Dose adjustment flag
# add SUPPEX.EXPLDOS to demonstrate the derivation for dose intensity
# The CDISC pilot EX datasets, contains exposure data for daily dosing. Care should be taken when
# the dosing frequency is different.


ex <- ex %>%
  mutate(
    EXADJ = case_when(
      USUBJID == "01-701-1034" & VISIT %in% c("WEEK 2", "WEEK 24") ~ "ADVERSE EVENT",
      USUBJID == "01-701-1148" & VISIT %in% c("WEEK 24") ~ "MEDICATION ERROR",
      TRUE ~ NA_character_
    ),
    EXDOSE = case_when(
      USUBJID == "01-701-1034" & VISIT %in% c("WEEK 2", "WEEK 24") ~ 0,
      USUBJID == "01-701-1148" & VISIT %in% c("WEEK 24") ~ 0,
      TRUE ~ EXDOSE
    )
  ) %>%
  # add SUPPEX.EXPLDOS to test for dose intensity
  mutate(EXPLDOS = if_else(EXTRT == "PLACEBO", 0, 54))


# Derivations ----

# Get list of ADSL vars required for derivations
adsl_vars <- exprs(TRTSDT, TRTSDTM, TRTEDTM)

# Part 1
# Join ADSL with ex and derive required dates, variables
adex0 <- ex %>%
  # Join ADSL with EX (only ADSL vars required for derivations)
  derive_vars_merged(
    dataset_add = adsl,
    new_vars = adsl_vars,
    by_vars = exprs(STUDYID, USUBJID)
  ) %>%
  ## Calculate ASTDTM, AENDTM using `derive_vars_dtm()` ----
  derive_vars_dtm(
    dtc = EXSTDTC,
    highest_imputation = "M",
    new_vars_prefix = "AST"
  ) %>%
  derive_vars_dtm(
    dtc = EXENDTC,
    highest_imputation = "M",
    date_imputation = "last",
    new_vars_prefix = "AEN"
  ) %>%
  ## Calculate ASTDY, AENDY ----
  derive_vars_dy(
    reference_date = TRTSDTM,
    source_vars = exprs(ASTDTM, AENDTM)
  ) %>%
  ## Add EXDUR, the duration of trt for each record ----
  derive_vars_duration(
    new_var = EXDURD,
    start_date = ASTDTM,
    end_date = AENDTM
  ) %>%
  ## Derive analysis end/start date ----
  derive_vars_dtm_to_dt(exprs(ASTDTM, AENDTM)) %>%
  mutate(
    # Compute the cumulative dose
    DOSEO = EXDOSE * EXDURD,
    PDOSEO = EXPLDOS * EXDURD
  )

# Part 2: 1:1 mapping ----

adex <- bind_rows(
  adex0 %>% mutate(PARAMCD = "DURD", AVAL = EXDURD),
  adex0 %>% mutate(PARAMCD = "DOSE", AVAL = DOSEO),
  adex0 %>% mutate(PARAMCD = "PLDOSE", AVAL = PDOSEO),
  adex0 %>% mutate(PARAMCD = "ADJ", AVALC = if_else(!is.na(EXADJ), "Y", NA_character_)),
  adex0 %>% mutate(PARAMCD = "ADJAE", AVALC = if_else(EXADJ == "ADVERSE EVENT", "Y", NA_character_))
) %>%
  mutate(PARCAT1 = "INDIVIDUAL")

# Part 3: Derive summary parameters ----
# Note that, for the functions `derive_param_exposure()`,
# `derive_param_doseint()` and `derive_param_computed()`, only the variables
# specified in `by_vars` will be populated in the newly created records.

adex <- adex %>%
  # Overall exposure
  call_derivation(
    derivation = derive_param_exposure,
    variable_params = list(
      params(
        set_values_to = exprs(
          PARAMCD = "TDOSE",
          PARCAT1 = "OVERALL",
          AVAL = sum(AVAL, na.rm = TRUE)
        ),
        input_code = "DOSE"
      ),
      params(
        set_values_to = exprs(
          PARAMCD = "TPDOSE",
          PARCAT1 = "OVERALL",
          AVAL = sum(AVAL, na.rm = TRUE)
        ),
        input_code = "PLDOSE"
      ),
      params(
        set_values_to = exprs(
          PARAMCD = "TDURD",
          PARCAT1 = "OVERALL",
          AVAL = sum(AVAL, na.rm = TRUE)
        ),
        input_code = "DURD"
      ),
      params(
        set_values_to = exprs(
          PARAMCD = "TADJ",
          PARCAT1 = "OVERALL",
          AVALC = if_else(sum(!is.na(AVALC)) > 0, "Y", NA_character_)
        ),
        input_code = "ADJ"
      ),
      params(
        set_values_to = exprs(
          PARAMCD = "TADJAE",
          PARCAT1 = "OVERALL",
          AVALC = if_else(sum(!is.na(AVALC)) > 0, "Y", NA_character_)
        ),
        input_code = "ADJAE"
      )
    ),
    dataset_add = adex,
    by_vars = exprs(STUDYID, USUBJID, !!!adsl_vars)
  ) %>%
  # W2-W24 exposure
  call_derivation(
    derivation = derive_param_exposure,
    variable_params = list(
      params(
        set_values_to = exprs(
          PARAMCD = "PDOSE",
          PARCAT1 = "WEEK 2-24",
          AVAL = sum(AVAL, na.rm = TRUE)
        ),
        input_code = "DOSE"
      ),
      params(
        set_values_to = exprs(
          PARAMCD = "PPDOSE",
          PARCAT1 = "WEEK 2-24",
          AVAL = sum(AVAL, na.rm = TRUE)
        ),
        input_code = "PLDOSE"
      ),
      params(
        set_values_to = exprs(
          PARAMCD = "PDURD",
          PARCAT1 = "WEEK 2-24",
          AVAL = sum(AVAL, na.rm = TRUE)
        ),
        input_code = "DURD"
      ),
      params(
        set_values_to = exprs(
          PARAMCD = "PADJ",
          PARCAT1 = "WEEK 2-24",
          AVALC = if_else(sum(!is.na(AVALC)) > 0, "Y", NA_character_)
        ),
        input_code = "ADJ"
      ),
      params(
        set_values_to = exprs(
          PARAMCD = "PADJAE",
          PARCAT1 = "WEEK 2-24",
          AVALC = if_else(sum(!is.na(AVALC)) > 0, "Y", NA_character_)
        ),
        input_code = "ADJAE"
      )
    ),
    dataset_add = adex,
    filter_add = VISIT %in% c("WEEK 2", "WEEK 24"),
    by_vars = exprs(STUDYID, USUBJID, !!!adsl_vars)
  ) %>%
  # Overall Dose intensity and W2-24 dose intensity
  call_derivation(
    derivation = derive_param_doseint,
    variable_params = list(
      params(
        set_values_to = exprs(PARAMCD = "TDOSINT"),
        tadm_code = "TDOSE",
        tpadm_code = "TPDOSE"
      ),
      params(
        set_values_to = exprs(PARAMCD = "PDOSINT"),
        tadm_code = "PDOSE",
        tpadm_code = "PPDOSE"
      )
    ),
    by_vars = exprs(
      STUDYID, USUBJID, !!!adsl_vars, PARCAT1, ASTDTM, ASTDT, AENDTM, AENDT
    )
  ) %>%
  # Overall/W2-24 Average daily dose
  call_derivation(
    derivation = derive_param_computed,
    variable_params = list(
      params(
        parameters = c("TDOSE", "TDURD"),
        set_values_to = exprs(
          AVAL = (AVAL.TDOSE / AVAL.TDURD),
          PARAMCD = "AVDDSE"
        )
      ),
      params(
        parameters = c("PDOSE", "PDURD"),
        set_values_to = exprs(
          AVAL = (AVAL.PDOSE / AVAL.PDURD),
          PARAMCD = "PAVDDSE"
        )
      )
    ),
    by_vars = exprs(
      STUDYID, USUBJID, !!!adsl_vars, PARCAT1, ASTDTM, ASTDT, AENDTM, AENDT
    )
  )

# Part 4: Derive/Assign the last required variables ----

# Assign PARAMCD, PARAM, and PARAMN
# ---- Lookup tables ----
param_lookup <- tibble::tribble(
  ~PARAMCD,                                                     ~PARAM, ~PARAMN,
  "DURD", "Study drug duration during constant dosing interval (days)",       1,
  "DOSE",     "Dose administered during constant dosing interval (mg)",       2,
  "PLDOSE",        "Planned dose during constant dosing interval (mg)",       3,
  "ADJ",               "Dose adjusted during constant dosing interval",       4,
  "ADJAE",  "Dose adjusted  due to AE during constant dosing interval",       5,
  "TDURD",                                   "Overall duration (days)",       7,
  "TDOSE",                              "Total dose administered (mg)",       8,
  "AVDDSE",                  "Average daily dose administered (mg/mg)",      10,
  "TPDOSE",                                  "Total planned dose (mg)",      11,
  "TADJ",                                 "Dose adjusted during study",      13,
  "TADJAE",                     "Dose adjusted during study due to AE",      14,
  "PDURD",                         "Overall duration in W2-W24 (days)",      19,
  "PDOSE",                    "Total dose administered in W2-W2 (mg)4",      20,
  "PPDOSE",                        "Total planned dose in W2-W24 (mg)",      21,
  "PAVDDSE",          "Average daily dose administered in W2-W24 (mg)",      23,
  "PADJ",                                "Dose adjusted during W2-W24",      24,
  "PADJAE",                       "Dose adjusted  in W2-W24 due to AE",      25,
  "TDOSINT",                              "Overall dose intensity (%)",      90,
  "PDOSINT",                                "W2-24 dose intensity (%)",      91
)

# Assign AVALCATx
avalcax_lookup <- exprs(
  ~PARAMCD,            ~condition,             ~AVALCAT1,
  "TDURD",             AVAL >= 90,          ">= 90 days",
  "TDURD", AVAL >= 30 & AVAL < 90, ">= 30 and < 90 days",
  "TDURD",              AVAL < 30,           "< 30 days",
  "PDURD",             AVAL >= 90,          ">= 90 days",
  "PDURD", AVAL >= 30 & AVAL < 90, ">= 30 and < 90 days",
  "PDURD",              AVAL < 30,           "< 30 days",
  "TDOSE",             AVAL < 100,            "< 100 mg",
  "TDOSE",            AVAL >= 100,           ">= 100 mg",
  "PDOSE",             AVAL < 100,            "< 100 mg",
  "PDOSE",            AVAL >= 100,           ">= 100 mg"
)

adex <- adex %>%
  # Add PARAMN and PARAM, AVALU
  derive_vars_merged(
    dataset_add = param_lookup,
    by_vars = exprs(PARAMCD)
  ) %>%
  # Derive AVALCATx
  derive_vars_cat(
    definition = avalcax_lookup,
    by_vars = exprs(PARAMCD)
  ) %>%
  # Calculate ASEQ
  derive_var_obs_number(
    new_var = ASEQ,
    by_vars = exprs(STUDYID, USUBJID),
    order = exprs(PARCAT1, ASTDT, VISIT, VISITNUM, EXSEQ, PARAMN),
    check_type = "error"
  )

# Join all ADSL with EX
adex <- adex %>%
  derive_vars_merged(
    dataset_add = select(adsl, !!!negate_vars(adsl_vars)),
    by_vars = exprs(STUDYID, USUBJID)
  )

# Final Steps, Select final variables and Add labels
# This process will be based on your metadata, no example given for this reason
# ...

# Save output ----

# Change to whichever directory you want to save the dataset in
dir <- tools::R_user_dir("admiral_templates_data", which = "cache")
if (!file.exists(dir)) {
  # Create the folder
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
}
save(adex, file = file.path(dir, "adex.rda"), compress = "bzip2")
```

ad_adlb.R

```
# Name: ADLB
#
# Label: Lab Analysis Dataset
#
# Input: adsl, lb
library(admiral)
library(pharmaversesdtm) # Contains example datasets from the CDISC pilot project
library(dplyr)
library(lubridate)
library(stringr)

# Load source datasets ----

# Use e.g. haven::read_sas to read in .sas7bdat, or other suitable functions
# as needed and assign to the variables below.
# For illustration purposes read in admiral test data

lb <- pharmaversesdtm::lb
adsl <- admiral::admiral_adsl

# When SAS datasets are imported into R using haven::read_sas(), missing
# character values from SAS appear as "" characters in R, instead of appearing
# as NA values. Further details can be obtained via the following link:
# https://pharmaverse.github.io/admiral/cran-release/articles/admiral.html#handling-of-missing-values # nolint

lb <- convert_blanks_to_na(lb)

# Look-up tables ----

# Assign PARAMCD, PARAM, and PARAMN
param_lookup <- tibble::tribble(
  ~LBTESTCD, ~PARAMCD,                                             ~PARAM, ~PARAMN,
  "ALB",        "ALB",                                    "Albumin (g/L)",       1,
  "ALP",      "ALKPH",                       "Alkaline Phosphatase (U/L)",       2,
  "ALT",        "ALT",                   "Alanine Aminotransferase (U/L)",       3,
  "ANISO",    "ANISO",                                       "Anisocytes",       4,
  "AST",        "AST",                 "Aspartate Aminotransferase (U/L)",       5,
  "BASO",      "BASO",                           "Basophils Abs (10^9/L)",       6,
  "BASOLE",  "BASOLE",                  "Basophils/Leukocytes (FRACTION)",       7,
  "BILI",      "BILI",                               "Bilirubin (umol/L)",       8,
  "BUN",        "BUN",                     "Blood Urea Nitrogen (mmol/L)",       9,
  "CA",          "CA",                                 "Calcium (mmol/L)",      10,
  "CHOL",    "CHOLES",                             "Cholesterol (mmol/L)",      11,
  "CK",          "CK",                          "Creatinine Kinase (U/L)",      12,
  "CL",          "CL",                                "Chloride (mmol/L)",      13,
  "COLOR",    "COLOR",                                            "Color",      14,
  "CREAT",    "CREAT",                              "Creatinine (umol/L)",      15,
  "EOS",        "EOS",                             "Eosinophils (10^9/L)",      16,
  "EOSLE",    "EOSLE",                "Eosinophils/Leukocytes (FRACTION)",      17,
  "GGT",        "GGT",                 "Gamma Glutamyl Transferase (U/L)",      18,
  "GLUC",      "GLUC",                                 "Glucose (mmol/L)",      19,
  "HBA1C",    "HBA1C",                               "Hemoglobin A1C (1)",      20,
  "HCT",        "HCT",                                   "Hematocrit (1)",      21,
  "HGB",        "HGB",                              "Hemoglobin (mmol/L)",      22,
  "K",        "POTAS",                               "Potassium (mmol/L)",      23,
  "KETONES",  "KETON",                                          "Ketones",      24,
  "LYM",      "LYMPH",                         "Lymphocytes Abs (10^9/L)",      25,
  "LYMLE",  "LYMPHLE",                "Lymphocytes/Leukocytes (FRACTION)",      26,
  "MACROCY", "MACROC",                                       "Macrocytes",      27,
  "MCH",        "MCH",      "Ery. Mean Corpuscular Hemoglobin (fmol(Fe))",      28,
  "MCHC",      "MCHC", "Ery. Mean Corpuscular HGB Concentration (mmol/L)",      29,
  "MCV",        "MCV",               "Ery. Mean Corpuscular Volume (f/L)",      30,
  "MICROCY", "MICROC",                                       "Microcytes",      31,
  "MONO",      "MONO",                               "Monocytes (10^9/L)",      32,
  "MONOLE",  "MONOLE",                  "Monocytes/Leukocytes (FRACTION)",      33,
  "PH",          "PH",                                               "pH",      34,
  "PHOS",      "PHOS",                               "Phosphate (mmol/L)",      35,
  "PLAT",      "PLAT",                                "Platelet (10^9/L)",      36,
  "POIKILO", "POIKIL",                                     "Poikilocytes",      37,
  "POLYCHR", "POLYCH",                                    "Polychromasia",      38,
  "PROT",      "PROT",                                    "Protein (g/L)",      39,
  "RBC",        "RBC",                              "Erythrocytes (TI/L)",      40,
  "SODIUM",  "SODIUM",                                  "Sodium (mmol/L)",      41,
  "SPGRAV",  "SPGRAV",                                 "Specific Gravity",      42,
  "TSH",        "TSH",                               "Thyrotropin (mU/L)",      43,
  "URATE",    "URATE",                                   "Urate (umol/L)",      44,
  "UROBIL",  "UROBIL",                                     "Urobilinogen",      45,
  "VITB12",  "VITB12",                             "Vitamin B12 (pmol/L)",      46,
  "WBC",        "WBC",                              "Leukocytes (10^9/L)",      47
)


# Derivations ----

# Get list of ADSL vars required for derivations
adsl_vars <- exprs(TRTSDT, TRTEDT, TRT01A, TRT01P)

adlb <- lb %>%
  # Join ADSL with LB (need TRTSDT for ADY derivation)
  derive_vars_merged(
    dataset_add = adsl,
    new_vars = adsl_vars,
    by_vars = exprs(STUDYID, USUBJID)
  ) %>%
  ## Calculate ADT, ADY ----
  derive_vars_dt(
    new_vars_prefix = "A",
    dtc = LBDTC
  ) %>%
  derive_vars_dy(reference_date = TRTSDT, source_vars = exprs(ADT))

adlb <- adlb %>%
  ## Add PARAMCD PARAM and PARAMN - from LOOK-UP table ----
  # Replace with PARAMCD lookup function
  derive_vars_merged_lookup(
    dataset_add = param_lookup,
    new_vars = exprs(PARAMCD, PARAM, PARAMN),
    by_vars = exprs(LBTESTCD),
    check_type = "none",
    print_not_mapped = FALSE
  ) %>%
  ## Calculate PARCAT1 AVAL AVALC ANRLO ANRHI ----
  mutate(
    PARCAT1 = LBCAT,
    AVAL = LBSTRESN,
    # Only populate AVALC if character value is non-redundant with AVAL
    AVALC = ifelse(
      is.na(LBSTRESN) | as.character(LBSTRESN) != LBSTRESC,
      LBSTRESC,
      NA
    ),
    ANRLO = LBSTNRLO,
    ANRHI = LBSTNRHI
  )

# Derive Absolute values from fractional Differentials using WBC
# Only derive where absolute values do not already exist
# Need to populate ANRLO and ANRHI for newly created records
adlb <- adlb %>%
  # Derive absolute Basophils
  derive_param_wbc_abs(
    by_vars = exprs(STUDYID, USUBJID, !!!adsl_vars, DOMAIN, VISIT, VISITNUM, ADT, ADY),
    set_values_to = exprs(
      PARAMCD = "BASO",
      PARAM = "Basophils Abs (10^9/L)",
      PARAMN = 6,
      DTYPE = "CALCULATION",
      PARCAT1 = "HEMATOLOGY"
    ),
    get_unit_expr = extract_unit(PARAM),
    diff_code = "BASOLE"
  ) %>%
  # Derive absolute Lymphocytes
  derive_param_wbc_abs(
    by_vars = exprs(STUDYID, USUBJID, !!!adsl_vars, DOMAIN, VISIT, VISITNUM, ADT, ADY),
    set_values_to = exprs(
      PARAMCD = "LYMPH",
      PARAM = "Lymphocytes Abs (10^9/L)",
      PARAMN = 25,
      DTYPE = "CALCULATION",
      PARCAT1 = "HEMATOLOGY"
    ),
    get_unit_expr = extract_unit(PARAM),
    diff_code = "LYMPHLE"
  )

## Get Visit Info ----
# See also the "Visit and Period Variables" vignette
# (https://pharmaverse.github.io/admiral/cran-release/articles/visits_periods.html#visits)
adlb <- adlb %>%
  # Derive Timing
  mutate(
    AVISIT = case_when(
      str_detect(VISIT, "SCREEN") ~ "Baseline",
      !is.na(VISIT) ~ str_to_title(VISIT),
      TRUE ~ NA_character_
    ),
    AVISITN = case_when(
      AVISIT == "Baseline" ~ 0,
      !is.na(VISITNUM) ~ VISITNUM
    )
  )

adlb <- adlb %>%
  ## Calculate ONTRTFL ----
  derive_var_ontrtfl(
    start_date = ADT,
    ref_start_date = TRTSDT,
    ref_end_date = TRTEDT,
    filter_pre_timepoint = AVISIT == "Baseline"
  )

## Calculate ANRIND : requires the reference ranges ANRLO, ANRHI ----
adlb <- adlb %>%
  derive_var_anrind()

## Derive baseline flags ----
adlb <- adlb %>%
  # Calculate BASETYPE
  mutate(
    BASETYPE = "LAST"
  ) %>%
  # Calculate ABLFL
  restrict_derivation(
    derivation = derive_var_extreme_flag,
    args = params(
      by_vars = exprs(STUDYID, USUBJID, BASETYPE, PARAMCD),
      order = exprs(ADT, VISITNUM, LBSEQ),
      new_var = ABLFL,
      mode = "last"
    ),
    filter = (!is.na(AVAL) & ADT <= TRTSDT & !is.na(BASETYPE))
  )

## Derive baseline information ----
adlb <- adlb %>%
  # Calculate BASE
  derive_var_base(
    by_vars = exprs(STUDYID, USUBJID, PARAMCD, BASETYPE),
    source_var = AVAL,
    new_var = BASE
  ) %>%
  # Calculate BASEC
  derive_var_base(
    by_vars = exprs(STUDYID, USUBJID, PARAMCD, BASETYPE),
    source_var = AVALC,
    new_var = BASEC
  ) %>%
  # Calculate BNRIND
  derive_var_base(
    by_vars = exprs(STUDYID, USUBJID, PARAMCD, BASETYPE),
    source_var = ANRIND,
    new_var = BNRIND
  ) %>%
  # Calculate CHG for post-baseline records
  # The decision on how to populate pre-baseline and baseline values of CHG is left to producer choice
  restrict_derivation(
    derivation = derive_var_chg,
    filter = AVISITN > 0
  ) %>%
  # Calculate PCHG for post-baseline records
  # The decision on how to populate pre-baseline and baseline values of PCHG is left to producer choice
  restrict_derivation(
    derivation = derive_var_pchg,
    filter = AVISITN > 0
  )

## Calculate lab grading ----

# Assign ATOXDSCL and ATOXDSCH to hold lab grading terms
# ATOXDSCL and ATOXDSCH hold terms defined by NCI-CTCAEv4.
# See (https://pharmaverse.github.io/admiral/cran-release/articles/lab_grading.html#implement_ctcv4)
grade_lookup <- tibble::tribble(
  ~PARAMCD,                  ~ATOXDSCL,                              ~ATOXDSCH,
  "ALB",             "Hypoalbuminemia",                          NA_character_,
  "ALKPH",               NA_character_,       "Alkaline phosphatase increased",
  "ALT",                 NA_character_,   "Alanine aminotransferase increased",
  "AST",                 NA_character_, "Aspartate aminotransferase increased",
  "BILI",                NA_character_,            "Blood bilirubin increased",
  "CA",                 "Hypocalcemia",                        "Hypercalcemia",
  "CHOLES",              NA_character_,                     "Cholesterol high",
  "CK",                  NA_character_,                        "CPK increased",
  "CREAT",               NA_character_,                 "Creatinine increased",
  "GGT",                 NA_character_,                        "GGT increased",
  "GLUC",               "Hypoglycemia",                        "Hyperglycemia",
  "HGB",                      "Anemia",                 "Hemoglobin increased",
  "POTAS",               "Hypokalemia",                         "Hyperkalemia",
  "LYMPH", "CD4 lymphocytes decreased",                          NA_character_,
  "PHOS",           "Hypophosphatemia",                          NA_character_,
  "PLAT",   "Platelet count decreased",                          NA_character_,
  "SODIUM",             "Hyponatremia",                        "Hypernatremia",
  "WBC",  "White blood cell decreased",                         "Leukocytosis",
)

# Assign grade criteria
# metadata atoxgr_criteria_ctcv4 used to implement NCI-CTCAEv4
# user could change to atoxgr_criteria_ctcv5 to implement NCI-CTCAEv5
# Note: Hyperglycemia and Hypophosphatemia not defined in NCI-CTCAEv5 so
# user would need to amend look-up table grade_lookup
# See (https://pharmaverse.github.io/admiral/cran-release/articles/lab_grading.html#implement_ctcv5)
grade_crit <- atoxgr_criteria_ctcv4


# Add ATOXDSCL and ATOXDSCH
adlb <- adlb %>%
  derive_vars_merged(
    dataset_add = grade_lookup,
    by_vars = exprs(PARAMCD)
  ) %>%
  # if implementing NCI-CTCAEv5 or NCI-CTCAEv6 then arguments `high_indicator` and `low_indicator`
  # should be assigned in all calls to function `derive_var_atoxgr_dir` below
  # Derive toxicity grade for low values ATOXGRL
  derive_var_atoxgr_dir(
    meta_criteria = grade_crit,
    new_var = ATOXGRL,
    tox_description_var = ATOXDSCL,
    criteria_direction = "L",
    get_unit_expr = extract_unit(PARAM)
  ) %>%
  # Derive toxicity grade for high values ATOXGRH
  # default metadata atoxgr_criteria_ctcv4 used
  derive_var_atoxgr_dir(
    meta_criteria = grade_crit,
    new_var = ATOXGRH,
    tox_description_var = ATOXDSCH,
    criteria_direction = "H",
    get_unit_expr = extract_unit(PARAM)
  ) %>%
  # (Optional) derive overall grade ATOXGR (combining ATOXGRL and ATOXGRH)
  derive_var_atoxgr() %>%
  # Derive baseline toxicity grade for low values BTOXGRL
  derive_var_base(
    by_vars = exprs(STUDYID, USUBJID, PARAMCD, BASETYPE),
    source_var = ATOXGRL,
    new_var = BTOXGRL
  ) %>%
  # Derive baseline toxicity grade for high values BTOXGRH
  derive_var_base(
    by_vars = exprs(STUDYID, USUBJID, PARAMCD, BASETYPE),
    source_var = ATOXGRH,
    new_var = BTOXGRH
  ) %>%
  # Derive baseline toxicity grade for for overall grade BTOXGR
  derive_var_base(
    by_vars = exprs(STUDYID, USUBJID, PARAMCD, BASETYPE),
    source_var = ATOXGR,
    new_var = BTOXGR
  )


## Calculate R2BASE, R2ANRLO and R2ANRHI ----
adlb <- adlb %>%
  derive_var_analysis_ratio(
    numer_var = AVAL,
    denom_var = BASE
  ) %>%
  derive_var_analysis_ratio(
    numer_var = AVAL,
    denom_var = ANRLO
  ) %>%
  derive_var_analysis_ratio(
    numer_var = AVAL,
    denom_var = ANRHI
  )

## SHIFT derivation ----
adlb <- adlb %>%
  # Derive shift from baseline for analysis indicator
  derive_var_shift(
    new_var = SHIFT1,
    from_var = BNRIND,
    to_var = ANRIND
  ) %>%
  # Derive shift from baseline for overall grade
  restrict_derivation(
    derivation = derive_var_shift,
    args = params(
      new_var = SHIFT2,
      from_var = BTOXGR,
      to_var = ATOXGR
    ),
    filter = !is.na(ATOXDSCL) | !is.na(ATOXDSCH)
  )

## Flag variables (ANL01FL, LVOTFL) ----
# ANL01FL: Flag last result within an AVISIT for post-baseline records
# LVOTFL: Flag last valid on-treatment record
adlb <- adlb %>%
  restrict_derivation(
    derivation = derive_var_extreme_flag,
    args = params(
      by_vars = exprs(USUBJID, PARAMCD, AVISIT),
      order = exprs(ADT, AVAL),
      new_var = ANL01FL,
      mode = "last"
    ),
    filter = !is.na(AVISITN) & ONTRTFL == "Y"
  ) %>%
  restrict_derivation(
    derivation = derive_var_extreme_flag,
    args = params(
      by_vars = exprs(USUBJID, PARAMCD),
      order = exprs(ADT, AVAL),
      new_var = LVOTFL,
      mode = "last"
    ),
    filter = ONTRTFL == "Y"
  )

## Get treatment information ----
# See also the "Visit and Period Variables" vignette
# (https://pharmaverse.github.io/admiral/cran-release/articles/visits_periods.html#treatment_bds)
adlb <- adlb %>%
  # Assign TRTA, TRTP
  mutate(
    TRTP = TRT01P,
    TRTA = TRT01A
  )

## Get extreme values ----
adlb <- adlb %>%
  # get MINIMUM value
  derive_extreme_records(
    dataset_add = adlb,
    by_vars = exprs(STUDYID, USUBJID, PARAMCD, BASETYPE),
    order = exprs(AVAL, ADT, AVISITN),
    mode = "first",
    filter_add = (!is.na(AVAL) & ONTRTFL == "Y"),
    set_values_to = exprs(
      AVISITN = 9997,
      AVISIT = "POST-BASELINE MINIMUM",
      DTYPE = "MINIMUM"
    )
  ) %>%
  # get MAXIMUM value
  derive_extreme_records(
    dataset_add = adlb,
    by_vars = exprs(STUDYID, USUBJID, PARAMCD, BASETYPE),
    order = exprs(desc(AVAL), ADT, AVISITN),
    mode = "first",
    filter_add = (!is.na(AVAL) & ONTRTFL == "Y"),
    set_values_to = exprs(
      AVISITN = 9998,
      AVISIT = "POST-BASELINE MAXIMUM",
      DTYPE = "MAXIMUM"
    )
  ) %>%
  # get LOV value
  derive_extreme_records(
    dataset_add = adlb,
    by_vars = exprs(STUDYID, USUBJID, PARAMCD, BASETYPE),
    order = exprs(ADT, AVISITN),
    mode = "last",
    filter_add = (ONTRTFL == "Y"),
    set_values_to = exprs(
      AVISITN = 9999,
      AVISIT = "POST-BASELINE LAST",
      DTYPE = "LOV"
    )
  )

## Get ASEQ ----
adlb <- adlb %>%
  # Calculate ASEQ
  derive_var_obs_number(
    new_var = ASEQ,
    by_vars = exprs(STUDYID, USUBJID),
    order = exprs(PARAMCD, ADT, AVISITN, VISITNUM),
    check_type = "error"
  )

# Add all ADSL variables
adlb <- adlb %>%
  derive_vars_merged(
    dataset_add = select(adsl, !!!negate_vars(adsl_vars)),
    by_vars = exprs(STUDYID, USUBJID)
  )

# Final Steps, Select final variables and Add labels
# This process will be based on your metadata, no example given for this reason
# ...

# Save output ----

# Change to whichever directory you want to save the dataset in
dir <- tools::R_user_dir("admiral_templates_data", which = "cache")
if (!file.exists(dir)) {
  # Create the folder
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
}
save(adlb, file = file.path(dir, "adlb.rda"), compress = "bzip2")
```

ad_adlbhy.R

```
# Name: ADLBHY
#
# Label: Lab Analysis Dataset for Hy's Law
#
# Input: adlb
library(admiral)
library(dplyr)
library(lubridate)

# Load source datasets ----

# Use e.g. haven::read_sas to read in .sas7bdat, or other suitable functions
# as needed and assign to the variables below.
# For illustration purposes read in admiral test data
# Using use_ad_template("adlb") and assigning the end object as admiral_adlb
# ADLBHY is a special dataset specifically used to check for potential drug induced liver injuries
# Please see "Hy's Law Implementation Guide" on the admiral website for additional information

adlb <- admiral::admiral_adlb

adlb_annotated <- adlb %>%
  filter(PARAMCD %in% c("AST", "ALT", "BILI") & is.na(DTYPE)) %>%
  # Assign flags for contribution to potential Hy's Law event
  slice_derivation(
    derive_vars_crit_flag,
    args = params(
      values_yn = TRUE
    ),
    derivation_slice(
      filter = PARAMCD %in% c("AST", "ALT"),
      args = params(
        condition = AVAL / ANRHI >= 3,
        description = paste(PARAMCD, ">=3xULN")
      )
    ),
    derivation_slice(
      filter = PARAMCD == "BILI",
      args = params(
        condition = AVAL / ANRHI >= 2,
        description = "BILI >= 2xULN"
      )
    )
  ) %>%
  select(STUDYID, USUBJID, TRT01A, PARAMCD, PARAM, LBSEQ, ADT, AVISIT, ADY, AVAL, ANRHI, CRIT1, CRIT1FL)

# Subset Datasets
altast_records <- adlb_annotated %>%
  filter(PARAMCD %in% c("AST", "ALT"))

bili_records <- adlb_annotated %>%
  filter(PARAMCD %in% c("BILI"))

# Flag elevated AST/ALT values with a BILI elevation within 14 days
hylaw_records <- derive_vars_joined(
  dataset = altast_records,
  dataset_add = bili_records,
  by_vars = exprs(STUDYID, USUBJID),
  order = exprs(ADY),
  join_type = "all",
  filter_join = 0 <= ADT.join - ADT & ADT.join - ADT <= 14 & CRIT1FL == "Y" & CRIT1FL.join == "Y",
  new_vars = exprs(BILI_LBSEQ = LBSEQ, BILI_DT = ADT, BILI_CRITFL = CRIT1FL),
  mode = "first"
)

hylaw_records_pts_visits <- hylaw_records %>%
  select(STUDYID, USUBJID, TRT01A) %>% # add AVISIT, ADT for by visit
  distinct()

hylaw_records_fls <- hylaw_records %>%
  select(STUDYID, USUBJID, TRT01A, CRIT1FL, BILI_CRITFL) %>% # add AVISIT, ADT for by visit
  distinct()

# Create new parameters based on records that present potential case
hylaw_params <- derive_param_exist_flag(
  dataset_ref = hylaw_records_pts_visits,
  dataset_add = hylaw_records_fls,
  condition = CRIT1FL == "Y" & BILI_CRITFL == "Y",
  false_value = "N",
  missing_value = "N",
  by_vars = exprs(STUDYID, USUBJID, TRT01A), # add AVISIT, ADT for by visit
  set_values_to = exprs(
    PARAMCD = "HYSLAW",
    PARAM = "ALT/AST >= 3xULN and BILI >= 2xULN"
  )
)

# Row bind back to relevant adlb-like dataset
adlbhy <- adlb_annotated %>%
  bind_rows(hylaw_params)

# Change to whichever directory you want to save the dataset in
dir <- tools::R_user_dir("admiral_templates_data", which = "cache")
if (!file.exists(dir)) {
  # Create the folder
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
}
save(adlbhy, file = file.path(dir, "adlbhy.rda"), compress = "bzip2")
```

ad_admh.R

```
# Name: ADMH
#
# Label: Medical History Analysis Dataset
#
# Input: mh, adsl
library(admiral)
library(pharmaversesdtm) # Contains example datasets from the CDISC pilot project
library(dplyr)
library(lubridate)

# Load source datasets ----

# Use e.g. haven::read_sas to read in .sas7bdat, or other suitable functions
# as needed and assign to the variables below.
# For illustration purposes read in admiral test data

mh <- pharmaversesdtm::mh
queries_mh <- admiral::queries_mh
adsl <- admiral::admiral_adsl

# When SAS datasets are imported into R using haven::read_sas(), missing
# character values from SAS appear as "" characters in R, instead of appearing
# as NA values. Further details can be obtained via the following link:
# https://pharmaverse.github.io/admiral/cran-release/articles/admiral.html#handling-of-missing-values # nolint

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
adsl_vars <- exprs(TRTSDT, TRTEDT, TRT01A, TRT01P, DTHDT, EOSDT)

admh <- mh %>%
  # join ADSL with MH
  derive_vars_merged(
    dataset_add = adsl,
    new_vars = adsl_vars,
    by_vars = exprs(STUDYID, USUBJID)
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
    max_dates = exprs(DTHDT, EOSDT)
  ) %>%
  # Derive analysis start relative day and analysis end relative day
  derive_vars_dy(
    reference_date = TRTSDT,
    source_vars = exprs(ASTDT, AENDT)
  ) %>%
  # Derive analysis date of medical history collection - ADT (company specific variable derivation)
  derive_vars_dt(
    dtc = MHDTC,
    new_vars_prefix = "A"
  ) %>%
  # Derive analysis relative day - ADY (company specific variable derivation)
  derive_vars_dy(
    reference_date = TRTSDT,
    source_vars = exprs(ADT)
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
    by_vars = exprs(USUBJID),
    order = exprs(ASTDT, MHSEQ),
    new_var = AOCCFL,
    mode = "first"
  ) %>%
  derive_var_extreme_flag(
    by_vars = exprs(USUBJID, MHBODSYS),
    order = exprs(USUBJID, MHBODSYS, MHCAT, MHDECOD, MHTERM, ASTDT, MHSEQ),
    new_var = AOCCSFL,
    mode = "first"
  ) %>%
  derive_var_extreme_flag(
    by_vars = exprs(USUBJID, MHDECOD),
    order = exprs(USUBJID, MHBODSYS, MHCAT, MHDECOD, MHTERM, ASTDT, MHSEQ),
    new_var = AOCCPFL,
    mode = "first"
  ) %>%
  # (company specific occurrence flag variables derivation)
  derive_var_extreme_flag(
    by_vars = exprs(USUBJID),
    order = exprs(USUBJID, AHIST, ASTDT, MHSEQ),
    new_var = AOCPFL,
    mode = "first"
  ) %>%
  derive_var_extreme_flag(
    by_vars = exprs(USUBJID, MHBODSYS),
    order = exprs(USUBJID, AHIST, MHBODSYS, MHCAT, ASTDT, MHSEQ),
    new_var = AOCPSFL,
    mode = "first"
  ) %>%
  derive_var_extreme_flag(
    by_vars = exprs(USUBJID, MHDECOD),
    order = exprs(USUBJID, AHIST, MHBODSYS, MHCAT, MHDECOD, MHTERM, ASTDT, MHSEQ),
    new_var = AOCPPFL,
    mode = "first"
  ) %>%
  ## Derive analysis flag (company specific variable derivation) ----
  mutate(ANL01FL = ifelse(MHOCCUR != "N", "Y", NA_character_)) %>%
  ## Assign TRTA, TRTP (company specific variables derivation) ----
  # See also the "Visit and Period Variables" vignette
  # (https://pharmaverse.github.io/admiral/cran-release/articles/visits_periods.html#treatment_bds)
  mutate(
    TRTP = TRT01P,
    TRTA = TRT01A
  ) %>%
  ## Assign APHASE and APHASEN Variable (company specific variable derivation) ----
  # See also the "Visit and Period Variables" vignette
  # (https://pharmaverse.github.io/admiral/cran-release/articles/visits_periods.html#periods_bds)
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
    by_vars = exprs(MHTERM),
    new_vars = exprs(MHTERMN)
  ),
  filter = (MHPRESP == "Y")
)

## Add all ADSL variables ----
admh <- admh %>%
  derive_vars_merged(
    dataset_add = select(adsl, !!!negate_vars(adsl_vars)),
    by_vars = exprs(STUDYID, USUBJID)
  )


# Final Steps, Select final variables and Add labels
# This process will be based on your metadata, no example given for this reason
# ...

# Save output ----

# Change to whichever directory you want to save the dataset in
dir <- tools::R_user_dir("admiral_templates_data", which = "cache")
if (!file.exists(dir)) {
  # Create the folder
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
}
save(admh, file = file.path(dir, "admh.rda"), compress = "bzip2")
```

ad_adpc.R

```
# Name: ADPC
#
# Label: Pharmacokinetics Concentrations Analysis Dataset
#
# Description: Based on simulated data, create ADPC analysis dataset
#   The dataset format is also suitable for Non-compartmental analysis (ADNCA)
#
# Input: pc, ex, vs, adsl
library(admiral)
library(dplyr)
library(lubridate)
library(stringr)

library(pharmaversesdtm) # Contains example datasets from the CDISC pilot project or simulated

# ---- Load source datasets ----

# Use e.g. haven::read_sas to read in .sas7bdat, or other suitable functions
# as needed and assign to the variables below.
# For illustration purposes read in admiral test data

# Load PC, EX, VS and ADSL
pc <- pharmaversesdtm::pc
ex <- pharmaversesdtm::ex
vs <- pharmaversesdtm::vs
adsl <- admiral::admiral_adsl

# When SAS datasets are imported into R using haven::read_sas(), missing
# character values from SAS appear as "" characters in R, instead of appearing
# as NA values. Further details can be obtained via the following link:
# https://pharmaverse.github.io/admiral/cran-release/articles/admiral.html#handling-of-missing-values # nolint

ex <- convert_blanks_to_na(ex)
pc <- convert_blanks_to_na(pc)
vs <- convert_blanks_to_na(vs)

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
    time_imputation = "00:00:00",
    ignore_seconds_flag = FALSE
  ) %>%
  # Derive dates and times from date/times
  derive_vars_dtm_to_dt(exprs(ADTM)) %>%
  derive_vars_dtm_to_tm(exprs(ADTM)) %>%
  derive_vars_dy(reference_date = TRTSDT, source_vars = exprs(ADT)) %>%
  # Derive event ID and nominal relative time from first dose (NFRLT)
  mutate(
    EVID = 0,
    DRUG = PCTEST,
  ) %>%
  derive_var_nfrlt(
    new_var = NFRLT,
    new_var_unit = FRLTU,
    out_unit = "HOURS",
    tpt_var = PCTPT,
    visit_day = VISITDY,
    treatment_duration = 0,
    range_method = "midpoint"
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
    EVID = 1
  ) %>%
  derive_var_nfrlt(
    new_var = NFRLT,
    new_var_unit = FRLTU,
    out_unit = "HOURS",
    visit_day = VISITDY
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
  derive_vars_dtm_to_tm(exprs(AENDTM)) %>%
  derive_vars_dy(reference_date = TRTSDT, source_vars = exprs(ADT))


# ---- Find first dose per treatment per subject ----
# ---- Join with ADPC data and keep only subjects with dosing ----

adpc_first_dose <- pc_dates %>%
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

adpc_prev <- adpc_first_dose %>%
  derive_vars_joined(
    dataset_add = ex_exp,
    by_vars = exprs(USUBJID),
    order = exprs(ADTM),
    new_vars = exprs(
      ADTM_prev = ADTM, EXDOSE_prev = EXDOSE, AVISIT_prev = AVISIT,
      AENDTM_prev = AENDTM
    ),
    join_vars = exprs(ADTM),
    join_type = "all",
    filter_add = NULL,
    filter_join = ADTM > ADTM.join,
    mode = "last",
    check_type = "none"
  )

# ---- Find next dose  ----

adpc_next <- adpc_prev %>%
  derive_vars_joined(
    dataset_add = ex_exp,
    by_vars = exprs(USUBJID),
    order = exprs(ADTM),
    new_vars = exprs(
      ADTM_next = ADTM, EXDOSE_next = EXDOSE, AVISIT_next = AVISIT,
      AENDTM_next = AENDTM
    ),
    join_vars = exprs(ADTM),
    join_type = "all",
    filter_add = NULL,
    filter_join = ADTM <= ADTM.join,
    mode = "first",
    check_type = "none"
  )

# ---- Find previous nominal time ----

adpc_nom_prev <- adpc_next %>%
  derive_vars_joined(
    dataset_add = ex_exp,
    by_vars = exprs(USUBJID),
    order = exprs(NFRLT),
    new_vars = exprs(NFRLT_prev = NFRLT),
    join_vars = exprs(NFRLT),
    join_type = "all",
    filter_add = NULL,
    filter_join = NFRLT > NFRLT.join,
    mode = "last",
    check_type = "none"
  )

# ---- Find next nominal time ----

adpc_nom_next <- adpc_nom_prev %>%
  derive_vars_joined(
    dataset_add = ex_exp,
    by_vars = exprs(USUBJID),
    order = exprs(NFRLT),
    new_vars = exprs(NFRLT_next = NFRLT),
    join_vars = exprs(NFRLT),
    join_type = "all",
    filter_add = NULL,
    filter_join = NFRLT <= NFRLT.join,
    mode = "first",
    check_type = "none"
  )

# ---- Combine ADPC and EX data ----
# Derive Relative Time Variables

adpc_arrlt <- bind_rows(adpc_nom_next, ex_exp) %>%
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
    out_unit = "HOURS",
    floor_in = FALSE,
    add_one = FALSE
  ) %>%
  # Derive Actual Relative Time from Reference Dose (ARRLT)
  derive_vars_duration(
    new_var = ARRLT,
    new_var_unit = RRLTU,
    start_date = ADTM_prev,
    end_date = ADTM,
    out_unit = "HOURS",
    floor_in = FALSE,
    add_one = FALSE
  ) %>%
  # Derive Actual Relative Time from Next Dose (AXRLT not kept)
  derive_vars_duration(
    new_var = AXRLT,
    start_date = ADTM_next,
    end_date = ADTM,
    out_unit = "HOURS",
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
  ) %>%
  # Derive dates and times from datetimes
  derive_vars_dtm_to_dt(exprs(FANLDTM)) %>%
  derive_vars_dtm_to_tm(exprs(FANLDTM)) %>%
  derive_vars_dtm_to_dt(exprs(PCRFTDTM)) %>%
  derive_vars_dtm_to_tm(exprs(PCRFTDTM))

# Derive Nominal Relative Time from Reference Dose (NRRLT)

adpc_nrrlt <- adpc_arrlt %>%
  mutate(
    NRRLT = case_when(
      EVID == 1 ~ 0,
      is.na(NFRLT_prev) ~ NFRLT - min_NFRLT,
      TRUE ~ NFRLT - NFRLT_prev
    ),
    NXRLT = case_when(
      EVID == 1 ~ 0,
      TRUE ~ NFRLT - NFRLT_next
    )
  )

# ---- Derive Analysis Variables ----
# Derive ATPTN, ATPT, ATPTREF, ABLFL and BASETYPE
# Derive planned dose DOSEP, actual dose DOSEA and units
# Derive PARAMCD and relative time units
# Derive AVAL, AVALU and AVALCAT1

adpc_aval <- adpc_nrrlt %>%
  mutate(
    PARCAT1 = PCSPEC,
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
    # Derive BASETYPE
    BASETYPE = paste(ATPTREF, "Baseline"),

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
    DOSEU = "mg",
  ) %>%
  # Derive relative time units
  mutate(
    # Derive PARAMCD
    PARAMCD = coalesce(PCTESTCD, "DOSE"),
    ALLOQ = PCLLOQ,
    # Derive AVAL
    AVAL = case_when(
      EVID == 1 ~ EXDOSE,
      TRUE ~ PCSTRESN
    ),
    AVALU = case_when(
      EVID == 1 ~ EXDOSU,
      TRUE ~ PCSTRESU
    ),
    AVALCAT1 = if_else(PCSTRESC == "<BLQ", PCSTRESC, prettyNum(signif(AVAL, digits = 3))),
  ) %>%
  # Add SRCSEQ
  mutate(
    SRCDOM = DOMAIN,
    SRCVAR = "SEQ",
    SRCSEQ = coalesce(PCSEQ, EXSEQ)
  )

# ---- Create DTYPE imputation

adpc_lloq <- adpc_aval %>%
  derive_extreme_records(
    dataset_add = adpc_aval,
    by_vars = exprs(USUBJID, PARAMCD, PARCAT1, AVISITN, AVISIT, ADTM, PCSEQ),
    order = exprs(ADTM, BASETYPE, EVID, ATPTN, PARCAT1),
    mode = "last",
    filter_add = PCSTRESC == "<BLQ" & is.na(AVAL),
    set_values_to = exprs(
      AVAL = ALLOQ * .5,
      DTYPE = "HALFLLOQ"
    )
  )

# ---- Create DTYPE copy records ----

dtype <- adpc_lloq %>%
  filter(NFRLT > 0 & NXRLT == 0 & EVID == 0 & !is.na(AVISIT_next)) %>%
  select(-PCRFTDT, -PCRFTTM) %>%
  # Re-derive variables in for DTYPE copy records
  mutate(
    ATPTREF = AVISIT_next,
    ARRLT = AXRLT,
    NRRLT = NXRLT,
    PCRFTDTM = ADTM_next,
    DOSEA = EXDOSE_next,
    BASETYPE = paste(AVISIT_next, "Baseline"),
    ATPT = "Pre-dose",
    ATPTN = -0.5,
    DTYPE = case_when(
      is.na(DTYPE) ~ "COPY",
      DTYPE == "HALFLLOQ" ~ "COPY/HALFLLOQ"
    )
  ) %>%
  derive_vars_dtm_to_dt(exprs(PCRFTDTM)) %>%
  derive_vars_dtm_to_tm(exprs(PCRFTDTM))

# ---- Combine original records and DTYPE copy records ----

adpc_dtype <- bind_rows(adpc_lloq, dtype) %>%
  arrange(STUDYID, USUBJID, BASETYPE, ADTM, NFRLT) %>%
  mutate(
    # Derive baseline flag for pre-dose records
    ABLFL = case_when(
      ATPT == "Pre-dose" & !is.na(AVAL) ~ "Y",
      TRUE ~ NA_character_
    ),
    # Derive MRRLT, ANL01FL and ANL02FL
    MRRLT = if_else(ARRLT < 0, 0, ARRLT),
    ANL01FL = "Y",
    ANL02FL = if_else(is.na(DTYPE), "Y", NA_character_)
  )

# ---- Derive BASE and Calculate Change from Baseline ----

adpc_base <- adpc_dtype %>%
  derive_var_base(
    by_vars = exprs(STUDYID, USUBJID, PARAMCD, PARCAT1, BASETYPE),
    source_var = AVAL,
    new_var = BASE,
    filter = ABLFL == "Y"
  )

# Calculate CHG for post-baseline records
# The decision on how to populate pre-baseline and baseline values of CHG is left to producer choice
adpc_chg <- restrict_derivation(
  adpc_base,
  derivation = derive_var_chg,
  filter = AVISITN > 0
)

# ---- Add ASEQ ----

adpc_aseq <- adpc_chg %>%
  # Calculate ASEQ
  derive_var_obs_number(
    new_var = ASEQ,
    by_vars = exprs(STUDYID, USUBJID),
    order = exprs(ADTM, BASETYPE, EVID, AVISITN, ATPTN, PARCAT1, DTYPE),
    check_type = "error"
  ) %>%
  # Remove temporary variables
  select(
    -DOMAIN, -PCSEQ, -starts_with("min"),
    -starts_with("max"), -starts_with("EX"), -ends_with("next"),
    -ends_with("prev"), -DRUG, -EVID, -AXRLT, -NXRLT
  ) %>%
  # Derive PARAM and PARAMN
  derive_vars_merged(dataset_add = select(param_lookup, -PCTESTCD), by_vars = exprs(PARAMCD))


#---- Derive additional baselines from VS ----

adpc_baselines <- adpc_aseq %>%
  derive_vars_merged(
    dataset_add = vs,
    filter_add = VSTESTCD == "HEIGHT",
    by_vars = exprs(STUDYID, USUBJID),
    new_vars = exprs(HTBL = VSSTRESN, HTBLU = VSSTRESU)
  ) %>%
  derive_vars_merged(
    dataset_add = vs,
    filter_add = VSTESTCD == "WEIGHT" & VSBLFL == "Y",
    by_vars = exprs(STUDYID, USUBJID),
    new_vars = exprs(WTBL = VSSTRESN, WTBLU = VSSTRESU)
  ) %>%
  mutate(
    BMIBL = compute_bmi(height = HTBL, weight = WTBL),
    BMIBLU = "kg/m^2"
  )

# ---- Add all ADSL variables ----

# Add all ADSL variables
adpc <- adpc_baselines %>%
  derive_vars_merged(
    dataset_add = select(adsl, !!!negate_vars(adsl_vars)),
    by_vars = exprs(STUDYID, USUBJID)
  )

# Final Steps, Select final variables and Add labels
# This process will be based on your metadata, no example given for this reason
# ...
# ---- Save output ----

# Change to whichever directory you want to save the dataset in
dir <- tools::R_user_dir("admiral_templates_data", which = "cache")
if (!file.exists(dir)) {
  # Create the folder
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
}
save(adpc, file = file.path(dir, "adpc.rda"), compress = "bzip2")
```

ad_adpp.R

```
# Name: ADPP
#
# Label: Pharmacokinetics Parameters Analysis Dataset
#
# Description: Based on simulated data, create ADPP analysis dataset
#
# Input: pp, adsl
library(admiral)
library(pharmaversesdtm) # Contains example datasets from the CDISC pilot project
library(dplyr)
library(lubridate)
library(stringr)


# Load source datasets ----

# Use e.g. haven::read_sas to read in .sas7bdat, or other suitable functions
# as needed and assign to the variables below.
# For illustration purposes read in admiral test data

# Load PP and Adsl
pp <- pharmaversesdtm::pp
admiral_adsl <- admiral::admiral_adsl

# When SAS datasets are imported into R using haven::read_sas(), missing
# character values from SAS appear as "" characters in R, instead of appearing
# as NA values. Further details can be obtained via the following link:
# https://pharmaverse.github.io/admiral/cran-release/articles/admiral.html#handling-of-missing-values # nolint

pp <- convert_blanks_to_na(pp)

# Lookup tables ----
param_lookup <- tibble::tribble(
  ~PPTESTCD,  ~PARAMCD,                                     ~PARAM, ~PARAMN,
  "AUCALL",   "AUCALL",                                  "AUC All",       1,
  "AUCIFO",   "AUCIFO",                         "AUC Infinity Obs",       2,
  "AUCIFOD", "AUCIFOD",            "AUC Infinity Obs Norm by Dose",       3,
  "AUCINT",   "AUCINT",                        "AUC from T1 to T2",       4,
  "AUCLST",   "AUCLST",                 "AUC to Last Nonzero Conc",       5,
  "AUCPEO",   "AUCPEO",                   "AUC %Extrapolation Obs",       6,
  "CEND",       "CEND", "Concentration at the end of the infusion",       7, # non Cdisc Term
  "CLO",         "CLO",                             "Total CL Obs",       8,
  "CLST",       "CLST",                        "Last Nonzero Conc",       9,
  "CMAX",       "CMAX",                                 "Max Conc",      10,
  "CMAXD",     "CMAXD",                    "Max Conc Norm by Dose",      11,
  "CSF",         "CSF",                      "CSF to Plasma Ratio",      12, # non Cdisc Term
  "LAMZ",       "LAMZ",                                 "Lambda z",      13,
  "LAMZHL",   "LAMZHL",                       "Half-Life Lambda z",      14,
  "LAMZNPT", "LAMZNPT",            "Number of Points for Lambda z",      15,
  "R2ADJ",     "R2ADJ",                       "R Squared Adjusted",      16,
  "TCEND",     "TCEND",                             "Time of CEND",      17, # non Cdisc Term
  "TLST",       "TLST",                "Time of Last Nonzero Conc",      18,
  "TMAX",       "TMAX",                             "Time of CMAX",      19,
  "VSSO",       "VSSO",                "Vol Dist Steady State Obs",      20,
  "RCAMINT", "RCAMINT",                                       "Ae",      21,
  "RENALCL", "RENALCL",                                      "CLR",      22
)

# Assign AVALCATx
avalcax_lookup <- exprs(
  ~PARAMCD, ~condition,       ~AVALCAT1, ~AVALCA1N,
  "AUCALL",  AVAL < 19,           "<19",         1,
  "AUCALL", AVAL >= 19,          ">=19",         2
)

attr(param_lookup$PPTESTCD, "label") <- "Parameter Short Name"

# Derivations ----

# Get list of ADSL vars required for derivations
adsl_vars <- exprs(TRTSDT, TRTEDT, DTHDT, EOSDT, TRT01P, TRT01A)

adpp_pp <- pp %>%
  # Join ADSL with PP (need TRTSDT for ADY derivation)
  derive_vars_merged(
    dataset_add = admiral_adsl,
    new_vars = adsl_vars,
    by_vars = exprs(STUDYID, USUBJID)
  ) %>%
  ## Calculate ADT, ADY ----
  derive_vars_dt(
    new_vars_prefix = "A",
    dtc = PPRFDTC
  ) %>%
  derive_vars_dy(reference_date = TRTSDT, source_vars = exprs(ADT))

adpp_aval <- adpp_pp %>%
  ## Add PARAMCD only - add PARAM etc later ----
  left_join(
    select(param_lookup, PPTESTCD, PARAMCD),
    by = "PPTESTCD"
  ) %>%
  ## Calculate PARCAT1, AVAL and AVALC ----
  mutate(
    PARCAT1 = PPCAT,
    AVAL = PPSTRESN,
    AVALU = PPSTRESU,
  ) %>%
  # Remove variables
  select(-PPSTRESN, -PPSTRESC) %>%
  # Add ASEQ
  mutate(
    SRCDOM = DOMAIN,
    SRCVAR = "SEQ",
    SRCSEQ = PPSEQ
  ) %>%
  select(-DOMAIN, -PPSEQ)

## Get visit info ----
# See also the "Visit and Period Variables" vignette
# (https://pharmaverse.github.io/admiral/cran-release/articles/visits_periods.html#visit_bds)
adpp_avisit <- adpp_aval %>%
  # Derive Timing
  mutate(
    AVISITN = ADY,
    AVISIT = paste("Day", ADY),
    VISITNUM = AVISITN,
    VISIT = AVISIT
  ) %>%
  ## Assign TRTA, TRTP ----
  # See also the "Visit and Period Variables" vignette
  # (https://pharmaverse.github.io/admiral/cran-release/articles/visits_periods.html#treatment_bds)
  mutate(
    TRTP = TRT01P,
    TRTA = TRT01A
  ) %>%
  ## Derive AVALCA1N and AVALCAT1 ----
  derive_vars_cat(
    definition = avalcax_lookup,
    by_vars = exprs(PARAMCD)
  )
# Add all ADSL variables
adpp <- adpp_avisit %>%
  derive_vars_merged(
    dataset_add = select(admiral_adsl, !!!negate_vars(adsl_vars)),
    by_vars = exprs(STUDYID, USUBJID)
  )

# Final Steps, Select final variables and Add labels
# This process will be based on your metadata, no example given for this reason
# ...
# Save output ----

# Change to whichever directory you want to save the dataset in
dir <- tools::R_user_dir("admiral_templates_data", which = "cache")
if (!file.exists(dir)) {
  # Create the folder
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
}
save(adpp, file = file.path(dir, "adpp.rda"), compress = "bzip2")
```

ad_adppk.R

```
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
library(pharmaversesdtm) # Contains example datasets from the CDISC pilot project or simulated

# ---- Load source datasets ----

# Use e.g. haven::read_sas to read in .sas7bdat, or other suitable functions
# as needed and assign to the variables below.
# For illustration purposes read in admiral test data

# Load PC, EX, VS, LB and ADSL
pc <- pharmaversesdtm::pc
ex <- pharmaversesdtm::ex
vs <- pharmaversesdtm::vs
lb <- pharmaversesdtm::lb
adsl <- admiral::admiral_adsl

# When SAS datasets are imported into R using haven::read_sas(), missing
# character values from SAS appear as "" characters in R, instead of appearing
# as NA values. Further details can be obtained via the following link:
# https://pharmaverse.github.io/admiral/cran-release/articles/admiral.html#handling-of-missing-values # nolint

ex <- convert_blanks_to_na(ex)
pc <- convert_blanks_to_na(pc)
vs <- convert_blanks_to_na(vs)
lb <- convert_blanks_to_na(lb)

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
    time_imputation = "00:00:00",
    ignore_seconds_flag = FALSE
  ) %>%
  # Derive dates and times from date/times
  derive_vars_dtm_to_dt(exprs(ADTM)) %>%
  derive_vars_dtm_to_tm(exprs(ADTM)) %>%
  # Derive event ID and nominal relative time from first dose (NFRLT)
  mutate(
    EVID = 0,
    DRUG = PCTEST
  ) %>%
  derive_var_nfrlt(
    new_var = NFRLT,
    out_unit = "HOURS",
    tpt_var = PCTPT,
    visit_day = VISITDY,
    treatment_duration = 0,
    range_method = "midpoint"
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
    EVID = 1
  ) %>%
  derive_var_nfrlt(
    new_var = NFRLT,
    out_unit = "HOURS",
    visit_day = VISITDY
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
    join_type = "all",
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
    join_type = "all",
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
    out_unit = "HOURS",
    floor_in = FALSE,
    add_one = FALSE
  ) %>%
  # Derive Actual Relative Time from Reference Dose (APRLT)
  derive_vars_duration(
    new_var = APRLT,
    start_date = ADTM_prev,
    end_date = ADTM,
    out_unit = "HOURS",
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
      PCSPEC == "PLASMA" ~ 2,
      TRUE ~ 3
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
    order = exprs(AFRLT, EVID, CMT),
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
    -starts_with("AVISIT"), -ends_with("TMF"), -starts_with("TRT"),
    -starts_with("ATPT"), -DRUG
  )

#---- Derive Covariates ----
# Include numeric values for STUDYIDN, USUBJIDN, SEXN, RACEN etc.

covar <- adsl %>%
  derive_vars_merged(
    dataset_add = country_code_lookup,
    new_vars = exprs(COUNTRYN = country_number, COUNTRYL = country_name),
    by_vars = exprs(COUNTRY = country_code),
  ) %>%
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
    )
  ) %>%
  select(
    STUDYID, STUDYIDN, SITEID, SITEIDN, USUBJID, USUBJIDN,
    SUBJID, SUBJIDN, AGE, SEX, SEXN, COHORT, COHORTC, ROUTE, ROUTEN,
    RACE, RACEN, ETHNIC, ETHNICN, FORM, FORMN, COUNTRY, COUNTRYN, COUNTRYL
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
      weight = WTBL,
      method = "Mosteller"
    ),
    CRCLBL = compute_egfr(
      creat = CREATBL, creatu = "SI", age = AGE, weight = WTBL, sex = SEX,
      method = "CRCL"
    ),
    EGFRBL = compute_egfr(
      creat = CREATBL, creatu = "SI", age = AGE, weight = WTBL, sex = SEX,
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

# Change to whichever directory you want to save the dataset in
dir <- tools::R_user_dir("admiral_templates_data", which = "cache")
if (!file.exists(dir)) {
  # Create the folder
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
}
save(adppk, file = file.path(dir, "adppk.rda"), compress = "bzip2")
```

ad_adsl.R

```
# Name: ADSL
#
# Label: Subject Level Analysis Dataset
#
# Input: dm, ds, ex, ae, lb
library(admiral)
library(pharmaversesdtm) # Contains example datasets from the CDISC pilot project
library(dplyr)
library(lubridate)
library(stringr)

# Load source datasets ----

# Use e.g. haven::read_sas to read in .sas7bdat, or other suitable functions
# as needed and assign to the variables below.
# For illustration purposes read in admiral test data.

dm <- pharmaversesdtm::dm
ds <- pharmaversesdtm::ds
ex <- pharmaversesdtm::ex
ae <- pharmaversesdtm::ae
lb <- pharmaversesdtm::lb

# When SAS datasets are imported into R using haven::read_sas(), missing
# character values from SAS appear as "" characters in R, instead of appearing
# as NA values. Further details can be obtained via the following link:
# https://pharmaverse.github.io/admiral/cran-release/articles/admiral.html#handling-of-missing-values # nolint


dm <- convert_blanks_to_na(dm)
ds <- convert_blanks_to_na(ds)
ex <- convert_blanks_to_na(ex)
ae <- convert_blanks_to_na(ae)
lb <- convert_blanks_to_na(lb)

# User defined functions ----

# Here are some examples of how you can create your own functions that
#  operates on vectors, which can be used in `mutate`.

# Grouping
format_racegr1 <- function(x) {
  case_when(
    x == "WHITE" ~ "White",
    x != "WHITE" ~ "Non-white",
    TRUE ~ "Missing"
  )
}

format_agegr1 <- function(x) {
  case_when(
    x < 18 ~ "<18",
    between(x, 18, 64) ~ "18-64",
    x > 64 ~ ">64",
    TRUE ~ "Missing"
  )
}

format_region1 <- function(x) {
  case_when(
    x %in% c("CAN", "USA") ~ "NA",
    !is.na(x) ~ "RoW",
    TRUE ~ "Missing"
  )
}

format_lddthgr1 <- function(x) {
  case_when(
    x <= 30 ~ "<= 30",
    x > 30 ~ "> 30",
    TRUE ~ NA_character_
  )
}

# EOSSTT mapping
format_eosstt <- function(x) {
  case_when(
    x %in% c("COMPLETED") ~ "COMPLETED",
    x %in% c("SCREEN FAILURE") ~ NA_character_,
    !is.na(x) ~ "DISCONTINUED",
    TRUE ~ "ONGOING"
  )
}

# Derivations ----
# impute start and end time of exposure to first and last respectively, do not impute date
ex_ext <- ex %>%
  derive_vars_dtm(
    dtc = EXSTDTC,
    new_vars_prefix = "EXST"
  ) %>%
  derive_vars_dtm(
    dtc = EXENDTC,
    new_vars_prefix = "EXEN",
    time_imputation = "last"
  )

adsl <- dm %>%
  ## derive treatment variables (TRT01P, TRT01A) ----
  # See also the "Visit and Period Variables" vignette
  # (https://pharmaverse.github.io/admiral/cran-release/articles/visits_periods.html#treatment_adsl)
  mutate(TRT01P = ARM, TRT01A = ACTARM) %>%
  ## derive treatment start date (TRTSDTM) ----
  derive_vars_merged(
    dataset_add = ex_ext,
    filter_add = (EXDOSE > 0 |
      (EXDOSE == 0 &
        str_detect(EXTRT, "PLACEBO"))) &
      !is.na(EXSTDTM),
    new_vars = exprs(TRTSDTM = EXSTDTM, TRTSTMF = EXSTTMF),
    order = exprs(EXSTDTM, EXSEQ),
    mode = "first",
    by_vars = exprs(STUDYID, USUBJID)
  ) %>%
  ## derive treatment end date (TRTEDTM) ----
  derive_vars_merged(
    dataset_add = ex_ext,
    filter_add = (EXDOSE > 0 |
      (EXDOSE == 0 &
        str_detect(EXTRT, "PLACEBO"))) & !is.na(EXENDTM),
    new_vars = exprs(TRTEDTM = EXENDTM, TRTETMF = EXENTMF),
    order = exprs(EXENDTM, EXSEQ),
    mode = "last",
    by_vars = exprs(STUDYID, USUBJID)
  ) %>%
  ## Derive treatment end/start date TRTSDT/TRTEDT ----
  derive_vars_dtm_to_dt(source_vars = exprs(TRTSDTM, TRTEDTM)) %>%
  ## derive treatment duration (TRTDURD) ----
  derive_var_trtdurd()

## Disposition dates, status ----
# convert character date to numeric date without imputation
ds_ext <- derive_vars_dt(
  ds,
  dtc = DSSTDTC,
  new_vars_prefix = "DSST"
)

# Screen fail date
adsl <- adsl %>%
  derive_vars_merged(
    dataset_add = ds_ext,
    by_vars = exprs(STUDYID, USUBJID),
    new_vars = exprs(SCRFDT = DSSTDT),
    filter_add = DSCAT == "DISPOSITION EVENT" & DSDECOD == "SCREEN FAILURE"
  ) %>%
  derive_vars_merged(
    dataset_add = ds_ext,
    by_vars = exprs(STUDYID, USUBJID),
    new_vars = exprs(EOSDT = DSSTDT),
    filter_add = DSCAT == "DISPOSITION EVENT" & DSDECOD != "SCREEN FAILURE"
  ) %>%
  # EOS status
  derive_vars_merged(
    dataset_add = ds_ext,
    by_vars = exprs(STUDYID, USUBJID),
    filter_add = DSCAT == "DISPOSITION EVENT",
    new_vars = exprs(EOSSTT = format_eosstt(DSDECOD)),
    missing_values = exprs(EOSSTT = "ONGOING")
  ) %>%
  # Last retrieval date
  derive_vars_merged(
    dataset_add = ds_ext,
    by_vars = exprs(STUDYID, USUBJID),
    new_vars = exprs(FRVDT = DSSTDT),
    filter_add = DSCAT == "OTHER EVENT" & DSDECOD == "FINAL RETRIEVAL VISIT"
  ) %>%
  # Derive Randomization Date
  derive_vars_merged(
    dataset_add = ds_ext,
    filter_add = DSDECOD == "RANDOMIZED",
    by_vars = exprs(STUDYID, USUBJID),
    new_vars = exprs(RANDDT = DSSTDT)
  ) %>%
  # Death date - impute partial date to first day/month
  derive_vars_dt(
    new_vars_prefix = "DTH",
    dtc = DTHDTC,
    highest_imputation = "M",
    date_imputation = "first"
  ) %>%
  # Relative Day of Death
  derive_vars_duration(
    new_var = DTHADY,
    start_date = TRTSDT,
    end_date = DTHDT
  ) %>%
  # Elapsed Days from Last Dose to Death
  derive_vars_duration(
    new_var = LDDTHELD,
    start_date = TRTEDT,
    end_date = DTHDT,
    add_one = FALSE
  ) %>%
  # Cause of Death and Traceability Variables
  derive_vars_extreme_event(
    by_vars = exprs(STUDYID, USUBJID),
    events = list(
      event(
        dataset_name = "ae",
        condition = AEOUT == "FATAL",
        set_values_to = exprs(DTHCAUS = AEDECOD, DTHDOM = DOMAIN),
      ),
      event(
        dataset_name = "ds",
        condition = DSDECOD == "DEATH" & grepl("DEATH DUE TO", DSTERM),
        set_values_to = exprs(DTHCAUS = DSTERM, DTHDOM = DOMAIN),
      )
    ),
    source_datasets = list(ae = ae, ds = ds),
    tmp_event_nr_var = event_nr,
    order = exprs(event_nr),
    mode = "first",
    new_vars = exprs(DTHCAUS = DTHCAUS, DTHDOM = DTHDOM)
  ) %>%
  # Death Cause Category
  mutate(DTHCGR1 = case_when(
    is.na(DTHDOM) ~ NA_character_,
    DTHDOM == "AE" ~ "ADVERSE EVENT",
    str_detect(DTHCAUS, "(PROGRESSIVE DISEASE|DISEASE RELAPSE)") ~ "PROGRESSIVE DISEASE",
    TRUE ~ "OTHER"
  ))

## Last known alive date ----
## DTC variables are converted to numeric dates imputing missing day and month
## to the first

adsl <- adsl %>%
  derive_vars_extreme_event(
    by_vars = exprs(STUDYID, USUBJID),
    events = list(
      event(
        dataset_name = "ae",
        order = exprs(AESTDTC, AESEQ),
        condition = !is.na(AESTDTC),
        set_values_to = exprs(
          LSTALVDT = convert_dtc_to_dt(AESTDTC, highest_imputation = "M"),
          seq = AESEQ
        ),
      ),
      event(
        dataset_name = "ae",
        order = exprs(AEENDTC, AESEQ),
        condition = !is.na(AEENDTC),
        set_values_to = exprs(
          LSTALVDT = convert_dtc_to_dt(AEENDTC, highest_imputation = "M"),
          seq = AESEQ
        ),
      ),
      event(
        dataset_name = "lb",
        order = exprs(LBDTC, LBSEQ),
        condition = !is.na(LBDTC),
        set_values_to = exprs(
          LSTALVDT = convert_dtc_to_dt(LBDTC, highest_imputation = "M"),
          seq = LBSEQ
        ),
      ),
      event(
        dataset_name = "adsl",
        condition = !is.na(TRTEDT),
        set_values_to = exprs(LSTALVDT = TRTEDT, seq = NA_integer_),
      )
    ),
    source_datasets = list(ae = ae, lb = lb, adsl = adsl),
    tmp_event_nr_var = event_nr,
    order = exprs(LSTALVDT, seq, event_nr),
    mode = "last",
    new_vars = exprs(LSTALVDT)
  ) %>%
  derive_var_merged_exist_flag(
    dataset_add = ex,
    by_vars = exprs(STUDYID, USUBJID),
    new_var = SAFFL,
    false_value = "N",
    missing_value = "N",
    condition = (EXDOSE > 0 | (EXDOSE == 0 & str_detect(EXTRT, "PLACEBO")))
  ) %>%
  ## Groupings and others variables ----
  mutate(
    RACEGR1 = format_racegr1(RACE),
    AGEGR1 = format_agegr1(AGE),
    REGION1 = format_region1(COUNTRY),
    LDDTHGR1 = format_lddthgr1(LDDTHELD),
    DTH30FL = if_else(LDDTHGR1 == "<= 30", "Y", NA_character_),
    DTHA30FL = if_else(LDDTHGR1 == "> 30", "Y", NA_character_),
    DTHB30FL = if_else(DTHDT <= TRTSDT + 30, "Y", NA_character_),
    DOMAIN = NULL
  )


# Save output ----

# Change to whichever directory you want to save the dataset in
dir <- tools::R_user_dir("admiral_templates_data", which = "cache")
if (!file.exists(dir)) {
  # Create the folder
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
}
save(adsl, file = file.path(dir, "adsl.rda"), compress = "bzip2")
```

ad_advs.R

```
# Name: ADVS
#
# Label: Vital Signs Analysis Dataset
#
# Input: adsl, vs
library(admiral)
library(pharmaversesdtm) # Contains example datasets from the CDISC pilot project
library(dplyr)
library(lubridate)
library(stringr)

# Load source datasets ----

# Use e.g. `haven::read_sas()` to read in .sas7bdat, or other suitable functions
# as needed and assign to the variables below.
# For illustration purposes read in admiral test data

vs <- pharmaversesdtm::vs
adsl <- admiral::admiral_adsl

# When SAS datasets are imported into R using haven::read_sas(), missing
# character values from SAS appear as "" characters in R, instead of appearing
# as NA values. Further details can be obtained via the following link:
# https://pharmaverse.github.io/admiral/cran-release/articles/admiral.html#handling-of-missing-values # nolint

vs <- convert_blanks_to_na(vs)

# Lookup tables ----

# Assign PARAMCD, PARAM, and PARAMN
param_lookup <- tibble::tribble(
  ~VSTESTCD, ~PARAMCD,                            ~PARAM, ~PARAMN,
  "SYSBP",    "SYSBP",  "Systolic Blood Pressure (mmHg)",       1,
  "DIABP",    "DIABP", "Diastolic Blood Pressure (mmHg)",       2,
  "PULSE",    "PULSE",          "Pulse Rate (beats/min)",       3,
  "WEIGHT",  "WEIGHT",                     "Weight (kg)",       4,
  "HEIGHT",  "HEIGHT",                     "Height (cm)",       5,
  "TEMP",      "TEMP",                 "Temperature (C)",       6,
  "MAP",        "MAP",   "Mean Arterial Pressure (mmHg)",       7,
  "BMI",        "BMI",         "Body Mass Index(kg/m^2)",       8,
  "BSA",        "BSA",          "Body Surface Area(m^2)",       9
)
attr(param_lookup$VSTESTCD, "label") <- "Vital Signs Test Short Name"


# Assign ANRLO/HI, A1LO/HI
range_lookup <- tibble::tribble(
  ~PARAMCD, ~ANRLO, ~ANRHI, ~A1LO, ~A1HI,
  "SYSBP",      90,    130,    70,   140,
  "DIABP",      60,     80,    40,    90,
  "PULSE",      60,    100,    40,   110,
  "TEMP",     36.5,   37.5,    35,    38
)

# Assign AVALCATx
avalcax_lookup <- exprs(
  ~PARAMCD,  ~condition,  ~AVALCAT1, ~AVALCA1N,
  "HEIGHT",  AVAL > 100,  ">100 cm",         1,
  "HEIGHT", AVAL <= 100, "<=100 cm",         2
)

# Derivations ----

# Get list of ADSL vars required for derivations
adsl_vars <- exprs(TRTSDT, TRTEDT, TRT01A, TRT01P)

advs <- vs %>%
  # Join ADSL with VS (need TRTSDT for ADY derivation)
  derive_vars_merged(
    dataset_add = adsl,
    new_vars = adsl_vars,
    by_vars = exprs(STUDYID, USUBJID)
  ) %>%
  ## Calculate ADT, ADY ----
  derive_vars_dt(
    new_vars_prefix = "A",
    dtc = VSDTC
  ) %>%
  derive_vars_dy(reference_date = TRTSDT, source_vars = exprs(ADT))

advs <- advs %>%
  ## Add PARAMCD only - add PARAM etc later ----
  derive_vars_merged_lookup(
    dataset_add = param_lookup,
    new_vars = exprs(PARAMCD),
    by_vars = exprs(VSTESTCD)
  ) %>%
  ## Calculate AVAL and AVALC ----
  # AVALC should only be mapped if it contains non-redundant information.
  mutate(
    # AVALC = VSSTRESC,
    AVAL = VSSTRESN
  ) %>%
  ## Derive new parameters based on existing records ----
  # Note that, for the following three `derive_param_*()` functions, only the
  # variables specified in `by_vars` will be populated in the newly created
  # records.

  # Derive Mean Arterial Pressure
  derive_param_map(
    by_vars = exprs(STUDYID, USUBJID, !!!adsl_vars, VISIT, VISITNUM, ADT, ADY, VSTPT, VSTPTNUM),
    set_values_to = exprs(PARAMCD = "MAP"),
    get_unit_expr = VSSTRESU,
    filter = VSSTAT != "NOT DONE" | is.na(VSSTAT)
  ) %>%
  # Derive Body Surface Area
  derive_param_bsa(
    by_vars = exprs(STUDYID, USUBJID, !!!adsl_vars, VISIT, VISITNUM, ADT, ADY, VSTPT, VSTPTNUM),
    method = "Mosteller",
    set_values_to = exprs(PARAMCD = "BSA"),
    get_unit_expr = VSSTRESU,
    filter = VSSTAT != "NOT DONE" | is.na(VSSTAT),
    constant_by_vars = exprs(USUBJID)
  ) %>%
  # Derive Body Mass Index
  derive_param_bmi(
    by_vars = exprs(STUDYID, USUBJID, !!!adsl_vars, VISIT, VISITNUM, ADT, ADY, VSTPT, VSTPTNUM),
    set_values_to = exprs(PARAMCD = "BMI"),
    get_unit_expr = VSSTRESU,
    filter = VSSTAT != "NOT DONE" | is.na(VSSTAT),
    constant_by_vars = exprs(USUBJID)
  )


## Get visit info ----
# See also the "Visit and Period Variables" vignette
# (https://pharmaverse.github.io/admiral/cran-release/articles/visits_periods.html#visits)
advs <- advs %>%
  # Derive Timing
  mutate(
    ATPTN = VSTPTNUM,
    ATPT = VSTPT,
    AVISIT = case_when(
      str_detect(VISIT, "SCREEN|UNSCHED|RETRIEVAL|AMBUL") ~ NA_character_,
      !is.na(VISIT) ~ str_to_title(VISIT),
      TRUE ~ NA_character_
    ),
    AVISITN = as.numeric(case_when(
      VISIT == "BASELINE" ~ "0",
      str_detect(VISIT, "WEEK") ~ str_trim(str_replace(VISIT, "WEEK", "")),
      TRUE ~ NA_character_
    ))
  )

## Derive a new record as a summary record (e.g. mean of the triplicates at each time point) ----
advs <- advs %>%
  derive_summary_records(
    dataset_add = advs,
    by_vars = exprs(STUDYID, USUBJID, !!!adsl_vars, PARAMCD, AVISITN, AVISIT, ADT, ADY),
    filter_add = !is.na(AVAL),
    set_values_to = exprs(
      AVAL = mean(AVAL),
      DTYPE = "AVERAGE"
    )
  )

advs <- advs %>%
  ## Calculate ONTRTFL ----
  derive_var_ontrtfl(
    start_date = ADT,
    ref_start_date = TRTSDT,
    ref_end_date = TRTEDT,
    filter_pre_timepoint = AVISIT == "Baseline"
  )

## Calculate ANRIND : requires the reference ranges ANRLO, ANRHI ----
# Also accommodates the ranges A1LO, A1HI
advs <- advs %>%
  derive_vars_merged(dataset_add = range_lookup, by_vars = exprs(PARAMCD)) %>%
  # Calculate ANRIND
  derive_var_anrind()

## Derive baseline flags ----
advs <- advs %>%
  # Calculate BASETYPE
  derive_basetype_records(
    basetypes = exprs(
      "LAST: AFTER LYING DOWN FOR 5 MINUTES" = ATPTN == 815,
      "LAST: AFTER STANDING FOR 1 MINUTE" = ATPTN == 816,
      "LAST: AFTER STANDING FOR 3 MINUTES" = ATPTN == 817,
      "LAST" = is.na(ATPTN)
    )
  ) %>%
  # Calculate ABLFL
  restrict_derivation(
    derivation = derive_var_extreme_flag,
    args = params(
      by_vars = exprs(STUDYID, USUBJID, BASETYPE, PARAMCD),
      order = exprs(ADT, VISITNUM, VSSEQ),
      new_var = ABLFL,
      mode = "last"
    ),
    filter = (!is.na(AVAL) &
      ADT <= TRTSDT & !is.na(BASETYPE) & is.na(DTYPE))
  )

## Derive baseline information ----
advs <- advs %>%
  # Calculate BASE
  derive_var_base(
    by_vars = exprs(STUDYID, USUBJID, PARAMCD, BASETYPE),
    source_var = AVAL,
    new_var = BASE
  ) %>%
  # Calculate BASEC
  # only if AVALC is mapped
  # derive_var_base(
  #   by_vars = exprs(STUDYID, USUBJID, PARAMCD, BASETYPE),
  #   source_var = AVALC,
  #   new_var = BASEC
  # ) %>%
  # Calculate BNRIND
  derive_var_base(
    by_vars = exprs(STUDYID, USUBJID, PARAMCD, BASETYPE),
    source_var = ANRIND,
    new_var = BNRIND
  ) %>%
  # Calculate CHG for post-baseline records
  # The decision on how to populate pre-baseline and baseline values of CHG is left to producer choice
  restrict_derivation(
    derivation = derive_var_chg,
    filter = AVISITN > 0
  ) %>%
  # Calculate PCHG for post-baseline records
  # The decision on how to populate pre-baseline and baseline values of PCHG is left to producer choice
  restrict_derivation(
    derivation = derive_var_pchg,
    filter = AVISITN > 0
  )

## ANL01FL: Flag last result within an AVISIT and ATPT for post-baseline records ----
advs <- advs %>%
  restrict_derivation(
    derivation = derive_var_extreme_flag,
    args = params(
      new_var = ANL01FL,
      by_vars = exprs(USUBJID, PARAMCD, AVISIT, ATPT, DTYPE),
      order = exprs(ADT, AVAL),
      mode = "last"
    ),
    filter = !is.na(AVISITN) & ONTRTFL == "Y"
  )

## Get treatment information ----
# See also the "Visit and Period Variables" vignette
# (https://pharmaverse.github.io/admiral/cran-release/articles/visits_periods.html#treatment_bds)
advs <- advs %>%
  # Assign TRTA, TRTP
  # Create End of Treatment Record
  derive_extreme_records(
    dataset_add = advs,
    by_vars = exprs(STUDYID, USUBJID, PARAMCD, ATPTN),
    order = exprs(ADT, AVISITN, AVAL),
    mode = "last",
    filter_add = (4 < AVISITN & AVISITN <= 13 & ANL01FL == "Y" & is.na(DTYPE)),
    set_values_to = exprs(
      AVISIT = "End of Treatment",
      AVISITN = 99,
      DTYPE = "LOV"
    )
  ) %>%
  mutate(
    TRTP = TRT01P,
    TRTA = TRT01A
  )

## Get ASEQ and AVALCATx and add PARAM/PARAMN ----
advs <- advs %>%
  # Calculate ASEQ
  derive_var_obs_number(
    new_var = ASEQ,
    by_vars = exprs(STUDYID, USUBJID),
    order = exprs(PARAMCD, ADT, AVISITN, VISITNUM, ATPTN, DTYPE),
    check_type = "error"
  ) %>%
  # Define condition and categories using derive_vars_cat
  derive_vars_cat(
    definition = avalcax_lookup,
    by_vars = exprs(PARAMCD)
  ) %>%
  # Derive PARAM and PARAMN
  derive_vars_merged(dataset_add = select(param_lookup, -VSTESTCD), by_vars = exprs(PARAMCD))



# Add all ADSL variables
advs <- advs %>%
  derive_vars_merged(
    dataset_add = select(adsl, !!!negate_vars(adsl_vars)),
    by_vars = exprs(STUDYID, USUBJID)
  )

# Final Steps, Select final variables and Add labels
# This process will be based on your metadata, no example given for this reason
# ...

# Save output ----

# Change to whichever directory you want to save the dataset in
dir <- tools::R_user_dir("admiral_templates_data", which = "cache")
if (!file.exists(dir)) {
  # Create the folder
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
}
save(advs, file = file.path(dir, "advs.rda"), compress = "bzip2")
```
