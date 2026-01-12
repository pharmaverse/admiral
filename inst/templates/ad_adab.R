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
