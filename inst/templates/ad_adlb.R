# Name: ADLB
#
# Label: Lab Analysis Dataset
#
# Input: adsl, lb
library(admiral)
library(admiral.test) # Contains example datasets from the CDISC pilot project
library(dplyr)
library(lubridate)
library(stringr)

# ---- Load source datasets ----

# Use e.g. haven::read_sas to read in .sas7bdat, or other suitable functions
# as needed and assign to the variables below.
# For illustration purposes read in admiral test data

data("lb")
data("adsl")

lb <- convert_blanks_to_na(lb)

# The CDISC Pilot Data contains no SUPPVS data
# If you have a SUPPLB then uncomment function below

# lb <- derive_vars_suppqual(lb, supplb) %>%

# ---- Derivations ----

# Get list of ADSL vars required for derivations
adsl_vars <- vars(TRTSDT, TRTEDT, TRT01A, TRT01P)

adlb <- lb %>%

  # Join ADSL with LB (need TRTSDT for ADY derivation)
  derive_vars_merged(
    dataset_add = adsl,
    new_vars = adsl_vars,
    by_vars = vars(STUDYID, USUBJID)
  ) %>%

  # Calculate ADT, ADY
  derive_vars_dt(
    new_vars_prefix = "A",
    dtc = LBDTC,
    flag_imputation = FALSE
  ) %>%

  derive_var_ady(reference_date = TRTSDT, date = ADT)

adlb <- adlb %>%
  # For labs in SI units
  # Add PARAMCD only - can create LOOK-UP table
  # PARAMCD currently only 8 chars
  # Calculate AVAL and AVALC
  mutate(
    PARCAT2 = "SI",
    PARAMCD = paste0(LBTESTCD, PARCAT2),
    PARAM = paste0(LBTEST, " (", LBSTRESU, ")"),
    AVAL  = LBSTRESN,
    AVALC = LBSTRESC,
    ANRLO = LBSTNRLO,
    ANRHI = LBSTNRHI
  )

# CSDISC data does not hold CV unit variables
# adlb_cv <- adlb %>%
#   For labs in CV (or US) units
#   Add PARAMCD only - can create LOOK-UP table later
#   PARAMCD currently only 8 chars
#   Calculate AVAL and AVALC
#  mutate(
#    PARCAT2 = "CV",
#    PARAMCD = paste0(LBTESTCD,PARCAT2),
#    PARAM = paste0(LBTEST, ' (',LBCVRESU,')'),
#    AVAL = LBCVRESN,
#    AVALC = LBCVRESC,
#    ANRLO = LBCVNRLO,
#    ANRHI = LBCVNRHI
#  )

# get visit info
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

  # Calculate ONTRTFL
  derive_var_ontrtfl(
    start_date = ADT,
    ref_start_date = TRTSDT,
    ref_end_date = TRTEDT,
    filter_pre_timepoint = AVISIT == "Baseline"
  )

# Calculate ANRIND : requires the reference ranges ANRLO, ANRHI
# Also accommodates the ranges A1LO, A1HI
adlb <- adlb %>%
  derive_var_anrind()

# Derive baseline flags
adlb <- adlb %>%
  # Calculate BASETYPE
  mutate(
    BASETYPE = "LAST"
  ) %>%

  # Calculate ABLFL
  derive_var_extreme_flag(
    by_vars = vars(STUDYID, USUBJID, BASETYPE, PARAMCD),
    order = vars(ADT, VISITNUM, LBSEQ),
    new_var = ABLFL,
    mode = "last",
    filter = (!is.na(AVAL) & ADT <= TRTSDT & !is.na(BASETYPE))
  )

# Derive baseline information
adlb <- adlb %>%

  # Calculate BASE
  derive_var_base(
    by_vars = vars(STUDYID, USUBJID, PARAMCD, BASETYPE),
    source_var = AVAL,
    new_var = BASE
  ) %>%

  # Calculate BASEC
  derive_var_base(
    by_vars = vars(STUDYID, USUBJID, PARAMCD, BASETYPE),
    source_var = AVALC,
    new_var = BASEC
  ) %>%

  # Calculate BNRIND
  derive_var_base(
    by_vars = vars(STUDYID, USUBJID, PARAMCD, BASETYPE),
    source_var = ANRIND,
    new_var = BNRIND
  ) %>%

  # Calculate CHG
  derive_var_chg() %>%

  # Calculate PCHG
  derive_var_pchg()


# Add RATIO derivations with new function??
# R2BASE, R2ANRLO R2ANRHI
adlb <- adlb %>%
  mutate(
    R2BASE = if_else(
      !is.na(AVAL) & !is.na(BASE) & BASE != 0, AVAL / BASE, NA_real_
    ),
    R2ANRLO = if_else(
      !is.na(AVAL) & !is.na(ANRLO) & ANRLO != 0, AVAL / ANRLO, NA_real_
    ),
    R2ANRHI = if_else(
      !is.na(AVAL) & !is.na(ANRHI) & ANRHI != 0, AVAL / ANRLO, NA_real_
    )
  )

# SHIFT derivation
adlb <- adlb %>%
  derive_var_shift(
    new_var = SHIFT1,
    from_var = BNRIND,
    to_var = ANRIND
    )

# ANL01FL: Flag last result within an AVISIT for post-baseline records
# LVOTFL: Flag last valid on-treatment record
adlb <- adlb %>%
  derive_var_extreme_flag(
    new_var = ANL01FL,
    by_vars = vars(USUBJID, PARAMCD, AVISIT),
    order = vars(ADT, AVAL),
    mode = "last",
    filter = !is.na(AVISITN) & ONTRTFL == "Y"
  ) %>%

  derive_var_extreme_flag(
    new_var = LVOTFL,
    by_vars = vars(USUBJID, PARAMCD),
    order = vars(ADT, AVAL),
    mode = "last",
    filter = ONTRTFL == "Y"
  )

# Get treatment information
adlb <- adlb %>%

  # Assign TRTA, TRTP
  mutate(
    TRTP = TRT01P,
    TRTA = TRT01A
  )

# Get ASEQ
adlb <- adlb %>%
  # Calculate ASEQ
  derive_var_obs_number(
    new_var = ASEQ,
    by_vars = vars(STUDYID, USUBJID),
    order = vars(PARAMCD, ADT, AVISITN, VISITNUM),
    check_type = "error"
  )

# Add all ADSL variables
adlb <- adlb %>%
  derive_vars_merged(
    dataset_add = select(adsl, !!!negate_vars(adsl_vars)),
    by_vars = vars(STUDYID, USUBJID)
  )

# Final Steps, Select final variables and Add labels
# This process will be based on your metadata, no example given for this reason
# ...

# ---- Save output ----

dir <- tempdir() # Change to whichever directory you want to save the dataset in
save(adlb, file = file.path(dir, "adlb.rda"), compress = "bzip2")
