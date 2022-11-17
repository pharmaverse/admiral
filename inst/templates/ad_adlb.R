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

# Load source datasets ----

# Use e.g. haven::read_sas to read in .sas7bdat, or other suitable functions
# as needed and assign to the variables below.
# For illustration purposes read in admiral test data

data("admiral_lb")
data("admiral_adsl")

lb <- admiral_lb
adsl <- admiral_adsl

# When SAS datasets are imported into R using haven::read_sas(), missing
# character values from SAS appear as "" characters in R, instead of appearing
# as NA values. Further details can be obtained via the following link:
# https://pharmaverse.github.io/admiral/articles/admiral.html#handling-of-missing-values

lb <- convert_blanks_to_na(lb)

# Look-up tables ----

# Assign PARAMCD, PARAM, and PARAMN
param_lookup <- tibble::tribble(
  ~LBTESTCD, ~PARAMCD, ~PARAM, ~PARAMN,
  "ALB", "ALB", "Albumin (g/L)", 1,
  "ALP", "ALKPH", "Alkaline Phosphatase (U/L)", 2,
  "ALT", "ALT", "Alanine Aminotransferase (U/L)", 3,
  "ANISO", "ANISO", "Anisocytes", 4,
  "AST", "AST", "Aspartate Aminotransferase (U/L)", 5,
  "BASO", "BASO", "Basophils (10^9/L)", 6,
  "BASOLE", "BASOLE", "Basophils/Leukocytes (FRACTION)", 7,
  "BILI", "BILI", "Bilirubin (umol/L)", 8,
  "BUN", "BUN", "Blood Urea Nitrogen (mmol/L)", 9,
  "CA", "CA", "Calcium (mmol/L)", 10,
  "CHOL", "CHOLES", "Cholesterol (mmol/L)", 11,
  "CK", "CK", "Creatinine Kinase (U/L)", 12,
  "CL", "CL", "Chloride (mmol/L)", 13,
  "COLOR", "COLOR", "Color", 14,
  "CREAT", "CREAT", "Creatinine (umol/L)", 15,
  "EOS", "EOS", "Eosinophils (10^9/L)", 16,
  "EOSLE", "EOSLE", "Eosinophils/Leukocytes (FRACTION)", 17,
  "GGT", "GGT", "Gamma Glutamyl Transferase (U/L)", 18,
  "GLUC", "GLUC", "Glucose (mmol/L)", 19,
  "HBA1C", "HBA1C", "Hemoglobin A1C (1)", 20,
  "HCT", "HCT", "Hematocrit (1)", 21,
  "HGB", "HGB", "Hemoglobin (mmol/L)", 22,
  "K", "POTAS", "Potassium (mmol/L)", 23,
  "KETONES", "KETON", "Ketones", 24,
  "LYM", "LYMPH", "Lymphocytes (10^9/L)", 25,
  "LYMLE", "LYMPHLE", "Lymphocytes/Leukocytes (FRACTION)", 26,
  "MACROCY", "MACROC", "Macrocytes", 27,
  "MCH", "MCH", "Ery. Mean Corpuscular Hemoglobin (fmol(Fe))", 28,
  "MCHC", "MCHC", "Ery. Mean Corpuscular HGB Concentration (mmol/L)", 29,
  "MCV", "MCV", "Ery. Mean Corpuscular Volume (f/L)", 30,
  "MICROCY", "MICROC", "Microcytes", 31,
  "MONO", "MONO", "Monocytes (10^9/L)", 32,
  "MONOLE", "MONOLE", "Monocytes/Leukocytes (FRACTION)", 33,
  "PH", "PH", "pH", 34,
  "PHOS", "PHOS", "Phosphate (mmol/L)", 35,
  "PLAT", "PLAT", "Platelet (10^9/L)", 36,
  "POIKILO", "POIKIL", "Poikilocytes", 37,
  "POLYCHR", "POLYCH", "Polychromasia", 38,
  "PROT", "PROT", "Protein (g/L)", 39,
  "RBC", "RBC", "Erythrocytes (TI/L)", 40,
  "SODIUM", "SODIUM", "Sodium (mmol/L)", 41,
  "SPGRAV", "SPGRAV", "Specific Gravity", 42,
  "TSH", "TSH", "Thyrotropin (mU/L)", 43,
  "URATE", "URATE", "Urate (umol/L)", 44,
  "UROBIL", "UROBIL", "Urobilinogen", 45,
  "VITB12", "VITB12", "Vitamin B12 (pmol/L)", 46,
  "WBC", "WBC", "Leukocytes (10^9/L)", 47
)


# Derivations ----

# Get list of ADSL vars required for derivations
adsl_vars <- vars(TRTSDT, TRTEDT, TRT01A, TRT01P)

adlb <- lb %>%
  # Join ADSL with LB (need TRTSDT for ADY derivation)
  derive_vars_merged(
    dataset_add = adsl,
    new_vars = adsl_vars,
    by_vars = vars(STUDYID, USUBJID)
  ) %>%
  ## Calculate ADT, ADY ----
  derive_vars_dt(
    new_vars_prefix = "A",
    dtc = LBDTC
  ) %>%
  derive_vars_dy(reference_date = TRTSDT, source_vars = vars(ADT))

adlb <- adlb %>%
  ## Add PARAMCD PARAM and PARAMN - from LOOK-UP table ----
  # Replace with PARAMCD lookup function
  derive_vars_merged_lookup(
    dataset_add = param_lookup,
    new_vars = vars(PARAMCD, PARAM, PARAMN),
    by_vars = vars(LBTESTCD),
    check_type = "none",
    print_not_mapped = FALSE
  ) %>%
  ## Calculate PARCAT1 AVAL AVALC ANRLO ANRHI ----
  mutate(
    PARCAT1 = LBCAT,
    AVAL = LBSTRESN,
    AVALC = LBSTRESC,
    ANRLO = LBSTNRLO,
    ANRHI = LBSTNRHI
  )

# Derive Absolute values from fractional Differentials using WBC
# Only derive where absolute values do not already exist
# Need to populate ANRLO and ANRHI for newly created records
adlb <- adlb %>%
  # Derive absolute Basophils
  derive_param_wbc_abs(
    by_vars = vars(STUDYID, USUBJID, !!!adsl_vars, DOMAIN, VISIT, VISITNUM, ADT, ADY),
    set_values_to = vars(
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
    by_vars = vars(STUDYID, USUBJID, !!!adsl_vars, DOMAIN, VISIT, VISITNUM, ADT, ADY),
    set_values_to = vars(
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
# (https://pharmaverse.github.io/admiral/articles/visits_periods.html#visits)
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
      by_vars = vars(STUDYID, USUBJID, BASETYPE, PARAMCD),
      order = vars(ADT, VISITNUM, LBSEQ),
      new_var = ABLFL,
      mode = "last"
    ),
    filter = (!is.na(AVAL) & ADT <= TRTSDT & !is.na(BASETYPE))
  )

## Derive baseline information ----
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


## Calculate lab grading ----

# Assign ATOXDSCL and ATOXDSCH to hold lab grading terms
# ATOXDSCL and ATOXDSCH hold terms defined by NCI-CTCAEv4.
# See (https://pharmaverse.github.io/admiral/articles/lab_grading.html#implement_ctcv4)
grade_lookup <- tibble::tribble(
  ~PARAMCD, ~ATOXDSCL, ~ATOXDSCH,
  "ALB", "Hypoalbuminemia", NA_character_,
  "ALKPH", NA_character_, "Alkaline phosphatase increased",
  "ALT", NA_character_, "Alanine aminotransferase increased",
  "AST", NA_character_, "Aspartate aminotransferase increased",
  "BILI", NA_character_, "Blood bilirubin increased",
  "CA", "Hypocalcemia", "Hypercalcemia",
  "CHOLES", NA_character_, "Cholesterol high",
  "CK", NA_character_, "CPK increased",
  "CREAT", NA_character_, "Creatinine increased",
  "GGT", NA_character_, "GGT increased",
  "GLUC", "Hypoglycemia", "Hyperglycemia",
  "HGB", "Anemia", "Hemoglobin increased",
  "POTAS", "Hypokalemia", "Hyperkalemia",
  "LYMPH", "CD4 lymphocytes decreased", NA_character_,
  "PHOS", "Hypophosphatemia", NA_character_,
  "PLAT", "Platelet count decreased", NA_character_,
  "SODIUM", "Hyponatremia", "Hypernatremia",
  "WBC", "White blood cell decreased", "Leukocytosis",
)

# Assign grade criteria
# metadata atoxgr_criteria_ctcv4 used to implement NCI-CTCAEv4
# user could change to atoxgr_criteria_ctcv5 to implement NCI-CTCAEv5
# Note: Hyperglycemia and Hypophosphatemia not defined in NCI-CTCAEv5 so
# user would need to amend look-up table grade_lookup
# See (https://pharmaverse.github.io/admiral/articles/lab_grading.html#implement_ctcv5)
grade_crit <- atoxgr_criteria_ctcv4


# Add ATOXDSCL and ATOXDSCH
adlb <- adlb %>%
  derive_vars_merged(
    dataset_add = grade_lookup,
    by_vars = vars(PARAMCD)
  ) %>%
  # Derive toxicity grade for low values ATOXGRL

  derive_var_atoxgr_dir(
    meta_criteria = grade_crit,
    new_var = ATOXGRL,
    tox_description_var = ATOXDSCL,
    criteria_direction = "L",
    get_unit_expr = extract_unit(PARAM)
  ) %>%
  # Derive toxicity grade for low values ATOXGRH
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
    by_vars = vars(STUDYID, USUBJID, PARAMCD, BASETYPE),
    source_var = ATOXGRL,
    new_var = BTOXGRL
  ) %>%
  # Derive baseline toxicity grade for high values BTOXGRH
  derive_var_base(
    by_vars = vars(STUDYID, USUBJID, PARAMCD, BASETYPE),
    source_var = ATOXGRH,
    new_var = BTOXGRH
  ) %>%
  # Derive baseline toxicity grade for for overall grade BTOXGR
  derive_var_base(
    by_vars = vars(STUDYID, USUBJID, PARAMCD, BASETYPE),
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
      by_vars = vars(USUBJID, PARAMCD, AVISIT),
      order = vars(ADT, AVAL),
      new_var = ANL01FL,
      mode = "last"
    ),
    filter = !is.na(AVISITN) & ONTRTFL == "Y"
  ) %>%
  restrict_derivation(
    derivation = derive_var_extreme_flag,
    args = params(
      by_vars = vars(USUBJID, PARAMCD),
      order = vars(ADT, AVAL),
      new_var = LVOTFL,
      mode = "last"
    ),
    filter = ONTRTFL == "Y"
  )

## Get treatment information ----
# See also the "Visit and Period Variables" vignette
# (https://pharmaverse.github.io/admiral/articles/visits_periods.html#treatment_bds)
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
    by_vars = vars(STUDYID, USUBJID, PARAMCD, BASETYPE),
    order = vars(AVAL, ADT, AVISITN),
    mode = "first",
    # "AVISITN < 9997" to evaluate only real visits
    filter = (!is.na(AVAL) & ONTRTFL == "Y" & AVISITN < 9997),
    set_values_to = vars(
      AVISITN = 9997,
      AVISIT = "POST-BASELINE MINIMUM",
      DTYPE = "MINIMUM"
    )
  ) %>%
  # get MAXIMUM value
  derive_extreme_records(
    by_vars = vars(STUDYID, USUBJID, PARAMCD, BASETYPE),
    order = vars(desc(AVAL), ADT, AVISITN),
    mode = "first",
    # "AVISITN < 9997" to evaluate only real visits
    filter = (!is.na(AVAL) & ONTRTFL == "Y" & AVISITN < 9997),
    set_values_to = vars(
      AVISITN = 9998,
      AVISIT = "POST-BASELINE MAXIMUM",
      DTYPE = "MAXIMUM"
    )
  ) %>%
  # get LOV value
  derive_extreme_records(
    by_vars = vars(STUDYID, USUBJID, PARAMCD, BASETYPE),
    order = vars(ADT, AVISITN),
    mode = "last",
    # "AVISITN < 9997" to evaluate only real visits
    filter = (ONTRTFL == "Y" & AVISITN < 9997),
    set_values_to = vars(
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

# Save output ----

dir <- tempdir() # Change to whichever directory you want to save the dataset in
saveRDS(adlb, file = file.path(dir, "adlb.rds"), compress = "bzip2")
