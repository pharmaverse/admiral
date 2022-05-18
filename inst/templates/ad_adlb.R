# Name: ADLB
#
# Label: Lab Analysis Dataset
#
# Input: adsl, lb
library(admiral)
library(admiraltest) # Contains example datasets from the CDISC pilot project
library(dplyr)
library(lubridate)
library(stringr)

# ---- Load source datasets ----

# Use e.g. haven::read_sas to read in .sas7bdat, or other suitable functions
# as needed and assign to the variables below.
# For illustration purposes read in admiral test data

data("admiral_lb")
data("admiral_adsl")

lb <- admiral_lb
adsl <- admiral_adsl

lb <- convert_blanks_to_na(lb)

# ---- Look-up tables ----

# Assign PARAMCD, PARAM, and PARAMN
param_lookup <- tibble::tribble(
  ~LBTESTCD, ~PARAMCD, ~PARAM, ~PARAMN,
  "ALB", "ALB", "Albumin (g/L)", 1,
  "ALP", "ALKPH", "Alkaline Phosphatase (U/L)", 2,
  "ALT", "ALT", "Alanine Aminotransferase (U/L)", 3,
  "ANISO", "ANISO", "Anisocytes", 4,
  "AST", "AST", "Aspartate Aminotransferase (U/L)", 5,
  "BASO", "BASO", "Basophils (GI/L)", 6,
  "BILI", "BILI", "Bilirubin (umol/L)", 7,
  "BUN", "BUN", "Blood Urea Nitrogen (mmol/L)", 8,
  "CA", "CA", "Calcium (mmol/L)", 9,
  "CHOLES", "CHOLES", "Cholesterol (mmol/L)", 10,
  "CK", "CK", "Creatinine Kinase (U/L)", 11,
  "CL", "CL", "Chloride (mmol/L)", 12,
  "COLOR", "COLOR", "Color", 13,
  "CREAT", "CREAT", "Creatinine (umol/L)", 14,
  "EOS", "EOS", "Eosinophils (GI/L)", 15,
  "GGT", "GGT", "Gamma Glutamyl Transferase (U/L)", 16,
  "GLUC", "GLUC", "Glucose (mmol/L)", 17,
  "HBA1C", "HBA1C", "Hemoglobin A1C (1)", 18,
  "HCT", "HCT", "Hematocrit (1)", 19,
  "HGB", "HGB", "Hemoglobin (mmol/L)", 20,
  "K", "POTAS", "Potassium (mmol/L)", 21,
  "KETONES", "KETON", "Ketones", 22,
  "LYM", "LYMPH", "Lymphocytes (GI/L)", 23,
  "MACROCY", "MACROC", "Macrocytes", 24,
  "MCH", "MCH", "Ery. Mean Corpuscular Hemoglobin (fmol(Fe))", 25,
  "MCHC", "MCHC", "Ery. Mean Corpuscular HGB Concentration (mmol/L)", 26,
  "MCV", "MCV", "Ery. Mean Corpuscular Volume (f/L)", 27,
  "MICROCY", "MICROC", "Microcytes", 28,
  "MONO", "MONO", "Monocytes (GI/L)", 29,
  "PH", "PH", "pH", 30,
  "PHOS", "PHOS", "Phosphate (mmol/L)", 31,
  "PLAT", "PLAT", "Platelet (GI/L)", 32,
  "POIKILO", "POIKIL", "Poikilocytes", 33,
  "POLYCHR", "POLYCH", "Polychromasia", 34,
  "PROT", "PROT", "Protein (g/L)", 35,
  "RBC", "RBC", "Erythrocytes (TI/L)", 36,
  "SODIUM", "SODIUM", "Sodium (mmol/L)", 37,
  "SPGRAV", "SPGRAV", "Specific Gravity", 38,
  "TSH", "TSH", "Thyrotropin (mU/L)", 39,
  "URATE", "URATE", "Urate (umol/L)", 40,
  "UROBIL", "UROBIL", "Urobilinogen", 41,
  "VITB12", "VITB12", "Vitamin B12 (pmol/L)", 42,
  "WBC", "WBC", "Leukocytes (GI/L)", 43
)


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
    flag_imputation = "none"
  ) %>%
  derive_vars_dy(reference_date = TRTSDT, source_vars = vars(ADT))

adlb <- adlb %>%
  # Add PARAMCD PARAM and PARAMN - from LOOK-UP table
  # Replace with PARAMCD lookup function
  derive_vars_merged(
    dataset_add = param_lookup,
    new_vars = vars(PARAMCD, PARAM, PARAMN),
    by_vars = vars(LBTESTCD)
  ) %>%
  # Calculate PARCAT1 AVAL AVALC ANRLO ANRHI
  mutate(
    PARCAT1 = LBCAT,
    AVAL = LBSTRESN,
    AVALC = LBSTRESC,
    ANRLO = LBSTNRLO,
    ANRHI = LBSTNRHI
  )

# Get Visit Info
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
adlb <- adlb %>%
  derive_var_anrind()

# Derive baseline flags
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


# Calculate R2BASE, R2ANRLO and R2ANRHI
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

# Get treatment information
adlb <- adlb %>%
  # Assign TRTA, TRTP
  mutate(
    TRTP = TRT01P,
    TRTA = TRT01A
  )

# Get extreme values
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
