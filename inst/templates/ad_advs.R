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
# https://pharmaverse.github.io/admiral/articles/admiral.html#handling-of-missing-values # nolint

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
# (https://pharmaverse.github.io/admiral/articles/visits_periods.html#visits)
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
# (https://pharmaverse.github.io/admiral/articles/visits_periods.html#treatment_bds)
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
