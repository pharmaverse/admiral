# Name: ADVS
#
# Label: Vital Signs Analysis Dataset
#
# Input: adsl, vs
library(admiral)
library(dplyr)
library(lubridate)
library(stringr)

# ---- Load source datasets ----

# Use e.g. haven::read_sas to read in .sas7bdat, or other suitable functions
# as needed and assign to the variables below.
# For illustration purposes read in admiral test data

data("vs")
data("adsl")

vs <- convert_blanks_to_na(vs)

# The CDISC Pilot Data contains no SUPPVS data
# If you have a SUPPVS then uncomment function below

# vs <- derive_vars_suppqual(vs, suppvs) %>%


# ---- Lookup tables ----

# Assign PARAMCD, PARAM, and PARAMN
param_lookup <- tibble::tribble(
  ~VSTESTCD, ~PARAMCD, ~PARAM, ~PARAMN,
  "SYSBP", "SYSBP", "Systolic Blood Pressure (mmHg)", 1,
  "DIABP", "DIABP", "Diastolic Blood Pressure (mmHg)", 2,
  "PULSE", "PULSE", "Pulse Rate (beats/min)", 3,
  "WEIGHT", "WEIGHT", "Weight (kg)", 4,
  "HEIGHT", "HEIGHT", "Height (cm)", 5,
  "TEMP", "TEMP", "Temperature (C)", 6,
  "MAP", "MAP", "Mean Arterial Pressure (mmHg)", 7,
  "BMI", "BMI", "Body Mass Index(kg/m^2)", 8,
  "BSA", "BSA", "Body Surface Area(m^2)", 9
)
attr(param_lookup$VSTESTCD, "label") <- "Vital Signs Test Short Name"


# Assign ANRLO/HI, A1LO/HI
range_lookup <- tibble::tribble(
  ~PARAMCD, ~ANRLO, ~ANRHI, ~A1LO, ~A1HI,
  "SYSBP", 90, 130, 70, 140,
  "DIABP", 60, 80, 40, 90,
  "PULSE", 60, 100, 40, 110,
  "TEMP", 36.5, 37.5, 35, 38
)
# ASSIGN AVALCAT1
avalcat_lookup <- tibble::tribble(
  ~PARAMCD, ~AVALCA1N, ~AVALCAT1,
  "HEIGHT", 1, ">100 cm",
  "HEIGHT", 2, "<= 100 cm"
)

# ---- User defined functions ----

# Here are some examples of how you can create your own functions that
#  operates on vectors, which can be used in `mutate`.
format_avalcat1n <- function(param, aval) {
  case_when(
    param == "HEIGHT" & aval > 140 ~ 1,
    param == "HEIGHT" & aval <= 140 ~ 2
  )
}

# ---- Derivations ----

advs <- vs %>%

  # Join ADSL with VS (need TRTSDT for ADY derivation)
  left_join(
    select(adsl, STUDYID, USUBJID, TRTSDT),
    by = c("STUDYID", "USUBJID")
  ) %>%

  # Calculate ADT, ADY
  derive_vars_dt(
    new_vars_prefix = "A",
    dtc = VSDTC,
    flag_imputation = FALSE
  ) %>%

  derive_var_ady(reference_date = TRTSDT, date = ADT) %>%

  select(-TRTSDT)

advs <- advs %>%
  # Add PARAMCD only - add PARAM etc later
  left_join(
    select(param_lookup, VSTESTCD, PARAMCD),
    by = "VSTESTCD"
  ) %>%

  # Calculate AVAL and AVALC
  mutate(
    AVAL = VSSTRESN,
    AVALC = VSSTRESC
  ) %>%

  # Derive new parameters based on existing records.
  # Derive Mean Arterial Pressure
  derive_param_map(
    by_vars = vars(STUDYID, USUBJID, VISIT, VISITNUM, ADT, ADY, VSTPT, VSTPTNUM),
    set_values_to = vars(PARAMCD = "MAP"),
    get_unit_expr = VSSTRESU,
    filter = VSSTAT != "NOT DONE" | is.na(VSSTAT)
  ) %>%

  # Derive Body Surface Area
  derive_param_bsa(
    by_vars = vars(STUDYID, USUBJID, VISIT, VISITNUM, ADT, ADY, VSTPT, VSTPTNUM),
    method = "Mosteller",
    set_values_to = vars(PARAMCD = "BSA"),
    get_unit_expr = VSSTRESU,
    filter = VSSTAT != "NOT DONE" | is.na(VSSTAT)
  ) %>%

  # Derive Body Surface Area
  derive_param_bmi(
    by_vars = vars(STUDYID, USUBJID, VISIT, VISITNUM, ADT, ADY, VSTPT, VSTPTNUM),
    set_values_to = vars(PARAMCD = "BMI"),
    get_unit_expr = VSSTRESU,
    filter = VSSTAT != "NOT DONE" | is.na(VSSTAT)
  )


# get visit info
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

# Derive a new record as a summary record (e.g. mean of the triplicates at each time point)
advs <- advs %>%
  derive_summary_records(
    by_vars = vars(STUDYID, USUBJID, PARAMCD, AVISITN, AVISIT, ADT, ADY),
    filter = !is.na(AVAL),
    analysis_var = AVAL,
    summary_fun = mean,
    set_values_to = vars(DTYPE = "AVERAGE")
  )

advs <- advs %>%

  # Join ADSL with VS (need TRTSDT and TRTEDT for ONTRTFL derivation)
  left_join(
    select(adsl, STUDYID, USUBJID, TRTSDT, TRTEDT),
    by = c("STUDYID", "USUBJID")
  ) %>%

  # Calculate ONTRTFL
  derive_var_ontrtfl(
    start_date = ADT,
    ref_start_date = TRTSDT,
    ref_end_date = TRTEDT,
    filter_pre_timepoint = AVISIT == "Baseline"
  ) %>%
  select(-TRTEDT)

# Calculate ANRIND : requires the reference ranges ANRLO, ANRHI
# Also accommodates the ranges A1LO, A1HI
advs <- advs %>%
  left_join(range_lookup, by = "PARAMCD") %>%
  # Calculate ANRIND
  derive_var_anrind()

# Derive baseline flags
advs <- advs %>%
  # Calculate BASETYPE
  derive_var_basetype(
    basetypes = exprs(
      "LAST: AFTER LYING DOWN FOR 5 MINUTES" = ATPTN == 815,
      "LAST: AFTER STANDING FOR 1 MINUTE" = ATPTN == 816,
      "LAST: AFTER STANDING FOR 3 MINUTES" = ATPTN == 817,
      "LAST" = is.na(ATPTN)
    )
  ) %>%

  # Calculate ABLFL
  derive_extreme_flag(
    by_vars = vars(STUDYID, USUBJID, BASETYPE, PARAMCD),
    order = vars(ADT, VISITNUM, VSSEQ),
    new_var = ABLFL,
    mode = "last",
    filter = (!is.na(AVAL) & ADT <= TRTSDT & !is.na(BASETYPE) & is.na(DTYPE))
  ) %>%
  select(-TRTSDT)

# Derive baseline information
advs <- advs %>%

  # Calculate BASE
  derive_var_base(
    by_vars = vars(STUDYID, USUBJID, PARAMCD, BASETYPE)
  ) %>%

  # Calculate BASEC
  derive_var_basec(
    by_vars = vars(STUDYID, USUBJID, PARAMCD, BASETYPE)
  ) %>%

  # Calculate BNRIND
  derive_baseline(
    by_vars = vars(STUDYID, USUBJID, PARAMCD, BASETYPE),
    source_var = ANRIND,
    new_var = BNRIND
  ) %>%

  # Calculate CHG
  derive_var_chg() %>%

  # Calculate PCHG
  derive_var_pchg()


# ANL01FL: Flag last result within an AVISIT and ATPT for post-baseline records
advs <- advs %>%
  derive_extreme_flag(
    new_var = ANL01FL,
    by_vars = vars(USUBJID, PARAMCD, AVISIT, ATPT, DTYPE),
    order = vars(ADT, AVAL),
    mode = "last",
    filter = !is.na(AVISITN) & ONTRTFL == "Y"
  )

# Get treatment information
advs <- advs %>%

  # Join ADSL with VS (need TRT01P/TRT01A for TRTA/TRTP derivation)
  left_join(
    select(adsl, STUDYID, USUBJID, TRT01A, TRT01P),
    by = c("STUDYID", "USUBJID")
  ) %>%

  # Assign TRTA, TRTP
  mutate(
    TRTP = TRT01P,
    TRTA = TRT01A
  ) %>%

  # Create End of Treatment Record
  derive_extreme_flag(
    by_vars = vars(STUDYID, USUBJID, PARAMCD, ATPTN),
    order = vars(ADT),
    new_var = EOTFL,
    mode = "last",
    filter = (4 < VISITNUM & VISITNUM <= 13 & ANL01FL == "Y" & is.na(DTYPE))
  ) %>%
  filter(EOTFL == "Y") %>%
  mutate(
    AVISIT = "End of Treatment",
    AVISITN = 99
  ) %>%
  union_all(advs) %>%
  select(-EOTFL, -TRT01A, -TRT01P)

# Get ASEQ and AVALCATx and add PARAM/PARAMN
advs <- advs %>%
  # Calculate ASEQ
  derive_obs_number(
    new_var = ASEQ,
    by_vars = vars(STUDYID, USUBJID),
    order = vars(PARAMCD, ADT, AVISITN, VISITNUM, ATPTN, DTYPE),
    check_type = "error"
  ) %>%

  # Derive AVALCA1N and AVALCAT1
  mutate(AVALCA1N = format_avalcat1n(param = PARAMCD, aval = AVAL)) %>%
  left_join(avalcat_lookup, by = c("PARAMCD", "AVALCA1N")) %>%

  # Derive PARAM and PARAMN
  left_join(select(param_lookup, -VSTESTCD), by = "PARAMCD")


# Add all ADSL variables
advs <- advs %>%

  left_join(adsl, by = c("STUDYID", "USUBJID"))

# Final Steps, Select final variables and Add labels
# This process will be based on your metadata, no example given for this reason
# ...

# ---- Save output ----

save(advs, file = "data/advs.rda", compress = "bzip2")
