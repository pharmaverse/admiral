# Name: ADVS
#
# Label: Vital Signs Analysis Dataset
#
# Description: Based on CDISC Pilot data, create ADVS analysis dataset
#
# Input: dm, vs
#

library(dplyr)
library(lubridate)
library(stringr)
library(admiral)

# Read in Data
# The CDISC Pilot Data contains no SUPPVS data
data("vs")
data("adsl")

# Join ADSL
advs <- left_join(vs,
                  adsl,
                  by = c("STUDYID", "USUBJID"))

# Calculate ADT
advs <- derive_vars_dt(advs,
                       new_vars_prefix = "A",
                       dtc = VSDTC,
                       flag_imputation = FALSE)

# Calculate ADY
advs <- derive_var_ady(advs, reference_date = TRTSDT, date = ADT)

# Calculate AVAL, AVALC, AVALU
advs <- mutate(advs,
               AVAL = VSSTRESN,
               AVALC = VSSTRESC,
               AVALU = VSSTRESU)

# Assign PARAMCD, PARAM, and PARAMN
param_lookup <- tibble::tribble(
  ~VSTESTCD, ~PARAMCD, ~PARAM, ~PARAMN,
  "SYSBP",   "SYSBP",  "Systolic Blood Pressure (mmHg)",  1,
  "DIABP",   "DIABP",  "Diastolic Blood Pressure (mmHg)", 2,
  "PULSE",   "PULSE",  "Pulse Rate (beats/min)",          3,
  "WEIGHT",  "WEIGHT", "Weight (kg)",                     4,
  "HEIGHT",  "HEIGHT", "Height (cm)",                     5,
  "TEMP",    "TEMP",   "Temperature (C)",                 6,
  NA,        "MAP",    "Mean Arterial Pressure (mmHg)",   7,
  NA,        "BMI",    "Body Mass Index (kg/m^2)",        8)

advs <- left_join(advs, select(param_lookup, VSTESTCD, PARAMCD), by = "VSTESTCD")

# Derive MAP (Mean Arterial Pressure from Systolic and Diastolic Pressure)
# This is an example of deriving a new record based on existing records.
# Note: this PARAMCD is not derived in the CDISC pilot and is presented
#       for demonstration purposes.
advs <- derive_param_map(
  advs,
  diabp_code = "DIABP",
  sysbp_code = "SYSBP",
  new_param = MAP,
  by_vars = vars(USUBJID, VISITNUM, VSDTC, VSTPT),
  set_values_to = vars(PARAMCD = "MAP")
)

# Derive BMI
advs <- derive_derived_param(
  advs,
  parameters = c("WEIGHT"),
  by_vars = vars(USUBJID, VISITNUM, VSDTC, VSTPT),
  constant_parameters = c("HEIGHT"),
  constant_by_vars = vars(USUBJID),
  analysis_value = AVAL.WEIGHT / (AVAL.HEIGHT / 100)^2,
  set_values_to = vars(PARAMCD = "BMI",
                       AVALU = "kg/m2")
)

# Assign PARAM and PARAMN
advs <- advs %>%
  left_join(select(param_lookup, -VSTESTCD), by = "PARAMCD")

# Derive Timing
advs <- mutate(advs,
               ATPTN = VSTPTNUM,
               ATPT = VSTPT,
               AVISIT = case_when(str_detect(VISIT, "SCREEN") |
                                    str_detect(VISIT, "UNSCHED") |
                                    str_detect(VISIT, "RETRIEVAL") |
                                    str_detect(VISIT, "AMBUL") ~ NA_character_,
                                    !is.na(VISIT) ~ str_to_title(VISIT),
                                  TRUE ~ NA_character_),
               AVISITN = case_when(VISIT == "BASELINE" ~ 0,
                                   str_detect(VISIT, "WEEK") ~
                                     as.numeric(str_sub(VISIT, start = 5)),
                                   TRUE ~ NA_real_))

# A summary records for each Visit. This could be applicable if triplicates
# are collected and the average will be summarized.
# This is an example of deriving a new summary record based on existing records.
# Note: This is not derived in the CDISC pilot and is presented
#       for demonstration purposes. Care should be taken to specify the
#       by_vars appropriately taking into consideration unscheduled visits and
#       the study's schedule of assessment.
advs <- derive_summary_records(
  advs,
  by_vars = vars(STUDYID, USUBJID, PARAMCD, VISITNUM, ADT),
  fns = list(AVAL ~ mean),
  set_values_to = vars(DTYPE = "AVERAGE"))

# ANL01FL: Flag last (and highest) results within an AVISIT and ATPT
advs <- derive_extreme_flag(
  advs,
  new_var = ANL01FL,
  by_vars = vars(USUBJID, PARAMCD, AVISIT, ATPT, DTYPE),
  order = vars(ADT, AVAL),
  mode = "last",
  flag_filter = (!is.na(AVISITN)))

# Calculate ONTRTFL
# Note: ONTRTFL is not calculated in the CDISC pilot
advs <- derive_var_ontrtfl(advs,
                           date = ADT,
                           ref_start_date = TRTSDT,
                           ref_end_date = TRTEDT)

# Calculate ANRIND
# Note: ANRIND along with ANRLO and ANRHI are not included in CDISC pilot
range_lookup <- tibble::tribble(
  ~PARAMCD, ~ANRLO, ~ANRHI, ~A1LO, ~A1HI,
  "SYSBP", 90, 130, 70, 140,
  "DIABP", 60, 80, 40, 90,
  "PULSE", 60, 100, 40, 110,
  "TEMP", 36.5, 37.5, 35, 38)

advs <- left_join(advs, range_lookup, by = "PARAMCD") %>%
  derive_var_anrind()

# Calculate BASETYPE
advs <- derive_var_basetype(
  dataset = advs,
  basetypes = exprs("AFTER LYING DOWN FOR 5 MINUTES" = ATPTN == 815,
                    "AFTER STANDING FOR 1 MINUTE" = ATPTN == 816,
                    "AFTER STANDING FOR 3 MINUTES" = ATPTN == 817))

# Calculate ABLFL
advs <- derive_extreme_flag(
  advs,
  new_var = ABLFL,
  by_vars = vars(STUDYID, USUBJID, BASETYPE, PARAMCD),
  order = vars(ADT, VISITNUM),
  mode = "last",
  flag_filter = (!is.na(AVAL) & ADT <= TRTSDT & !is.na(BASETYPE)))

# Calculate BASE & BASEC
advs <- derive_var_base(
  advs,
  by_vars = vars(STUDYID, USUBJID, PARAMCD, BASETYPE))

advs <- derive_var_basec(
  advs,
  by_vars = vars(STUDYID, USUBJID, PARAMCD, BASETYPE))

# Calculate CHG, PCHG
advs <- derive_var_chg(advs)

advs <- derive_var_pchg(advs)

# Assign TRTA, TRTP
advs <- mutate(advs,
               TRTP = TRT01P,
               TRTA = TRT01A)

# Create End of Treatment Record
advs <-
  derive_extreme_flag(
    advs,
    new_var = EOTFL,
    by_vars = vars(STUDYID, USUBJID, PARAMCD, ATPTN),
    order = vars(ADT),
    mode = "last",
    flag_filter = (4 < VISITNUM & VISITNUM <= 13 & ANL01FL == "Y")) %>%
    filter(EOTFL == "Y") %>%
  mutate(AVISIT = "End of Treatment",
         AVISITN = 99) %>%
  union_all(advs) %>%
  select(-EOTFL)

# Calculate ASEQ
advs <- derive_obs_number(
  advs,
  new_var = ASEQ,
  by_vars = vars(STUDYID, USUBJID),
  order = vars(PARAMCD, ADT, AVISITN, VISITNUM, ATPTN, DTYPE),
  check_type = "warning")

# Derive AVALCATs
# Note: Derivation of AVALCAT is not represented in the CDISC Pilot. It is
#       presented for demonstration purposes.
avalcat_lookup <- tibble::tribble(
  ~PARAMCD, ~AVALCA1N, ~AVALCAT1,
  "HEIGHT", 1, ">100 cm",
  "HEIGHT", 2, "<= 100 cm")

advs <- mutate(advs,
               AVALCA1N = case_when(PARAMCD == "HEIGHT" & AVAL > 100 ~ 1,
                                    PARAMCD == "HEIGHT" & AVAL <= 100 ~ 2)) %>%
  left_join(avalcat_lookup, by = c("PARAMCD", "AVALCA1N"))


# Final Steps, Select final variables and Add labels
# This process will be based on your metadata, no example given for this reason
advs <- advs

save(advs, file = "data/advs.rda", compress = TRUE)
