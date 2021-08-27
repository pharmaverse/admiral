# Name: ADVS
#
# Label: Vital Signs Analysis Dataset
#
# Input: adsl, vs
#

library(dplyr)
library(lubridate)
library(stringr)
library(admiral)

# Use e.g. haven::read_sas to read in .sas7bdat, or other suitable functions
#  as needed and assign to the variables below.
adsl <- NULL
vs <- NULL


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
    param == "HEIGHT" & aval > 100 ~ 1,
    param == "HEIGHT" & aval <= 100 ~ 2
  )
}

# ---- Derivations ----

# Join with ADSL immediately to pick up TRTSDT only (for ADY derivation)
# Join with ADSL again after derived parameters
advs_init <- vs %>%
  left_join(
    adsl %>% select(STUDYID, USUBJID, TRTSDT),
    by = c("STUDYID", "USUBJID")
  ) %>%
  # Calculate ADT
  derive_vars_dt(
    new_vars_prefix = "A",
    dtc = VSDTC,
    flag_imputation = FALSE
  ) %>%
  # Calculate ADY
  derive_var_ady(reference_date = TRTSDT, date = ADT) %>%
  # Add PARAMCD and PARAM
  left_join(param_lookup, by = "VSTESTCD") %>%
  # Calculate AVAL, AVALC, AVALU
  mutate(
    AVAL = VSSTRESN,
    AVALC = VSSTRESC,
    AVALU = VSSTRESU
  )

# Derive new parameters based on existing records.
advs_param <- advs_init %>%
  # Derive Mean Arterial Pressure
  derive_param_map(
    unit_var = NULL,
    by_vars = vars(USUBJID, VISIT, ADT, ADY, VSTPT, VSTPTNUM),
    set_values_to = vars(PARAMCD = "MAP", AVALU = "mmHg")
  ) %>%
  # Derive Body Surface Area
  derive_param_bsa(
    unit_var = NULL,
    by_vars = vars(USUBJID, VISIT, ADT, ADY, VSTPT, VSTPTNUM),
    set_values_to = vars(PARAMCD = "BSA", AVALU = "m^2")
  ) %>%
  # Derive Body Surface Area
  derive_param_bmi(
    unit_var = NULL,
    by_vars = vars(USUBJID, VISIT, ADT, ADY, VSTPT, VSTPTNUM),
    set_values_to = vars(PARAMCD = "BMI", AVALU = "kg/m^2")
  ) %>%
  # Add parameter details for derived parameters
  select(filter(PARAMCD %in% c("MAP", "BSA")),
         -PARAM, -PARAMN) %>%
  left_join(select(param_lookup, -VSTESTCD), by = "PARAMCD") %>%
  union_all(advs_init,filter(!(PARAMCD %in% c("MAP", "BSA"))))

# join ADSL vars and get visit info
advs_visit <- advs_param %>%
  select(-DOMAIN, -TRTSDT)
  left_join(adsl,by = c("STUDYID", "USUBJID")) %>%
  # Derive Timing
  mutate(
    ATPTN = VSTPTNUM,
    ATPT = VSTPT,
    AVISIT = case_when(
      str_detect(VISIT, "SCREEN") |
        str_detect(VISIT, "UNSCHED") |
        str_detect(VISIT, "RETRIEVAL") |
        str_detect(VISIT, "AMBUL") ~ NA_character_,
      !is.na(VISIT) ~ str_to_title(VISIT),
      TRUE ~ NA_character_
    ),
    AVISITN = case_when(
      VISIT == "BASELINE" ~ 0,
      str_detect(VISIT, "WEEK") ~ as.numeric(str_sub(VISIT, start = 5)),
      TRUE ~ NA_real_
    )
  ) %>%
  # Derive a new a new record as a summary record (e.g. mean of the triplictaes at each time point)
  # Note: This is not derived in the CDISC pilot and is presented
  #       for demonstration purposes.
  derive_summary_records(
    by_vars = vars(STUDYID, USUBJID, PARAMCD, VISITNUM, ADT),
    fns = list(AVAL ~ mean),
    set_values_to = vars(DTYPE = "AVERAGE")
  )

# Derive flag derivations
advs_flags <- advs_visit %>%
  # ANL01FL: Flag last (and highest) results within an AVISIT and ATPT
  derive_extreme_flag(
    new_var = ANL01FL,
    by_vars = vars(USUBJID, PARAMCD, AVISIT, ATPT, DTYPE),
    order = vars(ADT, AVAL),
    mode = "last",
    flag_filter = (!is.na(AVISITN))
  ) %>%

  # Calculate ONTRTFL
  # Note: ONTRTFL is not calculated in the CDISC pilot
  derive_var_ontrtfl(
    date = ADT,
    ref_start_date = TRTSDT,
    ref_end_date = TRTEDT
  ) %>%

  # Calculate ANRIND
  # Note: ANRIND along with ANRLO and ANRHI are not included in CDISC pilot
  left_join(range_lookup, by = "PARAMCD") %>%
  derive_var_anrind()

# Derive baseline derivations
advs_base <- advs_flags %>%
  # Calculate BASETYPE
  derive_var_basetype(
    basetypes = exprs(
      "AFTER LYING DOWN FOR 5 MINUTES" = ATPTN == 815,
      "AFTER STANDING FOR 1 MINUTE" = ATPTN == 816,
      "AFTER STANDING FOR 3 MINUTES" = ATPTN == 817
    )
  ) %>%
  # Calculate ABLFL
  derive_extreme_flag(
    new_var = ABLFL,
    by_vars = vars(STUDYID, USUBJID, BASETYPE, PARAMCD),
    order = vars(ADT, VSSEQ),
    mode = "last",
    flag_filter = (!is.na(AVAL) & ADT <= TRTSDT & !is.na(BASETYPE))
  ) %>%

  # Calculate BASE & BASEC
  derive_var_base(
    by_vars = vars(STUDYID, USUBJID, PARAMCD, BASETYPE)
  ) %>%
  derive_var_basec(
    by_vars = vars(STUDYID, USUBJID, PARAMCD, BASETYPE)
  ) %>%

  # Calculate CHG, PCHG
  derive_var_chg() %>%
  derive_var_pchg()

# Get treatment information
advs_trt <- advs_base %>%
  # Assign TRTA, TRTP
  mutate(
    TRTP = TRT01P,
    TRTA = TRT01A
  ) %>%

  # Create End of Treatment Record
  derive_extreme_flag(
    new_var = EOTFL,
    by_vars = vars(STUDYID, USUBJID, PARAMCD, ATPTN),
    order = vars(ADT),
    mode = "last",
    flag_filter = (4 < VISITNUM & VISITNUM <= 13 & ANL01FL == "Y")
  ) %>%
  filter(EOTFL == "Y") %>%
  mutate(
    AVISIT = "End of Treatment",
    AVISITN = 99
  ) %>%
  # union_all(advs) %>%
  select(-EOTFL)

# Get ASEQ and AVALCATx
advs_final <- advs_trt %>%
  # Calculate ASEQ
  derive_obs_number(
    new_var = ASEQ,
    by_vars = vars(STUDYID, USUBJID),
    order = vars(PARAMCD, ADT, AVISITN, VISITNUM, ATPTN, DTYPE),
    check_type = "warning"
  ) %>%

  # Derive AVALCATx
  # Note: Derivation of AVALCAT is not represented in the CDISC Pilot. It is
  #       presented for demonstration purposes.
  mutate(AVALCA1N = format_avalcat1n(param = PARAMCD, aval = AVAL)) %>%
  left_join(avalcat_lookup, by = "PARAMCD")


# Final Steps, Select final variables and Add labels
# This process will be based on your metadata, no example given for this reason
# ...

# ---- Save output ----

saveRDS(advs, file = "./ADVS.rds", compress = TRUE)
