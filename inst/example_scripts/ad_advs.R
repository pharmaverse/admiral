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
data ("adsl")

# Join ADSL
advs_adsl <- left_join(vs,
                       select(adsl, -DOMAIN),
                       by = c("STUDYID", "USUBJID"))

# Calculate ADT
advs_dt <- derive_vars_dt(advs_adsl,
                          new_vars_prefix = "A",
                          dtc = VSDTC,
                          flag_imputation = FALSE)

# Calculate ADY
advs_ady <- derive_var_ady(advs_dt, reference_date = TRTSDT, date = ADT)

# Calculate AVAL, AVALC, AVALU
advs_aval <- mutate(advs_dy,
                    AVAL = VSSTRESN,
                    AVALC = VSSTRESC,
                    AVALU = VSSTRESU)

# Derive MAP (Mean Arterial Pressure from Systolic and Diastolic Pressure)
# This is an example of deriving a new record based on existing records.
# Note: this PARAMCD is not derived in the CDISC pilot and is presented
#       for demonstration purposes.
sysbp <- filter(advs_aval, VSTESTCD == "SYSBP") %>%
  select(-VSSEQ, -VSTEST, -VSORRES, -VSORRESU, -VSSTRESC, -VSSTRESN,
         -VSSTRESU, -VSLOC, -VSBLFL, -VSSTAT, -AVALC)

diabp <- filter(advs_aval, VSTESTCD == "DIABP") %>%
  select(STUDYID, USUBJID, VISITNUM, VSDTC, VSTPT, AVAL) %>%
  rename(DBPAVAL = AVAL)

advs_map <- left_join(sysbp, diabp,
                      by = c("STUDYID", "USUBJID", "VISITNUM", "VSDTC",
                             "VSDTC", "VSTPT")) %>%
  mutate(AVAL = ((2 * DBPAVAL) + AVAL) / 3,
         VSTESTCD = "MAP") %>%
  filter(!is.na(AVAL)) %>%
  union_all(advs_aval) %>%
  select(-DBPAVAL)

# A summary records for each Visit. This could be applicable if triplicates
# are collected and the average will be summmarized.
# This is an example of deriving a new summary record based on existing records.
# Note: This is not derived in the CDISC pilot and is presented
#       for demonstration purposes.
advs_avg <- derive_summary_records(
  advs_map,
  by_vars = vars(STUDYID, USUBJID, VSTESTCD, VISITNUM, ADT),
  fns = list(AVAL ~ mean),
  set_values_to = vars(DTYPE = "AVERAGE"))

# Assign PARAMCD, PARAM, and PARAMN
advs_paramcd <-
  mutate(advs_avg,
         PARAMCD = VSTESTCD,
         PARAM = case_when(VSTESTCD %in% c("SYSBP", "DIABP") ~  str_c(VSTEST, "(mmHg)", sep = " "),
                           VSTESTCD == "HEIGHT" ~ str_c(VSTEST, "(cm)", sep = " "),
                           VSTESTCD == "WEIGHT" ~ str_c(VSTEST, "(kg)", sep = " "),
                           VSTESTCD == "TEMP" ~ str_c(VSTEST, "(C)", sep = " "),
                           VSTESTCD == "PULSE" ~ str_c(VSTEST, "(beats/min)", sep = " "),
                           VSTESTCD == "MAP" ~ "Mean Arterial Pressure (mmHg)"),
         PARAMN = case_when(VSTESTCD == "SYSBP" ~ 1,
                            VSTESTCD == "DIABP" ~ 2,
                            VSTESTCD == "PULSE" ~ 3,
                            VSTESTCD == "WEIGHT" ~ 4,
                            VSTESTCD == "HEIGHT" ~ 5,
                            VSTESTCD == "TEMP" ~ 6,
                            VSTESTCD == "MAP" ~ 7)) %>%
  mutate(VSTESTCD = ifelse(VSTESTCD == "MAP", NA_character_, VSTESTCD))

# Alternative Assignment example
# param_lookup <- tibble::tribble(
#   ~VSTESTCD, ~PARAMCD, ~PARAM, ~PARAMN,
#   "SYSBP",   "SYSBP",  "Systolic Blood Pressure (mmHg)",  1,
#   "DIABP",   "DIABP",  "Diastolic Blood Pressure (mmHg)", 2,
#   "PULSE",   "PULSE",  "Pulse Rate (beats/min)",          3,
#   "WEIGHT",  "WEIGHT", "Weight (kg)",                     4,
#   "HEIGHT",  "HEIGHT", "Height (cm)",                     5,
#   "TEMP",    "TEMP",   "Temperature (C)",                 6,
#   "MAP",     "MAP",    "Mean Arterial Pressure (mmHg)",   7)
#
# advs_paramcd_alt <- left_join(advs_aval, param_lookup, by = "VSTESTCD") %>%
#   mutate(VSTESTCD = ifelse(VSTESTCD == "MAP", NA_character_, VSTESTCD))

# Derive Timing
advs_timing <- mutate(advs_paramcd,
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

# ANL01FL: Choose non-missing results within an AVISIT
advs_anl01fl <- derive_extreme_flag(
  advs_timing,
  new_var = ANL01FL,
  by_vars = vars(USUBJID, PARAMCD, AVISIT, ATPT, DTYPE),
  order = vars(ADT, AVAL),
  mode = "last",
  flag_filter = (!is.na(AVISITN)))

# Calculate ONTRTFL
# Note: ONTRTFL is not calculated in the CDISC pilot
advs_ontrtfl <- derive_var_ontrtfl(advs_anl01fl,
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

advs_anrind <- left_join(advs_ontrtfl, range_lookup, by = "PARAMCD") %>%
  derive_var_anrind()

# Calculate BASETYPE
advs_basetype <- derive_var_basetype(
  dataset = advs_anrind,
  basetypes = exprs("AFTER LYING DOWN FOR 5 MINUTES" = ATPTN == 815,
                    "AFTER STANDING FOR 1 MINUTE" = ATPTN == 816,
                    "AFTER STANDING FOR 3 MINUTES" = ATPTN == 817,
                    NA_character = is.na(ATPTN)))

# Calculate ABLFL
advs_ablfl <- derive_extreme_flag(
  advs_basetype,
  new_var = ABLFL,
  by_vars = vars(STUDYID, USUBJID, BASETYPE, PARAMCD),
  order = vars(ADT),
  mode = "last",
  flag_filter = (ANL01FL == "Y" & AVISIT == "Baseline"))

# Calculate BASE & BASEC
advs_base <- derive_var_base(
  advs_ablfl,
  by_vars = vars(STUDYID, USUBJID, PARAMCD, BASETYPE))

advs_basec <- derive_var_basec(
  advs_base,
  by_vars = vars(STUDYID, USUBJID, PARAMCD, BASETYPE))

# Calculate CHG, PCHG
advs_chg <- derive_var_chg(advs_basec)

advs_pchg <- derive_var_pchg(advs_chg)

# Assign TRTA, TRTP
advs_trta_trtp <- mutate(advs_pchg,
                         TRTP = TRT01P,
                         TRTA = TRT01A)

# Derive the Average AVAL per PARAMCD, AVISIT

# Create End of Treatment Record
advs_eot <-
  derive_extreme_flag(
    advs_trta_trtp,
    new_var = EOTFL,
    by_vars = vars(STUDYID, USUBJID, PARAMCD, ATPTN),
    order = vars(ADT),
    mode = "last",
    flag_filter = (4 < VISITNUM & VISITNUM <= 13 & ANL01FL == "Y")) %>%
    filter(EOTFL == "Y") %>%
  mutate(AVISIT = "End of Treatment",
         AVISITN = 99) %>%
  union_all(advs_trta_trtp) %>%
  select(-EOTFL)

# Calculate ASEQ
advs_aseq <- derive_obs_number(
  advs_eot,
  new_var = ASEQ,
  by_vars = vars(STUDYID, USUBJID),
  order = vars(PARAMCD, ADT, AVISITN, ATPTN))

# Derive AVALCATs
# Note: Derivation of AVALCAT is not represented in the CDISC Pilot. It is
#       presented for demonstration purposes.
avalcat_lookup <- tibble::tribble(
  ~PARAMCD, ~AVALCA1N, ~AVALCAT1,
  "HEIGHT", 1, ">100 cm",
  "HEIGHT", 2, "<= 100 cm")

advs_avalcat <- mutate(advs_aseq,
                       AVALCA1N = case_when(PARAMCD == "HEIGHT" & AVAL > 100 ~ 1,
                                            PARAMCD == "HEIGHT" & AVAL <= 100 ~ 2)) %>%
  left_join(avalcat_lookup, by = "PARAMCD")


# Final Steps, Select final variables and Add labels
# This process will be based on your metadata, no example given for this reason
advs <- advs_avalcat

save(advs, file = "data/advs.rda", compress = TRUE)
