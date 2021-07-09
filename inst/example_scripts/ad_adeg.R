# Name: ADEG
#
# Label: Electrocardiogram Analysis Dataset
#
# Description: Based on CDISC Pilot data, create ADEG analysis dataset
#
# Input: dm, eg
#

library(dplyr)
library(lubridate)
library(stringr)
library(admiral)

# Read in Data
# The CDISC Pilot Data contains no EG data
data("eg")
data ("adsl")
data("vs")
# test adsl
adsl <- adsl %>% filter(., as.character(USUBJID) %in% c("01-701-1015", "01-701-1028"))

### Fake EG
eg <- tibble::tribble(
  ~STUDYID,       ~USUBJID,  ~DOMAIN, ~EGSEQ, ~EGTESTCD, ~EGTEST,         ~EGORRES,  ~EGSTRESC,   ~EGORRESU,    ~EGDTC,               ~EGTPT,                       ~EGTPTNUM, ~VISIT,     ~VISITNUM,
  "CDISCPILOT01", "01-701-1015", "EG", 1, "EGINTP", "ECG Interpretation", "NORMAL",   "NORMAL",   "",          "2014-01-02T15:25:40", "AFTER LYING DOWN FOR 5 MINUTES", 1,     "BASELINE", 0,
  "CDISCPILOT01", "01-701-1015", "EG", 2, "HR",     "Heart Rate",         70.14,      "70.14",    "beats/min", "2014-01-02T15:25:40", "AFTER STANDING FOR 1 MINUTE",    2,     "BASELINE", 0,
  "CDISCPILOT01", "01-701-1015", "EG", 3, "EGINTP", "ECG Interpretation", "ABNORMAL", "ABNORMAL", "",          "2014-01-16T15:25:40", "AFTER LYING DOWN FOR 5 MINUTES", 1,     "WEEK 2",   2,
  "CDISCPILOT01", "01-701-1015", "EG", 4, "QT",     "QT Duration",        370,        "370",      "msec",      "2014-01-16T15:25:40", "AFTER STANDING FOR 3 MINUTES",   3,     "WEEK 2",   2,
  "CDISCPILOT01", "01-701-1015", "EG", 5, "HR",     "Heart Rate",         62.66,      "62.66",    "beats/min", "2014-01-09T15:25:40", "AFTER STANDING FOR 3 MINUTES",   3,     "WEEK 1",   1,
  "CDISCPILOT01", "01-701-1015", "EG", 6, "EGINTP", "ECG Interpretation", "NORMAL",   "NORMAL",   "",          "2014-01-16T15:25:40", "AFTER STANDING FOR 3 MINUTES",   3,     "WEEK 2",   2,
  "CDISCPILOT01", "01-701-1015", "EG", 7, "RR",     "RR Duration",        710,        "710",      "msec",      "2014-01-16T15:25:40", "AFTER LYING DOWN FOR 5 MINUTES", 1,     "WEEK 2",   2,
  "CDISCPILOT01", "01-701-1028", "EG", 1, "HR",     "Heart Rate",         85.45,      "85.45",    "beats/min", "2013-07-19T15:25:40", "AFTER LYING DOWN FOR 5 MINUTES", 1,     "BASELINE", 0,
  "CDISCPILOT01", "01-701-1028", "EG", 2, "QT",     "QT Duration",        480,        "480",      "msec",      "2013-08-03T15:26:40", "AFTER STANDING FOR 3 MINUTES",   3,     "WEEK 2",   2,
  "CDISCPILOT01", "01-701-1028", "EG", 5, "QT",     "QT Duration",        350,        "350",      "msec",      "2013-08-03T15:27:40", "AFTER STANDING FOR 3 MINUTES",   3,     "WEEK 2",   2,
  "CDISCPILOT01", "01-701-1028", "EG", 5, "QT",     "QT Duration",        410,        "410",      "msec",      "2013-08-03T15:25:40", "AFTER STANDING FOR 3 MINUTES",   3,     "WEEK 2",   2,
  "CDISCPILOT01", "01-701-1028", "EG", 3, "HR",     "Heart Rate",         56.54,      "56.54",    "beats/min", "2013-08-10T15:25:40", "AFTER LYING DOWN FOR 5 MINUTES", 1,     "WEEK 3",   3,
  "CDISCPILOT01", "01-701-1028", "EG", 4, "RR",     "RR Duration",        842,        "842",      "msec",      "2013-07-19T15:25:40", "AFTER STANDING FOR 1 MINUTE",    2,     "BASELINE", 0
)


# Join ADSL
adeg <- left_join(select(adsl, -DOMAIN),
                  eg,
                  by = c("STUDYID", "USUBJID"))

# Calculate ADT
adeg <- derive_vars_dt(adeg,
                       new_vars_prefix = "A",
                       dtc = EGDTC,
                       flag_imputation = FALSE)

# Calculate ADY
adeg <- derive_var_ady(adeg, reference_date = TRTSDT, date = ADT)

# Calculate AVAL, AVALC, AVALU
adeg <- mutate(adeg,
               AVAL = EGORRES,
               AVALC = EGSTRESC,
               AVALU = EGORRESU)

# Assign PARAMCD, PARAM, and PARAMN
param_lookup <- tibble::tribble(
  ~EGTESTCD, ~PARAMCD, ~PARAM,                  ~PARAMN,
  "EGINTP", "EGINTP", "ECG Interpretation",     1,
  "HR",     "HR",     "Heart Rate (beats/min)", 2,
  "QT",     "QT",     "QT Duration (msec)",     3,
  "RR",     "RR",     "RR Duration (msec)",     4
)


adeg <- left_join(adeg, param_lookup, by = "EGTESTCD")

# Derive Timing
adeg <- mutate(adeg,
               ATPTN = EGTPTNUM,
               ATPT = EGTPT,
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
adeg <- derive_summary_records(
  adeg,
  by_vars = vars(STUDYID, USUBJID, PARAMCD, VISITNUM, ADT),
  fns = list(AVAL ~ mean),
  set_values_to = vars(DTYPE = "AVERAGE"))

# ANL01FL: Flag last (and highest) results within an AVISIT and ATPT
adeg <- derive_extreme_flag(
  adeg,
  new_var = ANL01FL,
  by_vars = vars(USUBJID, PARAMCD, AVISIT, ATPT, DTYPE),
  order = vars(ADT, AVAL),
  mode = "last",
  flag_filter = (!is.na(AVISITN)))

# Calculate ONTRTFL
# Note: ONTRTFL is not calculated in the CDISC pilot
adeg <- derive_var_ontrtfl(adeg,
                           date = ADT,
                           ref_start_date = TRTSDT,
                           ref_end_date = TRTEDT)

# Calculate ANRIND
# Note: ANRIND along with ANRLO and ANRHI are not included in CDISC pilot
range_lookup <- tibble::tribble(
  ~PARAMCD, ~ANRLO, ~ANRHI,
  "EGINTP", NA,     NA,
  "HR",     40,     100,
  "QT",     200,    500,
  "RR",     600,    1500
)

adeg <- left_join(adeg, range_lookup, by = "PARAMCD") %>%
  derive_var_anrind()

# Calculate BASETYPE
adeg <- mutate(adeg, BASETYPE = "LAST")

# Calculate ABLFL
adeg <- derive_extreme_flag(
  adeg,
  new_var = ABLFL,
  by_vars = vars(STUDYID, USUBJID, BASETYPE, PARAMCD),
  order = vars(ADT),
  mode = "last",
  flag_filter = (!is.na(AVAL) & ADT <= TRTSDT & !is.na(BASETYPE)))

# Calculate BASE & BASEC
adeg <- derive_var_base(
  adeg,
  by_vars = vars(STUDYID, USUBJID, PARAMCD, BASETYPE))

adeg <- derive_var_basec(
  adeg,
  by_vars = vars(STUDYID, USUBJID, PARAMCD, BASETYPE))

# Assign TRTA, TRTP
adeg <- mutate(adeg,
               TRTP = TRT01P,
               TRTA = TRT01A)

# Create End of Treatment Record
adeg <-
  derive_extreme_flag(
    adeg,
    new_var = EOTFL,
    by_vars = vars(STUDYID, USUBJID, PARAMCD, ATPTN),
    order = vars(ADT),
    mode = "last",
    flag_filter = (4 < VISITNUM & VISITNUM <= 13 & ANL01FL == "Y")) %>%
  filter(EOTFL == "Y") %>%
  mutate(AVISIT = "End of Treatment",
         AVISITN = 99) %>%
  union_all(adeg) %>%
  select(-EOTFL)

# Calculate ASEQ
adeg <- derive_obs_number(
  adeg,
  new_var = ASEQ,
  by_vars = vars(STUDYID, USUBJID),
  order = vars(PARAMCD, ADT, AVISITN, VISITNUM, ATPTN, DTYPE),
  check_type = "warning")

# Derive AVALCATs
# Note: Derivation of AVALCAT is not represented in the CDISC Pilot. It is
#       presented for demonstration purposes.
adeg <- mutate(adeg,
               AVALCAT = case_when(PARAMCD == "HR" & AVAL > 70 ~ ">70 beats/min",
                                    PARAMCD == "HR" & AVAL <= 70 ~ "<= 70 beats/min"))

# Final Steps, Select final variables and Add labels
# This process will be based on your metadata, no example given for this reason
adeg <- adeg

save(adeg, file = "data/adeg.rda", compress = TRUE)
