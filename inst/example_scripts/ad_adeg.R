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
# library(admiral)

devtools::load_all()

# Read in Data
data("adsl")
# The CDISC Pilot Data contains no EG data
# data("eg")
### Fake EG
data("vs")
eg <- vs %>%
  filter(VSTESTCD %in% c("PULSE", "DIABP", "SYSBP", "WEIGHT")) %>%
  set_names(~ str_to_upper(.) %>% str_replace_all("VS", "EG")) %>%
  mutate(
    EGTESTCD = case_when(
      EGTESTCD == "PULSE" ~ "HR",
      EGTESTCD == "SYSBP" ~ "RR",
      EGTESTCD == "DIABP" ~ "QT",
      EGTESTCD == "WEIGHT" ~ "EGINTP",
      TRUE ~ NA_character_
    ),
    EGTEST = case_when(
      EGTESTCD == "HR" ~ "Heart Rate",
      EGTESTCD == "RR" ~ "RR Duration",
      EGTESTCD == "QT" ~ "QT Duration",
      EGTESTCD == "EGINTP" ~ "ECG Interpretation",
      TRUE ~ NA_character_
    ),
    EGORRES = case_when(
      EGTESTCD == "RR" ~ as.character(EGSTRESN * 4),
      EGTESTCD == "QT" ~ as.character(EGSTRESN * 6),
      EGTESTCD == "EGINTP" & EGSTRESN > 85 ~ "ABNORMAL",
      EGTESTCD == "ECGINT" ~ "NORMAL",
      TRUE ~ EGORRES
    ),
    temp = if_else(EGTESTCD == "EGINTP", NA_character_, EGORRES),
    EGSTRESN = as.numeric(temp),
    EGSTRESC = EGORRES,
    EGSTRESU = case_when(
      EGTESTCD %in% c("RR", "QT") & EGSTAT != "NOT DONE" ~ "msec",
      EGTESTCD != "EGINTP" ~ EGSTRESU,
      TRUE ~ NA_character_
    )
  ) %>%
  select(-EGPOS)

# Assign PARAMCD, PARAM, and PARAMN
param_lookup <- tibble::tribble(
  ~EGTESTCD, ~PARAMCD, ~PARAM, ~PARAMN,
  "EGINTP", "EGINTP", "ECG Interpretation", 1,
  "HR", "HR", "Heart Rate (beats/min)", 2,
  "RR", "RR", "RR Duration (msec)", 3,
  NA_character_, "RRR", "RR Duration derived (msec)", 4,
  "QT", "QT", "QT Duration (msec)", 10,
  NA_character_, "QTCB", "QTcB Duration (msec)", 11,
  NA_character_, "QTCF", "QTcF Duration (msec)", 12,
  NA_character_, "QTCL", "QTlc Duration (msec)", 13,
)


# Start

# Join ADSL & EG
adeg <- adsl %>%
  left_join(eg,
    by = c("STUDYID", "USUBJID")
  ) %>%
  # Calculate ADT, ADY
  derive_vars_dtm(
    new_vars_prefix = "A",
    dtc = EGDTC,
    flag_imputation = "time"
  ) %>%
  derive_var_ady(reference_date = TRTSDT, date = ADTM) %>%
  # Calculate AVAL, AVALC, AVALU
  mutate(
    AVAL = EGSTRESN,
    AVALC = EGSTRESC,
    AVALU = tolower(EGSTRESU)
  ) %>%
  # Add PARAMCD to be able to use derive_param_rr...
  left_join(select(param_lookup, EGTESTCD, PARAMCD), by = "EGTESTCD") %>%
  ##############################################
  # Derived Parameters: QTcf, QTcB, RRd go there
  ##############################################
  # TEMORARY FIX to avoid
  # Error: It is expected that QT is measured in msec.
  # In the input dataset it is measured in `msec` and ``.
  filter(EGSTAT != "NOT DONE") %>%
  # END FIX
  derive_param_rr(
    filter = EGSTAT != "NOT DONE",
    new_param = "RRR",
    hr_code = "HR",
    by_vars = vars(STUDYID, USUBJID, VISIT, VISITNUM, EGTPT, EGTPTNUM, ADTM),
    # drop_values_from = vars(ends_with("RESU"))
  ) %>%
  derive_param_qtcb(
    filter = EGSTAT != "NOT DONE",
    new_param = "QTCB",
    qt_code = "QT",
    rr_code = "RR",
    by_vars = vars(STUDYID, USUBJID, VISIT, VISITNUM, EGTPT, EGTPTNUM, ADTM),
    # drop_values_from = vars(ends_with("RESU"))
  ) %>%
  derive_param_qtcf(
    filter = EGSTAT != "NOT DONE",
    new_param = "QTCF",
    qt_code = "QT",
    rr_code = "RR",
    by_vars = vars(STUDYID, USUBJID, VISIT, VISITNUM, EGTPT, EGTPTNUM, ADTM),
    # drop_values_from = vars(ends_with("RESU"))
  ) %>%
  left_join(select(param_lookup, -EGTESTCD), by = "PARAMCD") %>%
  # Derive Timing
  mutate(
    ADT=date(ADTM),
    ATPTN = EGTPTNUM,
    ATPT = EGTPT,
    AVISIT = case_when(
      str_detect(VISIT, "SCREEN") |str_detect(VISIT, "UNSCHED") |
        str_detect(VISIT, "RETRIEVAL") |str_detect(VISIT, "AMBUL") ~ NA_character_,
      !is.na(VISIT) ~ str_to_title(VISIT),
      TRUE ~ NA_character_
    ),
    #using the same code as advs : ie
    #AVISITN = case_when(VISIT == "BASELINE" ~ 0,
    #                    str_detect(VISIT, "WEEK") ~ as.numeric(str_sub(VISIT, start = 5)),
    #                    TRUE ~ NA_real_)
    #generates a warning: In eval_tidy(pair$rhs, env = default_env) : NAs introduced by coercion
    AVISITN = as.numeric(
      case_when(
        AVISIT == "Baseline" ~ "0",
        !is.na(AVISIT) ~ str_sub(AVISIT, start = 5),
        TRUE ~ NA_character_
      )
    ),
    BASETYPE="LAST"

  )%>%
  # Derive a summary records for each Visit. This could be applicable if triplicates
  # are collected and the average will be summarized.
  # This is an example of deriving a new summary record based on existing records.
  # We want to compute the mean of the triplicates at each visit, if we have at least 2 records
  # available
  derive_summary_records(
    by_vars = vars(STUDYID, USUBJID, PARAMCD, VISITNUM, VISIT, ADT),
    fns = list(AVAL ~ mean(., na.rm = TRUE)),
    filter_rows = (dplyr::n() >= 2 & PARAMCD !="EGINTP"),
    set_values_to = vars(DTYPE = "AVERAGE"),
    drop_values_from = vars(EGBLFL,EGORRESU, EGSTRESU)
  ) %>%
  # Calculate ABLFL
  derive_extreme_flag(
  new_var = ABLFL,
  by_vars = vars(STUDYID, USUBJID, BASETYPE, PARAMCD),
  order = vars(ADT,VISITNUM, ATPTN),
  mode = "last",
  flag_filter = (!is.na(AVAL) & ADT <= TRTSDT & !is.na(BASETYPE) & DTYPE=="AVERAGE")
)

# Calculate BASE & BASEC
adeg <- derive_var_base(
  adeg,
  by_vars = vars(STUDYID, USUBJID, PARAMCD, BASETYPE)
)

adeg <- derive_var_basec(
  adeg,
  by_vars = vars(STUDYID, USUBJID, PARAMCD, BASETYPE)
)

# ANL01FL: Flag last (and highest) results within an AVISIT and ATPT
advs <- derive_extreme_flag(
  advs,
  new_var = ANL01FL,
  by_vars = vars(USUBJID, PARAMCD, AVISIT, ATPT, DTYPE),
  order = vars(ADT, AVAL),
  mode = "last",
  flag_filter = (!is.na(AVISITN))
)

# Calculate ONTRTFL
# Note: ONTRTFL is not calculated in the CDISC pilot
advs <- derive_var_ontrtfl(advs,
                           date = ADT,
                           ref_start_date = TRTSDT,
                           ref_end_date = TRTEDT
)

# ANL01FL: Flag last (and highest) results within an AVISIT and ATPT
adeg <- derive_extreme_flag(
  adeg,
  new_var = ANL01FL,
  by_vars = vars(USUBJID, PARAMCD, AVISIT, ATPT, DTYPE),
  order = vars(ADT, AVAL),
  mode = "last",
  flag_filter = (!is.na(AVISITN))
)

# Calculate ONTRTFL
# Note: ONTRTFL is not calculated in the CDISC pilot
adeg <- derive_var_ontrtfl(adeg,
  date = ADT,
  ref_start_date = TRTSDT,
  ref_end_date = TRTEDT
)

# Calculate ANRIND
# Note: ANRIND along with ANRLO and ANRHI are not included in CDISC pilot
range_lookup <- tibble::tribble(
  ~PARAMCD, ~ANRLO, ~ANRHI,
  "EGINTP", NA, NA,
  "HR", 40, 100,
  "QT", 200, 500,
  "RR", 600, 1500
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
  flag_filter = (!is.na(AVAL) & ADT <= TRTSDT & !is.na(BASETYPE))
)

# Calculate BASE & BASEC
adeg <- derive_var_base(
  adeg,
  by_vars = vars(STUDYID, USUBJID, PARAMCD, BASETYPE)
)

adeg <- derive_var_basec(
  adeg,
  by_vars = vars(STUDYID, USUBJID, PARAMCD, BASETYPE)
)

# Assign TRTA, TRTP
adeg <- mutate(adeg,
  TRTP = TRT01P,
  TRTA = TRT01A
)

# Create End of Treatment Record
adeg <-
  derive_extreme_flag(
    adeg,
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
  union_all(adeg) %>%
  select(-EOTFL)

# Calculate ASEQ
adeg <- derive_obs_number(
  adeg,
  new_var = ASEQ,
  by_vars = vars(STUDYID, USUBJID),
  order = vars(PARAMCD, ADT, AVISITN, VISITNUM, ATPTN, DTYPE),
  check_type = "warning"
)

# Derive AVALCATs
# Note: Derivation of AVALCAT is not represented in the CDISC Pilot. It is
#       presented for demonstration purposes.
adeg <- mutate(adeg,
  AVALCAT = case_when(
    PARAMCD == "HR" & AVAL > 70 ~ ">70 beats/min",
    PARAMCD == "HR" & AVAL <= 70 ~ "<= 70 beats/min"
  )
)

# Final Steps, Select final variables and Add labels
# This process will be based on your metadata, no example given for this reason
adeg <- adeg

save(adeg, file = "data/adeg.rda", compress = TRUE)
