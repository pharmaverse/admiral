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
#library(admiral)

devtools::load_all()

# Read in Data
data("adsl")
# The CDISC Pilot Data contains no EG data
### Fake EG from VS
data("eg")

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
  NA_character_, "QTLC", "QTlc Duration (msec)", 13,
)
range_lookup <- tibble::tribble(
  ~PARAMCD, ~ANRLO, ~ANRHI,
  "EGINTP", NA, NA,
  "HR", 40, 100,
  "RR", 600, 1500,
  "QT", 350, 450,
  "QTCB", 350, 450,
  "QTCF", 350, 450,
  "QTLC", 350, 450,
)

# Start

# Join ADSL & EG
adeg <- adsl %>%
  select(STUDYID, USUBJID, starts_with("TRT")) %>%
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

  # Add required derived Parameters: QTcf, QTcB, RRd
  derive_param_rr(
    filter = EGSTAT != "NOT DONE",
    new_param = "RRR",
    hr_code = "HR",
    by_vars = vars(STUDYID, USUBJID, VISIT, VISITNUM, EGTPT, EGTPTNUM, ADTM)
  ) %>%
  derive_param_qtcb(
    filter = EGSTAT != "NOT DONE",
    new_param = "QTCB",
    qt_code = "QT",
    rr_code = "RR",
    by_vars = vars(STUDYID, USUBJID, VISIT, VISITNUM, EGTPT, EGTPTNUM, ADTM),
    drop_values_from = vars(EGSEQ, EGBLFL, ends_with("RESU"))
  ) %>%
  derive_param_qtcf(
    filter = EGSTAT != "NOT DONE",
    new_param = "QTCF",
    qt_code = "QT",
    rr_code = "RR",
    by_vars = vars(STUDYID, USUBJID, VISIT, VISITNUM, EGTPT, EGTPTNUM, ADTM),
    drop_values_from = vars(EGSEQ, EGBLFL, ends_with("RESU"))
  ) %>%
  derive_param_qtlc(
    filter = EGSTAT != "NOT DONE",
    new_param = "QTLC",
    qt_code = "QT",
    rr_code = "RR",
    by_vars = vars(STUDYID, USUBJID, VISIT, VISITNUM, EGTPT, EGTPTNUM, ADTM),
    drop_values_from = vars(EGSEQ, EGBLFL, ends_with("RESU"))
  ) %>%

  # add PARAM/PARAMN
  left_join(select(param_lookup, -EGTESTCD), by = "PARAMCD") %>%

  # Derive Timing, Assign BASETYPE, TRTA/P
  mutate(
    ADT = date(ADTM),
    ATPTN = EGTPTNUM,
    ATPT = EGTPT,
    AVISIT = case_when(
      str_detect(VISIT, "SCREEN") | str_detect(VISIT, "UNSCHED") |
        str_detect(VISIT, "RETRIEVAL") | str_detect(VISIT, "AMBUL") ~ NA_character_,
      !is.na(VISIT) ~ str_to_title(VISIT),
      TRUE ~ NA_character_
    ),
    AVISITN = as.numeric(
      case_when(
        AVISIT == "Baseline" ~ "0",
        !is.na(AVISIT) ~ str_sub(AVISIT, start = 5),
        TRUE ~ NA_character_
      )
    ),
    BASETYPE = "LAST",
    TRTP = TRT01P,
    TRTA = TRT01A
  ) %>%

  # Derive a summary records representing the mean of the triplicates at each visit (if least 2
  # records available) for all parameter except EGINTP
  derive_summary_records(
    by_vars = vars(STUDYID, USUBJID, PARAMCD, VISITNUM, VISIT, ADT),
    fns = list(AVAL ~ mean(., na.rm = TRUE)),
    filter_rows = (dplyr::n() >= 2 & PARAMCD != "EGINTP"),
    set_values_to = vars(DTYPE = "AVERAGE"),
    drop_values_from = vars(EGBLFL, EGORRESU, EGSTRESU)
  ) %>%

  # Calculate ABLFL, BASE & BASEC, CHG, PCHG
  derive_extreme_flag(
    new_var = ABLFL,
    by_vars = vars(STUDYID, USUBJID, BASETYPE, PARAMCD),
    order = vars(ADT, VISITNUM, ATPTN),
    mode = "last",
    flag_filter = ((!is.na(AVAL) | !is.na(AVALC)) &
      ADT <= TRTSDT &
      (DTYPE == "AVERAGE" | PARAMCD == "EGINTP"))
  ) %>%
  derive_var_base(
    by_vars = vars(STUDYID, USUBJID, PARAMCD, BASETYPE)
  ) %>%
  derive_var_basec(
    by_vars = vars(STUDYID, USUBJID, PARAMCD, BASETYPE)
  ) %>%
  derive_var_chg() %>%
  derive_var_pchg() %>%

  # Calculate ONTRTFL: from trt start up to 30 days after trt ends.
  derive_var_ontrtfl(
    date = ADT,
    ref_start_date = TRTSDT,
    ref_end_date = TRTEDT,
    ref_end_window = 30
  ) %>%

  # ANL01FL: Flag last (and highest) results within an AVISIT and ATPT
  derive_extreme_flag(
    new_var = ANL01FL,
    by_vars = vars(USUBJID, PARAMCD, AVISIT, ATPT, DTYPE),
    order = vars(ADT, AVAL),
    mode = "last",
    flag_filter = (!is.na(AVISITN))
  ) %>%

  # Calculate ANRIND
  left_join(range_lookup, by = "PARAMCD") %>%
  derive_var_anrind() %>%

  #Derive AVALCTx, CHGCATx
  mutate(
    AVALCAT1 = case_when(
      str_detect(PARAMCD, "QT") & AVAL <= 450 ~ "<= 450 msec",
      str_detect(PARAMCD, "QT") & AVAL > 450 & AVAL <= 480 ~ ">450<=480 msec",
      str_detect(PARAMCD, "QT") & AVAL > 480 & AVAL <= 500 ~ ">480<=500 msec",
      str_detect(PARAMCD, "QT") & AVAL > 500 ~ ">500 msec"
    ),
    CHGCAT1 = case_when(
      str_detect(PARAMCD, "QT") & CHG <= 30 ~ "<= 30 msec",
      str_detect(PARAMCD, "QT") & CHG > 30 & CHG <= 60 ~ ">30<=60 msec",
      str_detect(PARAMCD, "QT") & CHG > 60 ~ ">60 msec"
    ),
  ) %>%

  # Calculate ASEQ
  derive_obs_number(
    new_var = ASEQ,
    by_vars = vars(STUDYID, USUBJID),
    order = vars(PARAMCD, ADT, AVISITN, VISITNUM, ATPTN, DTYPE),
    check_type = "warning"
  )


# save(adeg, file = "data/adeg.rda", compress = TRUE)
