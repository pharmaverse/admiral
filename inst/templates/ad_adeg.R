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
data("adsl")
# The CDISC Pilot Data contains no EG data
### Fake EG for demonstration
data("eg")

# Assign PARAMCD, PARAM, and PARAMN
param_lookup <- tibble::tribble(
  ~EGTESTCD, ~PARAMCD, ~PARAM, ~PARAMN,
  "EGINTP", "EGINTP", "ECG Interpretation", 1,
  "HR", "HR", "Heart Rate (beats/min)", 2,
  "RR", "RR", "RR Duration (msec)", 3,
  "QT", "QT", "QT Duration (msec)", 10,
)
range_lookup <- tibble::tribble(
  ~PARAMCD, ~ANRLO, ~ANRHI,
  "EGINTP", NA, NA,
  "HR", 40, 100,
  "RR", 600, 1500,
  "QT", 350, 450,
  "RRR", 600, 1500,
  "QTCBR", 350, 450,
  "QTCFR", 350, 450,
  "QTLCR", 350, 450,
)

# Start

# Join ADSL & EG
adeg <- eg %>%

  left_join(select(adsl, STUDYID, USUBJID, starts_with("TRT")),
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
  left_join(param_lookup, by = "EGTESTCD") %>%

  # Add required derived Parameters: QTcf, QTcB, RRd
  derive_param_rr(
    by_vars = vars(STUDYID, USUBJID, VISIT, VISITNUM, EGTPT, EGTPTNUM, ADTM),
    set_values_to = vars(PARAMCD = "RRR",
                         PARAM = "RR Duration Rederived (msec)",
                         PARAMN = 4
                         ),
    hr_code = "HR",
    filter = EGSTAT != "NOT DONE"
  ) %>%

  derive_param_qtcb(
    by_vars = vars(STUDYID, USUBJID, VISIT, VISITNUM, EGTPT, EGTPTNUM, ADTM),
    set_values_to = vars(PARAMCD = "QTCBR",
                         PARAM = "QTcB - Bazett's Correction Formula Rederived (msec)",
                         PARAMN = 11
                         ),
    qt_code = "QT",
    rr_code = "RR",
    filter = EGSTAT != "NOT DONE"
  ) %>%

  derive_param_qtcf(
    by_vars = vars(STUDYID, USUBJID, VISIT, VISITNUM, EGTPT, EGTPTNUM, ADTM),
    set_values_to = vars(PARAMCD = "QTCFR",
                         PARAM = "QTcF - Fridericia's Correction Formula Rederived (msec)",
                         PARAMN = 12
                         ),
    qt_code = "QT",
    rr_code = "RR",
    filter = EGSTAT != "NOT DONE"
  ) %>%

  derive_param_qtlc(
    by_vars = vars(STUDYID, USUBJID, VISIT, VISITNUM, EGTPT, EGTPTNUM, ADTM),
    set_values_to = vars(PARAMCD = "QTLCR",
                         PARAM = "QTlc - Sagie's Correction Formula Rederived (msec)",
                         PARAMN = 13
                         ),
    qt_code = "QT",
    rr_code = "RR",
    filter = EGSTAT != "NOT DONE"
  ) %>%

  # Derive Timing, Assign BASETYPE, TRTA/P
  mutate(
    ADT = date(ADTM),
    ATPTN = EGTPTNUM,
    ATPT = EGTPT,
    AVISIT = case_when(
      str_detect(VISIT, "SCREEN|UNSCHED|RETRIEVAL|AMBUL") ~ NA_character_,
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
    analysis_var = AVAL,
    summary_fun = function(x) mean(x, na.rm = TRUE),
    filter_rows = (dplyr::n() >= 2 & PARAMCD != "EGINTP"),
    set_values_to = vars(DTYPE = "AVERAGE")
  ) %>%

  # Calculate ABLFL, BASE & BASEC, CHG, PCHG
  derive_extreme_flag(
    by_vars = vars(STUDYID, USUBJID, BASETYPE, PARAMCD),
    order = vars(ADT, VISITNUM, ATPTN),
    new_var = ABLFL,
    mode = "last",
    filter = ((!is.na(AVAL) | !is.na(AVALC)) &
      ADT <= TRTSDT &
      (DTYPE != "AVERAGE" | PARAMCD != "EGINTP"))
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
    by_vars = vars(USUBJID, PARAMCD, AVISIT, ATPT, DTYPE),
    order = vars(ADT, AVAL),
    new_var = ANL01FL,
    mode = "last",
    filter = (!is.na(AVISITN))
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


saveRDS(adeg, file = "./ADEG.rds", compress = TRUE)
