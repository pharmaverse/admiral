# Name: ADEG
#
# Label: Electrocardiogram Analysis Dataset
#
# Description: Based on CDISC Pilot data, create ADEG analysis dataset
#
# Input: adsl, eg
library(admiral)
library(admiral.test) # Contains example datasets from the CDISC pilot project
library(dplyr)
library(lubridate)
library(stringr)

# Load source datasets ----

# Use e.g. `haven::read_sas()` to read in .sas7bdat, or other suitable functions
# as needed and assign to the variables below.
# For illustration purposes read in admiral test data

data("admiral_adsl")
data("admiral_eg")

adsl <- admiral_adsl
eg <- admiral_eg

# When SAS datasets are imported into R using haven::read_sas(), missing
# character values from SAS appear as "" characters in R, instead of appearing
# as NA values. Further details can be obtained via the following link:
# https://pharmaverse.github.io/admiral/articles/admiral.html#handling-of-missing-values

eg <- convert_blanks_to_na(eg)

# Lookup tables ----

# Assign PARAMCD, PARAM, and PARAMN
param_lookup <- tibble::tribble(
  ~EGTESTCD, ~PARAMCD, ~PARAM, ~PARAMN,
  "ECGINT", "EGINTP", "ECG Interpretation", 1,
  "HR", "HR", "Heart Rate (beats/min)", 2,
  "RR", "RR", "RR Duration (msec)", 3,
  "RRR", "RRR", "RR Duration Rederived (msec)", 4,
  "QT", "QT", "QT Duration (msec)", 10,
  "QTCBR", "QTCBR", "QTcB - Bazett's Correction Formula Rederived (msec)", 11,
  "QTCFR", "QTCFR", "QTcF - Fridericia's Correction Formula Rederived (msec)", 12,
  "QTLCR", "QTLCR", "QTlc - Sagie's Correction Formula Rederived (msec)", 13,
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

# ASSIGN AVALCAT1
avalcat_lookup <- tibble::tribble(
  ~AVALCA1N, ~AVALCAT1,
  1, "<= 450 msec",
  2, ">450<=480 msec",
  3, ">480<=500 msec",
  4, ">500 msec"
)

# ASSIGN CHGCAT1
chgcat_lookup <- tibble::tribble(
  ~CHGCAT1N, ~CHGCAT1,
  1, "<= 30 msec",
  2, ">30<=60 msec",
  3, ">60 msec"
)

# Here are some examples of how you can create your own functions that
#  operates on vectors, which can be used in `mutate()`. Info then used for
# lookup table
format_avalca1n <- function(paramcd, aval) {
  case_when(
    str_detect(paramcd, "QT") & aval <= 450 ~ 1,
    str_detect(paramcd, "QT") & aval > 450 & aval <= 480 ~ 2,
    str_detect(paramcd, "QT") & aval > 480 & aval <= 500 ~ 3,
    str_detect(paramcd, "QT") & aval > 500 ~ 4
  )
}

format_chgcat1n <- function(paramcd, chg) {
  case_when(
    str_detect(paramcd, "QT") & chg <= 30 ~ 1,
    str_detect(paramcd, "QT") & chg > 30 & chg <= 60 ~ 2,
    str_detect(paramcd, "QT") & chg > 60 ~ 3
  )
}


# Derivations ----

# Get list of ADSL vars required for derivations
adsl_vars <- vars(TRTSDT, TRTEDT, TRT01A, TRT01P)

adeg <- eg %>%
  # Join ADSL & EG (need TRTSDT for ADY derivation)
  derive_vars_merged(
    dataset_add = adsl,
    new_vars = adsl_vars,
    by_vars = vars(STUDYID, USUBJID)
  ) %>%
  ## Calculate ADTM, ADY ----
  derive_vars_dtm(
    new_vars_prefix = "A",
    dtc = EGDTC,
  ) %>%
  derive_vars_dy(reference_date = TRTSDT, source_vars = vars(ADTM))

adeg <- adeg %>%
  ## Add PARAMCD only (add PARAM, etc later) ----
  derive_vars_merged_lookup(
    dataset_add = param_lookup,
    new_vars = vars(PARAMCD),
    by_vars = vars(EGTESTCD)
  ) %>%
  ## Calculate AVAL and AVALC ----
  mutate(
    AVAL = EGSTRESN,
    AVALC = EGSTRESC
  ) %>%
  ## Derive new parameters based on existing records ----
  # Note that, for the following four `derive_param_*()` functions, only the
  # variables specified in `by_vars` will be populated in the newly created
  # records.

  # Derive RRR
  derive_param_rr(
    by_vars = vars(STUDYID, USUBJID, !!!adsl_vars, VISIT, VISITNUM, EGTPT, EGTPTNUM, ADTM, ADY),
    set_values_to = vars(PARAMCD = "RRR"),
    hr_code = "HR",
    get_unit_expr = tolower(EGSTRESU),
    filter = EGSTAT != "NOT DONE" | is.na(EGSTAT)
  ) %>%
  # Derive QTCBR
  derive_param_qtc(
    by_vars = vars(STUDYID, USUBJID, !!!adsl_vars, VISIT, VISITNUM, EGTPT, EGTPTNUM, ADTM, ADY),
    method = "Bazett",
    set_values_to = vars(PARAMCD = "QTCBR"),
    qt_code = "QT",
    rr_code = "RR",
    get_unit_expr = EGSTRESU,
    filter = EGSTAT != "NOT DONE" | is.na(EGSTAT)
  ) %>%
  # Derive QTCFR
  derive_param_qtc(
    by_vars = vars(STUDYID, USUBJID, !!!adsl_vars, VISIT, VISITNUM, EGTPT, EGTPTNUM, ADTM, ADY),
    method = "Fridericia",
    set_values_to = vars(PARAMCD = "QTCFR"),
    qt_code = "QT",
    rr_code = "RR",
    get_unit_expr = EGSTRESU,
    filter = EGSTAT != "NOT DONE" | is.na(EGSTAT)
  ) %>%
  # Derive QTLCR
  derive_param_qtc(
    by_vars = vars(STUDYID, USUBJID, !!!adsl_vars, VISIT, VISITNUM, EGTPT, EGTPTNUM, ADTM, ADY),
    method = "Sagie",
    set_values_to = vars(PARAMCD = "QTLCR"),
    qt_code = "QT",
    rr_code = "RR",
    get_unit_expr = EGSTRESU,
    filter = EGSTAT != "NOT DONE" | is.na(EGSTAT)
  )

## Get visit info ----
# See also the "Visit and Period Variables" vignette
# (https://pharmaverse.github.io/admiral/articles/visits_periods.html#visits)
adeg <- adeg %>%
  # Derive Timing
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
        str_detect(VISIT, "WEEK") ~ str_trim(str_replace(VISIT, "WEEK", "")),
        TRUE ~ NA_character_
      )
    ),
  )

## Derive a summary records representing the mean of the triplicates at each visit ----
# (if least 2 records available) for all parameter except EGINTP
adeg <- adeg %>%
  derive_summary_records(
    by_vars = vars(STUDYID, USUBJID, !!!adsl_vars, PARAMCD, AVISITN, AVISIT, ADT),
    analysis_var = AVAL,
    summary_fun = function(x) mean(x, na.rm = TRUE),
    filter = dplyr::n() >= 2 & PARAMCD != "EGINTP",
    set_values_to = vars(DTYPE = "AVERAGE")
  )

adeg <- adeg %>%
  ## Calculate ONTRTFL: from trt start up to 30 days after trt ends ----
  derive_var_ontrtfl(
    start_date = ADT,
    ref_start_date = TRTSDT,
    ref_end_date = TRTEDT,
    ref_end_window = 30,
    filter_pre_timepoint = AVISIT == "Baseline"
  )

## Calculate ANRIND: requires the reference ranges ANRLO, ANRHI ----
# Also accommodates the ranges A1LO, A1HI
adeg <- adeg %>%
  derive_vars_merged(
    dataset_add = range_lookup,
    by_vars = vars(PARAMCD)
  ) %>%
  # Calculate ANRIND
  derive_var_anrind()

## Derive baseline flags ----
adeg <- adeg %>%
  # Calculate BASETYPE
  derive_var_basetype(
    basetypes = rlang::exprs(
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
      by_vars = vars(STUDYID, USUBJID, BASETYPE, PARAMCD),
      order = vars(ADT, VISITNUM, EGSEQ),
      new_var = ABLFL,
      mode = "last"
    ),
    filter = ((!is.na(AVAL) | !is.na(AVALC)) &
      ADT <= TRTSDT & !is.na(BASETYPE) & is.na(DTYPE) &
      PARAMCD != "EGINTP"
    )
  )

## Derive baseline information ----
adeg <- adeg %>%
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

## ANL01FL: Flag last result within an AVISIT and ATPT for post-baseline records ----
adeg <- adeg %>%
  restrict_derivation(
    derivation = derive_var_extreme_flag,
    args = params(
      by_vars = vars(USUBJID, PARAMCD, AVISIT, ATPT, DTYPE),
      order = vars(ADT, AVAL),
      new_var = ANL01FL,
      mode = "last"
    ),
    filter = !is.na(AVISITN) & ONTRTFL == "Y"
  )

## Get treatment information ----
# See also the "Visit and Period Variables" vignette
# (https://pharmaverse.github.io/admiral/articles/visits_periods.html#treatment_bds)
adeg <- adeg %>%
  # Assign TRTA, TRTP
  mutate(TRTP = TRT01P, TRTA = TRT01A)

## Get ASEQ and AVALCAT1/CHGCAT1 and add PARAM/PARAMN ----
adeg <- adeg %>%
  # Calculate ASEQ
  derive_var_obs_number(
    new_var = ASEQ,
    by_vars = vars(STUDYID, USUBJID),
    order = vars(PARAMCD, ADT, AVISITN, VISITNUM, ATPTN, DTYPE),
    check_type = "error"
  ) %>%
  # Derive AVALCA1N and AVALCAT1
  mutate(AVALCA1N = format_avalca1n(param = PARAMCD, aval = AVAL)) %>%
  derive_vars_merged(
    dataset_add = avalcat_lookup,
    by_vars = vars(AVALCA1N)
  ) %>%
  # Derive CHGCAT1N and CHGCAT1
  mutate(CHGCAT1N = format_chgcat1n(param = PARAMCD, chg = CHG)) %>%
  derive_vars_merged(dataset_add = chgcat_lookup, by_vars = vars(CHGCAT1N)) %>%
  # Derive PARAM and PARAMN
  derive_vars_merged(
    dataset_add = select(param_lookup, -EGTESTCD),
    by_vars = vars(PARAMCD)
  )

# Add all ADSL variables
adeg <- adeg %>%
  derive_vars_merged(
    dataset_add = select(adsl, !!!negate_vars(adsl_vars)),
    by_vars = vars(STUDYID, USUBJID)
  )


# Save output ----

dir <- tempdir() # Change to whichever directory you want to save the dataset in
saveRDS(adeg, file = file.path(dir, "adeg.rds"), compress = "bzip2")
