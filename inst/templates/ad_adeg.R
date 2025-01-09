# Name: ADEG
#
# Label: Electrocardiogram Analysis Dataset
#
# Description: Based on CDISC Pilot data, create ADEG analysis dataset
#
# Input: adsl, eg
library(admiral)
library(pharmaversesdtm) # Contains example datasets from the CDISC pilot project
library(dplyr)
library(lubridate)
library(stringr)

# Load source datasets ----

# Use e.g. `haven::read_sas()` to read in .sas7bdat, or other suitable functions
# as needed and assign to the variables below.
# For illustration purposes read in admiral test data

eg <- pharmaversesdtm::eg
adsl <- admiral::admiral_adsl

# When SAS datasets are imported into R using haven::read_sas(), missing
# character values from SAS appear as "" characters in R, instead of appearing
# as NA values. Further details can be obtained via the following link:
# https://pharmaverse.github.io/admiral/articles/admiral.html#handling-of-missing-values # nolint

eg <- convert_blanks_to_na(eg)

# Lookup tables ----

# Assign PARAMCD, PARAM, and PARAMN
param_lookup <- tibble::tribble(
  ~EGTESTCD, ~PARAMCD,                                                  ~PARAM, ~PARAMN,
  "ECGINT",  "EGINTP",                                    "ECG Interpretation",       1,
  "HR",          "HR",                                "Heart Rate (beats/min)",       2,
  "RR",          "RR",                                      "RR Duration (ms)",       3,
  "RRR",        "RRR",                            "RR Duration Rederived (ms)",       4,
  "QT",          "QT",                                      "QT Duration (ms)",      10,
  "QTCBR",    "QTCBR",     "QTcB - Bazett's Correction Formula Rederived (ms)",      11,
  "QTCFR",    "QTCFR", "QTcF - Fridericia's Correction Formula Rederived (ms)",      12,
  "QTLCR",    "QTLCR",      "QTlc - Sagie's Correction Formula Rederived (ms)",      13,
)

range_lookup <- tibble::tribble(
  ~PARAMCD, ~ANRLO, ~ANRHI,
  "EGINTP",     NA,     NA,
  "HR",         40,    100,
  "RR",        600,   1500,
  "QT",        350,    450,
  "RRR",       600,   1500,
  "QTCBR",     350,    450,
  "QTCFR",     350,    450,
  "QTLCR",     350,    450
)

# Assign AVALCAx
avalcax_lookup <- exprs(
  ~condition,                                                ~AVALCAT1, ~AVALCA1N,
  startsWith(PARAMCD, "QT") & AVAL <= 450,                 "<= 450 ms",         1,
  startsWith(PARAMCD, "QT") & AVAL > 450 & AVAL <= 480, ">450<=480 ms",         2,
  startsWith(PARAMCD, "QT") & AVAL > 480 & AVAL <= 500, ">480<=500 ms",         3,
  startsWith(PARAMCD, "QT") & AVAL > 500,                    ">500 ms",         4
)
# Assign CHGCAx
chgcax_lookup <- exprs(
  ~condition,                                           ~CHGCAT1, ~CHGCAT1N,
  startsWith(PARAMCD, "QT") & CHG <= 30,              "<= 30 ms",         1,
  startsWith(PARAMCD, "QT") & CHG > 30 & CHG <= 60, ">30<=60 ms",         2,
  startsWith(PARAMCD, "QT") & CHG > 60,                 ">60 ms",         3
)

# Derivations ----

# Get list of ADSL vars required for derivations
adsl_vars <- exprs(TRTSDT, TRTEDT, TRT01A, TRT01P)

adeg <- eg %>%
  # Join ADSL & EG (need TRTSDT for ADY derivation)
  derive_vars_merged(
    dataset_add = adsl,
    new_vars = adsl_vars,
    by_vars = exprs(STUDYID, USUBJID)
  ) %>%
  ## Calculate ADTM, ADY ----
  derive_vars_dtm(
    new_vars_prefix = "A",
    dtc = EGDTC,
  ) %>%
  derive_vars_dy(reference_date = TRTSDT, source_vars = exprs(ADTM))

adeg <- adeg %>%
  ## Add PARAMCD only (add PARAM, etc later) ----
  derive_vars_merged_lookup(
    dataset_add = param_lookup,
    new_vars = exprs(PARAMCD),
    by_vars = exprs(EGTESTCD)
  ) %>%
  ## Calculate AVAL and AVALC ----
  mutate(
    AVAL = EGSTRESN,
    AVALC = ifelse(
      is.na(EGSTRESN) | as.character(EGSTRESN) != EGSTRESC,
      EGSTRESC,
      NA
    )
  ) %>%
  ## Derive new parameters based on existing records ----
  # Note that, for the following four `derive_param_*()` functions, only the
  # variables specified in `by_vars` will be populated in the newly created
  # records.

  # Derive RRR
  derive_param_rr(
    by_vars = exprs(STUDYID, USUBJID, !!!adsl_vars, VISIT, VISITNUM, EGTPT, EGTPTNUM, ADTM, ADY),
    set_values_to = exprs(PARAMCD = "RRR"),
    hr_code = "HR",
    get_unit_expr = tolower(EGSTRESU),
    filter = EGSTAT != "NOT DONE" | is.na(EGSTAT)
  ) %>%
  # Derive QTCBR
  derive_param_qtc(
    by_vars = exprs(STUDYID, USUBJID, !!!adsl_vars, VISIT, VISITNUM, EGTPT, EGTPTNUM, ADTM, ADY),
    method = "Bazett",
    set_values_to = exprs(PARAMCD = "QTCBR"),
    qt_code = "QT",
    rr_code = "RR",
    get_unit_expr = EGSTRESU,
    filter = EGSTAT != "NOT DONE" | is.na(EGSTAT)
  ) %>%
  # Derive QTCFR
  derive_param_qtc(
    by_vars = exprs(STUDYID, USUBJID, !!!adsl_vars, VISIT, VISITNUM, EGTPT, EGTPTNUM, ADTM, ADY),
    method = "Fridericia",
    set_values_to = exprs(PARAMCD = "QTCFR"),
    qt_code = "QT",
    rr_code = "RR",
    get_unit_expr = EGSTRESU,
    filter = EGSTAT != "NOT DONE" | is.na(EGSTAT)
  ) %>%
  # Derive QTLCR
  derive_param_qtc(
    by_vars = exprs(STUDYID, USUBJID, !!!adsl_vars, VISIT, VISITNUM, EGTPT, EGTPTNUM, ADTM, ADY),
    method = "Sagie",
    set_values_to = exprs(PARAMCD = "QTLCR"),
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
    dataset_add = adeg,
    by_vars = exprs(STUDYID, USUBJID, !!!adsl_vars, PARAMCD, AVISITN, AVISIT, ADT, ADY),
    filter_add = dplyr::n() >= 2 & PARAMCD != "EGINTP",
    set_values_to = exprs(
      AVAL = mean(AVAL, na.rm = TRUE),
      DTYPE = "AVERAGE"
    )
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
    by_vars = exprs(PARAMCD)
  ) %>%
  # Calculate ANRIND
  derive_var_anrind()

## Derive baseline flags ----
adeg <- adeg %>%
  # Calculate BASETYPE
  derive_basetype_records(
    basetypes = exprs(
      "BASELINE DAY 1" = TRUE
    )
  ) %>%
  # Calculate ABLFL
  restrict_derivation(
    derivation = derive_var_extreme_flag,
    args = params(
      by_vars = exprs(STUDYID, USUBJID, BASETYPE, PARAMCD),
      order = exprs(ADT, VISITNUM, EGSEQ),
      new_var = ABLFL,
      mode = "last"
    ),
    filter = ((!is.na(AVAL) | !is.na(AVALC)) &
      ADT <= TRTSDT & !is.na(BASETYPE) & DTYPE == "AVERAGE" &
      PARAMCD != "EGINTP"
    )
  )

## Derive baseline information ----
adeg <- adeg %>%
  # Calculate BASE
  derive_var_base(
    by_vars = exprs(STUDYID, USUBJID, PARAMCD, BASETYPE),
    source_var = AVAL,
    new_var = BASE
  ) %>%
  # Calculate BASEC
  derive_var_base(
    by_vars = exprs(STUDYID, USUBJID, PARAMCD, BASETYPE),
    source_var = AVALC,
    new_var = BASEC
  ) %>%
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
adeg <- adeg %>%
  restrict_derivation(
    derivation = derive_var_extreme_flag,
    args = params(
      by_vars = exprs(USUBJID, PARAMCD, AVISIT, ATPT, DTYPE),
      order = exprs(ADT, AVAL),
      new_var = ANL01FL,
      mode = "last"
    ),
    filter = !is.na(AVISITN) & (ONTRTFL == "Y" | ABLFL == "Y") & DTYPE == "AVERAGE"
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
    by_vars = exprs(STUDYID, USUBJID),
    order = exprs(PARAMCD, ADT, AVISITN, VISITNUM, ATPTN, DTYPE),
    check_type = "error"
  ) %>%
  # Derive AVALCA1N and AVALCAT1
  derive_vars_cat(
    definition = avalcax_lookup
  ) %>%
  # Derive CHGCAT1N and CHGCAT1
  derive_vars_cat(
    definition = chgcax_lookup
  ) %>%
  # Derive PARAM and PARAMN
  derive_vars_merged(
    dataset_add = select(param_lookup, -EGTESTCD),
    by_vars = exprs(PARAMCD)
  )

# Add all ADSL variables
adeg <- adeg %>%
  derive_vars_merged(
    dataset_add = select(adsl, !!!negate_vars(adsl_vars)),
    by_vars = exprs(STUDYID, USUBJID)
  )


# Save output ----

# Change to whichever directory you want to save the dataset in
dir <- tools::R_user_dir("admiral_templates_data", which = "cache")
if (!file.exists(dir)) {
  # Create the folder
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
}
save(adeg, file = file.path(dir, "adeg.rda"), compress = "bzip2")
