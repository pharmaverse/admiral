# Name: ADPC
#
# Label: Pharmacokinetics Concentrations Analysis Dataset
#
# Description: Based on simulated data, create ADPC analysis dataset
#
# Input: pc, ex,  adsl
library(admiral)
library(dplyr)
library(lubridate)
library(stringr)

library(admiral.test) # Contains example datasets from the CDISC pilot project or simulated

# ---- Load source datasets ----

# Use e.g. haven::read_sas to read in .sas7bdat, or other suitable functions
# as needed and assign to the variables below.
# For illustration purposes read in admiral test data


# Load PC, EX and ADSL
data("admiral_pc")
data("admiral_ex")
data("admiral_adsl")

# When SAS datasets are imported into R using haven::read_sas(), missing
# character values from SAS appear as "" characters in R, instead of appearing
# as NA values. Further details can be obtained via the following link:
# https://pharmaverse.github.io/admiral/articles/admiral.html#handling-of-missing-values

# Load EX

ex <- convert_blanks_to_na(admiral_ex)

# Load PC and fix format issue for PCDTC

pc <- convert_blanks_to_na(admiral_pc) %>%
  rename(PCDTM = PCDTC) %>%
  mutate(PCDTC = as.character.Date(PCDTM, format = "%Y-%m-%dT%H:%M:%S")) %>%
  select(-PCDTM)

# ---- Lookup tables ----
param_lookup <- tibble::tribble(
  ~PCTESTCD, ~PARAMCD, ~PARAM, ~PARAMN,
  "XAN", "XAN", "Pharmacokinetic concentration of Xanomeline", 1,
  "DOSE", "DOSE", "Xanomeline Patch Dose", 2,
)


attr(param_lookup$PCTESTCD, "label") <- "Pharmacokinetic Test Short Name"

# ---- User defined functions ----

# Here are some examples of how you can create your own functions that
#  operates on vectors, which can be used in `mutate`.
format_avalcat1n <- function(param, aval) {
  case_when(
    param == "PKCONC" & aval < 1 ~ 1,
    param == "PKCONC" & aval >= 1 ~ 2,
    T ~ NA_real_
  )
}

# Use this function to expand nominal times based on EXDOSFRQ
derive_var_expand_nfrlt <- function(dataset) {
  assert_data_frame(dataset, required_vars = vars(NFRLT))

  dataset <- dataset %>%
    mutate(orig_NFRLT = NFRLT)

  for (i in (seq_along(dataset$NFRLT))) {
    if (i > 1) {
      if (dataset$orig_EXDOSFRQ[i - 1] == "QD") {
        dose_int <- 24
      } else if (dataset$orig_EXDOSFRQ[i - 1] == "QOD") {
        dose_int <- 48
      } else if (dataset$orig_EXDOSFRQ[i - 1] == "BID") {
        dose_int <- 12
      }

      if (dataset$orig_NFRLT[i] == dataset$orig_NFRLT[i - 1] &
        dataset$USUBJID[i] == dataset$USUBJID[i - 1]) {
        dataset$NFRLT[i] <- dataset$NFRLT[i - 1] + dose_int
      } else {
        dataset$NFRLT[i] <- dataset$orig_NFRLT[i]
      }
    }
  }

  return(dataset)
}


# ---- Derivations ----

# Get list of ADSL vars required for derivations
adsl_vars <- vars(TRTSDT, TRTSDTM, TRTEDT, DTHDT, EOSDT, TRT01P, TRT01A)

adpc <- pc %>%
  # Join ADSL with PC (need TRTSDT for ADY derivation)
  left_join(
    select(admiral_adsl, STUDYID, USUBJID, !!!adsl_vars),
    by = c("STUDYID", "USUBJID")
  ) %>%
  # Calculate ADTM, ADT, ADY
  derive_vars_dtm(
    new_vars_prefix = "A",
    dtc = PCDTC,
    date_imputation = "FIRST",
    time_imputation = "00:00:00",
    ignore_seconds_flag = FALSE,
    flag_imputation = "auto"
  ) %>%
  # Derive dates and times from date/times
  derive_vars_dtm_to_dt(vars(ADTM)) %>%
  derive_vars_dtm_to_tm(vars(ADTM)) %>%
  derive_vars_dy(reference_date = TRTSDT, source_vars = vars(ADT)) %>%
  # Derive event ID and nominal relative time from first dose (NFRLT)
  mutate(
    EVID = 0,
    DRUG = PCTEST,
    NFRLT = if_else(PCTPTNUM < 0, 0, PCTPTNUM), .after = USUBJID
  )


# ---- Get dosing information ----

ex <- ex %>%
  left_join(select(admiral_adsl, STUDYID, USUBJID, !!!adsl_vars),
    by = c("STUDYID", "USUBJID")
  ) %>%
  # Keep records with nonzero dose
  filter(EXDOSE > 0) %>%
  # Add time and set missing end date to start date
  # Impute missing time to 00:00:00
  # Derive Analysis Start and End Dates
  derive_vars_dtm(
    new_vars_prefix = "AST",
    dtc = EXSTDTC,
    date_imputation = "FIRST",
    time_imputation = "00:00:00",
    ignore_seconds_flag = FALSE,
    flag_imputation = "auto"
  ) %>%
  derive_vars_dtm(
    new_vars_prefix = "AEN",
    dtc = EXENDTC,
    date_imputation = "LAST",
    time_imputation = "00:00:00",
    ignore_seconds_flag = FALSE,
    flag_imputation = "auto"
  ) %>%
  # Derive event ID and nominal relative time from first dose (NFRLT)
  mutate(
    EVID = 1,
    NFRLT = 24 * (VISITDY - 1), .after = USUBJID
  ) %>%
  # Set missing end dates to start date
  mutate(AENDTM = case_when(
    is.na(AENDTM) ~ ASTDTM,
    TRUE ~ AENDTM
  )) %>%
  # Derive dates from date/times
  derive_vars_dtm_to_dt(vars(ASTDTM)) %>%
  derive_vars_dtm_to_dt(vars(AENDTM))


# ---- Expand dosing records between start and end dates ----

ex_exp <- ex %>%
  mutate(orig_EXDOSFRQ = EXDOSFRQ) %>%
  create_single_dose_dataset(
    dose_freq = EXDOSFRQ,
    start_date = ASTDT,
    start_datetime = ASTDTM,
    end_date = AENDT,
    end_datetime = AENDTM,
    lookup_table = dose_freq_lookup,
    lookup_column = CDISC_VALUE,
    keep_source_vars = vars(
      STUDYID, USUBJID, EVID, EXDOSFRQ, orig_EXDOSFRQ, EXDOSFRM,
      NFRLT, EXDOSE, EXTRT, ASTDT, ASTDTM, AENDT, AENDTM,
      VISIT, VISITNUM, VISITDY,
      TRT01A, TRT01P, DOMAIN, EXSEQ, !!!adsl_vars
    )
  ) %>%
  mutate(
    ADTM = ASTDTM,
    DRUG = EXTRT
  ) %>%
  # Derive dates and times from datetimes
  derive_vars_dtm_to_dt(vars(ADTM)) %>%
  derive_vars_dtm_to_tm(vars(ADTM)) %>%
  derive_vars_dtm_to_tm(vars(ASTDTM)) %>%
  derive_vars_dtm_to_tm(vars(AENDTM)) %>%
  derive_vars_dy(reference_date = TRTSDT, source_vars = vars(ADT))

# ---- Add nominal time from first dose (NFRLT) to EX data ----

ex_expn <- ex_exp %>%
  derive_var_expand_nfrlt() %>%
  # Derive AVISIT based on nominal relative time
  mutate(
    AVISITN = NFRLT %/% 24 + 1,
    AVISIT = paste("Day", AVISITN)
  )

# ---- Find first dose per treatment per subject ----

timezero <- ex_expn %>%
  group_by(STUDYID, USUBJID, DRUG) %>%
  summarize(timezero = min(ADTM, na.rm = TRUE))

# ---- Join with ADPC data and keep only subjects with dosing ----

adpc <- left_join(timezero, adpc, by = c("STUDYID", "USUBJID", "DRUG")) %>%
  ungroup() %>%
  # Derive AVISIT from nominal relative time
  mutate(
    AVISITN = NFRLT %/% 24 + 1,
    AVISIT = paste("Day", AVISITN),
  )

# ---- Find last dose  ----
# Use derive_vars_joined for consistency with other variables
# This is equivalent to derive_vars_last_dose in this case

adpc <- adpc %>%
  derive_vars_joined(
    dataset_add = ex_expn,
    by_vars = vars(USUBJID),
    order = vars(ADTM),
    new_vars = vars(
      ADTM_prev = ADTM, EXDOSE_prev = EXDOSE, AVISIT_prev = AVISIT,
      AENDTM_prev = AENDTM
    ),
    join_vars = vars(ADTM),
    filter_add = NULL,
    filter_join = ADTM > ADTM.join,
    mode = "last",
    check_type = "none"
  )

# ---- Find next dose  ----

adpc <- adpc %>%
  derive_vars_joined(
    dataset_add = ex_expn,
    by_vars = vars(USUBJID),
    order = vars(ADTM),
    new_vars = vars(
      ADTM_next = ADTM, EXDOSE_next = EXDOSE, AVISIT_next = AVISIT,
      AENDTM_next = AENDTM
    ),
    join_vars = vars(ADTM),
    filter_add = NULL,
    filter_join = ADTM <= ADTM.join,
    mode = "first",
    check_type = "none"
  )

# ---- Find previous nominal time ----

adpc <- adpc %>%
  derive_vars_joined(
    dataset_add = ex_expn,
    by_vars = vars(USUBJID),
    order = vars(NFRLT),
    new_vars = vars(NFRLT_prev = NFRLT),
    join_vars = vars(NFRLT),
    filter_add = NULL,
    filter_join = NFRLT > NFRLT.join,
    mode = "last",
    check_type = "none"
  )

# ---- Find next nominal time ----

adpc <- adpc %>%
  derive_vars_joined(
    dataset_add = ex_expn,
    by_vars = vars(USUBJID),
    order = vars(NFRLT),
    new_vars = vars(NFRLT_next = NFRLT),
    join_vars = vars(NFRLT),
    filter_add = NULL,
    filter_join = NFRLT <= NFRLT.join,
    mode = "first",
    check_type = "none"
  )

# ---- Combine ADPC and EX data ----
# Derive Relative Time Variables


adpc <- bind_rows(adpc, ex_expn) %>%
  group_by(USUBJID) %>%
  mutate(
    timezero = min(timezero, na.rm = TRUE),
    FANLDTM = timezero,
    min_NFRLT = min(NFRLT_prev, na.rm = TRUE),
    maxdate = max(ADT[EVID == 0], na.rm = TRUE), .after = USUBJID
  ) %>%
  arrange(USUBJID, ADTM) %>%
  ungroup() %>%
  filter(ADT <= maxdate) %>%
  # Derive Actual Relative Time from First Dose (AFRLT)
  derive_vars_duration(
    new_var = AFRLT,
    start_date = FANLDTM,
    end_date = ADTM,
    out_unit = "hours",
    floor_in = FALSE,
    add_one = FALSE
  ) %>%
  # Derive Actual Relative Time from Reference Dose (ARRLT)
  derive_vars_duration(
    new_var = ARRLT,
    start_date = ADTM_prev,
    end_date = ADTM,
    out_unit = "hours",
    floor_in = FALSE,
    add_one = FALSE
  ) %>%
  # Derive Actual Relative Time from Next Dose (AXRLT not kept)
  derive_vars_duration(
    new_var = AXRLT,
    start_date = ADTM_next,
    end_date = ADTM,
    out_unit = "hours",
    floor_in = FALSE,
    add_one = FALSE
  ) %>%
  # Derive Nominal Relative Time from Reference Dose (NRRLT)
  mutate(
    NRRLT = case_when(
      EVID == 1 ~ 0,
      is.na(NFRLT_prev) ~ NFRLT - min_NFRLT,
      TRUE ~ NFRLT - NFRLT_prev
    ),
    NXRLT = case_when(
      EVID == 1 ~ 0,
      TRUE ~ NFRLT - NFRLT_next
    ),
    ARRLT = case_when(
      EVID == 1 ~ 0,
      is.na(ARRLT) ~ AXRLT,
      TRUE ~ ARRLT
    ),
    ATPTN = case_when(
      EVID == 1 ~ 0,
      TRUE ~ PCTPTNUM
    ),
    ATPT = case_when(
      EVID == 1 ~ "Dose",
      TRUE ~ PCTPT
    ),
    ATPTREF = case_when(
      EVID == 1 ~ AVISIT,
      is.na(AVISIT_prev) ~ AVISIT_next,
      TRUE ~ AVISIT_prev
    ),
    # Derive baseline flag for pre-dose records
    ABLFL = case_when(
      ATPT == "Pre-dose" ~ "Y",
      TRUE ~ NA_character_
    ),
    # Derive BASETYPE
    BASETYPE = paste(ATPTREF, "Baseline"),
    # Derive Reference Dose Date
    PCRFTDTM = case_when(
      EVID == 1 ~ ADTM,
      is.na(ADTM_prev) ~ ADTM_next,
      TRUE ~ ADTM_prev
    ),
    # Derive Actual Dose
    DOSEA = case_when(
      EVID == 1 ~ EXDOSE,
      is.na(EXDOSE_prev) ~ EXDOSE_next,
      TRUE ~ EXDOSE_next
    ),
    # Derive Planned Dose
    DOSEP = case_when(
      TRT01P == "Xanomeline High Dose" ~ 81,
      TRT01P == "Xanomeline Low Dose" ~ 54
    ),
    DOSEU = "mg",
  ) %>%
  # Derive dates and times from datetimes
  derive_vars_dtm_to_dt(vars(FANLDTM)) %>%
  derive_vars_dtm_to_tm(vars(FANLDTM)) %>%
  derive_vars_dtm_to_dt(vars(PCRFTDTM)) %>%
  derive_vars_dtm_to_tm(vars(PCRFTDTM)) %>%
  mutate(
    RFLTU = "h",
    RRLTU = "h",
    PARAMCD = coalesce(PCTESTCD, "DOSE"),
    ALLOQ = PCLLOQ,
    AVAL = case_when(
      EVID == 1 ~ EXDOSE,
      PCSTRESC == "<BLQ" & NFRLT == 0 ~ 0,
      PCSTRESC == "<BLQ" & NFRLT > 0 ~ 0.5 * ALLOQ,
      TRUE ~ PCSTRESN
    ),
    AVALC = case_when(
      EVID == 1 ~ paste(EXDOSE),
      TRUE ~ PCSTRESC
    ),
    AVALCAT1 = if_else(PCSTRESC == "<BLQ", PCSTRESC, prettyNum(signif(AVAL, digits = 3))),
  ) %>%
  # Add SRCSEQ
  mutate(
    SRCDOM = DOMAIN,
    SRCVAR = "SEQ",
    SRCSEQ = coalesce(PCSEQ, EXSEQ)
  )

# ---- Create DTYPE copy records ----

dtype <- adpc %>%
  filter(NFRLT > 0 & NXRLT == 0 & EVID == 0 & !is.na(AVISIT_next)) %>%
  select(-PCRFTDT, -PCRFTTM) %>%
  # Re-derive variables in for DTYPE copy records
  mutate(
    ABLFL = NA_character_,
    ATPTREF = AVISIT_next,
    ARRLT = AXRLT,
    NRRLT = NXRLT,
    PCRFTDTM = ADTM_next,
    DOSEA = EXDOSE_next,
    BASETYPE = paste(AVISIT_next, "Baseline"),
    ATPT = "Pre-dose",
    ATPTN = NFRLT,
    ABLFL = "Y",
    DTYPE = "COPY"
  ) %>%
  derive_vars_dtm_to_dt(vars(PCRFTDTM)) %>%
  derive_vars_dtm_to_tm(vars(PCRFTDTM))

# ---- Combine original records and DTYPE copy records ----

adpc <- bind_rows(adpc, dtype) %>%
  arrange(STUDYID, USUBJID, BASETYPE, ADTM, NFRLT) %>%
  mutate(
    # Derive MRRLT, ANL01FL and ANL02FL
    MRRLT = if_else(ARRLT < 0, 0, ARRLT),
    ANL01FL = "Y",
    ANL02FL = if_else(is.na(DTYPE), "Y", NA_character_),
  ) %>%
  # Derive BASE
  derive_var_base(
    by_vars = vars(STUDYID, USUBJID, PARAMCD, BASETYPE),
    source_var = AVAL,
    new_var = BASE,
    filter = ABLFL == "Y"
  )


# ---- Calculate Change from Baseline ----

adpc <- derive_var_chg(adpc)

# ---- Add ASEQ ----

adpc <- adpc %>%
  # Calculate ASEQ
  derive_var_obs_number(
    new_var = ASEQ,
    by_vars = vars(STUDYID, USUBJID),
    order = vars(ADTM, BASETYPE, EVID, AVISITN, ATPTN, DTYPE),
    check_type = "error"
  ) %>%
  # Remove temporary variables
  select(
    -DOMAIN, -PCSEQ, -timezero, -starts_with("orig"), -starts_with("min"),
    -starts_with("max"), -starts_with("EX"), -ends_with("next"),
    -ends_with("prev"), -DRUG, -EVID, -AXRLT, -NXRLT
  ) %>%
  # Derive PARAM and PARAMN
  derive_vars_merged(dataset_add = select(param_lookup, -PCTESTCD), by_vars = vars(PARAMCD))


# ---- Add all ADSL variables ----

adpc <- adpc %>%
  left_join(select(admiral_adsl, !!!admiral:::negate_vars(adsl_vars)),
    by = c("STUDYID", "USUBJID")
  )

# Final Steps, Select final variables and Add labels
# This process will be based on your metadata, no example given for this reason
# ...
# ---- Save output ----

dir <- tempdir() # Change to whichever directory you want to save the dataset in
saveRDS(adpc, file = file.path(dir, "adpc.rds"), compress = "bzip2")
