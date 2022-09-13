# Name: ADPC
#
# Label: Pharmacokinetics Concentrations Analysis Dataset
#
# Description: Based on simulated data, create ADPC analysis dataset
#
# Input: pc, adsl
library(admiral)
library(dplyr)
library(lubridate)
library(stringr)
library(admiral.test) # Contains example datasets from the CDISC pilot project

# to remove
devtools::load_all()

# ---- Load source datasets ----

# Use e.g. haven::read_sas to read in .sas7bdat, or other suitable functions
# as needed and assign to the variables below.
# For illustration purposes read in admiral test data

# Load PC and Adsl
data("admiral_pc")
data("admiral_adsl")

# When SAS datasets are imported into R using haven::read_sas(), missing
# character values from SAS appear as "" characters in R, instead of appearing
# as NA values. Further details can be obtained via the following link:
# https://pharmaverse.github.io/admiral/articles/admiral.html#handling-of-missing-values

pc <- convert_blanks_to_na(admiral_pc) %>%
  rename(PCDTM = PCDTC) %>%
  mutate(PCDTC = as.character(as.Date(PCDTM)))


# ---- Lookup tables ----
param_lookup <- tibble::tribble(
  ~PCTESTCD, ~PARAMCD, ~PARAM, ~PARAMN,
  "XAN", "PKCONC", "Pharmacokinetic concentration of Xanomeline", 1,
)

# ASSIGN AVALCAT1
avalcat_lookup <- tibble::tribble(
  ~PARAMCD, ~AVALCA1N, ~AVALCAT1,
  "PKCONC", 1, "< 1",
  "PKCONC", 2, ">= 1"
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

# ---- Derivations ----

# Get list of ADSL vars required for derivations
adsl_vars <- vars(TRTSDT, TRTEDT, DTHDT, EOSDT, TRT01P, TRT01A)

adpc1 <- pc %>%
  # Join ADSL with PC (need TRTSDT for ADY derivation)
  left_join(
    select(admiral_adsl, STUDYID, USUBJID, !!!adsl_vars),
    by = c("STUDYID", "USUBJID")
  ) %>%
  # Calculate ADT, ADY
  derive_vars_dt(
    new_vars_prefix = "A",
    dtc = PCDTC
  ) %>%
  derive_vars_dy(reference_date = TRTSDT, source_vars = vars(ADT))

adpc2 <- adpc1 %>%
  # Add PARAMCD only - add PARAM etc later
  left_join(
    select(param_lookup, PCTESTCD, PARAMCD),
    by = "PCTESTCD"
  ) %>%
  # Calculate AVAL and AVALC
  mutate(
    AVAL = PCSTRESN,
    AVALC = PCSTRESC
  ) %>%
  # Remove variables
  select(-PCSTRESN, -PCSTRESC) %>%
  # Add ASEQ
  mutate(
    SRCDOM = DOMAIN,
    SRCVAR = "SEQ",
    SRCSEQ = PCSEQ
  ) %>%
  select(-DOMAIN, -PCSEQ)

# get visit info
adpc3 <- adpc2 %>%
  # Derive Timing
  mutate(
    # VISIT = "", # /!\ To remove
    # VISITNUM = NA, # /!\ To remove
    AVISIT = case_when(
      str_detect(VISIT, "SCREEN|UNSCHED|RETRIEVAL|AMBUL") ~ NA_character_,
      !is.na(VISIT) ~ str_to_title(VISIT),
      TRUE ~ NA_character_
    ),
    AVISITN = VISITNUM
  ) %>%
  # Assign TRTA, TRTP
  mutate(
    TRTP = TRT01P,
    TRTA = TRT01A
  ) %>%
  # Derive AVALCA1N and AVALCAT1
  mutate(AVALCA1N = format_avalcat1n(param = PARAMCD, aval = AVAL)) %>%
  derive_vars_merged(dataset_add = avalcat_lookup, by_vars = vars(PARAMCD, AVALCA1N))

# Add all ADSL variables
adpc <- adpc3 %>%
  left_join(select(admiral_adsl, !!!admiral:::negate_vars(adsl_vars)),
    by = c("STUDYID", "USUBJID")
  )

# Final Steps, Select final variables and Add labels
# This process will be based on your metadata, no example given for this reason
# ...
# ---- Save output ----

dir <- tempdir() # Change to whichever directory you want to save the dataset in
save(adpc, file = file.path(dir, "adpc.rda"), compress = "bzip2")
