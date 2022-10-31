# Name: ADPP
#
# Label: Pharmacokinetics Parameters Analysis Dataset
#
# Description: Based on simulated data, create ADPP analysis dataset
#
# Input: pp, adsl
library(admiral)
library(admiral.test) # Contains example datasets from the CDISC pilot project
library(dplyr)
library(lubridate)
library(stringr)


# Load source datasets ----

# Use e.g. haven::read_sas to read in .sas7bdat, or other suitable functions
# as needed and assign to the variables below.
# For illustration purposes read in admiral test data

# Load PP and Adsl
data("admiral_pp")
data("admiral_adsl")

# When SAS datasets are imported into R using haven::read_sas(), missing
# character values from SAS appear as "" characters in R, instead of appearing
# as NA values. Further details can be obtained via the following link:
# https://pharmaverse.github.io/admiral/articles/admiral.html#handling-of-missing-values

pp <- convert_blanks_to_na(admiral_pp)

# Lookup tables ----
param_lookup <- tibble::tribble(
  ~PPTESTCD, ~PARAMCD, ~PARAM, ~PARAMN,
  "AUCALL", "AUCALL", "AUC All", 1,
  "AUCIFO", "AUCIFO", "AUC Infinity Obs", 2,
  "AUCIFOD", "AUCIFOD", "AUC Infinity Obs Norm by Dose", 3,
  "AUCINT", "AUCINT", "AUC from T1 to T2", 4,
  "AUCLST", "AUCLST", "AUC to Last Nonzero Conc", 5,
  "AUCPEO", "AUCPEO", "AUC %Extrapolation Obs", 6,
  "CEND", "CEND", "Concentration at the end of the infusion", 7, # non Cdisc Term
  "CLO", "CLO", "Total CL Obs", 8,
  "CLST", "CLST", "Last Nonzero Conc", 9,
  "CMAX", "CMAX", "Max Conc", 10,
  "CMAXD", "CMAXD", "Max Conc Norm by Dose", 11,
  "CSF", "CSF", "CSF to Plasma Ratio", 12, # non Cdisc Term
  "LAMZ", "LAMZ", "Lambda z", 13,
  "LAMZHL", "LAMZHL", "Half-Life Lambda z", 14,
  "LAMZNPT", "LAMZNPT", "Number of Points for Lambda z", 15,
  "R2ADJ", "R2ADJ", "R Squared Adjusted", 16,
  "TCEND", "TCEND", "Time of CEND", 17, # non Cdisc Term
  "TLST", "TLST", "Time of Last Nonzero Conc", 18,
  "TMAX", "TMAX", "Time of CMAX", 19,
  "VSSO", "VSSO", "Vol Dist Steady State Obs", 20
)

# ASSIGN AVALCAT1
avalcat_lookup <- tibble::tribble(
  ~PARAMCD, ~AVALCA1N, ~AVALCAT1,
  "AUCALL", 1, "< 19",
  "AUCALL", 2, ">= 19"
)

attr(param_lookup$PPTESTCD, "label") <- "Parameter Short Name"

# User defined functions ----

# Here are some examples of how you can create your own functions that
#  operates on vectors, which can be used in `mutate`.
format_avalcat1n <- function(param, aval) {
  case_when(
    param == "AUCALL" & aval < 19 ~ 1,
    param == "AUCALL" & aval >= 19 ~ 2,
    T ~ NA_real_
  )
}

# Derivations ----

# Get list of ADSL vars required for derivations
adsl_vars <- vars(TRTSDT, TRTEDT, DTHDT, EOSDT, TRT01P, TRT01A)

adpp <- pp %>%
  # Join ADSL with PP (need TRTSDT for ADY derivation)
  left_join(
    select(admiral_adsl, STUDYID, USUBJID, !!!adsl_vars),
    by = c("STUDYID", "USUBJID")
  ) %>%
  ## Calculate ADT, ADY ----
  derive_vars_dt(
    new_vars_prefix = "A",
    dtc = PPRFDTC
  ) %>%
  derive_vars_dy(reference_date = TRTSDT, source_vars = vars(ADT))

adpp <- adpp %>%
  ## Add PARAMCD only - add PARAM etc later ----
  left_join(
    select(param_lookup, PPTESTCD, PARAMCD),
    by = "PPTESTCD"
  ) %>%
  ## Calculate AVAL and AVALC ----
  mutate(
    AVAL = PPSTRESN,
    AVALC = PPSTRESC
  ) %>%
  # Remove variables
  select(-PPSTRESN, -PPSTRESC) %>%
  # Add ASEQ
  mutate(
    SRCDOM = DOMAIN,
    SRCVAR = "SEQ",
    SRCSEQ = PPSEQ
  ) %>%
  select(-DOMAIN, -PPSEQ)

## Get visit info ----
# See also the "Visit and Period Variables" vignette
# (https://pharmaverse.github.io/admiral/articles/visits_periods.html#visit_bds)
adpp <- adpp %>%
  # Derive Timing
  mutate(
    VISIT = "", # /!\ To remove
    VISITNUM = NA, # /!\ To remove
    AVISIT = case_when(
      str_detect(VISIT, "SCREEN|UNSCHED|RETRIEVAL|AMBUL") ~ NA_character_,
      !is.na(VISIT) ~ str_to_title(VISIT),
      TRUE ~ NA_character_
    ),
    AVISITN = VISITNUM
  ) %>%
  ## Assign TRTA, TRTP ----
  # See also the "Visit and Period Variables" vignette
  # (https://pharmaverse.github.io/admiral/articles/visits_periods.html#treatment_bds)
  mutate(
    TRTP = TRT01P,
    TRTA = TRT01A
  ) %>%
  ## Derive AVALCA1N and AVALCAT1 ----
  mutate(AVALCA1N = format_avalcat1n(param = PARAMCD, aval = AVAL)) %>%
  derive_vars_merged(dataset_add = avalcat_lookup, by_vars = vars(PARAMCD, AVALCA1N))

# Add all ADSL variables
adpp <- adpp %>%
  left_join(admiral_adsl,
    by = c("STUDYID", "USUBJID")
  )

# Final Steps, Select final variables and Add labels
# This process will be based on your metadata, no example given for this reason
# ...
# Save output ----

dir <- tempdir() # Change to whichever directory you want to save the dataset in
saveRDS(adpp, file = file.path(dir, "adpp.rds"), compress = "bzip2")
