# Name: ADPP
#
# Label: Pharmacokinetics Parameters Analysis Dataset
#
# Description: Based on simulated data, create ADPP analysis dataset
#
# Input: pp, adsl
library(admiral)
library(pharmaversesdtm) # Contains example datasets from the CDISC pilot project
library(dplyr)
library(lubridate)
library(stringr)


# Load source datasets ----

# Use e.g. haven::read_sas to read in .sas7bdat, or other suitable functions
# as needed and assign to the variables below.
# For illustration purposes read in admiral test data

# Load PP and Adsl
pp <- pharmaversesdtm::pp
admiral_adsl <- admiral::admiral_adsl

# When SAS datasets are imported into R using haven::read_sas(), missing
# character values from SAS appear as "" characters in R, instead of appearing
# as NA values. Further details can be obtained via the following link:
# https://pharmaverse.github.io/admiral/cran-release/articles/admiral.html#handling-of-missing-values # nolint

pp <- convert_blanks_to_na(pp)

# Lookup tables ----
param_lookup <- tibble::tribble(
  ~PPTESTCD,  ~PARAMCD,                                     ~PARAM, ~PARAMN,
  "AUCALL",   "AUCALL",                                  "AUC All",       1,
  "AUCIFO",   "AUCIFO",                         "AUC Infinity Obs",       2,
  "AUCIFOD", "AUCIFOD",            "AUC Infinity Obs Norm by Dose",       3,
  "AUCINT",   "AUCINT",                        "AUC from T1 to T2",       4,
  "AUCLST",   "AUCLST",                 "AUC to Last Nonzero Conc",       5,
  "AUCPEO",   "AUCPEO",                   "AUC %Extrapolation Obs",       6,
  "CEND",       "CEND", "Concentration at the end of the infusion",       7, # non Cdisc Term
  "CLO",         "CLO",                             "Total CL Obs",       8,
  "CLST",       "CLST",                        "Last Nonzero Conc",       9,
  "CMAX",       "CMAX",                                 "Max Conc",      10,
  "CMAXD",     "CMAXD",                    "Max Conc Norm by Dose",      11,
  "CSF",         "CSF",                      "CSF to Plasma Ratio",      12, # non Cdisc Term
  "LAMZ",       "LAMZ",                                 "Lambda z",      13,
  "LAMZHL",   "LAMZHL",                       "Half-Life Lambda z",      14,
  "LAMZNPT", "LAMZNPT",            "Number of Points for Lambda z",      15,
  "R2ADJ",     "R2ADJ",                       "R Squared Adjusted",      16,
  "TCEND",     "TCEND",                             "Time of CEND",      17, # non Cdisc Term
  "TLST",       "TLST",                "Time of Last Nonzero Conc",      18,
  "TMAX",       "TMAX",                             "Time of CMAX",      19,
  "VSSO",       "VSSO",                "Vol Dist Steady State Obs",      20,
  "RCAMINT", "RCAMINT",                                       "Ae",      21,
  "RENALCL", "RENALCL",                                      "CLR",      22
)

# Assign AVALCATx
avalcax_lookup <- exprs(
  ~PARAMCD, ~condition,       ~AVALCAT1, ~AVALCA1N,
  "AUCALL",  AVAL < 19,           "<19",         1,
  "AUCALL", AVAL >= 19,          ">=19",         2
)

attr(param_lookup$PPTESTCD, "label") <- "Parameter Short Name"

# Derivations ----

# Get list of ADSL vars required for derivations
adsl_vars <- exprs(TRTSDT, TRTEDT, DTHDT, EOSDT, TRT01P, TRT01A)

adpp_pp <- pp %>%
  # Join ADSL with PP (need TRTSDT for ADY derivation)
  derive_vars_merged(
    dataset_add = admiral_adsl,
    new_vars = adsl_vars,
    by_vars = exprs(STUDYID, USUBJID)
  ) %>%
  ## Calculate ADT, ADY ----
  derive_vars_dt(
    new_vars_prefix = "A",
    dtc = PPRFDTC
  ) %>%
  derive_vars_dy(reference_date = TRTSDT, source_vars = exprs(ADT))

adpp_aval <- adpp_pp %>%
  ## Add PARAMCD only - add PARAM etc later ----
  left_join(
    select(param_lookup, PPTESTCD, PARAMCD),
    by = "PPTESTCD"
  ) %>%
  ## Calculate PARCAT1, AVAL and AVALC ----
  mutate(
    PARCAT1 = PPCAT,
    AVAL = PPSTRESN,
    AVALU = PPSTRESU,
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
# (https://pharmaverse.github.io/admiral/cran-release/articles/visits_periods.html#visit_bds)
adpp_avisit <- adpp_aval %>%
  # Derive Timing
  mutate(
    AVISITN = ADY,
    AVISIT = paste("Day", ADY),
    VISITNUM = AVISITN,
    VISIT = AVISIT
  ) %>%
  ## Assign TRTA, TRTP ----
  # See also the "Visit and Period Variables" vignette
  # (https://pharmaverse.github.io/admiral/cran-release/articles/visits_periods.html#treatment_bds)
  mutate(
    TRTP = TRT01P,
    TRTA = TRT01A
  ) %>%
  ## Derive AVALCA1N and AVALCAT1 ----
  derive_vars_cat(
    definition = avalcax_lookup,
    by_vars = exprs(PARAMCD)
  )
# Add all ADSL variables
adpp <- adpp_avisit %>%
  derive_vars_merged(
    dataset_add = select(admiral_adsl, !!!negate_vars(adsl_vars)),
    by_vars = exprs(STUDYID, USUBJID)
  )

# Final Steps, Select final variables and Add labels
# This process will be based on your metadata, no example given for this reason
# ...
# Save output ----

# Change to whichever directory you want to save the dataset in
dir <- tools::R_user_dir("admiral_templates_data", which = "cache")
if (!file.exists(dir)) {
  # Create the folder
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
}
save(adpp, file = file.path(dir, "adpp.rda"), compress = "bzip2")
