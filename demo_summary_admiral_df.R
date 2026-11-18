# Demo: `summary()` for `admiral_df` datasets
#
# Shows the summary for every ADaM class it tailors its report to:
#   - the subject classes  ADSL, BDS, OCCDS, TTE   (pharmaverseadam data)
#   - the device classes   ADDL, MDBDS, MDOCCDS, MDTTE   (built here, as there
#     is no example device data)
# and what it adds with `adsl =`, what its checks catch, and how it copes with a
# dataset still being derived.
#
# Run from the package root:
#   source("demo_summary_admiral_df.R", echo = TRUE)
#
# Not part of the package: do not commit (or add to .Rbuildignore).

devtools::load_all(quiet = TRUE)
library(tibble)
library(dplyr, warn.conflicts = FALSE)

adsl <- pharmaverseadam::adsl
advs <- pharmaverseadam::advs
adae <- pharmaverseadam::adae
adtte <- pharmaverseadam::adtte_onco

# 1. Subject-level classes ----------------------------------------------------

## ADSL: treatment arms, populations, subject flow, deaths, demographics
adsl_df <- admiral_df(adsl, dataset_keys = exprs(STUDYID, USUBJID), dataset_name = "ADSL")
summary(adsl_df)

## BDS: per parameter records, subjects, visits and AVAL; the visits in AVISITN
## order; records added by derivations (DTYPE); and the baseline check.
## These keys leave out the derived AVERAGE records and the unscheduled visits,
## so the structure check reports duplicates ...
advs_df <- admiral_df(
  advs,
  dataset_keys = exprs(STUDYID, USUBJID, PARAMCD, AVISIT, ATPT),
  dataset_name = "ADVS"
)
summary(advs_df)

## ... most of which go once the date and DTYPE are keys. The 39 left are
## repeated measurements on the same date in the pharmaverseadam data, which
## only ASEQ tells apart
advs_df <- admiral_df(
  advs_df,
  dataset_keys = exprs(STUDYID, USUBJID, PARAMCD, AVISIT, ATPT, ADT, DTYPE)
)
summary(advs_df)

## OCCDS: distinct terms, treatment-emergent, severity, seriousness, and the
## occurrence flag check
adae_df <- admiral_df(adae, dataset_keys = exprs(STUDYID, USUBJID, AESEQ), dataset_name = "ADAE")
summary(adae_df)

## TTE: events and censored per parameter. No keys set, so the structure is
## not checked, and the name is taken from the variable
adtte_df <- admiral_df(adtte)
summary(adtte_df)

# 2. Comparing with ADSL ------------------------------------------------------

## Subjects with records out of ADSL (overall, safety population, per arm),
## and that every subject with records is in ADSL
summary(advs_df, adsl = adsl)

## For OCCDS, also the treatment-emergent subjects out of the safety population
summary(adae_df, adsl = adsl)

## Subjects missing from ADSL are listed (the first and last few)
adae_orphans <- adae %>%
  mutate(USUBJID = if_else(row_number() <= 20, paste0(USUBJID, "-X"), USUBJID)) %>%
  admiral_df(dataset_name = "ADAE")
summary(adae_orphans, adsl = adsl)

# 3. Medical device classes ---------------------------------------------------
# Counted and checked per device (SPDEVID) rather than per subject. USUBJID is
# optional: where it is present, subjects are counted too and `adsl` is
# matched by it.

## ADDL: one record per device. Devices and subjects only
addl <- tribble(
  ~STUDYID, ~USUBJID,      ~SPDEVID, ~DEVTYPE,
  "DEV01",  "01-701-1015", "D-001",  "INSULIN PUMP",
  "DEV01",  "01-701-1015", "D-002",  "CGM SENSOR",
  "DEV01",  "01-701-1023", "D-003",  "INSULIN PUMP",
  "DEV01",  "01-701-1028", "D-004",  "INSULIN PUMP"
) %>%
  admiral_df(dataset_keys = exprs(STUDYID, SPDEVID), dataset_name = "ADDL")
get_dataset_class(addl) # inferred
summary(addl, adsl = adsl)

## MDBDS: device findings; the per-parameter table counts devices, and the
## baseline check is per device and parameter. No USUBJID here
mdbds <- tribble(
  ~STUDYID, ~SPDEVID, ~PARAMCD, ~PARAM,                ~AVISITN, ~AVISIT,    ~ABLFL, ~AVAL,
  "DEV01",  "D-001",  "BATTV",  "Battery Voltage (V)",        0, "BASELINE", "Y",      3.10,
  "DEV01",  "D-001",  "BATTV",  "Battery Voltage (V)",        4, "WEEK 4",   NA,       2.95,
  "DEV01",  "D-002",  "BATTV",  "Battery Voltage (V)",        0, "BASELINE", "Y",      3.08,
  "DEV01",  "D-002",  "BATTV",  "Battery Voltage (V)",        4, "WEEK 4",   NA,       2.71,
  "DEV01",  "D-001",  "OCCLN",  "Occlusion Count",            0, "BASELINE", "Y",         0,
  "DEV01",  "D-001",  "OCCLN",  "Occlusion Count",            4, "WEEK 4",   NA,          2,
  "DEV01",  "D-002",  "OCCLN",  "Occlusion Count",            0, "BASELINE", "Y",         0,
  "DEV01",  "D-002",  "OCCLN",  "Occlusion Count",            4, "WEEK 4",   NA,          1
) %>%
  admiral_df(dataset_keys = exprs(STUDYID, SPDEVID, PARAMCD, AVISITN), dataset_name = "ADDEVFND")
get_dataset_class(mdbds) # inferred
summary(mdbds)

## A derived record which leaves SPDEVID empty belongs to no device: it is not
## counted. The class must then be declared, as it is only inferred as device
## data when every record has a device
mdbds_avg <- bind_rows(
  mdbds,
  tibble(STUDYID = "DEV01", PARAMCD = "BATTV", AVISIT = "WEEK 4", AVISITN = 4, AVAL = 2.83)
)
get_dataset_class(mdbds_avg) # no longer inferred as MDBDS
mdbds_avg <- admiral_df(mdbds_avg, dataset_class = "MDBDS")
summary(mdbds_avg)

## MDOCCDS: device events; terms, severity, and the occurrence flags per device
mdoccds <- tribble(
  ~STUDYID, ~SPDEVID, ~DETERM,                     ~DEDECOD,           ~DESEV,     ~AOCCFL, ~AOCCPFL,
  "DEV01",  "D-001",  "PUMP OCCLUSION ALARM",      "OCCLUSION",        "MILD",     "Y",     "Y",
  "DEV01",  "D-001",  "PUMP OCCLUSION ALARM",      "OCCLUSION",        "MODERATE", NA,      NA,
  "DEV01",  "D-002",  "CRACKED RESERVOIR HOUSING", "MATERIAL FAILURE", "SEVERE",   "Y",     "Y",
  "DEV01",  "D-002",  "PUMP OCCLUSION ALARM",      "OCCLUSION",        "MILD",     NA,      "Y"
) %>%
  mutate(ASEQ = row_number(), .by = SPDEVID) %>%
  admiral_df(dataset_keys = exprs(STUDYID, SPDEVID, ASEQ), dataset_name = "ADDE")
get_dataset_class(mdoccds) # inferred
summary(mdoccds)

## MDTTE: time to device failure; events and censored per parameter. CNSR and
## STARTDT are what make it inferred as time-to-event data
mdtte <- tribble(
  ~STUDYID, ~SPDEVID, ~PARAMCD, ~PARAM,                         ~STARTDT,             ~AVAL, ~CNSR,
  "DEV01",  "D-001",  "TTDF",   "Time to Device Failure (days)", as.Date("2026-01-05"), 180,     1,
  "DEV01",  "D-002",  "TTDF",   "Time to Device Failure (days)", as.Date("2026-01-12"),  42,     0,
  "DEV01",  "D-003",  "TTDF",   "Time to Device Failure (days)", as.Date("2026-01-19"), 180,     1,
  "DEV01",  "D-004",  "TTDF",   "Time to Device Failure (days)", as.Date("2026-02-02"),  97,     0
) %>%
  admiral_df(dataset_keys = exprs(STUDYID, SPDEVID, PARAMCD), dataset_name = "ADDEVTTE")
get_dataset_class(mdtte) # inferred
summary(mdtte)

# 4. What the checks catch ----------------------------------------------------

## A second baseline per subject and parameter, which corrupts BASE and CHG
## but which the structure check cannot see (the visits differ)
advs_two_bl <- advs %>%
  mutate(ABLFL = if_else(PARAMCD == "SYSBP" & AVISIT == "Week 2", "Y", ABLFL)) %>%
  admiral_df(dataset_name = "ADVS")
summary(advs_two_bl)

## A second first occurrence per device and term
mdoccds_dup <- mutate(mdoccds, AOCCFL = "Y")
summary(mdoccds_dup)

## A merge which fans out duplicates subject records; the declared class keeps
## the ADSL details, and the structure check reports the extra records
adsl_dup <- bind_rows(adsl, head(adsl, 5)) %>%
  admiral_df(dataset_keys = exprs(STUDYID, USUBJID), dataset_name = "ADSL", dataset_class = "ADSL")
summary(adsl_dup)

## Keys renamed since they were set are reported, not checked on part of the key
advs_df %>%
  rename(ANLDT = ADT, ANLTPT = ATPT) %>%
  summary()

# 5. Datasets still being derived ---------------------------------------------

## Every item is optional: only what the variables allow is reported
adsl_wip <- admiral_df(select(adsl, STUDYID, USUBJID, AGE, SEX), dataset_class = "ADSL")
summary(adsl_wip)

## Typical use: check a derivation added what was expected. Here MAP is added
## for both subjects at the visits with both SYSBP and DIABP: 3 records
advs_small <- tribble(
  ~STUDYID,  ~USUBJID,      ~PARAMCD, ~AVISITN, ~AVISIT,    ~AVAL,
  "PILOT01", "01-701-1015", "DIABP",         0, "BASELINE",    51,
  "PILOT01", "01-701-1015", "SYSBP",         0, "BASELINE",   121,
  "PILOT01", "01-701-1015", "DIABP",         2, "WEEK 2",      55,
  "PILOT01", "01-701-1015", "SYSBP",         2, "WEEK 2",     125,
  "PILOT01", "01-701-1028", "DIABP",         0, "BASELINE",    79,
  "PILOT01", "01-701-1028", "SYSBP",         0, "BASELINE",   130,
  "PILOT01", "01-701-1028", "SYSBP",         2, "WEEK 2",     128
)
advs_map <- advs_small %>%
  derive_param_computed(
    by_vars = exprs(STUDYID, USUBJID, AVISITN, AVISIT),
    parameters = c("SYSBP", "DIABP"),
    set_values_to = exprs(
      AVAL = (AVAL.SYSBP + 2 * AVAL.DIABP) / 3,
      PARAMCD = "MAP"
    )
  ) %>%
  admiral_df(dataset_keys = exprs(STUDYID, USUBJID, PARAMCD, AVISITN), dataset_name = "ADVS")
summary(advs_map)

# 6. The summary object -------------------------------------------------------

## A list, so each piece can be used programmatically
s <- summary(adae_df, adsl = adsl)
names(s)
s$structure
s$details$terms
s$vs_adsl$incidence
