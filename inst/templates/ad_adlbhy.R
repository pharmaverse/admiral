# Name: ADLBHY
#
# Label: Lab Analysis Dataset for Hy's Law
#
# Input: adlb
library(admiral)
library(admiral.test) # Contains example datasets from the CDISC pilot project
library(dplyr)
library(lubridate)

# Load source datasets ----

# Use e.g. haven::read_sas to read in .sas7bdat, or other suitable functions
# as needed and assign to the variables below.
# For illustration purposes read in admiral test data
# Using use_ad_template("adlb") and assigning the end object as admiral_adlb
# ADLBHY is a special dataset specifically used to check for potential drug induced liver injuries
# Please see "Hy's Law Implementation Guide" on the admiral website for additional information


data("admiral_adlb")
adlb <- admiral_adlb

adlb_annotated <- adlb %>%
  filter(PARAMCD %in% c("AST", "ALT", "BILI") & is.na(DTYPE)) %>%
  mutate(CRIT1 = case_when(
    PARAMCD == "AST" ~ "AST >=3xULN",
    PARAMCD == "ALT" ~ "ALT >=3xULN",
    PARAMCD == "BILI" ~ "BILI >=2xULN"
  )) %>%
  # Assign flags for contribution to potential Hy's Law event
  call_derivation(
    .,
    derivation = derive_var_merged_exist_flag,
    dataset_add = .,
    by_vars = exprs(USUBJID, LBSEQ, PARAMCD, ADT),
    variable_params = list(
      params(
        new_var = CRIT1FL,
        condition = (
          (AVAL / ANRHI >= 3 & PARAMCD %in% c("AST", "ALT")) |
            (AVAL / ANRHI >= 2 & PARAMCD == "BILI")
        )
      )
    )
  ) %>%
  select(STUDYID, USUBJID, TRT01A, PARAMCD, LBSEQ, ADT, AVISIT, ADY, AVAL, ANRHI, CRIT1, CRIT1FL)

# Subset Datasets
altast_records <- adlb_annotated %>%
  filter(PARAMCD %in% c("AST", "ALT"))

bili_records <- adlb_annotated %>%
  filter(PARAMCD %in% c("BILI"))

# Use a join and filter to accomplish time-window search
hylaw_records <- derive_vars_joined(
  dataset = altast_records,
  dataset_add = bili_records,
  by_vars = exprs(STUDYID, USUBJID, ADY),
  order = exprs(ADY),
  filter_join = ADT.join - ADT <= 14 & CRIT1FL == "Y" & CRIT1FL.join == "Y",
  new_vars = exprs(BILI_LBSEQ = LBSEQ, BILI_DT = ADT, BILI_CRITFL = CRIT1FL),
  mode = "first"
)

hylaw_records_pts_visits <- hylaw_records %>%
  select(STUDYID, USUBJID, TRT01A) %>% # add AVISIT, ADT for by visit
  distinct()

hylaw_records_fls <- hylaw_records %>%
  select(STUDYID, USUBJID, TRT01A, CRIT1FL, BILI_CRITFL) %>% # add AVISIT, ADT for by visit
  distinct()

# Create new parameters based on records that present potential case
hylaw_params <- derive_param_exist_flag(
  dataset_adsl = hylaw_records_pts_visits,
  dataset_add = hylaw_records_fls,
  condition = CRIT1FL == "Y" & BILI_CRITFL == "Y",
  false_value = "N",
  missing_value = "N",
  subject_keys = exprs(STUDYID, USUBJID, TRT01A), # add AVISIT, ADT for by visit
  set_values_to = exprs(
    PARAMCD = "HYSLAW",
    PARAM = "ALT/AST >= 3xULN and BILI >= 2xULN"
  )
)

# Row bind back to relevant adlb-like dataset
adlbhy <- adlb_annotated %>%
  bind_rows(hylaw_params)
