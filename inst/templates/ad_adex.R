# Name: ADEX
#
# Label: Exposure Analysis Dataset
#
# Input: adsl, ex
#

library(admiral)
library(pharmaversesdtm) # Contains example datasets from the CDISC pilot project
library(dplyr)
library(lubridate)
library(stringr)

# Load source datasets ----
# Use e.g. haven::read_sas to read in .sas7bdat, or other suitable functions
#  as needed and assign to the variables below.
# The CDISC pilot datasets are used for demonstration purpose.

ex <- pharmaversesdtm::ex
adsl <- admiral::admiral_adsl

# When SAS datasets are imported into R using haven::read_sas(), missing
# character values from SAS appear as "" characters in R, instead of appearing
# as NA values. Further details can be obtained via the following link:
# https://pharmaverse.github.io/admiral/cran-release/articles/admiral.html#handling-of-missing-values # nolint

ex <- convert_blanks_to_na(ex)

# The CDISC pilot data does not contain EXADJ,nor a SUPPEX dataset
# add a fake EXADJ to demonstrate the derivation for Dose adjustment flag
# add SUPPEX.EXPLDOS to demonstrate the derivation for dose intensity
# The CDISC pilot EX datasets, contains exposure data for daily dosing. Care should be taken when
# the dosing frequency is different.


ex <- ex %>%
  mutate(
    EXADJ = case_when(
      USUBJID == "01-701-1034" & VISIT %in% c("WEEK 2", "WEEK 24") ~ "ADVERSE EVENT",
      USUBJID == "01-701-1148" & VISIT %in% c("WEEK 24") ~ "MEDICATION ERROR",
      TRUE ~ NA_character_
    ),
    EXDOSE = case_when(
      USUBJID == "01-701-1034" & VISIT %in% c("WEEK 2", "WEEK 24") ~ 0,
      USUBJID == "01-701-1148" & VISIT %in% c("WEEK 24") ~ 0,
      TRUE ~ EXDOSE
    )
  ) %>%
  # add SUPPEX.EXPLDOS to test for dose intensity
  mutate(EXPLDOS = if_else(EXTRT == "PLACEBO", 0, 54))


# Derivations ----

# Get list of ADSL vars required for derivations
adsl_vars <- exprs(TRTSDT, TRTSDTM, TRTEDTM)

# Part 1
# Join ADSL with ex and derive required dates, variables
adex0 <- ex %>%
  # Join ADSL with EX (only ADSL vars required for derivations)
  derive_vars_merged(
    dataset_add = adsl,
    new_vars = adsl_vars,
    by_vars = exprs(STUDYID, USUBJID)
  ) %>%
  ## Calculate ASTDTM, AENDTM using `derive_vars_dtm()` ----
  derive_vars_dtm(
    dtc = EXSTDTC,
    highest_imputation = "M",
    new_vars_prefix = "AST"
  ) %>%
  derive_vars_dtm(
    dtc = EXENDTC,
    highest_imputation = "M",
    date_imputation = "last",
    new_vars_prefix = "AEN"
  ) %>%
  ## Calculate ASTDY, AENDY ----
  derive_vars_dy(
    reference_date = TRTSDTM,
    source_vars = exprs(ASTDTM, AENDTM)
  ) %>%
  ## Add EXDUR, the duration of trt for each record ----
  derive_vars_duration(
    new_var = EXDURD,
    start_date = ASTDTM,
    end_date = AENDTM
  ) %>%
  ## Derive analysis end/start date ----
  derive_vars_dtm_to_dt(exprs(ASTDTM, AENDTM)) %>%
  mutate(
    # Compute the cumulative dose
    DOSEO = EXDOSE * EXDURD,
    PDOSEO = EXPLDOS * EXDURD
  )

# Part 2: 1:1 mapping ----

adex <- bind_rows(
  adex0 %>% mutate(PARAMCD = "DURD", AVAL = EXDURD),
  adex0 %>% mutate(PARAMCD = "DOSE", AVAL = DOSEO),
  adex0 %>% mutate(PARAMCD = "PLDOSE", AVAL = PDOSEO),
  adex0 %>% mutate(PARAMCD = "ADJ", AVALC = if_else(!is.na(EXADJ), "Y", NA_character_)),
  adex0 %>% mutate(PARAMCD = "ADJAE", AVALC = if_else(EXADJ == "ADVERSE EVENT", "Y", NA_character_))
) %>%
  mutate(PARCAT1 = "INDIVIDUAL")

# Part 3: Derive summary parameters ----
# Note that, for the functions `derive_param_exposure()`,
# `derive_param_doseint()` and `derive_param_computed()`, only the variables
# specified in `by_vars` will be populated in the newly created records.

adex <- adex %>%
  # Overall exposure
  call_derivation(
    derivation = derive_param_exposure,
    variable_params = list(
      params(
        set_values_to = exprs(
          PARAMCD = "TDOSE",
          PARCAT1 = "OVERALL",
          AVAL = sum(AVAL, na.rm = TRUE)
        ),
        input_code = "DOSE"
      ),
      params(
        set_values_to = exprs(
          PARAMCD = "TPDOSE",
          PARCAT1 = "OVERALL",
          AVAL = sum(AVAL, na.rm = TRUE)
        ),
        input_code = "PLDOSE"
      ),
      params(
        set_values_to = exprs(
          PARAMCD = "TDURD",
          PARCAT1 = "OVERALL",
          AVAL = sum(AVAL, na.rm = TRUE)
        ),
        input_code = "DURD"
      ),
      params(
        set_values_to = exprs(
          PARAMCD = "TADJ",
          PARCAT1 = "OVERALL",
          AVALC = if_else(sum(!is.na(AVALC)) > 0, "Y", NA_character_)
        ),
        input_code = "ADJ"
      ),
      params(
        set_values_to = exprs(
          PARAMCD = "TADJAE",
          PARCAT1 = "OVERALL",
          AVALC = if_else(sum(!is.na(AVALC)) > 0, "Y", NA_character_)
        ),
        input_code = "ADJAE"
      )
    ),
    dataset_add = adex,
    by_vars = exprs(STUDYID, USUBJID, !!!adsl_vars)
  ) %>%
  # W2-W24 exposure
  call_derivation(
    derivation = derive_param_exposure,
    variable_params = list(
      params(
        set_values_to = exprs(
          PARAMCD = "PDOSE",
          PARCAT1 = "WEEK 2-24",
          AVAL = sum(AVAL, na.rm = TRUE)
        ),
        input_code = "DOSE"
      ),
      params(
        set_values_to = exprs(
          PARAMCD = "PPDOSE",
          PARCAT1 = "WEEK 2-24",
          AVAL = sum(AVAL, na.rm = TRUE)
        ),
        input_code = "PLDOSE"
      ),
      params(
        set_values_to = exprs(
          PARAMCD = "PDURD",
          PARCAT1 = "WEEK 2-24",
          AVAL = sum(AVAL, na.rm = TRUE)
        ),
        input_code = "DURD"
      ),
      params(
        set_values_to = exprs(
          PARAMCD = "PADJ",
          PARCAT1 = "WEEK 2-24",
          AVALC = if_else(sum(!is.na(AVALC)) > 0, "Y", NA_character_)
        ),
        input_code = "ADJ"
      ),
      params(
        set_values_to = exprs(
          PARAMCD = "PADJAE",
          PARCAT1 = "WEEK 2-24",
          AVALC = if_else(sum(!is.na(AVALC)) > 0, "Y", NA_character_)
        ),
        input_code = "ADJAE"
      )
    ),
    dataset_add = adex,
    filter_add = VISIT %in% c("WEEK 2", "WEEK 24"),
    by_vars = exprs(STUDYID, USUBJID, !!!adsl_vars)
  ) %>%
  # Overall Dose intensity and W2-24 dose intensity
  call_derivation(
    derivation = derive_param_doseint,
    variable_params = list(
      params(
        set_values_to = exprs(PARAMCD = "TDOSINT"),
        tadm_code = "TDOSE",
        tpadm_code = "TPDOSE"
      ),
      params(
        set_values_to = exprs(PARAMCD = "PDOSINT"),
        tadm_code = "PDOSE",
        tpadm_code = "PPDOSE"
      )
    ),
    by_vars = exprs(
      STUDYID, USUBJID, !!!adsl_vars, PARCAT1, ASTDTM, ASTDT, AENDTM, AENDT
    )
  ) %>%
  # Overall/W2-24 Average daily dose
  call_derivation(
    derivation = derive_param_computed,
    variable_params = list(
      params(
        parameters = c("TDOSE", "TDURD"),
        set_values_to = exprs(
          AVAL = (AVAL.TDOSE / AVAL.TDURD),
          PARAMCD = "AVDDSE"
        )
      ),
      params(
        parameters = c("PDOSE", "PDURD"),
        set_values_to = exprs(
          AVAL = (AVAL.PDOSE / AVAL.PDURD),
          PARAMCD = "PAVDDSE"
        )
      )
    ),
    by_vars = exprs(
      STUDYID, USUBJID, !!!adsl_vars, PARCAT1, ASTDTM, ASTDT, AENDTM, AENDT
    )
  )

# Part 4: Derive/Assign the last required variables ----

# Assign PARAMCD, PARAM, and PARAMN
# ---- Lookup tables ----
param_lookup <- tibble::tribble(
  ~PARAMCD,                                                     ~PARAM, ~PARAMN,
  "DURD", "Study drug duration during constant dosing interval (days)",       1,
  "DOSE",     "Dose administered during constant dosing interval (mg)",       2,
  "PLDOSE",        "Planned dose during constant dosing interval (mg)",       3,
  "ADJ",               "Dose adjusted during constant dosing interval",       4,
  "ADJAE",  "Dose adjusted  due to AE during constant dosing interval",       5,
  "TDURD",                                   "Overall duration (days)",       7,
  "TDOSE",                              "Total dose administered (mg)",       8,
  "AVDDSE",                  "Average daily dose administered (mg/mg)",      10,
  "TPDOSE",                                  "Total planned dose (mg)",      11,
  "TADJ",                                 "Dose adjusted during study",      13,
  "TADJAE",                     "Dose adjusted during study due to AE",      14,
  "PDURD",                         "Overall duration in W2-W24 (days)",      19,
  "PDOSE",                    "Total dose administered in W2-W2 (mg)4",      20,
  "PPDOSE",                        "Total planned dose in W2-W24 (mg)",      21,
  "PAVDDSE",          "Average daily dose administered in W2-W24 (mg)",      23,
  "PADJ",                                "Dose adjusted during W2-W24",      24,
  "PADJAE",                       "Dose adjusted  in W2-W24 due to AE",      25,
  "TDOSINT",                              "Overall dose intensity (%)",      90,
  "PDOSINT",                                "W2-24 dose intensity (%)",      91
)

# Assign AVALCATx
avalcax_lookup <- exprs(
  ~PARAMCD,            ~condition,             ~AVALCAT1,
  "TDURD",             AVAL >= 90,          ">= 90 days",
  "TDURD", AVAL >= 30 & AVAL < 90, ">= 30 and < 90 days",
  "TDURD",              AVAL < 30,           "< 30 days",
  "PDURD",             AVAL >= 90,          ">= 90 days",
  "PDURD", AVAL >= 30 & AVAL < 90, ">= 30 and < 90 days",
  "PDURD",              AVAL < 30,           "< 30 days",
  "TDOSE",             AVAL < 100,            "< 100 mg",
  "TDOSE",            AVAL >= 100,           ">= 100 mg",
  "PDOSE",             AVAL < 100,            "< 100 mg",
  "PDOSE",            AVAL >= 100,           ">= 100 mg"
)

adex <- adex %>%
  # Add PARAMN and PARAM, AVALU
  derive_vars_merged(
    dataset_add = param_lookup,
    by_vars = exprs(PARAMCD)
  ) %>%
  # Derive AVALCATx
  derive_vars_cat(
    definition = avalcax_lookup,
    by_vars = exprs(PARAMCD)
  ) %>%
  # Calculate ASEQ
  derive_var_obs_number(
    new_var = ASEQ,
    by_vars = exprs(STUDYID, USUBJID),
    order = exprs(PARCAT1, ASTDT, VISIT, VISITNUM, EXSEQ, PARAMN),
    check_type = "error"
  )

# Join all ADSL with EX
adex <- adex %>%
  derive_vars_merged(
    dataset_add = select(adsl, !!!negate_vars(adsl_vars)),
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
save(adex, file = file.path(dir, "adex.rda"), compress = "bzip2")
