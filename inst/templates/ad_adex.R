# Name: ADEX
#
# Label: Exposure Analysis Dataset
#
# Input: adsl, ex
#

library(admiral)
library(admiral.test) # Contains example datasets from the CDISC pilot project
library(dplyr)
library(lubridate)
library(stringr)

# Use e.g. haven::read_sas to read in .sas7bdat, or other suitable functions
#  as needed and assign to the variables below.
# The CDISC pilot datasets are used for demonstration purpose.
data("adsl")
data("ex")

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


# ---- Derivations ----

# Get list of ADSL vars required for derivations
adsl_vars <- vars(TRTSDT, TRTSDTM, TRTEDTM)

# Part 1
# Join ADSL with ex and derive required dates, variables
adex0 <- ex %>%
  # Join ADSL with EX (only ADSL vars required for derivations)
  left_join(
    adsl %>% select(STUDYID, USUBJID, !!!adsl_vars),
    by = c("STUDYID", "USUBJID")
  ) %>%

  # Calculate ASTDTM, AENDTM using derive_vars_dtm
  derive_vars_dtm(dtc = EXSTDTC, date_imputation = "first", new_vars_prefix = "AST") %>%
  derive_vars_dtm(dtc = EXENDTC, date_imputation = "last", new_vars_prefix = "AEN") %>%

  # Calculate ASTDY, AENDY
  derive_var_astdy(date = ASTDTM) %>%
  derive_var_aendy(date = AENDTM) %>%

  # add EXDUR, the duration of trt for each record
  derive_vars_duration(
    new_var = EXDURD,
    start_date = ASTDTM,
    end_date = AENDTM
  ) %>%

  # Derive analysis end/start date
  derive_vars_dtm_to_dt(vars(ASTDTM, AENDTM)) %>%

  mutate(
    # Compute the cumulative dose
    DOSEO = EXDOSE * EXDURD,
    PDOSEO = EXPLDOS * EXDURD
  )

# Part 2
# 1:1 mapping

adex <- bind_rows(
  adex0 %>% mutate(PARAMCD = "DURD", AVAL = EXDURD),
  adex0 %>% mutate(PARAMCD = "DOSE", AVAL = DOSEO),
  adex0 %>% mutate(PARAMCD = "PLDOSE", AVAL = PDOSEO),
  adex0 %>% mutate(PARAMCD = "ADJ", AVALC = if_else(!is.na(EXADJ), "Y", NA_character_)),
  adex0 %>% mutate(PARAMCD = "ADJAE", AVALC = if_else(EXADJ == "ADVERSE EVENT", "Y", NA_character_))
) %>%
  mutate(PARCAT1 = "INDIVIDUAL")

  # Part 3
  # Derive summary parameters
adex <- adex %>%
  # Overall exposure
  call_derivation(
    derivation = derive_param_exposure,
    variable_params = list(
      params(
        set_values_to = vars(PARAMCD = "TDOSE", PARCAT1 = "OVERALL"),
        input_code = "DOSE",
        analysis_var = AVAL,
        summary_fun = function(x) sum(x, na.rm = TRUE)
      ),
      params(
        set_values_to = vars(PARAMCD = "TPDOSE", PARCAT1 = "OVERALL"),
        input_code = "PLDOSE",
        analysis_var = AVAL,
        summary_fun = function(x) sum(x, na.rm = TRUE)
      ),
      params(
        set_values_to = vars(PARAMCD = "TDURD", PARCAT1 = "OVERALL"),
        input_code = "DURD",
        analysis_var = AVAL,
        summary_fun = function(x) sum(x, na.rm = TRUE)
      ),
      params(
        set_values_to = vars(PARAMCD = "TADJ", PARCAT1 = "OVERALL"),
        input_code = "ADJ",
        analysis_var = AVALC,
        summary_fun = function(x) if_else(sum(!is.na(x)) > 0, "Y", NA_character_)
      ),
      params(
        set_values_to = vars(PARAMCD = "TADJAE", PARCAT1 = "OVERALL"),
        input_code = "ADJAE",
        analysis_var = AVALC,
        summary_fun = function(x) if_else(sum(!is.na(x)) > 0, "Y", NA_character_)
      )
    ),
    by_vars = vars(STUDYID, USUBJID, !!!adsl_vars)
  ) %>%

  # W2-W24 exposure
  call_derivation(
    derivation = derive_param_exposure,
    variable_params = list(
      params(
        set_values_to = vars(PARAMCD = "PDOSE", PARCAT1 = "WEEK 2-24"),
        input_code = "DOSE",
        analysis_var = AVAL,
        summary_fun = function(x) sum(x, na.rm = TRUE)
      ),
      params(
        set_values_to = vars(PARAMCD = "PPDOSE", PARCAT1 = "WEEK 2-24"),
        input_code = "PLDOSE",
        analysis_var = AVAL,
        summary_fun = function(x) sum(x, na.rm = TRUE)
      ),
      params(
        set_values_to = vars(PARAMCD = "PDURD", PARCAT1 = "WEEK 2-24"),
        input_code = "DURD",
        analysis_var = AVAL,
        summary_fun = function(x) sum(x, na.rm = TRUE)
      ),
      params(
        set_values_to = vars(PARAMCD = "PADJ", PARCAT1 = "WEEK 2-24"),
        input_code = "ADJ",
        analysis_var = AVALC,
        summary_fun = function(x) if_else(sum(!is.na(x)) > 0, "Y", NA_character_)
      ),
      params(
        set_values_to = vars(PARAMCD = "PADJAE", PARCAT1 = "WEEK 2-24"),
        input_code = "ADJAE",
        analysis_var = AVALC,
        summary_fun = function(x) if_else(sum(!is.na(x)) > 0, "Y", NA_character_)
      )
    ),
    filter = VISIT %in% c("WEEK 2", "WEEK 24"),
    by_vars = vars(STUDYID, USUBJID, !!!adsl_vars)
  ) %>%

  # Overall Dose intensity and W2-24 dose intensity
  call_derivation(
    derivation = derive_param_doseint,
    variable_params = list(
      params(set_values_to = vars(PARAMCD = "TDOSINT"), tadm_code = "TDOSE", tpadm_code = "TPDOSE"),
      params(set_values_to = vars(PARAMCD = "PDOSINT"), tadm_code = "PDOSE", tpadm_code = "PPDOSE")
    ),
    by_vars = vars(
      STUDYID, USUBJID, !!!adsl_vars, PARCAT1, ASTDTM, ASTDT, AENDTM, AENDT
    )
  ) %>%
  # Overall/W2-24 Average daily dose
  call_derivation(
    derivation = derive_derived_param,
    variable_params = list(
      params(
        parameters = c("TDOSE", "TDURD"),
        analysis_value = (AVAL.TDOSE / AVAL.TDURD),
        set_values_to = vars(PARAMCD = "AVDDSE")
      ),
      params(
        parameters = c("PDOSE", "PDURD"),
        analysis_value = (AVAL.PDOSE / AVAL.PDURD),
        set_values_to = vars(PARAMCD = "PAVDDSE")
      )
    ),
    by_vars = vars(
      STUDYID, USUBJID, !!!adsl_vars, PARCAT1, ASTDTM, ASTDT, AENDTM, AENDT
    )
  )

# Part 4
# Derive/Assign the last required variables

# Assign PARAMCD, PARAM, and PARAMN
# ---- Lookup tables ----
param_lookup <- tibble::tribble(
  ~PARAMCD, ~PARAM, ~PARAMN,
  "DURD", "Study drug duration during constant dosing interval (days)", 1,
  "DOSE", "Dose administered during constant dosing interval (mg)", 2,
  "PLDOSE", "Planned dose during constant dosing interval (mg)", 3,
  "ADJ", "Dose adjusted during constant dosing interval", 4,
  "ADJAE", "Dose adjusted  due to AE during constant dosing interval", 5,
  "TDURD", "Overall duration (days)", 7,
  "TDOSE", "Total dose administered (mg)", 8,
  "AVDDSE", "Average daily dose administered (mg/mg)", 10,
  "TPDOSE", "Total planned dose (mg)", 11,
  "TADJ", "Dose adjusted during study", 13,
  "TADJAE", "Dose adjusted during study due to AE", 14,
  "PDURD", "Overall duration in W2-W24 (days)", 19,
  "PDOSE", "Total dose administered in W2-W2 (mg)4", 20,
  "PPDOSE", "Total planned dose in W2-W24 (mg)", 21,
  "PAVDDSE", "Average daily dose administered in W2-W24 (mg)", 23,
  "PADJ", "Dose adjusted during W2-W24", 24,
  "PADJAE", "Dose adjusted  in W2-W24 due to AE", 25,
  "TDOSINT", "Overall dose intensity (%)", 90,
  "PDOSINT", "W2-24 dose intensity (%)", 91
)


# ---- User defined functions ----
# Derive AVALCAT1
# Here are some examples of how you can create your own functions that
#  operates on vectors, which can be used in `mutate`.
format_avalcat1 <- function(param, aval) {
  case_when(
    param %in% c("TDURD", "PDURD") & aval < 30 & !is.na(aval) ~ "< 30 days",
    param %in% c("TDURD", "PDURD") & aval >= 30 & aval < 90 ~ ">= 30 and < 90 days",
    param %in% c("TDURD", "PDURD") & aval >= 90 ~ ">=90 days",
    param %in% c("TDOSE", "PDOSE") & aval < 100 & !is.na(aval) ~ "< 100 mg",
    param %in% c("TDOSE", "PDOSE") & aval >= 100 ~ ">= 100 mg",
    TRUE ~ NA_character_
  )
}

adex <- adex %>%
  # Add PARAMN and PARAM, AVALU
  left_join(param_lookup, by = "PARAMCD") %>%

  # Derive AVALCATx
  mutate(AVALCAT1 = format_avalcat1(param = PARAMCD, aval = AVAL)) %>%

  # Calculate ASEQ
  derive_var_obs_number(
    new_var = ASEQ,
    by_vars = vars(STUDYID, USUBJID),
    order = vars(PARCAT1, ASTDT, VISIT, VISITNUM, EXSEQ, PARAMN),
    check_type = "error"
  )

# Join all ADSL with EX
adex <- adex %>%

  left_join(select(adsl, !!!admiral:::negate_vars(adsl_vars)),
            by = c("STUDYID", "USUBJID")
  )

# Final Steps, Select final variables and Add labels
# This process will be based on your metadata, no example given for this reason
# ...

# ---- Save output ----

save(adex, file = "data/adex.rda", compress = "bzip2")
