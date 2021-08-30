# Name: ADEX
#
# Label: Exposure Analysis Dataset
#
# Input: adsl, ex
#

library(dplyr)
library(lubridate)
library(stringr)
library(admiral)


# Use e.g. haven::read_sas to read in .sas7bdat, or other suitable functions
#  as needed and assign to the variables below.
# The CDISC pilot datasets are used for demonstration purpose.
data("adsl")
data("ex")

# The CDISC pilot data does not contain EXADJ,nor a SUPPVS dataset
# add a fake EXADJ to test for Dose adjustment flag
# add SUPPEX.EXPLDOS to test for dose intensity
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


# ---- Lookup tables ----

# Assign PARAMCD, PARAM, and PARAMN
param_lookup <- tibble::tribble(
  ~PARAMCD, ~PARAM, ~PARAMN,
  "DURD", "Study drug duration during constant dosing interval (days)", 1,
  "DOSE", "Dose administered during constant dosing interval (mg)", 2,
  "PLDOSE", "Planned dose during constant dosing interval (mg)", 3,
  "ADJ", "Dose adjusted during constant dosing interval", 4,
  "ADJAE", "Dose adjusted  due to AE during constant dosing interval", 5,
  "TDURD", "Overall duration (days)", 8,
  "TDOSE", "Total dose administered (mg)", 9,
  "ADOSE", "Average dose administered (mg)", 10,
  "TPDOSE", "Total planned dose (mg)", 11,
  "TADJ", "Dose adjusted during study", 13,
  "TADJAE", "Dose adjusted during study due to AE", 14,
  "PDURD", "Overall duration in W2-W24 (days)", 19,
  "PDOSE", "Total dose administered in W2-W2 (mg)4", 20,
  "PPDOSE", "Total planned dose in W2-W24 (mg)", 21,
  "PADOSE", "Average dose administered in W2-W24 (mg)", 22,
  "PADJ", "Dose adjusted during W2-W24", 24,
  "PADJAE", "Dose adjusted  in W2-W24 due to AE", 25,
  "TDOSINT", "Overall dose intensity (%)", 90,
  "PDOSINT", "W2-24 dose intensity (%)", 91
)


# ---- User defined functions ----

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

# ---- Derivations ----
adsl0 <- adsl %>%
  select(STUDYID, USUBJID, starts_with("TRT"))

# Part 1
# Join ADSL with ex and derive required dates, variables
adex0 <- ex %>%
  left_join(adsl0, by = c("STUDYID", "USUBJID")) %>%

  # Calculate ASTDTM, AENDTM using derive_vars_dtm
  call_derivation(
    derivation = derive_vars_dtm,
    variable_params = list(
      params(dtc = EXSTDTC, date_imputation = "first", new_vars_prefix = "AST"),
      params(dtc = EXENDTC, date_imputation = "last", new_vars_prefix = "AEN")
    ),
    flag_imputation = "none",
    min_dates = list(TRTSDT),
    max_dates = list(TRTEDT)
  ) %>%

  # Calculate ASTDY, AENDY
  derive_var_astdy(reference_date = TRTSDT, date = ASTDTM) %>%
  derive_var_aendy(reference_date = TRTSDT, date = AENDTM) %>%

  # add EXDUR, the duration of trt for each record
  derive_duration(
    new_var = EXDUR,
    new_var_unit = EXDURU,
    start_date = ASTDTM,
    end_date = AENDTM
  ) %>%
  mutate(ASTDT = date(ASTDTM), AENDT = date(AENDTM))

# Part 2
# 1:1 mapping

adex <- bind_rows(
  adex0 %>% mutate(PARAMCD = "DURD", AVAL = EXDUR),
  adex0 %>% mutate(PARAMCD = "DOSE", AVAL = EXDOSE),
  adex0 %>% mutate(PARAMCD = "PLDOSE", AVAL = EXPLDOS),
  adex0 %>% mutate(PARAMCD = "ADJ", AVALC = if_else(!is.na(EXADJ), "Y", NA_character_)),
  adex0 %>% mutate(PARAMCD = "ADJAE", AVALC = if_else(EXADJ == "ADVERSE EVENT", "Y", NA_character_))
) %>%
  mutate(PARCAT1 = "INDIVIDUAL") %>%

  # Part 3
  # Derive summary parameters

  # Overall exposure
  call_derivation(
    derivation = derive_exposure_params,
    variable_params = list(
      params(
        set_values_to = vars(PARAMCD = "TDOSE", PARCAT1 = "OVERALL"),
        input_param = "DOSE",
        fns = AVAL ~ sum(., na.rm = TRUE)
      ),
      params(
        set_values_to = vars(PARAMCD = "TPDOSE", PARCAT1 = "OVERALL"),
        input_param = "PLDOSE",
        fns = AVAL ~ sum(., na.rm = TRUE)
      ),
      params(
        set_values_to = vars(PARAMCD = "TDURD", PARCAT1 = "OVERALL"),
        input_param = "DURD",
        fns = AVAL ~ sum(., na.rm = TRUE)
      ),
      params(
        set_values_to = vars(PARAMCD = "AVDOSE", PARCAT1 = "OVERALL"),
        input_param = "DOSE",
        fns = AVAL ~ mean(., na.rm = TRUE)
      ),
      params(
        set_values_to = vars(PARAMCD = "TADJ", PARCAT1 = "OVERALL"),
        input_param = "ADJ",
        fns = AVALC ~ if_else(sum(!is.na(.)) > 0, "Y", NA_character_)
      ),
      params(
        set_values_to = vars(PARAMCD = "TADJAE", PARCAT1 = "OVERALL"),
        input_param = "ADJAE",
        fns = AVALC ~ if_else(sum(!is.na(.)) > 0, "Y", NA_character_)
      )
    ),
    by_vars = vars(STUDYID, USUBJID, TRT01P, TRT01A, TRTSDTM, TRTEDTM, TRTSDT, TRTEDT)
  ) %>%

  # W2-W24 exposure
  call_derivation(
    derivation = derive_exposure_params,
    variable_params = list(
      params(
        set_values_to = vars(PARAMCD = "PDOSE", PARCAT1 = "WEEK 2-24"),
        input_param = "DOSE",
        fns = AVAL ~ sum(., na.rm = TRUE)
      ),
      params(
        set_values_to = vars(PARAMCD = "PPDOSE", PARCAT1 = "WEEK 2-24"),
        input_param = "PLDOSE",
        fns = AVAL ~ sum(., na.rm = TRUE)
      ),
      params(
        set_values_to = vars(PARAMCD = "PDURD", PARCAT1 = "WEEK 2-24"),
        input_param = "DURD",
        fns = AVAL ~ sum(., na.rm = TRUE)
      ),
      params(
        set_values_to = vars(PARAMCD = "PAVDOSE", PARCAT1 = "WEEK 2-24"),
        input_param = "DOSE",
        fns = AVAL ~ mean(., na.rm = TRUE)
      ),
      params(
        set_values_to = vars(PARAMCD = "PADJ", PARCAT1 = "WEEK 2-24"),
        input_param = "ADJ",
        fns = AVALC ~ if_else(sum(!is.na(.)) > 0, "Y", NA_character_)
      ),
      params(
        set_values_to = vars(PARAMCD = "PADJAE", PARCAT1 = "WEEK 2-24"),
        input_param = "ADJAE",
        fns = AVALC ~ if_else(sum(!is.na(.)) > 0, "Y", NA_character_)
      )
    ),
    filter = VISIT %in% c("WEEK 2", "WEEK 24"),
    by_vars = vars(STUDYID, USUBJID, TRT01P, TRT01A, TRTSDTM, TRTEDTM, TRTSDT, TRTEDT)
  ) %>%

  # Overall Dose intensity and W2-24 dose intensity
  call_derivation(
    derivation = derive_param_doseint,
    variable_params = list(
      params(new_param = "TDOSINT", tadm_code = "TDOSE", tpadm_code = "TPDOSE"),
      params(new_param = "PDOSINT", tadm_code = "PDOSE", tpadm_code = "PPDOSE")
    ),
    by_vars = vars(STUDYID, USUBJID, TRT01P, TRT01A, TRTSDTM, TRTEDTM, TRTSDT, TRTEDT, PARCAT1)
  ) %>%


  # Part 4
  # Derive/Assign the last required variables
  # Add PARAMN and PARAM, AVALU
  left_join(param_lookup, by = "PARAMCD") %>%

  # Derive AVALCATx
  mutate(AVALCAT1 = format_avalcat1(param = PARAMCD, aval = AVAL)) %>%

  # Calculate ASEQ
  derive_obs_number(
    new_var = ASEQ,
    by_vars = vars(STUDYID, USUBJID),
    order = vars(ASTDT, VISIT, VISITNUM, EXSEQ, PARAMCD),
    check_type = "warning"
  )


# Final Steps, Select final variables and Add labels
# This process will be based on your metadata, no example given for this reason
# ...

# ---- Save output ----

saveRDS(adex, file = "./ADEX.rds", compress = TRUE)
