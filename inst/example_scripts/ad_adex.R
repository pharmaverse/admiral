# Name: ADEX
#
# Label: Vital Signs Analysis Dataset
#
# Input: adsl, ex
#

library(dplyr)
library(lubridate)
library(rlang)
library(styler)
library(roxygen2)
library(devtools)
library(diffdf)
library(assertthat)
library(testthat)
library(stringr)
library(styler)
library(purrr)
# library(haven)
library(admiral)
devtools::load_all()

# Use e.g. haven::read_sas to read in .sas7bdat, or other suitable functions
#  as needed and assign to the variables below.
# adsl <- NULL
# ex <- NULL

data("adsl")
data("ex")

# The CDISC pilot data does not contain EXADJ,
# add a fake one to test for Dose adjustment flag
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
  mutate(EXPLDOS = 54)


# ---- Lookup tables ----

# Assign PARAMCD, PARAM, and PARAMN
param_lookup <- tibble::tribble(
  ~PARAMCD, ~PARAM, ~PARAMN,
  "DURD", "Study drug duration during constant dosing interval (days)", 1,
  "DOSE", "Dose administered during constant dosing interval", 2,
  "PLDOSE", "Planned dose during constant dosing interval", 3,
  "ADJ", "Dose adjusted during constant dosing interval", 4,
  "ADJAE", "Dose adjusted  due to AE during constant dosing interval", 5,
  "TDURD", "Overall duration (days)", 8,
  "TDOSE", "Total dose administered", 9,
  "ADOSE", "Average dose administered", 10,
  "TPDOSE", "Total planned dose", 11,
  "TNDOSMIS", "Total number of missed doses", 12,
  "TADJ", "Dose adjusted during study", 13,
  "TADJAE", "Dose adjusted during study due to AE", 14,
  "PDURD", "Overall duration in W2-W24 (days)", 19,
  "PDOSE", "Total dose administered in W2-W24", 20,
  "PADOSE", "Average dose administered in W2-W24", 21,
  "PNDOSMIS", "Total number of missed doses in W2-W24", 22,
  "PADJ", "Dose adjusted during W2-W24", 23,
  "PADJAE", "Dose adjusted  in W2-W24 due to AE", 24,
  "TDOSINT", "Overall dose intensity", 90
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

# Part 1
# Join ADSL with ex and derive required dates, variables
adex0 <- adsl %>%
  select(STUDYID, USUBJID, TRTSDT) %>%
  right_join(ex, by = c("STUDYID", "USUBJID")) %>%
  # Calculate ASTDTM, AENDTM
  derive_vars_dtm(
    new_vars_prefix = "AST",
    dtc = EXSTDTC,
    flag_imputation = FALSE
  ) %>%
  derive_vars_dtm(
    new_vars_prefix = "AEN",
    dtc = EXENDTC,
    flag_imputation = FALSE
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
# will be appended once the derived param are available
adex1_1 <- bind_rows(
  adex0 %>% mutate(PARAMCD = "DURD", AVAL = EXDUR),
  adex0 %>% mutate(PARAMCD = "DOSE", AVAL = EXDOSE),
  adex0 %>% mutate(PARAMCD = "PLDOSE", AVAL = EXPLDOS),
  adex0 %>% mutate(PARAMCD = "ADJ", AVALC = if_else(!is.na(EXADJ), "Y", NA_character_)),
  adex0 %>% mutate(PARAMCD = "ADJAE", AVALC = if_else(EXADJ == "ADVERSE EVENT", "Y", NA_character_))
) %>%
  mutate(PARCAT1 = "INDIVIDUAL")

# Part 3 - summary parameters

# 3.1 - Define exposure parameters specification
# Overall
# Total dose and total duration
tdose <- exposure_var_spec(
  new_parameter = "TDOSE",
  input_parameter = "DOSE",
  fns = list(AVAL ~ sum(., na.rm = TRUE))
)
# Total duration
tdurd <- exposure_var_spec(
  new_parameter = "TDURD",
  input_parameter = "DURD",
  fns = list(AVAL ~ sum(., na.rm = TRUE))
)
# Average overall dose
adose <- exposure_var_spec(
  new_parameter = "ADOSE",
  input_parameter = "DOSE",
  fns = list(AVAL ~ mean(., na.rm = TRUE))
)
# Any Dose adjustment?
tadj <- exposure_var_spec(
  new_parameter = "TADJ",
  input_parameter = "ADJ",
  fns = list(AVALC ~ if_else(sum(!is.na(.)) > 0, "Y", NA_character_))
)
# Any Dose adjustment due to AE?
tadjae <- exposure_var_spec(
  new_parameter = "TADJAE",
  input_parameter = "ADJAE",
  fns = list(AVALC ~ if_else(sum(!is.na(. == "ADVERSE EVENT")) > 0, "Y", NA_character_))
)

# W2-24 expousre
pdose <- exposure_var_spec(
  new_parameter = "PDOSE",
  input_parameter = "DOSE",
  fns = list(AVAL ~ sum(., na.rm = TRUE))
)
pdurd <- exposure_var_spec(
  new_parameter = "PDURD",
  input_parameter = "DURD",
  fns = list(AVAL ~ sum(., na.rm = TRUE))
)
padose <- exposure_var_spec(
  new_parameter = "PADOSE",
  input_parameter = "DOSE",
  fns = list(AVAL ~ mean(., na.rm = TRUE))
)
padj <- exposure_var_spec(
  new_parameter = "PADJ",
  input_parameter = "ADJ",
  fns = list(AVALC ~ if_else(sum(!is.na(.)) > 0, "Y", NA_character_))
)
padjae <- exposure_var_spec(
  new_parameter = "PADJAE",
  input_parameter = "ADJAE",
  fns = list(AVALC ~ if_else(sum(!is.na(. == "ADVERSE EVENT")) > 0, "Y", NA_character_))
)
# 3.2 Derive summary parameters
adex <- adex1_1 %>%
  # add overall exposure
  derive_exposure_parameters(
    filter_rows = NULL,
    by_vars = vars(USUBJID, EXTRT),
    set_values_to = vars(PARCAT1 = "OVERALL"),
    drop_values_from = vars(EXPLDOS, EXDOSU, EXDOSFRM, EXDOSFRQ, EXROUTE, EXDURU),
    tdose , tdurd, adose, tadj, tadjae
  ) %>%
  # add W2-24 exposure
  derive_exposure_parameters(
    filter_rows = VISIT %in% c("WEEK 2", "WEEK 24"),
    by_vars = vars(USUBJID, EXTRT),
    set_values_to = vars(PARCAT1 = "W2-24"),
    drop_values_from = vars(EXPLDOS, EXDOSU, EXDOSFRM, EXDOSFRQ, EXROUTE, EXDURU),
    sum_w2_24, avg_w2_24, adj_w2_24, adjae_w2_24
  )


# Part 4
# Dose intensity goes there
# ...


# Part 5 Set 1:1 mapping and exposure parameters together
# Derive/Assign the last required variables

adex <- adex %>%
  # Add PARAMN and PARAM
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

save(adex, file = "/PATH/TO/SAVE/ADex", compress = TRUE)
