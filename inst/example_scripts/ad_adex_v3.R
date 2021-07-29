# Name: ADEX
#
# Label: Vital Signs Analysis Dataset
#
# Input: adsl, ex
#

library(dplyr)
library(lubridate)
library(stringr)
library(admiral)

# Use e.g. haven::read_sas to read in .sas7bdat, or other suitable functions
#  as needed and assign to the variables below.
adsl <- NULL
ex <- NULL

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
# Assign PARAMCD, PARAM, and PARAMN
param_lookup <- tibble::tribble(
  ~PARAMCD, ~PARAM, ~PARAMN,
  "DURD", "Study drug duration during constant dosing interval (days)", 1,
  "DOSE", "Dose administered during constant dosing interval", 2,
  "PLDOSE", "Planned dose during constant dosing interval", 3,
  "ADJ", "Dose adjusted during constant dosing interval", 4,
  "ADJAE", "Dose adjusted  due to AE during constant dosing interval", 5,
  "TDURD", "Overall duration (days)", 8,
  "AVDOSE", "Average dose", 9,
  "TDOSE", "Total dose administered", 10,
  "TPDOSE", "Total planned dose", 11,
  "TNDOSMIS", "Total number of missed doses", 12,
  "TADJ", "Dose adjusted during study", 13,
  "TADJAE", "Dose adjusted during study due to AE", 14,
  "PDURD", "Overall duration in W2-W24 (days)", 19,
  "AVDOSE", "Average dose in W2-W24", 20,
  "PDOSE", "Total dose administered in W2-W24", 21,
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
    TRUE~NA_character_
  )
}

# ---- Derivations ----

#Part 1
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
  # add EXDUR, the duration of trt for each record
  derive_duration(
    new_var = EXDUR,
    new_var_unit = EXDURU,
    start_date = ASTDTM,
    end_date = AENDTM
  )


#Part 2
#1:1 mapping -TEMPORARY
# will be appended once the derived param are avaialble
adex1_1<- bind_rows(
  adex0 %>% mutate(PARAMCD="DURD", AVAL=EXDUR),
  adex0 %>% mutate(PARAMCD="DOSE", AVAL=EXDOSE),
  adex0 %>% mutate(PARAMCD="AVDOSE", AVAL=EXDOSE),
  adex0 %>% mutate(PARAMCD="PLDOSE", AVAL=EXPLDOS),
  adex0 %>% mutate(PARAMCD="ADJ",  AVALC=if_else(!is.na(EXADJ), "Y",NA_character_)),
  adex0 %>% mutate(PARAMCD="ADJAE",  AVALC=if_else(EXADJ=="ADVERSE EVENT", "Y",NA_character_))
) %>%
  mutate(PARCAT1="INDIVIDUAL") %>%

#Part - Summary parameter derivation
  # derive summary records for OVERALL exposure
  #keep the dates
  derive_summary_records(
    by_vars = vars(STUDYID, USUBJID, EXTRT, PARAMCD),
    filter_rows = PARCAT1=="INDIVIDUAL",
    fns = list(AVAL ~ sum(., na.rm = TRUE),
               AVAL~ mean(na.omit(.), na.rm = TRUE),
               AVALC ~ if_else(sum(!is.na(.)) > 0, "Y", NA_character_),
               ASTDTM ~ min(., na.rm=TRUE),
               AENDTM ~ max(coalesce(AENDTM, ASTDTM), na.rm=TRUE)),
    set_values_to = vars(PARCAT1 = rep("OVERALL",5),
                         #PARAMCD=paste0("T",PARAMCD)
                         DTYPE= c("SUM", "AVERAGE", "", "", "")
                         ),
    drop_values_from = vars(EXPLDOS, EXDOSU, EXDOSFRM, EXDOSFRQ, EXROUTE, EXDURU)

  ) %>%
  #Remove non required record/values
  filter(!(PARAMCD %in% c("ADJ", "ADJAE") & DTYPE =="AVERAGE"),
         )%>%
  mutate(AVAL=if_else(PARAMCD %in% c("ADJ", "ADJAE"), NA_real_,AVAL ) )%>%
  group_by(STUDYID, USUBJID, EXSEQ, EXTRT,PARAMCD) %>%
  #Coalesce all rows by group
  summarise_all(~first(na.omit(.))) %>%
  #Re label PARAMCD for overall parameters
  mutate(PARCAT1=coalesce(PARCAT1, PARCAT1.y),
         PARAMCD= case_when(
           PARCAT1=="OVERALL" & DTYPE == "AVERAGE" ~ "AVDOSE",
           PARCAT1=="OVERALL" & PARAMCD !="PLDOSE" ~ paste0("T",PARAMCD),
           PARCAT1=="OVERALL" & PARAMCD =="PLDOSE" ~ "TPDOSE",
           TRUE ~ PARAMCD
         )
         ) %>%
  #remocve temp var
  select(-ends_with(".x"), -ends_with(".y")) %>%
  #ungroup() %>%
  # derive summary records for exposure in W2-W24
  derive_summary_records(
    by_vars = vars(STUDYID, USUBJID, EXTRT, PARAMCD),
    filter_rows = VISIT %in% c("WEEK 2", "WEEK 24") & PARCAT1=="INDIVIDUAL",
    fns = list(AVAL ~ sum(., na.rm = TRUE),
               AVALC ~ if_else(sum(!is.na(.)) > 0, "Y", NA_character_),
               ASTDTM ~ min(., na.rm=TRUE),
               AENDTM ~ max(coalesce(AENDTM, ASTDTM), na.rm=TRUE)
               ),
    set_values_to = vars(PARCAT1 = rep("W2-24",4)),
    drop_values_from = vars(EXPLDOS, EXDOSU, EXDOSFRM, EXDOSFRQ, EXROUTE, EXDURU)
  ) %>%

   #Re label PARAMCD for W2-W24 parameters
  mutate(PARCAT1=coalesce(PARCAT1, PARCAT1.y),
         PARAMCD= case_when(
           PARCAT1=="W2-24" & PARAMCD !="PLDOSE" ~ paste0("T",PARAMCD),
           PARCAT1=="W2-24" & PARAMCD =="PLDOSE" ~ "PDOSE",
           TRUE ~ PARAMCD
         )
  ) %>%
  #remocve temp var
  select(-ends_with(".x"), -ends_with(".y")) %>%
  group_by(STUDYID, USUBJID, EXTRT, PARAMCD) %>%
  #Coalesce all row by by group
  summarise_all(~first(na.omit(.)))%>%


    mutate(ASTDT=date(ASTDTM), AENDT=date(AENDTM)          )
  # Calculate ASTDY, AENDY
  derive_var_astdy(reference_date = TRTSDT, date = ASTDTM) %>%
  derive_var_aendy(reference_date = TRTSDT, date = AENDTM) %>%


      # Add PARAMCD and PARAM
  left_join(param_lookup, by = "PARAMCD") %>%
  # Derive AVALCATx
  mutate(AVALCAT1 = format_avalcat1(param = PARAMCD, aval = AVAL)) %>%
  # Calculate ASEQ
derive_obs_number(
  new_var = ASEQ,
  by_vars = vars(STUDYID, USUBJID),
  order = vars(PARAMCD, ASTDT, VISIT, VISITNUM, PARAMCD),
  check_type = "warning"
)



  # Final Steps, Select final variables and Add labels
  # This process will be based on your metadata, no example given for this reason
  # ...

  # ---- Save output ----

  save(adex, file = "/PATH/TO/SAVE/ADex", compress = TRUE)
