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
param_lookup <- tibble::tribble(
  ~PARAMCD, ~PARAM, ~PARAMN,
  "DURD", "Study drug duration during constant dosing interval (days)", 1,
  "DOSE", "Dose administered during constant dosing interval", 2,
  "PLDOSE", "Planned dose during constant dosing interval", 3,
  "ADJ", "Dose adjusted during constant dosing interval", 4,
  "ADJAE", "Dose adjusted  due to AE during constant dosing interval", 5,
  "TDURD", "Overall duration (days)", 9,
  "TDOSE", "Total dose administered", 10,
  "TPDOSE", "Total planned dose", 11,
  "TNDOSMIS", "Total number of missed doses", 12,
  "TADJ", "Dose adjusted during study", 13,
  "TADJAE", "Dose adjusted during study due to AE", 14,
  "PDURD", "Overall duration in W2-W24 (days)", 20,
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
  # Calculate ASTDY, AENDY
  derive_var_astdy(reference_date = TRTSDT, date = ASTDTM) %>%
  derive_var_aendy(reference_date = TRTSDT, date = AENDTM) %>%
  # add EXDUR, the duration of trt for each record
  derive_duration(
    new_var = EXDUR,
    new_var_unit = EXDURU,
    start_date = ASTDTM,
    end_date = AENDTM
  )%>%
  mutate(ASTDT=date(ASTDTM), AENDT=date(AENDTM))

#Part 2
#1:1 mapping
adex1<- rbind(
  adex0 %>% mutate(PARAMCD="DURD", AVAL=EXDUR, AVALC=NA_character_),
  adex0 %>% mutate(PARAMCD="DOSE", AVAL=EXDOSE, AVALC=NA_character_),
  adex0 %>% mutate(PARAMCD="PLDOSE", AVAL=EXPLDOS, AVALC=NA_character_),
  adex0 %>% mutate(PARAMCD="ADJ", AVAL=NA, AVALC=if_else(!is.na(EXADJ), "Y",NA_character_)),
  adex0 %>% mutate(PARAMCD="ADJAE", AVAL=NA, AVALC=if_else(EXADJ=="ADVERSE EVENT", "Y",NA_character_))
) %>%

#Part 3 - summary parameters
  # derive summary records for OVERALL exposure:
  #duration
  derive_summary_records(
    by_vars = vars(STUDYID, USUBJID, EXTRT),
    filter_rows = PARAMCD=="DURD",
    #fns = list(AVAL ~ sum(., na.rm = TRUE) ),
    fns = list(AVAL ~ date(AENDTM)-date(ASTDTM)+1),
    set_values_to = vars(
      PARCAT1 = rep("OVERALL", 1),
      PARAMCD = c("TDURD")
      #does not work
      #,ASTDTM=min(ASTDTM, na.rm=TRUE)
    ),
    drop_values_from = vars(EXDOSU, EXDOSFRM, EXDOSFRQ, EXROUTE, EXDURU)
  ) %>%
  #add start date/end dates
  group_by(STUDYID, USUBJID, EXTRT) %>%
  mutate(
    ASTDTM= case_when(PARCAT1=="OVERALL" ~ min(ASTDTM, na.rm=TRUE),
                      TRUE ~ ASTDTM),
    AENDTM= case_when(PARCAT1=="OVERALL" ~ max(coalesce(AENDTM,ASTDTM), na.rm=TRUE),
                      TRUE~AENDTM),
    ASTDT=date(ASTDTM),
    AENDT-date(AENDTM)
  )

  #Dose
  derive_summary_records(
    by_vars = vars(STUDYID, USUBJID, EXTRT),
    filter_rows = PARAMCD=="DOSE",
    fns = list(AVAL ~ sum(., na.rm = TRUE)),
    set_values_to = vars(
      PARCAT1 = rep("OVERALL", 1),
      PARAMCD = c("TDOSE")
    ),
    drop_values_from = vars(EXDOSU, EXDOSFRM, EXDOSFRQ, EXROUTE, EXDURU)
  )  %>%
  #Dose ad?
  derive_summary_records(
    by_vars = vars(STUDYID, USUBJID, EXTRT),
    filter_rows = PARAMCD=="ADJ",
    fns = list(AVALC ~ if_else(sum(!is.na(.)) > 0, "Y", NA_character_)),
    set_values_to = vars(
      PARCAT1 = rep("OVERALL", 1),
      PARAMCD = c("TADJ")
    ),
    drop_values_from = vars(EXDOSU, EXDOSFRM, EXDOSFRQ, EXROUTE, EXDURU)
  )   %>%
  #Dose adj due to AE
  derive_summary_records(
    by_vars = vars(STUDYID, USUBJID, EXTRT),
    filter_rows = PARAMCD=="ADJAE",
    fns = list(AVALC ~ if_else(sum(!is.na(.)) > 0, "Y", NA_character_)),
    set_values_to = vars(
      PARCAT1 = rep("OVERALL", 1),
      PARAMCD = c("TADJAE")
    ),
    drop_values_from = vars(EXDOSU, EXDOSFRM, EXDOSFRQ, EXROUTE, EXDURU)
  )   %>%
  #Exposure at W2-W24
  #Duration
  derive_summary_records(
    by_vars = vars(STUDYID, USUBJID, EXTRT),
    filter_rows = PARAMCD=="DOSE" & VISIT %in% c("WEEK 2", "WEEK 24"),
    fns = list(AVAL ~ sum(., na.rm = TRUE) ),
    set_values_to = vars(
      PARCAT1 = rep("W2-24", 1),
      PARAMCD = c("PDOSE")
    ),
    drop_values_from = vars(EXDOSU, EXDOSFRM, EXDOSFRQ, EXROUTE, EXDURU)
  )
  # repeat for other parameters
  # ...


#Part 4
  #Derive ASTDT/AENDT for summary param
  #date for overall param
  dates_overall<-adex1 %>%
    group_by(STUDYID, USUBJID, EXTRT) %>%
    summarise(
      ASTDTM = min(ASTDTM, na.rm = TRUE),
      AENDTM = max(coalesce(AENDTM, ASTDTM), na.rm = TRUE),
      ASTDT=date(ASTDTM),
      AENDT=date(AENDTM),
      PARCAT1="OVERALL"
  )
  #Dates for W2-W24 param
  dates_w224<-adex1 %>%
    filter(VISIT %in%  c("WEEK 2", "WEEK 24"))%>%
    group_by(STUDYID, USUBJID, EXTRT) %>%
    summarise(
      ASTDTM = min(ASTDTM, na.rm = TRUE),
      AENDTM = max(coalesce(AENDTM, ASTDTM), na.rm = TRUE),
      ASTDT=date(ASTDTM),
      AENDT=date(AENDTM),
      PARCAT1="W2-24"
    )


#Add the date to the summary param
adex<-adex1 %>%
  left_join(dates_overall, by= c("STUDYID", "USUBJID", "EXTRT","PARCAT1")) %>%
  mutate(ASTDTM=coalesce(ASTDTM.x, ASTDTM.y),
         ASTDT=coalesce(ASTDT.x, ASTDT.y),
         AENDTM=coalesce(AENDTM.x, AENDTM.y),
         AENDT=coalesce(AENDT.x, AENDT.y)
         ) %>%
  select(-ends_with(".x"), -ends_with(".y")) %>%
  left_join(dates_w224, by= c("STUDYID", "USUBJID", "EXTRT","PARCAT1")) %>%
    mutate(ASTDTM=coalesce(ASTDTM.x, ASTDTM.y),
           ASTDT=coalesce(ASTDT.x, ASTDT.y),
           AENDTM=coalesce(AENDTM.x, AENDTM.y),
           AENDT=coalesce(AENDT.x, AENDT.y),
           PARCAT1.x=PARCAT1,
           PARCAT1=if_else(is.na( PARCAT1.x), "INDIVIDUAL",  PARCAT1.x)
    ) %>%
    select(-ends_with(".x"), -ends_with(".y"))%>%

  # Add PARAMN and PARAM
  left_join(param_lookup, by = "PARAMCD") %>%
  # Derive AVALCATx
  # Note: Derivation of AVALCAT is not represented in the CDISC Pilot. It is
  #       presented for demonstration purposes.
  mutate(AVALCAT1 = format_avalcat1(param = PARAMCD, aval = AVAL)) %>%
  # Calculate ASEQ
  derive_obs_number(
    new_var = ASEQ,
    by_vars = vars(STUDYID, USUBJID),
    order = vars( ASTDT, VISIT, VISITNUM,EXSEQ,PARAMN),
    check_type = "warning"
  )



# Final Steps, Select final variables and Add labels
# This process will be based on your metadata, no example given for this reason
# ...

# ---- Save output ----

save(adex, file = "/PATH/TO/SAVE/ADex", compress = TRUE)
