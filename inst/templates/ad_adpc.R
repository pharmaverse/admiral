# Name: ADPC
#
# Label: Pharmacokinetics Concentrations Analysis Dataset
#
# Description: Based on simulated data, create ADPC analysis dataset
#
# Input: pc, adsl

# Clean environment
rm(list=ls(all=TRUE))

library(admiral)
library(dplyr)
library(lubridate)
library(stringr)
library(tidylog)
library(admiral.test) # Contains example datasets from the CDISC pilot project or simulated

# ---- Load source datasets ----

# Use e.g. haven::read_sas to read in .sas7bdat, or other suitable functions
# as needed and assign to the variables below.
# For illustration purposes read in admiral test data

# Load PC Adsl and EX
data("admiral_pc")
data("admiral_ex")
data("admiral_adsl")

# When SAS datasets are imported into R using haven::read_sas(), missing
# character values from SAS appear as "" characters in R, instead of appearing
# as NA values. Further details can be obtained via the following link:
# https://pharmaverse.github.io/admiral/articles/admiral.html#handling-of-missing-values

# add missing T in PCDTC
admiral_pc$PCDTC <-sub(' ','T',admiral_pc$PCDTC)

pc <- convert_blanks_to_na(admiral_pc)


# ---- Lookup tables ----
param_lookup <- tibble::tribble(
  ~PCTESTCD, ~PARAMCD, ~PARAM, ~PARAMN,
  "XAN", "PKCONC", "Pharmacokinetic concentration of Xanomeline", 1,
)

attr(param_lookup$PCTESTCD, "label") <- "Pharmacokinetic Test Short Name"

#remove placebo dose administration
ex <- filter(admiral_ex,EXDOSE>0)

#Check different treatment name
distinct(select(ex,EXTRT))

# Patient with Dose data
doseid = distinct(select(ex,USUBJID)) #168 patients

# Patient with PK data
pkid = distinct(select(pc,USUBJID)) #254 patients

# Restrict PK data for patients with Dose data
pcdose=dplyr::inner_join(pc,doseid,by='USUBJID')
pkdoseid = distinct(select(pcdose,USUBJID))  #168 patients

# Restrict PK data for patients without Dose data. Check that all PK samples are BLQ
pcnodose=dplyr::anti_join(pc,doseid,by='USUBJID')
pknodoseid = distinct(select(pcnodose,USUBJID)) #86 patients
qcnopk = distinct(select(pcnodose,PCSTRESC)) # check if we have non BLQ data

#replace missing end date by start date
ex <-ex %>% mutate (EXENDTC=case_when(EXENDTC !="" ~ EXENDTC, EXENDTC =="" ~ EXSTDTC))

# Convert EX start date to numeric
ex2<-derive_vars_dt(dataset=ex,new_vars_prefix="AST",dtc=EXSTDTC)

# Convert EX end date to numeric
ex3<-derive_vars_dt(dataset=ex2,new_vars_prefix="AEN",dtc=EXENDTC)

# add missing dose clock time
ex4<-mutate(ex3,ASTDTC=paste(ex3$ASTDT,"T01:00:00"))
ex4a<-mutate(ex4,ASTDTC=str_remove(ASTDTC," "))
ex5<-mutate(ex4a,AENDTC=paste(ex3$AENDT,"T01:00:00"))
ex5a<-mutate(ex5,AENDTC=str_remove(AENDTC," "))

#Convert dose date/time to numeric
ex6<-derive_vars_dtm(dataset =ex5a,dtc = ASTDTC,new_vars_prefix = "AST",)
ex6a<-derive_vars_dtm(dataset =ex6,dtc = AENDTC,new_vars_prefix = "AEN",)

# Create 1 dose per date
ex7=create_single_dose_dataset(ex6a,start_date=ASTDT,end_date=AENDT,keep_source_vars = vars(STUDYID,USUBJID, EXDOSFRQ, ASTDT, ASTDTM, AENDT, AENDTM,EXSTDTC,EXENDTC,EXDOSE,EXTRT))
#ex7=create_single_dose_dataset(ex6a,start_date=ASTDT,end_date=AENDT,keep_source_vars = all_vars())

#Convert PK sampling date/time to numeric
pcdose2<-derive_vars_dtm(dataset = pcdose,dtc = PCDTC,new_vars_prefix = "PC",)

#Assign to each PK sample the corresponding dose date/time
NCA <- derive_var_last_dose_date(
  dataset=pcdose2,
  dataset_ex=ex7,
  dose_date=ASTDTM,
  analysis_date=PCDTM,
  new_var=PCRFDTM,
  single_dose_condition = (EXDOSFRQ == "ONCE"))

# Get list of ADSL
adsl_vars=select(admiral_adsl,USUBJID,AGE,AGEU,SEX,RACE,TRT01P)

NCA1=dplyr::left_join(NCA,adsl_vars,by=c('STUDYID','USUBJID'))

# reshape data
NCA2 = NCA %>%
  rename(AVALC=PCSTRESC) %>%
  rename(AVAL=PCSTRESN) %>%
  rename(VISITDY=PCDY) %>%
  rename(AVALU=PCSTRESU) %>%
  rename(PARAM=PCTEST) %>%
  rename(PARAMCD=PCTESTCD)

NCA2$AVISIT<-NCA2$VISIT
NCA2$AVISITN<-NCA2$VISITNUM
NCA2$RELTMU<-'Hrs'
NCA2$NRELTM2<-NCA2$PCTPTNUM
NCA2$ASTDTM<-NCA2$PCRFDTM

#Do not get why ARELTM2 disply sec
NCA2$ARELTM2<-(NCA2$PCDTM-NCA2$PCRFDTM)/60


# Get list of ADSL
adsl_vars=select(admiral_adsl,USUBJID,AGE,AGEU,SEX,RACE,TRT01P)


# Get list of ADSL
dose_cov=select(ex7,USUBJID,ASTDTM,EXDOSE)
dose_cov1=dose_cov %>%
  rename(ACTDOSE=EXDOSE)


NCA3=dplyr::left_join(NCA2,dose_cov1,by=c('USUBJID','ASTDTM'))

FDOSE=arrange(ex7,USUBJID,ASTDTM)
FDOSE1 <- FDOSE %>%
  group_by(USUBJID) %>%
  filter(row_number(USUBJID) == 1)
FDOSE2=select(FDOSE1,USUBJID,ASTDTM) %>%
  rename(FIRSTDT=ASTDTM)

NCA4=dplyr::left_join(NCA3,FDOSE2,by='USUBJID')

#Do not get why ARELTM1 disply sec
NCA4$ARELTM1<-(NCA4$PCDTM-NCA4$FIRSTDT)/60

NCA5=dplyr::left_join(NCA4,adsl_vars,by='USUBJID')


ADPC=select(NCA5,STUDYID,USUBJID,AGE,AGEU,SEX,RACE,TRT01P,VISITNUM,VISIT,VISITDY,PCDTC,PCTPT,PCTPTNUM,AVAL,AVALC,AVALU,PARAM,PARAMCD,AVISIT,AVISITN,RELTMU,ARELTM2,NRELTM2,ARELTM1)
ADPC1<-ADPC[order(ADPC$USUBJID)]




########################################################################################################################################################################

# ---- Derivations ----

# Get list of ADSL
adsl_vars <- vars(STUDYID, USUBJID, AGE, AGEU, SEX, RACE, TRT01P)

adpc1 <- pc %>%
  # Join ADSL with PC (need TRTSDT for ADY derivation)
  left_join(
    select(admiral_adsl, STUDYID, USUBJID, !!!adsl_vars),
    by = c("STUDYID", "USUBJID")
  ) %>%
  # Calculate ADT, ADY
  derive_vars_dt(
    new_vars_prefix = "A",
    dtc = PCDTC
  ) %>%
  derive_vars_dy(reference_date = TRTSDT, source_vars = vars(ADT))

adpc2 <- adpc1 %>%
  # Add PARAMCD only - add PARAM etc later
  left_join(
    select(param_lookup, PCTESTCD, PARAMCD),
    by = "PCTESTCD"
  ) %>%
  # Calculate AVAL and AVALC
  mutate(
    AVAL = PCSTRESN,
    AVALC = PCSTRESC
  ) %>%
  # Remove variables
  select(-PCSTRESN, -PCSTRESC) %>%
  # Add ASEQ
  mutate(
    SRCDOM = DOMAIN,
    SRCVAR = "SEQ",
    SRCSEQ = PCSEQ
  ) %>%
  select(-DOMAIN, -PCSEQ)

# get visit info
adpc3 <- adpc2 %>%
  # Derive Timing
  mutate(
    AVISIT = case_when(
      str_detect(VISIT, "SCREEN|UNSCHED|RETRIEVAL|AMBUL") ~ NA_character_,
      !is.na(VISIT) ~ str_to_title(VISIT),
      TRUE ~ NA_character_
    ),
    AVISITN = VISITNUM
  ) %>%
  # Assign TRTA, TRTP
  mutate(
    TRTP = TRT01P,
    TRTA = TRT01A
  ) %>%
  # Derive AVALCA1N and AVALCAT1
  mutate(AVALCA1N = format_avalcat1n(param = PARAMCD, aval = AVAL)) %>%
  derive_vars_merged(dataset_add = avalcat_lookup, by_vars = vars(PARAMCD, AVALCA1N))

# Add all ADSL variables
adpc <- adpc3 %>%
  left_join(select(admiral_adsl, !!!admiral:::negate_vars(adsl_vars)),
    by = c("STUDYID", "USUBJID")
  )

# Final Steps, Select final variables and Add labels
# This process will be based on your metadata, no example given for this reason
# ...
# ---- Save output ----

dir <- tempdir() # Change to whichever directory you want to save the dataset in
save(adpc, file = file.path(dir, "adpc.rda"), compress = "bzip2")


# ---- User defined functions ----

# Here are some examples of how you can create your own functions that
#  operates on vectors, which can be used in `mutate`.
format_avalcat1n <- function(param, aval) {
  case_when(
    param == "PKCONC" & aval < 1 ~ 1,
    param == "PKCONC" & aval >= 1 ~ 2,
    T ~ NA_real_
  )
}
# ASSIGN AVALCAT1
avalcat_lookup <- tibble::tribble(
  ~PARAMCD, ~AVALCA1N, ~AVALCAT1,
  "PKCONC", 1, "< 1",
  "PKCONC", 2, ">= 1"
)
