# Name: ADAB
#
# Label: Anti-Drug Antibody Analysis Dataset
#
# Description: Based on simulated data, create ADAB analysis dataset
#
# Input: is_ada, ex, adsl
library(admiral)
library(dplyr)
library(lubridate)
library(stringr)
library(pharmaversesdtm) # Contains example datasets from the CDISC pilot project or simulated

# ---- Load source datasets ----

# Use e.g. haven::read_sas to read in .sas7bdat, or other suitable functions
# as needed and assign to the variables below.
# For illustration purposes read in admiral test data

# Load IS, EX and ADSL
# Note: these two lines will be udpated when pharmaversesdtm::is_ada is online
#is <- pharmaversesdtm::is_ada
is <- arrow::read_parquet("data/pharmaverse/is_ada.parquet")
ex <- pharmaversesdtm::ex
adsl <- admiral::admiral_adsl

# When SAS datasets are imported into R using haven::read_sas(), missing
# character values from SAS appear as "" characters in R, instead of appearing
# as NA values. Further details can be obtained via the following link:
# https://pharmaverse.github.io/admiral/cran-release/articles/admiral.html#handling-of-missing-values # nolint

ex <- convert_blanks_to_na(ex)
is <- convert_blanks_to_na(is)

# Define values for records with overall values
# Suggested are AVISIT=Overall, AVISITN=11111
overall_avisit <- "OVERALL"
overall_avisitn <- 11111

# ---- Derivations ----

is_dates <- is %>%
  # Filter as needed
  filter( toupper(ISBDAGNT) == "XANOMELINE") %>%
  # Join ADSL with is (need TRTSDT for ADY derivation)
  derive_vars_merged(
    dataset_add = adsl,
    new_vars =  exprs(TRTSDT),
    by_vars = exprs(STUDYID, USUBJID)
  ) %>%
  # Derive analysis date/time
  # Impute missing time to 00:00:00
  derive_vars_dtm(
    new_vars_prefix = "A",
    highest_imputation = "s",
    dtc = ISDTC,
    time_imputation = "00:00:00"
  ) %>%
  # Derive dates and times from date/times
  derive_vars_dtm_to_dt(exprs(ADTM)) %>%
  derive_vars_dtm_to_tm(exprs(ADTM)) %>%
  derive_vars_dy(reference_date = TRTSDT, source_vars = exprs(ADT)) %>%
  # Derive DRUG, ADATYPE nominal relative time from first dose (NFRLT)
  mutate(
    # Set AVISIT and AVISITN based on VISIT and VISITNUM
    AVISIT = VISIT,
    AVISITN = VISITNUM,
    # Map the analyte test to matching DRUG based on EX
    #  This is especially critical when multipel analytes and EX.EXTRT instances
    DRUG = case_when(
      toupper(ISBDAGNT) == "XANOMELINE" ~ "XANOMELINE",
    ),
    # Set ADATYPE based on SDTM V2.0 ISTESTCD, can also assign using other
    #   values in local version of IS.
    ADATYPE = ISTESTCD,

    # Set ADAPARM based on ISBDAGNT or other values
    #   ADAPARM will later become PARAMCD
    ADAPARM =  ISBDAGNT,

    # Assing NFRLT basedon how pharmaversesdtm.IS_ADA is structured or study specific
    #   This example maps any discontinuation visits as 77777 and unscheduled as 99999
    #   Units for this ADAB sample will be Days.  ISTPTNUM is in hours, VISITDY is days
    NFRLT = case_when (
      grepl("EARLY DISC", toupper(VISIT)) ~ 77777,
      grepl("TREATMENT DISC", toupper(VISIT)) ~ 77777,
      grepl("UNSCHEDULED", toupper(VISIT)) ~ 99999,
      grepl("BASELINE", toupper(VISIT)) & !is.na(ISTPTNUM) ~ ISTPTNUM  / 24,
      !is.na(VISITDY) ~ (VISITDY - 1)
    ),
    FRLTU = "DAYS"
  )


# ---- Get dosing information ----

ex_dates <- ex %>%
  # Add ADSL Vaiables for ADY
  derive_vars_merged(
    dataset_add = adsl,
    new_vars = exprs(TRTSDT),
    by_vars = exprs(STUDYID, USUBJID)
  ) %>%
  # Keep appliable desired records based on EXTRT and/or dose values (>=0, >0, etc.)
  filter(grepl("XANOMELINE", toupper(EXTRT)) |  grepl("PLACEBO", toupper(EXTRT)), EXDOSE >= 0) %>%
  # Add analysis datetime variables and set missing end date to start date
  # Impute missing time to 00:00:00
  # Note all times are missing for dosing records in this example data
  derive_vars_dtm(
    new_vars_prefix = "AST",
    dtc = EXSTDTC,
    time_imputation = "00:00:00"
  ) %>%
  derive_vars_dtm(
    new_vars_prefix = "AEN",
    dtc = EXENDTC,
    time_imputation = "00:00:00"
  ) %>%
  # Derive dates from date/times
  derive_vars_dtm_to_dt(exprs(ASTDTM)) %>%
  derive_vars_dtm_to_dt(exprs(AENDTM)) %>%
  # Set missing end dates to start date
  mutate(
    AENDTM = case_when(
    is.na(AENDTM) ~ ASTDTM,
    TRUE ~ AENDTM
   )
  ) %>%
  derive_vars_dtm_to_dt(exprs(ASTDTM)) %>%
  derive_vars_dy(reference_date = TRTSDT, source_vars = exprs(ASTDT))

 # This template computes only AFRLT, see ADPC template if need EX dose expansion example.

ex_exp <- ex_dates %>%
  mutate(
    # Map EXTRT to DRUG to match DRUG values in is_dates above
    #  This will be used to merge first dose into IS working data.
    DRUG = case_when(
      grepl("XANOMELINE", toupper(EXTRT)) ~  "XANOMELINE",
      grepl("PLACEBO", toupper(EXTRT)) ~  "XANOMELINE",
    ),

    # Compute Nominal Time
    NFRLT = case_when(
      VISITDY == 1 ~ 0,
      TRUE ~ VISITDY
    )
  )

# Derive AFRLT in IS data

is_afrlt <- is_dates %>%
  derive_vars_merged(
    dataset_add = ex_exp,
    filter_add = (EXDOSE >= 0 & !is.na(ASTDTM)),
    new_vars = exprs(FANLDTM = ASTDTM, FANLTMF = ASTTMF),
    order = exprs(ASTDTM, EXSEQ),
    mode = "first",
    by_vars = exprs(STUDYID, USUBJID, DRUG)
  )  %>%
  derive_vars_dtm_to_dt(exprs(FANLDTM)) %>%
  derive_vars_dtm_to_tm(exprs(FANLDTM)) %>%
  derive_vars_duration(
    new_var = AFRLT,
    start_date = FANLDTM,
    end_date = ADTM,
    out_unit = "days",
    floor_in = FALSE,
    add_one = FALSE
  ) %>%
  mutate (
    ACTTMFDS = as.numeric(difftime(ADTM, FANLDTM, units = "hours")),
    ACTTMFDS = case_when(
      FRLTU == "DAYS" ~ ACTTMFDS / 24,
      FRLTU == "HOURS" ~ ACTTMFDS
    )
  )  %>%
  arrange(STUDYID, USUBJID, DRUG, ADATYPE, ADAPARM, ADTM, NFRLT)


# Compute or assign BASETYPE, APERIOD and APHASE ----------------------------------
# Add study specific code as applicable using ADEX or ADSL APx and PHx vars

is_basetype <- is_afrlt %>%
  mutate(
    APERIOD = 1,
    APERIODC = "Period 01",
    APHASE = NA_character_,
    APHASEN = NA_integer_,
    BASETYPE = "DOUBLE_BLINDED"
  )

  # Review Unique AAB VISITs and NOMTIMES data
  is_basetype  %>%
    group_by(DRUG, ISTESTCD, ADATYPE, ISBDAGNT, VISIT, AVISIT, AVISITN, ISTPT, NFRLT) %>%
    summarize (n = n()) %>%
    print(n = 100)

  # Assign AVAL, AVALC, AVALU and DTYPE for each ISTESTCD and ISBDAGNT

  is_aval <- is_basetype %>%
    mutate(
      MRT = case_when(
        ADATYPE == "ADA_BAB" & ISBDAGNT == "XANOMELINE" ~ 1.4,
        TRUE ~ NA_real_
      ),
      DTL = case_when(
        ADATYPE == "ADA_BAB" & ISBDAGNT == "XANOMELINE" ~ 999,
        TRUE ~ NA_real_
      ),
      RESULTC = case_when(
        toupper(ISSTRESC) %in% c(
          "NEGATIVE", "NEGATIVE SCREEN", "NEGATIVE IMMUNODEPLETION",
          "NEGATIVE CONFIRMATION"
        ) ~ "NEGATIVE",
        toupper(ISSTRESC) %in% c(
          "NEGATIVE TITER", "<1.70", "< 1.70", "<1.30", "< 1.30",  "<1.40", "POSITIVE IMMUNODEPLETION", "POSITIVE CONFIRMATION",
          "POSITIVE"
        ) ~ "POSITIVE",
        ISSTRESN > 0 ~ "POSITIVE",
        TRUE ~ NA_character_
      ),
      RESULTN = case_when(
        toupper(RESULTC) == "POSITIVE" ~ 1,
        toupper(RESULTC) == "NEGATIVE" ~ 0,
        TRUE ~ NA_integer_
      ),
      AVAL = case_when(
        ADATYPE == "ADA_BAB" & toupper(RESULTC) == "POSITIVE" & !is.na(ISSTRESN) ~ ISSTRESN,
        ADATYPE == "ADA_BAB" & toupper(RESULTC) == "POSITIVE" & is.na(ISSTRESN) & !is.na(MRT) ~ MRT,
      ),
      AVALC = case_when(
        # NABSTAT gets ISSTRESC, Standard ADA is special rules when non NEGATIVE results
        ADATYPE == "ADA_NAB" ~ ISSTRESC,
        toupper(RESULTC) == "POSITIVE" & ADATYPE == "ADA_BAB" & !is.na(ISSTRESN) ~
          as.character(ISSTRESN),
        toupper(RESULTC) == "POSITIVE" & ADATYPE == "ADA_BAB" & is.na(ISSTRESN) &
          (grepl("<", ISSTRESC) | grepl("NEGATIVE TITER", toupper(ISSTRESC))) & !is.na(MRT) ~
          paste("<", as.character(MRT), sep = ""),
        # If positive with no numeric ISSTRESN set to N.C.
        toupper(RESULTC) == "POSITIVE" & ADATYPE == "ADA_BAB" ~ "N.C.",
        toupper(RESULTC) == "NEGATIVE" & ADATYPE == "ADA_BAB" ~ NA_character_,
        TRUE ~ NA_character_
      ),
      AVALU = case_when(
        ADATYPE == "ADA_NAB" ~ ISSTRESU,
        ADATYPE == "ADA_BAB" & toupper(RESULTC) == "POSITIVE" & !is.na(ISSTRESN) ~ ISSTRESU,
        ADATYPE == "ADA_BAB" & toupper(RESULTC) == "POSITIVE" & is.na(ISSTRESN) & !is.na(MRT) ~ "titer",
        TRUE ~ NA_character_
      ),
      DTYPE = case_when(
        ADATYPE == "ADA_BAB" & toupper(RESULTC) == "POSITIVE" & is.na(ISSTRESN) & !is.na(MRT) ~ "MRT",
        TRUE ~ NA_character_
      )
    )

# Begin computation of parameters -----------------------------------------

  # Identify Best Baseline for each analyte and parameter type
  # Baseline is NFRLT <= 0 or Unscheduled AND the ADA Date is on or before the date of first dose.

  is_baseline <- is_aval %>%
    # Calculate ABLFL, if more than one record for the 'order' and 'filter', will throw a duplicate record warning
    #   User can decide how to adjust the order, filter or prior data step.
    restrict_derivation(
      derivation = derive_var_extreme_flag,
      args = params(
        by_vars = exprs(STUDYID, USUBJID, BASETYPE, ADATYPE, ADAPARM),
        order = exprs(ADT, NFRLT),
        new_var = ABLFL,
        mode = "last"
      ),
      # flag baseline based on time not values due to ADA parameters can in some cases permit missing.
      filter = (NFRLT <= 0 & (ADT <= FANLDT) | (NFRLT > 60000 & (ADT <= FANLDT)) & !is.na(BASETYPE))
    ) %>%
    mutate(
      # VALID flags for use later:
      # VALIDBASE flags non-missing values on baseline (by each ADATYPE and ADAPARM)
      #   Note: VALIDBASE is not used as thie templates allows a baseline to be valid as
      #        as long as its present (can be missing), adapt as needed.
      # VALIDPOST flags non-missing values on post-baseline (by each ADATYPE and ADAPARM)
      VALIDBASE = case_when(
        ABLFL == 'Y' & (!is.na(AVALC) | !is.na(RESULTC) | !is.na(AVAL)) ~ "Y",
      ),
      VALIDPOST = case_when(
        ADTM > FANLDTM & is.na(ABLFL) & (!is.na(RESULTC) | !is.na(AVAL)) ~ "Y",
        ADTM > FANLDTM & is.na(ABLFL) & (is.na(RESULTC) & is.na(AVAL)) ~ "N",
        TRUE ~ NA_character_
      )
    ) %>%
    arrange(STUDYID, USUBJID, DRUG, ADATYPE, ADAPARM, ADTM, NFRLT)

   check_baseline <- is_baseline %>%
     select (STUDYID, USUBJID, BASETYPE, ADATYPE, ADAPARM, AVISIT, ADT, NFRLT, AVAL, AVAL, AVALC, RESULTC, RESULTN, ABLFL, VALIDBASE, VALIDPOST)

  # Compute BASE and CHG
  is_aval_change <- is_baseline %>%
    derive_var_base(
      by_vars = exprs(!!!get_admiral_option("subject_keys"), BASETYPE, ADATYPE, ADAPARM),
      source_var = AVAL,
      new_var = BASE,
      filter = ABLFL == "Y"
    )

  # If not ABLFL set CHG to NA
  is_aval_change <- derive_var_chg(is_aval_change) %>%
    mutate(
      CHG = case_when(
        ABLFL == "Y" ~ NA_integer_,
        TRUE ~ CHG
      )
    )

  check_aval_change <- is_aval_change %>%
    select (USUBJID, BASETYPE, ADATYPE, ADAPARM, AVISIT, ADT, NFRLT, AVAL, AVAL, AVALC, RESULTC, RESULTN, ABLFL, VALIDBASE, VALIDPOST, AVAL, BASE, CHG)


  # Interpreted RESULT BASE
  is_result_change <- is_aval_change %>%
    derive_var_base(
      by_vars = exprs(!!!get_admiral_option("subject_keys"), BASETYPE, ADATYPE, ADAPARM),
      source_var = RESULTN,
      new_var = BASE_RESULT,
      filter = ABLFL == "Y"
    )

  check_result_change <- is_result_change %>%
    select (USUBJID, BASETYPE, ADATYPE, ADAPARM, AVISIT, ADT, NFRLT, AVAL, AVAL, AVALC, RESULTC, RESULTN, ABLFL, VALIDBASE, VALIDPOST, BASE_RESULT, BASE, CHG)

  # Get base only data for use later
  base_data <- is_result_change %>%
    filter(ABLFL == "Y") %>%
    select(
      !!!get_admiral_option("subject_keys"), DRUG, BASETYPE, ADATYPE, ADAPARM, BASE_RESULT, BASE,
      ABLFL
    )

  # ADABLPFL will flag subjects that had a baseline record on main ADA_BAB analytes
  # ADPBLPFL will flag subjects that had a valid (non missing) post baseline record on main ADA_BAB analytes


  # For later use (if just need to merge in the adablpfl flag by paramcd)
  adablpfl <- is_result_change %>%
    filter(ABLFL == "Y") %>%
    distinct(!!!get_admiral_option("subject_keys"), DRUG, BASETYPE, ADATYPE, ADAPARM, .keep_all = TRUE) %>%
    select(!!!get_admiral_option("subject_keys"), DRUG, BASETYPE, ADATYPE, ADAPARM, BASE_RESULT, BASE, ABLFL) %>%
    rename(ADABLPFL = ABLFL)


  # Add Base_data then Compute CHG and By Visit ADA Flags
  # Note: sample assumes one BASETYPE, if multiple APERIOD,
  #    Add additional code per specs as needed.

  is_visit_flags <- is_result_change %>%
    mutate(
      TFLAGV = case_when(

        VALIDPOST == "Y" & is.na(ABLFL) & ADATYPE == "ADA_BAB" & (BASE_RESULT == 1 & (CHG >= 0.6)) ~ 2,

        VALIDPOST == "Y" & is.na(ABLFL) & ADATYPE == "ADA_BAB" &  BASE_RESULT == 1 & (((CHG < 0.6) & !is.na(AVAL)) | (RESULTN == 0)) ~ 3,

        VALIDPOST == "Y" & is.na(ABLFL) & ADATYPE == "ADA_BAB" & ((BASE_RESULT == 0 | is.na(BASE_RESULT)) & RESULTN == 0) ~ 0,

        VALIDPOST == "Y" & is.na(ABLFL) & ADATYPE == "ADA_BAB" & ((BASE_RESULT == 0 | is.na(BASE_RESULT)) & RESULTN == 1) ~ 1
      ),
      PBFLAGV = case_when(
        !is.na(TFLAGV) & TFLAGV %in% c(1, 2) ~ 1,
        !is.na(TFLAGV) & TFLAGV %in% c(0, 3) ~ 0
      ),
      ADASTATV = case_when(
        !is.na(PBFLAGV) & PBFLAGV == 1 ~ "ADA+",
        !is.na(PBFLAGV) & PBFLAGV == 0 ~ "ADA-",
        VALIDPOST == "Y" & is.na(ABLFL)  & ADATYPE == "ADA_BAB" & is.na(PBFLAGV) ~ "MISSING",
        TRUE ~ "MISSING"
      ),
    )


  # Postbaseline must be valid post data (result not missing)
  post_data <- is_visit_flags %>%
    filter(VALIDPOST == "Y") %>%
    select(
      !!!get_admiral_option("subject_keys"), DRUG, BASETYPE, ADATYPE, ADAPARM, RESULTN, AVAL, ADTM,
      CHG, VALIDPOST
    ) %>%
    rename(AVAL_P = AVAL) %>%
    rename(RESULT_P = RESULTN)



  adpblpfl <- post_data %>%
    distinct(!!!get_admiral_option("subject_keys"), BASETYPE, ADATYPE, ADAPARM, .keep_all = TRUE) %>%
    select(!!!get_admiral_option("subject_keys"), BASETYPE, ADATYPE, ADAPARM, RESULT_P, AVAL_P, VALIDPOST) %>%
    rename(ADPBLPFL = VALIDPOST)

  post_data <- post_data %>% select(-VALIDPOST)



  # Compute BFLAG, TFLAG, PBFLAG ------------------


  most_post_result <- post_data %>%
    group_by(!!!get_admiral_option("subject_keys"), DRUG, BASETYPE, ADATYPE, ADAPARM) %>%
    summarize(RESULT_P = max(RESULT_P)) %>%
    ungroup()


  most_post_aval <- post_data %>%
    group_by(!!!get_admiral_option("subject_keys"), BASETYPE, ADATYPE, ADAPARM) %>%
    filter(!is.na(AVAL_P)) %>%
    summarize(AVAL_P = max(AVAL_P)) %>%
    ungroup()

  most_post_chg <- post_data %>%
    group_by(!!!get_admiral_option("subject_keys"), BASETYPE, ADATYPE, ADAPARM) %>%
    filter(!is.na(AVAL_P)) %>%
    summarize(MAXCHG = max(CHG)) %>%
    ungroup()


  most_post <- most_post_result %>%
    derive_vars_merged(
      dataset_add = most_post_aval,
      by_vars = exprs(!!!get_admiral_option("subject_keys"), BASETYPE, ADATYPE, ADAPARM)
    )

  most_post <- most_post %>%
    derive_vars_merged(
      dataset_add = most_post_chg,
      by_vars = exprs(!!!get_admiral_option("subject_keys"), BASETYPE, ADATYPE, ADAPARM)
    )



  # Use outer Join to combine baseline with most post results
  flag_data <- full_join(base_data, most_post, by = c("DRUG", "STUDYID", "USUBJID", "BASETYPE", "ADATYPE", "ADAPARM")
  ) %>%
    arrange(DRUG, BASETYPE, ADATYPE, ADAPARM, !!!get_admiral_option("subject_keys"), )


  # Set the values for AVAL usage (numeric) then later convert to AVALC when creating parameters
  flag_data <- flag_data %>%
    mutate(
      BFLAG =
        case_when(
          BASE_RESULT == 0 ~ 0,
          BASE_RESULT == 1 ~ 1
        ),
      TFLAG =
        case_when(
          (BFLAG == 0 | is.na(BFLAG)) & RESULT_P == 0 ~ 0,
          (BFLAG == 0 | is.na(BFLAG)) & RESULT_P == 1 ~ 1,
          BFLAG == 1 & (MAXCHG >= 0.6) ~ 2,
          BFLAG == 1 & !is.na(MAXCHG) & (MAXCHG < 0.6) ~ 3,
          BFLAG == 1 & is.na(MAXCHG) & RESULT_P == 0  ~ 3
        ),
      PBFLAG =
        case_when(
          ADATYPE == "ADA_BAB" & (TFLAG == 1 | TFLAG == 2) ~ 1,
          ADATYPE == "ADA_BAB" & (TFLAG == 0 | TFLAG == 3) ~ 0,
          ADATYPE == "ADA_NAB" & RESULT_P == 0 ~ 0,
          ADATYPE == "ADA_NAB" & RESULT_P == 1 ~ 1
        ),
      ADASTAT = case_when(
        ADATYPE == "ADA_BAB" & PBFLAG == 1 ~ 1,
        ADATYPE == "ADA_BAB" & PBFLAG == 0 ~ 0
      )
    )


  # Pull out Main ADA Titer ADASTAT status for NABSTAT option 1 as ADASTAT_MAIN

  flag_data <- flag_data %>%
    derive_vars_merged(
      dataset_add = flag_data,
      filter_add = ADATYPE == "ADA_BAB",
      new_vars = exprs(ADASTAT_MAIN = ADASTAT),
      by_vars = exprs(!!!get_admiral_option("subject_keys"), DRUG)
    )


  # Compute NABPOSTMISS onto flag_data for NABSTAT

  flag_data <- flag_data %>%
    derive_var_merged_exist_flag(
      dataset_add = is_visit_flags,
      new_var = NABPOSTMISS,
      condition = is.na(ABLFL) & VALIDPOST == "N" & ADATYPE == "ADA_NAB",
      by_vars = exprs(!!!get_admiral_option("subject_keys"), BASETYPE, ADATYPE, ADAPARM)
    )


  # Compute NAB Stat using both methods
  # Note: Options 2 added as a placekeeper, starter method, adjust as needed for
  #   current or study specific specs.

  flag_data <- flag_data %>%
    mutate(

      # For Option 1, if any post baseline NAB were blank without a Positive (NABPOSTMISS = Y), NABSTAT = missing
      nabstat_opt1 = case_when(
        # Based on ADASTAT (EMERNEG/EMERPOS) and NAB results
        ADATYPE == "ADA_NAB" & ADASTAT_MAIN == 1 & RESULT_P == 1 ~ 1,
        ADATYPE == "ADA_NAB" & ADASTAT_MAIN == 1 & RESULT_P == 0 &
                       (is.na(NABPOSTMISS) | NABPOSTMISS == "N")~ 0
      ),
      nabstat_opt2 = case_when(
        # Based on ADASTAT (EMERNEG/EMERPOS) Only.
        ADATYPE == "ADA_NAB" & ADASTAT_MAIN == 1 ~ 1,
        ADATYPE == "ADA_NAB" & ADASTAT_MAIN == 0 ~ 0,
      ),
    )



  # Drop the vars no longer neded from flag_data before merging with main ADAB
  flag_data <- flag_data %>% select(-BASE, -BASE_RESULT, -AVAL_P, -RESULT_P, -DRUG, -ABLFL)

  # Put TFLAG, BFLAG and PBFLAG onto the main ADAB dataset

  is_flagdata <- is_visit_flags %>%
   # main_aab_flagdata
    derive_vars_merged(
      dataset_add = flag_data,
      by_vars = exprs(!!!get_admiral_option("subject_keys"), BASETYPE, ADATYPE, ADAPARM)
    )

  # PERSADA and TRANADA

  per_tran_pre <- is_flagdata %>%
    filter(VALIDPOST == "Y") %>%
    select(!!!get_admiral_option("subject_keys"), BASETYPE, ADATYPE, ADAPARM, ADTM, FANLDTM, FANLDT, TFLAGV, ISSEQ) %>%
    group_by(!!!get_admiral_option("subject_keys"), BASETYPE, ADATYPE, ADAPARM,) %>%
    mutate(
      MaxADTM = max(ADTM)
    ) %>%
    ungroup() %>%
    mutate(
      LFLAG = case_when(
        ADTM == MaxADTM ~ 1
      ),
      LFLAGPOS = case_when(
        ADTM == MaxADTM & (TFLAGV == 1 | TFLAGV == 2) ~ 1
      )
    )


  # KD Update to template code: If last record has two or more dupes, later  keep the one with the higher ISSEQ
  per_tran <- per_tran_pre %>%
    filter(TFLAGV == 1 | TFLAGV == 2)   %>%
    group_by(!!!get_admiral_option("subject_keys"), BASETYPE, ADATYPE, ADAPARM, ADTM) %>%
    mutate(
      LastSEQ = max(ISSEQ),
    ) %>%
    ungroup()  %>%
    group_by(!!!get_admiral_option("subject_keys"), BASETYPE, ADATYPE, ADAPARM, ADTM) %>%
    mutate(
      PARM_N = n()
    ) %>%
    ungroup()  %>%
    mutate(
      drop_dupe = case_when(
        LFLAG  == 1 & ISSEQ != LastSEQ ~ TRUE,
        TRUE ~ FALSE
      )
    )

  # Potential last OBS dupe dates
  per_tran_dupes <- per_tran  %>%
    filter (PARM_N > 1 & LFLAG == 1)

  # if get dupe signal error below, uncomment this and evaluate the dupes
  per_tran <- per_tran %>%
    #filter (drop_dupe == FALSE)  %>%
    select(-MaxADTM, -LFLAG, -LastSEQ, -drop_dupe, - ISSEQ) %>%
    group_by(!!!get_admiral_option("subject_keys"), BASETYPE, ADATYPE, ADAPARM,) %>%
    mutate(
      FPPDTM = min(ADTM),
      LPPDTM = max(ADTM)
    ) %>%
    ungroup()

  per_tran_inc_last <- per_tran %>%
    group_by(!!!get_admiral_option("subject_keys"), BASETYPE, ADATYPE, ADAPARM,) %>%
    summarize(COUNT_INC_LAST = n()) %>%
    ungroup()


  per_tran_exc_last <- per_tran %>%
    filter(is.na(LFLAGPOS)) %>%
    group_by(!!!get_admiral_option("subject_keys"), BASETYPE, ADATYPE, ADAPARM,) %>%
    summarize(COUNT_EXC_LAST = n()) %>%
    ungroup()


  # Reduce per_tran to one record per PARAMCD, use ADTM = LPPDTM to get best last record

  per_tran <- per_tran %>%
    filter(ADTM == LPPDTM)


  per_tran <- per_tran %>%
    derive_vars_merged(
      dataset_add = per_tran_inc_last,
      by_vars = exprs(!!!get_admiral_option("subject_keys"), BASETYPE, ADATYPE, ADAPARM)
    )


  per_tran <- per_tran %>%
    derive_vars_merged(
      dataset_add = per_tran_exc_last,
      by_vars = exprs(!!!get_admiral_option("subject_keys"), BASETYPE, ADATYPE, ADAPARM)
    )


  per_tran <- per_tran %>%
    mutate(
      FPPDT = as_date(FPPDTM),
      LPPDT = as_date(LPPDTM),
      ADADUR = case_when(
        LPPDTM - FPPDTM == 0 & (!is.na(LPPDTM) & !is.na(FPPDTM)) ~ 1 / 7,
        !is.na(LPPDTM) & !is.na(FPPDTM)
        ~ ((as.numeric(difftime(LPPDTM, FPPDTM, units = "secs")) / (60 * 60 * 24)) + 1) / 7,
        !is.na(LPPDT) & !is.na(LPPDT)
        ~ (as.numeric(LPPDT - FPPDT) + 1) / 7
      ),
      TIMADA = case_when(
        !is.na(FPPDTM) & !is.na(FANLDTM) ~ as.numeric(difftime(FPPDTM, FANLDTM, units = "weeks")),
        !is.na(FPPDT) & !is.na(FANLDT) ~ (as.numeric(FPPDT - FANLDT) + 1) / 7
      ),
      tdur = case_when(
        !is.na(LPPDTM) & !is.na(LPPDTM) ~ as.numeric(difftime(LPPDTM, FPPDTM, units = "secs")) / (7 * 3600 * 24),
        !is.na(LPPDT) & !is.na(LPPDT) ~ as.numeric(LPPDT - FPPDT + 1) / 7
      ),
      TRANADA = case_when(
        TFLAGV == 1 & ((COUNT_EXC_LAST == 1 | (COUNT_INC_LAST >= 2 & tdur < 16)) & is.na(LFLAGPOS)) ~ 1
      ),
      PERSADA = case_when(
        TFLAGV == 1 & (COUNT_INC_LAST == 1 & (COUNT_EXC_LAST <= 0 | is.na(COUNT_EXC_LAST)) |
                         (COUNT_INC_LAST >= 2 & tdur >= 16) | LFLAGPOS == 1) ~ 1
      ),
      TRANADAE = case_when(
        TFLAGV >= 1 & ((COUNT_EXC_LAST == 1 | (COUNT_INC_LAST >= 2 & tdur < 16)) & is.na(LFLAGPOS)) ~ 1
      ),
      PERSADAE = case_when(
        TFLAGV >= 1 & (COUNT_INC_LAST == 1 & (COUNT_EXC_LAST <= 0 | is.na(COUNT_EXC_LAST)) |
                         (COUNT_INC_LAST >= 2 & tdur >= 16) | LFLAGPOS == 1) ~ 1
      ),
      INDUCED = case_when(
        TFLAGV == 1 ~ "Y",
        TRUE ~ "N"
      ),
      ENHANCED = case_when(
        TFLAGV == 2 ~ "Y",
        TRUE ~ "N"
      )
    )

  # Drop temporary vars that do not need to be merged into main ADAB
  per_tran <- per_tran %>%
    select(-ADTM, -TFLAGV, -FANLDTM, -FANLDT, -LFLAGPOS, -FPPDT, -LPPDT)

  # Put PERSADA, TRANADA, TDUR, ADADUR onto main ADAB

  main_aab_pertran <- is_flagdata %>%
    derive_vars_merged(
      dataset_add = per_tran,
      by_vars = exprs(!!!get_admiral_option("subject_keys"), BASETYPE, ADATYPE, ADAPARM)
    )


  main_aab_rtimes <- main_aab_pertran %>%
    mutate(
      ATPT = ISTPT,
      ADADUR = round(ADADUR, digits = 4),
      TIMADA = round(TIMADA, digits = 4),
      PERSADA = case_when(
        is.na(PERSADA) ~ 0,
        TRUE ~ PERSADA
      ),
      PERSADAE = case_when(
        is.na(PERSADAE) ~ 0,
        TRUE ~ PERSADAE
      ),
      TRANADA = case_when(
        is.na(TRANADA) ~ 0,
        TRUE ~ TRANADA
      ),
      TRANADAE = case_when(
        is.na(TRANADAE) ~ 0,
        TRUE ~ TRANADAE
      ),
      NABSTAT = nabstat_opt1
    )


  # Put ADABLPFL and ADPBFPFL onto main dataset.  If These are not Y set to N
  # NOTE: These were computed by ISTESTCD, we only need them for and on the main
  # Analytes (Original and Interpeted reuslts observations by ISTESTCD and EXTRT)


  main_aab <- main_aab_rtimes %>%
    derive_vars_merged(
      dataset_add = adablpfl,
      new_vars = exprs(ADABLPFL),
      by_vars = exprs(!!!get_admiral_option("subject_keys"), BASETYPE, ADATYPE, ADAPARM)
    ) %>%
    derive_vars_merged(
      dataset_add = adpblpfl,
      new_vars = exprs(ADPBLPFL),
      by_vars = exprs(!!!get_admiral_option("subject_keys"), BASETYPE, ADATYPE, ADAPARM)
    ) %>%
    mutate(
      ADABLPFL = case_when(
        ADABLPFL == "Y" ~ "Y",
        TRUE ~ NA_character_
      ),
      ADPBLPFL = case_when(
        ADPBLPFL == "Y" ~ "Y",
        TRUE ~ NA_character_
      )
    )

  view_keys <- main_aab %>%
    distinct (ISTESTCD, ISTEST,  ISBDAGNT)  %>%
    arrange (ISTESTCD,  ISTEST,  ISBDAGNT)


  main_aab_review <- main_aab %>%
    select(
      !!!get_admiral_option("subject_keys"), VISIT, ISTESTCD, ISBDAGNT, DRUG, ADTM, ADT, ADY, BASETYPE, ADATYPE, ADAPARM,
       MRT, DTL, ABLFL, ADABLPFL, ADPBLPFL, VALIDBASE, VALIDPOST, NFRLT, AFRLT, RESULTC, RESULTN,
      ISSTRESC, ISSTRESN,AVALC,  AVAL, AVALU, BASE, BFLAG, TFLAG, PBFLAG, ADASTAT, ADASTAT_MAIN, NABSTAT, nabstat_opt1,
      tdur, ADADUR, TIMADA, TRANADA, PERSADA
    )  %>%
    arrange(USUBJID, ISTESTCD, ADY)   %>%
    filter (1==1) %>%
    filter (ADATYPE=="ADA_NAB" & NABSTAT==0)
    #filter (USUBJID == "01-703-1439" | USUBJID == "01-704-1017" | USUBJID == "01-716-1063"  )





  # Begin Creation of each PARAM for the final ADAB format using main_aab ---------

  # Note:  Specs for assigning AVALC from AVAL and vice versa may change or vary
  #   Adjust AVAL and AVALC for each below PARAMCD, PARAM
  #   as needed for current or study specific requirements.

  # Save as core_aab to be the input for all the param assemblies

  core_aab <- main_aab %>%
    mutate(
      # Assign original PARAM (--TEST) to PARCAT1 before customizing PARAM for parameters
      PARAM = case_when(
        ADATYPE == "ADA_BAB" ~ paste ("Anti-", ADAPARM, " Antibody", sep = ""),
        ADATYPE == "ADA_NAB" ~ paste ("Anti-", ADAPARM, " Neutralizing Antibody", sep = ""),
      ),
      PARAMCD = ADAPARM,
      PARCAT1 = PARAM
    )

  core_aab %>%
    distinct(ADATYPE, PARAMCD, PARAM, PARCAT1)

  # By Visit Primary ADA ISTESTCD Titer Results
  adab_titer <- core_aab %>%
    filter(ADATYPE == "ADA_BAB") %>%
    mutate(
      # For ADASTAT, append "Titer Units" to PARAM
      PARAM = paste(PARAM, "Titer Units", sep = " "),
    )


  # By Visit NAB ISTESTCD Results
  adab_nabvis <- core_aab %>%
    filter(ADATYPE == "ADA_NAB") %>%
    # These two flags, BASE and CHG are only kept on the primary ADA ISTESTCD test
    select(-BASE, -CHG, -ADABLPFL, -ADPBLPFL)


  # By Visit Main ADA Titer Interpreted RESULT data
  adab_result <- core_aab %>%
    filter(ADATYPE == "ADA_BAB") %>%
    mutate(
      PARAMCD = "RESULTy",
      PARAM = paste("ADA interpreted per sample result,", PARCAT1, sep = " "),
      AVALC = toupper(RESULTC),
      AVAL = case_when(
        AVALC == "NEGATIVE" ~ 0,
        AVALC == "POSITIVE" ~ 1
      ),
      AVALU = NA_character_
    ) %>%
    # DTYPE is only kept on the parent ISTESTCD param, recompute CHG and BASE for RESULTy
    select(-BASE, -CHG, -DTYPE)

  adab_result <- adab_result %>%
    derive_var_base(
      by_vars = exprs(!!!get_admiral_option("subject_keys"), BASETYPE, ADATYPE, ADAPARM),
      source_var = AVAL,
      new_var = BASE,
      filter = ABLFL == "Y"
    )

  adab_result <- derive_var_chg(adab_result) %>%
    mutate(
      CHG = case_when(
        ABLFL == "Y" ~ NA_integer_,
        TRUE ~ CHG
      )
    )

  adab_result_review <- adab_result %>%
    arrange (USUBJID, BASETYPE, ISTESTCD, PARAMCD, AVISIT, NFRLT)  %>%
    select(USUBJID, BASETYPE, ISTESTCD, PARAMCD, AVISIT, NFRLT, AVALC, AVAL, ABLFL, BASE, CHG)



  # By Visit NAB Interpreted RESULT data (keep all IS vars plus ____)
  adab_nabres <- core_aab %>%
    filter(ADATYPE == "ADA_NAB") %>%
    mutate(
      PARAMCD = "RESULTy",
      PARAM = paste("NAB interpreted per sample result,", PARCAT1, sep = " "),
      AVALC = toupper(RESULTC),
      AVAL = case_when(
        AVALC == "NEGATIVE" ~ 0,
        AVALC == "POSITIVE" ~ 1
      ),
      AVALU = NA_character_
    ) %>%
    # These two flags, BASE and CHG are only kept on the primary ADA test
    select(-BASE, -CHG, -ADABLPFL, -ADPBLPFL)


  # ---- Derive BASE and Calculate Change from Baseline on NAB RESULT records  ----

  adab_nabres <- adab_nabres %>%
    derive_var_base(
      by_vars = exprs(!!!get_admiral_option("subject_keys"), BASETYPE, ADATYPE, ADAPARM),
      source_var = AVAL,
      new_var = BASE,
      filter = ABLFL == "Y"
    )

  adab_nabres <- derive_var_chg(adab_nabres) %>%
    mutate(
      CHG = case_when(
        ABLFL == "Y" ~ NA_integer_,
        TRUE ~ CHG
      )
    )

  # By Visit Titer TFLAGV data -----------------------

  # Note: By Visit parameters do not keep CHG, MRT and DTL
  adab_tflagv <- core_aab %>%
    filter(ADATYPE == "ADA_BAB") %>%
    mutate(
      PARAMCD = "TFLAGV",
      PARAM = paste("Treatment related ADA by Visit,", PARCAT1, sep = " "),
      AVAL = TFLAGV,
      AVALC = case_when(
        is.na(AVAL) ~ NA_character_,
        TRUE ~ as.character(AVAL)
      ),
      AVALU = NA_character_
    ) %>%
    # Drop BASE, CHG, MRT and DTL and the two --FL flags (are only kept on the primary ADA test)
    select(-BASE, -CHG, -MRT, -DTL, -ADABLPFL, -ADPBLPFL)


  # By Visit Titer PBFLAGV data ---------------------
  # Note: By Visit parameters do not keep CHG, MRT and DTL
  adab_pbflagv <- core_aab %>%
    filter(ADATYPE == "ADA_BAB") %>%
    mutate(
      PARAMCD = "PBFLAGV",
      PARAM = paste("Post Baseline Pos/Neg by Visit,", PARCAT1, sep = " "),
      AVALC = case_when(
        PBFLAGV == "1" | PBFLAGV == "2" ~ "POSITIVE",
        PBFLAGV == "0" | PBFLAGV == "3" ~ "NEGATIVE",
        TRUE ~ "MISSING"
      ),
      AVAL = case_when(
        AVALC == "POSITIVE" ~ 1,
        AVALC == "NEGATIVE" ~ 0
      ),
      AVALU = NA_character_
    ) %>%
    # Drop BASE, CHG, MRT and DTL and the two --FL flags (are only kept on the primary ADA test)
    select(-BASE, -CHG, -MRT, -DTL, -ADABLPFL, -ADPBLPFL)


  # By Visit Titer ADASTATV data
  # Note: By Visit parameters do not keep CHG, MRT and DTL
  adab_adastatv <- core_aab %>%
    filter(ADATYPE == "ADA_BAB") %>%
    mutate(
      PARAMCD = "ADASTATV",
      PARAM = paste("ADA Status of a patient by Visit,", PARCAT1, sep = " "),
      AVALC = ADASTATV,
      AVAL = case_when(
        AVALC == "ADA+" ~ 1,
        AVALC == "ADA-" ~ 0
      ),
      AVALU = NA_character_
    ) %>%
    # Drop BASE, CHG, MRT and DTL and the two --FL flags (are only kept on the primary ADA test)
    select(-BASE, -CHG, -MRT, -DTL, -ADABLPFL, -ADPBLPFL)



  # Next below are the individual params
  # assign the AVISIT and AVISITN for individual params

  core_aab <- core_aab %>%
    mutate(
      AVISIT = overall_avisit,
      AVISITN = overall_avisitn
    )


  # Get Patient flag BFLAG
  adab_bflag <- core_aab %>%
    filter(ADATYPE == "ADA_BAB") %>%
    distinct(
      !!!get_admiral_option("subject_keys"), BASETYPE, PARAM, ADATYPE, ADAPARM, ISSPEC, ISTESTCD, ISTEST, ISBDAGNT,
      AVISITN, AVISIT, PARCAT1, BFLAG
    ) %>%
    mutate(
      PARAMCD = "BFLAGy",
      PARAM = paste("Baseline Pos/Neg,", PARCAT1, sep = " "),
      AVAL = BFLAG,
      AVALC = case_when(
        AVAL == 1 ~ "Y",
        AVAL == 0 ~ "N",
        TRUE ~ NA_character_
      ),
      AVALU = NA_character_
    )


  # Get Patient flag INDUCD
  adab_incucd <- core_aab %>%
    filter(ADATYPE == "ADA_BAB") %>%
    distinct(
      !!!get_admiral_option("subject_keys"), BASETYPE, PARAM, ADATYPE, ADAPARM, ISSPEC, ISTESTCD, ISTEST, ISBDAGNT,
      AVISITN, AVISIT, PARCAT1, TFLAG
    ) %>%
    mutate(
      PARAMCD = "INDUCDy",
      PARAM = paste("Treatment induced ADA,", PARCAT1, sep = " "),
      AVAL = case_when(
        TFLAG == 1 ~ 1,
        TRUE ~ 0
      ),
      AVALC = case_when(
        AVAL == 1 ~ "Y",
        AVAL == 0 ~ "N"
      ),
      AVALU = NA_character_
    )


  # Get Patient flag ENHANC
  adab_enhanc <- core_aab %>%
    filter(ADATYPE == "ADA_BAB") %>%
    distinct(
      !!!get_admiral_option("subject_keys"), BASETYPE, PARAM, ADATYPE, ADAPARM, ISSPEC, ISTESTCD, ISTEST, ISBDAGNT,
      AVISITN, AVISIT, PARCAT1, TFLAG
    ) %>%
    mutate(
      PARAMCD = "ENHANCy",
      PARAM = paste("Treatment enhanced ADA,", PARCAT1, sep = " "),
      AVAL = case_when(
        TFLAG == 2 ~ 1,
        TRUE ~ 0
      ),
      AVALC = case_when(
        AVAL == 1 ~ "Y",
        AVAL == 0 ~ "N"
      ),
      AVALU = NA_character_
    )


  # Get Patient flag EMERPOS
  adab_emerpos <- core_aab %>%
    filter(ADATYPE == "ADA_BAB") %>%
    distinct(
      !!!get_admiral_option("subject_keys"), BASETYPE, PARAM, ADATYPE, ADAPARM, ISSPEC, ISTESTCD, ISTEST, ISBDAGNT,
      AVISITN, AVISIT, PARCAT1, PBFLAG
    ) %>%
    mutate(
      PARAMCD = "EMERPOSy",
      PARAM = paste("Treatment Emergent - Positive,", PARCAT1, sep = " "),
      AVAL = case_when(
        PBFLAG == 1 ~ 1,
        TRUE ~ 0
      ),
      AVALC = case_when(
        AVAL == 1 ~ "Y",
        AVAL == 0 ~ "N"
      ),
      AVALU = NA_character_
    )


  # Get Patient flag TRUNAFF
  adab_trunaff <- core_aab %>%
    filter(ADATYPE == "ADA_BAB") %>%
    distinct(
      !!!get_admiral_option("subject_keys"), BASETYPE, PARAM, ADATYPE, ADAPARM, ISSPEC, ISTESTCD, ISTEST, ISBDAGNT,
      AVISITN, AVISIT, PARCAT1, TFLAG
    ) %>%
    mutate(
      PARAMCD = "TRUNAFFy",
      PARAM = paste("Treatment unaffected,", PARCAT1, sep = " "),
      AVAL = case_when(
        TFLAG == 3 ~ 1,
        TRUE ~ 0
      ),
      AVALC = case_when(
        AVAL == 1 ~ "Y",
        AVAL == 0 ~ "N"
      ),
      AVALU = NA_character_
    )


  # Get Patient flag EMERNEG
  adab_emerneg <- core_aab %>%
    filter(ADATYPE == "ADA_BAB") %>%
    distinct(
      !!!get_admiral_option("subject_keys"), BASETYPE, PARAM, ADATYPE, ADAPARM, ISSPEC, ISTESTCD, ISTEST, ISBDAGNT,
      AVISITN, AVISIT, PARCAT1, PBFLAG
    ) %>%
    mutate(
      PARAMCD = "EMERNEGy",
      PARAM = paste("Treatment Emergent - Negative,", PARCAT1, sep = " "),
      AVAL = case_when(
        PBFLAG == 0 ~ 1,
        TRUE ~ 0
      ),
      AVALC = case_when(
        AVAL == 1 ~ "Y",
        AVAL == 0 ~ "N"
      ),
      AVALU = NA_character_
    )


  # Get Patient flag NOTRREL
  adab_notrrel <- core_aab %>%
    filter(ADATYPE == "ADA_BAB") %>%
    distinct(
      !!!get_admiral_option("subject_keys"), BASETYPE, PARAM, ADATYPE, ADAPARM, ISSPEC, ISTESTCD, ISTEST, ISBDAGNT,
      AVISITN, AVISIT, PARCAT1, PBFLAG
    ) %>%
    mutate(
      PARAMCD = "NOTRRELy",
      PARAM = paste("No treatment related ADA,", PARCAT1, sep = " "),
      AVAL = case_when(
        is.na(PBFLAG) ~ 1,
        TRUE ~ 0
      ),
      AVALC = case_when(
        AVAL == 1 ~ "Y",
        AVAL == 0 ~ "N"
      ),
      AVALU = NA_character_
    )


  # Get Patient flag ADASTAT
  adab_adastat <- core_aab %>%
    filter(ADATYPE == "ADA_BAB") %>%
    distinct(
      !!!get_admiral_option("subject_keys"), BASETYPE, PARAM, ADATYPE, ADAPARM, ISSPEC, ISTESTCD, ISTEST, ISBDAGNT,
      AVISITN, AVISIT, PARCAT1, ADASTAT
    ) %>%
    mutate(
      PARAMCD = "ADASTATy",
      PARAM = paste("ADA Status of a patient,", PARCAT1, sep = " "),
      AVAL = ADASTAT,
      AVALC = case_when(
        AVAL == 1 ~ "POSITIVE",
        AVAL == 0 ~ "NEGATIVE",
        TRUE ~ "MISSING"
      ),
      AVALU = NA_character_
    )


  # Get Patient flag TIMADA
  adab_timada <- core_aab %>%
    filter(ADATYPE == "ADA_BAB") %>%
    distinct(
      !!!get_admiral_option("subject_keys"), BASETYPE, PARAM, ADATYPE, ADAPARM, ISSPEC, ISTESTCD, ISTEST, ISBDAGNT,
      AVISITN, AVISIT, PARCAT1, TIMADA
    ) %>%
    mutate(
      PARAMCD = "TIMADAy",
      PARAM = paste("Time to onset of ADA (Weeks),", PARCAT1, sep = " "),
      AVAL = TIMADA,
      AVALC = NA_character_,
      AVALU = case_when(
        !is.na(AVAL) ~ "WEEKS",
        TRUE ~ NA_character_
      )
    )


  # Get Patient flag PERSADA
  adab_persada <- core_aab %>%
    filter(ADATYPE == "ADA_BAB") %>%
    distinct(
      !!!get_admiral_option("subject_keys"), BASETYPE, PARAM, ADATYPE, ADAPARM, ISSPEC, ISTESTCD, ISTEST, ISBDAGNT,
      AVISITN, AVISIT, PARCAT1, PERSADA
    ) %>%
    mutate(
      PARAMCD = "PERSADAy",
      PARAM = paste("Persistent ADA,", PARCAT1, sep = " "),
      AVAL = PERSADA,
      AVALC = case_when(
        AVAL == 1 ~ "Y",
        AVAL == 0 ~ "N",
        TRUE ~ NA_character_
      ),
      AVALU = NA_character_
    )


  # Get Patient flag TRANADA
  adab_tranada <- core_aab %>%
    filter(ADATYPE == "ADA_BAB") %>%
    distinct(
      !!!get_admiral_option("subject_keys"), BASETYPE, PARAM, ADATYPE, ADAPARM, ISSPEC, ISTESTCD, ISTEST, ISBDAGNT,
      AVISITN, AVISIT, PARCAT1, TRANADA
    ) %>%
    mutate(
      PARAMCD = "TRANADAy",
      PARAM = paste("Transient ADA,", PARCAT1, sep = " "),
      AVAL = TRANADA,
      AVALC = case_when(
        AVAL == 1 ~ "Y",
        AVAL == 0 ~ "N",
        TRUE ~ NA_character_
      ),
      AVALU = NA_character_
    )


  # Get Patient flag NABSTAT
  adab_nabstat <- core_aab %>%
    filter(ADATYPE == "ADA_NAB") %>%
    distinct(
      !!!get_admiral_option("subject_keys"), BASETYPE, PARAM, ADATYPE, ADAPARM, ISSPEC, ISTESTCD, ISTEST, ISBDAGNT,
      AVISITN, AVISIT, PARCAT1, NABSTAT
    ) %>%
    mutate(
      PARAMCD = "NABSTATy",
      PARAM = paste("NAB Status of a patient,", PARCAT1, sep = " "),
      AVAL = NABSTAT,
      AVALC = case_when(
        AVAL == 1 ~ "POSITIVE",
        AVAL == 0 ~ "NEGATIVE",
        TRUE ~ "MISSING"
      ),
      AVALU = NA_character_
    )


  # Get Patient flag ADADUR
  adab_adadur <- core_aab %>%
    filter(ADATYPE == "ADA_BAB") %>%
    distinct(
      !!!get_admiral_option("subject_keys"), BASETYPE, PARAM, ADATYPE, ADAPARM, ISSPEC, ISTESTCD, ISTEST, ISBDAGNT,
      AVISITN, AVISIT, PARCAT1, ADADUR
    ) %>%
    mutate(
      PARAMCD = "ADADURy",
      PARAM = paste("ADA Duration (Weeks),", PARCAT1, sep = " "),
      AVAL = ADADUR,
      AVALC = NA_character_,
      AVALU = case_when(
        !is.na(AVAL) ~ "WEEKS",
        TRUE ~ NA_character_
      )
    )


  # Get Patient flag FPPDTM
  adab_fppdtm <- core_aab %>%
    filter(ADATYPE == "ADA_BAB") %>%
    distinct(
      !!!get_admiral_option("subject_keys"), BASETYPE, PARAM, ADATYPE, ADAPARM, ISSPEC, ISTESTCD, ISTEST, ISBDAGNT,
      AVISITN, AVISIT, PARCAT1, FPPDTM
    ) %>%
    mutate(
      PARAMCD = "FPPDTMy",
      PARAM = paste("First Post Dose Positive Datetime,", PARCAT1, sep = " "),
      AVAL = (as.numeric(FPPDTM)),
      AVALC = case_when(
        is.na(FPPDTM) ~ "",
        TRUE ~ toupper(as.character(format(FPPDTM, "%d%b%Y:%H:%M:%S")))
      ),
      AVALU = NA_character_
    )


  # Get Patient flag LPPDTM
  adab_lppdtm <- core_aab %>%
    filter(ADATYPE == "ADA_BAB") %>%
    distinct(
      !!!get_admiral_option("subject_keys"), BASETYPE, PARAM, ADATYPE, ADAPARM, ISSPEC, ISTESTCD, ISTEST, ISBDAGNT,
      AVISITN, AVISIT, PARCAT1, LPPDTM
    ) %>%
    mutate(
      PARAMCD = "LPPDTMy",
      PARAM = paste("Last Post Dose Positive Datetime,", PARCAT1, sep = " "),
      AVAL = as.numeric(LPPDTM),
      AVALC = case_when(
        is.na(LPPDTM) ~ "",
        TRUE ~ toupper(as.character(format(LPPDTM, "%d%b%Y:%H:%M:%S")))
      ),
      AVALU = NA_character_
    )



  # Set all the standard PARAM components together -----------------------------------

  adab_paramcds <- bind_rows(
    adab_titer, adab_nabvis, adab_result, adab_nabres, adab_bflag,
    adab_incucd, adab_enhanc, adab_emerpos, adab_trunaff, adab_emerneg, adab_notrrel,
    adab_adastat, adab_timada, adab_persada, adab_tranada, adab_nabstat
  )

  # In this sample, also have BY VISIT parameters,
  adab_visits <- bind_rows(adab_tflagv, adab_pbflagv, adab_adastatv)


  # Create the final parameterized dataset --------------------------------------------

  # Standard parameters:
  #adab_study <- adab_paramcds
  # To include BY VISIT parameters
  adab_study <- bind_rows(adab_paramcds,adab_visits)


  # Drop Temp vars and ADA Flag vars that are now parameterized
  adab_study <- adab_study %>%
    select(
      -TIMADA, -ADADUR, -TRANADA, -PERSADA, -TRANADAE, -PERSADAE, -INDUCED, -ENHANCED, -RESULTC, -RESULTN,
      -ADASTAT, -BFLAG, -TFLAG, -PBFLAG, -FPPDTM, -LPPDTM, -TFLAGV, -PBFLAGV, -ADASTATV,
      -nabstat_opt1, -nabstat_opt2, -NABSTAT, -MAXCHG, -VALIDBASE, -VALIDPOST,
       -tdur, -ADASTAT_MAIN, -NABPOSTMISS, -TRTSDT
    )


  # Merge in ADSL Static and Computed Values --------------------------------

  adab_adsl <- adab_study %>%
    derive_vars_merged(
      dataset_add = adsl,
      by_vars = exprs(!!!get_admiral_option("subject_keys"))
    )


  # Compute COHORT From ADSL, other source could be ADSUB
  adab_cohort <- adab_adsl %>%
    derive_vars_merged(
      dataset_add = adsl,
      new_vars = exprs(COHORT = ARMCD),
      by_vars = exprs(!!!get_admiral_option("subject_keys"))
    )


  # Compute ADAFL
  # Other method could be ADSL.SAFFL, ADAB.ADPBLPFL by ISTESTCD, etc.
  adab_adafl <- adab_cohort %>%
    derive_vars_merged(
      dataset_add = adsl,
      new_vars = exprs(ADAFL = SAFFL),
      by_vars = exprs(!!!get_admiral_option("subject_keys"))
    )


  view_keys1  <- adab_adafl %>%
    distinct (ISTESTCD, ISTEST, ISBDAGNT, PARCAT1, ADATYPE, ADAPARM, PARAMCD, PARAM)  %>%
    arrange (ISTESTCD,  ISTEST, ISBDAGNT, PARCAT1, ADAPARM, ADATYPE, PARAMCD, PARAM)



  # Study Specific Specs Post-Processing ------------------------------------

  # Create a Tibble to map above processed PARAMCD to the study specs\
  #  PARAM and PARAMCD


  adab_param_data <- tribble(
    ~PARAMCD,   ~ADATYPE,  ~ADAPARM, ~PARAMCD_NEW, ~PARAM_SUFFIX,
    "RESULTy",   "ADA_BAB", "XANOMELINE", "RESULT",  NA_character_,
    "RESULTy",   "ADA_NAB", "XANOMELINE", "RESULT2", NA_character_,
    "BFLAGy",   "ADA_BAB", "XANOMELINE", "BFLAG",   NA_character_,
    "INDUCDy",  "ADA_BAB", "XANOMELINE", "INDUCD",  NA_character_,
    "ENHANCy",  "ADA_BAB", "XANOMELINE", "ENHANC",  NA_character_,
    "EMERPOSy", "ADA_BAB", "XANOMELINE", "EMERPOS", NA_character_,
    "TRUNAFFy", "ADA_BAB", "XANOMELINE", "TRUNAFF", NA_character_,
    "EMERNEGy", "ADA_BAB", "XANOMELINE", "EMERNEG", NA_character_,
    "NOTRRELy", "ADA_BAB", "XANOMELINE", "NOTRREL", NA_character_,
    "ADASTATy", "ADA_BAB", "XANOMELINE", "ADASTAT", NA_character_,
    "TIMADAy",  "ADA_BAB", "XANOMELINE", "TIMADA",  NA_character_,
    "PERSADAy", "ADA_BAB", "XANOMELINE", "PERSADA", NA_character_,
    "TRANADAy", "ADA_BAB", "XANOMELINE", "TRANADA", NA_character_,
    "NABSTATy", "ADA_NAB", "XANOMELINE", "NABSTAT", NA_character_,
    "ADASTATV", "ADA_BAB", "XANOMELINE", "ADASTTV",  NA_character_,
    "TFLAGV",   "ADA_BAB", "XANOMELINE", "TFLAGV",   NA_character_,
    "PBFLAGV",  "ADA_BAB", "XANOMELINE", "PBFLAGV",  NA_character_,
  )


  # Merge the Paramter dataset into the main data
  adab_params <- adab_adafl %>%
    derive_vars_merged(
      dataset_add = adab_param_data,
      by_vars = exprs(PARAMCD, ADAPARM,ADATYPE )
    ) %>%
    # for the original data, use original Test Code full text vars
    mutate(
      PARAMCD = case_when(
        PARAMCD == ADAPARM ~ PARAMCD,
        TRUE ~ PARAMCD_NEW
      ),
      # Assign PARAM_NEW based on PARAM and PARAMCD_NEW
      PARAM = case_when(
        !is.na(PARAM_SUFFIX) ~ paste(PARAM, PARAM_SUFFIX, sep = " "),
        TRUE ~ PARAM
      )
    )


  view_keys2  <- adab_params %>%
    distinct (ISTESTCD, ISTEST, ISBDAGNT, PARCAT1, ADAPARM, PARAMCD, PARAM)  %>%
    arrange (ISTESTCD,  ISTEST, ISBDAGNT, PARCAT1, ADAPARM, PARAMCD, PARAM)


  # Sort by the standard Key then Compute ASEQ


  adab_prefinal <- adab_params %>%
    # Calculate ASEQ
    derive_var_obs_number(
      new_var = ASEQ,
      by_vars = exprs(!!!get_admiral_option("subject_keys")),
      order = exprs(PARCAT1, PARAMCD, BASETYPE, NFRLT, AFRLT, ISSEQ),
      check_type = "error"
    )



  review_adab_prefinal <- adab_prefinal  %>%
    select(USUBJID, SUBJID, ISDTC, ADTM, FANLDTM, NFRLT, AFRLT,FRLTU, ASEQ, PARAMCD, PARAM, ISSTRESN, ISSTRESC, AVAL, AVALC, BASE, CHG)   %>%
    filter ( 1==1) %>%
    filter (PARAMCD == "RESULT2")


  adab <- adab_prefinal  %>%
    select(
      STUDYID, USUBJID, SUBJID, SITEID, ASEQ,
      REGION1,  COUNTRY, ETHNIC,
      AGE, AGEU, SEX, RACE,
      SAFFL, TRT01P, TRT01A,
      TRTSDTM, TRTSDT, TRTEDTM, TRTEDT,
      ISSEQ, ISTESTCD, ISTEST, ISCAT, ISBDAGNT,
      ISSTRESC, ISSTRESN, ISSTRESU,
      ISSTAT, ISREASND, ISSPEC,
      DTL, MRT,
      VISITNUM, VISIT, VISITDY,
      EPOCH, ISDTC, ISDY, ISTPT, ISTPTNUM,
      PARAM, PARAMCD, PARCAT1,
      AVAL,AVALC, AVALU,
      BASETYPE, BASE, CHG, DTYPE,
      ADTM, ADT, ADY, ATMF,
      AVISIT, AVISITN, ATPT,
      APHASE, APHASEN, APERIOD, APERIODC,
      FANLDTM, FANLDT, FANLTM, FANLTMF,
      NFRLT, AFRLT, FRLTU,
      ABLFL, ADABLPFL, ADPBLPFL, ADAFL
    )

  # Save output ----

  # Save as PARQUET file
  arrow::write_parquet(adab, file.path("data/pharmaverse/adab.parquet"))







