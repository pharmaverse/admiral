# Name: ADLB
#
# Label: Lab Analysis Dataset
#
# Input: adsl, lb
library(admiral)
library(pharmaversesdtm) # Contains example datasets from the CDISC pilot project
library(dplyr)
library(lubridate)
library(stringr)

# Load source datasets ----

# Use e.g. haven::read_sas to read in .sas7bdat, or other suitable functions
# as needed and assign to the variables below.
# For illustration purposes read in admiral test data

lb <- pharmaversesdtm::lb
adsl <- admiral::admiral_adsl

# When SAS datasets are imported into R using haven::read_sas(), missing
# character values from SAS appear as "" characters in R, instead of appearing
# as NA values. Further details can be obtained via the following link:
# https://pharmaverse.github.io/admiral/articles/admiral.html#handling-of-missing-values # nolint

lb <- convert_blanks_to_na(lb)

# Look-up tables ----

# Assign PARAMCD, PARAM, and PARAMN
param_lookup <- tibble::tribble(
  ~LBTESTCD, ~PARAMCD,                                             ~PARAM, ~PARAMN,
  "ALB",        "ALB",                                    "Albumin (g/L)",       1,
  "ALP",      "ALKPH",                       "Alkaline Phosphatase (U/L)",       2,
  "ALT",        "ALT",                   "Alanine Aminotransferase (U/L)",       3,
  "ANISO",    "ANISO",                                       "Anisocytes",       4,
  "AST",        "AST",                 "Aspartate Aminotransferase (U/L)",       5,
  "BASO",      "BASO",                           "Basophils Abs (10^9/L)",       6,
  "BASOLE",  "BASOLE",                  "Basophils/Leukocytes (FRACTION)",       7,
  "BILI",      "BILI",                               "Bilirubin (umol/L)",       8,
  "BUN",        "BUN",                     "Blood Urea Nitrogen (mmol/L)",       9,
  "CA",          "CA",                                 "Calcium (mmol/L)",      10,
  "CHOL",    "CHOLES",                             "Cholesterol (mmol/L)",      11,
  "CK",          "CK",                          "Creatinine Kinase (U/L)",      12,
  "CL",          "CL",                                "Chloride (mmol/L)",      13,
  "COLOR",    "COLOR",                                            "Color",      14,
  "CREAT",    "CREAT",                              "Creatinine (umol/L)",      15,
  "EOS",        "EOS",                             "Eosinophils (10^9/L)",      16,
  "EOSLE",    "EOSLE",                "Eosinophils/Leukocytes (FRACTION)",      17,
  "GGT",        "GGT",                 "Gamma Glutamyl Transferase (U/L)",      18,
  "GLUC",      "GLUC",                                 "Glucose (mmol/L)",      19,
  "HBA1C",    "HBA1C",                               "Hemoglobin A1C (1)",      20,
  "HCT",        "HCT",                                   "Hematocrit (1)",      21,
  "HGB",        "HGB",                              "Hemoglobin (mmol/L)",      22,
  "K",        "POTAS",                               "Potassium (mmol/L)",      23,
  "KETONES",  "KETON",                                          "Ketones",      24,
  "LYM",      "LYMPH",                         "Lymphocytes Abs (10^9/L)",      25,
  "LYMLE",  "LYMPHLE",                "Lymphocytes/Leukocytes (FRACTION)",      26,
  "MACROCY", "MACROC",                                       "Macrocytes",      27,
  "MCH",        "MCH",      "Ery. Mean Corpuscular Hemoglobin (fmol(Fe))",      28,
  "MCHC",      "MCHC", "Ery. Mean Corpuscular HGB Concentration (mmol/L)",      29,
  "MCV",        "MCV",               "Ery. Mean Corpuscular Volume (f/L)",      30,
  "MICROCY", "MICROC",                                       "Microcytes",      31,
  "MONO",      "MONO",                               "Monocytes (10^9/L)",      32,
  "MONOLE",  "MONOLE",                  "Monocytes/Leukocytes (FRACTION)",      33,
  "PH",          "PH",                                               "pH",      34,
  "PHOS",      "PHOS",                               "Phosphate (mmol/L)",      35,
  "PLAT",      "PLAT",                                "Platelet (10^9/L)",      36,
  "POIKILO", "POIKIL",                                     "Poikilocytes",      37,
  "POLYCHR", "POLYCH",                                    "Polychromasia",      38,
  "PROT",      "PROT",                                    "Protein (g/L)",      39,
  "RBC",        "RBC",                              "Erythrocytes (TI/L)",      40,
  "SODIUM",  "SODIUM",                                  "Sodium (mmol/L)",      41,
  "SPGRAV",  "SPGRAV",                                 "Specific Gravity",      42,
  "TSH",        "TSH",                               "Thyrotropin (mU/L)",      43,
  "URATE",    "URATE",                                   "Urate (umol/L)",      44,
  "UROBIL",  "UROBIL",                                     "Urobilinogen",      45,
  "VITB12",  "VITB12",                             "Vitamin B12 (pmol/L)",      46,
  "WBC",        "WBC",                              "Leukocytes (10^9/L)",      47
)


# Derivations ----

# Get list of ADSL vars required for derivations
adsl_vars <- exprs(TRTSDT, TRTEDT, TRT01A, TRT01P)

adlb <- lb %>%
  # Join ADSL with LB (need TRTSDT for ADY derivation)
  derive_vars_merged(
    dataset_add = adsl,
    new_vars = adsl_vars,
    by_vars = exprs(STUDYID, USUBJID)
  ) %>%
  ## Calculate ADT, ADY ----
  derive_vars_dt(
    new_vars_prefix = "A",
    dtc = LBDTC
  ) %>%
  derive_vars_dy(reference_date = TRTSDT, source_vars = exprs(ADT))

adlb <- adlb %>%
  ## Add PARAMCD PARAM and PARAMN - from LOOK-UP table ----
  # Replace with PARAMCD lookup function
  derive_vars_merged_lookup(
    dataset_add = param_lookup,
    new_vars = exprs(PARAMCD, PARAM, PARAMN),
    by_vars = exprs(LBTESTCD),
    check_type = "none",
    print_not_mapped = FALSE
  ) %>%
  ## Calculate PARCAT1 AVAL AVALC ANRLO ANRHI ----
  mutate(
    PARCAT1 = LBCAT,
    AVAL = LBSTRESN,
    # Only populate AVALC if character value is non-redundant with AVAL
    AVALC = ifelse(
      is.na(LBSTRESN) | as.character(LBSTRESN) != LBSTRESC,
      LBSTRESC,
      NA
    ),
    ANRLO = LBSTNRLO,
    ANRHI = LBSTNRHI
  )

# Derive Absolute values from fractional Differentials using WBC
# Only derive where absolute values do not already exist
# Need to populate ANRLO and ANRHI for newly created records
adlb <- adlb %>%
  # Derive absolute Basophils
  derive_param_wbc_abs(
    by_vars = exprs(STUDYID, USUBJID, !!!adsl_vars, DOMAIN, VISIT, VISITNUM, ADT, ADY),
    set_values_to = exprs(
      PARAMCD = "BASO",
      PARAM = "Basophils Abs (10^9/L)",
      PARAMN = 6,
      DTYPE = "CALCULATION",
      PARCAT1 = "HEMATOLOGY"
    ),
    get_unit_expr = extract_unit(PARAM),
    diff_code = "BASOLE"
  ) %>%
  # Derive absolute Lymphocytes
  derive_param_wbc_abs(
    by_vars = exprs(STUDYID, USUBJID, !!!adsl_vars, DOMAIN, VISIT, VISITNUM, ADT, ADY),
    set_values_to = exprs(
      PARAMCD = "LYMPH",
      PARAM = "Lymphocytes Abs (10^9/L)",
      PARAMN = 25,
      DTYPE = "CALCULATION",
      PARCAT1 = "HEMATOLOGY"
    ),
    get_unit_expr = extract_unit(PARAM),
    diff_code = "LYMPHLE"
  )

## Get Visit Info ----
# See also the "Visit and Period Variables" vignette
# (https://pharmaverse.github.io/admiral/articles/visits_periods.html#visits)
adlb <- adlb %>%
  # Derive Timing
  mutate(
    AVISIT = case_when(
      str_detect(VISIT, "SCREEN") ~ "Baseline",
      !is.na(VISIT) ~ str_to_title(VISIT),
      TRUE ~ NA_character_
    ),
    AVISITN = case_when(
      AVISIT == "Baseline" ~ 0,
      !is.na(VISITNUM) ~ VISITNUM
    )
  )

adlb <- adlb %>%
  ## Calculate ONTRTFL ----
  derive_var_ontrtfl(
    start_date = ADT,
    ref_start_date = TRTSDT,
    ref_end_date = TRTEDT,
    filter_pre_timepoint = AVISIT == "Baseline"
  )

## Calculate ANRIND : requires the reference ranges ANRLO, ANRHI ----
adlb <- adlb %>%
  derive_var_anrind()

## Derive baseline flags ----
adlb <- adlb %>%
  # Calculate BASETYPE
  mutate(
    BASETYPE = "LAST"
  ) %>%
  # Calculate ABLFL
  restrict_derivation(
    derivation = derive_var_extreme_flag,
    args = params(
      by_vars = exprs(STUDYID, USUBJID, BASETYPE, PARAMCD),
      order = exprs(ADT, VISITNUM, LBSEQ),
      new_var = ABLFL,
      mode = "last"
    ),
    filter = (!is.na(AVAL) & ADT <= TRTSDT & !is.na(BASETYPE))
  )

## Derive baseline information ----
adlb <- adlb %>%
  # Calculate BASE
  derive_var_base(
    by_vars = exprs(STUDYID, USUBJID, PARAMCD, BASETYPE),
    source_var = AVAL,
    new_var = BASE
  ) %>%
  # Calculate BASEC
  derive_var_base(
    by_vars = exprs(STUDYID, USUBJID, PARAMCD, BASETYPE),
    source_var = AVALC,
    new_var = BASEC
  ) %>%
  # Calculate BNRIND
  derive_var_base(
    by_vars = exprs(STUDYID, USUBJID, PARAMCD, BASETYPE),
    source_var = ANRIND,
    new_var = BNRIND
  ) %>%
  # Calculate CHG for post-baseline records
  # The decision on how to populate pre-baseline and baseline values of CHG is left to producer choice
  restrict_derivation(
    derivation = derive_var_chg,
    filter = AVISITN > 0
  ) %>%
  # Calculate PCHG for post-baseline records
  # The decision on how to populate pre-baseline and baseline values of PCHG is left to producer choice
  restrict_derivation(
    derivation = derive_var_pchg,
    filter = AVISITN > 0
  )

## Calculate lab grading ----

# Assign ATOXDSCL and ATOXDSCH to hold lab grading terms
# ATOXDSCL and ATOXDSCH hold terms defined by NCI-CTCAEv4.
# See (https://pharmaverse.github.io/admiral/articles/lab_grading.html#implement_ctcv4)
grade_lookup <- tibble::tribble(
  ~PARAMCD,                  ~ATOXDSCL,                              ~ATOXDSCH,
  "ALB",             "Hypoalbuminemia",                          NA_character_,
  "ALKPH",               NA_character_,       "Alkaline phosphatase increased",
  "ALT",                 NA_character_,   "Alanine aminotransferase increased",
  "AST",                 NA_character_, "Aspartate aminotransferase increased",
  "BILI",                NA_character_,            "Blood bilirubin increased",
  "CA",                 "Hypocalcemia",                        "Hypercalcemia",
  "CHOLES",              NA_character_,                     "Cholesterol high",
  "CK",                  NA_character_,                        "CPK increased",
  "CREAT",               NA_character_,                 "Creatinine increased",
  "GGT",                 NA_character_,                        "GGT increased",
  "GLUC",               "Hypoglycemia",                        "Hyperglycemia",
  "HGB",                      "Anemia",                 "Hemoglobin increased",
  "POTAS",               "Hypokalemia",                         "Hyperkalemia",
  "LYMPH", "CD4 lymphocytes decreased",                          NA_character_,
  "PHOS",           "Hypophosphatemia",                          NA_character_,
  "PLAT",   "Platelet count decreased",                          NA_character_,
  "SODIUM",             "Hyponatremia",                        "Hypernatremia",
  "WBC",  "White blood cell decreased",                         "Leukocytosis",
)

adlb_test <- adlb %>%
  filter(PARAMCD == "ALB") %>%
  mutate(PARAM = "Albumin (g/dL)",
         AVAL = AVAL/10,
         ANRLO = ANRLO/10
         )

adlb_all <- adlb %>%
  bind_rows(adlb_test)

# Assign grade criteria
# metadata atoxgr_criteria_ctcv4 used to implement NCI-CTCAEv4
# user could change to atoxgr_criteria_ctcv5 to implement NCI-CTCAEv5
# Note: Hyperglycemia and Hypophosphatemia not defined in NCI-CTCAEv5 so
# user would need to amend look-up table grade_lookup
# See (https://pharmaverse.github.io/admiral/articles/lab_grading.html#implement_ctcv5)
grade_crit <- atoxgr_criteria_ctcv5

case_when(  is.na(AVAL) ~ NA_character_,  signif(AVAL, signif_dig) < 2 ~ "3",  signif(AVAL, signif_dig) < 3 ~ "2",  is.na(ANRLO) ~ NA_character_,  signif(AVAL, signif_dig) < signif(ANRLO, signif_dig) ~ "1",  signif(AVAL, signif_dig) >= signif(ANRLO, signif_dig) ~ "0"  )


grade_crit_alb <- grade_crit %>%
  filter(TERM == "Hypoalbuminemia") %>%
  mutate(SI_UNIT_CHECK = "g/dL",
         GRADE_CRITERIA_CODE = 'case_when(  is.na(AVAL) ~ NA_character_,  signif(AVAL, signif_dig) < 2 ~ "3",  signif(AVAL, signif_dig) < 3 ~ "2",  is.na(ANRLO) ~ NA_character_,  signif(AVAL, signif_dig) < signif(ANRLO, signif_dig) ~ "1",  signif(AVAL, signif_dig) >= signif(ANRLO, signif_dig) ~ "0"  )',
  )

grade_crit_all <- grade_crit %>%
  bind_rows(grade_crit_alb)


adlb_alb <- adlb_chk2 %>%
  filter(PARAMCD == "ALB" & ATOXGRL > 1)

terms_in_vad <- input1_dbili_daids %>%
  filter(!is.na(ATOXDSCH)) %>%
  distinct(ATOXDSCH) %>%
  mutate(
    TERM = ATOXDSCH,
    TERM_UPPER = toupper(TERM)
  )
atoxgr_dir <- meta %>%
  filter(!is.na(GRADE_CRITERIA_CODE) & toupper(DIRECTION) == toupper("H")) %>%
  select(TERM, DIRECTION, SI_UNIT_CHECK, FILTER, GRADE_CRITERIA_CODE, VAR_CHECK) %>%
  mutate(
    TERM_UPPER = toupper(TERM),
    SI_UNIT_UPPER = toupper(SI_UNIT_CHECK)
  )


# only keep terms that exist in both ADLB data and criteria metadata
list_of_terms <- terms_in_vad %>%
  semi_join(atoxgr_dir, by = "TERM_UPPER") %>%
  arrange(TERM)

out_data <- input_dbili_daids %>%
  filter(ATOXDSCH %notin% (list_of_terms$TERM) | is.na(ATOXDSCH)) %>%
  mutate(ATOXGRH := NA_character_)
to_be_graded <- input_dbili_daids %>%
  filter(ATOXDSCH %in% (list_of_terms$TERM))


for (i in seq_along(list_of_terms$TERM)) {
  # filter metadata on a term
  meta_this_term <- atoxgr_dir %>%
    filter(TERM_UPPER == list_of_terms$TERM_UPPER[i])

  grade_this_term <- to_be_graded %>%
    filter(ATOXDSCH == list_of_terms$TERM[i])

  # get unique list of FILTERS (possibly more than one unit)
  meta_this_filter_uniq <- meta_this_term %>%
    select(FILTER) %>%
    distinct()
  signif_dig <- get_admiral_option("signif_digits")
  # Within each TERM check if there are FILTERs to be applied
  # if FILTER not missing then loop through each FILTER for the TERM already specified
  for (j in seq_along(meta_this_filter_uniq$FILTER)) {
    # subset using FILTER if its not empty
    if (!is.na(meta_this_filter_uniq$FILTER[j])) {
      meta_this_filter <- meta_this_term %>%
        filter(FILTER == meta_this_term$FILTER[j])

      # filter lab data using FILTER from metadata
      grade_this_filter <- grade_this_term %>%
        filter(eval(parse(text = meta_this_filter_uniq$FILTER[j])))

    } else {
      meta_this_filter <- meta_this_term

      grade_this_filter <- grade_this_term
    }


    for (x in seq_along(meta_this_filter$SI_UNIT_CHECK)) {

      if (!is.na(meta_this_filter$SI_UNIT_CHECK[x])) {

        meta_this_filter_unit <- meta_this_filter %>%
          filter(SI_UNIT_CHECK == meta_this_filter$SI_UNIT_CHECK[x])

        grade_this_filter_unit <- grade_this_filter %>%
          filter(meta_this_filter_unit$SI_UNIT_UPPER == toupper(AVALU) |
                   is.na(toupper(AVALU)))
      } else {
        meta_this_filter_unit <- meta_this_filter

        grade_this_filter_unit <- grade_this_filter
      }


      # Put list of variables required for criteria in a vector
      list_of_vars <- gsub("\\s+", "", unlist(strsplit(meta_this_filter_unit$VAR_CHECK, ",")))

      # check variables required in criteria exist on data
      assert_data_frame(grade_this_filter_unit, required_vars = exprs(!!!syms(list_of_vars)))

      if ("BNRIND" %in% list_of_vars) {
        # check input parameter is character value
        assert_character_vector(abnormal_indicator, optional = FALSE)
      }

      # apply criteria when SI unit matches
      grade_this_filter_unit <- grade_this_filter_unit %>%
        mutate(
          temp_flag = meta_this_filter_unit$SI_UNIT_UPPER == toupper(AVALU) |
            is.na(meta_this_filter_unit$SI_UNIT_UPPER),
          ATOXGRH = if_else(
            temp_flag, eval(parse(text = meta_this_filter_unit$GRADE_CRITERIA_CODE)), NA_character_
          )
        ) %>%
        select(-temp_flag)

      # add data just graded to data already processed
      out_data <- bind_rows(out_data, grade_this_filter_unit)

      if (!is.na(meta_this_filter$SI_UNIT_CHECK[x])) {
        grade_this_filter <- grade_this_filter %>%
          filter(meta_this_filter_unit$SI_UNIT_UPPER != toupper(AVALU) &
          !is.na(meta_this_filter_unit$SI_UNIT_UPPER))

        if (x == length(meta_this_filter$SI_UNIT_CHECK)) {
          out_data <- bind_rows(out_data, grade_this_filter)
        }
      }

      if (j == 2 & x == 1L){
        filter_start21 <- grade_this_filter_unit
        filter_rem12 <- grade_this_filter
        filter_meta2 <- meta_this_filter
        out_data2 <-  out_data
      }

    }
    # remove lab data just graded from data still to be graded for the specified TERM
    grade_this_term <- grade_this_term %>%
      filter(!(eval(parse(text = meta_this_filter_unit$FILTER))))

  }


  # remove lab data with TERM just graded from data still to be graded
  to_be_graded <- to_be_graded %>%
    filter(ATOXDSCH != list_of_terms$TERM[i])
}


NCICTC5 <- grade_crit %>%
  filter(!is.na(SI_UNIT_CHECK))



























for (i in seq_along(list_of_terms$TERM)) {
 if (i == length(list_of_terms$TERM)) {
   x <- "Y"
 }
}

test <- length(list_of_terms$TERM)

# Add ATOXDSCL and ATOXDSCH
adlb_chk3 <- adlb_all %>%
  derive_vars_merged(
    dataset_add = grade_lookup,
    by_vars = exprs(PARAMCD)
  ) %>%
  # Derive toxicity grade for low values ATOXGRL

  derive_var_atoxgr_dir(
    meta_criteria = grade_crit_all,
    new_var = ATOXGRL,
    tox_description_var = ATOXDSCL,
    criteria_direction = "L",
    get_unit_expr = extract_unit(PARAM)
  ) %>%
  # Derive toxicity grade for low values ATOXGRH
  # default metadata atoxgr_criteria_ctcv4 used
  derive_var_atoxgr_dir(
    meta_criteria = grade_crit,
    new_var = ATOXGRH,
    tox_description_var = ATOXDSCH,
    criteria_direction = "H",
    abnormal_indicator = "HIGH",
    get_unit_expr = extract_unit(PARAM)
  ) %>%
  # (Optional) derive overall grade ATOXGR (combining ATOXGRL and ATOXGRH)
  derive_var_atoxgr() %>%
  # Derive baseline toxicity grade for low values BTOXGRL
  derive_var_base(
    by_vars = exprs(STUDYID, USUBJID, PARAMCD, BASETYPE),
    source_var = ATOXGRL,
    new_var = BTOXGRL
  ) %>%
  # Derive baseline toxicity grade for high values BTOXGRH
  derive_var_base(
    by_vars = exprs(STUDYID, USUBJID, PARAMCD, BASETYPE),
    source_var = ATOXGRH,
    new_var = BTOXGRH
  ) %>%
  # Derive baseline toxicity grade for for overall grade BTOXGR
  derive_var_base(
    by_vars = exprs(STUDYID, USUBJID, PARAMCD, BASETYPE),
    source_var = ATOXGR,
    new_var = BTOXGR
  )


## Calculate R2BASE, R2ANRLO and R2ANRHI ----
adlb <- adlb %>%
  derive_var_analysis_ratio(
    numer_var = AVAL,
    denom_var = BASE
  ) %>%
  derive_var_analysis_ratio(
    numer_var = AVAL,
    denom_var = ANRLO
  ) %>%
  derive_var_analysis_ratio(
    numer_var = AVAL,
    denom_var = ANRHI
  )

## SHIFT derivation ----
adlb <- adlb %>%
  # Derive shift from baseline for analysis indicator
  derive_var_shift(
    new_var = SHIFT1,
    from_var = BNRIND,
    to_var = ANRIND
  ) %>%
  # Derive shift from baseline for overall grade
  restrict_derivation(
    derivation = derive_var_shift,
    args = params(
      new_var = SHIFT2,
      from_var = BTOXGR,
      to_var = ATOXGR
    ),
    filter = !is.na(ATOXDSCL) | !is.na(ATOXDSCH)
  )

## Flag variables (ANL01FL, LVOTFL) ----
# ANL01FL: Flag last result within an AVISIT for post-baseline records
# LVOTFL: Flag last valid on-treatment record
adlb <- adlb %>%
  restrict_derivation(
    derivation = derive_var_extreme_flag,
    args = params(
      by_vars = exprs(USUBJID, PARAMCD, AVISIT),
      order = exprs(ADT, AVAL),
      new_var = ANL01FL,
      mode = "last"
    ),
    filter = !is.na(AVISITN) & ONTRTFL == "Y"
  ) %>%
  restrict_derivation(
    derivation = derive_var_extreme_flag,
    args = params(
      by_vars = exprs(USUBJID, PARAMCD),
      order = exprs(ADT, AVAL),
      new_var = LVOTFL,
      mode = "last"
    ),
    filter = ONTRTFL == "Y"
  )

## Get treatment information ----
# See also the "Visit and Period Variables" vignette
# (https://pharmaverse.github.io/admiral/articles/visits_periods.html#treatment_bds)
adlb <- adlb %>%
  # Assign TRTA, TRTP
  mutate(
    TRTP = TRT01P,
    TRTA = TRT01A
  )

## Get extreme values ----
adlb <- adlb %>%
  # get MINIMUM value
  derive_extreme_records(
    dataset_add = adlb,
    by_vars = exprs(STUDYID, USUBJID, PARAMCD, BASETYPE),
    order = exprs(AVAL, ADT, AVISITN),
    mode = "first",
    filter_add = (!is.na(AVAL) & ONTRTFL == "Y"),
    set_values_to = exprs(
      AVISITN = 9997,
      AVISIT = "POST-BASELINE MINIMUM",
      DTYPE = "MINIMUM"
    )
  ) %>%
  # get MAXIMUM value
  derive_extreme_records(
    dataset_add = adlb,
    by_vars = exprs(STUDYID, USUBJID, PARAMCD, BASETYPE),
    order = exprs(desc(AVAL), ADT, AVISITN),
    mode = "first",
    filter_add = (!is.na(AVAL) & ONTRTFL == "Y"),
    set_values_to = exprs(
      AVISITN = 9998,
      AVISIT = "POST-BASELINE MAXIMUM",
      DTYPE = "MAXIMUM"
    )
  ) %>%
  # get LOV value
  derive_extreme_records(
    dataset_add = adlb,
    by_vars = exprs(STUDYID, USUBJID, PARAMCD, BASETYPE),
    order = exprs(ADT, AVISITN),
    mode = "last",
    filter_add = (ONTRTFL == "Y"),
    set_values_to = exprs(
      AVISITN = 9999,
      AVISIT = "POST-BASELINE LAST",
      DTYPE = "LOV"
    )
  )

## Get ASEQ ----
adlb <- adlb %>%
  # Calculate ASEQ
  derive_var_obs_number(
    new_var = ASEQ,
    by_vars = exprs(STUDYID, USUBJID),
    order = exprs(PARAMCD, ADT, AVISITN, VISITNUM),
    check_type = "error"
  )

# Add all ADSL variables
adlb <- adlb %>%
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
save(adlb, file = file.path(dir, "adlb.rda"), compress = "bzip2")
