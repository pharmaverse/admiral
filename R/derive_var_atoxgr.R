#' Derive Lab Toxicity Grade 0 - 4
#'
#' @description
#' Derives a character lab grade based on severity/toxicity criteria.
#'
#' @param dataset Input data set
#'
#'   The columns specified by `tox_description_var` parameter is expected.
#'
#' @param new_var Name of the character grade variable to create, for example, `ATOXGRH`
#' or `ATOXGRL`.
#'
#' @param tox_description_var Variable containing the description of the grading
#' criteria. For example: "Anemia" or "INR Increased".
#'
#' @param meta_criteria Metadata data set holding the criteria (normally a case statement)
#'
#'   Default: `atoxgr_criteria_ctcv4`
#'
#'   {admiral} metadata data set `atoxgr_criteria_ctcv4` implements
#'   [Common Terminology Criteria for Adverse Events (CTCAE)
#'    v4.0](https://ctep.cancer.gov/protocoldevelopment/electronic_applications/ctc.htm)
#'
#'   The metadata should have the following variables:
#'
#' - `TERM`: variable to hold the term describing the criteria applied to a particular lab test,
#'   eg. "Anemia" or "INR Increased". Note: the variable is case insensitive.
#' - `DIRECTION`: variable to hold the direction of the abnormality of a particular lab test
#'   value. "L" is for LOW values, "H" is for HIGH values. Note: the variable is case insensitive.
#' - `SI_UNIT_CHECK`: variable to hold unit of particular lab test. Used to check against input data
#'   if criteria is based on absolute values.
#' - `VAR_CHECK`: variable to hold comma separated list of variables used in criteria. Used to check
#'   against input data that variables exist.
#' - `GRADE_CRITERIA_CODE`: variable to hold code that creates grade based on defined criteria.
#'
#' @param criteria_direction Direction (L= Low, H = High) of toxicity grade.
#'
#' Permitted Values: "L", "H"
#'
#' @param get_unit_expr An expression providing the unit of the parameter
#'
#'   The result is used to check the units of the input parameters. Compared with
#'   `SI_UNIT_CHECK` in metadata (see `meta_criteria` parameter).
#'
#'   Permitted Values: A variable containing unit from the input dataset, or a function call,
#'   for example, `get_unit_expr = extract_unit(PARAM)`.
#'
#'
#' @details
#' `new_var` is derived with values NA, "0", "1", "2", "3", "4", where "4" is the most
#' severe grade
#' - "4" is where the lab value satisfies the criteria for grade 4.
#' - "3" is where the lab value satisfies the criteria for grade 3.
#' - "2" is where the lab value satisfies the criteria for grade 2.
#' - "1" is where the lab value satisfies the criteria for grade 1.
#' - "0" is where a grade can be derived and is not grade "1", "2", "3" or "4".
#' - NA is where a grade cannot be derived.
#'
#' @author Gordon Miller
#'
#' @return The input dataset with the character variable added
#'
#' @keywords der_bds_findings
#'
#' @family der_bds_findings
#'
#' @export
#'
#' @examples
#' library(dplyr, warn.conflicts = FALSE)
#'
#' data <- tibble::tribble(
#'   ~ATOXDSCL,                     ~AVAL,  ~ANRLO,   ~ANRHI, ~PARAM,
#'   "Hypoglycemia",                119,    4,        7,      "Glucose (mmol/L)",
#'   "Hypoglycemia",                120,    4,        7,      "Glucose (mmol/L)",
#'   "Anemia",                      129,    120,      180,    "Hemoglobin (g/L)",
#'   "White blood cell decreased",  10,     5,        20,     "White blood cell (10^9/L)",
#'   "White blood cell decreased",  15,     5,        20,     "White blood cell (10^9/L)",
#'   "Anemia",                      140,    120,      180,    "Hemoglobin (g/L)"
#' )
#'
#' derive_var_atoxgr_dir(data,
#'   new_var = ATOXGRL,
#'   tox_description_var = ATOXDSCL,
#'   meta_criteria = atoxgr_criteria_ctcv4,
#'   criteria_direction = "L",
#'   get_unit_expr = extract_unit(PARAM)
#' )
#'
#' data <- tibble::tribble(
#'   ~ATOXDSCH,                     ~AVAL,  ~ANRLO,   ~ANRHI, ~PARAM,
#'   "Hyperglycemia",               119,    4,        7,      "Glucose (mmol/L)",
#'   "Hyperglycemia",               120,    4,        7,      "Glucose  (mmol/L)",
#'   "GGT increased",               129,    0,        30,     "Gamma Glutamyl Transferase (U/L)",
#'   "Lymphocyte count increased",  4,      1,        4,      "Lymphocytes Abs (10^9/L)",
#'   "Lymphocyte count increased",  2,      1,        4,      "Lymphocytes Abs (10^9/L)",
#'   "GGT increased",               140,    120,      180,    "Gamma Glutamyl Transferase (U/L)"
#' )
#'
#' derive_var_atoxgr_dir(data,
#'   new_var = ATOXGRH,
#'   tox_description_var = ATOXDSCH,
#'   meta_criteria = atoxgr_criteria_ctcv4,
#'   criteria_direction = "H",
#'   get_unit_expr = extract_unit(PARAM)
#' )
derive_var_atoxgr_dir <- function(dataset,
                                  new_var,
                                  tox_description_var,
                                  meta_criteria = atoxgr_criteria_ctcv4,
                                  criteria_direction,
                                  get_unit_expr) {
  new_var <- assert_symbol(enquo(new_var))
  tox_description_var <- assert_symbol(enquo(tox_description_var))
  get_unit_expr <- assert_expr(enquo(get_unit_expr))

  # check input parameter has correct value
  assert_character_scalar(criteria_direction, values = c("L", "H"))

  # Check Grade description variable exists on input data set
  assert_data_frame(dataset, required_vars = vars(!!tox_description_var))

  # Check metadata data set has required variables
  assert_data_frame(
    meta_criteria,
    required_vars = vars(TERM, GRADE_CRITERIA_CODE, DIRECTION, SI_UNIT_CHECK, VAR_CHECK)
  )
  # check DIRECTION has expected values L or H
  assert_character_vector(meta_criteria$DIRECTION, values = c("L", "H"))


  # Get list of terms from criteria metadata with particular direction
  # L = low (Hypo) H = high (Hyper)
  atoxgr_dir <- meta_criteria %>%
    filter(!is.na(GRADE_CRITERIA_CODE) & toupper(DIRECTION) == toupper(criteria_direction)) %>%
    select(TERM, DIRECTION, SI_UNIT_CHECK, GRADE_CRITERIA_CODE, VAR_CHECK) %>%
    mutate(
      TERM_UPPER = toupper(TERM),
      SI_UNIT_UPPER = toupper(SI_UNIT_CHECK)
    )

  # from ADLB VAD get distinct list of terms to be graded
  terms_in_vad <- dataset %>%
    filter(!is.na(!!tox_description_var)) %>%
    distinct(!!tox_description_var) %>%
    mutate(
      TERM = !!tox_description_var,
      TERM_UPPER = toupper(TERM)
    )

  # only keep terms that exist in both ADLB data and criteria metadata
  list_of_terms <- terms_in_vad %>%
    semi_join(atoxgr_dir, by = "TERM_UPPER") %>%
    arrange(TERM)

  # output lab data not to be graded
  # this will be appended to in for loop after each term is graded
  out_data <- dataset %>%
    filter(!!tox_description_var %notin% (list_of_terms$TERM) | is.na(!!tox_description_var)) %>%
    mutate(!!new_var := NA_character_)

  # get lab data to be graded
  to_be_graded <- dataset %>%
    filter(!!tox_description_var %in% (list_of_terms$TERM))

  # for each TERM apply criteria and create grade derivation
  for (i in seq_along(list_of_terms$TERM)) {

    # filter metadata on a term
    meta_this_term <- atoxgr_dir %>%
      filter(TERM_UPPER == list_of_terms$TERM_UPPER[i])

    # Put list of variables required for criteria in a vector
    list_of_vars <- gsub("\\s+", "", unlist(strsplit(meta_this_term$VAR_CHECK, ",")))

    # filter lab data on term and apply criteria to derive grade
    grade_this_term <- to_be_graded %>%
      filter(!!tox_description_var == list_of_terms$TERM[i])

    # check variables required in criteria exist on data
    assert_data_frame(grade_this_term, required_vars = vars(!!!syms(list_of_vars)))

    # apply criteria when SI unit matches
    grade_this_term <- grade_this_term %>%
      mutate(
        temp_flag = meta_this_term$SI_UNIT_UPPER == toupper(!!get_unit_expr) |
          is.na(meta_this_term$SI_UNIT_UPPER),
        !!new_var := if_else(
          temp_flag, eval(parse(text = meta_this_term$GRADE_CRITERIA_CODE)), NA_character_
        )
      ) %>%
      select(-temp_flag)

    # remove lab data just graded from data still to be graded
    to_be_graded <- to_be_graded %>%
      filter(!!tox_description_var != list_of_terms$TERM[i])

    # append lab data just graded to output data
    out_data <- bind_rows(out_data, grade_this_term)
  }

  out_data
}


#' Derive Lab High toxicity Grade 0 - 4 and Low Toxicity Grades 0 - (-4)
#'
#' @description
#'
#' Derives character lab grade based on high and low severity/toxicity grade(s).
#'
#' @param dataset Input data set
#'
#'   The columns `ATOXGRL`, `ATOXGRH` and specified by `lotox_description_var`,
#'   and `hitox_description_var` parameters are expected.
#'
#' @param lotox_description_var Variable containing the toxicity grade description
#' for low values, eg. "Anemia"
#'
#' @param hitox_description_var Variable containing the toxicity grade description
#' for low values, eg. "Hemoglobin Increased".
#'
#' @details
#' Created variable `ATOXGR` will contain values "-4", "-3", "-2", "-1" for low values
#' and "1", "2", "3", "4" for high values, and will contain "0" if value is gradable
#' and does not satisfy any of the criteria for high or low values. ATOXGR is set to
#' missing if information not available to give a grade.
#'
#' Function applies the following rules:
#' - High and low missing - overall missing
#' - Low grade not missing and > 0 - overall holds low grade
#' - High grade not missing and > 0 - overall holds high grade
#' - (Only high direction OR low direction is NORMAL) and high grade normal - overall NORMAL
#' - (Only low direction OR high direction is NORMAL) and low grade normal - overall NORMAL
#' - otherwise set to missing
#'
#' @author Gordon Miller
#'
#' @return The input data set with the character variable added
#'
#' @keywords der_bds_findings
#'
#' @family der_bds_findings
#'
#' @export
#'
#' @examples
#' library(dplyr, warn.conflicts = FALSE)
#'
#' adlb <- tibble::tribble(
#'   ~ATOXDSCL,          ~ATOXDSCH,        ~ATOXGRL,      ~ATOXGRH,
#'   "Hypoglycemia",     "Hyperglycemia",  NA_character_, "0",
#'   "Hypoglycemia",     "Hyperglycemia",  "0",           "1",
#'   "Hypoglycemia",     "Hyperglycemia",  "0",           "0",
#'   NA_character_,      "INR Increased",  NA_character_, "0",
#'   "Hypophosphatemia", NA_character_,    "1",           NA_character_
#' )
#'
#' derive_var_atoxgr(adlb)
derive_var_atoxgr <- function(dataset,
                              lotox_description_var = ATOXDSCL,
                              hitox_description_var = ATOXDSCH) {
  lotox_description_var <- assert_symbol(enquo(lotox_description_var))
  hitox_description_var <- assert_symbol(enquo(hitox_description_var))

  assert_data_frame(
    dataset,
    required_vars = vars(
      !!lotox_description_var,
      ATOXGRL,
      !!hitox_description_var,
      ATOXGRH
    )
  )

  lowgrade_char <- unique(dataset$ATOXGRL)
  assert_character_vector(lowgrade_char, values = c("0", "1", "2", "3", "4", NA_character_))

  highgrade_char <- unique(dataset$ATOXGRH)
  assert_character_vector(highgrade_char, values = c("0", "1", "2", "3", "4", NA_character_))


  # High and low missing - overall missing
  # Low grade not missing and > 0 - overall holds low grade
  # High grade not missing and > 0 - overall holds high grade
  # (Only high direction OR low direction is NORMAL) and high grade normal - overall NORMAL
  # (Only low direction OR high direction is NORMAL) and low grade normal - overall NORMAL
  # otherwise set to missing

  dataset %>%
    mutate(ATOXGR = case_when(
      is.na(ATOXGRL) & is.na(ATOXGRH) ~ NA_character_,
      !is.na(ATOXGRL) & ATOXGRL >= "1" ~ paste0("-", ATOXGRL),
      !is.na(ATOXGRH) & ATOXGRH >= "1" ~ ATOXGRH,
      (ATOXGRL == "0" | is.na(!!lotox_description_var)) & ATOXGRH == "0" ~ "0",
      (ATOXGRH == "0" | is.na(!!hitox_description_var)) & ATOXGRL == "0" ~ "0",
      TRUE ~ NA_character_
    ))
}
