#' Derive Lab toxicity grade 0 - 4
#'
#' Derives a character variable containing lab grade based on severity/toxicity criteria.
#'
#' @param dataset Input data set
#'
#'   The columns specified by `tox_description_var` parameter is expected.
#'
#' @param new_var Name of the character grade variable to create.
#'
#' @param tox_description_var Variable containing the description of the grading
#' criteria eg. Anemia.
#'
#' @param meta_criteria metadata data set holding the criteria (normally a case statement)
#'
#'   The metadata should have the following variables:
#'
#' - `TERM`: variable to hold the term describing the criteria applied to a particular lab test,
#'   eg. 'Anemia' or 'INR Increased'. Note: the variable is case insensitive.
#' - `DIRECTION`: variable to hold the direction of the abnormality of a particular lab test
#'   value. 'L' is for LOW values, 'H' is for HIGH values. Note: the variable is case insensitive.
#' - `SI_UNIT_CHECK`: variable to hold unit of particular lab test. Used to check against input data
#'   if criteria is based on absolute values.
#' - `IMPLEMENTATION`: variable to hold case statement to create grade based on defined criteria.
#'
#' @param criteria_direction Direction (L= Low, H = High) of toxicity grade.
#'
#' @details
#' `new_var` is derived with values NA, '0', '1', '2', '3', '4', where '4' is the most
#' severe grade
#' - '4' is where the lab value satisfies the criteria for grade 4.
#' - '3' is where the lab value satisfies the criteria for grade 3.
#' - '2' is where the lab value satisfies the criteria for grade 2.
#' - '1' is where the lab value satisfies the criteria for grade 1.
#' - '0' is where a grade can be derived and is not grade '1', '2', '3' or '4'.
#' - NA is where a grade cannot be derived.
#'
#' @author Gordon Miller
#'
#' @return The input dataset with the character grade added
#'
#' @keywords adam bds adlb derivation
#'
#' @export
#'
#' @examples
#' library(dplyr, warn.conflicts = FALSE)
#'
#' data <- tibble::tribble(
#'   ~ATOXDSCL,                     ~AVAL,  ~ANRLO,   ~ANRHI,
#'   "Hypoglycemia",                119,    4,        7,
#'   "Hypoglycemia",                120,    4,        7,
#'   "Anemia",                      129,    120,      180,
#'   "White blood cell decreased",  10,     5,        20,
#'   "White blood cell decreased",   15,     5,        20,
#'   "Anemia",                      140,    120,      180,
#' )
#'
#' derive_var_atoxgr_dir(adlb,
#'   new_var = ATOXGRL,
#'   tox_description_var = ATOXDSCL,
#'   meta_criteria = atoxgr_criteria_ctcv4,
#'   criteria_direction = "L"
#' )
#'
#' data <- tibble::tribble(
#'   ~ATOXDSCH,                     ~AVAL,  ~ANRLO,   ~ANRHI,
#'   "Hyperglycemia",               119,    4,        7,
#'   "Hyperglycemia",               120,    4,        7,
#'   "Hemoglobin increased",        129,    120,      180,
#'   "Lymphocyte count increased",  4,      1,        4,
#'   "Lymphocyte count increased",  2,      1,        4,
#'   "Hemoglobin increased",        140,    120,      180,
#' )
#'
#' derive_var_atoxgr_dir(adlb,
#'   new_var = ATOXGRH,
#'   tox_description_var = ATOXDSCH,
#'   meta_criteria = atoxgr_criteria_ctcv4,
#'   criteria_direction = "H"
#' )


derive_var_atoxgr_dir <- function(dataset,
                                  new_var,
                                  tox_description_var,
                                  meta_criteria,
                                  criteria_direction) {

  new_var <- assert_symbol(enquo(new_var))
  tox_description_var <- assert_symbol(enquo(tox_description_var))


  # Get list of terms from criteria metadata with particular direction
  # L = low (Hypo) H = high (Hyper)
  atoxgr_dir <- meta_criteria %>%
    filter(!is.na(IMPLEMENTATION) & toupper(DIRECTION) == toupper(criteria_direction)
    ) %>%
    select(TERM, DIRECTION, SI_UNIT_CHECK, IMPLEMENTATION
    ) %>%
    mutate(TERM_UPPER = toupper(TERM),
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

  # Output terms in ADLB VAD that are not in criteria metadata
  # put out a WARNING with list of terms - to be done

  # output lab data not to be graded
  # this will be appended to in for loop after each term is graded
  out_data <- dataset %>%
    filter(!!tox_description_var %notin% (list_of_terms$TERM) | is.na(!!tox_description_var)
    )

  # get lab data to be graded
  to_be_graded <- dataset %>%
    filter(!!tox_description_var %in% (list_of_terms$TERM)
    )

  # for each TERM apply criteria and create grade derivation
  for (i in seq_along(list_of_terms$TERM)) {

    # filter metadata on a term
    meta_this_term <-  atoxgr_dir %>%
      filter(TERM_UPPER == list_of_terms$TERM_UPPER[i]
      )

    # filter lab data on term and apply criteria to derive grade
    grade_this_term <- to_be_graded %>%
      filter(!!tox_description_var == list_of_terms$TERM[i]) %>%
      mutate(
        TEMP_VAR = meta_this_term$SI_UNIT_UPPER == toupper(AVALU) |
          is.na(meta_this_term$SI_UNIT_UPPER),
        !!new_var := if_else(
          TEMP_VAR, eval(parse(text = meta_this_term$IMPLEMENTATION)), NA_character_
        )
      ) %>%
      select(-TEMP_VAR)

    # remove lab data just graded from data still to be graded
    to_be_graded <- to_be_graded %>%
      filter(!!tox_description_var != list_of_terms$TERM[i])

    # append lab data just graded to output data
    out_data <- bind_rows(out_data, grade_this_term)

  }

  all_data <- out_data

}


#' Derive Lab High toxicity grade 0 - 4 and Low toxicity grades 0 - (-4)
#'
#' Derives character variable containing lab grade based on high and low
#' severity/toxicity grade(s).
#'
#' @param dataset Input data set
#'
#'   The columns specified by `lotox_description_var`, `lotox_var`,
#'   `hitox_description_var` and `hitox_var` parameters are expected.
#'
#' @param new_var Name of the character grade variable to create. Will contain
#' values "-4", "-3", "-2", "-1" for low values and "1", "2", "3", "4" for high values.
#' Contains 0 if value is gradable and does not satisfy any of the criteria for high or low
#' values. Is set to missing if information not available to give a grade.
#'
#' @param lotox_description_var Variable containing the description of the grading
#' criteria eg. Anemia. Grade description for low values.
#'
#' @param lotox_var Variable containing the grading for low values.
#'
#' @param hitox_description_var Variable containing the description of the grading
#' criteria eg. Hemoglobin Increased. Grade description for high values.
#'
#' @param hitox_var Variable containing the grading for high values.
#'
#' @author Gordon Miller
#'
#' @return The input data set with the character grades (low and high) added
#'
#' @keywords adam bds adlb derivation
#'
#' @export
#'
#' @examples
#' library(dplyr, warn.conflicts = FALSE)
#'
#' data <- tibble::tribble(
#'   ~ATOXDSCL,          ~ATOXDSCH,        ~ATOXGRL,      ~ATOXGRH,
#'   "Hypoglycemia",     "Hyperglycemia",  NA_character_, "0",
#'   "Hypoglycemia",     "Hyperglycemia",  "0",           "1",
#'   "Hypoglycemia",     "Hyperglycemia",  "0",           "0",
#'   NA_character_,      "INR Increased",  NA_character_, "0",
#'   "Hypophosphatemia", NA_character_,    "1",           NA_character_,
#' )
#'
#' derive_var_atoxgr(adlb)


derive_vars_atoxgr <- function(dataset,
                               new_var = ATOXGR,
                               lotox_description_var = ATOXDSCL,
                               lotox_var = ATOXGRL,
                               hitox_description_var = ATOXDSCH,
                               hitox_var = ATOXGRH
                               ) {

  new_var <- assert_symbol(enquo(new_var))
  lotox_description_var <- assert_symbol(enquo(lotox_description_var))
  lotox_var <- assert_symbol(enquo(lotox_var))
  hitox_description_var <- assert_symbol(enquo(hitox_description_var))
  hitox_var <- assert_symbol(enquo(hitox_var))

  # High and low missing - overall missing
  # Low grade not missing and > 0 - overall holds low grade
  # High grade not missing and > 0 - overall holds high grade
  # (Only high direction OR low direction is NORMAL) and high grade normal - overall NORMAL
  # (Only low direction OR high direction is NORMAL) and low grade normal - overall NORMAL
  # otherwise set to missing

  dataset %>%
    mutate(!!new_var := case_when(
      is.na(!!lotox_var) & is.na(!!hitox_var) ~ NA_character_,
      !is.na(!!lotox_var) & !!lotox_var >= "1" ~ paste0('-', !!ATOXGRL),
      !is.na(!!hitox_var) & !!hitox_var >= "1" ~ !!hitox_var,
      (!!lotox_var == "0" | is.na(!!lotox_description_var)) & !!hitox_var == "0" ~ "0",
      (!!hitox_var == "0" | is.na(!!hitox_description_var)) & !!lotox_var == "0" ~ "0",
      TRUE ~ NA_character_)
    )

}
