#' Create Metadata Holding Grading Criteria for NCI-CTCAEv4
#'
#' @description
#' Reads in excel file and creates metadata data set
#'
#' @details
#' Reads in excel file containing grading criteria, the spreadsheet has the following columns:
#' - `SOC`: variable to hold the SOC of the lab test criteria.
#' - `TERM`: variable to hold the term describing the criteria applied to a particular lab test,
#'   eg. 'Anemia' or 'INR Increased'. Note: the variable is case insensitive.
#' - `Grade 1`: Criteria defining lab value as Grade 1.
#' - `Grade 2`: Criteria defining lab value as Grade 2.
#' - `Grade 3`: Criteria defining lab value as Grade 3.
#' - `Grade 4`: Criteria defining lab value as Grade 4.
#' - `Grade 5`: Criteria defining lab value as Grade 5.
#' - `Definition`: Holds the definition of the lab test abnormality.
#' - `GRADE_CRITERIA_CODE`: variable to hold code that creates grade based on defined criteria.
#' - `SI_UNIT_CHECK`: variable to hold unit of particular lab test. Used to check against input data
#'   if criteria is based on absolute values.
#' - `VAR_CHECK`: List of variables required to implement lab grade criteria. Use to check against
#'   input data.
#' - `DIRECTION`: variable to hold the direction of the abnormality of a particular lab test
#'   value. 'L' is for LOW values, 'H' is for HIGH values. Note: the variable is case insensitive.
#' - `COMMENT`: Holds any information regarding rationale behind implementation of grading criteria.
#'
#' Note: Variables `SOC`, `TERM`, `Grade 1`, `Grade 2`,`Grade 3`,`Grade 4`,`Grade 5`, `Definition`
#' are from the source document on NCI-CTC website defining the grading criteria.
#' From these variables only 'TERM' is used in the {admiral} code, the rest are for information and
#' tracability only.
#'
#' @author Gordon Miller
#'
#' @return Dataset with grading criteria
#'
#' @keywords der_bds_findings
#'
#' @family der_bds_findings
#'
#' @rdname atoxgr_sources
#'
#' @export

atoxgr_criteria_ctcv4 <- system.file("adlb_grading_spec.xlsx", package = "admiral") %>%
  read_excel(sheet = "NCICTCAEv4") %>%
  mutate(GRADE_CRITERIA_CODE = gsub("[\r\n]", " ", GRADE_CRITERIA_CODE))
