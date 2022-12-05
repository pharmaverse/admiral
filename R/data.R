#' Single Dose Exposure Dataset
#'
#' A derived dataset with single dose per date.
#' @keywords datasets
#' @family datasets
#' @source
#' Derived from the `ex` dataset using `{admiral}` and `{dplyr}` (\url{https://github.com/pharmaverse/admiral/blob/main/inst/example_scripts/derive_single_dose.R})
"ex_single"

#' Queries Dataset
#' @keywords datasets
#' @family datasets
#' @source
#' An example of standard query dataset to be used in deriving variables in ADAE and ADCM
#'
"queries"

#' Queries MH Dataset
#' @keywords datasets
#' @family datasets
#' @source
#' An example of standard query MH dataset to be used in deriving variables in ADMH
#'
"queries_mh"

#' Subject Level Analysis Dataset
#'
#' An example subject level analysis dataset
#' @keywords datasets
#' @family datasets
#' @source
#' Derived from the `dm` and `ds` datasets using `{admiral}` (\url{https://github.com/pharmaverse/admiral/blob/main/inst/templates/ad_adsl.R})
#'
"admiral_adsl"

#' Metadata Holding Grading Criteria for NCI-CTCAEv4
#'
#' @details
#' This metadata has its origin in the ADLB Grading Spec Excel file which ships with `{admiral}`
#' and can be accessed using `system.file("adlb_grading/adlb_grading_spec.xlsx", package = "admiral")`
#' in sheet = "NCICTCAEv4".
#' The dataset contained in there has the following columns:
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
#' [**Common Terminology Criteria for Adverse Events (CTCAE)v4.0**](https://ctep.cancer.gov/protocoldevelopment/electronic_applications/ctc.htm#ctc_40)
#' From these variables only 'TERM' is used in the {admiral} code, the rest are for information and
#' tracability only.
#'
#' @author Gordon Miller
#'
#' @keywords metadata
#' @family metadata
"atoxgr_criteria_ctcv4"

#' Metadata Holding Grading Criteria for NCI-CTCAEv5
#'
#' @details
#' This metadata has its origin in the ADLB Grading Spec Excel file which ships with `{admiral}`
#' and can be accessed using `system.file("adlb_grading/adlb_grading_spec.xlsx", package = "admiral")`
#' in sheet = "NCICTCAEv5".
#' The dataset contained in there has the following columns:
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
#' [**Common Terminology Criteria for Adverse Events (CTCAE)v5.0**](https://ctep.cancer.gov/protocoldevelopment/electronic_applications/ctc.htm#ctc_50)
#' From these variables only 'TERM' is used in the {admiral} code, the rest are for information and
#' traceability only.
#'
#' @author Gordon Miller
#'
#' @keywords metadata
#' @family metadata
"atoxgr_criteria_ctcv5"
