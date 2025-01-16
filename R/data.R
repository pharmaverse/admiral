#' Single Dose Exposure Dataset
#'
#' A derived dataset with single dose per date.
#' @keywords datasets
#' @family datasets
#' @source
#' Derived from the `ex` dataset using `{admiral}` and `{dplyr}` (\url{https://github.com/pharmaverse/admiral/blob/main/data-raw/create_ex_single.R})
"ex_single"

#' Example `QS` Dataset
#'
#' An example `QS` dataset based on the examples from the CDISC ADaM Supplements
#' [Generalized Anxiety Disorder 7-Item Version 2
#' (GAD-7)](https://www.cdisc.org/standards/foundational/qrs/generalized-anxiety-disorder-7-item-version-2-0)
#' and [Geriatric Depression Scale Short Form
#' (GDS-SF)](https://www.cdisc.org/standards/foundational/qrs/geriatric-depression-scale-short-form-0).
#' @keywords datasets
#' @family datasets
#' @source
#' Created by (\url{https://github.com/pharmaverse/admiral/blob/main/data-raw/create_example_qs.R})
"example_qs"

#' Queries Dataset
#' @keywords datasets
#' @family datasets
#' @source
#' An example of standard query dataset to be used in deriving Standardized
#' MedDRA Query variables in ADAE
#'
"queries"

#' Queries MH Dataset
#' @keywords datasets
#' @family datasets
#' @source
#' An example of standard query MH dataset to be used in deriving Standardized
#' MedDRA Query variables in ADMH
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

#' Lab Analysis Dataset
#'
#' An example of lab analysis dataset
#' @keywords datasets
#' @family datasets
#' @source
#' Derived from the `adlb` template, then further filtered
#' due to dataset size by the following USUBJIDs:
#' 01-701-1015, 01-701-1023, 01-701-1028, 01-701-1033,
#' 01-701-1034, 01-701-1047, 01-701-1097, 01-705-1186,
#' 01-705-1292, 01-705-1310, 01-708-1286
#'
"admiral_adlb"

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
#' From these variables only 'TERM' is used in the `{admiral}` code, the rest are for information and
#' traceability only.
#'
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
#' From these variables only 'TERM' is used in the `{admiral}` code, the rest are for information and
#' traceability only.
#'
#'
#' @keywords metadata
#' @family metadata
"atoxgr_criteria_ctcv5"

#' Metadata Holding Grading Criteria for DAIDs
#'
#' @details
#' This metadata has its origin in the ADLB Grading Spec Excel file which ships with `{admiral}`
#' and can be accessed using `system.file("adlb_grading/adlb_grading_spec.xlsx", package = "admiral")`
#' in sheet = "DAIDS".
#' The dataset contained in there has the following columns:
#' - `SOC`: variable to hold the SOC of the lab test criteria.
#' - `TERM`: variable to hold the term describing the criteria applied to a particular lab test,
#'   eg. 'Anemia' or 'INR Increased'. Note: the variable is case insensitive.
#' - `SUBGROUP` : Description of sub-group of subjects were grading will be applied (i.e. >= 18 years)
#' - `Grade 1`: Criteria defining lab value as Grade 1.
#' - `Grade 2`: Criteria defining lab value as Grade 2.
#' - `Grade 3`: Criteria defining lab value as Grade 3.
#' - `Grade 4`: Criteria defining lab value as Grade 4.
#' - `Grade 5`: Criteria defining lab value as Grade 5.
#' - `Definition`: Holds the definition of the lab test abnormality.
#' - `FILTER` : `admiral` code to apply the filter based on SUBGROUP column.
#' - `GRADE_CRITERIA_CODE`: variable to hold code that creates grade based on defined criteria.
#' - `SI_UNIT_CHECK`: variable to hold unit of particular lab test. Used to check against input data
#'   if criteria is based on absolute values.
#' - `VAR_CHECK`: List of variables required to implement lab grade criteria. Use to check against
#'   input data.
#' - `DIRECTION`: variable to hold the direction of the abnormality of a particular lab test
#'   value. 'L' is for LOW values, 'H' is for HIGH values. Note: the variable is case insensitive.
#' - `COMMENT`: Holds any information regarding rationale behind implementation of grading criteria.
#'
#' Note: Variables `SOC`, `TERM`, `SUBGROUP`, `Grade 1`, `Grade 2`,`Grade 3`,`Grade 4`,`Grade 5`, `Definition`
#' are from the source document on DAIDS website defining the grading criteria.
#' [Division of AIDS (DAIDS) Table for Grading the Severity of Adult and Pediatric Adverse Events
#' From these variables only 'TERM' is used in the `{admiral}` code, the rest are for information and
#' traceability only.
#'
#'
#' @keywords metadata
#' @family metadata
"atoxgr_criteria_daids"

#' Country Code Lookup
#'
#' @description
#' These pre-defined country codes are sourced from
#' [ISO 3166 Standards](https://www.iso.org/iso-3166-country-codes.html).
#' See also [Wikipedia](https://en.wikipedia.org/wiki/ISO_3166-1_alpha-3).
#'
#' @details
#'
#' `country_code` is the 3-letter ISO 3166-1 county code commonly found in the
#' ADSL `COUNTRY` variable.
#' `country_name` is the country long name corresponding to to the 3-letter code.
#' `country_number` is the numeric code corresponding to an alphabetic sorting of
#' the 3-letter codes.
#'
#' To see the entire table in the console, run `print(country_code_lookup)`.
#'
#' @seealso [dose_freq_lookup]
#'
#' @keywords metadata
#'
#' @family metadata
#'
#' @examples
#' library(tibble)
#' library(dplyr, warn.conflicts = FALSE)
#'
#' # Create reference dataset for periods
#' adsl <- tribble(
#'   ~USUBJID, ~SEX, ~COUNTRY,
#'   "ST01-01", "F", "AUT",
#'   "ST01-02", "M", "MWI",
#'   "ST01-03", "F", "GBR",
#'   "ST01-04", "M", "CHE",
#'   "ST01-05", "M", "NOR",
#'   "ST01-06", "F", "JPN",
#'   "ST01-07", "F", "USA"
#' )
#'
#' adsl %>%
#'   derive_vars_merged(
#'     dataset_add = country_code_lookup,
#'     new_vars = exprs(COUNTRYN = country_number, COUNTRYL = country_name),
#'     by_vars = exprs(COUNTRY = country_code)
#'   )
#'
#' @rdname country_code_lookup
"country_code_lookup"

#' Pre-Defined Dose Frequencies
#'
#' @description
#' These pre-defined dose frequencies are sourced from
#' [CDISC](https://evs.nci.nih.gov/ftp1/CDISC/SDTM/SDTM%20Terminology.pdf). The
#' number of rows to generate using `create_single_dose_dataset()` arguments
#' `start_date` and `end_date` is derived from `DOSE_COUNT`, `DOSE_WINDOW`, and
#' `CONVERSION_FACTOR` with appropriate functions from `lubridate`.
#'
#' @details
#' `NCI_CODE` and `CDISC_VALUE` are included from the CDISC source for
#' traceability.
#'
#' `DOSE_COUNT` represents the number of doses received in one single unit of
#' `DOSE_WINDOW`. For example, for `CDISC_VALUE=="10 DAYS PER MONTH"`,
#' `DOSE_WINDOW=="MONTH"` and `DOSE_COUNT==10`. Similarly, for
#' `CDISC_VALUE=="EVERY 2 WEEKS"`, `DOSE_WINDOW=="WEEK"` and
#' `DOSE_COUNT==0.5` (to yield one dose every two weeks).
#'
#' `CONVERSION_FACTOR` is used to convert `DOSE_WINDOW` units `"WEEK"`,
#'  `"MONTH"`, and `"YEAR"` to the unit `"DAY"`.
#'
#' For example, for `CDISC_VALUE=="10 DAYS PER MONTH"`, `CONVERSION_FACTOR`
#' is `0.0329`. One day of a month is assumed to be `1 / 30.4375` of a month (one
#' day is assumed to be `1/365.25` of a year).
#' Given only `start_date` and `end_date` in the aggregate dataset, `CONVERSION_FACTOR`
#' is used to calculate specific dates for`start_date` and `end_date` in the
#' resulting single dose dataset for the doses that occur. In such cases, doses
#' are assumed to occur at evenly spaced increments over the interval.
#'
#' To see the entire table in the console, run `print(dose_freq_lookup)`.
#'
#' @seealso [create_single_dose_dataset()]
#'
#' @keywords metadata
#' @family metadata
#'
#' @rdname dose_freq_lookup
"dose_freq_lookup"
