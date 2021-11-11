#' Single Dose Exposure Dataset
#'
#' A derived dataset with single dose per date.
#'
#' @source
#' Derived from the [ex] dataset using `{admiral}` and `{dplyr}` (\url{https://github.com/Roche-GSK/admiral/blob/master/inst/example_scripts/derive_single_dose.R})
"ex_single"

#' Queries Dataset
#'
#' An example of standard query dataset to be used in deriving variables in ADAE and ADCM
#'
"queries"

#' Subject Level Analysis Dataset
#'
#' An example subject level analysis dataset
#'
#' @source
#' Derived from the [dm] and [ds] datasets using `{admiral}` (\url{https://github.com/Roche-GSK/admiral/blob/master/inst/example_scripts/ad_adsl.R})
#'
"adsl"

#' Adverse Event Analysis Dataset
#'
#' An example adverse event analysis dataset
#'
#' @source
#' Derived from the [adsl] and [ae] datasets using `{admiral}` (\url{https://github.com/Roche-GSK/admiral/blob/master/inst/example_scripts/ad_adae.R})
#'
"adae"

#' Concomitant Medication Analysis Dataset
#'
#' An example concomitant medication analysis dataset
#'
#' @source
#' Derived from the [adsl] and [cm] datasets using `{admiral}` (\url{https://github.com/Roche-GSK/admiral/blob/master/inst/example_scripts/ad_adcm.R})
#'
"adcm"

#' Exposure Analysis Dataset
#'
#' An example exposure analysis dataset
#'
#' @source
#' Derived from the [adsl] and [ex] datasets using `{admiral}` (\url{https://github.com/Roche-GSK/admiral/blob/master/inst/example_scripts/ad_adex.R})
#'
"adex"

#' Vital Signs Analysis Dataset
#'
#' An example vital signs analysis dataset
#'
#' @source
#' Derived from the [adsl] and [vs] datasets using `{admiral}` (\url{https://github.com/Roche-GSK/admiral/blob/master/inst/example_scripts/ad_advs.R})
#'
"advs"
