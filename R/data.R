#' Adverse Events Dataset
#'
#' A SDTM AE dataset from the CDISC pilot project
#'
#' @source \url{https://bitbucket.cdisc.org/projects/CED/repos/sdtm-adam-pilot-project/browse/updated-pilot-submission-package/900172/m5/datasets/cdiscpilot01/tabulations/sdtm/ae.xpt}
"ae"

#' Concomitant Medication Dataset
#'
#' A SDTM CM dataset from the CDISC pilot project
#'
#' @source \url{https://bitbucket.cdisc.org/projects/CED/repos/sdtm-adam-pilot-project/browse/updated-pilot-submission-package/900172/m5/datasets/cdiscpilot01/tabulations/sdtm/cm.xpt}
"cm"

#' Demography Dataset
#'
#' A SDTM DM dataset from the CDISC pilot project
#'
#' @source \url{https://bitbucket.cdisc.org/projects/CED/repos/sdtm-adam-pilot-project/browse/updated-pilot-submission-package/900172/m5/datasets/cdiscpilot01/tabulations/sdtm/dm.xpt}
"dm"

#' Disposition Dataset
#'
#' A SDTM DS dataset from the CDISC pilot project
#'
#' @source \url{https://bitbucket.cdisc.org/projects/CED/repos/sdtm-adam-pilot-project/browse/updated-pilot-submission-package/900172/m5/datasets/cdiscpilot01/tabulations/sdtm/ds.xpt}
"ds"

#' Exposure Dataset
#'
#' A SDTM EX dataset from the CDISC pilot project
#'
#' @source \url{https://bitbucket.cdisc.org/projects/CED/repos/sdtm-adam-pilot-project/browse/updated-pilot-submission-package/900172/m5/datasets/cdiscpilot01/tabulations/sdtm/ex.xpt}
"ex"

#' Single Dose Exposure Dataset
#'
#' A derived dataset with single dose per date.
#'
#' @source
#' Derived from the [ex] dataset using `{admiral}` and `{dplyr}` (\url{https://github.com/Roche-GSK/admiral/blob/master/inst/example_scripts/derive_single_dose.R})
"ex_single"

#' Laboratory Measurements Dataset
#'
#' A SDTM LB dataset from the CDISC pilot project
#'
#' @source \url{https://bitbucket.cdisc.org/projects/CED/repos/sdtm-adam-pilot-project/browse/updated-pilot-submission-package/900172/m5/datasets/cdiscpilot01/tabulations/sdtm/lb.xpt}
"lb"

#' Medical History Dataset
#'
#' A SDTM MH dataset from the CDISC pilot project
#'
#' @source \url{https://bitbucket.cdisc.org/projects/CED/repos/sdtm-adam-pilot-project/browse/updated-pilot-submission-package/900172/m5/datasets/cdiscpilot01/tabulations/sdtm/mh.xpt}
"mh"

#' Questionnaire Dataset
#'
#' A SDTM QS dataset from the CDISC pilot project
#'
#' @source \url{https://bitbucket.cdisc.org/projects/CED/repos/sdtm-adam-pilot-project/browse/updated-pilot-submission-package/900172/m5/datasets/cdiscpilot01/tabulations/sdtm/qs.xpt}
"qs"

#' Supplemental Adverse Events Dataset
#'
#' A SDTM SUPPAE dataset from the CDISC pilot project
#'
#' @source \url{https://bitbucket.cdisc.org/projects/CED/repos/sdtm-adam-pilot-project/browse/updated-pilot-submission-package/900172/m5/datasets/cdiscpilot01/tabulations/sdtm/suppae.xpt}
"suppae"

#' Supplemental Disposition Dataset
#'
#' A SDTM SUPPDS dataset from the CDISC pilot project
#'
#' @source \url{https://bitbucket.cdisc.org/projects/CED/repos/sdtm-adam-pilot-project/browse/updated-pilot-submission-package/900172/m5/datasets/cdiscpilot01/tabulations/sdtm/suppds.xpt}
"suppds"

#' Supplemental Demography Dataset
#'
#' A SDTM SUPPDM dataset from the CDISC pilot project
#'
#' @source \url{https://bitbucket.cdisc.org/projects/CED/repos/sdtm-adam-pilot-project/browse/updated-pilot-submission-package/900172/m5/datasets/cdiscpilot01/tabulations/sdtm/suppdm.xpt}
"suppdm"

#' Trial Design Dataset
#'
#' A SDTM TS dataset from the CDISC pilot project
#'
#' @source \url{https://bitbucket.cdisc.org/projects/CED/repos/sdtm-adam-pilot-project/browse/updated-pilot-submission-package/900172/m5/datasets/cdiscpilot01/tabulations/sdtm/ts.xpt}
"ts"

#' Vital Signs Dataset
#'
#' A SDTM VS dataset from the CDISC pilot project
#'
#' @source \url{https://bitbucket.cdisc.org/projects/CED/repos/sdtm-adam-pilot-project/browse/updated-pilot-submission-package/900172/m5/datasets/cdiscpilot01/tabulations/sdtm/vs.xpt}
"vs"

#' Subject Level Analysis Dataset
#'
#' An example subject level analysis dataset
#'
#' @source
#' Derived from the [dm] and [ds] datasets using `{admiral}` (\url{https://github.com/Roche-GSK/admiral/blob/master/inst/example_scripts/ad_adsl.R})
#'
"adsl"

#' Vital Signs Analysis Dataset
#'
#' An example vital signs analysis dataset
#'
#' @source
#' Derived from the [adsl] and [vs] datasets using `{admiral}` (\url{https://github.com/Roche-GSK/admiral/blob/master/inst/example_scripts/ad_advs.R})
#'
"advs"

#' Queries Dataset
#'
#' An example of standard query dataset to be used in deriving variables in ADAE and ADCM
#'
"queries"
