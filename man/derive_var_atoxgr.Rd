% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/derive_var_atoxgr.R
\name{derive_var_atoxgr}
\alias{derive_var_atoxgr}
\title{Derive Lab High toxicity Grade 0 - 4 and Low Toxicity Grades 0 - (-4)}
\usage{
derive_var_atoxgr(
  dataset,
  lotox_description_var = ATOXDSCL,
  hitox_description_var = ATOXDSCH
)
}
\arguments{
\item{dataset}{Input dataset

The variables specified by the \code{lotox_description_var} and \code{hitox_description_var} arguments are expected to be in the dataset.
\code{ATOXGRL}, and \code{ATOXGRH} are expected as well.

\describe{
\item{Default value}{none}
}}

\item{lotox_description_var}{Variable containing the toxicity grade description
for low values, eg. "Anemia"

\describe{
\item{Default value}{\code{ATOXDSCL}}
}}

\item{hitox_description_var}{Variable containing the toxicity grade description
for high values, eg. "Hemoglobin Increased".

\describe{
\item{Default value}{\code{ATOXDSCH}}
}}
}
\value{
The input data set with the character variable added
}
\description{
Derives character lab grade based on high and low severity/toxicity grade(s).
}
\details{
Created variable \code{ATOXGR} will contain values "-4", "-3", "-2", "-1" for low values
and "1", "2", "3", "4" for high values, and will contain "0" if value is gradable
and does not satisfy any of the criteria for high or low values. ATOXGR is set to
missing if information not available to give a grade.

Function applies the following rules:
\itemize{
\item High and low missing - overall missing
\item Low grade not missing and > 0 - overall holds low grade
\item High grade not missing and > 0 - overall holds high grade
\item (Only high direction OR low direction is NORMAL) and high grade normal - overall NORMAL
\item (Only low direction OR high direction is NORMAL) and low grade normal - overall NORMAL
\item otherwise set to missing
}
}
\examples{
library(tibble)

adlb <- tribble(
  ~ATOXDSCL,          ~ATOXDSCH,        ~ATOXGRL,      ~ATOXGRH,
  "Hypoglycemia",     "Hyperglycemia",  NA_character_, "0",
  "Hypoglycemia",     "Hyperglycemia",  "0",           "1",
  "Hypoglycemia",     "Hyperglycemia",  "0",           "0",
  NA_character_,      "INR Increased",  NA_character_, "0",
  "Hypophosphatemia", NA_character_,    "1",           NA_character_
)

derive_var_atoxgr(adlb)
}
\seealso{
BDS-Findings Functions that returns variable appended to dataset: 
\code{\link{derive_basetype_records}()},
\code{\link{derive_var_analysis_ratio}()},
\code{\link{derive_var_anrind}()},
\code{\link{derive_var_atoxgr_dir}()},
\code{\link{derive_var_base}()},
\code{\link{derive_var_chg}()},
\code{\link{derive_var_ontrtfl}()},
\code{\link{derive_var_pchg}()},
\code{\link{derive_var_shift}()},
\code{\link{derive_vars_crit_flag}()}
}
\concept{der_bds_findings}
\keyword{der_bds_findings}
