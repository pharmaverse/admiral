% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/derive_var_chg.R
\name{derive_var_chg}
\alias{derive_var_chg}
\title{Derive Change from Baseline}
\usage{
derive_var_chg(dataset)
}
\arguments{
\item{dataset}{Input dataset \code{AVAL} and \code{BASE} are expected.

\describe{
\item{Default value}{none}
}}
}
\value{
The input dataset with an additional column named \code{CHG}
}
\description{
Derive change from baseline (\code{CHG}) in a BDS dataset
}
\details{
Change from baseline is calculated by subtracting the baseline value
from the analysis value.
}
\examples{
library(tibble)

advs <- tribble(
  ~USUBJID, ~PARAMCD, ~AVAL, ~ABLFL, ~BASE,
  "P01",    "WEIGHT", 80,    "Y",    80,
  "P01",    "WEIGHT", 80.8,  NA,     80,
  "P01",    "WEIGHT", 81.4,  NA,     80,
  "P02",    "WEIGHT", 75.3,  "Y",    75.3,
  "P02",    "WEIGHT", 76,    NA,     75.3
)
derive_var_chg(advs)
}
\seealso{
BDS-Findings Functions that returns variable appended to dataset: 
\code{\link{derive_basetype_records}()},
\code{\link{derive_var_analysis_ratio}()},
\code{\link{derive_var_anrind}()},
\code{\link{derive_var_atoxgr}()},
\code{\link{derive_var_atoxgr_dir}()},
\code{\link{derive_var_base}()},
\code{\link{derive_var_ontrtfl}()},
\code{\link{derive_var_pchg}()},
\code{\link{derive_var_shift}()},
\code{\link{derive_vars_crit_flag}()}
}
\concept{der_bds_findings}
\keyword{der_bds_findings}
