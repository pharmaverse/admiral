% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/derive_var_anrind.R
\name{derive_var_anrind}
\alias{derive_var_anrind}
\title{Derive Reference Range Indicator}
\usage{
derive_var_anrind(
  dataset,
  signif_dig = get_admiral_option("signif_digits"),
  use_a1hia1lo = FALSE
)
}
\arguments{
\item{dataset}{Input dataset
\code{ANRLO}, \code{ANRHI}, and \code{AVAL} are expected and if \code{use_a1hia1lo} is set to \code{TRUE},
\code{A1LO} and \code{A1H1} are expected as well.

\describe{
\item{Default value}{none}
}}

\item{signif_dig}{Number of significant digits to use when comparing values.

Significant digits used to avoid floating point discrepancies when comparing numeric values.
See blog: \href{https://pharmaverse.github.io/blog/posts/2023-10-30_floating_point/floating_point.html}{How admiral handles floating points}

\describe{
\item{Default value}{\code{get_admiral_option("signif_digits")}}
}}

\item{use_a1hia1lo}{A logical value indicating whether to use \code{A1H1} and \code{A1LO} in
the derivation of \code{ANRIND}.

\describe{
\item{Default value}{\code{FALSE}}
}}
}
\value{
The input dataset with additional column \code{ANRIND}
}
\description{
Derive Reference Range Indicator
}
\details{
In the case that \code{A1H1} and \code{A1LO} are to be used, \code{ANRIND} is set to:
\itemize{
\item \code{"NORMAL"} if \code{AVAL} is greater or equal \code{ANRLO} and less than
or equal \code{ANRHI}; or if \code{AVAL} is greater than or equal \code{ANRLO} and \code{ANRHI}
is missing; or if \code{AVAL} is less than or equal \code{ANRHI} and \code{ANRLO} is
missing
\item \code{"LOW"} if \code{AVAL} is less than \code{ANRLO} and either \code{A1LO} is missing or \code{AVAL}
is greater than or equal \code{A1LO}
\item \code{"HIGH"} if \code{AVAL} is greater than \code{ANRHI} and either \code{A1HI} is missing or \code{AVAL}
is less than or equal \code{A1HI}
\item \code{"LOW LOW"} if \code{AVAL} is less than \code{A1LO}
\item \code{"HIGH HIGH"} if \code{AVAL} is greater than \code{A1HI}
}

In the case that \code{A1H1} and \code{A1LO} are not to be used, \code{ANRIND} is set to:
\itemize{
\item \code{"NORMAL"} if \code{AVAL} is greater or equal \code{ANRLO} and less than
or equal \code{ANRHI}; or if \code{AVAL} is greater than or equal \code{ANRLO} and \code{ANRHI}
is missing; or if \code{AVAL} is less than or equal \code{ANRHI} and \code{ANRLO} is
missing
\item \code{"LOW"} if \code{AVAL} is less than \code{ANRLO}
\item \code{"HIGH"} if \code{AVAL} is greater than \code{ANRHI}
}
}
\examples{
library(tibble)
library(dplyr, warn.conflicts = FALSE)

vs <- tibble::tribble(
  ~USUBJID, ~PARAMCD, ~AVAL, ~ANRLO, ~ANRHI, ~A1LO, ~A1HI,
  "P01",       "PUL",    70,     60,    100,    40,   110,
  "P01",       "PUL",    57,     60,    100,    40,   110,
  "P01",       "PUL",    60,     60,    100,    40,   110,
  "P01",     "DIABP",   102,     60,     80,    40,    90,
  "P02",       "PUL",   109,     60,    100,    40,   110,
  "P02",       "PUL",   100,     60,    100,    40,   110,
  "P02",     "DIABP",    80,     60,     80,    40,    90,
  "P03",       "PUL",    39,     60,    100,    40,   110,
  "P03",       "PUL",    40,     60,    100,    40,   110
)

vs \%>\% derive_var_anrind(use_a1hia1lo = TRUE)
vs \%>\% derive_var_anrind(use_a1hia1lo = FALSE)

}
\seealso{
BDS-Findings Functions that returns variable appended to dataset: 
\code{\link{derive_basetype_records}()},
\code{\link{derive_var_analysis_ratio}()},
\code{\link{derive_var_atoxgr}()},
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
