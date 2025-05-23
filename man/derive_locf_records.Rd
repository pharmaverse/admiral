% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/derive_locf_records.R
\name{derive_locf_records}
\alias{derive_locf_records}
\title{Derive LOCF (Last Observation Carried Forward) Records}
\usage{
derive_locf_records(
  dataset,
  dataset_ref,
  by_vars,
  analysis_var = AVAL,
  order,
  keep_vars = NULL
)
}
\arguments{
\item{dataset}{Input dataset

The variables specified by the \code{by_vars}, \code{analysis_var}, \code{order}, and \code{keep_vars} arguments are expected to be in the dataset.

\describe{
\item{Default value}{none}
}}

\item{dataset_ref}{Expected observations dataset

Data frame with all the combinations of \code{PARAMCD}, \code{PARAM}, \code{AVISIT},
\code{AVISITN}, ... which are expected in the dataset is expected.

\describe{
\item{Default value}{none}
}}

\item{by_vars}{Grouping variables

For each group defined by \code{by_vars} those observations from \code{dataset_ref}
are added to the output dataset which do not have a corresponding observation
in the input dataset or for which \code{analysis_var} is \code{NA} for the corresponding observation
in the input dataset.

\describe{
\item{Default value}{none}
}}

\item{analysis_var}{Analysis variable.

\emph{Default}: \code{AVAL}

\describe{
\item{Permitted values}{a variable}
\item{Default value}{\code{AVAL}}
}}

\item{order}{Sort order

The dataset is sorted by \code{order} before carrying the last observation
forward (e.g. \code{AVAL}) within each \code{by_vars}.

For handling of \code{NA}s in sorting variables see \href{../articles/generic.html#sort_order}{Sort Order}.

\describe{
\item{Default value}{none}
}}

\item{keep_vars}{Variables that need carrying the last observation forward

Keep variables that need carrying the last observation forward other than \code{analysis_var}
(e.g., \code{PARAMN}, \code{VISITNUM}). If by default \code{NULL}, only variables specified in
\code{by_vars} and \code{analysis_var} will be populated in the newly created records.

\describe{
\item{Default value}{\code{NULL}}
}}
}
\value{
The input dataset with the new "LOCF" observations added for each
\code{by_vars}. Note, a variable will only be populated in the new parameter rows
if it is specified in \code{by_vars}.
}
\description{
Adds LOCF records as new observations for each 'by group' when the dataset
does not contain observations for missed visits/time points.
}
\details{
For each group (with respect to the variables specified for the
by_vars parameter) those observations from \code{dataset_ref} are added to
the output dataset
\itemize{
\item which do not have a corresponding observation in the input dataset or
\item for which \code{analysis_var} is NA for the corresponding observation in the input dataset.

For the new observations, \code{analysis_var} is set to the non-missing \code{analysis_var} of the
previous observation in the input dataset (when sorted by \code{order}) and
\code{DTYPE} is set to "LOCF".
}
}
\examples{

library(dplyr)
library(tibble)

advs <- tribble(
  ~STUDYID,  ~USUBJID,      ~PARAMCD, ~PARAMN, ~AVAL, ~AVISITN, ~AVISIT,
  "CDISC01", "01-701-1015", "PULSE",        1,    65,        0, "BASELINE",
  "CDISC01", "01-701-1015", "DIABP",        2,    79,        0, "BASELINE",
  "CDISC01", "01-701-1015", "DIABP",        2,    80,        2, "WEEK 2",
  "CDISC01", "01-701-1015", "DIABP",        2,    NA,        4, "WEEK 4",
  "CDISC01", "01-701-1015", "DIABP",        2,    NA,        6, "WEEK 6",
  "CDISC01", "01-701-1015", "SYSBP",        3,   130,        0, "BASELINE",
  "CDISC01", "01-701-1015", "SYSBP",        3,   132,        2, "WEEK 2",
  "CDISC01", "01-701-1028", "PULSE",        1,    61,        0, "BASELINE",
  "CDISC01", "01-701-1028", "PULSE",        1,    60,        6, "WEEK 6",
  "CDISC01", "01-701-1028", "DIABP",        2,    51,        0, "BASELINE",
  "CDISC01", "01-701-1028", "DIABP",        2,    50,        2, "WEEK 2",
  "CDISC01", "01-701-1028", "DIABP",        2,    51,        4, "WEEK 4",
  "CDISC01", "01-701-1028", "DIABP",        2,    50,        6, "WEEK 6",
  "CDISC01", "01-701-1028", "SYSBP",        3,   121,        0, "BASELINE",
  "CDISC01", "01-701-1028", "SYSBP",        3,   121,        2, "WEEK 2",
  "CDISC01", "01-701-1028", "SYSBP",        3,   121,        4, "WEEK 4",
  "CDISC01", "01-701-1028", "SYSBP",        3,   121,        6, "WEEK 6"
)

# A dataset with all the combinations of PARAMCD, PARAM, AVISIT, AVISITN, ... which are expected.
advs_expected_obsv <- tribble(
  ~PARAMCD, ~AVISITN, ~AVISIT,
  "PULSE",         0, "BASELINE",
  "PULSE",         6, "WEEK 6",
  "DIABP",         0, "BASELINE",
  "DIABP",         2, "WEEK 2",
  "DIABP",         4, "WEEK 4",
  "DIABP",         6, "WEEK 6",
  "SYSBP",         0, "BASELINE",
  "SYSBP",         2, "WEEK 2",
  "SYSBP",         4, "WEEK 4",
  "SYSBP",         6, "WEEK 6"
)

derive_locf_records(
  dataset = advs,
  dataset_ref = advs_expected_obsv,
  by_vars = exprs(STUDYID, USUBJID, PARAMCD),
  order = exprs(AVISITN, AVISIT),
  keep_vars = exprs(PARAMN)
) |>
  arrange(USUBJID, PARAMCD, AVISIT)

}
\seealso{
BDS-Findings Functions for adding Parameters/Records: 
\code{\link{default_qtc_paramcd}()},
\code{\link{derive_expected_records}()},
\code{\link{derive_extreme_event}()},
\code{\link{derive_extreme_records}()},
\code{\link{derive_param_bmi}()},
\code{\link{derive_param_bsa}()},
\code{\link{derive_param_computed}()},
\code{\link{derive_param_doseint}()},
\code{\link{derive_param_exist_flag}()},
\code{\link{derive_param_exposure}()},
\code{\link{derive_param_framingham}()},
\code{\link{derive_param_map}()},
\code{\link{derive_param_qtc}()},
\code{\link{derive_param_rr}()},
\code{\link{derive_param_wbc_abs}()},
\code{\link{derive_summary_records}()}
}
\author{
G Gayatri
}
\concept{der_prm_bds_findings}
\keyword{der_prm_bds_findings}
