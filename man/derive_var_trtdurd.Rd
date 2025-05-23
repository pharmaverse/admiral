% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/derive_var_trtdurd.R
\name{derive_var_trtdurd}
\alias{derive_var_trtdurd}
\title{Derive Total Treatment Duration (Days)}
\usage{
derive_var_trtdurd(dataset, start_date = TRTSDT, end_date = TRTEDT)
}
\arguments{
\item{dataset}{Input dataset

The variables specified by the \code{start_date} and \code{end_date} arguments are expected to be in the dataset.

\describe{
\item{Default value}{none}
}}

\item{start_date}{The start date

A date or date-time object is expected.

Refer to \code{derive_vars_dt()} to impute and derive a date from a date
character vector to a date object.

\describe{
\item{Default value}{\code{TRTSDT}}
}}

\item{end_date}{The end date

A date or date-time object is expected.

Refer to \code{derive_vars_dt()} to impute and derive a date from a date
character vector to a date object.

\describe{
\item{Default value}{\code{TRTEDT}}
}}
}
\value{
The input dataset with \code{TRTDURD} added
}
\description{
Derives total treatment duration (days) (\code{TRTDURD}).

\strong{Note:} This is a wrapper function for the more generic \code{derive_vars_duration()}.
}
\details{
The total treatment duration is derived as the number of days from
start to end date plus one.
}
\examples{
library(tibble)
library(lubridate)

data <- tribble(
  ~TRTSDT, ~TRTEDT,
  ymd("2020-01-01"), ymd("2020-02-24")
)

derive_var_trtdurd(data)
}
\seealso{
\code{\link[=derive_vars_duration]{derive_vars_duration()}}

Date/Time Derivation Functions that returns variable appended to dataset: 
\code{\link{derive_vars_dt}()},
\code{\link{derive_vars_dtm}()},
\code{\link{derive_vars_dtm_to_dt}()},
\code{\link{derive_vars_dtm_to_tm}()},
\code{\link{derive_vars_duration}()},
\code{\link{derive_vars_dy}()}
}
\concept{der_date_time}
\keyword{der_date_time}
\keyword{der_gen}
