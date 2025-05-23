% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/derive_param_exposure.R
\name{derive_param_exposure}
\alias{derive_param_exposure}
\title{Add an Aggregated Parameter and Derive the Associated Start and End Dates}
\usage{
derive_param_exposure(
  dataset = NULL,
  dataset_add,
  by_vars,
  input_code,
  filter_add = NULL,
  set_values_to = NULL
)
}
\arguments{
\item{dataset}{Input dataset

The variables specified by the \code{by_vars} argument are expected to be in the dataset.

\describe{
\item{Default value}{\code{NULL}}
}}

\item{dataset_add}{Additional dataset

The variables specified for \code{by_vars}, \code{analysis_var}, \code{PARAMCD},
alongside either \code{ASTDTM} and \code{AENDTM} or \code{ASTDT} and \code{AENDT} are also expected.
Observations from the specified dataset are going to be used to calculate and added
as new records to the input dataset (\code{dataset}).

\describe{
\item{Default value}{none}
}}

\item{by_vars}{Grouping variables

For each group defined by \code{by_vars} an observation is added to the output
dataset. Only variables specified in \code{by_vars} will be populated
in the newly created records.

\describe{
\item{Default value}{none}
}}

\item{input_code}{Required parameter code

The observations where \code{PARAMCD} equals the specified value are considered to compute the
summary record.

\describe{
\item{Permitted values}{A character of \code{PARAMCD} value}
\item{Default value}{none}
}}

\item{filter_add}{Filter condition as logical expression to apply during
summary calculation. By default, filtering expressions are computed within
\code{by_vars} as this will help when an aggregating, lagging, or ranking
function is involved.

For example,
\itemize{
\item \code{filter_add = (AVAL > mean(AVAL, na.rm = TRUE))} will filter all \code{AVAL}
values greater than mean of \code{AVAL} with in \code{by_vars}.
\item \code{filter_add = (dplyr::n() > 2)} will filter n count of \code{by_vars} greater
than 2.
}

\describe{
\item{Default value}{\code{NULL}}
}}

\item{set_values_to}{Variable-value pairs

Set a list of variables to some specified value for the new observation(s)
\itemize{
\item LHS refer to a variable. It is expected that at least \code{PARAMCD} is defined.
\item RHS refers to the values to set to the variable. This can be a string, a symbol, a numeric
value, \code{NA}, or an expression.
(e.g.  \code{exprs(PARAMCD = "TDOSE",PARCAT1 = "OVERALL")}).
}

\describe{
\item{Permitted values}{List of variable-value pairs}
\item{Default value}{\code{NULL}}
}}
}
\value{
The input dataset with a new record added for each group (as defined
by \code{by_vars} parameter). That is, a variable will only
be populated in this new record if it is specified in \code{by_vars}.
For each new record,
\itemize{
\item \code{set_values_to} lists each specified variable and computes its value,
\item the variable(s) specified on the LHS of \code{set_values_to} are set to their paired value (RHS).
In addition, the start and end date are computed as the minimum/maximum dates by \code{by_vars}.
}

If the input datasets contains
\itemize{
\item both \code{AxxDTM} and \code{AxxDT} then all \code{ASTDTM},\code{AENDTM}, \code{ASTDT}, \code{AENDT} are computed
\item only \code{AxxDTM} then \code{ASTDTM},\code{AENDTM} are computed
\item only \code{AxxDT} then \code{ASTDT},\code{AENDT} are computed.
}
}
\description{
Add a record computed from the aggregated analysis value of another parameter and compute the
start (\code{ASTDT(M)})and end date (\code{AENDT(M)}) as the minimum and maximum date by \code{by_vars}.
}
\details{
For each group (with respect to the variables specified for the \code{by_vars} parameter),
an observation is added to the output dataset and the defined values are set to the defined
variables
}
\examples{
library(tibble)
library(dplyr, warn.conflicts = FALSE)
library(lubridate, warn.conflicts = FALSE)
library(stringr, warn.conflicts = FALSE)
adex <- tribble(
  ~USUBJID, ~PARAMCD, ~AVAL, ~AVALC, ~VISIT, ~ASTDT, ~AENDT,
  "1015", "DOSE", 80, NA_character_, "BASELINE", ymd("2014-01-02"), ymd("2014-01-16"),
  "1015", "DOSE", 85, NA_character_, "WEEK 2", ymd("2014-01-17"), ymd("2014-06-18"),
  "1015", "DOSE", 82, NA_character_, "WEEK 24", ymd("2014-06-19"), ymd("2014-07-02"),
  "1015", "ADJ", NA, NA_character_, "BASELINE", ymd("2014-01-02"), ymd("2014-01-16"),
  "1015", "ADJ", NA, NA_character_, "WEEK 2", ymd("2014-01-17"), ymd("2014-06-18"),
  "1015", "ADJ", NA, NA_character_, "WEEK 24", ymd("2014-06-19"), ymd("2014-07-02"),
  "1017", "DOSE", 80, NA_character_, "BASELINE", ymd("2014-01-05"), ymd("2014-01-19"),
  "1017", "DOSE", 50, NA_character_, "WEEK 2", ymd("2014-01-20"), ymd("2014-05-10"),
  "1017", "DOSE", 65, NA_character_, "WEEK 24", ymd("2014-05-10"), ymd("2014-07-02"),
  "1017", "ADJ", NA, NA_character_, "BASELINE", ymd("2014-01-05"), ymd("2014-01-19"),
  "1017", "ADJ", NA, "ADVERSE EVENT", "WEEK 2", ymd("2014-01-20"), ymd("2014-05-10"),
  "1017", "ADJ", NA, NA_character_, "WEEK 24", ymd("2014-05-10"), ymd("2014-07-02")
) \%>\%
  mutate(ASTDTM = ymd_hms(paste(ASTDT, "00:00:00")), AENDTM = ymd_hms(paste(AENDT, "00:00:00")))

# Cumulative dose
adex \%>\%
  derive_param_exposure(
    dataset_add = adex,
    by_vars = exprs(USUBJID),
    set_values_to = exprs(
      PARAMCD = "TDOSE",
      PARCAT1 = "OVERALL",
      AVAL = sum(AVAL, na.rm = TRUE)
    ),
    input_code = "DOSE"
  ) \%>\%
  select(-ASTDTM, -AENDTM)

# average dose in w2-24
adex \%>\%
  derive_param_exposure(
    dataset_add = adex,
    by_vars = exprs(USUBJID),
    filter_add = VISIT \%in\% c("WEEK 2", "WEEK 24"),
    set_values_to = exprs(
      PARAMCD = "AVDW224",
      PARCAT1 = "WEEK2-24",
      AVAL = mean(AVAL, na.rm = TRUE)
    ),
    input_code = "DOSE"
  ) \%>\%
  select(-ASTDTM, -AENDTM)

# Any dose adjustment?
adex \%>\%
  derive_param_exposure(
    dataset_add = adex,
    by_vars = exprs(USUBJID),
    set_values_to = exprs(
      PARAMCD = "TADJ",
      PARCAT1 = "OVERALL",
      AVALC = if_else(sum(!is.na(AVALC)) > 0, "Y", NA_character_)
    ),
    input_code = "ADJ"
  ) \%>\%
  select(-ASTDTM, -AENDTM)
}
\seealso{
BDS-Findings Functions for adding Parameters/Records: 
\code{\link{default_qtc_paramcd}()},
\code{\link{derive_expected_records}()},
\code{\link{derive_extreme_event}()},
\code{\link{derive_extreme_records}()},
\code{\link{derive_locf_records}()},
\code{\link{derive_param_bmi}()},
\code{\link{derive_param_bsa}()},
\code{\link{derive_param_computed}()},
\code{\link{derive_param_doseint}()},
\code{\link{derive_param_exist_flag}()},
\code{\link{derive_param_framingham}()},
\code{\link{derive_param_map}()},
\code{\link{derive_param_qtc}()},
\code{\link{derive_param_rr}()},
\code{\link{derive_param_wbc_abs}()},
\code{\link{derive_summary_records}()}
}
\concept{der_prm_bds_findings}
\keyword{der_prm_bds_findings}
