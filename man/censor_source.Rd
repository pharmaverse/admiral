% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/derive_param_tte.R
\name{censor_source}
\alias{censor_source}
\title{Create a \code{censor_source} Object}
\usage{
censor_source(
  dataset_name,
  filter = NULL,
  date,
  censor = 1,
  set_values_to = NULL,
  order = NULL
)
}
\arguments{
\item{dataset_name}{The name of the source dataset

The name refers to the dataset provided by the \code{source_datasets} parameter
of \code{derive_param_tte()}.

\describe{
\item{Default value}{none}
}}

\item{filter}{An unquoted condition for selecting the observations from
\code{dataset} which are events or possible censoring time points.

\describe{
\item{Default value}{\code{NULL}}
}}

\item{date}{A variable or expression providing the date of the event or
censoring. A date, or a datetime can be specified. An unquoted symbol or
expression is expected.

Refer to \code{derive_vars_dt()} or \code{convert_dtc_to_dt()} to impute and derive a
date from a date character vector to a date object.

\describe{
\item{Default value}{none}
}}

\item{censor}{Censoring value

CDISC strongly recommends using \code{0} for events and positive integers for
censoring.

\describe{
\item{Default value}{\code{0}}
}}

\item{set_values_to}{A named list returned by \code{exprs()} defining the variables
to be set for the event or censoring, e.g. \code{exprs(EVENTDESC = "DEATH", SRCDOM = "ADSL", SRCVAR = "DTHDT")}. The values must be a symbol, a
character string, a numeric value, an expression, or \code{NA}.

\describe{
\item{Default value}{\code{NULL}}
}}

\item{order}{Sort order

An optional named list returned by \code{exprs()} defining additional variables
that the source dataset is sorted on after \code{date}.

\describe{
\item{Permitted values}{list of variables created by \code{exprs()} e.g. \code{exprs(ASEQ)}.}
\item{Default value}{\code{order}}
}}
}
\value{
An object of class \code{censor_source}, inheriting from class \code{tte_source}
}
\description{
\code{censor_source} objects are used to define censorings as input for the
\code{derive_param_tte()} function.

\strong{Note:} This is a wrapper function for the more generic \code{tte_source()}.
}
\examples{
# Last study date known alive censor

censor_source(
  dataset_name = "adsl",
  date = LSTALVDT,
  set_values_to = exprs(
    EVNTDESC = "ALIVE",
    SRCDOM = "ADSL",
    SRCVAR = "LSTALVDT"
  )
)
}
\seealso{
\code{\link[=derive_param_tte]{derive_param_tte()}}, \code{\link[=event_source]{event_source()}}

Source Objects: 
\code{\link{basket_select}()},
\code{\link{death_event}},
\code{\link{event}()},
\code{\link{event_joined}()},
\code{\link{event_source}()},
\code{\link{flag_event}()},
\code{\link{query}()},
\code{\link{records_source}()},
\code{\link{tte_source}()}
}
\concept{source_specifications}
\keyword{source_specifications}
