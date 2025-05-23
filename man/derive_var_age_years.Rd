% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/derive_vars_aage.R
\name{derive_var_age_years}
\alias{derive_var_age_years}
\title{Derive Age in Years}
\usage{
derive_var_age_years(dataset, age_var, age_unit = NULL, new_var)
}
\arguments{
\item{dataset}{Input dataset

The variables specified by the \code{age_var} argument are expected to be in the dataset.

\describe{
\item{Default value}{none}
}}

\item{age_var}{Age variable.

A numeric object is expected.

\describe{
\item{Default value}{none}
}}

\item{age_unit}{Age unit.

The \code{age_unit} argument is only expected when there is NOT a variable \code{age_var+U}
in \code{dataset}. This gives the unit of the \code{age_var} variable and is used to convert
AGE to 'years' so that grouping can occur.

\describe{
\item{Permitted values}{'years', 'months', 'weeks', 'days', 'hours', 'minutes', 'seconds'}
\item{Default value}{\code{NULL}}
}}

\item{new_var}{New age variable to be created in years. The returned values are
doubles and NOT integers.
'

\describe{
\item{Default value}{none}
}}
}
\value{
The input dataset (\code{dataset}) with \code{new_var} variable added in years.
}
\description{
Converts the given age variable (\code{age_var}) to the unit 'years' from the current
units given in the \code{age_var+U} variable or \code{age_unit} argument and stores
in a new variable (\code{new_var}).
}
\details{
This function is used to convert an age variable into the unit 'years'
which can then be used to create age groups. The resulting column contains the
equivalent years as a double. Note, underlying computations assume an equal number
of days in each year (365.25).
}
\examples{
library(tibble)

# Derive age with age units specified
data <- tribble(
  ~AGE, ~AGEU,
  27, "days",
  24, "months",
  3, "years",
  4, "weeks",
  1, "years"
)

derive_var_age_years(data, AGE, new_var = AAGE)

# Derive age without age units variable specified
data <- tribble(
  ~AGE,
  12,
  24,
  36,
  48
)
derive_var_age_years(data, AGE, age_unit = "months", new_var = AAGE)
}
\seealso{
\code{\link[=derive_vars_duration]{derive_vars_duration()}}

ADSL Functions that returns variable appended to dataset: 
\code{\link{derive_vars_aage}()},
\code{\link{derive_vars_extreme_event}()},
\code{\link{derive_vars_period}()}
}
\concept{der_adsl}
\keyword{der_adsl}
