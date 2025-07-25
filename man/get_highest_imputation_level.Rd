% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/derive_vars_dt_dtm_utils.R
\name{get_highest_imputation_level}
\alias{get_highest_imputation_level}
\title{Get Highest Imputation Level}
\usage{
get_highest_imputation_level(highest_imputation, create_datetime)
}
\arguments{
\item{highest_imputation}{A character indicating the highest imputation level.

\describe{
\item{Default value}{none}
}}

\item{create_datetime}{A logical indicating whether datetime factors levels are required.

\describe{
\item{Default value}{none}
}}
}
\value{
A \code{dt_level()} or \code{dtm_level()} object representing the highest imputation level.
}
\description{
Returns the \code{dt_level()} or \code{dtm_level()} representation of the \code{highest_imputation}
character value. The level object allows comparisons of levels.
}
\examples{
# Get highest imputation level for date
highest_level_date <- admiral:::get_highest_imputation_level(
  highest_imputation = "Y",
  create_datetime = FALSE
)
print(highest_level_date)

# Get highest imputation level for datetime
highest_level_datetime <- admiral:::get_highest_imputation_level(
  highest_imputation = "Y",
  create_datetime = TRUE
)
print(highest_level_datetime)

# Get highest imputation level for date with month level
highest_level_month_date <- admiral:::get_highest_imputation_level(
  highest_imputation = "M",
  create_datetime = FALSE
)
print(highest_level_month_date)

# Get highest imputation level for datetime with hour level
highest_level_hour_datetime <- admiral:::get_highest_imputation_level(
  highest_imputation = "h",
  create_datetime = TRUE
)
print(highest_level_hour_datetime)

}
\keyword{internal}
