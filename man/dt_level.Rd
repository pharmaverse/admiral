% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/derive_vars_dt_dtm_utils.R
\name{dt_level}
\alias{dt_level}
\title{Create a \code{dt_level} object}
\usage{
dt_level(level)
}
\arguments{
\item{level}{Date level

\describe{
\item{Permitted values}{\code{"Y"} (year, highest level), \code{"M"} (month), \code{"D"}
(day), \code{"n"} (none, lowest level)}
\item{Default value}{none}
}}
}
\value{
A \code{dt_level} object
}
\description{
Create a \code{dt_level} object
}
\details{
A \code{dt_level} object is an ordered factor, i.e., two objects can be
compared.
}
\seealso{
Utilities used for date imputation: 
\code{\link{dtm_level}()},
\code{\link{get_imputation_target_date}()},
\code{\link{get_imputation_target_time}()},
\code{\link{get_partialdatetime}()},
\code{\link{restrict_imputed_dtc_dt}()},
\code{\link{restrict_imputed_dtc_dtm}()}
}
\concept{utils_impute}
\keyword{internal}
