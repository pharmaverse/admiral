% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/compute_qual_imputation.R
\name{compute_qual_imputation}
\alias{compute_qual_imputation}
\title{Function to Impute Values When Qualifier Exists in Character Result}
\usage{
compute_qual_imputation(character_value, imputation_type = 1, factor = 0)
}
\arguments{
\item{character_value}{Character version of value to be imputed

\describe{
\item{Default value}{none}
}}

\item{imputation_type}{(default value=1)
Valid Values:
1: Strip <, >, = and convert to numeric.
2: imputation_type=1 and if the character value contains a < or >, the number of
of decimals associated with the character value is found and then a factor of
1/10^(number of decimals + 1) will be added/subtracted from the numeric value.
If no decimals exists, a factor of 1/10 will be added/subtracted from the value.

\describe{
\item{Default value}{\code{1}}
}}

\item{factor}{Numeric value (default=0), when using \code{imputation_type} = 1, this
value can be added or subtracted when the qualifier is removed.

\describe{
\item{Default value}{\code{0}}
}}
}
\value{
The imputed value
}
\description{
Derive an imputed value
}
\examples{
compute_qual_imputation("<40")
}
\seealso{
BDS-Findings Functions that returns a vector: 
\code{\link{compute_bmi}()},
\code{\link{compute_bsa}()},
\code{\link{compute_egfr}()},
\code{\link{compute_framingham}()},
\code{\link{compute_map}()},
\code{\link{compute_qtc}()},
\code{\link{compute_qual_imputation_dec}()},
\code{\link{compute_rr}()},
\code{\link{compute_scale}()},
\code{\link{transform_range}()}
}
\concept{com_bds_findings}
\keyword{com_bds_findings}
