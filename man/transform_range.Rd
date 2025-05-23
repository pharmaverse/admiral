% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/transform_range.R
\name{transform_range}
\alias{transform_range}
\title{Transform Range}
\usage{
transform_range(
  source,
  source_range,
  target_range,
  flip_direction = FALSE,
  outside_range = "NA"
)
}
\arguments{
\item{source}{A vector of values to be transformed

A numeric vector is expected.

\describe{
\item{Default value}{none}
}}

\item{source_range}{The permitted source range

A numeric vector containing two elements is expected, representing the
lower and upper bounds of the permitted source range.

\describe{
\item{Default value}{none}
}}

\item{target_range}{The target range

A numeric vector containing two elements is expected, representing the
lower and upper bounds of the target range.

\describe{
\item{Default value}{none}
}}

\item{flip_direction}{Flip direction of the range?

The transformed values will be reversed within the target range, e.g.
within the range 0 to 100, 25 would be reversed to 75.

\describe{
\item{Permitted values}{\code{TRUE}, \code{FALSE}}
\item{Default value}{\code{FALSE}}
}}

\item{outside_range}{Handling of values outside the source range

Values outside the source range (\code{source_range}) are transformed to \code{NA}.

If \code{"warning"} or \code{"error"} is specified, a warning or error is issued if
\code{source} includes any values outside the source range.

\describe{
\item{Permitted values}{\code{"NA"}, \code{"warning"}, \code{"error"}}
\item{Default value}{\code{"NA"}}
}}
}
\value{
The source linearly transformed to the target range
}
\description{
Transforms results from the source range to the target range. For example,
for transforming source values 1, 2, 3, 4, 5 to 0, 25, 50, 75, 100.
}
\details{
Returns the values of \code{source} linearly transformed from the source
range (\code{source_range}) to the target range (\code{target_range}). Values outside
the source range are set to \code{NA}.
}
\examples{
transform_range(
  source = c(1, 4, 3, 6, 5),
  source_range = c(1, 5),
  target_range = c(0, 100)
)

transform_range(
  source = c(1, 4, 3, 6, 5),
  source_range = c(1, 5),
  target_range = c(0, 100),
  flip_direction = TRUE
)
}
\seealso{
BDS-Findings Functions that returns a vector: 
\code{\link{compute_bmi}()},
\code{\link{compute_bsa}()},
\code{\link{compute_egfr}()},
\code{\link{compute_framingham}()},
\code{\link{compute_map}()},
\code{\link{compute_qtc}()},
\code{\link{compute_qual_imputation}()},
\code{\link{compute_qual_imputation_dec}()},
\code{\link{compute_rr}()},
\code{\link{compute_scale}()}
}
\concept{com_bds_findings}
\keyword{com_bds_findings}
