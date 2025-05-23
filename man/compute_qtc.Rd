% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/derive_param_qtc.R
\name{compute_qtc}
\alias{compute_qtc}
\title{Compute Corrected QT}
\usage{
compute_qtc(qt, rr, method)
}
\arguments{
\item{qt}{QT interval

A numeric vector is expected. It is expected that QT is measured in ms or
msec.

\describe{
\item{Default value}{none}
}}

\item{rr}{RR interval

A numeric vector is expected. It is expected that RR is measured in ms or
msec.

\describe{
\item{Default value}{none}
}}

\item{method}{Method used to QT correction

\describe{
\item{Permitted values}{\code{"Bazett"}, \code{"Fridericia"}, \code{"Sagie"}}
\item{Default value}{none}
}}
}
\value{
QT interval in ms
}
\description{
Computes corrected QT using Bazett's, Fridericia's or Sagie's formula.
}
\details{
Depending on the chosen \code{method} one of the following formulae is used.

\emph{Bazett}: \deqn{\frac{QT}{\sqrt{\frac{RR}{1000}}}}{QT/\sqrt(RR/1000)}

\emph{Fridericia}: \deqn{\frac{QT}{\sqrt[3]{\frac{RR}{1000}}}}{QT/(RR/1000)^(1/3)}

\emph{Sagie}: \deqn{1000\left(\frac{QT}{1000} + 0.154\left(1 - \frac{RR}{1000}\right)\right)}{
1000(QT/1000 + 0.154(1 - RR/1000))}

Usually this computation function can not be used with \verb{\%>\%}.
}
\examples{
compute_qtc(qt = 350, rr = 857, method = "Bazett")

compute_qtc(qt = 350, rr = 857, method = "Fridericia")

compute_qtc(qt = 350, rr = 857, method = "Sagie")
}
\seealso{
\code{\link[=derive_param_qtc]{derive_param_qtc()}}

BDS-Findings Functions that returns a vector: 
\code{\link{compute_bmi}()},
\code{\link{compute_bsa}()},
\code{\link{compute_egfr}()},
\code{\link{compute_framingham}()},
\code{\link{compute_map}()},
\code{\link{compute_qual_imputation}()},
\code{\link{compute_qual_imputation_dec}()},
\code{\link{compute_rr}()},
\code{\link{compute_scale}()},
\code{\link{transform_range}()}
}
\concept{com_bds_findings}
\keyword{com_bds_findings}
