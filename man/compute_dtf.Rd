% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/derive_vars_dt.R
\name{compute_dtf}
\alias{compute_dtf}
\title{Derive the Date Imputation Flag}
\usage{
compute_dtf(dtc, dt)
}
\arguments{
\item{dtc}{The date character vector (\code{--DTC}).

A character date is expected in a format like \code{yyyy-mm-ddThh:mm:ss} (partial or complete).

\describe{
\item{Default value}{none}
}}

\item{dt}{The  Date vector to compare.

A date object is expected.

\describe{
\item{Default value}{none}
}}
}
\value{
The date imputation flag (\verb{*DTF}) (character value of \code{"D"}, \code{"M"} , \code{"Y"} or \code{NA})
}
\description{
Derive the date imputation flag (\verb{*DTF}) comparing a date character vector
(\code{--DTC}) with a Date vector (\verb{*DT}).
}
\details{
Usually this computation function can not be used with \verb{\%>\%}.
}
\examples{
compute_dtf(dtc = "2019-07", dt = as.Date("2019-07-18"))
compute_dtf(dtc = "2019", dt = as.Date("2019-07-18"))
compute_dtf(dtc = "--06-01T00:00", dt = as.Date("2022-06-01"))
compute_dtf(dtc = "2022-06--T00:00", dt = as.Date("2022-06-01"))
compute_dtf(dtc = "2022---01T00:00", dt = as.Date("2022-06-01"))
compute_dtf(dtc = "2022----T00:00", dt = as.Date("2022-06-01"))
}
\seealso{
Date/Time Computation Functions that returns a vector: 
\code{\link{compute_age_years}()},
\code{\link{compute_duration}()},
\code{\link{compute_tmf}()},
\code{\link{convert_date_to_dtm}()},
\code{\link{convert_dtc_to_dt}()},
\code{\link{convert_dtc_to_dtm}()},
\code{\link{impute_dtc_dt}()},
\code{\link{impute_dtc_dtm}()}
}
\concept{com_date_time}
\keyword{com_date_time}
