% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/slice_derivation.R
\name{slice_derivation}
\alias{slice_derivation}
\title{Execute a Derivation with Different Arguments for Subsets of the Input Dataset}
\usage{
slice_derivation(dataset, derivation, args = NULL, ...)
}
\arguments{
\item{dataset}{Input dataset

\describe{
\item{Default value}{none}
}}

\item{derivation}{Derivation

A function that performs a specific derivation is expected. A derivation
adds variables or observations to a dataset. The first argument of a
derivation must expect a dataset and the derivation must return a dataset.
All expected arguments for the derivation function must be provided through
the \code{params()} object passed to the \code{args} argument or be provided in \emph{every}
\code{derivation_slice()}.

\describe{
\item{Default value}{none}
}}

\item{args}{Arguments of the derivation

A \code{params()} object is expected.

\describe{
\item{Default value}{\code{NULL}}
}}

\item{...}{A \code{derivation_slice()} object is expected

Each slice defines a subset of the input dataset and some of the parameters
for the derivation. The derivation is called on the subset with the
parameters specified by the \code{args} parameter and the \code{args} field of the
\code{derivation_slice()} object. If a parameter is specified for both, the
value in \code{derivation_slice()} overwrites the one in \code{args}.

\describe{
\item{Default value}{none}
}}
}
\value{
The input dataset with the variables derived by the derivation added
}
\description{
The input dataset is split into slices (subsets) and for each slice the
derivation is called separately. Some or all arguments of the derivation
may vary depending on the slice.
}
\details{
For each slice the derivation is called on the subset defined by the
\code{filter} field of the \code{derivation_slice()} object and with the parameters
specified by the \code{args} parameter and the \code{args} field of the
\code{derivation_slice()} object. If a parameter is specified for both, the
value in \code{derivation_slice()} overwrites the one in \code{args}.
\itemize{
\item Observations that match with more than one slice are only considered for
the first matching slice.
\item The derivation is called for slices with no observations.
\item Observations with no match to any of the slices are included in the
output dataset but the derivation is not called for them.
}

It is also possible to pass functions from outside the \code{{admiral}} package
to \code{slice_derivation()}, e.g. an extension package function, or
\code{dplyr::mutate()}. The only requirement for a function being passed to \code{derivation} is that
it must take a dataset as its first argument and return a dataset.
}
\examples{
library(tibble)
library(stringr)
advs <- tribble(
  ~USUBJID, ~VSDTC,       ~VSTPT,
  "1",      "2020-04-16", NA_character_,
  "1",      "2020-04-16", "BEFORE TREATMENT"
)

# For the second slice filter is set to TRUE. Thus derive_vars_dtm is called
# with time_imputation = "last" for all observations which do not match for the
# first slice.
slice_derivation(
  advs,
  derivation = derive_vars_dtm,
  args = params(
    dtc = VSDTC,
    new_vars_prefix = "A"
  ),
  derivation_slice(
    filter = str_detect(VSTPT, "PRE|BEFORE"),
    args = params(time_imputation = "first")
  ),
  derivation_slice(
    filter = TRUE,
    args = params(time_imputation = "last")
  )
)

}
\seealso{
\code{\link[=params]{params()}} \code{\link[=restrict_derivation]{restrict_derivation()}} \code{\link[=call_derivation]{call_derivation()}}

Higher Order Functions: 
\code{\link{call_derivation}()},
\code{\link{derivation_slice}()},
\code{\link{restrict_derivation}()}
}
\concept{high_order_function}
\keyword{high_order_function}
