% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/derive_vars_query.R
\name{derive_vars_query}
\alias{derive_vars_query}
\title{Derive Query Variables}
\usage{
derive_vars_query(dataset, dataset_queries)
}
\arguments{
\item{dataset}{Input dataset

\describe{
\item{Default value}{none}
}}

\item{dataset_queries}{A dataset containing required columns \code{PREFIX},
\code{GRPNAME}, \code{SRCVAR}, \code{TERMCHAR} and/or \code{TERMNUM}, and optional columns
\code{GRPID}, \code{SCOPE}, \code{SCOPEN}.

\code{create_query_data()} can be used to create the dataset.

\describe{
\item{Default value}{none}
}}
}
\value{
The input dataset with query variables derived.
}
\description{
Derive Query Variables
}
\details{
This function can be used to derive CDISC variables such as
\code{SMQzzNAM}, \code{SMQzzCD}, \code{SMQzzSC}, \code{SMQzzSCN}, and \code{CQzzNAM} in ADAE and
ADMH, and variables such as \code{SDGzzNAM}, \code{SDGzzCD}, and \code{SDGzzSC} in ADCM.
An example usage of this function can be found in the \code{vignette("occds")}.

A query dataset is expected as an input to this function. See the
\code{vignette("queries_dataset")} for descriptions, or call \code{data("queries")}
for an example of a query dataset.

For each unique element in \code{PREFIX}, the corresponding "NAM"
variable will be created. For each unique \code{PREFIX}, if \code{GRPID} is
not "" or NA, then the corresponding "CD" variable is created; similarly,
if \code{SCOPE} is not "" or NA, then the corresponding "SC" variable will
be created; if \code{SCOPEN} is not "" or NA, then the corresponding
"SCN" variable will be created.

For each record in \code{dataset}, the "NAM" variable takes the value of
\code{GRPNAME} if the value of \code{TERMCHAR} or \code{TERMNUM} in \code{dataset_queries} matches
the value of the respective SRCVAR in \code{dataset}.
Note that \code{TERMCHAR} in \code{dataset_queries} dataset may be NA only when \code{TERMNUM}
is non-NA and vice versa. The matching is case insensitive.
The "CD", "SC", and "SCN" variables are derived accordingly based on
\code{GRPID}, \code{SCOPE}, and \code{SCOPEN} respectively,
whenever not missing.
}
\examples{
library(tibble)
data("queries")
adae <- tribble(
  ~USUBJID, ~ASTDTM, ~AETERM, ~AESEQ, ~AEDECOD, ~AELLT, ~AELLTCD,
  "01", "2020-06-02 23:59:59", "ALANINE AMINOTRANSFERASE ABNORMAL",
  3, "Alanine aminotransferase abnormal", NA_character_, NA_integer_,
  "02", "2020-06-05 23:59:59", "BASEDOW'S DISEASE",
  5, "Basedow's disease", NA_character_, 1L,
  "03", "2020-06-07 23:59:59", "SOME TERM",
  2, "Some query", "Some term", NA_integer_,
  "05", "2020-06-09 23:59:59", "ALVEOLAR PROTEINOSIS",
  7, "Alveolar proteinosis", NA_character_, NA_integer_
)
derive_vars_query(adae, queries)
}
\seealso{
\code{\link[=create_query_data]{create_query_data()}}

OCCDS Functions: 
\code{\link{derive_var_trtemfl}()},
\code{\link{derive_vars_atc}()}
}
\concept{der_occds}
\keyword{der_occds}
