% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/user_utils.R
\name{print_named_list}
\alias{print_named_list}
\title{Print Named List}
\usage{
print_named_list(list, indent = 0)
}
\arguments{
\item{list}{A named list

\describe{
\item{Default value}{none}
}}

\item{indent}{Indent

The output is indented by the specified number of characters.

\describe{
\item{Default value}{\code{0}}
}}
}
\value{
No return value, called for side effects
}
\description{
Print Named List
}
\examples{
print_named_list(death_event)
}
\seealso{
Utilities for printing: 
\code{\link{print.adam_templates}()},
\code{\link{print.duplicates}()},
\code{\link{print.source}()}
}
\concept{utils_print}
\keyword{internal}
