<!-- Please do not edit the README.md file as it is auto-generated. Only edit the README.Rmd file -->

# admiraldev <img src="man/figures/logo.png" align="right" alt="" width="120" />

ADaM in R Asset Library Development Utilities

<!-- badges: start -->

<!-- badges: end -->

## Purpose

Tools for developing functions and maintaining a healthy code base
within the family of admiral R packages. `{admiraldev}` is intended to
be used when developing `{admiral}` or `{admiral}` extension packages.

__NOTE:__ This package is not intended for standalone use but rather as
a central dependency for all developer utilities of `{admiral}` and its
extension packages

## Installation

The package is available from CRAN and can be installed by running
`install.packages("admiraldev")`.

To install the latest development version of the package directly from
GitHub use the following code:

```
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}

remotes::install_github("pharmaverse/admiraldev", ref = "devel")
```

## Release Schedule

`{admiraldev}` is to be officially released to CRAN one week before an
official release of `{admiral}`. You can find the release schedule for
`{admiral}` packages
[here](https://pharmaverse.github.io/admiral/#release-schedule).
