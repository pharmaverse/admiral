# admiral <img src="man/figures/logo.png" align="right" width="200" style="margin-left:50px;"/>

<!-- badges: start -->
[![CRAN status](https://www.r-pkg.org/badges/version/admiral)](https://CRAN.R-project.org/package=admiral)
[![R-CMD-check](https://github.com/pharmaverse/admiral/workflows/R-CMD-check/badge.svg)](https://github.com/pharmaverse/admiral/actions)
[![Codecov test coverage](https://codecov.io/gh/Roche-GSK/admiral/branch/main/graph/badge.svg)](https://codecov.io/gh/Roche-GSK/admiral?branch=main)
<!-- badges: end -->

ADaM in R Asset Library

## Purpose

To provide an open source, modularized toolbox that enables the pharmaceutical programming community to develop ADaM datasets in R.

## Installation

Once the package is available from CRAN you'll be able to install it using `install.packages("admiral")`.

In the meantime you can install the latest release of the package directly from GitHub.

```r
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}
remotes::install_github("Roche-GSK/admiral.test") # This is a required dependency of {admiral}
remotes::install_github("Roche-GSK/admiral")
```

## Scope

* A toolbox of re-usable functions and utilities to create ADaM datasets using R scripts in a modular manner (an "opinionated" design strategy)
* Pharmaceutical communities and companies are encouraged to contribute to admiral following the provided programming strategy and modular approach
* All functions are documented, tested, include examples and are listed in the
[Reference section](https://roche-gsk.github.io/admiral/reference/index.html)
* Vignettes how to create ADSL, BDS and OCCUR datasets, including example scripts
* Vignettes for ADaM dataset specific functionality (i.e. dictionary coding, date imputation, SMQs ...)

## Usage

* Think of admiral as a toolbox of modular blocks (R functions) to create analysis derivations
    * Each block has a stand alone purpose, i.e. each function provides a specific functionality
    * Data Scientists can create their own blocks, i.e. create own R functions
* Constructing ADaM dataset should become like building out of blocks that are based on admiral modular functions and user created modular functions

## Expectations

* admiral will never cover 100% of eventualities that could be needed to produce ADaMs across each and every company/disease area/study - ADaM is infinite
* Some flexibility can be added to the functions, but only where there is an agreed common need across the industry as this has to balance vs ease of usage and testing
* One of our principle design decisions in creating admiral was to prioritise transparency and simplicity for our users - and not to let this ever become a "black-box" toolkit
* We hope admiral offers a chance for users to be programmers - this is not a "run 1 line and an ADaM appears" solution or an attempt to automate ADaM
* It is expected for companies to adopt admiral that a company-specific extension package would likely be needed (e.g. `{admiral.roche}` or `{admiral.gsk}`)
    * Consider this an opportunity to influence company-specific standards, as admiral offers a bridge towards a more industry-aligned implementation of ADaM
* From the core `{admiral}` package covering common functions and mostly safety templates, there will be further package extensions dedicated to certain disease area endpoints (e.g. `{admiral.onco}` or `{admiral.hiv}`)

## References and Documentation

* Please go to [Get Started](https://roche-gsk.github.io/admiral/articles/admiral.html) section to start using admiral
* Please see the [Programming Strategy](https://roche-gsk.github.io/admiral/articles/programming_strategy.html) to understand how functions are created
* Please see the [FAQ](https://roche-gsk.github.io/admiral/articles/faq.html) for the most frequent questions
* Please see the [Contribution Model](https://roche-gsk.github.io/admiral/articles/contribution_model.html) for how to get involved with making contributions

## Contact 

We use the following for support and communications between user and developer community:
* [Slack](https://app.slack.com/client/T028PB489D3/C02M8KN8269) - for informal discussions, Q&A and building our user community. If you don't have access, use this [link](https://join.slack.com/t/pharmaverse/shared_invite/zt-yv5atkr4-Np2ytJ6W_QKz_4Olo7Jo9A) to join the pharmaverse Slack workspace
* [GitHub Issues](https://github.com/Roche-GSK/admiral/issues) - for direct feedback, enhancement requests or raising bugs
