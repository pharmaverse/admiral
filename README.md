# admiral <img src="man/figures/logo.png" align="right" width="200" style="margin-left:50px;"/>

<!-- badges: start -->
[![Codecov test coverage](https://codecov.io/gh/Roche-GSK/admiral/branch/master/graph/badge.svg)](https://codecov.io/gh/Roche-GSK/admiral?branch=master)
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
remotes::install_github("Roche-GSK/admiral")
```

## Scope

* Build a toolbox of re-usable functions and utilities to create ADaM datasets using R scripts in a modular manner
* Pharmaceutical communities and companies are encouraged to contribute to admiral following the provided programming strategy and modular approach
* All functions are documented, tested, include examples and are listed [in the Reference section](https://roche-gsk.github.io/admiral/reference/index.html) 
* Vignettes how to create ADSL, BDS and OCCUR datasets, including example scripts
* Vignettes for ADaM dataset specific functionality (i.e. dictionary coding, date imputation, SMQs ...)

## Usage

* Think of admiral as a toolbox of modular blocks (R functions) to create analysis derivations
    * Each block has a stand alone purpose, i.e. each function provides a specific functionality
    * Data Scientists can create their own blocks, i.e. create own R functions

* Constructing ADaM dataset should become like building out of blocks that are based on admiral modular functions and user created modular functions


## References and Documentation

* Please go to [Get Started](https://roche-gsk.github.io/admiral/articles/admiral.html) section to start using admiral
* Please see the [Programming Strategy](https://roche-gsk.github.io/admiral/articles/programming_strategy.html) to understand how functions are created
* Please see the [FAQ](https://roche-gsk.github.io/admiral/articles/faq.html) for the most frequent questions

## Contact 

In case your questions are not answered please contact us on [Slack](https://app.slack.com/client/T028PB489D3/C028SJ83KM1)
