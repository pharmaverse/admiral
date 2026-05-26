# Exploratory mutation testing with `muttest`

This folder keeps mutation-testing setup out of the package code paths.

## Scope

The script currently targets simpler `admiral` functions:

- `compute_age_years()` (`R/compute_age_years.R`)
- `compute_scale()` (`R/compute_scale.R`)
- `filter_exist()` / `filter_not_exist()` (`R/filter_exist.R`)

Matching tests are expected in:

- `tests/testthat/test-compute_age_years.R`
- `tests/testthat/test-compute_scale.R`
- `tests/testthat/test-filter_exist.R`

## Run

From the repository root:

```r
install.packages("muttest")
source("dev/muttest/run_mutation_checks.R")
```

The script creates a mutation plan for the selected source files and runs `muttest`
with `FileTestStrategy`, so only matching test files are executed for each mutant.
