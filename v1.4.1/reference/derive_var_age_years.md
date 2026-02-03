# Derive Age in Years

Converts the given age variable (`age_var`) to the unit 'years' from the
current units given in the `age_var+U` variable or `age_unit` argument
and stores in a new variable (`new_var`).

## Usage

``` r
derive_var_age_years(dataset, age_var, age_unit = NULL, new_var)
```

## Arguments

- dataset:

  Input dataset

  The variables specified by the `age_var` argument are expected to be
  in the dataset.

  Default value

  :   none

- age_var:

  Age variable.

  A numeric object is expected.

  Default value

  :   none

- age_unit:

  Age unit.

  The `age_unit` argument is only expected when there is NOT a variable
  `age_var+U` in `dataset`. This gives the unit of the `age_var`
  variable and is used to convert AGE to 'years' so that grouping can
  occur.

  Permitted values

  :   'years', 'months', 'weeks', 'days', 'hours', 'minutes', 'seconds'

  Default value

  :   `NULL`

- new_var:

  New age variable to be created in years. The returned values are
  doubles and NOT integers. '

  Default value

  :   none

## Value

The input dataset (`dataset`) with `new_var` variable added in years.

## Details

This function is used to convert an age variable into the unit 'years'
which can then be used to create age groups. The resulting column
contains the equivalent years as a double. Note, underlying computations
assume an equal number of days in each year (365.25).

## See also

[`derive_vars_duration()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/derive_vars_duration.md)

ADSL Functions that returns variable appended to dataset:
[`derive_vars_aage()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/derive_vars_aage.md),
[`derive_vars_extreme_event()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/derive_vars_extreme_event.md),
[`derive_vars_period()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/derive_vars_period.md)

## Examples

``` r
library(tibble)

# Derive age with age units specified
data <- tribble(
  ~AGE, ~AGEU,
  27, "days",
  24, "months",
  3, "years",
  4, "weeks",
  1, "years"
)

derive_var_age_years(data, AGE, new_var = AAGE)
#> # A tibble: 5 × 3
#>     AGE AGEU     AAGE
#>   <dbl> <chr>   <dbl>
#> 1    27 days   0.0739
#> 2    24 months 2     
#> 3     3 years  3     
#> 4     4 weeks  0.0767
#> 5     1 years  1     

# Derive age without age units variable specified
data <- tribble(
  ~AGE,
  12,
  24,
  36,
  48
)
derive_var_age_years(data, AGE, age_unit = "months", new_var = AAGE)
#> # A tibble: 4 × 2
#>     AGE  AAGE
#>   <dbl> <dbl>
#> 1    12     1
#> 2    24     2
#> 3    36     3
#> 4    48     4
```
