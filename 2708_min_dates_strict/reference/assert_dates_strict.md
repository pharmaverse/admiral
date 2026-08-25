# Assert Strict Date Ranges

Asserts that the minimum dates specified in `min_dates_strict` are not
greater than the maximum dates specified in `max_dates_strict`.

## Usage

``` r
assert_dates_strict(min_dates_strict, max_dates_strict)
```

## Arguments

- min_dates_strict:

  A list of minimum dates to check.

  Default value

  :   none

- max_dates_strict:

  A list of maximum dates to check.

  Default value

  :   none

## Value

Invisibly returns `NULL` if the assertion passes.

## Details

If any minimum date is greater than its corresponding maximum date, an
error is thrown with the invalid combinations of minimum and maximum
dates.
