# Asserts Requirements for Terms for Queries

The function checks the requirements for terms for queries provided by
the user. The terms could have been provided directly in the query
definition or via a user provided function for accessing a SMQ or SDG
database.

## Usage

``` r
assert_terms(terms, expect_grpname = FALSE, expect_grpid = FALSE, source_text)
```

## Arguments

- terms:

  Terms provided by user

  Default value

  :   none

- expect_grpname:

  Is the `GRPNAME` column expected?

  Default value

  :   `FALSE`

- expect_grpid:

  Is the `GRPID` column expected?

  Default value

  :   `FALSE`

- source_text:

  Text describing the source of the terms, e.g.,
  `"the data frame provided for the `definition` element"`.

  Default value

  :   none

## Value

An error is issued if

- `terms` is not a data frame,

- `terms` has zero observations,

- the `SRCVAR` variable is not in `terms`,

- neither the `TERMCHAR` nor the `TERMNUM` variable is in `terms`,

- `expect_grpname == TRUE` and the `GRPNAME` variable is not in `terms`,

- `expect_grpid == TRUE` and the `GRPID` variable is not in `terms`,

## Examples

``` r
try(
  assert_terms(
    terms = 42,
    source_text = "object provided by the `definition` element"
  )
)
#> Error in assert_terms(terms = 42, source_text = "object provided by the `definition` element") : 
#>   could not find function "assert_terms"
```
