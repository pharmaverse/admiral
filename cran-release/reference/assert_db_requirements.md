# Check required parameters for a basket

If a basket (SMQ, SDG, ....) are requested, the version and a function
to access the database must be provided. The function checks these
requirements.

## Usage

``` r
assert_db_requirements(
  version,
  version_arg_name,
  fun,
  fun_arg_name,
  queries,
  i
)
```

## Arguments

- version:

  Version provided by user

  Default value

  :   none

- version_arg_name:

  Name of the argument providing the version

  Default value

  :   none

- fun:

  Function provided by user

  Default value

  :   none

- fun_arg_name:

  Name of the argument providing the function

  Default value

  :   none

- queries:

  Queries provide by user

  Default value

  :   none

- i:

  Index of query being checked

  Default value

  :   none

## Value

An error is issued if `version` or `fun` is null.
