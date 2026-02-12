# Get Terms from the Queries Database

The function checks if all requirements to access the database are
fulfilled (version and access function are available, see
[`assert_db_requirements()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/assert_db_requirements.md)),
reads the terms from the database, and checks if the dataset with the
terms is in the expected format (see
[`assert_terms()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/assert_terms.md)).

## Usage

``` r
get_terms_from_db(
  version,
  fun,
  fun_name,
  queries,
  definition,
  expect_grpname = FALSE,
  expect_grpid = FALSE,
  i,
  temp_env
)
```

## Arguments

- version:

  Version

  The version must be non null. Otherwise, an error is issued. The value
  is passed to the access function (`fun`).

  Default value

  :   none

- fun:

  Access function

  The access function must be non null. Otherwise, an error is issued.
  The function is called to retrieve the terms.

  Default value

  :   none

- fun_name:

  Name of access function

  The character name of the access function, usually created with
  `deparse(substitute(fun))`. This must be non null. Otherwise, an error
  is issued.

  Default value

  :   none

- queries:

  Queries

  List of all queries passed to
  [`create_query_data()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/create_query_data.md).
  It is used for error messages.

  Default value

  :   none

- definition:

  Definition of the query

  The definition is passed to the access function. It defines which
  terms are returned.

  Default value

  :   none

- expect_grpname:

  Is `GRPNAME` expected in the output dataset?

  Default value

  :   `FALSE`

- expect_grpid:

  Is `GRPID` expected in the output dataset?

  Default value

  :   `FALSE`

- i:

  Index of `definition` in `queries`

  The value is used for error messages.

  Default value

  :   none

- temp_env:

  Temporary environment

  The value is passed to the access function.

  Default value

  :   none

## Value

Output dataset of the access function
