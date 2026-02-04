# Create an `query` object

A `query` object defines a query, e.g., a Standard MedDRA Query (SMQ), a
Standardized Drug Grouping (SDG), or a customized query (CQ). It is used
as input to
[`create_query_data()`](https:/pharmaverse.github.io/admiral/main/reference/create_query_data.md).

## Usage

``` r
query(prefix, name = auto, id = NULL, add_scope_num = FALSE, definition = NULL)
```

## Arguments

- prefix:

  The value is used to populate `PREFIX` in the output dataset of
  [`create_query_data()`](https:/pharmaverse.github.io/admiral/main/reference/create_query_data.md),
  e.g., `"SMQ03"`

  Default value

  :   none

- name:

  The value is used to populate `GRPNAME` in the output dataset of
  [`create_query_data()`](https:/pharmaverse.github.io/admiral/main/reference/create_query_data.md).
  If the `auto` keyword is specified, the variable is set to the name of
  the query in the SMQ/SDG database.

  Permitted values

  :   A character scalar or the `auto` keyword. The `auto` keyword is
      permitted only for queries which are defined by an
      [`basket_select()`](https:/pharmaverse.github.io/admiral/main/reference/basket_select.md)
      object.

  Default value

  :   `auto`

- id:

  The value is used to populate `GRPID` in the output dataset of
  [`create_query_data()`](https:/pharmaverse.github.io/admiral/main/reference/create_query_data.md).
  If the `auto` keyword is specified, the variable is set to the id of
  the query in the SMQ/SDG database.

  Permitted values

  :   A integer scalar or the `auto` keyword. The `auto` keyword is
      permitted only for queries which are defined by an
      [`basket_select()`](https:/pharmaverse.github.io/admiral/main/reference/basket_select.md)
      object.

  Default value

  :   `NULL`

- add_scope_num:

  Determines if `SCOPEN` in the output dataset of
  [`create_query_data()`](https:/pharmaverse.github.io/admiral/main/reference/create_query_data.md)
  is populated

  If the parameter is set to `TRUE`, the definition must be an
  [`basket_select()`](https:/pharmaverse.github.io/admiral/main/reference/basket_select.md)
  object.

  *Default*: `FALSE`

  Permitted values

  :   `TRUE`, `FALSE`

  Default value

  :   `FALSE`

- definition:

  Definition of terms belonging to the query

  There are three different ways to define the terms:

  - An
    [`basket_select()`](https:/pharmaverse.github.io/admiral/main/reference/basket_select.md)
    object is specified to select a query from the SMQ database.

  - A data frame with columns `SRCVAR` and `TERMCHAR` or `TERMNUM` can
    be specified to define the terms of a customized query. The `SRCVAR`
    should be set to the name of the variable which should be used to
    select the terms, e.g., `"AEDECOD"` or `"AELLTCD"`. `SRCVAR` does
    not need to be constant within a query. For example a query can be
    based on `AEDECOD` and `AELLT`.

    If `SRCVAR` refers to a character variable, `TERMCHAR` should be set
    to the value the variable. If it refers to a numeric variable,
    `TERMNUM` should be set to the value of the variable. If only
    character variables or only numeric variables are used, `TERMNUM` or
    `TERMCHAR` respectively can be omitted.

  - A list of data frames and
    [`basket_select()`](https:/pharmaverse.github.io/admiral/main/reference/basket_select.md)
    objects can be specified to define a customized query based on
    custom terms and SMQs. The data frames must have the same structure
    as described for the previous item.

  Permitted values

  :   an
      [`basket_select()`](https:/pharmaverse.github.io/admiral/main/reference/basket_select.md)
      object, a data frame, or a list of data frames and
      [`basket_select()`](https:/pharmaverse.github.io/admiral/main/reference/basket_select.md)
      objects.

  Default value

  :   `NULL`

## Value

An object of class `query`.

## See also

[`create_query_data()`](https:/pharmaverse.github.io/admiral/main/reference/create_query_data.md),
[`basket_select()`](https:/pharmaverse.github.io/admiral/main/reference/basket_select.md),
[`vignette("queries_dataset")`](https:/pharmaverse.github.io/admiral/main/articles/queries_dataset.md)

Source Objects:
[`basket_select()`](https:/pharmaverse.github.io/admiral/main/reference/basket_select.md),
[`censor_source()`](https:/pharmaverse.github.io/admiral/main/reference/censor_source.md),
[`death_event`](https:/pharmaverse.github.io/admiral/main/reference/tte_source_objects.md),
[`event()`](https:/pharmaverse.github.io/admiral/main/reference/event.md),
[`event_joined()`](https:/pharmaverse.github.io/admiral/main/reference/event_joined.md),
[`event_source()`](https:/pharmaverse.github.io/admiral/main/reference/event_source.md),
[`flag_event()`](https:/pharmaverse.github.io/admiral/main/reference/flag_event.md),
[`records_source()`](https:/pharmaverse.github.io/admiral/main/reference/records_source.md),
[`tte_source()`](https:/pharmaverse.github.io/admiral/main/reference/tte_source.md)

## Examples

``` r
# create a query for an SMQ
library(tibble)
library(dplyr, warn.conflicts = FALSE)

# create a query for a SMQ
query(
  prefix = "SMQ02",
  id = auto,
  definition = basket_select(
    name = "Pregnancy and neonatal topics (SMQ)",
    scope = "NARROW",
    type = "smq"
  )
)
#> <query> object
#> prefix: "SMQ02"
#> name: auto
#> id: auto
#> add_scope_num: FALSE
#> definition:
#>   <basket_select> object
#>   name: "Pregnancy and neonatal topics (SMQ)"
#>   id: NULL
#>   scope: "NARROW"
#>   type: "smq"

# create a query for an SDG
query(
  prefix = "SDG01",
  id = auto,
  definition = basket_select(
    name = "5-aminosalicylates for ulcerative colitis",
    scope = NA_character_,
    type = "sdg"
  )
)
#> <query> object
#> prefix: "SDG01"
#> name: auto
#> id: auto
#> add_scope_num: FALSE
#> definition:
#>   <basket_select> object
#>   name: "5-aminosalicylates for ulcerative colitis"
#>   id: NULL
#>   scope: "NA"
#>   type: "sdg"

# creating a query for a customized query
cqterms <- tribble(
  ~TERMCHAR, ~TERMNUM,
  "APPLICATION SITE ERYTHEMA", 10003041L,
  "APPLICATION SITE PRURITUS", 10003053L
) %>%
  mutate(SRCVAR = "AEDECOD")

query(
  prefix = "CQ01",
  name = "Application Site Issues",
  definition = cqterms
)
#> <query> object
#> prefix: "CQ01"
#> name: "Application Site Issues"
#> add_scope_num: FALSE
#> definition:
#> # A tibble: 2 × 3
#>   TERMCHAR                   TERMNUM SRCVAR 
#>   <chr>                        <int> <chr>  
#> 1 APPLICATION SITE ERYTHEMA 10003041 AEDECOD
#> 2 APPLICATION SITE PRURITUS 10003053 AEDECOD

# creating a customized query based on SMQs and additional terms
query(
  prefix = "CQ03",
  name = "Special issues of interest",
  definition = list(
    cqterms,
    basket_select(
      name = "Pregnancy and neonatal topics (SMQ)",
      scope = "NARROW",
      type = "smq"
    ),
    basket_select(
      id = 8050L,
      scope = "BROAD",
      type = "smq"
    )
  )
)
#> <query> object
#> prefix: "CQ03"
#> name: "Special issues of interest"
#> add_scope_num: FALSE
#> definition:
#>   c("APPLICATION SITE ERYTHEMA", "APPLICATION SITE PRURITUS")
#> c(10003041, 10003053)
#> c("AEDECOD", "AEDECOD")
#>   Pregnancy and neonatal topics (SMQ)
#> NULL
#> NARROW
#> smq
#>   NULL
#> 8050
#> BROAD
#> smq
```
