# Creates a queries dataset as input dataset to the `dataset_queries` argument in `derive_vars_query()`

Creates a queries dataset as input dataset to the `dataset_queries`
argument in the
[`derive_vars_query()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/derive_vars_query.md)
function as defined in the
[`vignette("queries_dataset")`](https:/pharmaverse.github.io/admiral/v1.4.1/articles/queries_dataset.md).

## Usage

``` r
create_query_data(queries, version = NULL, get_terms_fun = NULL)
```

## Arguments

- queries:

  List of queries

  A list of
  [`query()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/query.md)
  objects is expected.

  Default value

  :   none

- version:

  Dictionary version

  The dictionary version used for coding the terms should be specified.
  If any of the queries is a basket (SMQ, SDG, ....) or a customized
  query including a basket, the parameter needs to be specified.

  Permitted values

  :   A character string (the expected format is company-specific)

  Default value

  :   `NULL`

- get_terms_fun:

  Function which returns the terms

  For each query specified for the `queries` parameter referring to a
  basket (i.e., those where the `definition` field is set to a
  [`basket_select()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/basket_select.md)
  object or a list which contains at least one
  [`basket_select()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/basket_select.md)
  object) the specified function is called to retrieve the terms
  defining the query. This function is not provided by admiral as it is
  company specific, i.e., it has to be implemented at company level.

  The function must return a dataset with all the terms defining the
  basket. The output dataset must contain the following variables.

  - `SRCVAR`: the variable to be used for defining a term of the basket,
    e.g., `AEDECOD`

  - `TERMCHAR`: the name of the term if the variable `SRCVAR` is
    referring to is character

  - `TERMNUM` the numeric id of the term if the variable `SRCVAR` is
    referring to is numeric

  - `GRPNAME`: the name of the basket. The values must be the same for
    all observations.

  The function must provide the following parameters

  - `basket_select`: A
    [`basket_select()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/basket_select.md)
    object.

  - `version`: The dictionary version. The value specified for the
    `version` in the `create_query_data()` call is passed to this
    parameter.

  - `keep_id`: If set to `TRUE`, the output dataset must contain the
    `GRPID` variable. The variable must be set to the numeric id of the
    basket.

  - `temp_env`: A temporary environment is passed to this parameter. It
    can be used to store data which is used for all baskets in the
    `create_query_data()` call. For example if SMQs need to be read from
    a database all SMQs can be read and stored in the environment when
    the first SMQ is handled. For the other SMQs the terms can be
    retrieved from the environment instead of accessing the database
    again.

  Default value

  :   `NULL`

## Value

A dataset to be used as input dataset to the `dataset_queries` argument
in
[`derive_vars_query()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/derive_vars_query.md)

## Details

For each
[`query()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/query.md)
object listed in the `queries` argument, the terms belonging to the
query (`SRCVAR`, `TERMCHAR`, `TERMNUM`) are determined with respect to
the `definition` field of the query: if the definition field of the
[`query()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/query.md)
object is

- a
  [`basket_select()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/basket_select.md)
  object, the terms are read from the basket database by calling the
  function specified for the `get_terms_fun` parameter.

- a data frame, the terms stored in the data frame are used.

- a list of data frames and
  [`basket_select()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/basket_select.md)
  objects, all terms from the data frames and all terms read from the
  basket database referenced by the
  [`basket_select()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/basket_select.md)
  objects are collated.

The following variables (as described in
[`vignette("queries_dataset")`](https:/pharmaverse.github.io/admiral/v1.4.1/articles/queries_dataset.md)
are created:

- `PREFIX`: Prefix of the variables to be created by
  [`derive_vars_query()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/derive_vars_query.md)
  as specified by the `prefix` element.

- `GRPNAME`: Name of the query as specified by the `name` element.

- `GRPID`: Id of the query as specified by the `id` element. If the `id`
  element is not specified for a query, the variable is set to `NA`. If
  the `id` element is not specified for any query, the variable is not
  created.

- `SCOPE`: scope of the query as specified by the `scope` element of the
  [`basket_select()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/basket_select.md)
  object. For queries not defined by a
  [`basket_select()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/basket_select.md)
  object, the variable is set to `NA`. If none of the queries is defined
  by a
  [`basket_select()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/basket_select.md)
  object, the variable is not created.

- `SCOPEN`: numeric scope of the query. It is set to `1` if the scope is
  broad. Otherwise it is set to `2`. If the `add_scope_num` element
  equals `FALSE`, the variable is set to `NA`. If the `add_scope_num`
  element equals `FALSE` for all baskets or none of the queries is an
  basket , the variable is not created.

- `SRCVAR`: Name of the variable used to identify the terms.

- `TERMCHAR`: Value of the term variable if it is a character variable.

- `TERMNUM`: Value of the term variable if it is a numeric variable.

- `VERSION`: Set to the value of the `version` argument. If it is not
  specified, the variable is not created.

## See also

[`derive_vars_query()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/derive_vars_query.md),
[`query()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/query.md),
[`basket_select()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/basket_select.md),
[`vignette("queries_dataset")`](https:/pharmaverse.github.io/admiral/v1.4.1/articles/queries_dataset.md)

Creating auxiliary datasets:
[`consolidate_metadata()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/consolidate_metadata.md),
[`create_period_dataset()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/create_period_dataset.md),
[`create_single_dose_dataset()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/create_single_dose_dataset.md)

## Examples

``` r
library(tibble)
library(dplyr, warn.conflicts = FALSE)
library(pharmaversesdtm)
library(admiral)

# creating a query dataset for a customized query
cqterms <- tribble(
  ~TERMCHAR, ~TERMNUM,
  "APPLICATION SITE ERYTHEMA", 10003041L,
  "APPLICATION SITE PRURITUS", 10003053L
) %>%
  mutate(SRCVAR = "AEDECOD")

cq <- query(
  prefix = "CQ01",
  name = "Application Site Issues",
  definition = cqterms
)

create_query_data(queries = list(cq))
#> # A tibble: 2 × 5
#>   TERMCHAR                   TERMNUM SRCVAR  PREFIX GRPNAME                
#>   <chr>                        <int> <chr>   <chr>  <chr>                  
#> 1 APPLICATION SITE ERYTHEMA 10003041 AEDECOD CQ01   Application Site Issues
#> 2 APPLICATION SITE PRURITUS 10003053 AEDECOD CQ01   Application Site Issues

# create a query dataset for SMQs
pregsmq <- query(
  prefix = "SMQ02",
  id = auto,
  definition = basket_select(
    name = "Pregnancy and neonatal topics (SMQ)",
    scope = "NARROW",
    type = "smq"
  )
)

bilismq <- query(
  prefix = "SMQ04",
  definition = basket_select(
    id = 20000121L,
    scope = "BROAD",
    type = "smq"
  )
)

# The get_terms function from pharmaversesdtm is used for this example.
# In a real application a company-specific function must be used.
create_query_data(
  queries = list(pregsmq, bilismq),
  get_terms_fun = pharmaversesdtm:::get_terms,
  version = "20.1"
)
#> # A tibble: 43 × 7
#>    TERMCHAR                           SRCVAR GRPNAME  GRPID SCOPE PREFIX VERSION
#>    <chr>                              <chr>  <chr>    <int> <chr> <chr>  <chr>  
#>  1 Achromotrichia congenital          AEDEC… Pregna… 2.00e7 NARR… SMQ02  20.1   
#>  2 Craniosynostosis                   AEDEC… Pregna… 2.00e7 NARR… SMQ02  20.1   
#>  3 Hypophosphatasia                   AEDEC… Pregna… 2.00e7 NARR… SMQ02  20.1   
#>  4 Congenital pyelocaliectasis        AEDEC… Pregna… 2.00e7 NARR… SMQ02  20.1   
#>  5 Uterine contractions during pregn… AEDEC… Pregna… 2.00e7 NARR… SMQ02  20.1   
#>  6 Ductus arteriosus premature closu… AEDEC… Pregna… 2.00e7 NARR… SMQ02  20.1   
#>  7 Pseudotruncus arteriosus           AEDEC… Pregna… 2.00e7 NARR… SMQ02  20.1   
#>  8 Lipomeningocele                    AEDEC… Pregna… 2.00e7 NARR… SMQ02  20.1   
#>  9 Macrocephaly                       AEDEC… Pregna… 2.00e7 NARR… SMQ02  20.1   
#> 10 Carnitine palmitoyltransferase de… AEDEC… Pregna… 2.00e7 NARR… SMQ02  20.1   
#> # ℹ 33 more rows

# create a query dataset for SDGs
sdg <- query(
  prefix = "SDG01",
  id = auto,
  definition = basket_select(
    name = "5-aminosalicylates for ulcerative colitis",
    scope = NA_character_,
    type = "sdg"
  )
)

# The get_terms function from pharmaversesdtm is used for this example.
# In a real application a company-specific function must be used.
create_query_data(
  queries = list(sdg),
  get_terms_fun = pharmaversesdtm:::get_terms,
  version = "2019-09"
)
#> # A tibble: 16 × 7
#>    TERMCHAR                            SRCVAR GRPNAME GRPID SCOPE PREFIX VERSION
#>    <chr>                               <chr>  <chr>   <int> <chr> <chr>  <chr>  
#>  1 AMINOSALICYLIC ACID                 CMDEC… 5-amin…   220 NA    SDG01  2019-09
#>  2 AMINOSALICYLATE CALCIUM             CMDEC… 5-amin…   220 NA    SDG01  2019-09
#>  3 AMINOSALICYLATE CALCIUM ALUMINIUM   CMDEC… 5-amin…   220 NA    SDG01  2019-09
#>  4 AMINOSALICYLATE SODIUM              CMDEC… 5-amin…   220 NA    SDG01  2019-09
#>  5 SODIUM AMINOSALICYLATE DIHYDRATE    CMDEC… 5-amin…   220 NA    SDG01  2019-09
#>  6 AMINOSALICYLATE SODIUM;AMINOSALICY… CMDEC… 5-amin…   220 NA    SDG01  2019-09
#>  7 SULFASALAZINE                       CMDEC… 5-amin…   220 NA    SDG01  2019-09
#>  8 CALCIUM BENZAMIDOSALICYLATE         CMDEC… 5-amin…   220 NA    SDG01  2019-09
#>  9 OLSALAZINE                          CMDEC… 5-amin…   220 NA    SDG01  2019-09
#> 10 OLSALAZINE SODIUM                   CMDEC… 5-amin…   220 NA    SDG01  2019-09
#> 11 MESALAZINE                          CMDEC… 5-amin…   220 NA    SDG01  2019-09
#> 12 BALSALAZIDE                         CMDEC… 5-amin…   220 NA    SDG01  2019-09
#> 13 BALSALAZIDE SODIUM                  CMDEC… 5-amin…   220 NA    SDG01  2019-09
#> 14 BALSALAZIDE DISODIUM DIHYDRATE      CMDEC… 5-amin…   220 NA    SDG01  2019-09
#> 15 DERSALAZINE                         CMDEC… 5-amin…   220 NA    SDG01  2019-09
#> 16 DERSALAZINE SODIUM                  CMDEC… 5-amin…   220 NA    SDG01  2019-09

# creating a query dataset for a customized query including SMQs
# The get_terms function from pharmaversesdtm is used for this example.
# In a real application a company-specific function must be used.
create_query_data(
  queries = list(
    query(
      prefix = "CQ03",
      name = "Special issues of interest",
      definition = list(
        basket_select(
          name = "Pregnancy and neonatal topics (SMQ)",
          scope = "NARROW",
          type = "smq"
        ),
        cqterms
      )
    )
  ),
  get_terms_fun = pharmaversesdtm:::get_terms,
  version = "20.1"
)
#> # A tibble: 23 × 6
#>    TERMCHAR                                SRCVAR GRPNAME TERMNUM PREFIX VERSION
#>    <chr>                                   <chr>  <chr>     <int> <chr>  <chr>  
#>  1 Achromotrichia congenital               AEDEC… Specia…      NA CQ03   20.1   
#>  2 Craniosynostosis                        AEDEC… Specia…      NA CQ03   20.1   
#>  3 Hypophosphatasia                        AEDEC… Specia…      NA CQ03   20.1   
#>  4 Congenital pyelocaliectasis             AEDEC… Specia…      NA CQ03   20.1   
#>  5 Uterine contractions during pregnancy   AEDEC… Specia…      NA CQ03   20.1   
#>  6 Ductus arteriosus premature closure     AEDEC… Specia…      NA CQ03   20.1   
#>  7 Pseudotruncus arteriosus                AEDEC… Specia…      NA CQ03   20.1   
#>  8 Lipomeningocele                         AEDEC… Specia…      NA CQ03   20.1   
#>  9 Macrocephaly                            AEDEC… Specia…      NA CQ03   20.1   
#> 10 Carnitine palmitoyltransferase deficie… AEDEC… Specia…      NA CQ03   20.1   
#> # ℹ 13 more rows
```
