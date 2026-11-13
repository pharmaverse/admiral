# Create an `admiral_df` Dataset

**\[experimental\]**

Tags a data frame with the `admiral_df` class and stores its metadata:
the key variables – the variables which uniquely identify its records –
the dataset name, and the ADaM dataset class. Call it on the finished
dataset, once it has been derived, e.g. just before summarizing it: the
metadata is not guaranteed to survive the derivations which build the
dataset (see Details). The metadata is read back with
[`get_dataset_keys()`](https:/pharmaverse.github.io/admiral/3193-test-refactor-if-else/reference/get_dataset_keys.md),
[`get_dataset_name()`](https:/pharmaverse.github.io/admiral/3193-test-refactor-if-else/reference/get_dataset_name.md),
and
[`get_dataset_class()`](https:/pharmaverse.github.io/admiral/3193-test-refactor-if-else/reference/get_dataset_class.md).

Called on a dataset which is already an `admiral_df`, it updates the
metadata rather than rebuilding it: an argument left as `NULL` keeps the
value stored by the previous call, so each item can be set on its own.

The keys describe the record structure the dataset is intended to have;
they are not checked against its records when they are stored.
`{admiral}` never infers the keys of a dataset: if they have not been
set with this function, functions which need them ask for them to be
set. The class, by contrast, is inferred from the dataset if it has not
been set (see
[`infer_dataset_class()`](https:/pharmaverse.github.io/admiral/3193-test-refactor-if-else/reference/infer_dataset_class.md)),
but a class set here always takes precedence.

## Usage

``` r
admiral_df(
  dataset,
  dataset_keys = NULL,
  dataset_name = NULL,
  dataset_class = NULL
)
```

## Arguments

- dataset:

  A data frame. It does not need to have the `admiral_df` class already;
  the class is added if it is missing.

  Permitted values

  :   a data frame

  Default value

  :   none

- dataset_keys:

  The key variables of `dataset`, in key order, e.g.
  `exprs(USUBJID, PARAMCD, AVISIT)`, stored in the
  `"admiral_dataset_keys"` attribute. They must be variables of
  `dataset`: if any is not, a warning is issued and the keys are not
  set, so any keys set by a previous call are kept. An empty list,
  [`exprs()`](https:/pharmaverse.github.io/admiral/3193-test-refactor-if-else/reference/reexport-exprs.md),
  records that the dataset has no key variables defined, which is
  distinct from its keys never having been set. `NULL` leaves keys
  stored by a previous call untouched.

  Permitted values

  :   list of variables created by
      [`exprs()`](https:/pharmaverse.github.io/admiral/3193-test-refactor-if-else/reference/reexport-exprs.md)
      or `NULL`

  Default value

  :   `NULL`

- dataset_name:

  The name of the dataset, e.g. `"ADVS"`, stored in the
  `"admiral_dataset_name"` attribute. There are no restrictions on its
  casing: it is stored exactly as given. `NULL` leaves a name stored by
  a previous call untouched, so re-keying a dataset does not lose its
  name.

  Permitted values

  :   character scalar or `NULL`

  Default value

  :   `NULL`

- dataset_class:

  The ADaM dataset class, stored in the `"admiral_dataset_class"`
  attribute. This is the ADaM class of the dataset, not its R class
  (which is `admiral_df`). It must be given in upper case, as listed
  below. `"TTE"` is the time-to-event subclass of `"BDS"`. `"ADDL"`,
  `"MDBDS"`, `"MDTTE"`, and `"MDOCCDS"` are the medical device classes
  (device-level, and the device versions of `"BDS"`, `"TTE"`, and
  `"OCCDS"`). The class is not checked against the variables of
  `dataset`. `NULL` leaves a class stored by a previous call untouched.

  Permitted values

  :   `"ADSL"`, `"BDS"`, `"OCCDS"`, `"TTE"`, `"ADDL"`, `"MDBDS"`,
      `"MDTTE"`, `"MDOCCDS"`, or `NULL`

  Default value

  :   `NULL`

## Value

`dataset` with the `admiral_df` class, and the `"admiral_dataset_keys"`,
`"admiral_dataset_name"`, and `"admiral_dataset_class"` attributes for
those of `dataset_keys`, `dataset_name`, and `dataset_class` which were
supplied (or had been stored by a previous call).

## Details

The keys can come from anywhere: the `by_vars` of the derivation which
produced the dataset, a dataset specification, or a hand-written list. A
specification provides them as a character vector, which
[`chr2vars()`](https:/pharmaverse.github.io/admiral/3193-test-refactor-if-else/reference/chr2vars.md)
converts. `{admiral}` does not read specification objects itself; to use
the keys defined in a `{metacore}` specification, extract them in the
calling code, e.g.
`chr2vars(metacore::get_keys(spec, "ADVS")$variable)`.

All `{admiral}` dataset functions tag their output with the `admiral_df`
class, but they do not carry the metadata through. The core `{dplyr}`
verbs (`mutate()`, [`filter()`](https://rdrr.io/r/stats/filter.html),
`arrange()`, `select()`, `rename()`, `slice()`, `distinct()`, joins
where `dataset` is `x`, ...) restore data frame attributes, so the
metadata survives those. It is *not* preserved by operations which build
a new data frame, in particular:

- `summarise()`,

- [`tidyr::pivot_longer()`](https://tidyr.tidyverse.org/reference/pivot_longer.html)/`pivot_wider()`,
  `nest()`/`unnest()`, `complete()`,

- `bind_rows()` when `dataset` is not the *first* argument,

- joins where `dataset` is the `y` (rather than the `x`) argument,

- `ungroup()`, after `group_by()` or `rowwise()` (which themselves keep
  the attributes but drop the `admiral_df` class).

Many derivations use these internally, so metadata set at the start of a
derivation pipeline is lost along the way. Set it once the dataset is
finished, or re-apply this function if the dataset changes afterwards.

The keys are stored as plain variable names, so they are not updated by
`rename()` or removed by `select()`: a key which is no longer in the
dataset means it has changed since its keys were set.

## See also

Utilities used within Derivation functions:
[`extract_unit()`](https:/pharmaverse.github.io/admiral/3193-test-refactor-if-else/reference/extract_unit.md),
[`get_dataset_class()`](https:/pharmaverse.github.io/admiral/3193-test-refactor-if-else/reference/get_dataset_class.md),
[`get_dataset_keys()`](https:/pharmaverse.github.io/admiral/3193-test-refactor-if-else/reference/get_dataset_keys.md),
[`get_dataset_name()`](https:/pharmaverse.github.io/admiral/3193-test-refactor-if-else/reference/get_dataset_name.md),
[`get_flagged_records()`](https:/pharmaverse.github.io/admiral/3193-test-refactor-if-else/reference/get_flagged_records.md),
[`get_not_mapped()`](https:/pharmaverse.github.io/admiral/3193-test-refactor-if-else/reference/get_not_mapped.md),
[`get_vars_query()`](https:/pharmaverse.github.io/admiral/3193-test-refactor-if-else/reference/get_vars_query.md)

## Examples

### Tagging a dataset with its metadata

    library(tibble)

    advs <- tribble(
      ~USUBJID, ~PARAMCD, ~AVISIT,    ~AVAL,
      "1",      "DIABP",  "BASELINE",    51,
      "1",      "SYSBP",  "BASELINE",   121,
      "2",      "DIABP",  "BASELINE",    79
    )

    advs <- admiral_df(
      advs,
      dataset_keys = exprs(USUBJID, PARAMCD, AVISIT),
      dataset_name = "ADVS",
      dataset_class = "BDS"
    )

    get_dataset_keys(advs)
    #> [[1]]
    #> USUBJID
    #>
    #> [[2]]
    #> PARAMCD
    #>
    #> [[3]]
    #> AVISIT
    #>
    get_dataset_name(advs)
    #> [1] "ADVS"
    get_dataset_class(advs)
    #> [1] "BDS"
    class(advs)
    #> [1] "admiral_df" "tbl_df"     "tbl"        "data.frame"

### Keys from a dataset specification

A dataset specification provides the keys as a character vector, which
[`chr2vars()`](https:/pharmaverse.github.io/admiral/3193-test-refactor-if-else/reference/chr2vars.md)
converts. `{admiral}` does not read specification objects itself, so
they are extracted in the calling code; e.g. with a `{metacore}`
specification that is `metacore::get_keys()`:

    admiral_df(
      advs,
      dataset_keys = chr2vars(metacore::get_keys(spec, "ADVS")$variable),
      dataset_name = "ADVS"
    )

    advs <- admiral_df(
      advs,
      dataset_keys = chr2vars(c("USUBJID", "PARAMCD", "AVISIT"))
    )

    # `dataset_name` was not repeated, so the name from the previous call stands
    get_dataset_name(advs)
    #> [1] "ADVS"

### Declaring that a dataset has no key variables

An empty `dataset_keys` is stored all the same. It is not the same as
leaving the attribute unset: it records that no key variables are
defined.

    advs_no_keys <- admiral_df(advs, dataset_keys = exprs())

    get_dataset_keys(advs_no_keys)
    #> list()
