# Higher Order Functions

## Introduction

This vignette explains some of the more advanced options of
[admiral](https://pharmaverse.github.io/admiral/) related to higher
order functions. A higher order function is a function that takes
another function as input. By introducing these higher order functions,
we intend to give the user greater power over our derivation functions,
whilst trying to negate the need for adding additional
[admiral](https://pharmaverse.github.io/admiral/) functions or
arguments, or the user needing many separate steps.

The functions covered here are:

- [`call_derivation()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/call_derivation.md):
  Call a single derivation multiple times with some arguments being
  fixed across iterations and others varying
- [`restrict_derivation()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/restrict_derivation.md):
  Execute a single derivation on a subset of the input dataset
- [`slice_derivation()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/slice_derivation.md):
  The input dataset is split into slices (subsets) and for each slice a
  single derivation is called separately. Some or all arguments of the
  derivation may vary depending on the slice.

The examples in this vignette all showcase the higher order functions
taking other [admiral](https://pharmaverse.github.io/admiral/) functions
as input, but in principle any functions that modify a dataset (e.g. an
extension package function, or
[`dplyr::mutate()`](https://dplyr.tidyverse.org/reference/mutate.html),
etc.) can be passed too. The only requirement for a function being
passed to `derivation` is that it must take a dataset as its first
argument and return a dataset.

### Required Packages

The examples of this vignette require the following packages.

For example purpose, the ADSL dataset—which is included in
[admiral](https://pharmaverse.github.io/admiral/)—and the SDTM datasets
from [pharmaversesdtm](https://pharmaverse.github.io/pharmaversesdtm/)
are used.

``` r
library(admiral)
library(pharmaversesdtm)
library(dplyr, warn.conflicts = FALSE)

ae <- pharmaversesdtm::ae
vs <- pharmaversesdtm::vs
adsl <- admiral::admiral_adsl

ae <- convert_blanks_to_na(ae)
vs <- convert_blanks_to_na(vs)
```

The following code creates a minimally viable ADAE dataset to be used
where needed in the following examples.

``` r
adae <- ae %>%
  left_join(adsl, by = c("STUDYID", "USUBJID")) %>%
  derive_vars_dt(
    new_vars_prefix = "AST",
    dtc = AESTDTC,
    highest_imputation = "M"
  ) %>%
  mutate(TRTEMFL = if_else(ASTDT >= TRTSDT, "Y", NA_character_))
```

## Call Derivation

This function exists purely for convenience to save the user repeating
numerous similar derivation function calls. It is best used when
multiple derived variables have very similar specifications with only
slight variations.

As an example, imagine the case where all the parameters in a BDS ADaM
required both a highest value flag and a lowest value flag.

Here is an example of how to achieve this **without** using
[`call_derivation()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/call_derivation.md):

``` r
vs_without <- vs %>%
  derive_var_extreme_flag(
    by_vars = exprs(USUBJID, VSTESTCD),
    order = exprs(VSORRES, VSSEQ),
    new_var = AHIFL,
    mode = "last"
  ) %>%
  derive_var_extreme_flag(
    by_vars = exprs(USUBJID, VSTESTCD),
    order = exprs(VSORRES, VSSEQ),
    new_var = ALOFL,
    mode = "first"
  )
```

Here is an example of how to achieve the same **with** using
[`call_derivation()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/call_derivation.md),
where any different arguments are passed using
[`params()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/params.md):

``` r
vs_with <- vs %>%
  call_derivation(
    derivation = derive_var_extreme_flag,
    variable_params = list(
      params(new_var = AHIFL, mode = "last"),
      params(new_var = ALOFL, mode = "first")
    ),
    by_vars = exprs(USUBJID, VSTESTCD),
    order = exprs(VSORRES, VSSEQ)
  )
```

In the example, you can see how in these higher order functions,
`derivation` is where the user supplies the name of the derivation
function to apply, with no trailing parentheses required. Then
`variable_params` is used to pass a list of the different arguments
needed for each derived variable.

The advantage of this higher order function would be further highlighted
with examples where more than two variable derivations had similar
needs, such as the below case where multiple time to AE parameters are
derived in one call. Note that this example relies on pre-defined
`tte_source` objects, as explained at [Creating a BDS Time-to-Event
ADaM](https:/pharmaverse.github.io/admiral/v1.4.1/articles/bds_tte.md).

``` r
adaette <- call_derivation(
  derivation = derive_param_tte,
  variable_params = list(
    params(
      event_conditions = list(ae_event),
      set_values_to = exprs(PARAMCD = "TTAE")
    ),
    params(
      event_conditions = list(ae_ser_event),
      set_values_to = exprs(PARAMCD = "TTSERAE")
    ),
    params(
      event_conditions = list(ae_sev_event),
      set_values_to = exprs(PARAMCD = "TTSEVAE")
    ),
    params(
      event_conditions = list(ae_wd_event),
      set_values_to = exprs(PARAMCD = "TTWDAE")
    )
  ),
  dataset_adsl = adsl,
  source_datasets = list(adsl = adsl, adae = adae),
  censor_conditions = list(lastalive_censor)
)
#> Warning: Dataset "adae" contains duplicate records with respect to `STUDYID`, `USUBJID`,
#> and `ASTDT`
#> ℹ Run `admiral::get_duplicates_dataset()` to access the duplicate records
```

Developing your ADaM scripts this way using
[`call_derivation()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/call_derivation.md)
could give the following benefits:

- code becomes more efficient and readable
- maintenance would be eased in case of specification changes
- downstream quality checking would require less effort

## Restrict Derivation

The idea behind this function is that sometimes you want to apply a
derivation only for certain records from the input dataset. Introducing
[`restrict_derivation()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/restrict_derivation.md)
therefore gives the users the ability to achieve this across any
function, without each function needing to have such an argument to
allow for this.

An example would be if you wanted to flag the first occurring AE with
the highest severity for each patient, but you only wanted to do this
for records occurring on or after study day 1.

Here is how you could achieve this using
[`restrict_derivation()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/restrict_derivation.md),
where the function arguments are passed using
[`params()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/params.md)
and the restriction criteria is given using `filter`:

``` r
ae <- ae %>%
  mutate(TEMP_AESEVN = as.integer(factor(AESEV, levels = c("SEVERE", "MODERATE", "MILD")))) %>%
  restrict_derivation(
    derivation = derive_var_extreme_flag,
    args = params(
      new_var = AHSEVFL,
      by_vars = exprs(USUBJID),
      order = exprs(TEMP_AESEVN, AESTDY, AESEQ),
      mode = "first"
    ),
    filter = AESTDY >= 1
  )
```

## Slice Derivation

This function in a way combines the features of the above two. It allows
a single derivation to be applied with different arguments for different
slices (subsets) of records from the input dataset. You could do this
with separate
[`restrict_derivation()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/restrict_derivation.md)
calls for each different set of records, but
[`slice_derivation()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/slice_derivation.md)
allows to achieve this in one call.

An example would be if you wanted to achieve the same derivation as
above for records occurring on or after study day 1, but for
pre-treatment AEs you wanted to flag only the last occurring AE.

Here is how you could achieve this using
[`slice_derivation()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/slice_derivation.md),
where the function arguments are passed using
[`params()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/params.md)
and via the different slices controlled by `filter`:

``` r
ae <- ae %>%
  slice_derivation(
    derivation = derive_var_extreme_flag,
    args = params(
      new_var = AHSEV2FL,
      by_vars = exprs(USUBJID)
    ),
    derivation_slice(
      filter = AESTDY >= 1,
      args = params(order = exprs(TEMP_AESEVN, AESTDY, AESEQ), mode = "first")
    ),
    derivation_slice(
      filter = TRUE,
      args = params(order = exprs(AESTDY, AESEQ), mode = "last")
    )
  )
```

As you can see in the example, the `derivation_slice` ordering is
important. Here we addressed all the AEs on or after study day 1 first,
and then we used `filter = TRUE` option to catch all remaining records
(in this case pre-treatment AEs).

The ordering is perhaps shown even more when we look at the below
example where three slices are taken. Remember that observations that
match with more than one slice are only considered for the first
matching slice. So in this case we’re creating a flag for each patient
for the record with the first severe AE, and then the first moderate AE,
and finally flagging the last occurring AE where not severe or moderate.

``` r
ae <- ae %>%
  slice_derivation(
    derivation = derive_var_extreme_flag,
    args = params(
      new_var = AHSEV3FL,
      by_vars = exprs(USUBJID)
    ),
    derivation_slice(
      filter = AESEV == "SEVERE",
      args = params(order = exprs(AESTDY, AESEQ), mode = "first")
    ),
    derivation_slice(
      filter = AESEV == "MODERATE",
      args = params(order = exprs(AESTDY, AESEQ), mode = "first")
    ),
    derivation_slice(
      filter = TRUE,
      args = params(order = exprs(AESTDY, AESEQ), mode = "last")
    )
  )
```

The order is only important when the slices are not mutually exclusive,
so in the above case the moderate AE slice could have been above the
severe AE slice, for example, and there would have been no difference to
the result. However the third slice had to come last to check all
remaining (i.e. not severe or moderate) records only.

## Combining Higher Order Functions

It is also possible to use two or more of the higher order functions in
combination. The most likely use case is employing
[`call_derivation()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/call_derivation.md)
to repeat similar derivations multiple times, but within a call to
[`restrict_derivation()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/restrict_derivation.md)
to execute the derivations on a subset of the input dataset.

For instance, let’s pick up VS from the
[`call_derivation()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/call_derivation.md)
example and suppose that we wish to derive the highest value and lowest
value flags for just the temperature parameter. We can do this by
ensuring that the derivation passed to
[`restrict_derivation()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/restrict_derivation.md)
is
[`call_derivation()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/call_derivation.md),
and that the arguments for the latter are all passed through `args`.

``` r
vs_hilotemp <- vs %>%
  restrict_derivation(
    derivation = call_derivation,
    args = params(
      derivation = derive_var_extreme_flag,
      variable_params = list(
        params(new_var = ATMPHIFL, mode = "last"),
        params(new_var = ATMPLOFL, mode = "first")
      ),
      by_vars = exprs(USUBJID, VSTESTCD),
      order = exprs(VSORRES, VSSEQ)
    ),
    filter = VSTESTCD == "TEMP"
  )
```
