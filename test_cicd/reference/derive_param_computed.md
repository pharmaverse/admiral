# Adds a Parameter Computed from the Analysis Value of Other Parameters

Adds a parameter computed from the analysis value of other parameters.
It is expected that the analysis value of the new parameter is defined
by an expression using the analysis values of other parameters, such as
addition/sum, subtraction/difference, multiplication/product,
division/ratio, exponentiation/logarithm, or by formula.  
  
For example mean arterial pressure (MAP) can be derived from systolic
(SYSBP) and diastolic blood pressure (DIABP) with the formula \$\$MAP =
\frac{SYSBP + 2 DIABP}{3}\$\$

## Usage

``` r
derive_param_computed(
  dataset = NULL,
  dataset_add = NULL,
  by_vars,
  parameters,
  set_values_to,
  filter = NULL,
  constant_by_vars = NULL,
  constant_parameters = NULL,
  keep_nas = FALSE
)
```

## Arguments

- dataset:

  Input dataset

  The variables specified by the `by_vars` argument are expected to be
  in the dataset. `PARAMCD` is expected as well.

  The variable specified by `by_vars` and `PARAMCD` must be a unique key
  of the input dataset after restricting it by the filter condition
  (`filter` parameter) and to the parameters specified by `parameters`.

  Permitted values

  :   a dataset, i.e., a `data.frame` or tibble

  Default value

  :   `NULL`

- dataset_add:

  Additional dataset

  The variables specified by the `by_vars` parameter are expected.

  The variable specified by `by_vars` and `PARAMCD` must be a unique key
  of the additional dataset after restricting it to the parameters
  specified by `parameters`.

  If the argument is specified, the observations of the additional
  dataset are considered in addition to the observations from the input
  dataset (`dataset` restricted by `filter`).

  Permitted values

  :   a dataset, i.e., a `data.frame` or tibble

  Default value

  :   `NULL`

- by_vars:

  Grouping variables

  For each group defined by `by_vars` an observation is added to the
  output dataset. Only variables specified in `by_vars` will be
  populated in the newly created records.

  Permitted values

  :   list of variables created by
      [`exprs()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/reexport-exprs.md),
      e.g., `exprs(USUBJID, VISIT)`

  Default value

  :   none

- parameters:

  Required parameter codes

  It is expected that all parameter codes (`PARAMCD`) which are required
  to derive the new parameter are specified for this parameter or the
  `constant_parameters` parameter.

  If observations should be considered which do not have a parameter
  code, e.g., if an SDTM dataset is used, temporary parameter codes can
  be derived by specifying a list of expressions. The name of the
  element defines the temporary parameter code and the expression the
  condition for selecting the records. For example
  `parameters = exprs(HGHT = VSTESTCD == "HEIGHT")` selects the
  observations with `VSTESTCD == "HEIGHT"` from the input data
  (`dataset` and `dataset_add`), sets `PARAMCD = "HGHT"` for these
  observations, and adds them to the observations to consider.

  Unnamed elements in the list of expressions are considered as
  parameter codes. For example,
  `parameters = exprs(WEIGHT, HGHT = VSTESTCD == "HEIGHT")` uses the
  parameter code `"WEIGHT"` and creates a temporary parameter code
  `"HGHT"`.

  Permitted values

  :   A character vector of `PARAMCD` values or a list of expressions

  Default value

  :   none

- set_values_to:

  Variables to be set

  The specified variables are set to the specified values for the new
  observations. The values of variables of the parameters specified by
  `parameters` can be accessed using `<variable name>.<parameter code>`.
  For example

      exprs(
        AVAL = (AVAL.SYSBP + 2 * AVAL.DIABP) / 3,
        PARAMCD = "MAP"
      )

  defines the analysis value and parameter code for the new parameter.

  Variable names in the expression must not contain more than one dot.

  Note that `dplyr` helper functions such as
  [`dplyr::starts_with()`](https://tidyselect.r-lib.org/reference/starts_with.html)
  should be avoided unless the list of variable-value pairs is clearly
  specified in a statement via the `set_values_to` argument.

  Permitted values

  :   list of named expressions created by a formula using
      [`exprs()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/reexport-exprs.md),
      e.g., `exprs(AVALC = VSSTRESC, AVAL = yn_to_numeric(AVALC))`

  Default value

  :   none

- filter:

  Filter condition

  The specified condition is applied to the input dataset before
  deriving the new parameter, i.e., only observations fulfilling the
  condition are taken into account.

  Permitted values

  :   an unquoted condition, e.g., `AVISIT == "BASELINE"`

  Default value

  :   `NULL`

- constant_by_vars:

  By variables for constant parameters

  The constant parameters (parameters that are measured only once) are
  merged to the other parameters using the specified variables. (Refer
  to Example 2)

  Permitted values

  :   list of variables created by
      [`exprs()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/reexport-exprs.md),
      e.g., `exprs(USUBJID, VISIT)`

  Default value

  :   `NULL`

- constant_parameters:

  Required constant parameter codes

  It is expected that all the parameter codes (`PARAMCD`) which are
  required to derive the new parameter and are measured only once are
  specified here. For example if BMI should be derived and height is
  measured only once while weight is measured at each visit. Height
  could be specified in the `constant_parameters` parameter. (Refer to
  Example 2)

  If observations should be considered which do not have a parameter
  code, e.g., if an SDTM dataset is used, temporary parameter codes can
  be derived by specifying a list of expressions. The name of the
  element defines the temporary parameter code and the expression the
  condition for selecting the records. For example
  `constant_parameters = exprs(HGHT = VSTESTCD == "HEIGHT")` selects the
  observations with `VSTESTCD == "HEIGHT"` from the input data
  (`dataset` and `dataset_add`), sets `PARAMCD = "HGHT"` for these
  observations, and adds them to the observations to consider.

  Unnamed elements in the list of expressions are considered as
  parameter codes. For example,
  `constant_parameters = exprs(WEIGHT, HGHT = VSTESTCD == "HEIGHT")`
  uses the parameter code `"WEIGHT"` and creates a temporary parameter
  code `"HGHT"`.

  Permitted values

  :   A character vector of `PARAMCD` values or a list of expressions

  Default value

  :   `NULL`

- keep_nas:

  Keep observations with `NA`s

  If the argument is set to `TRUE`, observations are added even if some
  of the values contributing to the computed value are `NA` (see Example
  1b).

  If the argument is set to a list of variables, observations are added
  even if some of specified variables are `NA` (see Example 1c).

  Permitted values

  :   `TRUE`, `FALSE`, or a list of variables created by
      [`exprs()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/reexport-exprs.md)
      e.g. `exprs(ADTF, ATMF)`

  Default value

  :   `FALSE`

## Value

The input dataset with the new parameter added. Note, a variable will
only be populated in the new parameter rows if it is specified in
`by_vars`.

## Details

For each group (with respect to the variables specified for the
`by_vars` parameter) an observation is added to the output dataset if
the filtered input dataset (`dataset`) or the additional dataset
(`dataset_add`) contains exactly one observation for each parameter code
specified for `parameters` and all contributing values like `AVAL.SYSBP`
are not `NA`. The `keep_nas` can be used to specify variables for which
`NA`s are acceptable. See also Example 1b and 1c.

For the new observations the variables specified for `set_values_to` are
set to the provided values. The values of the other variables of the
input dataset are set to `NA`.

## See also

BDS-Findings Functions for adding Parameters/Records:
[`default_qtc_paramcd()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/default_qtc_paramcd.md),
[`derive_expected_records()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/derive_expected_records.md),
[`derive_extreme_event()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/derive_extreme_event.md),
[`derive_extreme_records()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/derive_extreme_records.md),
[`derive_locf_records()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/derive_locf_records.md),
[`derive_param_bmi()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/derive_param_bmi.md),
[`derive_param_bsa()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/derive_param_bsa.md),
[`derive_param_doseint()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/derive_param_doseint.md),
[`derive_param_exist_flag()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/derive_param_exist_flag.md),
[`derive_param_exposure()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/derive_param_exposure.md),
[`derive_param_framingham()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/derive_param_framingham.md),
[`derive_param_map()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/derive_param_map.md),
[`derive_param_qtc()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/derive_param_qtc.md),
[`derive_param_rr()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/derive_param_rr.md),
[`derive_param_wbc_abs()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/derive_param_wbc_abs.md),
[`derive_summary_records()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/derive_summary_records.md)

## Examples

### Example 1 - Data setup

Examples 1a, 1b, and 1c use the following `ADVS` data.

    ADVS <- tribble(
      ~USUBJID,      ~PARAMCD, ~PARAM,                            ~AVAL, ~VISIT,
      "01-701-1015", "DIABP",  "Diastolic Blood Pressure (mmHg)",    51, "BASELINE",
      "01-701-1015", "DIABP",  "Diastolic Blood Pressure (mmHg)",    50, "WEEK 2",
      "01-701-1015", "SYSBP",  "Systolic Blood Pressure (mmHg)",    121, "BASELINE",
      "01-701-1015", "SYSBP",  "Systolic Blood Pressure (mmHg)",    121, "WEEK 2",
      "01-701-1028", "DIABP",  "Diastolic Blood Pressure (mmHg)",    79, "BASELINE",
      "01-701-1028", "DIABP",  "Diastolic Blood Pressure (mmHg)",    80, "WEEK 2",
      "01-701-1028", "SYSBP",  "Systolic Blood Pressure (mmHg)",    130, "BASELINE",
      "01-701-1028", "SYSBP",  "Systolic Blood Pressure (mmHg)",     NA, "WEEK 2"
    ) %>%
      mutate(
        AVALU = "mmHg",
        ADT = case_when(
          VISIT == "BASELINE" ~ as.Date("2024-01-10"),
          VISIT == "WEEK 2" ~ as.Date("2024-01-24")
        ),
        ADTF = NA_character_
      )

### Example 1a - Adding a parameter computed from a formula (`parameters`, `set_values_to`)

Derive mean arterial pressure (MAP) from systolic (SYSBP) and diastolic
blood pressure (DIABP).

- Here, for each `USUBJID` and `VISIT` group (specified in `by_vars`),
  an observation is added to the output dataset when the filtered input
  dataset (`dataset`) contains exactly one observation for each
  parameter code specified for `parameters` and all contributing values
  (e.g., `AVAL.SYSBP` and `AVAL.DIABP`) are not `NA`. Indeed, patient
  `01-701-1028` does not get a `"WEEK 2"`-derived record as `AVAL` is
  `NA` for their `"WEEK 2"` systolic blood pressure.

    derive_param_computed(
      ADVS,
      by_vars = exprs(USUBJID, VISIT),
      parameters = c("SYSBP", "DIABP"),
      set_values_to = exprs(
        AVAL = (AVAL.SYSBP + 2 * AVAL.DIABP) / 3,
        PARAMCD = "MAP",
        PARAM = "Mean Arterial Pressure (mmHg)",
        AVALU = "mmHg",
        ADT = ADT.SYSBP
      )
    ) %>%
    select(-PARAM)
    #> # A tibble: 11 × 7
    #>    USUBJID     PARAMCD  AVAL VISIT    AVALU ADT        ADTF
    #>    <chr>       <chr>   <dbl> <chr>    <chr> <date>     <chr>
    #>  1 01-701-1015 DIABP    51   BASELINE mmHg  2024-01-10 <NA>
    #>  2 01-701-1015 DIABP    50   WEEK 2   mmHg  2024-01-24 <NA>
    #>  3 01-701-1015 SYSBP   121   BASELINE mmHg  2024-01-10 <NA>
    #>  4 01-701-1015 SYSBP   121   WEEK 2   mmHg  2024-01-24 <NA>
    #>  5 01-701-1028 DIABP    79   BASELINE mmHg  2024-01-10 <NA>
    #>  6 01-701-1028 DIABP    80   WEEK 2   mmHg  2024-01-24 <NA>
    #>  7 01-701-1028 SYSBP   130   BASELINE mmHg  2024-01-10 <NA>
    #>  8 01-701-1028 SYSBP    NA   WEEK 2   mmHg  2024-01-24 <NA>
    #>  9 01-701-1015 MAP      74.3 BASELINE mmHg  2024-01-10 <NA>
    #> 10 01-701-1015 MAP      73.7 WEEK 2   mmHg  2024-01-24 <NA>
    #> 11 01-701-1028 MAP      96   BASELINE mmHg  2024-01-10 <NA> 

### Example 1b - Keeping missing values for any source variables (`keep_nas = TRUE`)

Use option `keep_nas = TRUE` to derive MAP in the case where some/all
values of a variable used in the computation are missing.

- Note that observations will be added here even if some of the values
  contributing to the computed values are `NA`. In particular, patient
  `01-701-1028` does get a `"WEEK 2"`-derived record as compared to
  Example 1a, but with `AVAL = NA`.

    derive_param_computed(
      ADVS,
      by_vars = exprs(USUBJID, VISIT),
      parameters = c("SYSBP", "DIABP"),
      set_values_to = exprs(
        AVAL = (AVAL.SYSBP + 2 * AVAL.DIABP) / 3,
        PARAMCD = "MAP",
        PARAM = "Mean Arterial Pressure (mmHg)",
        AVALU = "mmHg",
        ADT = ADT.SYSBP,
        ADTF = ADTF.SYSBP
      ),
      keep_nas = TRUE
    )%>%
    select(-PARAM)
    #> # A tibble: 12 × 7
    #>    USUBJID     PARAMCD  AVAL VISIT    AVALU ADT        ADTF
    #>    <chr>       <chr>   <dbl> <chr>    <chr> <date>     <chr>
    #>  1 01-701-1015 DIABP    51   BASELINE mmHg  2024-01-10 <NA>
    #>  2 01-701-1015 DIABP    50   WEEK 2   mmHg  2024-01-24 <NA>
    #>  3 01-701-1015 SYSBP   121   BASELINE mmHg  2024-01-10 <NA>
    #>  4 01-701-1015 SYSBP   121   WEEK 2   mmHg  2024-01-24 <NA>
    #>  5 01-701-1028 DIABP    79   BASELINE mmHg  2024-01-10 <NA>
    #>  6 01-701-1028 DIABP    80   WEEK 2   mmHg  2024-01-24 <NA>
    #>  7 01-701-1028 SYSBP   130   BASELINE mmHg  2024-01-10 <NA>
    #>  8 01-701-1028 SYSBP    NA   WEEK 2   mmHg  2024-01-24 <NA>
    #>  9 01-701-1015 MAP      74.3 BASELINE mmHg  2024-01-10 <NA>
    #> 10 01-701-1015 MAP      73.7 WEEK 2   mmHg  2024-01-24 <NA>
    #> 11 01-701-1028 MAP      96   BASELINE mmHg  2024-01-10 <NA>
    #> 12 01-701-1028 MAP      NA   WEEK 2   mmHg  2024-01-24 <NA> 

### Example 1c - Keeping missing values for some source variables (`keep_nas = exprs()`)

Use option `keep_nas = exprs(ADTF)` to derive MAP in the case where
some/all values of a variable used in the computation are missing but
keeping `NA` values of `ADTF`.

- This is subtly distinct from Examples 1a and 1b. In 1a, we do not get
  new derived records if any of the source records have a value of `NA`
  for a variable that is included in `set_values_to`. In 1b, we do the
  opposite and allow the creation of new records regardless of how many
  `NA`s we encounter in the source variables.

- Here, we want to disregard `NA` values but only from the variables
  that are specified via `keep_na_values`.

- This is important because we have added `ADTF` in `set_values_to`, but
  all values of this variable are `NA`. As such, in order to get any
  derived records at all, but continue not getting one when `AVAL` is
  `NA` in any of the source records, (see patient `"01-701-1028"`
  again), we specify `keep_nas = exprs(ADTF)`.

    derive_param_computed(
      ADVS,
      by_vars = exprs(USUBJID, VISIT),
      parameters = c("SYSBP", "DIABP"),
      set_values_to = exprs(
        AVAL = (AVAL.SYSBP + 2 * AVAL.DIABP) / 3,
        PARAMCD = "MAP",
        PARAM = "Mean Arterial Pressure (mmHg)",
        AVALU = "mmHg",
        ADT = ADT.SYSBP,
        ADTF = ADTF.SYSBP
      ),
      keep_nas = exprs(ADTF)
    )
    #> # A tibble: 11 × 8
    #>    USUBJID     PARAMCD PARAM                   AVAL VISIT AVALU ADT        ADTF
    #>    <chr>       <chr>   <chr>                  <dbl> <chr> <chr> <date>     <chr>
    #>  1 01-701-1015 DIABP   Diastolic Blood Press…  51   BASE… mmHg  2024-01-10 <NA>
    #>  2 01-701-1015 DIABP   Diastolic Blood Press…  50   WEEK… mmHg  2024-01-24 <NA>
    #>  3 01-701-1015 SYSBP   Systolic Blood Pressu… 121   BASE… mmHg  2024-01-10 <NA>
    #>  4 01-701-1015 SYSBP   Systolic Blood Pressu… 121   WEEK… mmHg  2024-01-24 <NA>
    #>  5 01-701-1028 DIABP   Diastolic Blood Press…  79   BASE… mmHg  2024-01-10 <NA>
    #>  6 01-701-1028 DIABP   Diastolic Blood Press…  80   WEEK… mmHg  2024-01-24 <NA>
    #>  7 01-701-1028 SYSBP   Systolic Blood Pressu… 130   BASE… mmHg  2024-01-10 <NA>
    #>  8 01-701-1028 SYSBP   Systolic Blood Pressu…  NA   WEEK… mmHg  2024-01-24 <NA>
    #>  9 01-701-1015 MAP     Mean Arterial Pressur…  74.3 BASE… mmHg  2024-01-10 <NA>
    #> 10 01-701-1015 MAP     Mean Arterial Pressur…  73.7 WEEK… mmHg  2024-01-24 <NA>
    #> 11 01-701-1028 MAP     Mean Arterial Pressur…  96   BASE… mmHg  2024-01-10 <NA> 

### Example 2 - Derivations using parameters measured only once (`constant_parameters` and `constant_by_vars`)

Derive BMI where `HEIGHT` is measured only once.

- In the above examples, for each parameter specified in the
  `parameters` argument, we expect one record per by group, where the by
  group is specified in `by_vars`. However, if a parameter is only
  measured once, it can be specified in `constant_parameters` instead.

- A modified by group still needs to be provided for the constant
  parameters. This can be done via `constant_by_vars`.

- See the example below, where weight is measured for each patient at
  each visit (`by_vars = exprs(USUBJID, VISIT)`), while height is
  measured for each patient only at the first visit
  (`constant_parameters = "HEIGHT"`,
  `constant_by_vars = exprs(USUBJID`)).

    ADVS <- tribble(
      ~USUBJID,      ~PARAMCD, ~PARAM,        ~AVAL, ~AVALU, ~VISIT,
      "01-701-1015", "HEIGHT", "Height (cm)", 147.0, "cm",   "SCREENING",
      "01-701-1015", "WEIGHT", "Weight (kg)",  54.0, "kg",   "SCREENING",
      "01-701-1015", "WEIGHT", "Weight (kg)",  54.4, "kg",   "BASELINE",
      "01-701-1015", "WEIGHT", "Weight (kg)",  53.1, "kg",   "WEEK 2",
      "01-701-1028", "HEIGHT", "Height (cm)", 163.0, "cm",   "SCREENING",
      "01-701-1028", "WEIGHT", "Weight (kg)",  78.5, "kg",   "SCREENING",
      "01-701-1028", "WEIGHT", "Weight (kg)",  80.3, "kg",   "BASELINE",
      "01-701-1028", "WEIGHT", "Weight (kg)",  80.7, "kg",   "WEEK 2"
    )

    derive_param_computed(
      ADVS,
      by_vars = exprs(USUBJID, VISIT),
      parameters = "WEIGHT",
      set_values_to = exprs(
        AVAL = AVAL.WEIGHT / (AVAL.HEIGHT / 100)^2,
        PARAMCD = "BMI",
        PARAM = "Body Mass Index (kg/m^2)",
        AVALU = "kg/m^2"
      ),
      constant_parameters = c("HEIGHT"),
      constant_by_vars = exprs(USUBJID)
    )
    #> # A tibble: 14 × 6
    #>    USUBJID     PARAMCD PARAM                     AVAL AVALU  VISIT
    #>    <chr>       <chr>   <chr>                    <dbl> <chr>  <chr>
    #>  1 01-701-1015 HEIGHT  Height (cm)              147   cm     SCREENING
    #>  2 01-701-1015 WEIGHT  Weight (kg)               54   kg     SCREENING
    #>  3 01-701-1015 WEIGHT  Weight (kg)               54.4 kg     BASELINE
    #>  4 01-701-1015 WEIGHT  Weight (kg)               53.1 kg     WEEK 2
    #>  5 01-701-1028 HEIGHT  Height (cm)              163   cm     SCREENING
    #>  6 01-701-1028 WEIGHT  Weight (kg)               78.5 kg     SCREENING
    #>  7 01-701-1028 WEIGHT  Weight (kg)               80.3 kg     BASELINE
    #>  8 01-701-1028 WEIGHT  Weight (kg)               80.7 kg     WEEK 2
    #>  9 01-701-1015 BMI     Body Mass Index (kg/m^2)  25.0 kg/m^2 SCREENING
    #> 10 01-701-1015 BMI     Body Mass Index (kg/m^2)  25.2 kg/m^2 BASELINE
    #> 11 01-701-1015 BMI     Body Mass Index (kg/m^2)  24.6 kg/m^2 WEEK 2
    #> 12 01-701-1028 BMI     Body Mass Index (kg/m^2)  29.5 kg/m^2 SCREENING
    #> 13 01-701-1028 BMI     Body Mass Index (kg/m^2)  30.2 kg/m^2 BASELINE
    #> 14 01-701-1028 BMI     Body Mass Index (kg/m^2)  30.4 kg/m^2 WEEK 2   

### Example 3 - Derivations including data from an additional dataset (`dataset_add`) and non-`AVAL` variables

Use data from an additional dataset and other variables than `AVAL`.

- In this example, the dataset specified via `dataset_add` (e.g., `QS`)
  is an SDTM dataset. There is no parameter code in the dataset.

- The `parameters` argument is therefore used to specify a list of
  expressions to derive temporary parameter codes.

- Then, `set_values_to` is used to specify the values for the new
  observations of each variable, and variable-value pairs from both
  datasets are referenced via
  [`exprs()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/reexport-exprs.md).

    QS <- tribble(
      ~USUBJID, ~AVISIT,   ~QSTESTCD, ~QSORRES, ~QSSTRESN,
      "1",      "WEEK 2",  "CHSF112", NA,               1,
      "1",      "WEEK 2",  "CHSF113", "Yes",           NA,
      "1",      "WEEK 2",  "CHSF114", NA,               1,
      "1",      "WEEK 4",  "CHSF112", NA,               2,
      "1",      "WEEK 4",  "CHSF113", "No",            NA,
      "1",      "WEEK 4",  "CHSF114", NA,               1
    )

    ADCHSF <- tribble(
      ~USUBJID, ~AVISIT,  ~PARAMCD, ~QSSTRESN, ~AVAL,
      "1",      "WEEK 2", "CHSF12", 1,             6,
      "1",      "WEEK 2", "CHSF14", 1,             6,
      "1",      "WEEK 4", "CHSF12", 2,            12,
      "1",      "WEEK 4", "CHSF14", 1,             6
    ) %>%
      mutate(QSORRES = NA_character_)

    derive_param_computed(
      ADCHSF,
      dataset_add = QS,
      by_vars = exprs(USUBJID, AVISIT),
      parameters = exprs(CHSF12, CHSF13 = QSTESTCD %in% c("CHSF113"), CHSF14),
      set_values_to = exprs(
        AVAL = case_when(
          QSORRES.CHSF13 == "Not applicable" ~ 0,
          QSORRES.CHSF13 == "Yes" ~ 38,
          QSORRES.CHSF13 == "No" ~ if_else(
            QSSTRESN.CHSF12 > QSSTRESN.CHSF14,
            25,
            0
          )
        ),
        PARAMCD = "CHSF13"
      )
    )
    #> # A tibble: 6 × 6
    #>   USUBJID AVISIT PARAMCD QSSTRESN  AVAL QSORRES
    #>   <chr>   <chr>  <chr>      <dbl> <dbl> <chr>
    #> 1 1       WEEK 2 CHSF12         1     6 <NA>
    #> 2 1       WEEK 2 CHSF14         1     6 <NA>
    #> 3 1       WEEK 4 CHSF12         2    12 <NA>
    #> 4 1       WEEK 4 CHSF14         1     6 <NA>
    #> 5 1       WEEK 2 CHSF13        NA    38 <NA>
    #> 6 1       WEEK 4 CHSF13        NA    25 <NA>   

### Example 4 - Computing more than one variable

Specify more than one variable-value pair via `set_values_to`.

- In this example, the values of `AVALC`, `ADTM`, `ADTF`, `PARAMCD`, and
  `PARAM` are determined via distinctly defined analysis values and
  parameter codes.

- This is different from Example 3 as more than one variable is derived.

    ADLB_TBILIALK <- tribble(
      ~USUBJID, ~PARAMCD, ~AVALC, ~ADTM,        ~ADTF,
      "1",      "ALK2",   "Y",    "2021-05-13", NA_character_,
      "1",      "TBILI2", "Y",    "2021-06-30", "D",
      "2",      "ALK2",   "Y",    "2021-12-31", "M",
      "2",      "TBILI2", "N",    "2021-11-11", NA_character_,
      "3",      "ALK2",   "N",    "2021-04-03", NA_character_,
      "3",      "TBILI2", "N",    "2021-04-04", NA_character_
    ) %>%
      mutate(ADTM = ymd(ADTM))

    derive_param_computed(
      dataset_add = ADLB_TBILIALK,
      by_vars = exprs(USUBJID),
      parameters = c("ALK2", "TBILI2"),
      set_values_to = exprs(
        AVALC = if_else(AVALC.TBILI2 == "Y" & AVALC.ALK2 == "Y", "Y", "N"),
        ADTM = pmax(ADTM.TBILI2, ADTM.ALK2),
        ADTF = if_else(ADTM == ADTM.TBILI2, ADTF.TBILI2, ADTF.ALK2),
        PARAMCD = "TB2AK2",
        PARAM = "TBILI > 2 times ULN and ALKPH <= 2 times ULN"
      ),
      keep_nas = TRUE
    )
    #> # A tibble: 3 × 6
    #>   USUBJID AVALC ADTM       ADTF  PARAMCD PARAM
    #>   <chr>   <chr> <date>     <chr> <chr>   <chr>
    #> 1 1       Y     2021-06-30 D     TB2AK2  TBILI > 2 times ULN and ALKPH <= 2 tim…
    #> 2 2       N     2021-12-31 M     TB2AK2  TBILI > 2 times ULN and ALKPH <= 2 tim…
    #> 3 3       N     2021-04-04 <NA>  TB2AK2  TBILI > 2 times ULN and ALKPH <= 2 tim…
