# Adds a Parameter for Framingham Heart Study Cardiovascular Disease 10-Year Risk Score

Adds a record for framingham score (FCVD101) for each by group (e.g.,
subject and visit) where the source parameters are available.

## Usage

``` r
derive_param_framingham(
  dataset,
  by_vars,
  set_values_to = exprs(PARAMCD = "FCVD101"),
  sysbp_code = "SYSBP",
  chol_code = "CHOL",
  cholhdl_code = "CHOLHDL",
  age = AGE,
  sex = SEX,
  smokefl = SMOKEFL,
  diabetfl = DIABETFL,
  trthypfl = TRTHYPFL,
  get_unit_expr,
  filter = NULL
)
```

## Arguments

- dataset:

  Input dataset

  The variables specified by the `by_vars` argument are expected to be
  in the dataset. `PARAMCD`, and `AVAL` are expected as well.

  The variable specified by `by_vars` and `PARAMCD` must be a unique key
  of the input dataset after restricting it by the filter condition
  (`filter` parameter) and to the parameters specified by `sysbp_code`,
  `chol_code` and `hdl_code`.

  Default value

  :   none

- by_vars:

  Grouping variables

  Only variables specified in `by_vars` will be populated in the newly
  created records.

  Permitted values

  :   list of variables created by
      [`exprs()`](https:/pharmaverse.github.io/admiral/main/reference/reexport-exprs.md),
      e.g., `exprs(USUBJID, VISIT)`

  Default value

  :   none

- set_values_to:

  Variables to be set

  The specified variables are set to the specified values for the new
  observations. For example `exprs(PARAMCD = "MAP")` defines the
  parameter code for the new parameter.

  Permitted values

  :   List of variable-value pairs

  Default value

  :   `exprs(PARAMCD = "MAP")`

- sysbp_code:

  Systolic blood pressure parameter code

  The observations where `PARAMCD` equals the specified value are
  considered as the systolic blood pressure assessments.

  Permitted values

  :   a character scalar, i.e., a character vector of length one

  Default value

  :   `"SYSBP"`

- chol_code:

  Total serum cholesterol code

  The observations where `PARAMCD` equals the specified value are
  considered as the total cholesterol assessments. This must be measured
  in mg/dL.

  Permitted values

  :   a character scalar, i.e., a character vector of length one

  Default value

  :   `"CHOL"`

- cholhdl_code:

  HDL serum cholesterol code

  The observations where `PARAMCD` equals the specified value are
  considered as the HDL cholesterol assessments. This must be measured
  in mg/dL.

  Permitted values

  :   a character scalar, i.e., a character vector of length one

  Default value

  :   `"CHOLHDL"`

- age:

  Subject age

  A variable containing the subject's age.

  Permitted values

  :   A numeric variable name that refers to a subject age column of the
      input dataset

  Default value

  :   `AGE`

- sex:

  Subject sex

  A variable containing the subject's sex.

  Permitted values

  :   A character variable name that refers to a subject sex column of
      the input dataset

  Default value

  :   `SEX`

- smokefl:

  Smoking status flag

  A flag indicating smoking status.

  Permitted values

  :   A character variable name that refers to a smoking status column
      of the input dataset.

  Default value

  :   `SMOKEFL`

- diabetfl:

  Diabetic flag

  A flag indicating diabetic status.

  Permitted values

  :   A character variable name that refers to a diabetic status column
      of the input dataset

  Default value

  :   `DIABETFL`

- trthypfl:

  Treated with hypertension medication flag

  A flag indicating if a subject was treated with hypertension
  medication.

  Permitted values

  :   A character variable name that refers to a column that indicates
      whether a subject is treated for high blood pressure

  Default value

  :   `TRTHYPFL`

- get_unit_expr:

  An expression providing the unit of the parameter

  The result is used to check the units of the input parameters.

  Permitted values

  :   An expression which is evaluable in the input dataset and results
      in a character value

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

## Value

The input dataset with the new parameter added

## Details

The values of `age`, `sex`, `smokefl`, `diabetfl` and `trthypfl` will be
added to the `by_vars` list. The predicted probability of having
cardiovascular disease (CVD) within 10-years according to Framingham
formula. See AHA Journal article General Cardiovascular Risk Profile for
Use in Primary Care for reference.

**For Women:**

|                            |            |
|----------------------------|------------|
| **Factor**                 | **Amount** |
| Age                        | 2.32888    |
| Total Chol                 | 1.20904    |
| HDL Chol                   | -0.70833   |
| Sys BP                     | 2.76157    |
| Sys BP + Hypertension Meds | 2.82263    |
| Smoker                     | 0.52873    |
| Non-Smoker                 | 0          |
| Diabetic                   | 0.69154    |
| Not Diabetic               | 0          |
| Average Risk               | 26.1931    |
| Risk Period                | 0.95012    |

**For Men:**

|                            |            |
|----------------------------|------------|
| **Factor**                 | **Amount** |
| Age                        | 3.06117    |
| Total Chol                 | 1.12370    |
| HDL Chol                   | -0.93263   |
| Sys BP                     | 1.93303    |
| Sys BP + Hypertension Meds | 2.99881    |
| Smoker                     | .65451     |
| Non-Smoker                 | 0          |
| Diabetic                   | 0.57367    |
| Not Diabetic               | 0          |
| Average Risk               | 23.9802    |
| Risk Period                | 0.88936    |

**The equation for calculating risk:**

\$\$RiskFactors = (log(Age) \* AgeFactor) + (log(TotalChol) \*
TotalCholFactor) + (log(CholHDL) \* CholHDLFactor) \\ + (log(SysBP) \*
SysBPFactor) + Smoker + Diabetes Present - AvgRisk\$\$

\$\$Risk = 100 \* (1 - RiskPeriodFactor^{RiskFactors})\$\$

## See also

[`compute_framingham()`](https:/pharmaverse.github.io/admiral/main/reference/compute_framingham.md)

BDS-Findings Functions for adding Parameters/Records:
[`default_qtc_paramcd()`](https:/pharmaverse.github.io/admiral/main/reference/default_qtc_paramcd.md),
[`derive_expected_records()`](https:/pharmaverse.github.io/admiral/main/reference/derive_expected_records.md),
[`derive_extreme_event()`](https:/pharmaverse.github.io/admiral/main/reference/derive_extreme_event.md),
[`derive_extreme_records()`](https:/pharmaverse.github.io/admiral/main/reference/derive_extreme_records.md),
[`derive_locf_records()`](https:/pharmaverse.github.io/admiral/main/reference/derive_locf_records.md),
[`derive_param_bmi()`](https:/pharmaverse.github.io/admiral/main/reference/derive_param_bmi.md),
[`derive_param_bsa()`](https:/pharmaverse.github.io/admiral/main/reference/derive_param_bsa.md),
[`derive_param_computed()`](https:/pharmaverse.github.io/admiral/main/reference/derive_param_computed.md),
[`derive_param_doseint()`](https:/pharmaverse.github.io/admiral/main/reference/derive_param_doseint.md),
[`derive_param_exist_flag()`](https:/pharmaverse.github.io/admiral/main/reference/derive_param_exist_flag.md),
[`derive_param_exposure()`](https:/pharmaverse.github.io/admiral/main/reference/derive_param_exposure.md),
[`derive_param_map()`](https:/pharmaverse.github.io/admiral/main/reference/derive_param_map.md),
[`derive_param_qtc()`](https:/pharmaverse.github.io/admiral/main/reference/derive_param_qtc.md),
[`derive_param_rr()`](https:/pharmaverse.github.io/admiral/main/reference/derive_param_rr.md),
[`derive_param_wbc_abs()`](https:/pharmaverse.github.io/admiral/main/reference/derive_param_wbc_abs.md),
[`derive_summary_records()`](https:/pharmaverse.github.io/admiral/main/reference/derive_summary_records.md)

## Examples

``` r
library(tibble)

adcvrisk <- tribble(
  ~USUBJID, ~PARAMCD, ~PARAM, ~AVAL, ~AVALU,
  ~VISIT, ~AGE, ~SEX, ~SMOKEFL, ~DIABETFL, ~TRTHYPFL,
  "01-701-1015", "SYSBP", "Systolic Blood Pressure (mmHg)", 121,
  "mmHg", "BASELINE", 44, "F", "N", "N", "N",
  "01-701-1015", "SYSBP", "Systolic Blood Pressure (mmHg)", 115,
  "mmHg", "WEEK 2", 44, "F", "N", "N", "Y",
  "01-701-1015", "CHOL", "Total Cholesterol (mg/dL)", 216.16,
  "mg/dL", "BASELINE", 44, "F", "N", "N", "N",
  "01-701-1015", "CHOL", "Total Cholesterol (mg/dL)", 210.78,
  "mg/dL", "WEEK 2", 44, "F", "N", "N", "Y",
  "01-701-1015", "CHOLHDL", "Cholesterol/HDL-Cholesterol (mg/dL)", 54.91,
  "mg/dL", "BASELINE", 44, "F", "N", "N", "N",
  "01-701-1015", "CHOLHDL", "Cholesterol/HDL-Cholesterol (mg/dL)", 26.72,
  "mg/dL", "WEEK 2", 44, "F", "N", "N", "Y",
  "01-701-1028", "SYSBP", "Systolic Blood Pressure (mmHg)", 119,
  "mmHg", "BASELINE", 55, "M", "Y", "Y", "Y",
  "01-701-1028", "SYSBP", "Systolic Blood Pressure (mmHg)", 101,
  "mmHg", "WEEK 2", 55, "M", "Y", "Y", "Y",
  "01-701-1028", "CHOL", "Total Cholesterol (mg/dL)", 292.01,
  "mg/dL", "BASELINE", 55, "M", "Y", "Y", "Y",
  "01-701-1028", "CHOL", "Total Cholesterol (mg/dL)", 246.73,
  "mg/dL", "WEEK 2", 55, "M", "Y", "Y", "Y",
  "01-701-1028", "CHOLHDL", "Cholesterol/HDL-Cholesterol (mg/dL)", 65.55,
  "mg/dL", "BASELINE", 55, "M", "Y", "Y", "Y",
  "01-701-1028", "CHOLHDL", "Cholesterol/HDL-Cholesterol (mg/dL)", 44.62,
  "mg/dL", "WEEK 2", 55, "M", "Y", "Y", "Y"
)


adcvrisk %>%
  derive_param_framingham(
    by_vars = exprs(USUBJID, VISIT),
    set_values_to = exprs(
      PARAMCD = "FCVD101",
      PARAM = "FCVD1-Framingham CVD 10-Year Risk Score (%)"
    ),
    get_unit_expr = AVALU
  )
#> # A tibble: 16 × 11
#>    USUBJID     PARAMCD PARAM       AVAL AVALU VISIT   AGE SEX   SMOKEFL DIABETFL
#>    <chr>       <chr>   <chr>      <dbl> <chr> <chr> <dbl> <chr> <chr>   <chr>   
#>  1 01-701-1015 SYSBP   Systolic… 121    mmHg  BASE…    44 F     N       N       
#>  2 01-701-1015 SYSBP   Systolic… 115    mmHg  WEEK…    44 F     N       N       
#>  3 01-701-1015 CHOL    Total Ch… 216.   mg/dL BASE…    44 F     N       N       
#>  4 01-701-1015 CHOL    Total Ch… 211.   mg/dL WEEK…    44 F     N       N       
#>  5 01-701-1015 CHOLHDL Choleste…  54.9  mg/dL BASE…    44 F     N       N       
#>  6 01-701-1015 CHOLHDL Choleste…  26.7  mg/dL WEEK…    44 F     N       N       
#>  7 01-701-1028 SYSBP   Systolic… 119    mmHg  BASE…    55 M     Y       Y       
#>  8 01-701-1028 SYSBP   Systolic… 101    mmHg  WEEK…    55 M     Y       Y       
#>  9 01-701-1028 CHOL    Total Ch… 292.   mg/dL BASE…    55 M     Y       Y       
#> 10 01-701-1028 CHOL    Total Ch… 247.   mg/dL WEEK…    55 M     Y       Y       
#> 11 01-701-1028 CHOLHDL Choleste…  65.6  mg/dL BASE…    55 M     Y       Y       
#> 12 01-701-1028 CHOLHDL Choleste…  44.6  mg/dL WEEK…    55 M     Y       Y       
#> 13 01-701-1015 FCVD101 FCVD1-Fr…   3.14 NA    BASE…    44 F     N       N       
#> 14 01-701-1015 FCVD101 FCVD1-Fr…   5.80 NA    WEEK…    44 F     N       N       
#> 15 01-701-1028 FCVD101 FCVD1-Fr…  42.3  NA    BASE…    55 M     Y       Y       
#> 16 01-701-1028 FCVD101 FCVD1-Fr…  37.5  NA    WEEK…    55 M     Y       Y       
#> # ℹ 1 more variable: TRTHYPFL <chr>

derive_param_framingham(
  adcvrisk,
  by_vars = exprs(USUBJID, VISIT),
  set_values_to = exprs(
    PARAMCD = "FCVD101",
    PARAM = "FCVD1-Framingham CVD 10-Year Risk Score (%)"
  ),
  get_unit_expr = extract_unit(PARAM)
)
#> # A tibble: 16 × 11
#>    USUBJID     PARAMCD PARAM       AVAL AVALU VISIT   AGE SEX   SMOKEFL DIABETFL
#>    <chr>       <chr>   <chr>      <dbl> <chr> <chr> <dbl> <chr> <chr>   <chr>   
#>  1 01-701-1015 SYSBP   Systolic… 121    mmHg  BASE…    44 F     N       N       
#>  2 01-701-1015 SYSBP   Systolic… 115    mmHg  WEEK…    44 F     N       N       
#>  3 01-701-1015 CHOL    Total Ch… 216.   mg/dL BASE…    44 F     N       N       
#>  4 01-701-1015 CHOL    Total Ch… 211.   mg/dL WEEK…    44 F     N       N       
#>  5 01-701-1015 CHOLHDL Choleste…  54.9  mg/dL BASE…    44 F     N       N       
#>  6 01-701-1015 CHOLHDL Choleste…  26.7  mg/dL WEEK…    44 F     N       N       
#>  7 01-701-1028 SYSBP   Systolic… 119    mmHg  BASE…    55 M     Y       Y       
#>  8 01-701-1028 SYSBP   Systolic… 101    mmHg  WEEK…    55 M     Y       Y       
#>  9 01-701-1028 CHOL    Total Ch… 292.   mg/dL BASE…    55 M     Y       Y       
#> 10 01-701-1028 CHOL    Total Ch… 247.   mg/dL WEEK…    55 M     Y       Y       
#> 11 01-701-1028 CHOLHDL Choleste…  65.6  mg/dL BASE…    55 M     Y       Y       
#> 12 01-701-1028 CHOLHDL Choleste…  44.6  mg/dL WEEK…    55 M     Y       Y       
#> 13 01-701-1015 FCVD101 FCVD1-Fr…   3.14 NA    BASE…    44 F     N       N       
#> 14 01-701-1015 FCVD101 FCVD1-Fr…   5.80 NA    WEEK…    44 F     N       N       
#> 15 01-701-1028 FCVD101 FCVD1-Fr…  42.3  NA    BASE…    55 M     Y       Y       
#> 16 01-701-1028 FCVD101 FCVD1-Fr…  37.5  NA    WEEK…    55 M     Y       Y       
#> # ℹ 1 more variable: TRTHYPFL <chr>
```
