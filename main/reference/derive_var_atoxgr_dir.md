# Derive Lab Toxicity Grade 0 - 4

Derives a character lab grade based on severity/toxicity criteria.

## Usage

``` r
derive_var_atoxgr_dir(
  dataset,
  new_var,
  tox_description_var,
  meta_criteria,
  criteria_direction,
  abnormal_indicator = NULL,
  high_indicator = NULL,
  low_indicator = NULL,
  get_unit_expr,
  signif_dig = get_admiral_option("signif_digits")
)
```

## Arguments

- dataset:

  Input dataset

  The variables specified by the `tox_description_var` argument are
  expected to be in the dataset.

  Default value

  :   none

- new_var:

  Name of the character grade variable to create, for example, `ATOXGRH`
  or `ATOXGRL`.

  Default value

  :   none

- tox_description_var:

  Variable containing the description of the grading criteria. For
  example: "Anemia" or "INR Increased".

  Default value

  :   none

- meta_criteria:

  Metadata data set holding the criteria (normally a case statement)

  Permitted values

  :   `atoxgr_criteria_ctcv4`, `atoxgr_criteria_ctcv5`,
      `atoxgr_criteria_ctcv6`, `atoxgr_criteria_daids`

      - `atoxgr_criteria_ctcv4` implements [Common Terminology Criteria
        for Adverse Events (CTCAE)
        v4.0](https://dctd.cancer.gov/research/ctep-trials/trial-development#ctcae-and-ctep-codes)

      - `atoxgr_criteria_ctcv5` implements [Common Terminology Criteria
        for Adverse Events (CTCAE)
        v5.0](https://dctd.cancer.gov/research/ctep-trials/for-sites/adverse-events#ctep-ctcae)

      - `atoxgr_criteria_ctcv6` implements [Common Terminology Criteria
        for Adverse Events (CTCAE)
        v6.0](https://dctd.cancer.gov/research/ctep-trials/for-sites/adverse-events#ctep-ctcae)

      - `atoxgr_criteria_daids` implements [Division of AIDS (DAIDS)
        Table for Grading the Severity of Adult and Pediatric Adverse
        Events](https://rsc.niaid.nih.gov/sites/default/files/daidsgradingcorrectedv21.pdf)

      The metadata should have the following variables:

      - `TERM`: variable to hold the term describing the criteria
        applied to a particular lab test, eg. "Anemia" or "INR
        Increased". Note: the variable is case insensitive.

      - `DIRECTION`: variable to hold the direction of the abnormality
        of a particular lab test value. "L" is for LOW values, "H" is
        for HIGH values. Note: the variable is case insensitive.

      - `UNIT_CHECK`: variable to hold unit of particular lab test. Used
        to check against input data if criteria is based on absolute
        values.

      - `VAR_CHECK`: variable to hold comma separated list of variables
        used in criteria. Used to check against input data that
        variables exist.

      - `GRADE_CRITERIA_CODE`: variable to hold code that creates grade
        based on defined criteria.

      - `FILTER`: Required only for DAIDS grading, specifies `admiral`
        code to filter the lab data based on a subset of subjects (e.g.
        AGE \> 18 YEARS)

  Default value

  :   none

- criteria_direction:

  Direction (L= Low, H = High) of toxicity grade.

  Permitted values

  :   "L", "H"

  Default value

  :   none

- abnormal_indicator:

  **\[deprecated\]** Please use `low_indicator` and `high_indicator`
  instead.

  Default value

  :   `NULL`

- high_indicator:

  Value in `BNRIND` derivation to indicate an abnormal high value.
  Usually "HIGH" for `criteria_direction` = "H".

  This is only required when `meta_criteria = atoxgr_criteria_ctcv5` or
  `meta_criteria = atoxgr_criteria_ctcv6` and `BNRIND` is a required
  variable. Currently, for terms `"Alanine aminotransferase increased"`,
  `"Aspartate aminotransferase increased"`,
  `"Blood bilirubin increased"` and `"GGT increased"` for both sets of
  criteria. Also, term `"Alkaline phosphatase increased"` for
  `meta_criteria = atoxgr_criteria_ctcv5`.

  Default value

  :   `NULL`

- low_indicator:

  Value in `BNRIND` derivation to indicate an abnormal low value.
  Usually "LOW" for `criteria_direction` = "L".

  This is only required when `meta_criteria = atoxgr_criteria_ctcv6` and
  `BNRIND` is a required variable. Currently, only for term
  `"Creatinine increased"`.

  Default value

  :   `NULL`

- get_unit_expr:

  An expression providing the unit of the parameter

  The result is used to check the units of the input parameters.
  Compared with `UNIT_CHECK` in metadata (see `meta_criteria`
  parameter).

  Permitted values

  :   A variable containing unit from the input dataset, or a function
      call, for example, `get_unit_expr = extract_unit(PARAM)`.

  Default value

  :   none

- signif_dig:

  Number of significant digits to use when comparing a lab value against
  another value.

  Significant digits used to avoid floating point discrepancies when
  comparing numeric values. See blog: [How admiral handles floating
  points](https://pharmaverse.github.io/blog/posts/2023-10-30_floating_point/floating_point.html)

  Default value

  :   `get_admiral_option("signif_digits")`

## Value

The input dataset with the character variable added

## Details

`new_var` is derived with values NA, "0", "1", "2", "3", "4", where "4"
is the most severe grade

- "4" is where the lab value satisfies the criteria for grade 4.

- "3" is where the lab value satisfies the criteria for grade 3.

- "2" is where the lab value satisfies the criteria for grade 2.

- "1" is where the lab value satisfies the criteria for grade 1.

- "0" is where a grade can be derived and is not grade "1", "2", "3" or
  "4".

- NA is where a grade cannot be derived.

## See also

BDS-Findings Functions that returns variable appended to dataset:
[`derive_basetype_records()`](https:/pharmaverse.github.io/admiral/main/reference/derive_basetype_records.md),
[`derive_var_analysis_ratio()`](https:/pharmaverse.github.io/admiral/main/reference/derive_var_analysis_ratio.md),
[`derive_var_anrind()`](https:/pharmaverse.github.io/admiral/main/reference/derive_var_anrind.md),
[`derive_var_atoxgr()`](https:/pharmaverse.github.io/admiral/main/reference/derive_var_atoxgr.md),
[`derive_var_base()`](https:/pharmaverse.github.io/admiral/main/reference/derive_var_base.md),
[`derive_var_chg()`](https:/pharmaverse.github.io/admiral/main/reference/derive_var_chg.md),
[`derive_var_nfrlt()`](https:/pharmaverse.github.io/admiral/main/reference/derive_var_nfrlt.md),
[`derive_var_ontrtfl()`](https:/pharmaverse.github.io/admiral/main/reference/derive_var_ontrtfl.md),
[`derive_var_pchg()`](https:/pharmaverse.github.io/admiral/main/reference/derive_var_pchg.md),
[`derive_var_shift()`](https:/pharmaverse.github.io/admiral/main/reference/derive_var_shift.md),
[`derive_vars_crit_flag()`](https:/pharmaverse.github.io/admiral/main/reference/derive_vars_crit_flag.md)

## Examples

``` r
library(tibble)

data <- tribble(
  ~ATOXDSCL,                    ~AVAL, ~ANRLO, ~ANRHI, ~PARAM,
  "Hypoglycemia",               119,   4,      7,      "Glucose (mmol/L)",
  "Lymphocyte count decreased", 0.7,   1,      4,      "Lymphocytes Abs (10^9/L)",
  "Anemia",                     129,   120,    180,    "Hemoglobin (g/L)",
  "White blood cell decreased", 10,    5,      20,     "White blood cell (10^9/L)",
  "White blood cell decreased", 15,    5,      20,     "White blood cell (10^9/L)",
  "Anemia",                     140,   120,    180,    "Hemoglobin (g/L)"
)

derive_var_atoxgr_dir(data,
  new_var = ATOXGRL,
  tox_description_var = ATOXDSCL,
  meta_criteria = atoxgr_criteria_ctcv5,
  criteria_direction = "L",
  get_unit_expr = extract_unit(PARAM)
)
#> # A tibble: 6 × 6
#>   ATOXDSCL                    AVAL ANRLO ANRHI PARAM                     ATOXGRL
#>   <chr>                      <dbl> <dbl> <dbl> <chr>                     <chr>  
#> 1 Anemia                     129     120   180 Hemoglobin (g/L)          0      
#> 2 Anemia                     140     120   180 Hemoglobin (g/L)          0      
#> 3 Hypoglycemia               119       4     7 Glucose (mmol/L)          0      
#> 4 Lymphocyte count decreased   0.7     1     4 Lymphocytes Abs (10^9/L)  2      
#> 5 White blood cell decreased  10       5    20 White blood cell (10^9/L) 0      
#> 6 White blood cell decreased  15       5    20 White blood cell (10^9/L) 0      

data <- tribble(
  ~ATOXDSCH,                     ~AVAL,  ~ANRLO,   ~ANRHI, ~PARAM,
  "CPK increased",               129,    0,        30,     "Creatine Kinase (U/L)",
  "Lymphocyte count increased",  4,      1,        4,      "Lymphocytes Abs (10^9/L)",
  "Lymphocyte count increased",  2,      1,        4,      "Lymphocytes Abs (10^9/L)",
  "CPK increased",               140,    120,      180,    "Creatine Kinase (U/L)"
)

derive_var_atoxgr_dir(data,
  new_var = ATOXGRH,
  tox_description_var = ATOXDSCH,
  meta_criteria = atoxgr_criteria_ctcv5,
  criteria_direction = "H",
  get_unit_expr = extract_unit(PARAM)
)
#> # A tibble: 4 × 6
#>   ATOXDSCH                    AVAL ANRLO ANRHI PARAM                    ATOXGRH
#>   <chr>                      <dbl> <dbl> <dbl> <chr>                    <chr>  
#> 1 CPK increased                129     0    30 Creatine Kinase (U/L)    2      
#> 2 CPK increased                140   120   180 Creatine Kinase (U/L)    0      
#> 3 Lymphocyte count increased     4     1     4 Lymphocytes Abs (10^9/L) 0      
#> 4 Lymphocyte count increased     2     1     4 Lymphocytes Abs (10^9/L) 0      
```
