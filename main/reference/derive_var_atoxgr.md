# Derive Lab High toxicity Grade 0 - 4 and Low Toxicity Grades 0 - (-4)

Derives character lab grade based on high and low severity/toxicity
grade(s).

## Usage

``` r
derive_var_atoxgr(
  dataset,
  lotox_description_var = ATOXDSCL,
  hitox_description_var = ATOXDSCH
)
```

## Arguments

- dataset:

  Input dataset

  The variables specified by the `lotox_description_var` and
  `hitox_description_var` arguments are expected to be in the dataset.
  `ATOXGRL`, and `ATOXGRH` are expected as well.

  Default value

  :   none

- lotox_description_var:

  Variable containing the toxicity grade description for low values, eg.
  "Anemia"

  Default value

  :   `ATOXDSCL`

- hitox_description_var:

  Variable containing the toxicity grade description for high values,
  eg. "Hemoglobin Increased".

  Default value

  :   `ATOXDSCH`

## Value

The input data set with the character variable added

## Details

Created variable `ATOXGR` will contain values "-4", "-3", "-2", "-1" for
low values and "1", "2", "3", "4" for high values, and will contain "0"
if value is gradable and does not satisfy any of the criteria for high
or low values. ATOXGR is set to missing if information not available to
give a grade.

Function applies the following rules:

- High and low missing - overall missing

- Low grade not missing and \> 0 - overall holds low grade

- High grade not missing and \> 0 - overall holds high grade

- (Only high direction OR low direction is NORMAL) and high grade
  normal - overall NORMAL

- (Only low direction OR high direction is NORMAL) and low grade
  normal - overall NORMAL

- otherwise set to missing

## See also

BDS-Findings Functions that returns variable appended to dataset:
[`derive_basetype_records()`](https:/pharmaverse.github.io/admiral/main/reference/derive_basetype_records.md),
[`derive_var_analysis_ratio()`](https:/pharmaverse.github.io/admiral/main/reference/derive_var_analysis_ratio.md),
[`derive_var_anrind()`](https:/pharmaverse.github.io/admiral/main/reference/derive_var_anrind.md),
[`derive_var_atoxgr_dir()`](https:/pharmaverse.github.io/admiral/main/reference/derive_var_atoxgr_dir.md),
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

adlb <- tribble(
  ~ATOXDSCL,          ~ATOXDSCH,        ~ATOXGRL,      ~ATOXGRH,
  "Hypoglycemia",     "Hyperglycemia",  NA_character_, "0",
  "Hypoglycemia",     "Hyperglycemia",  "0",           "1",
  "Hypoglycemia",     "Hyperglycemia",  "0",           "0",
  NA_character_,      "INR Increased",  NA_character_, "0",
  "Hypophosphatemia", NA_character_,    "1",           NA_character_
)

derive_var_atoxgr(adlb)
#> # A tibble: 5 × 5
#>   ATOXDSCL         ATOXDSCH      ATOXGRL ATOXGRH ATOXGR
#>   <chr>            <chr>         <chr>   <chr>   <chr> 
#> 1 Hypoglycemia     Hyperglycemia NA      0       NA    
#> 2 Hypoglycemia     Hyperglycemia 0       1       1     
#> 3 Hypoglycemia     Hyperglycemia 0       0       0     
#> 4 NA               INR Increased NA      0       0     
#> 5 Hypophosphatemia NA            1       NA      -1    
```
