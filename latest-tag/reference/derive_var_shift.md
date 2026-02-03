# Derive Shift

Derives a character shift variable containing concatenated shift in
values based on user-defined pairing, e.g., shift from baseline to
analysis value, shift from baseline grade to analysis grade, ...

## Usage

``` r
derive_var_shift(
  dataset,
  new_var,
  from_var,
  to_var,
  missing_value = "NULL",
  sep_val = " to "
)
```

## Arguments

- dataset:

  Input dataset

  The variables specified by the `from_var` and `to_var` arguments are
  expected to be in the dataset.

  Default value

  :   none

- new_var:

  Name of the character shift variable to create.

  Default value

  :   none

- from_var:

  Variable containing value to shift from.

  Default value

  :   none

- to_var:

  Variable containing value to shift to.

  Default value

  :   none

- missing_value:

  Character string to replace missing values in `from_var` or `to_var`.

  Default value

  :   `"NULL"`

- sep_val:

  Character string to concatenate values of `from_var` and `to_var`.

  Default value

  :   `" to "`

## Value

The input dataset with the character shift variable added

## Details

`new_var` is derived by concatenating the values of `from_var` to values
of `to_var` (e.g. "NORMAL to HIGH"). When `from_var` or `to_var` has
missing value, the missing value is replaced by `missing_value` (e.g.
"NORMAL to NULL").

## See also

BDS-Findings Functions that returns variable appended to dataset:
[`derive_basetype_records()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/derive_basetype_records.md),
[`derive_var_analysis_ratio()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/derive_var_analysis_ratio.md),
[`derive_var_anrind()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/derive_var_anrind.md),
[`derive_var_atoxgr()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/derive_var_atoxgr.md),
[`derive_var_atoxgr_dir()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/derive_var_atoxgr_dir.md),
[`derive_var_base()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/derive_var_base.md),
[`derive_var_chg()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/derive_var_chg.md),
[`derive_var_nfrlt()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/derive_var_nfrlt.md),
[`derive_var_ontrtfl()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/derive_var_ontrtfl.md),
[`derive_var_pchg()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/derive_var_pchg.md),
[`derive_vars_crit_flag()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/derive_vars_crit_flag.md)

## Examples

``` r
library(tibble)

data <- tribble(
  ~USUBJID, ~PARAMCD, ~AVAL, ~ABLFL, ~BNRIND,  ~ANRIND,
  "P01",    "ALB",       33, "Y",    "LOW",    "LOW",
  "P01",    "ALB",       38, NA,     "LOW",    "NORMAL",
  "P01",    "ALB",       NA, NA,     "LOW",    NA,
  "P02",    "ALB",       37, "Y",    "NORMAL", "NORMAL",
  "P02",    "ALB",       49, NA,     "NORMAL", "HIGH",
  "P02",    "SODIUM",   147, "Y",    "HIGH",   "HIGH"
)

data %>%
  convert_blanks_to_na() %>%
  derive_var_shift(
    new_var = SHIFT1,
    from_var = BNRIND,
    to_var = ANRIND
  )
#> # A tibble: 6 × 7
#>   USUBJID PARAMCD  AVAL ABLFL BNRIND ANRIND SHIFT1          
#>   <chr>   <chr>   <dbl> <chr> <chr>  <chr>  <chr>           
#> 1 P01     ALB        33 Y     LOW    LOW    LOW to LOW      
#> 2 P01     ALB        38 NA    LOW    NORMAL LOW to NORMAL   
#> 3 P01     ALB        NA NA    LOW    NA     LOW to NULL     
#> 4 P02     ALB        37 Y     NORMAL NORMAL NORMAL to NORMAL
#> 5 P02     ALB        49 NA    NORMAL HIGH   NORMAL to HIGH  
#> 6 P02     SODIUM    147 Y     HIGH   HIGH   HIGH to HIGH    

# or only populate post-baseline records
data %>%
  convert_blanks_to_na() %>%
  restrict_derivation(
    derivation = derive_var_shift,
    args = params(
      new_var = SHIFT1,
      from_var = BNRIND,
      to_var = ANRIND
    ),
    filter = is.na(ABLFL)
  )
#> # A tibble: 6 × 7
#>   USUBJID PARAMCD  AVAL ABLFL BNRIND ANRIND SHIFT1        
#>   <chr>   <chr>   <dbl> <chr> <chr>  <chr>  <chr>         
#> 1 P01     ALB        38 NA    LOW    NORMAL LOW to NORMAL 
#> 2 P01     ALB        NA NA    LOW    NA     LOW to NULL   
#> 3 P02     ALB        49 NA    NORMAL HIGH   NORMAL to HIGH
#> 4 P01     ALB        33 Y     LOW    LOW    NA            
#> 5 P02     ALB        37 Y     NORMAL NORMAL NA            
#> 6 P02     SODIUM    147 Y     HIGH   HIGH   NA            
```
