# Country Code Lookup

These pre-defined country codes are sourced from [ISO 3166
Standards](https://www.iso.org/iso-3166-country-codes.html). See also
[Wikipedia](https://en.wikipedia.org/wiki/ISO_3166-1_alpha-3).

## Usage

``` r
country_code_lookup
```

## Format

An object of class `tbl_df` (inherits from `tbl`, `data.frame`) with 249
rows and 3 columns.

## Details

`country_code` is the 3-letter ISO 3166-1 county code commonly found in
the ADSL `COUNTRY` variable. `country_name` is the country long name
corresponding to to the 3-letter code. `country_number` is the numeric
code corresponding to an alphabetic sorting of the 3-letter codes.

To see the entire table in the console, run
`print(country_code_lookup)`.

## See also

[dose_freq_lookup](https:/pharmaverse.github.io/admiral/main/reference/dose_freq_lookup.md)

Other metadata:
[`atoxgr_criteria_ctcv4`](https:/pharmaverse.github.io/admiral/main/reference/atoxgr_criteria_ctcv4.md),
[`atoxgr_criteria_ctcv4_uscv`](https:/pharmaverse.github.io/admiral/main/reference/atoxgr_criteria_ctcv4_uscv.md),
[`atoxgr_criteria_ctcv5`](https:/pharmaverse.github.io/admiral/main/reference/atoxgr_criteria_ctcv5.md),
[`atoxgr_criteria_ctcv5_uscv`](https:/pharmaverse.github.io/admiral/main/reference/atoxgr_criteria_ctcv5_uscv.md),
[`atoxgr_criteria_ctcv6`](https:/pharmaverse.github.io/admiral/main/reference/atoxgr_criteria_ctcv6.md),
[`atoxgr_criteria_ctcv6_uscv`](https:/pharmaverse.github.io/admiral/main/reference/atoxgr_criteria_ctcv6_uscv.md),
[`atoxgr_criteria_daids`](https:/pharmaverse.github.io/admiral/main/reference/atoxgr_criteria_daids.md),
[`atoxgr_criteria_daids_uscv`](https:/pharmaverse.github.io/admiral/main/reference/atoxgr_criteria_daids_uscv.md),
[`dose_freq_lookup`](https:/pharmaverse.github.io/admiral/main/reference/dose_freq_lookup.md)

## Examples

``` r
library(tibble)
library(dplyr, warn.conflicts = FALSE)

# Create reference dataset for periods
adsl <- tribble(
  ~USUBJID, ~SEX, ~COUNTRY,
  "ST01-01", "F", "AUT",
  "ST01-02", "M", "MWI",
  "ST01-03", "F", "GBR",
  "ST01-04", "M", "CHE",
  "ST01-05", "M", "NOR",
  "ST01-06", "F", "JPN",
  "ST01-07", "F", "USA"
)

adsl %>%
  derive_vars_merged(
    dataset_add = country_code_lookup,
    new_vars = exprs(COUNTRYN = country_number, COUNTRYL = country_name),
    by_vars = exprs(COUNTRY = country_code)
  )
#> # A tibble: 7 × 5
#>   USUBJID SEX   COUNTRY COUNTRYN COUNTRYL                                       
#>   <chr>   <chr> <chr>      <dbl> <chr>                                          
#> 1 ST01-01 F     AUT           16 Austria                                        
#> 2 ST01-02 M     MWI          157 Malawi                                         
#> 3 ST01-03 F     GBR           80 United Kingdom of Great Britain and Northern I…
#> 4 ST01-04 M     CHE           42 Switzerland                                    
#> 5 ST01-05 M     NOR          168 Norway                                         
#> 6 ST01-06 F     JPN          116 Japan                                          
#> 7 ST01-07 F     USA          235 United States of America                       
```
