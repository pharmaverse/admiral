# Get Default Parameter Code for Corrected QT

Get Default Parameter Code for Corrected QT

## Usage

``` r
default_qtc_paramcd(method)
```

## Arguments

- method:

  Method used to QT correction

  Permitted values

  :   `"Bazett"`, `"Fridericia"`, `"Sagie"`

  Default value

  :   none

## Value

`"QTCBR"` if `method` is `"Bazett"`, `"QTCFR"` if it's `"Fridericia"` or
`"QTLCR"` if it's `"Sagie"`. An error otherwise.

## See also

[`derive_param_qtc()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/derive_param_qtc.md)

BDS-Findings Functions for adding Parameters/Records:
[`derive_expected_records()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/derive_expected_records.md),
[`derive_extreme_event()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/derive_extreme_event.md),
[`derive_extreme_records()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/derive_extreme_records.md),
[`derive_locf_records()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/derive_locf_records.md),
[`derive_param_bmi()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/derive_param_bmi.md),
[`derive_param_bsa()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/derive_param_bsa.md),
[`derive_param_computed()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/derive_param_computed.md),
[`derive_param_doseint()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/derive_param_doseint.md),
[`derive_param_exist_flag()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/derive_param_exist_flag.md),
[`derive_param_exposure()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/derive_param_exposure.md),
[`derive_param_framingham()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/derive_param_framingham.md),
[`derive_param_map()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/derive_param_map.md),
[`derive_param_qtc()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/derive_param_qtc.md),
[`derive_param_rr()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/derive_param_rr.md),
[`derive_param_wbc_abs()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/derive_param_wbc_abs.md),
[`derive_summary_records()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/derive_summary_records.md)

## Examples

``` r
default_qtc_paramcd("Sagie")
#> $PARAMCD
#> [1] "QTLCR"
#> 
```
