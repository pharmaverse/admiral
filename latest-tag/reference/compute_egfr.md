# Compute Estimated Glomerular Filtration Rate (eGFR) for Kidney Function

Compute Kidney Function Tests:

- Estimated Creatinine Clearance (CRCL) by Cockcroft-Gault equation

- Estimated Glomerular Filtration Rate (eGFR) by CKD-EPI or MDRD
  equations

## Usage

``` r
compute_egfr(creat, creatu = "SI", age, weight, sex, race = NULL, method)
```

## Arguments

- creat:

  Creatinine

  A numeric vector is expected.

  Default value

  :   none

- creatu:

  Creatinine Units

  A character vector is expected.

  Expected Values: `"SI"`, `"CV"`, `"umol/L"`, `"mg/dL"`

  Default value

  :   `"SI"`

- age:

  Age (years)

  A numeric vector is expected.

  Default value

  :   none

- weight:

  Weight (kg)

  A numeric vector is expected if `method = "CRCL"`

  Default value

  :   none

- sex:

  Gender

  A character vector is expected.

  Expected Values: `"M"`, `"F"`

  Default value

  :   none

- race:

  Race

  A character vector is expected if `method = "MDRD"`

  Expected Values: `"BLACK OR AFRICAN AMERICAN"` and others

  Default value

  :   `NULL`

- method:

  Method

  A character vector is expected.

  Expected Values: `"CRCL"`, `"CKD-EPI"`, `"MDRD"`

  Default value

  :   none

## Value

A numeric vector of egfr values

## Details

Calculates an estimate of Glomerular Filtration Rate (eGFR)

**CRCL Creatinine Clearance (Cockcroft-Gault)**

For Creatinine in umol/L:

\$\$\frac{(140 - age) \times weight(kg) \times
constant}{Serum\\Creatinine(\mu mol/L)}\$\$

\$\$Constant = 1.04\\for\\females, 1.23\\for\\males\$\$

For Creatinine in mg/dL:

\$\$\frac{(140 - age) \times weight(kg) \times (0.85\\if\\female)}{72
\times Serum\\Creatinine(mg/dL)}\$\$

units = mL/min

**CKD-EPI Chronic Kidney Disease Epidemiology Collaboration formula**

\$\$eGFR = 142 \times min(SCr/{\kappa}, 1)^{\alpha} \times
max(SCr/{\kappa}, 1)^{-1.200} \times 0.9938^{Age} \times 1.012
\[if\\female\]\$\$

SCr = standardized serum creatinine in mg/dL (Note SCr(mg/dL) =
Creat(umol/L) / 88.42)

\$\$\kappa\$\$ = 0.7 (females) or 0.9 (males) \$\$\alpha\$\$ = -0.241
(female) or -0.302 (male) units = mL/min/1.73 m2

**MDRD Modification of Diet in Renal Disease formula**

\$\$eGFR = 175 \times (SCr)^{-1.154} \times (age)^{-0.203} \times 0.742
\[if\\female\] \times 1.212 \[if\\Black\]\$\$

SCr = standardized serum creatinine in mg/dL (Note SCr(mg/dL) =
Creat(umol/L) / 88.42)

units = mL/min/1.73 m2

## See also

BDS-Findings Functions that returns a vector:
[`compute_bmi()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/compute_bmi.md),
[`compute_bsa()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/compute_bsa.md),
[`compute_framingham()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/compute_framingham.md),
[`compute_map()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/compute_map.md),
[`compute_qtc()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/compute_qtc.md),
[`compute_qual_imputation()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/compute_qual_imputation.md),
[`compute_qual_imputation_dec()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/compute_qual_imputation_dec.md),
[`compute_rr()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/compute_rr.md),
[`compute_scale()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/compute_scale.md),
[`transform_range()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/transform_range.md)

## Examples

``` r
compute_egfr(
  creat = 90, creatu = "umol/L", age = 53, weight = 85, sex = "M", method = "CRCL"
)
#> [1] 101.065

compute_egfr(
  creat = 90, creatu = "umol/L", age = 53, sex = "M", race = "ASIAN", method = "MDRD"
)
#> [1] 76.58319

compute_egfr(
  creat = 70, creatu = "umol/L", age = 52, sex = "F", race = "BLACK OR AFRICAN AMERICAN",
  method = "MDRD"
)
#> [1] 92.40002

compute_egfr(
  creat = 90, creatu = "umol/L", age = 53, sex = "M", method = "CKD-EPI"
)
#> [1] 88.10399


base <- tibble::tribble(
  ~STUDYID, ~USUBJID, ~AGE, ~SEX, ~RACE, ~WTBL, ~CREATBL, ~CREATBLU,
  "P01", "P01-1001", 55, "M", "WHITE", 90.7, 96.3, "umol/L",
  "P01", "P01-1002", 52, "F", "BLACK OR AFRICAN AMERICAN", 68.5, 70, "umol/L",
  "P01", "P01-1003", 67, "M", "BLACK OR AFRICAN AMERICAN", 85.0, 77, "umol/L",
  "P01", "P01-1004", 76, "F", "ASIAN", 60.7, 65, "umol/L",
)

base %>%
  dplyr::mutate(
    CRCL_CG = compute_egfr(
      creat = CREATBL, creatu = CREATBLU, age = AGE, weight = WTBL, sex = SEX,
      method = "CRCL"
    ),
    EGFR_EPI = compute_egfr(
      creat = CREATBL, creatu = CREATBLU, age = AGE, weight = WTBL, sex = SEX,
      method = "CKD-EPI"
    ),
    EGFR_MDRD = compute_egfr(
      creat = CREATBL, creatu = CREATBLU, age = AGE, weight = WTBL, sex = SEX,
      race = RACE, method = "MDRD"
    ),
  )
#> # A tibble: 4 × 11
#>   STUDYID USUBJID    AGE SEX   RACE       WTBL CREATBL CREATBLU CRCL_CG EGFR_EPI
#>   <chr>   <chr>    <dbl> <chr> <chr>     <dbl>   <dbl> <chr>      <dbl>    <dbl>
#> 1 P01     P01-1001    55 M     WHITE      90.7    96.3 umol/L      98.5     80.2
#> 2 P01     P01-1002    52 F     BLACK OR…  68.5    70   umol/L      89.6     89.7
#> 3 P01     P01-1003    67 M     BLACK OR…  85      77   umol/L      99.1     94.5
#> 4 P01     P01-1004    76 F     ASIAN      60.7    65   umol/L      62.2     84.5
#> # ℹ 1 more variable: EGFR_MDRD <dbl>
```
