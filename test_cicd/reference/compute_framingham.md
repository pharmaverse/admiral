# Compute Framingham Heart Study Cardiovascular Disease 10-Year Risk Score

Computes Framingham Heart Study Cardiovascular Disease 10-Year Risk
Score (FCVD101) based on systolic blood pressure, total serum
cholesterol (mg/dL), HDL serum cholesterol (mg/dL), sex, smoking status,
diabetic status, and treated for hypertension flag.

## Usage

``` r
compute_framingham(sysbp, chol, cholhdl, age, sex, smokefl, diabetfl, trthypfl)
```

## Arguments

- sysbp:

  Systolic blood pressure

  A numeric vector is expected.

  Default value

  :   none

- chol:

  Total serum cholesterol (mg/dL)

  A numeric vector is expected.

  Default value

  :   none

- cholhdl:

  HDL serum cholesterol (mg/dL)

  A numeric vector is expected.

  Default value

  :   none

- age:

  Age (years)

  A numeric vector is expected.

  Default value

  :   none

- sex:

  Gender

  A character vector is expected. Expected Values: 'M' 'F'

  Default value

  :   none

- smokefl:

  Smoking Status

  A character vector is expected. Expected Values: 'Y' 'N'

  Default value

  :   none

- diabetfl:

  Diabetic Status

  A character vector is expected. Expected Values: 'Y' 'N'

  Default value

  :   none

- trthypfl:

  Treated for hypertension status

  A character vector is expected. Expected Values: 'Y' 'N'

  Default value

  :   none

## Value

A numeric vector of Framingham values

## Details

The predicted probability of having cardiovascular disease (CVD) within
10-years according to Framingham formula. See AHA Journal article
General Cardiovascular Risk Profile for Use in Primary Care for
reference.

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

\$\$Risk = 100 \* (1 - RiskPeriodFactor ^ exp(RiskFactors))\$\$

## See also

[`derive_param_framingham()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/derive_param_framingham.md)

BDS-Findings Functions that returns a vector:
[`compute_bmi()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/compute_bmi.md),
[`compute_bsa()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/compute_bsa.md),
[`compute_egfr()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/compute_egfr.md),
[`compute_map()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/compute_map.md),
[`compute_qtc()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/compute_qtc.md),
[`compute_qual_imputation()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/compute_qual_imputation.md),
[`compute_qual_imputation_dec()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/compute_qual_imputation_dec.md),
[`compute_rr()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/compute_rr.md),
[`compute_scale()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/compute_scale.md),
[`transform_range()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/transform_range.md)

## Examples

``` r
compute_framingham(
  sysbp = 133, chol = 216.16, cholhdl = 54.91, age = 53,
  sex = "M", smokefl = "N", diabetfl = "N", trthypfl = "N"
)
#> [1] 10.37514

compute_framingham(
  sysbp = 161, chol = 186.39, cholhdl = 64.19, age = 52,
  sex = "F", smokefl = "Y", diabetfl = "N", trthypfl = "Y"
)
#> [1] 16.40353
```
