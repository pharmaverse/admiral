test_that("Test 1: All Ratio Variables are Created", {
 expected_data <- tibble::tribble(
  ~USUBJID, ~PARAMCD, ~AVAL, ~BASE, ~ANRLO, ~ANRHI, ~R2BASE, ~R2ANRLO, ~R2ANRHI,
  "P01",    "ALT",    27,    27,    6,      34,     1,        4.5,      0.7941176,
  "P01",    "ALT",    41,    27,    6,      34,     1.5185185,6.8333333,1.2058824,
  "P01",    "ALT",    17,    27,    6,      34,     0.6296296,2.8333333,0.5000000,
  "P02",    "ALB",    38,    38,    33,     49,     1.0000000,1.1515152,0.7755102,
  "P02",    "ALB",    39,    38,    33,     49,     1.0263158,1.1818182,0.7959184,
  "P02",    "ALB",    37,    38,    33,     49,     0.9736842,1.1212121,0.7551020
  )

 data <-tibble::tribble(
   ~USUBJID, ~PARAMCD, ~AVAL, ~BASE, ~ANRLO, ~ANRHI,
  "P01",    "ALT",     27,    27,    6,      34,
  "P01",    "ALT",     41,    27,    6,      34,
  "P01",    "ALT",     17,    27,    6,      34,
  "P02",    "ALB",     38,    38,    33,     49,
  "P02",    "ALB",     39,    38,    33,     49,
  "P02",    "ALB",     37,    38,    33,     49
  )

 data %>% derive_vars_analysis_ratio()

})
