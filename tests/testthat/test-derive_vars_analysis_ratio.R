test_that("Test 1: All Ratio Variables are Created", {
 expected_data <- tibble::tribble(
  ~USUBJID, ~PARAMCD, ~SEQ, ~AVAL, ~BASE, ~ANRLO, ~ANRHI, ~R2BASE, ~R2ANRLO, ~R2ANRHI,
  "P01",    "ALT",    1,    27,    27,    6,      34,     1,        4.5,     0.794,
  "P01",    "ALT",    2,    41,    27,    6,      34,     1.52,    6.83,    1.21,
  "P01",    "ALT",    3,    17,    27,    6,      34,     0.63,     2.83,    0.5,
  "P02",    "ALB",    1,    38,    38,    33,     49,     1,        1.15,    0.776,
  "P02",    "ALB",    2,    39,    38,    33,     49,     1.03,     1.18,    0.796,
  "P02",    "ALB",    3,    37,    38,    33,     49,     0.974,    1.12,    0.755
  )

 data <-tibble::tribble(
   ~USUBJID, ~PARAMCD, ~SEQ,  ~AVAL, ~BASE, ~ANRLO, ~ANRHI,
  "P01",    "ALT",     1,     27,    27,    6,      34,
  "P01",    "ALT",     2,     41,    27,    6,      34,
  "P01",    "ALT",     3,     17,    27,    6,      34,
  "P02",    "ALB",     1,     38,    38,    33,     49,
  "P02",    "ALB",     2,     39,    38,    33,     49,
  "P02",    "ALB",     3,     37,    38,    33,     49
  )

 actual_data <- data %>% derive_vars_analysis_ratio(ratio_vars = "high")

 expect_dfs_equal(expected_data,
                  actual_data,
                  keys = c("USUBJID", "PARAMCD", "SEQ"),
                  tolerance = 0.2)
})
