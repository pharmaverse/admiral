test_that("BMI is calculated correctly", {
  # Expected values are taken from the Center of Disease Control and Prevention's
  # (CDC) 'Adult BMI Calculator' at
  # https://www.cdc.gov/healthyweight/assessing/bmi/adult_bmi/metric_bmi_calculator/bmi_calculator.html
  expect_equal(round(compute_bmi(180.3, 74.3), 1L), 22.9)
  expect_equal(round(compute_bmi(169, 51.3), 1L), 18)
  expect_equal(round(compute_bmi(175.9, 94.5), 1L), 30.5)
})

test_that("BMI parameter is correctly added to input dataset", {
  input <- tibble::tribble(
    ~USUBJID, ~PARAMCD, ~AVAL, ~VSSTRESU,
    "P01",    "HEIGHT", 180.3, "cm",
    "P01",    "WEIGHT",  74.3, "kg",
    "P02",    "HEIGHT", 169.0, "cm",
    "P02",    "WEIGHT",  51.3, "kg",
    "P03",    "HEIGHT", 175.9, "cm",
    "P03",    "WEIGHT",  94.5, "kg"
  )
  output <- derive_param_bmi(input, by_vars = vars(USUBJID), get_unit_expr = VSSTRESU)

  expect_true(nrow(output) == nrow(input) + 3L)
  expect_true(nrow(filter(output, PARAMCD == "BMI")) == 3L)
  expect_identical(
    output %>% filter(PARAMCD == "BMI") %>% pull(AVAL) %>% round(1L),
    c(22.9, 18, 30.5)
  )
})
