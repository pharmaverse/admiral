# derive_vars_computed ----
## Test 1: new variable is derived correctly ----
test_that("derive_vars_computed Test 1: new variable is derived correctly", {
  adsl <- tribble(
    ~STUDYID,   ~USUBJID, ~AGE,   ~AGEU,
    "PILOT01", "01-1302",   61, "YEARS",
    "PILOT01", "17-1344",   64, "YEARS"
  )

  advs <- tribble(
    ~STUDYID, ~USUBJID, ~PARAMCD, ~PARAM, ~VISIT, ~AVAL, ~AVALU, ~ABLFL,
    "PILOT01", "01-1302", "HEIGHT", "Height (cm)", "SCREENING", 177.8, "cm", "Y",
    "PILOT01", "01-1302", "WEIGHT", "Weight (kg)", "SCREENING", 81.19, "kg", "N",
    "PILOT01", "01-1302", "WEIGHT", "Weight (kg)", "BASELINE", 82.1, "kg", "Y",
    "PILOT01", "01-1302", "WEIGHT", "Weight (kg)", "WEEK 2", 81.19, "kg", "N",
    "PILOT01", "01-1302", "WEIGHT", "Weight (kg)", "WEEK 4", 82.56, "kg", "N",
    "PILOT01", "01-1302", "WEIGHT", "Weight (kg)", "WEEK 6", 80.74, "kg", "N",
    "PILOT01", "17-1344", "HEIGHT", "Height (cm)", "SCREENING", 163.5, "cm", "Y",
    "PILOT01", "17-1344", "WEIGHT", "Weight (kg)", "SCREENING", 58.06, "kg", "N",
    "PILOT01", "17-1344", "WEIGHT", "Weight (kg)", "BASELINE", 58.06, "kg", "Y",
    "PILOT01", "17-1344", "WEIGHT", "Weight (kg)", "WEEK 2", 58.97, "kg", "N",
    "PILOT01", "17-1344", "WEIGHT", "Weight (kg)", "WEEK 4", 57.97, "kg", "N",
    "PILOT01", "17-1344", "WEIGHT", "Weight (kg)", "WEEK 6", 58.97, "kg", "N"
  )

  expected_output <- adsl %>%
    mutate(BMIBL = case_when(
      USUBJID == "01-1302" ~ 25.9704601042,
      USUBJID == "17-1344" ~ 21.7190846262
    ))


  expect_dfs_equal(
    derive_vars_computed(
      dataset = adsl,
      dataset_add = advs,
      by_vars = exprs(STUDYID, USUBJID),
      parameters = c("WEIGHT"),
      constant_by_vars = exprs(STUDYID, USUBJID),
      constant_parameters = c("HEIGHT"),
      new_vars = exprs(BMIBL = compute_bmi(height = AVAL.HEIGHT, weight = AVAL.WEIGHT)),
      filter_add = ABLFL == "Y"
    ),
    expected_output,
    keys = c("USUBJID")
  )
})

## Test 2: no new variables added if filtered dataset is empty ----
test_that("derive_vars_computed Test 2: no new variables added if filtered dataset is empty", {
  adsl <- tribble(
    ~STUDYID,   ~USUBJID, ~AGE,   ~AGEU, ~SAFFL,
    "PILOT01", "01-1302",   61, "YEARS", "N",
    "PILOT01", "17-1344",   64, "YEARS", "N"
  )

  advs <- tribble(
    ~STUDYID, ~USUBJID, ~PARAMCD, ~PARAM, ~VISIT, ~AVAL, ~AVALU, ~ABLFL,
    "PILOT01", "01-1302", "HEIGHT", "Height (cm)", "SCREENING", 177.8, "cm", "Y",
    "PILOT01", "01-1302", "WEIGHT", "Weight (kg)", "SCREENING", 81.19, "kg", "N",
    "PILOT01", "01-1302", "WEIGHT", "Weight (kg)", "BASELINE", 82.1, "kg", "Y",
    "PILOT01", "01-1302", "WEIGHT", "Weight (kg)", "WEEK 2", 81.19, "kg", "N",
    "PILOT01", "01-1302", "WEIGHT", "Weight (kg)", "WEEK 4", 82.56, "kg", "N",
    "PILOT01", "01-1302", "WEIGHT", "Weight (kg)", "WEEK 6", 80.74, "kg", "N",
    "PILOT01", "17-1344", "HEIGHT", "Height (cm)", "SCREENING", 163.5, "cm", "Y",
    "PILOT01", "17-1344", "WEIGHT", "Weight (kg)", "SCREENING", 58.06, "kg", "N",
    "PILOT01", "17-1344", "WEIGHT", "Weight (kg)", "BASELINE", 58.06, "kg", "Y",
    "PILOT01", "17-1344", "WEIGHT", "Weight (kg)", "WEEK 2", 58.97, "kg", "N",
    "PILOT01", "17-1344", "WEIGHT", "Weight (kg)", "WEEK 4", 57.97, "kg", "N",
    "PILOT01", "17-1344", "WEIGHT", "Weight (kg)", "WEEK 6", 58.97, "kg", "N"
  )

  expected <- adsl

  expect_warning(
    derive_vars_computed(
      dataset = adsl,
      dataset_add = advs,
      by_vars = exprs(USUBJID),
      parameters = c("WEIGHT"),
      constant_by_vars = exprs(USUBJID),
      constant_parameters = c("HEIGHT"),
      new_vars = exprs(BMIBL = compute_bmi(height = AVAL.HEIGHT, weight = AVAL.WEIGHT)),
      filter_add = ABLFL == "Y",
      filter = SAFFL == "Y"
    ) %>%
      expect_dfs_equal(expected,
        keys = c("USUBJID")
      ),
    "The dataset does not contain any observations fullfiling the filter condition .*"
  )
})

## Test 3: no new variables are added if a parameter is missing ----
test_that("derive_vars_computed Test 3: no new variables are added if a parameter is missing", {
  adsl <- tribble(
    ~STUDYID,   ~USUBJID, ~AGE,   ~AGEU,
    "PILOT01", "01-1302",   61, "YEARS",
    "PILOT01", "17-1344",   64, "YEARS"
  )

  advs <- tribble(
    ~STUDYID, ~USUBJID, ~PARAMCD, ~PARAM, ~VISIT, ~AVAL, ~AVALU, ~ABLFL,
    "PILOT01", "01-1302", "HEIGHT", "Height (cm)", "SCREENING", 177.8, "cm", "Y",
    "PILOT01", "01-1302", "WEIGHT", "Weight (kg)", "SCREENING", 81.19, "kg", "N",
    "PILOT01", "01-1302", "WEIGHT", "Weight (kg)", "BASELINE", 82.1, "kg", "Y",
    "PILOT01", "01-1302", "WEIGHT", "Weight (kg)", "WEEK 2", 81.19, "kg", "N",
    "PILOT01", "01-1302", "WEIGHT", "Weight (kg)", "WEEK 4", 82.56, "kg", "N",
    "PILOT01", "01-1302", "WEIGHT", "Weight (kg)", "WEEK 6", 80.74, "kg", "N",
    "PILOT01", "17-1344", "HEIGHT", "Height (cm)", "SCREENING", 163.5, "cm", "Y",
    "PILOT01", "17-1344", "WEIGHT", "Weight (kg)", "SCREENING", 58.06, "kg", "N",
    "PILOT01", "17-1344", "WEIGHT", "Weight (kg)", "BASELINE", 58.06, "kg", "Y",
    "PILOT01", "17-1344", "WEIGHT", "Weight (kg)", "WEEK 2", 58.97, "kg", "N",
    "PILOT01", "17-1344", "WEIGHT", "Weight (kg)", "WEEK 4", 57.97, "kg", "N",
    "PILOT01", "17-1344", "WEIGHT", "Weight (kg)", "WEEK 6", 58.97, "kg", "N"
  )

  expected <- adsl

  expect_warning(
    derive_vars_computed(
      dataset = adsl,
      dataset_add = advs,
      by_vars = exprs(STUDYID, USUBJID),
      parameters = c("WEIGHT"),
      constant_by_vars = exprs(STUDYID, USUBJID),
      constant_parameters = c("HEIGHT"),
      new_vars = exprs(BMIBL = compute_bmi(height = AVAL.HEIGHT, weight = AVAL.WEIGHT)),
      filter_add = ABLFL == "Y" & PARAMCD == "HEIGHT"
    )
    %>%
      expect_dfs_equal(expected,
        keys = c("USUBJID")
      ),
    "The dataset does not contain any observations fullfiling the filter condition .*"
  )
})

# assert_parameters_argument ----
## Test 4: error if argument is of wrong type ----
test_that("assert_parameters_argument Test 4: error if argument is of wrong type", {
  expect_error(
    assert_parameters_argument(myparameters <- c(1, 2, 3)),
    regexp = paste(
      "`myparameters` must be a character vector or a list of expressions",
      "but it is a double vector."
    ),
    fixed = TRUE
  )
})


# get_hori_data ----
## Test 5: error if variables with more than one dot ----
test_that("get_temp_data Test 5: error if variables with more than one dot", {
  input <- tribble(
    ~STUDYID, ~USUBJID, ~PARAMCD, ~PARAM, ~VISIT, ~AVAL, ~AVALU,
    "PILOT01", "01-1302", "HEIGHT", "Height (cm)", "SCREENING", 160.4, "cm",
    "PILOT01", "01-1302", "HEIGHT", "Height (cm)", "BASELINE", 177.81, "cm",
    "PILOT01", "01-1302", "WEIGHT", "Weight (kg)", "SCREENING", 81.19, "kg",
    "PILOT01", "01-1302", "WEIGHT", "Weight (kg)", "BASELINE", 82.145, "kg",
    "PILOT01", "17-1344", "HEIGHT", "Height (cm)", "SCREENING", 150.5, "cm",
    "PILOT01", "17-1344", "HEIGHT", "Height (cm)", "BASELINE", 163.5, "cm",
    "PILOT01", "17-1344", "WEIGHT", "Weight (kg)", "SCREENING", 58.06, "kg",
    "PILOT01", "17-1344", "WEIGHT", "Weight (kg)", "BASELINE", 58.06, "kg"
  )

  expect_error(
    get_temp_data(
      input,
      parameters = exprs(HEIGHT, WEIGHT),
      by_vars = exprs(USUBJID, VISIT),
      set_values_to = exprs(BMIBL = compute_bmi(height = AVAL.HEIGHT, weight = AVAL.WEI.GHT)),
      filter = NULL
    ),
    regexp = paste(
      "The `new_vars` argument contains variable names with more than one dot:",
      "`AVAL.WEI.GHT`",
      sep = "\n"
    ),
    fixed = TRUE
  )
})
