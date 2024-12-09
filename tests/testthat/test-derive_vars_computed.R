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
    "PILOT01", "01-1302", "WEIGHT", "Weight (kg)", "SCREENING", 81.19, "kg", NA_character_,
    "PILOT01", "01-1302", "WEIGHT", "Weight (kg)", "BASELINE", 82.1, "kg", "Y",
    "PILOT01", "01-1302", "WEIGHT", "Weight (kg)", "WEEK 2", 81.19, "kg", NA_character_,
    "PILOT01", "01-1302", "WEIGHT", "Weight (kg)", "WEEK 4", 82.56, "kg", NA_character_,
    "PILOT01", "01-1302", "WEIGHT", "Weight (kg)", "WEEK 6", 80.74, "kg", NA_character_,
    "PILOT01", "17-1344", "HEIGHT", "Height (cm)", "SCREENING", 163.5, "cm", "Y",
    "PILOT01", "17-1344", "WEIGHT", "Weight (kg)", "SCREENING", 58.06, "kg", NA_character_,
    "PILOT01", "17-1344", "WEIGHT", "Weight (kg)", "BASELINE", 58.06, "kg", "Y",
    "PILOT01", "17-1344", "WEIGHT", "Weight (kg)", "WEEK 2", 58.97, "kg", NA_character_,
    "PILOT01", "17-1344", "WEIGHT", "Weight (kg)", "WEEK 4", 57.97, "kg", NA_character_,
    "PILOT01", "17-1344", "WEIGHT", "Weight (kg)", "WEEK 6", 58.97, "kg", NA_character_
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
    ~STUDYID,   ~USUBJID, ~AGE,   ~AGEU,
    "PILOT01", "01-1302",   61, "YEARS",
    "PILOT01", "17-1344",   64, "YEARS"
  )

  advs <- tribble(
    ~STUDYID, ~USUBJID, ~PARAMCD, ~PARAM, ~VISIT, ~AVAL, ~AVALU, ~ABLFL,
    "PILOT01", "01-1302", "HEIGHT", "Height (cm)", "SCREENING", 177.8, "cm", "Y",
    "PILOT01", "01-1302", "WEIGHT", "Weight (kg)", "SCREENING", 81.19, "kg", NA_character_,
    "PILOT01", "01-1302", "WEIGHT", "Weight (kg)", "BASELINE", 82.1, "kg", "Y",
    "PILOT01", "01-1302", "WEIGHT", "Weight (kg)", "WEEK 2", 81.19, "kg", NA_character_,
    "PILOT01", "01-1302", "WEIGHT", "Weight (kg)", "WEEK 4", 82.56, "kg", NA_character_,
    "PILOT01", "01-1302", "WEIGHT", "Weight (kg)", "WEEK 6", 80.74, "kg", NA_character_,
    "PILOT01", "17-1344", "HEIGHT", "Height (cm)", "SCREENING", 163.5, "cm", "Y",
    "PILOT01", "17-1344", "WEIGHT", "Weight (kg)", "SCREENING", 58.06, "kg", NA_character_,
    "PILOT01", "17-1344", "WEIGHT", "Weight (kg)", "BASELINE", 58.06, "kg", "Y",
    "PILOT01", "17-1344", "WEIGHT", "Weight (kg)", "WEEK 2", 58.97, "kg", NA_character_,
    "PILOT01", "17-1344", "WEIGHT", "Weight (kg)", "WEEK 4", 57.97, "kg", NA_character_,
    "PILOT01", "17-1344", "WEIGHT", "Weight (kg)", "WEEK 6", 58.97, "kg", NA_character_
  )

  expected <- adsl %>%
    mutate(BMIBL = NA_integer_)

  expect_warning(
    derive_vars_computed(
      dataset = adsl,
      dataset_add = advs,
      by_vars = exprs(USUBJID),
      parameters = c("WEIGHT"),
      constant_by_vars = exprs(USUBJID),
      constant_parameters = c("HEIGHT"),
      new_vars = exprs(BMIBL = compute_bmi(height = AVAL.HEIGHT, weight = AVAL.WEIGHT)),
      filter_add = ABLFL == ""
    ) %>%
      expect_dfs_equal(expected,
        keys = c("USUBJID")
      ),
    "The input dataset does not contain any observations fullfiling the filter condition .*"
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
    "PILOT01", "01-1302", "WEIGHT", "Weight (kg)", "SCREENING", 81.19, "kg", NA_character_,
    "PILOT01", "01-1302", "WEIGHT", "Weight (kg)", "BASELINE", 82.1, "kg", "Y",
    "PILOT01", "01-1302", "WEIGHT", "Weight (kg)", "WEEK 2", 81.19, "kg", NA_character_,
    "PILOT01", "01-1302", "WEIGHT", "Weight (kg)", "WEEK 4", 82.56, "kg", NA_character_,
    "PILOT01", "01-1302", "WEIGHT", "Weight (kg)", "WEEK 6", 80.74, "kg", NA_character_,
    "PILOT01", "17-1344", "WEIGHT", "Weight (kg)", "SCREENING", 58.06, "kg", NA_character_,
    "PILOT01", "17-1344", "WEIGHT", "Weight (kg)", "BASELINE", 58.06, "kg", "Y",
    "PILOT01", "17-1344", "WEIGHT", "Weight (kg)", "WEEK 2", 58.97, "kg", NA_character_,
    "PILOT01", "17-1344", "WEIGHT", "Weight (kg)", "WEEK 4", 57.97, "kg", NA_character_,
    "PILOT01", "17-1344", "WEIGHT", "Weight (kg)", "WEEK 6", 58.97, "kg", NA_character_
  )

  expected <- adsl %>%
    mutate(BMIBL = NA_integer_)

  expect_warning(
    derive_vars_computed(
      dataset = adsl,
      dataset_add = advs,
      by_vars = exprs(STUDYID, USUBJID),
      parameters = c("WEIGHT"),
      constant_by_vars = exprs(STUDYID, USUBJID),
      constant_parameters = c("HEIGHT"),
      new_vars = exprs(BMIBL = compute_bmi(height = AVAL.HEIGHT, weight = AVAL.WEIGHT)),
      filter_add = ABLFL == "Y"
    )
    %>%
      expect_dfs_equal(expected,
        keys = c("USUBJID")
      ),
    "The input dataset does not contain any observations fullfiling the filter condition .*"
  )
})
