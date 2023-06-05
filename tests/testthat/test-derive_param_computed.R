# derive_param_computed ----
## Test 1: new observations are derived correctly ----
test_that("derive_param_computed Test 1: new observations are derived correctly", {
  input <- tibble::tribble(
    ~USUBJID, ~PARAMCD, ~PARAM, ~AVAL, ~AVALU, ~VISIT,
    "01-701-1015", "DIABP", "Diastolic Blood Pressure (mmHg)", 51, "mmHg", "BASELINE",
    "01-701-1015", "DIABP", "Diastolic Blood Pressure (mmHg)", 50, "mmHg", "WEEK 2",
    "01-701-1015", "SYSBP", "Systolic Blood Pressure (mmHg)", 121, "mmHg", "BASELINE",
    "01-701-1015", "SYSBP", "Systolic Blood Pressure (mmHg)", 121, "mmHg", "WEEK 2",
    "01-701-1028", "DIABP", "Diastolic Blood Pressure (mmHg)", 79, "mmHg", "BASELINE",
    "01-701-1028", "DIABP", "Diastolic Blood Pressure (mmHg)", 80, "mmHg", "WEEK 2",
    "01-701-1028", "SYSBP", "Systolic Blood Pressure (mmHg)", 130, "mmHg", "BASELINE",
    "01-701-1028", "SYSBP", "Systolic Blood Pressure (mmHg)", 132, "mmHg", "WEEK 2"
  )

  new_obs <-
    inner_join(input %>% filter(PARAMCD == "DIABP") %>% select(USUBJID, VISIT, AVAL),
      input %>% filter(PARAMCD == "SYSBP") %>% select(USUBJID, VISIT, AVAL),
      by = c("USUBJID", "VISIT"),
      suffix = c(".DIABP", ".SYSBP")
    ) %>%
    mutate(
      AVAL = (2 * AVAL.DIABP + AVAL.SYSBP) / 3,
      PARAMCD = "MAP",
      PARAM = "Mean arterial pressure (mmHg)",
      AVALU = "mmHg"
    ) %>%
    select(-AVAL.DIABP, -AVAL.SYSBP)
  expected_output <- bind_rows(input, new_obs)

  expect_dfs_equal(
    derive_param_computed(
      input,
      parameters = exprs(SYSBP, DIABP),
      by_vars = exprs(USUBJID, VISIT),
      analysis_value = (AVAL.SYSBP + 2 * AVAL.DIABP) / 3,
      set_values_to = exprs(
        PARAMCD = "MAP",
        PARAM = "Mean arterial pressure (mmHg)",
        AVALU = "mmHg"
      )
    ),
    expected_output,
    keys = c("USUBJID", "PARAMCD", "VISIT")
  )
})

## Test 2: new observations with constant parameters ----
test_that("derive_param_computed Test 2: new observations with constant parameters", {
  input <- tibble::tribble(
    ~USUBJID,      ~PARAMCD, ~PARAM,        ~AVAL, ~AVALU, ~VISIT,
    "01-701-1015", "HEIGHT", "Height (cm)", 147.0, "cm",   "SCREENING",
    "01-701-1015", "WEIGHT", "Weight (kg)",  54.0, "kg",   "SCREENING",
    "01-701-1015", "WEIGHT", "Weight (kg)",  54.4, "kg",   "BASELINE",
    "01-701-1015", "WEIGHT", "Weight (kg)",  53.1, "kg",   "WEEK 2",
    "01-701-1028", "HEIGHT", "Height (cm)", 163.0, "cm",   "SCREENING",
    "01-701-1028", "WEIGHT", "Weight (kg)",  78.5, "kg",   "SCREENING",
    "01-701-1028", "WEIGHT", "Weight (kg)",  80.3, "kg",   "BASELINE",
    "01-701-1028", "WEIGHT", "Weight (kg)",  80.7, "kg",   "WEEK 2"
  )

  new_obs <-
    inner_join(input %>% filter(PARAMCD == "HEIGHT") %>% select(USUBJID, AVAL),
      input %>% filter(PARAMCD == "WEIGHT") %>% select(USUBJID, VISIT, AVAL),
      by = c("USUBJID"),
      suffix = c(".HEIGHT", ".WEIGHT")
    ) %>%
    mutate(
      AVAL = AVAL.WEIGHT / (AVAL.HEIGHT / 100)^2,
      PARAMCD = "BMI",
      PARAM = "Body Mass Index (kg/m2)",
      AVALU = "kg/m2"
    ) %>%
    select(-AVAL.HEIGHT, -AVAL.WEIGHT)
  expected_output <- bind_rows(input, new_obs)

  expect_dfs_equal(
    derive_param_computed(
      input,
      parameters = c("WEIGHT"),
      by_vars = exprs(USUBJID, VISIT),
      constant_parameters = c("HEIGHT"),
      constant_by_vars = exprs(USUBJID),
      analysis_value = AVAL.WEIGHT / (AVAL.HEIGHT / 100)^2,
      set_values_to = exprs(
        PARAMCD = "BMI",
        PARAM = "Body Mass Index (kg/m2)",
        AVALU = "kg/m2"
      )
    ),
    expected_output,
    keys = c("USUBJID", "PARAMCD", "VISIT")
  )
})

## Test 3: no new observations if filtered dataset is empty ----
test_that("derive_param_computed Test 3: no new observations if filtered dataset is empty", {
  input <- tibble::tribble(
    ~USUBJID, ~PARAMCD, ~PARAM, ~AVAL, ~AVALU, ~VISIT,
    "01-701-1015", "DIABP", "Diastolic Blood Pressure (mmHg)", 51, "mmHg", "BASELINE",
    "01-701-1015", "DIABP", "Diastolic Blood Pressure (mmHg)", 50, "mmHg", "WEEK 2",
    "01-701-1015", "SYSBP", "Systolic Blood Pressure (mmHg)", 121, "mmHg", "BASELINE",
    "01-701-1015", "SYSBP", "Systolic Blood Pressure (mmHg)", 121, "mmHg", "WEEK 2",
    "01-701-1028", "DIABP", "Diastolic Blood Pressure (mmHg)", 79, "mmHg", "BASELINE",
    "01-701-1028", "DIABP", "Diastolic Blood Pressure (mmHg)", 80, "mmHg", "WEEK 2",
    "01-701-1028", "SYSBP", "Systolic Blood Pressure (mmHg)", 130, "mmHg", "BASELINE",
    "01-701-1028", "SYSBP", "Systolic Blood Pressure (mmHg)", 132, "mmHg", "WEEK 2"
  )

  expect_warning(
    derive_param_computed(
      input,
      filter = VISIT == "WEEK 24",
      parameters = c("SYSBP", "DIABP"),
      by_vars = exprs(USUBJID, VISIT),
      analysis_value = (AVAL.SYSBP + 2 * AVAL.DIABP) / 3,
      set_values_to = exprs(
        PARAMCD = "MAP",
        PARAM = "Mean arterial pressure (mmHg)",
        AVALU = "mmHg"
      )
    ) %>%
      expect_dfs_equal(input,
        keys = c("USUBJID", "PARAMCD", "VISIT")
      ),
    "The input dataset does not contain any observations fullfiling the filter condition .*"
  )
})

## Test 4: no new observations are added if a parameter is missing ----
test_that("derive_param_computed Test 4: no new observations are added if a parameter is missing", {
  input <- tibble::tribble(
    ~USUBJID, ~PARAMCD, ~PARAM, ~AVAL, ~AVALU, ~VISIT,
    "01-701-1015", "DIABP", "Diastolic Blood Pressure (mmHg)", 51, "mmHg", "BASELINE",
    "01-701-1015", "DIABP", "Diastolic Blood Pressure (mmHg)", 50, "mmHg", "WEEK 2",
    "01-701-1015", "SYSBP", "Systolic Blood Pressure (mmHg)", 121, "mmHg", "BASELINE",
    "01-701-1015", "SYSBP", "Systolic Blood Pressure (mmHg)", 121, "mmHg", "WEEK 2",
    "01-701-1028", "DIABP", "Diastolic Blood Pressure (mmHg)", 79, "mmHg", "BASELINE",
    "01-701-1028", "DIABP", "Diastolic Blood Pressure (mmHg)", 80, "mmHg", "WEEK 2",
    "01-701-1028", "SYSBP", "Systolic Blood Pressure (mmHg)", 130, "mmHg", "BASELINE",
    "01-701-1028", "SYSBP", "Systolic Blood Pressure (mmHg)", 132, "mmHg", "WEEK 2"
  )

  expect_warning(
    derive_param_computed(
      input,
      filter = PARAMCD == "DIABP",
      parameters = exprs(SYSBP, DIABP),
      by_vars = exprs(USUBJID, VISIT),
      analysis_value = (AVAL.SYSBP + 2 * AVAL.DIABP) / 3,
      set_values_to = exprs(
        PARAMCD = "MAP",
        PARAM = "Mean arterial pressure (mmHg)",
        AVALU = "mmHg"
      )
    )
    %>%
      expect_dfs_equal(input,
        keys = c("USUBJID", "PARAMCD", "VISIT")
      ),
    "The input dataset does not contain any observations fullfiling the filter condition .*"
  )
})


## Test 5: `dataset_add`, creating new parameters ----
test_that("derive_param_computed Test 5: `dataset_add`, creating new parameters", {
  qs <- tibble::tribble(
    ~USUBJID, ~AVISIT,   ~QSTESTCD, ~QSORRES, ~QSSTRESN,
    "1",      "WEEK 2",  "CHSF112", NA,               1,
    "1",      "WEEK 2",  "CHSF113", "Yes",           NA,
    "1",      "WEEK 2",  "CHSF114", NA,               1,
    "1",      "WEEK 4",  "CHSF112", NA,               2,
    "1",      "WEEK 4",  "CHSF113", "No",            NA,
    "1",      "WEEK 4",  "CHSF114", NA,               1
  )

  adchsf <- tibble::tribble(
    ~USUBJID, ~AVISIT,  ~PARAMCD, ~QSORRES, ~QSSTRESN, ~AVAL,
    "1",      "WEEK 2", "CHSF12", NA,       1,             6,
    "1",      "WEEK 2", "CHSF14", NA,       1,             6,
    "1",      "WEEK 4", "CHSF12", NA,       2,            12,
    "1",      "WEEK 4", "CHSF14", NA,       1,             6
  )

  expected <- bind_rows(
    adchsf,
    tibble::tribble(
      ~USUBJID, ~AVISIT,  ~PARAMCD, ~AVAL,
      "1",      "WEEK 2", "CHSF13",    38,
      "1",      "WEEK 4", "CHSF13",    25
    )
  )

  expect_dfs_equal(
    base = expected,
    compare = derive_param_computed(
      adchsf,
      dataset_add = qs,
      by_vars = exprs(USUBJID, AVISIT),
      parameters = exprs(CHSF12, CHSF13 = QSTESTCD %in% c("CHSF113", "CHSF213"), CHSF14),
      analysis_value = case_when(
        QSORRES.CHSF13 == "Not applicable" ~ 0,
        QSORRES.CHSF13 == "Yes" ~ 38,
        QSORRES.CHSF13 == "No" ~ if_else(
          QSSTRESN.CHSF12 > QSSTRESN.CHSF14,
          25,
          0
        )
      ),
      set_values_to = exprs(PARAMCD = "CHSF13")
    ),
    keys = c("USUBJID", "PARAMCD", "AVISIT")
  )
})

## Test 6: no input dataset ----
test_that("derive_param_computed Test 6: no input dataset", {
  qs <- tibble::tribble(
    ~USUBJID, ~AVISIT,   ~QSTESTCD, ~QSORRES, ~QSSTRESN,
    "1",      "WEEK 2",  "CHSF112", NA,               1,
    "1",      "WEEK 2",  "CHSF113", "Yes",           NA,
    "1",      "WEEK 2",  "CHSF114", NA,               1,
    "1",      "WEEK 4",  "CHSF112", NA,               2,
    "1",      "WEEK 4",  "CHSF213", "No",            NA,
    "1",      "WEEK 4",  "CHSF114", NA,               1
  )

  expected <- tibble::tribble(
    ~USUBJID, ~AVISIT,  ~PARAMCD, ~AVAL,
    "1",      "WEEK 2", "CHSF13",    38,
    "1",      "WEEK 4", "CHSF13",    25
  )

  expect_dfs_equal(
    base = expected,
    compare = derive_param_computed(
      dataset_add = qs,
      by_vars = exprs(USUBJID, AVISIT),
      parameters = exprs(
        CHSF12 = QSTESTCD == "CHSF112",
        CHSF13 = QSTESTCD %in% c("CHSF113", "CHSF213"),
        CHSF14 = QSTESTCD == "CHSF114"
      ),
      analysis_value = case_when(
        QSORRES.CHSF13 == "Not applicable" ~ 0,
        QSORRES.CHSF13 == "Yes" ~ 38,
        QSORRES.CHSF13 == "No" ~ if_else(
          QSSTRESN.CHSF12 > QSSTRESN.CHSF14,
          25,
          0
        )
      ),
      set_values_to = exprs(PARAMCD = "CHSF13")
    ),
    keys = c("USUBJID", "PARAMCD", "AVISIT")
  )
})

## Test 7: expression in constant_parameters ----
test_that("derive_param_computed Test 7: expression in constant_parameters", {
  input <- tibble::tribble(
    ~USUBJID,      ~PARAMCD, ~PARAM,        ~AVAL, ~AVALU, ~VISIT,
    "01-701-1015", "WEIGHT", "Weight (kg)",  54.0, "kg",   "SCREENING",
    "01-701-1015", "WEIGHT", "Weight (kg)",  54.4, "kg",   "BASELINE",
    "01-701-1015", "WEIGHT", "Weight (kg)",  53.1, "kg",   "WEEK 2",
    "01-701-1028", "HEIGHT", "Height (cm)", 163.0, "cm",   "SCREENING",
    "01-701-1028", "WEIGHT", "Weight (kg)",  78.5, "kg",   "SCREENING",
    "01-701-1028", "WEIGHT", "Weight (kg)",  80.3, "kg",   "BASELINE",
    "01-701-1028", "WEIGHT", "Weight (kg)",  80.7, "kg",   "WEEK 2"
  )

  vs <- tibble::tribble(
    ~USUBJID,      ~VSTESTCD, ~VSTEST,  ~VSSTRESN, ~VSSTRESU,
    "01-701-1015", "HGHT",    "Height",     147.0, "cm"
  )

  new_obs <-
    inner_join(vs %>% filter(VSTESTCD == "HGHT") %>% select(USUBJID, AVAL = VSSTRESN),
      input %>% filter(PARAMCD == "WEIGHT") %>% select(USUBJID, VISIT, AVAL),
      by = c("USUBJID"),
      suffix = c(".HEIGHT", ".WEIGHT")
    ) %>%
    mutate(
      AVAL = AVAL.WEIGHT / (AVAL.HEIGHT / 100)^2,
      PARAMCD = "BMI",
      PARAM = "Body Mass Index (kg/m2)",
      AVALU = "kg/m2"
    ) %>%
    select(-AVAL.HEIGHT, -AVAL.WEIGHT)
  expected_output <- bind_rows(input, new_obs)

  expect_dfs_equal(
    derive_param_computed(
      input,
      dataset_add = vs,
      parameters = exprs(WEIGHT),
      by_vars = exprs(USUBJID, VISIT),
      constant_parameters = exprs("HEIGHT" = VSTESTCD == "HGHT"),
      constant_by_vars = exprs(USUBJID),
      analysis_value = AVAL.WEIGHT / (VSSTRESN.HEIGHT / 100)^2,
      set_values_to = exprs(
        PARAMCD = "BMI",
        PARAM = "Body Mass Index (kg/m2)",
        AVALU = "kg/m2"
      )
    ),
    expected_output,
    keys = c("USUBJID", "PARAMCD", "VISIT")
  )
})

## Test 8: no new observations if a constant parameter is missing ----
test_that("derive_param_computed Test 8: no new observations if a constant parameter is missing", {
  input <- tibble::tribble(
    ~USUBJID,      ~PARAMCD, ~PARAM,        ~AVAL, ~AVALU, ~VISIT,
    "01-701-1015", "WEIGHT", "Weight (kg)",  54.0, "kg",   "SCREENING",
    "01-701-1015", "WEIGHT", "Weight (kg)",  54.4, "kg",   "BASELINE",
    "01-701-1015", "WEIGHT", "Weight (kg)",  53.1, "kg",   "WEEK 2",
    "01-701-1028", "WEIGHT", "Weight (kg)",  78.5, "kg",   "SCREENING",
    "01-701-1028", "WEIGHT", "Weight (kg)",  80.3, "kg",   "BASELINE",
    "01-701-1028", "WEIGHT", "Weight (kg)",  80.7, "kg",   "WEEK 2"
  )

  expect_warning(
    output <- derive_param_computed(
      input,
      parameters = c("WEIGHT"),
      by_vars = exprs(USUBJID, VISIT),
      constant_parameters = c("HEIGHT"),
      constant_by_vars = exprs(USUBJID),
      analysis_value = AVAL.WEIGHT / (AVAL.HEIGHT / 100)^2,
      set_values_to = exprs(
        PARAMCD = "BMI",
        PARAM = "Body Mass Index (kg/m2)",
        AVALU = "kg/m2"
      )
    ),
    regexp = paste(
      paste(
        "The input dataset does not contain any observations fullfiling the filter",
        "condition (NULL) for the parameter codes (PARAMCD) `HEIGHT`"
      ),
      "No new observations were added.",
      sep = "\n"
    ),
    fixed = TRUE
  )

  expect_dfs_equal(
    output,
    input,
    keys = c("USUBJID", "PARAMCD", "VISIT")
  )
})

# assert_parameters_argument ----
## Test 9: error if argument is of wrong type ----
test_that("assert_parameters_argument Test 9: error if argument is of wrong type", {
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
## Test 10: error if variables with more than one dot ----
test_that("get_hori_data Test 10: error if variables with more than one dot", {
  input <- tibble::tribble(
    ~USUBJID, ~PARAMCD, ~PARAM, ~AVAL, ~AVALU, ~VISIT,
    "01-701-1015", "DIABP", "Diastolic Blood Pressure (mmHg)", 51, "mmHg", "BASELINE",
    "01-701-1015", "DIABP", "Diastolic Blood Pressure (mmHg)", 50, "mmHg", "WEEK 2",
    "01-701-1015", "SYSBP", "Systolic Blood Pressure (mmHg)", 121, "mmHg", "BASELINE",
    "01-701-1015", "SYSBP", "Systolic Blood Pressure (mmHg)", 121, "mmHg", "WEEK 2",
    "01-701-1028", "DIABP", "Diastolic Blood Pressure (mmHg)", 79, "mmHg", "BASELINE",
    "01-701-1028", "DIABP", "Diastolic Blood Pressure (mmHg)", 80, "mmHg", "WEEK 2",
    "01-701-1028", "SYSBP", "Systolic Blood Pressure (mmHg)", 130, "mmHg", "BASELINE",
    "01-701-1028", "SYSBP", "Systolic Blood Pressure (mmHg)", 132, "mmHg", "WEEK 2"
  )

  expect_error(
    get_hori_data(
      input,
      parameters = exprs(SYSBP, DIABP),
      by_vars = exprs(USUBJID, VISIT),
      analysis_value = (AVAL.SYSBP + 2 * AVAL.DIA.BP) / 3,
      filter = NULL
    ),
    regexp = paste(
      "The `analysis_value` argument contains variable names with more than on dot:",
      "`AVAL.DIA.BP`",
      sep = "\n"
    ),
    fixed = TRUE
  )
})
