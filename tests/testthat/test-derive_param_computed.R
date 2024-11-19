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
    inner_join(
      input %>% filter(PARAMCD == "DIABP") %>% select(USUBJID, VISIT, AVAL),
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
      set_values_to = exprs(
        AVAL = (AVAL.SYSBP + 2 * AVAL.DIABP) / 3,
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
      set_values_to = exprs(
        AVAL = AVAL.WEIGHT / (AVAL.HEIGHT / 100)^2,
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
      set_values_to = exprs(
        AVAL = (AVAL.SYSBP + 2 * AVAL.DIABP) / 3,
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
      set_values_to = exprs(
        AVAL = (AVAL.SYSBP + 2 * AVAL.DIABP) / 3,
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
    ~USUBJID, ~AVISIT,  ~PARAMCD, ~QSSTRESN, ~AVAL,
    "1",      "WEEK 2", "CHSF12", 1,             6,
    "1",      "WEEK 2", "CHSF14", 1,             6,
    "1",      "WEEK 4", "CHSF12", 2,            12,
    "1",      "WEEK 4", "CHSF14", 1,             6
  ) %>% mutate(
    QSORRES = NA_character_
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
      set_values_to = exprs(
        AVAL = case_when(
          QSORRES.CHSF13 == "Not applicable" ~ 0,
          QSORRES.CHSF13 == "Yes" ~ 38,
          QSORRES.CHSF13 == "No" ~ if_else(
            QSSTRESN.CHSF12 > QSSTRESN.CHSF14,
            25,
            0
          )
        ),
        PARAMCD = "CHSF13"
      )
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
      set_values_to = exprs(
        AVAL = case_when(
          QSORRES.CHSF13 == "Not applicable" ~ 0,
          QSORRES.CHSF13 == "Yes" ~ 38,
          QSORRES.CHSF13 == "No" ~ if_else(
            QSSTRESN.CHSF12 > QSSTRESN.CHSF14,
            25,
            0
          )
        ),
        PARAMCD = "CHSF13"
      )
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
      set_values_to = exprs(
        AVAL = AVAL.WEIGHT / (VSSTRESN.HEIGHT / 100)^2,
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
      set_values_to = exprs(
        AVAL = AVAL.WEIGHT / (AVAL.HEIGHT / 100)^2,
        PARAMCD = "BMI",
        PARAM = "Body Mass Index (kg/m2)",
        AVALU = "kg/m2"
      )
    ),
    regexp = "The input dataset does not contain any observations fullfiling the filter"
  )

  expect_dfs_equal(
    output,
    input,
    keys = c("USUBJID", "PARAMCD", "VISIT")
  )
})

## Test 9: compute multiple variables, keep_nas ----
test_that("derive_param_computed Test 9: compute multiple variables, keep_nas", {
  adlb_tbilialk <- tibble::tribble(
    ~USUBJID, ~PARAMCD, ~AVALC, ~ADTM,        ~ADTF,
    "1",      "ALK2",   "Y",    "2021-05-13", NA_character_,
    "1",      "TBILI2", "Y",    "2021-06-30", "D",
    "2",      "ALK2",   "Y",    "2021-12-31", "M",
    "2",      "TBILI2", "N",    "2021-11-11", NA_character_,
    "3",      "ALK2",   "N",    "2021-04-03", NA_character_,
    "3",      "TBILI2", "N",    "2021-04-04", NA_character_
  ) %>%
    mutate(ADTM = lubridate::ymd(ADTM))

  expected <- tibble::tribble(
    ~USUBJID, ~AVALC, ~ADTM,        ~ADTF,
    "1",      "Y",    "2021-06-30", "D",
    "2",      "N",    "2021-12-31", "M",
    "3",      "N",    "2021-04-04", NA_character_
  ) %>%
    mutate(
      ADTM = lubridate::ymd(ADTM),
      PARAMCD = "TB2AK2",
      PARAM = "TBILI > 2 times ULN and ALKPH <= 2 times ULN"
    )

  actual <- derive_param_computed(
    dataset_add = adlb_tbilialk,
    by_vars = exprs(USUBJID),
    parameters = c("ALK2", "TBILI2"),
    set_values_to = exprs(
      AVALC = if_else(AVALC.TBILI2 == "Y" & AVALC.ALK2 == "Y", "Y", "N"),
      ADTM = pmax(ADTM.TBILI2, ADTM.ALK2),
      ADTF = if_else(ADTM == ADTM.TBILI2, ADTF.TBILI2, ADTF.ALK2),
      PARAMCD = "TB2AK2",
      PARAM = "TBILI > 2 times ULN and ALKPH <= 2 times ULN"
    ),
    keep_nas = TRUE
  )

  expect_dfs_equal(
    base = expected,
    compare = actual,
    keys = c("USUBJID")
  )
})

## Test 10: specify variables for keep_nas ----
test_that("derive_param_computed Test 10: specify variables for keep_nas", {
  advs <- tibble::tribble(
    ~USUBJID,      ~PARAMCD, ~PARAM,                            ~AVAL, ~AVISIT,
    "01-701-1015", "DIABP",  "Diastolic Blood Pressure (mmHg)",    51, "BASELINE",
    "01-701-1015", "DIABP",  "Diastolic Blood Pressure (mmHg)",    50, "WEEK 2",
    "01-701-1015", "SYSBP",  "Systolic Blood Pressure (mmHg)",    121, "BASELINE",
    "01-701-1015", "SYSBP",  "Systolic Blood Pressure (mmHg)",    121, "WEEK 2",
    "01-701-1028", "DIABP",  "Diastolic Blood Pressure (mmHg)",    79, "BASELINE",
    "01-701-1028", "DIABP",  "Diastolic Blood Pressure (mmHg)",    80, "WEEK 2",
    "01-701-1028", "SYSBP",  "Systolic Blood Pressure (mmHg)",    130, "BASELINE",
    "01-701-1028", "SYSBP",  "Systolic Blood Pressure (mmHg)",     NA, "WEEK 2"
  ) %>%
    mutate(
      ADT = case_when(
        AVISIT == "BASELINE" ~ as.Date("2024-01-10"),
        AVISIT == "WEEK 2" ~ as.Date("2024-01-24")
      ),
      ADTF = NA_character_
    )

  expected <- bind_rows(
    advs,
    tibble::tribble(
      ~USUBJID,      ~AVISIT,    ~AVAL,
      "01-701-1015", "BASELINE", (121 + 2 * 51) / 3,
      "01-701-1015", "WEEK 2",   (121 + 2 * 50) / 3,
      "01-701-1028", "BASELINE", (130 + 2 * 79) / 3
    ) %>%
      mutate(
        PARAMCD = "MAP",
        PARAM = "Mean Arterial Pressure (mmHg)",
        ADT = case_when(
          AVISIT == "BASELINE" ~ as.Date("2024-01-10"),
          AVISIT == "WEEK 2" ~ as.Date("2024-01-24")
        ),
        ADTF = NA_character_
      )
  )

  actual <- derive_param_computed(
    advs,
    by_vars = exprs(USUBJID, AVISIT),
    parameters = c("SYSBP", "DIABP"),
    set_values_to = exprs(
      AVAL = (AVAL.SYSBP + 2 * AVAL.DIABP) / 3,
      PARAMCD = "MAP",
      PARAM = "Mean Arterial Pressure (mmHg)",
      ADT = ADT.SYSBP,
      ADTF = ADTF.SYSBP
    ),
    keep_nas = exprs(ADTF)
  )

  expect_dfs_equal(
    base = expected,
    compare = actual,
    keys = c("USUBJID", "PARAMCD", "AVISIT")
  )
})

## Test 11: error if keep_nas is invalid ----
test_that("derive_param_computed Test 11: error if keep_nas is invalid", {
  advs <- tibble::tribble(
    ~USUBJID,      ~PARAMCD, ~PARAM,                            ~AVAL, ~AVISIT,
    "01-701-1015", "DIABP",  "Diastolic Blood Pressure (mmHg)",    51, "BASELINE",
    "01-701-1015", "DIABP",  "Diastolic Blood Pressure (mmHg)",    50, "WEEK 2",
    "01-701-1015", "SYSBP",  "Systolic Blood Pressure (mmHg)",    121, "BASELINE",
    "01-701-1015", "SYSBP",  "Systolic Blood Pressure (mmHg)",    121, "WEEK 2",
    "01-701-1028", "DIABP",  "Diastolic Blood Pressure (mmHg)",    79, "BASELINE",
    "01-701-1028", "DIABP",  "Diastolic Blood Pressure (mmHg)",    80, "WEEK 2",
    "01-701-1028", "SYSBP",  "Systolic Blood Pressure (mmHg)",    130, "BASELINE",
    "01-701-1028", "SYSBP",  "Systolic Blood Pressure (mmHg)",     NA, "WEEK 2"
  )

  expect_snapshot(
    derive_param_computed(
      advs,
      by_vars = exprs(USUBJID, AVISIT),
      parameters = c("SYSBP", "DIABP"),
      set_values_to = exprs(
        AVAL = (AVAL.SYSBP + 2 * AVAL.DIABP) / 3,
        PARAMCD = "MAP",
        PARAM = "Mean Arterial Pressure (mmHg)"
      ),
      keep_nas = 3
    ),
    error = TRUE
  )
})

## Test 12: inform if no new records due to NAs ----
test_that("derive_param_computed Test 12: inform if no new records due to NAs", {
  advs <- tibble::tribble(
    ~USUBJID,      ~PARAMCD, ~PARAM,                            ~AVAL, ~AVISIT,
    "01-701-1015", "DIABP",  "Diastolic Blood Pressure (mmHg)",    51, "BASELINE",
    "01-701-1015", "DIABP",  "Diastolic Blood Pressure (mmHg)",    50, "WEEK 2",
    "01-701-1015", "SYSBP",  "Systolic Blood Pressure (mmHg)",    121, "BASELINE",
    "01-701-1015", "SYSBP",  "Systolic Blood Pressure (mmHg)",    121, "WEEK 2",
    "01-701-1028", "DIABP",  "Diastolic Blood Pressure (mmHg)",    79, "BASELINE",
    "01-701-1028", "DIABP",  "Diastolic Blood Pressure (mmHg)",    80, "WEEK 2",
    "01-701-1028", "SYSBP",  "Systolic Blood Pressure (mmHg)",    130, "BASELINE",
    "01-701-1028", "SYSBP",  "Systolic Blood Pressure (mmHg)",     NA, "WEEK 2"
  ) %>%
    mutate(
      ADT = case_when(
        AVISIT == "BASELINE" ~ as.Date("2024-01-10"),
        AVISIT == "WEEK 2" ~ as.Date("2024-01-24")
      ),
      ADTF = NA_character_
    )

  expect_snapshot(derive_param_computed(
    advs,
    by_vars = exprs(USUBJID, AVISIT),
    parameters = c("SYSBP", "DIABP"),
    set_values_to = exprs(
      AVAL = (AVAL.SYSBP + 2 * AVAL.DIABP) / 3,
      PARAMCD = "MAP",
      PARAM = "Mean Arterial Pressure (mmHg)",
      ADT = ADT.SYSBP,
      ADTF = ADTF.SYSBP
    )
  ))
})

# assert_parameters_argument ----
## Test 13: error if argument is of wrong type ----
test_that("assert_parameters_argument Test 13: error if argument is of wrong type", {
  expect_snapshot(
    assert_parameters_argument(myparameters <- c(1, 2, 3)),
    error = TRUE
  )
})

# get_hori_data ----
## Test 14: error if variables with more than one dot ----
test_that("get_hori_data Test 14: error if variables with more than one dot", {
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

  expect_snapshot(
    get_hori_data(
      input,
      parameters = exprs(SYSBP, DIABP),
      by_vars = exprs(USUBJID, VISIT),
      set_values_to = exprs(AVAL = (AVAL.SYSBP + 2 * AVAL.DIA.BP) / 3),
      filter = NULL
    ),
    error = TRUE
  )
})
