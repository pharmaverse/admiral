test_that("new observations are derived correctly when zero_doses is NULL", {
  # nolint start
  input <- tibble::tribble(
    ~USUBJID, ~PARAMCD, ~AVAL, ~VISIT,
    "01-701-1015", "TNDOSE", 52, "WEEK 1",
    "01-701-1015", "TSNDOSE", 96, "WEEK 1",
    "01-701-1015", "TNDOSE", 52, "WEEK 2",
    "01-701-1015", "TSNDOSE", 0, "WEEK 2",
    "01-701-1028", "TNDOSE", 52, "WEEK 1",
    "01-701-1028", "TSNDOSE", NA, "WEEK 1",
    "01-701-1028", "TNDOSE", 0, "WEEK 2",
    "01-701-1028", "TSNDOSE", 99, "WEEK 2",
    "01-701-1028", "TNDOSE", 0, "WEEK 3",
    "01-701-1028", "TSNDOSE", NA, "WEEK 3",
    "01-701-1028", "TSNDOSE", 0, "WEEK 4",
    "01-701-1028", "TNDOSE", NA, "WEEK 5",
    "01-701-1028", "TSNDOSE", 0, "WEEK 6",
    "01-701-1028", "TSNDOSE", NA, "WEEK 7",
    "01-701-1028", "TNDOSE", 0, "WEEK 8",
    "01-701-1028", "TSNDOSE", 0, "WEEK 8",
  )
  # nolint end

  new_obs <-
    inner_join(
      input %>% filter(PARAMCD == "TNDOSE" & !is.na(AVAL)) %>% select(USUBJID, VISIT, AVAL),
      input %>% filter(PARAMCD == "TSNDOSE" & !is.na(AVAL)) %>% select(USUBJID, VISIT, AVAL),
      by = c("USUBJID", "VISIT"),
      suffix = c(".TNDOSE", ".TSNDOSE")
    ) %>%
    mutate(
      AVAL = AVAL.TNDOSE / AVAL.TSNDOSE * 100,
      PARAMCD = "TNDOSINT"
    ) %>%
    select(-AVAL.TSNDOSE, -AVAL.TNDOSE)
  expected_output <- bind_rows(input, new_obs)

  expect_dfs_equal(
    derive_param_doseint(input,
      by_vars = vars(USUBJID, VISIT)
    ),
    expected_output,
    keys = c("USUBJID", "PARAMCD", "VISIT")
  )
})

test_that("new observations are derived correctly when zero_doses is Y", {
  # nolint start
  input <- tibble::tribble(
    ~USUBJID, ~PARAMCD, ~AVAL, ~VISIT,
    "01-701-1015", "TNDOSE", 52, "WEEK 1",
    "01-701-1015", "TSNDOSE", 96, "WEEK 1",
    "01-701-1015", "TNDOSE", 52, "WEEK 2",
    "01-701-1015", "TSNDOSE", 0, "WEEK 2",
    "01-701-1015", "TNDOSE", 0, "WEEK 3",
    "01-701-1015", "TSNDOSE", 0, "WEEK 3",
    "01-701-1028", "TNDOSE", 52, "WEEK 1",
    "01-701-1028", "TSNDOSE", NA, "WEEK 1",
    "01-701-1028", "TNDOSE", 0, "WEEK 2",
    "01-701-1028", "TSNDOSE", 99, "WEEK 2",
    "01-701-1028", "TNDOSE", 0, "WEEK 3",
    "01-701-1028", "TSNDOSE", NA, "WEEK 3",
    "01-701-1028", "TSNDOSE", 0, "WEEK 4",
    "01-701-1028", "TNDOSE", NA, "WEEK 5",
    "01-701-1028", "TSNDOSE", 0, "WEEK 6",
    "01-701-1028", "TSNDOSE", NA, "WEEK 7",
  )
  # nolint end

  new_obs <-
    inner_join(
      input %>% filter(PARAMCD == "TNDOSE" & !is.na(AVAL)) %>% select(USUBJID, VISIT, AVAL),
      input %>% filter(PARAMCD == "TSNDOSE" & !is.na(AVAL)) %>% select(USUBJID, VISIT, AVAL),
      by = c("USUBJID", "VISIT"),
      suffix = c(".TNDOSE", ".TSNDOSE")
    ) %>%
    mutate(
      AVAL = case_when(
        AVAL.TSNDOSE == 0 & AVAL.TNDOSE > 0 ~ 100,
        AVAL.TSNDOSE == 0 & AVAL.TNDOSE == 0 ~ 0,
        TRUE ~ AVAL.TNDOSE / AVAL.TSNDOSE * 100
      ),
      PARAMCD = "TNDOSINT"
    ) %>%
    select(-AVAL.TSNDOSE, -AVAL.TNDOSE)
  expected_output <- bind_rows(input, new_obs)

  expect_dfs_equal(
    derive_param_doseint(input,
      by_vars = vars(USUBJID, VISIT),
      zero_doses = "100"
    ),
    expected_output,
    keys = c("USUBJID", "PARAMCD", "VISIT")
  )
})
