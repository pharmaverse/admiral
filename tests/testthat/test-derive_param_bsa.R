context("test-derive_advs_params")

input <- tibble::tribble(
  ~USUBJID,      ~PARAMCD,     ~PARAM,        ~AVAL, ~AVALU,  ~VISIT,
  "01-701-1015", "HEIGHT",     "Height (cm)", 170,   "cm",   "BASELINE",
  "01-701-1015", "WEIGHT",     "Weight (kg)",  75,   "kg",   "BASELINE",
  "01-701-1015", "WEIGHT",     "Weight (kg)",  78,   "kg",   "MONTH 1",
  "01-701-1015", "WEIGHT",     "Weight (kg)",  80,   "kg",   "MONTH 2",
  "01-701-1028", "HEIGHT",     "Height (cm)", 185,   "cm",   "BASELINE",
  "01-701-1028", "WEIGHT",     "Weight (kg)",  90,   "kg",   "BASELINE",
  "01-701-1028", "WEIGHT",     "Weight (kg)",  88,   "kg",   "MONTH 1",
  "01-701-1028", "WEIGHT",     "Weight (kg)",  85,   "kg",   "MONTH 2",
)

test_that("new observations are derived correctly with Mosteller method", {

  new_obs <-
    inner_join(input %>% filter(PARAMCD == "HEIGHT") %>% select(USUBJID, VISIT, AVAL, AVALU),
               input %>% filter(PARAMCD == "WEIGHT") %>% select(USUBJID, VISIT, AVAL),
               by = c("USUBJID", "VISIT"),
               suffix = c(".HEIGHT", ".WEIGHT")) %>%
    mutate(AVAL = sqrt(AVAL.HEIGHT * AVAL.WEIGHT / 3600),
           PARAMCD = "BSA",
           PARAM = "Body Surface Area",
           AVALU = "m^2") %>%
    select(-AVAL.HEIGHT, -AVAL.WEIGHT)
  expected_output <- bind_rows(input, new_obs)

  expect_dfs_equal(
    derive_param_bsa(input, by_vars = vars(USUBJID, VISIT), method = "Mosteller"),
    expected_output,
    keys = c("USUBJID", "PARAMCD", "VISIT")
  )
})

test_that("new observations are derived correctly with DuBois & DuBois method", {

new_obs <-
  inner_join(input %>% filter(PARAMCD == "HEIGHT") %>% select(USUBJID, VISIT, AVAL, AVALU),
             input %>% filter(PARAMCD == "WEIGHT") %>% select(USUBJID, VISIT, AVAL),
             by = c("USUBJID", "VISIT"),
             suffix = c(".HEIGHT", ".WEIGHT")) %>%
  mutate(AVAL = 0.20247 * (AVAL.HEIGHT/100) ^ 0.725 * AVAL.WEIGHT ^ 0.425,
         PARAMCD = "BSA",
         PARAM = "Body Surface Area",
         AVALU = "m^2") %>%
  select(-AVAL.HEIGHT, -AVAL.WEIGHT)
expected_output <- bind_rows(input, new_obs)

expect_dfs_equal(
  derive_param_bsa(input, by_vars = vars(USUBJID, VISIT), method = "DuBois-DuBois"),
  expected_output,
  keys = c("USUBJID", "PARAMCD", "VISIT")
)
})

test_that("new observations are derived correctly with Haycock method", {

  new_obs <-
    inner_join(input %>% filter(PARAMCD == "HEIGHT") %>% select(USUBJID, VISIT, AVAL, AVALU),
               input %>% filter(PARAMCD == "WEIGHT") %>% select(USUBJID, VISIT, AVAL),
               by = c("USUBJID", "VISIT"),
               suffix = c(".HEIGHT", ".WEIGHT")) %>%
    mutate(AVAL = 0.024265 * AVAL.HEIGHT ^ 0.3964 * AVAL.WEIGHT ^ 0.5378,
           PARAMCD = "BSA",
           PARAM = "Body Surface Area",
           AVALU = "m^2") %>%
    select(-AVAL.HEIGHT, -AVAL.WEIGHT)
  expected_output <- bind_rows(input, new_obs)

  expect_dfs_equal(
    derive_param_bsa(input, by_vars = vars(USUBJID, VISIT), method = "Haycock"),
    expected_output,
    keys = c("USUBJID", "PARAMCD", "VISIT")
  )
})

test_that("new observations are derived correctly with Gehan & George method", {

  new_obs <-
    inner_join(input %>% filter(PARAMCD == "HEIGHT") %>% select(USUBJID, VISIT, AVAL, AVALU),
               input %>% filter(PARAMCD == "WEIGHT") %>% select(USUBJID, VISIT, AVAL),
               by = c("USUBJID", "VISIT"),
               suffix = c(".HEIGHT", ".WEIGHT")) %>%
    mutate(AVAL = 0.0235 * AVAL.HEIGHT ^ 0.42246 * AVAL.WEIGHT ^ 0.51456,
           PARAMCD = "BSA",
           PARAM = "Body Surface Area",
           AVALU = "m^2") %>%
    select(-AVAL.HEIGHT, -AVAL.WEIGHT)
  expected_output <- bind_rows(input, new_obs)

  expect_dfs_equal(
    derive_param_bsa(input, by_vars = vars(USUBJID, VISIT), method = "Gehan-George"),
    expected_output,
    keys = c("USUBJID", "PARAMCD", "VISIT")
  )
})

test_that("new observations are derived correctly with Boyd method", {

  new_obs <-
    inner_join(input %>% filter(PARAMCD == "HEIGHT") %>% select(USUBJID, VISIT, AVAL, AVALU),
               input %>% filter(PARAMCD == "WEIGHT") %>% select(USUBJID, VISIT, AVAL),
               by = c("USUBJID", "VISIT"),
               suffix = c(".HEIGHT", ".WEIGHT")) %>%
    mutate(AVAL = 0.0003207 * (AVAL.HEIGHT ^ 0.3) *
             (1000 * AVAL.WEIGHT) ^ (0.7285 - (0.0188 * log10(1000 * AVAL.WEIGHT))),
           PARAMCD = "BSA",
           PARAM = "Body Surface Area",
           AVALU = "m^2") %>%
    select(-AVAL.HEIGHT, -AVAL.WEIGHT)
  expected_output <- bind_rows(input, new_obs)

  expect_dfs_equal(
    derive_param_bsa(input, by_vars = vars(USUBJID, VISIT), method = "Boyd"),
    expected_output,
    keys = c("USUBJID", "PARAMCD", "VISIT")
  )
})

test_that("new observations are derived correctly with Fujimoto method", {

  new_obs <-
    inner_join(input %>% filter(PARAMCD == "HEIGHT") %>% select(USUBJID, VISIT, AVAL, AVALU),
               input %>% filter(PARAMCD == "WEIGHT") %>% select(USUBJID, VISIT, AVAL),
               by = c("USUBJID", "VISIT"),
               suffix = c(".HEIGHT", ".WEIGHT")) %>%
    mutate(AVAL =  0.008883 * AVAL.HEIGHT ^ 0.663 * AVAL.WEIGHT ^ 0.444,
           PARAMCD = "BSA",
           PARAM = "Body Surface Area",
           AVALU = "m^2") %>%
    select(-AVAL.HEIGHT, -AVAL.WEIGHT)
  expected_output <- bind_rows(input, new_obs)

  expect_dfs_equal(
    derive_param_bsa(input, by_vars = vars(USUBJID, VISIT), method = "Fujimoto"),
    expected_output,
    keys = c("USUBJID", "PARAMCD", "VISIT")
  )
})

test_that("new observations are derived correctly with Takahira method", {

  new_obs <-
    inner_join(input %>% filter(PARAMCD == "HEIGHT") %>% select(USUBJID, VISIT, AVAL, AVALU),
               input %>% filter(PARAMCD == "WEIGHT") %>% select(USUBJID, VISIT, AVAL),
               by = c("USUBJID", "VISIT"),
               suffix = c(".HEIGHT", ".WEIGHT")) %>%
    mutate(AVAL =  0.007241 * AVAL.HEIGHT ^ 0.725 * AVAL.WEIGHT ^ 0.425,
           PARAMCD = "BSA",
           PARAM = "Body Surface Area",
           AVALU = "m^2") %>%
    select(-AVAL.HEIGHT, -AVAL.WEIGHT)
  expected_output <- bind_rows(input, new_obs)

  expect_dfs_equal(
    derive_param_bsa(input, by_vars = vars(USUBJID, VISIT), method = "Takahira"),
    expected_output,
    keys = c("USUBJID", "PARAMCD", "VISIT")
  )
})

test_that("an error is issued if an invalid method is specified", {
  expect_error(
    derive_param_bsa(input, by_vars = vars(USUBJID, VISIT), method = "unknown-method"),
    paste("`method` must be one of 'Mosteller', 'DuBois-DuBois', 'Haycock', 'Gehan-George',",
          "'Boyd', 'Fujimoto' or 'Takahira' but is 'unknown-method'")
  )
})

input <- tibble::tribble(
  ~USUBJID,      ~PARAMCD,     ~PARAM,        ~AVAL, ~AVALU,  ~VISIT,
  "01-701-1015", "HEIGHT",     "Height (cm)", 170,   "cm",   "BASELINE",
  "01-701-1015", "WEIGHT",     "Weight (kg)",  75,   "kg",   "BASELINE",
  "01-701-1015", "WEIGHT",     "Weight (kg)",  78,   "kg",   "MONTH 1",
  "01-701-1015", "HEIGHT",     "Height (cm)", 170,   "cm",   "MONTH 2",
  "01-701-1015", "WEIGHT",     "Weight (kg)",  80,   "kg",   "MONTH 2",
  "01-701-1028", "HEIGHT",     "Height (cm)", 185,   "cm",   "BASELINE",
  "01-701-1028", "WEIGHT",     "Weight (kg)",  90,   "kg",   "BASELINE",
  "01-701-1028", "WEIGHT",     "Weight (kg)",  88,   "kg",   "MONTH 1",
  "01-701-1028", "HEIGHT",     "Height (cm)", 185,   "cm",   "MONTH 2",
  "01-701-1028", "WEIGHT",     "Weight (kg)",  85,   "kg",   "MONTH 2",
)

test_that("new observations are derived correctly whenever HEIGHT and WEIGHT are available regardless of visit", {

  new_obs <-
    inner_join(input %>% filter(PARAMCD == "HEIGHT") %>% select(USUBJID, VISIT, AVAL, AVALU),
               input %>% filter(PARAMCD == "WEIGHT") %>% select(USUBJID, VISIT, AVAL),
               by = c("USUBJID", "VISIT"),
               suffix = c(".HEIGHT", ".WEIGHT")) %>%
    mutate(AVAL = sqrt(AVAL.HEIGHT * AVAL.WEIGHT / 3600),
           PARAMCD = "BSA",
           PARAM = "Body Surface Area",
           AVALU = "m^2") %>%
    select(-AVAL.HEIGHT, -AVAL.WEIGHT)
  expected_output <- bind_rows(input, new_obs)

  expect_dfs_equal(
    derive_param_bsa(input, by_vars = vars(USUBJID, VISIT), method = "Mosteller"),
    expected_output,
    keys = c("USUBJID", "PARAMCD", "VISIT")
  )
})
