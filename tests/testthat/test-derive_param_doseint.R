context("test-derive_param_doseint")


test_that("new observations are derived correctly when zero_doses is NULL", {
  input <- tibble::tribble(
    ~USUBJID,      ~PARAMCD, ~PARAM,         ~AVAL, ~VISIT,
    "01-701-1015", "TNDOSE", "Total Admin",  52,    "WEEK 1",
    "01-701-1015", "TSNDOSE","Planned Dose", 96,    "WEEK 1",
    "01-701-1015", "TNDOSE", "Total Admin",  52,    "WEEK 2",
    "01-701-1015", "TSNDOSE","Planned Dose", 0,     "WEEK 2",
    "01-701-1028", "TNDOSE", "Total Admin",  52,    "WEEK 1",
    "01-701-1028", "TSNDOSE","Planned Dose", NA,    "WEEK 1",
    "01-701-1028", "TNDOSE", "Total Admin",  0,     "WEEK 2",
    "01-701-1028", "TSNDOSE","Planned Dose", 99,    "WEEK 2",
    "01-701-1028", "TNDOSE", "Total Admin",  0,     "WEEK 3",
    "01-701-1028", "TSNDOSE","Planned Dose", NA,    "WEEK 3",
    "01-701-1028", "TSNDOSE","Total Admin",  0,     "WEEK 4",
    "01-701-1028", "TNDOSE", "Total Admin",  NA,    "WEEK 5",
    "01-701-1028", "TSNDOSE","Planned Dose", 0,     "WEEK 6",
    "01-701-1028", "TSNDOSE","Planned Dose", NA,    "WEEK 7",
  )

  new_obs <-
    inner_join(input %>% filter(PARAMCD == "TNDOSE" & !is.na(AVAL)) %>% select(USUBJID, VISIT, AVAL),
               input %>% filter(PARAMCD == "TSNDOSE" & !is.na(AVAL)) %>% select(USUBJID, VISIT, AVAL),
               by = c("USUBJID", "VISIT"),
               suffix = c(".TNDOSE", ".TSNDOSE")) %>%
    mutate(AVAL = AVAL.TNDOSE / AVAL.TSNDOSE * 100,
           PARAMCD = "TNDOSINT",
           PARAM = "Dose Intensity (%)") %>%
    select(-AVAL.TSNDOSE, -AVAL.TNDOSE)
  expected_output <- bind_rows(input, new_obs)

  expect_dfs_equal(derive_param_doseint(input,
                                     by_vars = vars(USUBJID, VISIT)),
                   expected_output,
                   keys = c("USUBJID", "PARAMCD", "VISIT"))
})

test_that("new observations are derived correctly when zero_doses is Y", {
  input <- tibble::tribble(
    ~USUBJID,      ~PARAMCD, ~PARAM,         ~AVAL, ~VISIT,
    "01-701-1015", "TNDOSE", "Total Admin",  52,    "WEEK 1",
    "01-701-1015", "TSNDOSE","Planned Dose", 96,    "WEEK 1",
    "01-701-1015", "TNDOSE", "Total Admin",  52,    "WEEK 2",
    "01-701-1015", "TSNDOSE","Planned Dose", 0,     "WEEK 2",
    "01-701-1015", "TNDOSE", "Total Admin",  0,     "WEEK 3",
    "01-701-1015", "TSNDOSE","Planned Dose", 0,     "WEEK 3",
    "01-701-1028", "TNDOSE", "Total Admin",  52,    "WEEK 1",
    "01-701-1028", "TSNDOSE","Planned Dose", NA,    "WEEK 1",
    "01-701-1028", "TNDOSE", "Total Admin",  0,     "WEEK 2",
    "01-701-1028", "TSNDOSE","Planned Dose", 99,    "WEEK 2",
    "01-701-1028", "TNDOSE", "Total Admin",  0,     "WEEK 3",
    "01-701-1028", "TSNDOSE","Planned Dose", NA,    "WEEK 3",
    "01-701-1028", "TSNDOSE","Total Admin",  0,     "WEEK 4",
    "01-701-1028", "TNDOSE", "Total Admin",  NA,    "WEEK 5",
    "01-701-1028", "TSNDOSE","Planned Dose", 0,     "WEEK 6",
    "01-701-1028", "TSNDOSE","Planned Dose", NA,    "WEEK 7",
  )

  new_obs <-
    inner_join(input %>% filter(PARAMCD == "TNDOSE" & !is.na(AVAL)) %>% select(USUBJID, VISIT, AVAL),
               input %>% filter(PARAMCD == "TSNDOSE" & !is.na(AVAL)) %>% select(USUBJID, VISIT, AVAL),
               by = c("USUBJID", "VISIT"),
               suffix = c(".TNDOSE", ".TSNDOSE")) %>%
    mutate(AVAL = case_when(AVAL.TSNDOSE == 0 & AVAL.TNDOSE > 0 ~ 100,
                            AVAL.TSNDOSE == 0 & AVAL.TNDOSE == 0 ~ 0,
                            TRUE ~ AVAL.TNDOSE / AVAL.TSNDOSE * 100),
           PARAMCD = "TNDOSINT",
           PARAM = "Dose Intensity (%)") %>%
    select(-AVAL.TSNDOSE, -AVAL.TNDOSE)
  expected_output <- bind_rows(input, new_obs)

  expect_dfs_equal(derive_param_doseint(input,
                                        by_vars = vars(USUBJID, VISIT),
                                        zero_doses = "Y"),
                   expected_output,
                   keys = c("USUBJID", "PARAMCD", "VISIT"))
})

