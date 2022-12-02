# create_period_dataset ----
## Test 1: periods ----
test_that("create_period_dataset Test 1: periods", {
  adsl <- tibble::tribble(
    ~USUBJID, ~AP01SDT,     ~AP01EDT,     ~AP02SDT,     ~AP02EDT,
    "1",      "2021-01-04", "2021-02-06", "2021-02-07", "2021-03-07",
    "2",      "2021-02-02", "2021-03-02", "2021-03-03", "2021-04-01"
  ) %>%
    mutate(
      dplyr::across(matches("AP\\d\\d[ES]DT"), ymd)
    ) %>%
    mutate(
      STUDYID = "xyz"
    )

  expected <- tibble::tribble(
    ~USUBJID, ~APERIOD, ~APERSDT,     ~APEREDT,
    "1",             1, "2021-01-04", "2021-02-06",
    "1",             2, "2021-02-07", "2021-03-07",
    "2",             1, "2021-02-02", "2021-03-02",
    "2",             2, "2021-03-03", "2021-04-01"
  ) %>%
    mutate(
      STUDYID = "xyz",
      APERIOD = as.integer(APERIOD),
      dplyr::across(matches("APER[ES]DT"), ymd)
    )

  expect_dfs_equal(
    base = expected,
    compare = create_period_dataset(
      adsl,
      new_vars = vars(APERSDT = APxxSDT, APEREDT = APxxEDT)
    ),
    keys = c("USUBJID", "APERIOD")
  )
})

## Test 2: phases ----
test_that("create_period_dataset Test 2: phases", {
  adsl <- tibble::tribble(
    ~USUBJID, ~PH1SDT,      ~PH1EDT,      ~PH2SDT,      ~PH2EDT,      ~APHASE1,    ~APHASE2,
    "1",      "2021-01-04", "2021-02-06", "2021-02-07", "2021-03-07", "TREATMENT", "FUP",
    "2",      "2021-02-02", "2021-03-02", NA,           NA,           "TREATMENT", NA
  ) %>%
    mutate(
      dplyr::across(matches("PH\\d[ES]DT"), ymd)
    ) %>%
    mutate(
      STUDYID = "xyz"
    )

  expected <- tibble::tribble(
    ~USUBJID, ~APHASEN, ~PHSDT,       ~PHEDT,       ~APHASE,
    "1",             1, "2021-01-04", "2021-02-06", "TREATMENT",
    "1",             2, "2021-02-07", "2021-03-07", "FUP",
    "2",             1, "2021-02-02", "2021-03-02", "TREATMENT"
  ) %>%
    mutate(
      STUDYID = "xyz",
      APHASEN = as.integer(APHASEN),
      dplyr::across(matches("PH[ES]DT"), ymd)
    )

  expect_dfs_equal(
    base = expected,
    compare = create_period_dataset(
      adsl,
      new_vars = vars(PHSDT = PHwSDT, PHEDT = PHwEDT, APHASE = APHASEw)
    ),
    keys = c("USUBJID", "APHASEN")
  )
})

## Test 3: subperiods ----
test_that("create_period_dataset Test 3: subperiods", {
  adsl <- tibble::tribble(
    ~USUBJID, ~P01S1SDT,    ~P01S1EDT,    ~P01S2SDT,    ~P01S2EDT,    ~P02S1SDT,    ~P02S1EDT,
    "1",      "2021-01-04", "2021-01-19", "2021-01-20", "2021-02-06", "2021-02-07", "2021-03-07",
    "2",      "2021-02-02", "2021-03-02", NA,           NA,           "2021-03-03", "2021-04-01"
  ) %>%
    mutate(
      dplyr::across(matches("PH\\d\\dS\\d[ES]DT"), ymd)
    ) %>%
    mutate(
      STUDYID = "xyz"
    )

  expected <- tibble::tribble(
    ~USUBJID, ~APERIOD, ~ASPER, ~ASPRSDT,     ~ASPREDT,
    "1",             1,      1, "2021-01-04", "2021-01-19",
    "1",             1,      2, "2021-01-20", "2021-02-06",
    "1",             2,      1, "2021-02-07", "2021-03-07",
    "2",             1,      1, "2021-02-02", "2021-03-02",
    "2",             2,      1, "2021-03-03", "2021-04-01"
  ) %>%
    mutate(
      STUDYID = "xyz",
      APERIOD = as.integer(APERIOD),
      ASPER = as.integer(ASPER),
      dplyr::across(matches("APER[ES]DT"), ymd)
    )

  expect_dfs_equal(
    base = expected,
    compare = create_period_dataset(
      adsl,
      new_vars = vars(ASPRSDT = PxxSwSDT, ASPREDT = PxxSwEDT)
    ),
    keys = c("USUBJID", "APERIOD", "ASPER")
  )
})

## Test 4: error if no period/phase variable on RHS ----
test_that("create_period_dataset Test 4: error if no period/phase variable on RHS", {
  adsl <- tibble::tribble(
    ~USUBJID, ~AP01SDT,     ~AP01EDT,     ~AP02SDT,     ~AP02EDT,
    "1",      "2021-01-04", "2021-02-06", "2021-02-07", "2021-03-07",
    "2",      "2021-02-02", "2021-03-02", "2021-03-03", "2021-04-01"
  ) %>%
    mutate(
      dplyr::across(matches("AP\\d\\d[ES]DT"), ymd)
    ) %>%
    mutate(
      STUDYID = "xyz"
    )

  expect_error(
    create_period_dataset(
      adsl,
      new_vars = vars(USUBJ = USUBJID)
    ),
    regexp = paste(
      paste0(
        "The right hand side values of `new_vars` have to be CDISC style ",
        "subperiod, period, or phase variables."
      ),
      "I.e., they must contain the xx or w fragment, e.g., APxxSDT, PxxSwSDT, or PHwSDT.",
      sep = "\n"
    ),
    fixed = TRUE
  )
})

## Test 5: error if different type of RHSs ----
test_that("create_period_dataset Test 5: error if different type of RHSs", {
  adsl <- tibble::tribble(
    ~USUBJID, ~AP01SDT,     ~AP01EDT,     ~AP02SDT,     ~AP02EDT,
    "1",      "2021-01-04", "2021-02-06", "2021-02-07", "2021-03-07",
    "2",      "2021-02-02", "2021-03-02", "2021-03-03", "2021-04-01"
  ) %>%
    mutate(
      dplyr::across(matches("AP\\d\\d[ES]DT"), ymd)
    ) %>%
    mutate(
      STUDYID = "xyz"
    )

  expect_error(
    create_period_dataset(
      adsl,
      new_vars = vars(APERSDT = APxxSDT, ASPRSDT = PxxSwSDT)
    ),
    regexp = paste(
      "More than one type of subperiod, period, or phase variables is specified for `new_vars`:",
      "subperiod: `PxxSwSDT`",
      "period: `APxxSDT`",
      sep = "\n"
    ),
    fixed = TRUE
  )
})

## Test 6: error if RHS variable not in input dataset ----
test_that("create_period_dataset Test 6: error if RHS variable not in input dataset", {
  adsl <- tibble::tribble(
    ~USUBJID, ~AP01SDT,     ~AP01EDT,     ~AP02SDT,     ~AP02EDT,
    "1",      "2021-01-04", "2021-02-06", "2021-02-07", "2021-03-07",
    "2",      "2021-02-02", "2021-03-02", "2021-03-03", "2021-04-01"
  ) %>%
    mutate(
      dplyr::across(matches("AP\\d\\d[ES]DT"), ymd)
    ) %>%
    mutate(
      STUDYID = "xyz"
    )

  expect_error(
    create_period_dataset(
      adsl,
      new_vars = vars(PHSDT = PHwSDT)
    ),
    regexp = "No variables of the form PHwSDT were found in the input dataset.",
    fixed = TRUE
  )
})

# derive_vars_period ----
## Test 7: periods ----
test_that("derive_vars_period Test 7: periods", {
  expected <- tibble::tribble(
    ~USUBJID, ~AP01SDT,     ~AP01EDT,     ~AP02SDT,     ~AP02EDT,
    "1",      "2021-01-04", "2021-02-06", "2021-02-07", "2021-03-07",
    "2",      "2021-02-02", "2021-03-02", "2021-03-03", "2021-04-01"
  ) %>%
    mutate(
      dplyr::across(matches("AP\\d\\d[ES]DT"), ymd)
    ) %>%
    mutate(
      STUDYID = "xyz"
    )

  period_ref <- tibble::tribble(
    ~USUBJID, ~APERIOD, ~APERSDT,     ~APEREDT,
    "1",             1, "2021-01-04", "2021-02-06",
    "1",             2, "2021-02-07", "2021-03-07",
    "2",             1, "2021-02-02", "2021-03-02",
    "2",             2, "2021-03-03", "2021-04-01"
  ) %>%
    mutate(
      STUDYID = "xyz",
      APERIOD = as.integer(APERIOD),
      dplyr::across(matches("APER[ES]DT"), ymd)
    )

  adsl <- tibble::tibble(STUDYID = "xyz", USUBJID = c("1", "2"))

  expect_dfs_equal(
    base = expected,
    compare = derive_vars_period(
      adsl,
      dataset_ref = period_ref,
      new_vars = vars(APxxSDT = APERSDT, APxxEDT = APEREDT)
    ),
    keys = c("USUBJID")
  )
})

## Test 8: phases ----
test_that("derive_vars_period Test 8: phases", {
  expected <- tibble::tribble(
    ~USUBJID, ~PH1SDT,      ~PH1EDT,      ~PH2SDT,      ~PH2EDT,      ~APHASE1,    ~APHASE2,
    "1",      "2021-01-04", "2021-02-06", "2021-02-07", "2021-03-07", "TREATMENT", "FUP",
    "2",      "2021-02-02", "2021-03-02", NA,           NA,           "TREATMENT", NA
  ) %>%
    mutate(
      dplyr::across(matches("PH\\d[ES]DT"), ymd)
    ) %>%
    mutate(
      STUDYID = "xyz"
    )

  phase_ref <- tibble::tribble(
    ~USUBJID, ~APHASEN, ~PHSDT,       ~PHEDT,       ~APHASE,
    "1",             1, "2021-01-04", "2021-02-06", "TREATMENT",
    "1",             2, "2021-02-07", "2021-03-07", "FUP",
    "2",             1, "2021-02-02", "2021-03-02", "TREATMENT"
  ) %>%
    mutate(
      STUDYID = "xyz",
      APHASEN = as.integer(APHASEN),
      dplyr::across(matches("PH[ES]DT"), ymd)
    )

  adsl <- tibble(STUDYID = "xyz", USUBJID = c("1", "2"))

  expect_dfs_equal(
    base = expected,
    compare = derive_vars_period(
      adsl,
      dataset_ref = phase_ref,
      new_vars = vars(PHwSDT = PHSDT, PHwEDT = PHEDT, APHASEw = APHASE)
    ),
    keys = c("USUBJID")
  )
})

## Test 9: subperiods ----
test_that("derive_vars_period Test 9: subperiods", {
  expected <- tibble::tribble(
    ~USUBJID, ~P01S1SDT,    ~P01S1EDT,    ~P01S2SDT,    ~P01S2EDT,    ~P02S1SDT,    ~P02S1EDT,
    "1",      "2021-01-04", "2021-01-19", "2021-01-20", "2021-02-06", "2021-02-07", "2021-03-07",
    "2",      "2021-02-02", "2021-03-02", NA,           NA,           "2021-03-03", "2021-04-01"
  ) %>%
    mutate(
      dplyr::across(matches("PH\\d\\dS\\d[ES]DT"), ymd)
    ) %>%
    mutate(
      STUDYID = "xyz"
    )

  subperiod_ref <- tibble::tribble(
    ~USUBJID, ~APERIOD, ~ASPER, ~ASPRSDT,     ~ASPREDT,
    "1",             1,      1, "2021-01-04", "2021-01-19",
    "1",             1,      2, "2021-01-20", "2021-02-06",
    "1",             2,      1, "2021-02-07", "2021-03-07",
    "2",             1,      1, "2021-02-02", "2021-03-02",
    "2",             2,      1, "2021-03-03", "2021-04-01"
  ) %>%
    mutate(
      STUDYID = "xyz",
      APERIOD = as.integer(APERIOD),
      ASPER = as.integer(ASPER),
      dplyr::across(matches("APER[ES]DT"), ymd)
    )

  adsl <- tibble(STUDYID = "xyz", USUBJID = c("1", "2"))

  expect_dfs_equal(
    base = expected,
    compare = derive_vars_period(
      adsl,
      dataset_ref = subperiod_ref,
      new_vars = vars(PxxSwSDT = ASPRSDT, PxxSwEDT = ASPREDT)
    ),
    keys = c("USUBJID")
  )
})
## Test 10: error if no period/phase variable on LHS ----
test_that("derive_vars_period Test 10: error if no period/phase variable on LHS", {
  period_ref <- tibble::tribble(
    ~USUBJID, ~APERIOD, ~APERSDT,     ~APEREDT,
    "1",             1, "2021-01-04", "2021-02-06",
    "1",             2, "2021-02-07", "2021-03-07",
    "2",             1, "2021-02-02", "2021-03-02",
    "2",             2, "2021-03-03", "2021-04-01"
  ) %>%
    mutate(
      STUDYID = "xyz",
      APERIOD = as.integer(APERIOD),
      dplyr::across(matches("APER[ES]DT"), ymd)
    )

  adsl <- tibble(STUDYID = "xyz", USUBJID = c("1", "2"))

  expect_error(
    derive_vars_period(
      adsl,
      dataset_ref = period_ref,
      new_vars = vars(USUBJ = USUBJID)
    ),
    regexp = paste(
      paste0(
        "The left hand side values of `new_vars` have to be CDISC style ",
        "subperiod, period, or phase variables."
      ),
      "I.e., they must contain the xx or w fragment, e.g., APxxSDT, PxxSwSDT, or PHwSDT.",
      sep = "\n"
    ),
    fixed = TRUE
  )
})

## Test 11: error if different type of LHSs ----
test_that("derive_vars_period Test 11: error if different type of LHSs", {
  period_ref <- tibble::tribble(
    ~USUBJID, ~APERIOD, ~APERSDT,     ~APEREDT,
    "1",             1, "2021-01-04", "2021-02-06",
    "1",             2, "2021-02-07", "2021-03-07",
    "2",             1, "2021-02-02", "2021-03-02",
    "2",             2, "2021-03-03", "2021-04-01"
  ) %>%
    mutate(
      STUDYID = "xyz",
      APERIOD = as.integer(APERIOD),
      dplyr::across(matches("APER[ES]DT"), ymd)
    )

  adsl <- tibble(STUDYID = "xyz", USUBJID = c("1", "2"))

  expect_error(
    derive_vars_period(
      adsl,
      dataset_ref = period_ref,
      new_vars = vars(APxxSDT = APERSDT, PxxSwSDT = ASPRSDT)
    ),
    regexp = paste(
      "More than one type of subperiod, period, or phase variables is specified for `new_vars`:",
      "subperiod: `PxxSwSDT`",
      "period: `APxxSDT`",
      sep = "\n"
    ),
    fixed = TRUE
  )
})
