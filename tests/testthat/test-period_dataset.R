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
      new_vars = exprs(APERSDT = APxxSDT, APEREDT = APxxEDT)
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
      new_vars = exprs(PHSDT = PHwSDT, PHEDT = PHwEDT, APHASE = APHASEw)
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
      new_vars = exprs(ASPRSDT = PxxSwSDT, ASPREDT = PxxSwEDT)
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

  expect_snapshot(
    create_period_dataset(
      adsl,
      new_vars = exprs(USUBJ = USUBJID)
    ),
    error = TRUE
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

  expect_snapshot(
    create_period_dataset(
      adsl,
      new_vars = exprs(APERSDT = APxxSDT, ASPRSDT = PxxSwSDT)
    ),
    error = TRUE
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

  expect_snapshot(
    create_period_dataset(
      adsl,
      new_vars = exprs(PHSDT = PHwSDT)
    ),
    error = TRUE
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
      new_vars = exprs(APxxSDT = APERSDT, APxxEDT = APEREDT)
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
      new_vars = exprs(PHwSDT = PHSDT, PHwEDT = PHEDT, APHASEw = APHASE)
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
      new_vars = exprs(PxxSwSDT = ASPRSDT, PxxSwEDT = ASPREDT)
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

  expect_snapshot(
    derive_vars_period(
      adsl,
      dataset_ref = period_ref,
      new_vars = exprs(USUBJ = USUBJID)
    ),
    error = TRUE
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

  expect_snapshot(
    derive_vars_period(
      adsl,
      dataset_ref = period_ref,
      new_vars = exprs(APxxSDT = APERSDT, PxxSwSDT = ASPRSDT)
    ),
    error = TRUE
  )
})

## Test 12: DT and DTM columns exist, pulls only one unique col ----
test_that("create_period_dataset Test 12: DT and DTM columns exist, pulls only one unique col", {
  adsl <- tibble::tribble(
    ~USUBJID, ~AP01SDT,     ~AP01SDTM,                 ~AP01EDT,     ~AP02SDT,     ~AP02EDT,
    "1",      "2021-01-04", "2021-01-04T12:00:00", "2021-02-06", "2021-02-07", "2021-03-07",
    "2",      "2021-02-02", "2021-02-02T12:00:00", "2021-03-02", "2021-03-03", "2021-04-01"
  ) %>%
    mutate(
      dplyr::across(matches("AP\\d\\d[ES]DT\\b"), ymd),
      dplyr::across(matches("AP\\d\\d[ES]DTM"), ymd_hms),
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
      new_vars = exprs(APERSDT = APxxSDT, APEREDT = APxxEDT)
    ),
    keys = c("USUBJID", "APERIOD")
  )
})

## Test 13: Periods derived even when some variables in dataset_ref are unused ----
test_that("derive_vars_period Test 13: periods", {
  expected <- tibble::tribble(
    ~USUBJID, ~AP01SDT,     ~AP01EDT,
    "1",      "2021-01-04", "2021-02-06",
    "2",      "2021-02-02", "2021-03-02"
  ) %>%
    mutate(
      across(matches("AP\\d\\d[ES]DT"), ymd)
    ) %>%
    mutate(
      STUDYID = "xyz"
    )

  period_ref <- tibble::tribble(
    ~USUBJID, ~APERIOD, ~APERSDT,     ~APEREDT,     ~APERSDTM,             ~APEREDTM,
    "1",             1, "2021-01-04", "2021-02-06", "2021-01-04 10:01:02", "2021-02-06 14:03:04",
    "2",             1, "2021-02-02", "2021-03-02", "2021-02-02 04:11:34", "2021-03-02 16:55:59"
  ) %>%
    mutate(
      STUDYID = "xyz",
      APERIOD = as.integer(APERIOD),
      APERSDT = ymd(APERSDT),
      APEREDT = ymd(APEREDT),
      APERSDTM = ymd_hms(APERSDTM),
      APEREDTM = ymd_hms(APEREDTM)
    )

  adsl <- tibble::tibble(STUDYID = "xyz", USUBJID = c("1", "2"))

  expect_dfs_equal(
    base = expected,
    compare = derive_vars_period(
      adsl,
      dataset_ref = period_ref,
      new_vars = exprs(APxxSDT = APERSDT, APxxEDT = APEREDT)
    ),
    keys = c("USUBJID")
  )
})

## Test 14: periods variable derived when only one variable needed ----
test_that("derive_vars_period Test 14: periods variable derived when only one variable needed", {
  expected <- tibble::tribble(
    ~USUBJID, ~AP01SDT,
    "1",      "2021-01-04",
    "2",      "2021-02-02",
  ) %>%
    mutate(
      across(matches("AP\\d\\d[ES]DT"), ymd)
    ) %>%
    mutate(
      STUDYID = "xyz"
    )

  period_ref <- tibble::tribble(
    ~USUBJID, ~APERIOD, ~APERSDT,     ~APEREDT,     ~APERSDTM,             ~APEREDTM,
    "1",             1, "2021-01-04", "2021-02-06", "2021-01-04 10:01:02", "2021-02-06 14:03:04",
    "2",             1, "2021-02-02", "2021-03-02", "2021-02-02 04:11:34", "2021-03-02 16:55:59"
  ) %>%
    mutate(
      STUDYID = "xyz",
      APERIOD = as.integer(APERIOD),
      APERSDT = ymd(APERSDT),
      APEREDT = ymd(APEREDT),
      APERSDTM = ymd_hms(APERSDTM),
      APEREDTM = ymd_hms(APEREDTM)
    )

  adsl <- tibble::tibble(STUDYID = "xyz", USUBJID = c("1", "2"))

  expect_dfs_equal(
    base = expected,
    compare = derive_vars_period(
      adsl,
      dataset_ref = period_ref,
      new_vars = exprs(APxxSDT = APERSDT)
    ),
    keys = c("USUBJID")
  )
})
