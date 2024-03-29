# slice_derivation ----
## Test 1: slice derivation ----
test_that("slice_derivation Test 1: slice derivation", {
  advs <- tibble::tribble(
    ~USUBJID, ~VSDTC,       ~VSTPT,             ~VSSEQ,
    "1",      "2020-04-16", NA_character_,      1,
    "1",      "2020-04-16", "BEFORE TREATMENT", 2
  )

  actual <- slice_derivation(
    advs,
    derivation = derive_vars_dtm,
    args = params(
      dtc = VSDTC,
      new_vars_prefix = "A"
    ),
    derivation_slice(
      filter = str_detect(VSTPT, "PRE|BEFORE"),
      args = params(time_imputation = "first")
    ),
    derivation_slice(
      filter = TRUE,
      args = params(time_imputation = "last")
    )
  )

  expected <- mutate(
    advs,
    ADTM = c(ymd_hms("2020-04-16 23:59:59"), ymd_hms("2020-04-16 00:00:00")),
    ATMF = "H"
  )

  expect_dfs_equal(
    base = expected,
    compare = actual,
    keys = c("USUBJID", "VSSEQ")
  )
})

## Test 2: non matching observations ----
test_that("slice_derivation Test 2: non matching observations", {
  advs <- tibble::tribble(
    ~USUBJID, ~VSDTC,       ~VSTPT,             ~VSSEQ,
    "1",      "2020-04-16", NA_character_,      1,
    "1",      "2020-04-16", "BEFORE TREATMENT", 2
  )

  actual <- slice_derivation(
    advs,
    derivation = derive_vars_dtm,
    args = params(
      dtc = VSDTC,
      new_vars_prefix = "A"
    ),
    derivation_slice(
      filter = str_detect(VSTPT, "PRE|BEFORE"),
      args = params(time_imputation = "first")
    )
  )

  expected <- mutate(
    advs,
    ADTM = c(ymd_hms(NA), ymd_hms("2020-04-16 00:00:00")),
    ATMF = c(NA_character_, "H")
  )

  expect_dfs_equal(
    base = expected,
    compare = actual,
    keys = c("USUBJID", "VSSEQ")
  )
})

## Test 3: empty slice ----
test_that("slice_derivation Test 3: empty slice", {
  advs <- tibble::tribble(
    ~USUBJID, ~VSDTC,       ~VSTPT,        ~VSSEQ,
    "1",      "2020-04-16", NA_character_, 1,
    "1",      "2020-04-16", NA_character_, 2
  )

  actual <- slice_derivation(
    advs,
    derivation = derive_vars_dtm,
    args = params(
      dtc = VSDTC,
      new_vars_prefix = "A"
    ),
    derivation_slice(
      filter = str_detect(VSTPT, "PRE|BEFORE"),
      args = params(time_imputation = "first")
    ),
    derivation_slice(
      filter = TRUE,
      args = params(time_imputation = "last")
    )
  )

  expected <- mutate(
    advs,
    ADTM = c(ymd_hms("2020-04-16 23:59:59"), ymd_hms("2020-04-16 23:59:59")),
    ATMF = "H"
  )

  expect_dfs_equal(
    base = expected,
    compare = actual,
    keys = c("USUBJID", "VSSEQ")
  )
})

## Test 4: slice without arguments ----
test_that("slice_derivation Test 4: slice without arguments", {
  advs <- tibble::tribble(
    ~USUBJID, ~VSDTC,       ~VSTPT,             ~VSSEQ,
    "1",      "2020-04-16", NA_character_,      1,
    "1",      "2020-04-16", "BEFORE TREATMENT", 2
  )

  actual <- slice_derivation(
    advs,
    derivation = derive_vars_dtm,
    args = params(
      dtc = VSDTC,
      new_vars_prefix = "A",
      time_imputation = "last"
    ),
    derivation_slice(
      filter = str_detect(VSTPT, "PRE|BEFORE"),
      args = params(time_imputation = "first")
    ),
    derivation_slice(
      filter = TRUE
    )
  )

  expected <- mutate(
    advs,
    ADTM = c(ymd_hms("2020-04-16 23:59:59"), ymd_hms("2020-04-16 00:00:00")),
    ATMF = "H"
  )

  expect_dfs_equal(
    base = expected,
    compare = actual,
    keys = c("USUBJID", "VSSEQ")
  )
})
