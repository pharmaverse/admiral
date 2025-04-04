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

## Test 5: calling it in a function ----
test_that("slice_derivation Test 5: calling it in a function", {
  expected <- tibble::tribble(
    ~USUBJID, ~ARM, ~MVAL,
    "1",      "A",  "A_val",
    "2",      "B",  "B_val"
  )

  a_data <- tibble::tribble(
    ~USUBJID, ~AVALC,
    "1",      "A_val"
  )

  b_data <- tibble::tribble(
    ~USUBJID, ~AVALC,
    "2",      "B_val"
  )

  my_merge <- function(dataset, my_b_data) {
    slice_derivation(
      dataset,
      derivation = derive_vars_merged,
      args = params(
        by_vars = exprs(USUBJID),
        new_vars = exprs(MVAL = AVALC)
      ),
      derivation_slice(
        filter = ARM == "A",
        args = params(dataset_add = a_data)
      ),
      derivation_slice(
        filter = ARM == "B",
        args = params(dataset_add = my_b_data)
      )
    )
  }

  expect_dfs_equal(
    base = expected,
    compare = my_merge(
      dataset = select(expected, -MVAL),
      my_b_data = b_data
    ),
    keys = "USUBJID"
  )
})

## Test 6: slice on 0-row dataset ----
test_that("slice_derivation Test 6: slice on 0-row dataset", {
  advs <- tibble::tribble(
    ~USUBJID, ~VSDTC,       ~VSTPT,             ~VSSEQ,
    "1",      "2020-04-16", NA_character_,      1,
    "1",      "2020-04-16", "BEFORE TREATMENT", 2
  )

  actual <- slice_derivation(
    advs[c(-1, -2), ],
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
    advs[c(-1, -2), ],
    ADTM = as.POSIXct(numeric(0), tz = "UTC"),
    ATMF = character(0)
  )

  expect_identical(
    expected,
    actual
  )
})

## Test 7: Error thrown if a mandatory argument is not in arg or all slices----
test_that("slice_derivation Test 7: Error if a mandatory argument is not in arg or all slices", {
  advs <- tibble::tribble(
    ~USUBJID, ~VSDTC,       ~VSTPT,             ~VSSEQ,
    "1",      "2020-04-16", NA_character_,      1,
    "1",      "2020-04-16", "BEFORE TREATMENT", 2
  )

  expect_snapshot(
    actual <- slice_derivation(
      advs,
      derivation = derive_vars_dtm,
      args = params(
        dtc = VSDTC
      ),
      derivation_slice(
        filter = str_detect(VSTPT, "PRE|BEFORE"),
        args = params(time_imputation = "first", new_vars_prefix = "A")
      ),
      derivation_slice(
        filter = TRUE,
        args = params(time_imputation = "last")
      )
    ),
    error = TRUE
  )
})
