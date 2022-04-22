# derive_vars_dy Test 1----
  test_that("derive_vars_dy Test1: Single --DT input when ref date is --DTM", {
    datain <- tibble::tribble(
      ~STUDYID, ~USUBJID, ~TRTSDTM, ~ASTDT,
      "TEST01", "PAT01", "2014-01-17T23:59:59", "2014-01-18",

    ) %>%
      mutate(TRTSDTM = lubridate::as_datetime(TRTSDTM),
             ASTDT = lubridate::ymd(ASTDT))

    expected_output <- tibble::tribble(
      ~STUDYID, ~USUBJID, ~TRTSDTM, ~ASTDT,  ~ASTDY,
      "TEST01", "PAT01", "2014-01-17T23:59:59", "2014-01-18", 2

    ) %>%
      mutate(TRTSDTM = lubridate::as_datetime(TRTSDTM),
             ASTDT = lubridate::ymd(ASTDT))

      actual_output <- derive_vars_dy(datain,
                                            reference_date = TRTSDTM,
                                            source_vars = vars(ASTDT))

            expect_dfs_equal(
              expected_output,
              actual_output,
              keys = c("STUDYID", "USUBJID")
            )
          })

# derive_vars_dy Test 2----
test_that("derive_vars_dy Test2: Multiple --DT input when ref date is --DTM,
          with name of --DY var specified in par call", {
  datain <- tibble::tribble(
    ~STUDYID, ~USUBJID, ~TRTSDTM, ~ASTDT, ~AENDT, ~DTHDT,
    "TEST01", "PAT01", "2014-01-17T23:59:59", "2014-01-18", "2014-01-20", "2014-02-01",
    "TEST01", "PAT02", "2014-01-17T23:59:59", "2014-01-17", "2014-01-20", "2014-03-01",
    "TEST01", "PAT03", "2014-01-17T23:59:59", "2014-01-15", "2014-01-20", "2014-05-01"

  ) %>%
    mutate(TRTSDTM = lubridate::as_datetime(TRTSDTM),
           ASTDT = lubridate::ymd(ASTDT),
           AENDT = lubridate::ymd(AENDT),
           DTHDT = lubridate::ymd(DTHDT))

  expected_output <- tibble::tribble(
  ~STUDYID, ~USUBJID, ~TRTSDTM, ~ASTDT, ~AENDT,  ~DTHDT, ~TRTSDY, ~ASTDY, ~AENDY, ~DEATHDY,
 "TEST01", "PAT01", "2014-01-17T23:59:59", "2014-01-18", "2014-01-20", "2014-02-01", 1, 2, 4, 16,
 "TEST01", "PAT02", "2014-01-17T23:59:59", "2014-01-17", "2014-01-20", "2014-03-01", 1, 1, 4, 44,
 "TEST01", "PAT03", "2014-01-17T23:59:59", "2014-01-15", "2014-01-20", "2014-05-01", 1, -2, 4, 105,

  ) %>%
    mutate(TRTSDTM = lubridate::as_datetime(TRTSDTM),
           ASTDT = lubridate::ymd(ASTDT),
           AENDT = lubridate::ymd(AENDT),
           DTHDT = lubridate::ymd(DTHDT))

  actual_output <- derive_vars_dy(datain,
                                  reference_date = TRTSDTM,
                                  source_vars = vars(TRTSDTM, ASTDT, AENDT, DEATHDY = DTHDT))

  expect_dfs_equal(
    expected_output,
    actual_output,
    keys = c("STUDYID", "USUBJID")
  )
})

# derive_vars_dy Test 3----
test_that("derive_vars_dy Test3: Combo of --DT/--DTM input when ref date is --DTM", {
  datain <- tibble::tribble(
    ~STUDYID, ~USUBJID, ~TRTSDTM, ~ASTDTM, ~AENDT,
    "TEST01", "PAT01", "2014-01-17T23:59:59", "2014-01-18T13:09:O9", "2014-01-20"

  ) %>%
    mutate(TRTSDTM = lubridate::as_datetime(TRTSDTM),
           ASTDTM = lubridate::as_datetime(ASTDTM),
           AENDT = lubridate::ymd(AENDT))

  expected_output <- tibble::tribble(
    ~STUDYID, ~USUBJID, ~TRTSDTM, ~ASTDTM, ~AENDT,  ~TRTSDY, ~ASTDY, ~AENDY,
    "TEST01", "PAT01", "2014-01-17T23:59:59", "2014-01-18T13:09:O9", "2014-01-20", 1, 2, 4

  ) %>%
    mutate(TRTSDTM = lubridate::as_datetime(TRTSDTM),
           ASTDTM = lubridate::as_datetime(ASTDTM),
           AENDT = lubridate::ymd(AENDT))

  actual_output <- derive_vars_dy(datain,
                                  reference_date = TRTSDTM,
                                  source_vars = vars(TRTSDTM, ASTDTM, AENDT))

  expect_dfs_equal(
    expected_output,
    actual_output,
    keys = c("STUDYID", "USUBJID")
  )
})

#derive_vars_dy Test 4----
test_that("derive_vars_dy Test4: Single --DT input when ref date is --DT", {
  datain <- tibble::tribble(
    ~STUDYID, ~USUBJID, ~TRTSDT, ~ASTDT,
    "TEST01", "PAT01", "2014-01-17", "2014-01-18",

  ) %>%
    mutate(TRTSDT = lubridate::ymd(TRTSDT),
           ASTDT = lubridate::ymd(ASTDT))

  expected_output <- tibble::tribble(
    ~STUDYID, ~USUBJID, ~TRTSDT, ~ASTDT,  ~ASTDY,
    "TEST01", "PAT01", "2014-01-17", "2014-01-18", 2

  ) %>%
    mutate(TRTSDT = lubridate::ymd(TRTSDT),
           ASTDT = lubridate::ymd(ASTDT))

  actual_output <- derive_vars_dy(datain,
                                  reference_date = TRTSDT,
                                  source_vars = vars(ASTDT))

  expect_dfs_equal(
    expected_output,
    actual_output,
    keys = c("STUDYID", "USUBJID")
  )
})

# derive_vars_dy Test 5----
test_that("derive_vars_dy Test5: Multiple --DT input when ref date is --DT", {
  datain <- tibble::tribble(
    ~STUDYID, ~USUBJID, ~TRTSDT, ~ASTDT, ~AENDT,
    "TEST01", "PAT01", "2014-01-17", "2014-01-18", "2014-01-20"

  ) %>%
    mutate(TRTSDT = lubridate::ymd(TRTSDT),
           ASTDT = lubridate::ymd(ASTDT),
           AENDT = lubridate::ymd(AENDT))

  expected_output <- tibble::tribble(
    ~STUDYID, ~USUBJID, ~TRTSDT, ~ASTDT, ~AENDT,  ~TRTSDY, ~ASTDY, ~AENDY,
    "TEST01", "PAT01", "2014-01-17", "2014-01-18", "2014-01-20", 1, 2, 4

  ) %>%
    mutate(TRTSDT = lubridate::ymd(TRTSDT),
           ASTDT = lubridate::ymd(ASTDT),
           AENDT = lubridate::ymd(AENDT))

  actual_output <- derive_vars_dy(datain,
                                  reference_date = TRTSDT,
                                  source_vars = vars(TRTSDT, ASTDT, AENDT))

  expect_dfs_equal(
    expected_output,
    actual_output,
    keys = c("STUDYID", "USUBJID")
  )
})

# derive_vars_dy Test 6----
test_that("derive_vars_dy Test6: Combo of --DT/--DTM input when ref date is --DT", {
  datain <- tibble::tribble(
    ~STUDYID, ~USUBJID, ~TRTSDT, ~ASTDTM, ~AENDT,
    "TEST01", "PAT01", "2014-01-17", "2014-01-18T13:09:O9", "2014-01-20"

  ) %>%
    mutate(TRTSDT = lubridate::ymd(TRTSDT),
           ASTDTM = lubridate::as_datetime(ASTDTM),
           AENDT = lubridate::ymd(AENDT))

  expected_output <- tibble::tribble(
    ~STUDYID, ~USUBJID, ~TRTSDT, ~ASTDTM, ~AENDT,  ~TRTSDY, ~ASTDY, ~AENDY,
    "TEST01", "PAT01", "2014-01-17", "2014-01-18T13:09:O9", "2014-01-20", 1, 2, 4

  ) %>%
    mutate(TRTSDT = lubridate::ymd(TRTSDT),
           ASTDTM = lubridate::as_datetime(ASTDTM),
           AENDT = lubridate::ymd(AENDT))

  actual_output <- derive_vars_dy(datain,
                                  reference_date = TRTSDT,
                                  source_vars = vars(TRTSDT, ASTDTM, AENDT))

  expect_dfs_equal(
    expected_output,
    actual_output,
    keys = c("STUDYID", "USUBJID")
  )
})

# derive_vars_dy Test 7----
test_that("derive_vars_dy Test7: All dates as --DTM", {
  datain <- tibble::tribble(
    ~STUDYID, ~USUBJID, ~TRTSDTM, ~ASTDTM, ~AENDTM,
    "TEST01", "PAT01", "2014-01-17T16:34:O9", "2014-01-18T13:09:O9", "2014-01-20T08:29:05"

  ) %>%
    mutate(TRTSDTM = lubridate::as_datetime(TRTSDTM),
           ASTDTM = lubridate::as_datetime(ASTDTM),
           AENDTM = lubridate::as_datetime(AENDTM))

  expected_output <- tibble::tribble(
    ~STUDYID, ~USUBJID, ~TRTSDTM, ~ASTDTM, ~AENDTM,  ~TRTSDY, ~ASTDY, ~AENDY,
    "TEST01", "PAT01", "2014-01-17T16:34:O9", "2014-01-18T13:09:O9", "2014-01-20T08:29:05", 1, 2, 4

  ) %>%
    mutate(TRTSDTM = lubridate::as_datetime(TRTSDTM),
           ASTDTM = lubridate::as_datetime(ASTDTM),
           AENDTM = lubridate::as_datetime(AENDTM))

  actual_output <- derive_vars_dy(datain,
                                  reference_date = TRTSDTM,
                                  source_vars = vars(TRTSDTM, ASTDTM, AENDTM))

  expect_dfs_equal(
    expected_output,
    actual_output,
    keys = c("STUDYID", "USUBJID")
  )
})

# derive_vars_dy Test 8----
test_that("derive_vars_dy Test8: An error is issued if source variables do not end in DT or DTM", { # nolint
    datain <- tibble::tribble(
      ~STUDYID, ~USUBJID, ~TRTSDTW, ~ASTDTW, ~AENDTW,
      "TEST01", "PAT01", "2014-01-17T16:34:O9", "2014-01-18T13:09:O9", "2014-01-20T08:29:05"
    )

  expect_error(
    derive_vars_dy(datain,
                   reference_date = TRTSDTW,
                   source_vars = vars(TRTSDTW, ASTDTW, AENDTW)),
    "source_vars must end in DT or DTM or be explicitly and uniquely named.\nPlease name or rename the following source_vars:\nTRTSDTW, ASTDTW, AENDTW" # nolint
  )
})
