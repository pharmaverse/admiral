test_that("arg_name works", {
  expect_equal(arg_name(sym("a")), "a")
  expect_equal(arg_name(call("enquo", sym("a"))), "a")
  expect_error(arg_name("a"), "Could not extract argument name from")
})

test_that("`convert_dtm_to_dtc` is in correct format", {
  expect_equal(
    convert_dtm_to_dtc(as.POSIXct("2022-04-05 15:34:07 UTC")),
    "2022-04-05T15:34:07"
  )
})

test_that("`convert_dtm_to_dtc` Error is thrown if dtm is not in correct format", {
  expect_error(
    convert_dtm_to_dtc("2022-04-05T15:26:14"),
    "lubridate::is.instant(dtm) is not TRUE",
    fixed = TRUE
  )
})

test_that("input is returned as is if filter is NULL", {
  input <- tibble::tribble(
    ~USUBJID, ~VSTESTCD, ~VSSTRESN,
    "P01", "WEIGHT", 80.9,
    "P01", "HEIGHT", 189.2
  )

  expect_dfs_equal(
    input,
    filter_if(input, quo(NULL)),
    keys = c("USUBJID", "VSTESTCD")
  )
})

test_that("input is filtered if filter is not NULL", {
  input <- tibble::tribble(
    ~USUBJID, ~VSTESTCD, ~VSSTRESN,
    "P01", "WEIGHT", 80.9,
    "P01", "HEIGHT", 189.2
  )

  expect_dfs_equal(
    input[1L, ],
    filter_if(input, quo(VSTESTCD == "WEIGHT")),
    keys = c("USUBJID", "VSTESTCD")
  )
})
