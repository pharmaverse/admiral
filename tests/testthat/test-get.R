test_that("get_constant_vars Test 1: without ignore_vars", {
  data <- tibble::tribble(
    ~USUBJID, ~AGE, ~AVISIT,
    "1",      26,   "BASELINE",
    "1",      26,   "WEEK 1",
    "2",      42,   "BASELINE",
    "2",      42,   "WEEK 1"
  )

  expect_equal(
    get_constant_vars(data, by_vars = vars(USUBJID)),
    vars(USUBJID, AGE)
  )
})

test_that("get_constant_vars Test 2: with ignore_vars", {
  data <- tibble::tribble(
    ~USUBJID, ~AGE, ~WGTBL, ~HGTBL, ~AVISIT,
    "1",      26,   61,     172,    "BASELINE",
    "1",      26,   61,     172,    "WEEK 1",
    "2",      42,   72,     183,    "BASELINE",
    "2",      42,   72,     183,    "WEEK 1"
  )

  expect_equal(
    get_constant_vars(data, by_vars = vars(USUBJID), ignore_vars = vars(WGTBL, HGTBL)),
    vars(USUBJID, AGE)
  )
})
