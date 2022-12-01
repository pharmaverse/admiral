# get_constant_vars ----
## Test 1: without ignore_vars ----
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

## Test 2: with ignore_vars ----
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

# get_duplicates ----
## Test 3: x atomic vector ----
test_that("get_duplicates Test 3: x atomic vector", {
  x <- c("a", "a", "b", "c", "d", "d", 1, 1, 4)

  expect_equal(
    get_duplicates(x),
    c("a", "d", 1)
  )
})

# get_source_vars ----
## Test 4: x is a list of quosures ----
test_that("get_source_vars Test 4: x is a list of quosures", {
  x <- vars(DTHDOM = "AE", DTHSEQ = AESEQ)

  expect_equal(
    get_source_vars(x),
    x[2]
  )
})

## Test 5: quosures is NULL ----
test_that("get_source_vars Test 5: quosures is NULL", {
  expect_equal(
    get_source_vars(NULL),
    quo_c(NULL)
  )
})
