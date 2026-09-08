# as_admiral_df ----
## Test 1: adds the admiral_df class to a data frame ----
test_that("as_admiral_df Test 1: adds the admiral_df class to a data frame", {
  input <- tibble::tribble(
    ~USUBJID, ~AVAL,
    "1",      10,
    "2",      20
  )

  actual <- as_admiral_df(input)

  expect_s3_class(actual, "admiral_df")
  expect_equal(class(actual), c("admiral_df", class(input)))
})

## Test 2: is idempotent when the class is already present ----
test_that("as_admiral_df Test 2: is idempotent when the class is already present", {
  input <- tibble::tribble(
    ~USUBJID, ~AVAL,
    "1",      10
  )

  once <- as_admiral_df(input)
  twice <- as_admiral_df(once)

  expect_equal(class(once), class(twice))
  expect_equal(sum(class(twice) == "admiral_df"), 1L)
})

## Test 3: preserves the existing classes of the data frame ----
test_that("as_admiral_df Test 3: preserves the existing classes of the data frame", {
  input <- tibble::tribble(
    ~USUBJID, ~AVAL,
    "1",      10
  )

  actual <- as_admiral_df(input)

  expect_true(all(c("tbl_df", "tbl", "data.frame") %in% class(actual)))
})

## Test 4: returns NULL unchanged ----
test_that("as_admiral_df Test 4: returns NULL unchanged", {
  expect_null(as_admiral_df(NULL))
})
