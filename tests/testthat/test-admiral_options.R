# get_admiral_option ----
## Test 1: get works ----
test_that("get_admiral_option Test 1: get works", {
  expect_equal(get_admiral_option("subject_keys"), exprs(STUDYID, USUBJID))
})

## Test 2: common typo gives error to select available options ----
test_that("get_admiral_option Test 2: common typo gives error to select available options", {
  expect_error(get_admiral_option("subject_key"))
})

## Test 3: non-character argument triggers assertion error ----
test_that("get_admiral_option Test 3: non-character argument triggers assertion error", {
  subject_keys <- 1
  expect_error(get_admiral_option(subject_keys),
    class = "assert_character_scalar"
  )
})

# set_admiral_options ----
## Test 4: set works ----
test_that("set_admiral_options Test 4: set works", {
  set_admiral_options(subject_keys = exprs(STUDYID, USUBJID2))
  expect_equal(get_admiral_option("subject_keys"), exprs(STUDYID, USUBJID2))
})

## Test 5: unexpected function input for set gives error ----
test_that("set_admiral_options Test 5: unexpected function input for set gives error", {
  expect_error(
    set_admiral_options(subject_keys = rlang::quos(STUDYID, USUBJID2)),
    class = "assert-admiraldev"
  )
  expect_error(
    set_admiral_options(subject_keys = STUDYID),
    class = "assert-admiraldev"
  )
})

## Test 6: non-positive numeric argument triggers assertion error ----
test_that("get_admiral_option Test 6: non-positive numeric argument triggers assertion error", {
  sigfigs <- 0
  expect_error(
    set_admiral_options(signif_digits = sigfigs),
    class = "assert_integer_scalar"
  )
  sigfigs <- -1
  expect_error(
    set_admiral_options(signif_digits = sigfigs),
    class = "assert_integer_scalar"
  )
})

## Test 7: set works ----
test_that("set_admiral_options Test 7: set works", {
  sigfigs <- get_admiral_option("signif_digits")
  set_admiral_options(signif_digits = sigfigs + 1)
  expect_equal(get_admiral_option("signif_digits"), sigfigs + 1)
  set_admiral_options(signif_digits = sigfigs)
})
set_admiral_options(subject_keys = exprs(STUDYID, USUBJID))
