# get_admiral_option ----
## Test 1: get works ----
test_that("get_admiral_option Test 1: get works", {
  expect_equal(get_admiral_option("subject_keys"), vars(STUDYID, USUBJID))
})

## Test 2: common typo gives error to select available options ----
test_that("get_admiral_option Test 2: common typo gives error to select available options", {
  expect_error(get_admiral_option("subject_key"))
})

## Test 3: non-character argument triggers assertion error ----
test_that("get_admiral_option Test 3: non-character argument triggers assertion error", {
  subject_keys <- 1
  expect_error(get_admiral_option(subject_keys), "`option` must be a character scalar but is `1`")
})

# set_admiral_options ----
## Test 4: set works ----
test_that("set_admiral_options Test 4: set works", {
  set_admiral_options(subject_keys = vars(STUDYID, USUBJID2))
  expect_equal(get_admiral_option("subject_keys"), vars(STUDYID, USUBJID2))
})

## Test 5: unexpected function input for set gives error ----
test_that("set_admiral_options Test 5: unexpected function input for set gives error", {
  expect_error(set_admiral_options(subject_keys = quo_c(STUDYID, USUBJID2)))
  expect_error(set_admiral_options(subject_keys = STUDYID))
})
set_admiral_options(subject_keys = vars(STUDYID, USUBJID))
