# get_admiral_options ----
## Test 1: get works ----
test_that("get_admiral_options,Test 1: get works", {
  expect_equal(get_admiral_options(subject_keys), vars(STUDYID, USUBJID))
})

## Test 2: common typo gives error ----
test_that("get_admiral_options, Test 2: common typo gives error", {
  expect_error(get_admiral_options(subject_key))
})

## Test 3: unexpected function input for get gives error ----
test_that("get_admiral_options, Test 3: unexpected function input for get gives error",{
  expect_error(get_admiral_options("subject_key"))
})

# set_admiral_options ----
## Test 4: set works ----
# test_that("set_admiral_options, Test 4: set works", {
#   set_admiral_options(subject_keys = vars(STUDYID, USUBJID2))
#   expect_equal(admiral_options$subject_keys, vars(STUDYID, USUBJID2))
# })

## Test 5: unexpected function input for set gives error ----
test_that("set_admiral_options, Test 5: unexpected function input for set gives error", {
  expect_error(set_admiral_options(subject_keys = quo_c(STUDYID, USUBJID2)))
  expect_error(set_admiral_options(subject_keys = STUDYID))
})
