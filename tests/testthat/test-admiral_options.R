# get_admiral_option ----
## Test 1: get works ----
test_that("get_admiral_option Test 1: get works", {
  expect_equal(get_admiral_option("subject_keys"), exprs(STUDYID, USUBJID))
})

## Test 2: common typo gives error to select available options ----
test_that("get_admiral_option Test 2: common typo gives error to select available options", {
  expect_snapshot(
    get_admiral_option("subject_key"),
    error = TRUE
  )
})

## Test 3: non-character argument triggers assertion error ----
test_that("get_admiral_option Test 3: non-character argument triggers assertion error", {
  subject_keys <- 1
  expect_error(
    get_admiral_option(subject_keys),
    class = "assert_character_scalar"
  )
})

# set_admiral_options ----
## Test 4: setting subject_keys ----
test_that("set_admiral_options Test 4: setting subject_keys", {
  subject_keys <- get_admiral_option("subject_keys")
  set_admiral_options(subject_keys = exprs(STUDYID, USUBJID2))
  expect_equal(get_admiral_option("subject_keys"), exprs(STUDYID, USUBJID2))
  set_admiral_options(subject_keys = subject_keys)
})

## Test 5: error if subject_keys is invalid ----
test_that("set_admiral_options Test 5: error if subject_keys is invalid", {
  expect_error(
    set_admiral_options(subject_keys = rlang::quos(STUDYID, USUBJID2)),
    class = "assert_vars"
  )
  expect_error(
    set_admiral_options(subject_keys = STUDYID),
    class = "assert_vars"
  )
})

## Test 6: error if signif_digits is invalid ----
test_that("set_admiral_options Test 6: error if signif_digits is invalid", {
  expect_error(
    set_admiral_options(signif_digits = 0),
    class = "assert_integer_scalar"
  )
  expect_error(
    set_admiral_options(signif_digits = -1),
    class = "assert_integer_scalar"
  )
})

## Test 7: setting signif_digits ----
test_that("set_admiral_options Test 7: setting signif_digits", {
  sigfigs <- get_admiral_option("signif_digits")
  set_admiral_options(signif_digits = sigfigs + 1)
  expect_equal(get_admiral_option("signif_digits"), sigfigs + 1)
  set_admiral_options(signif_digits = sigfigs)
})

## Test 8: setting save_memory ----
test_that("set_admiral_options Test 8: setting save_memory", {
  save_memory <- get_admiral_option("save_memory")
  set_admiral_options(save_memory = TRUE)
  expect_equal(get_admiral_option("save_memory"), TRUE)
  set_admiral_options(save_memory = save_memory)
})

## Test 9: error if save_memory is invalid ----
test_that("set_admiral_options Test 9: error if save_memory is invalid", {
  expect_error(
    set_admiral_options(save_memory = 0),
    class = "assert_logical_scalar"
  )
})
