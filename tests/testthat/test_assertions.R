library(cdiscpilot)

test_that("an error is thrown if a required variable is missing", {
  data(dm)

  expect_error(
    assert_has_variables(dm, "TRT01P"),
    "Required variable `TRT01P` is missing."
  )
})

test_that("no error is thrown if a required variable exists", {
  data(dm)

  expect_error(assert_has_variables(dm, "USUBJID"), NA)
})

test_that("assert_filter_cond works as expected", {
  fc <- quo(AGE == 64)
  expect_identical(
    assert_filter_cond(fc),
    fc
  )

  fc <- quo()
  expect_error(
    assert_filter_cond(arg = fc),
    "Argument `fc` is missing, with no default"
  )

  expect_identical(
    assert_filter_cond(arg = fc, optional = TRUE),
    fc
  )

  fc <- quo("string")
  expect_error(
    assert_filter_cond(arg = fc),
    "`fc` must be a filter condition but is `\"string\"`"
  )
})
