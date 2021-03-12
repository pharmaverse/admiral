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
