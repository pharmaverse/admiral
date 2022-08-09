

test_that("negate_vars returns list of negated variables", {
  expect_identical(negate_vars(vars(var1, var2)), rlang::exprs(-var1, -var2))
})

test_that("negate_vars returns NULL if input is NULL", {
  expect_identical(negate_vars(NULL), NULL)
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
