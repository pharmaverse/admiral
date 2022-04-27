test_that("Test 1: Error is issued to function call if function errors", {
  expect_error(
    call_user_fun(compute_bmi(height = 172, weight = "hallo"))
  )
})
