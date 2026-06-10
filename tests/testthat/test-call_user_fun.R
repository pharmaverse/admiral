test_that("Test 1: deprecation error", {
  expect_snapshot(
    call_user_fun(compute_bmi(height = 172, weight = "hallo")),
    error = TRUE
  )
})
