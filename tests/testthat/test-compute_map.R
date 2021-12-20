test_that("MAP based on diastolic and systolic blood pressure", {
  expect_equal(compute_map(
    diabp = 51,
    sysbp = 121
  ),
  (2 * 51 + 121) / 3)
})

test_that("MAP based on diastolic and systolic blood pressure", {
  expect_equal(compute_map(diabp = 51,
                           sysbp = 121,
                           hr = 59),
               51 + 0.01 * exp(4.14 - 40.74 / 59) * (121 - 51))
})
