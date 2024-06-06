# expect_dfs_equal ----
## Test 1: expect_dfs_equal works ----
test_that("expect_dfs_equal Test 1: expect_dfs_equal works", {
  a <- data.frame(x = 1:3, y = 4:6)
  b <- data.frame(x = 1:3, y = 5:7)

  expect_error(
    expect_dfs_equal(a, b, keys = "x"),
    regexp = "Differences found\\.*"
  )
})
