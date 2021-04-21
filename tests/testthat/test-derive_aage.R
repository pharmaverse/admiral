context("test-derive_aage")


test_that("duration and unit variable are added", {
  input <- tibble::tribble(
    ~BRTHDT, ~RANDDT,
    ymd("1999-09-09"), ymd("2020-02-20")
  )
  expected_output <- mutate(input, AAGE = 20, AAGEU = "YEARS")

  expect_equal(derive_aage(input), expected_output)
})
