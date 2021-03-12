context("test-derive_var_astdy")


test_that("ASTDY is added", {
  input <- tibble::tribble(
    ~TRTSDT, ~ASTDT,
    ymd('2020-01-01'), ymd('2020-02-24'),
    ymd('2020-01-01'), ymd('2020-01-01'),
    ymd('2020-02-24'), ymd('2020-01-01'))

  expected_output <- input %>% mutate(ASTDY := c(55, 1, -54))

  expect_equal(derive_var_astdy(input),
               expected_output)}
)