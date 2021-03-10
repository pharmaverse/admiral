context("test-derive_duration")


test_that("duration and unit variable are added", {
  input <- tibble::tribble(
    ~BRTHDT, ~RANDDT,
    ymd('1999-09-09'), ymd('2020-02-20'))

  expected_output <- input %>% mutate(AGE := 20, AGEU := 'YEARS')

  expect_equal(derive_duration(input,
                               new_col = AGE,
                               start_date = BRTHDT,
                               end_date = RANDDT,
                               unit_col = AGEU,
                               out_unit = 'years',
                               trunc_out = TRUE),
               expected_output)}
  )