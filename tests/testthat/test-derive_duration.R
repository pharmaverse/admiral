context("test-derive_duration")


test_that("duration and unit variable are added", {
  input <- tibble::tribble(
    ~BIRTHDT, ~RANDDT,
    ymd('1999-09-09'), ymd('2020-02-20'))

  expected_output <- input %>% mutate(AGE := 20, AGEU := 'YEARS')

  expect_equal(derive_duration(input,
                               newcol = AGE,
                               startdate = BIRTHDT,
                               enddate = RANDDT,
                               unitcol = AGEU,
                               out_unit = 'years',
                               trunc_out = TRUE),
               expected_output)}
  )