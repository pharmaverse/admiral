context("test-compute_bsa")

test_that("Mosteller method - single height and weight values", {
  expect_equal(
    compute_bsa(
      height = 170,
      weight = 75,
      method = "Mosteller"
    ),
    sqrt(170 * 75 / 3600)
  )
})

test_that("Mosteller method - height and weight vectors", {
  expect_equal(
    compute_bsa(
      height = c(170, 185),
      weight = c(75, 90),
      method = "Mosteller"
    ),
    sqrt(c(170, 185) * c(75, 90) / 3600)
  )
})

test_that("DuBois-DuBois method - height and weight vectors", {
  expect_equal(
    compute_bsa(
      height = c(170, 185),
      weight = c(75, 90),
      method = "DuBois-DuBois"
    ),
    0.20247 * (c(170, 185) / 100)^0.725 * c(75, 90)^0.425
  )
})

test_that("Haycock method - height and weight vectors", {
  expect_equal(
    compute_bsa(
      height = c(170, 185),
      weight = c(75, 90),
      method = "Haycock"
    ),
    0.024265 * c(170, 185)^0.3964 * c(75, 90)^0.5378
  )
})

test_that("Gehan-George - height and weight vectors", {
  expect_equal(
    compute_bsa(
      height = c(170, 185),
      weight = c(75, 90),
      method = "Gehan-George"
    ),
    0.0235 * c(170, 185)^0.42246 * c(75, 90)^0.51456
  )
})

test_that("Boyd - height and weight vectors", {
  expect_equal(
    compute_bsa(
      height = c(170, 185),
      weight = c(75, 90),
      method = "Boyd"
    ),
    0.0003207 * (c(170, 185)^0.3) *
      (1000 * c(75, 90)) ^ (0.7285 - (0.0188 * log10(1000 * c(75, 90))))
  )
})

test_that("Fujimoto - height and weight vectors", {
  expect_equal(
    compute_bsa(
      height = c(170, 185),
      weight = c(75, 90),
      method = "Fujimoto"
    ),
    0.008883 * c(170, 185)^0.663 * c(75, 90)^0.444
  )
})

test_that("Takahira - height and weight vectors", {
  expect_equal(
    compute_bsa(
      height = c(170, 185),
      weight = c(75, 90),
      method = "Takahira"
    ),
    0.007241 * c(170, 185)^0.725 * c(75, 90)^0.425
  )
})

test_that("an error is issued if an invalid method is specified", {
  expect_error(
    compute_bsa(
      height = c(170, 185),
      weight = c(75, 90),
      method = "unknown-method"
    ),
    paste(
      "`method` must be one of 'Mosteller', 'DuBois-DuBois', 'Haycock', 'Gehan-George',",
      "'Boyd', 'Fujimoto' or 'Takahira' but is 'unknown-method'"
    )
  )
})
