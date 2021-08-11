context("test-compute_bsa")

test_that("Mosteller method - single height and weight values", {
  expect_equal(compute_bsa(
    height = 170,
    weight = 75,
    method = "Mosteller"),
  1.88 )
})

test_that("Mosteller method - height and weight vectors", {
  expect_equal(compute_bsa(
    height = c(170, 185),
    weight = c(75, 90),
    method = "Mosteller"),
    c(1.88, 2.15))
})

test_that("DuBois-DuBois method - height and weight vectors", {
  expect_equal(compute_bsa(
    height = c(170, 185),
    weight = c(75, 90),
    method = "DuBois-DuBois"),
    c(1.86, 2.14))
})

test_that("Haycock method - height and weight vectors", {
  expect_equal(compute_bsa(
    height = c(170, 185),
    weight = c(75, 90),
    method = "Haycock"),
    c(1.89, 2.16))
})

test_that("Gehan-George - height and weight vectors", {
  expect_equal(compute_bsa(
    height = c(170, 185),
    weight = c(75, 90),
    method = "Gehan-George"),
    c(1.90, 2.16))
})

test_that("Boyd - height and weight vectors", {
  expect_equal(compute_bsa(
    height = c(170, 185),
    weight = c(75, 90),
    method = "Boyd"),
    c(1.91, 2.16))
})

test_that("Fujimoto - height and weight vectors", {
  expect_equal(compute_bsa(
    height = c(170, 185),
    weight = c(75, 90),
    method = "Fujimoto"),
    c(1.82, 2.09))
})

test_that("Takahira - height and weight vectors", {
  expect_equal(compute_bsa(
    height = c(170, 185),
    weight = c(75, 90),
    method = "Takahira"),
    c(1.88, 2.16))
})

