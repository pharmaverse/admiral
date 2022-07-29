
test_that("blank strings are turned into `NA`", {
  expect_identical(
    convert_blanks_to_na(c("a", "", "b")),
    c("a", NA, "b")
  )
})

test_that("attributes are preserved when converting blanks to `NA`", {
  input <- structure(letters, names = rev(letters), label = "Letters")
  input[c(1, 9, 23)] <- NA
  output <- convert_blanks_to_na(input)

  expect_identical(attr(output, "label"), "Letters")
  expect_identical(names(output), rev(letters))
})

test_that("blank strings are turned into `NA` inside data frames", {
  input <- tibble::tibble(
    a = structure(c("a", "b", "", "c"), label = "A"),
    b = structure(c(1, NA, 21, 9), label = "B"),
    c = structure(c(TRUE, FALSE, TRUE, TRUE), label = "C"),
    d = structure(c("", "", "s", "q"), label = "D")
  )
  expected_output <- tibble::tibble(
    a = structure(c("a", "b", NA, "c"), label = "A"),
    b = structure(c(1, NA, 21, 9), label = "B"),
    c = structure(c(TRUE, FALSE, TRUE, TRUE), label = "C"),
    d = structure(c(NA, NA, "s", "q"), label = "D")
  )

  expect_identical(convert_blanks_to_na(input), expected_output)
})

test_that("negate_vars returns list of negated variables", {
  expect_identical(negate_vars(vars(var1, var2)), rlang::exprs(-var1, -var2))
})

test_that("negate_vars returns NULL if input is NULL", {
  expect_identical(negate_vars(NULL), NULL)
})


test_that("`convert_dtm_to_dtc` is in correct format", {
  expect_equal(
    convert_dtm_to_dtc(as.POSIXct("2022-04-05 15:34:07 UTC")),
    "2022-04-05T15:34:07"
  )
})

test_that("`convert_dtm_to_dtc` Error is thrown if dtm is not in correct format", {
  expect_error(
    convert_dtm_to_dtc("2022-04-05T15:26:14"),
    "lubridate::is.instant(dtm) is not TRUE",
    fixed = TRUE
  )
})
