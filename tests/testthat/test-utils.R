test_that("atomic vectors of length 1", {
  expect_identical(what_is_it(NULL), "`NULL`")
  expect_identical(what_is_it(TRUE), "`TRUE`")
  expect_identical(what_is_it(NA), "`NA`")
  expect_identical(what_is_it("Text"), '`"Text"`')
  expect_identical(what_is_it("3"), '`"3"`')
  expect_identical(what_is_it(4L), "`4`")
  expect_identical(what_is_it(2.42), "`2.42`")
})

test_that("vectors", {
  expect_identical(what_is_it(letters), "a character vector")
  expect_identical(what_is_it(1:10), "an integer vector")
  expect_identical(what_is_it(c(1.2, 3)), "a double vector")
  expect_identical(what_is_it(c(TRUE, FALSE)), "a logical vector")
  expect_identical(what_is_it(list(1, letters, TRUE)), "a list")
})

test_that("S3 objects", {
  expect_identical(what_is_it(mtcars), "a data frame")
  expect_identical(what_is_it(factor(letters)), "a factor")
  expect_identical(what_is_it(lm(hp ~ mpg, data = mtcars)), "an object of class 'lm'")
  expect_identical(what_is_it(quo(4 / 1)), "an object of class 'quosure'")
})


test_that("S4 objects", {
  expect_identical(what_is_it(lubridate::days(1)), "a S4 object of class 'Period'")
})

test_that("symbols", {
  expect_identical(what_is_it(quote(USUBJID)), "a symbol")
})

test_that("input is returned as is if filter is NULL", {
  input <- tibble::tribble(
    ~USUBJID, ~VSTESTCD, ~VSSTRESN,
    "P01", "WEIGHT", 80.9,
    "P01", "HEIGHT", 189.2
  )

  expect_dfs_equal(
    input,
    filter_if(input, quo(NULL)),
    keys = c("USUBJID", "VSTESTCD")
  )
})

test_that("input is filtered if filter is not NULL", {
  input <- tibble::tribble(
    ~USUBJID, ~VSTESTCD, ~VSSTRESN,
    "P01", "WEIGHT", 80.9,
    "P01", "HEIGHT", 189.2
  )

  expect_dfs_equal(
    input[1L, ],
    filter_if(input, quo(VSTESTCD == "WEIGHT")),
    keys = c("USUBJID", "VSTESTCD")
  )
})

test_that("enumerate works", {
  expect_equal(enumerate(letters[1]), "`a`")
  expect_equal(enumerate(letters[1:3]), "`a`, `b` and `c`")
})

test_that("squote works", {
  expect_equal(squote(letters[1]), "'a'")
  expect_equal(squote(letters[1:3]), c("'a'", "'b'", "'c'"))
})

test_that("arg_name works", {
  expect_equal(arg_name(sym("a")), "a")
  expect_equal(arg_name(call("enquo", sym("a"))), "a")
  expect_error(arg_name("a"), "Could not extract argument name from")
})

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

test_that("`get_one_to_many_dataset()` returns a data frame after a previous error", {
  try(assert_one_to_one(admiral_adsl, vars(STUDYID), vars(SITEID)), silent = TRUE)

  expect_true(is.data.frame(get_one_to_many_dataset()))
})

test_that("`get_many_to_one_dataset()` returns a data frame after a previous error", {
  try(assert_one_to_one(admiral_adsl, vars(SITEID), vars(STUDYID)), silent = TRUE)

  expect_true(is.data.frame(get_many_to_one_dataset()))
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

test_that("get_constant_vars Test 1: without ignore_vars", {
  data <- tibble::tribble(
    ~USUBJID, ~AGE, ~AVISIT,
    "1",      26,   "BASELINE",
    "1",      26,   "WEEK 1",
    "2",      42,   "BASELINE",
    "2",      42,   "WEEK 1"
  )

  expect_equal(
    get_constant_vars(data, by_vars = vars(USUBJID)),
    vars(USUBJID, AGE)
  )
})

test_that("get_constant_vars Test 2: with ignore_vars", {
  data <- tibble::tribble(
    ~USUBJID, ~AGE, ~WGTBL, ~HGTBL, ~AVISIT,
    "1",      26,   61,     172,    "BASELINE",
    "1",      26,   61,     172,    "WEEK 1",
    "2",      42,   72,     183,    "BASELINE",
    "2",      42,   72,     183,    "WEEK 1"
  )

  expect_equal(
    get_constant_vars(data, by_vars = vars(USUBJID), ignore_vars = vars(WGTBL, HGTBL)),
    vars(USUBJID, AGE)
  )
})
