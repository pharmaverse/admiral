# convert_blanks_to_na ----
## Test 1: blank strings are turned into `NA` ----
test_that("convert_blanks_to_na Test 1: blank strings are turned into `NA`", {
  expect_identical(
    convert_blanks_to_na(c("a", "", "b")),
    c("a", NA, "b")
  )
})

## Test 2: attributes are preserved when converting blanks to `NA` ----
test_that("convert_blanks_to_na Test 2: attributes are preserved when converting blanks to `NA`", {
  input <- structure(letters, names = rev(letters), label = "Letters")
  input[c(1, 9, 23)] <- NA
  output <- convert_blanks_to_na(input)

  expect_identical(attr(output, "label"), "Letters")
  expect_identical(names(output), rev(letters))
})

## Test 3: blank strings are turned into `NA` inside data frames ----
test_that("convert_blanks_to_na Test 3: blank strings are turned into `NA` inside data frames", {
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

# convert_blanks_to_na.list ----
## Test 4: `convert_blanks_to_na.list` produces a lists ----
test_that("convert_blanks_to_na.list Test 4: `convert_blanks_to_na.list` produces a lists", {
  x <- c("", "", "")
  expected_output <- lapply(x, convert_blanks_to_na)
  actual_output <- convert_blanks_to_na.list(x)

  expect_equal(expected_output, actual_output)
})

# negate_vars ----
## Test 5: negate_vars returns list of negated variables ----
test_that("negate_vars Test 5: negate_vars returns list of negated variables", {
  expect_identical(negate_vars(vars(var1, var2)), rlang::exprs(-var1, -var2))
})

## Test 6: negate_vars returns NULL if input is NULL ----
test_that("negate_vars Test 6: negate_vars returns NULL if input is NULL", {
  expect_identical(negate_vars(NULL), NULL)
})

# get_one_to_many_dataset ----
## Test 7: `get_one_to_many_dataset()` returns a data frame after a previous error ----
test_that("get_one_to_many_dataset Test 7: `get_one_to_many_dataset()` returns a data frame after a previous error", {
  try(assert_one_to_one(admiral_adsl, vars(STUDYID), vars(SITEID)), silent = TRUE)

  expect_true(is.data.frame(get_one_to_many_dataset()))
})

# get_many_to_one_dataset ----
## Test 8: `get_many_to_one_dataset()` returns a data frame after a previous error ----
test_that("get_many_to_one_dataset Test 8: `get_many_to_one_dataset()` returns a data frame after a previous error", {
  try(assert_one_to_one(admiral_adsl, vars(SITEID), vars(STUDYID)), silent = TRUE)

  expect_true(is.data.frame(get_many_to_one_dataset()))
})

# print.source ----
## Test 9: `source` objects are printed as intended ----
test_that("print.source Test 9: `source` objects are printed as intended", {
  ttae <- event_source(
    dataset_name = "ae",
    date = AESTDTC,
    set_values_to = vars(
      EVENTDESC = "AE",
      SRCDOM = "AE",
      SRCVAR = "AESTDTC",
      SRCSEQ = AESEQ
    )
  )
  expected_print_output <- c(
    "<event_source> object",
    "dataset_name: \"ae\"",
    "filter: NULL",
    "date: AESTDTC",
    "censor: 0",
    "set_values_to:",
    "  EVENTDESC: \"AE\"",
    "  SRCDOM: \"AE\"",
    "  SRCVAR: \"AESTDTC\"",
    "  SRCSEQ: AESEQ"
  )
  expect_identical(capture.output(print(ttae)), expected_print_output)
})

## Test 10: `source` objects containing `source` objects ----
test_that("print.source Test 10: `source` objects containing `source` objects", {
  slice <-
    derivation_slice(
      filter = AVISITN > 0,
      args = params(new_var = CHG)
    )
  expected_print_output <- c(
    "<derivation_slice> object",
    "filter: AVISITN > 0",
    "args:",
    "  <params> object",
    "  new_var: CHG"
  )
  expect_identical(capture.output(print(slice)), expected_print_output)
})

## Test 11: `source` objects containing `data.frame` ----
test_that("print.source Test 11: `source` objects containing `data.frame`", {
  cqterms <- tribble(
    ~TERM_NAME,                  ~TERM_ID,
    "APPLICATION SITE ERYTHEMA", 10003041L,
    "APPLICATION SITE PRURITUS", 10003053L
  ) %>%
    mutate(TERM_LEVEL = "AEDECOD")

  cq <- query(
    prefix = "CQ01",
    name = "Application Site Issues",
    definition = cqterms
  )
  expected_print_output <- c(
    "<query> object",
    "prefix: \"CQ01\"",
    "name: \"Application Site Issues\"",
    "add_scope_num: FALSE",
    "definition:",
    "# A tibble: 2 x 3",
    "  TERM_NAME                  TERM_ID TERM_LEVEL",
    "  <chr>                        <int> <chr>     ",
    "1 APPLICATION SITE ERYTHEMA 10003041 AEDECOD   ",
    "2 APPLICATION SITE PRURITUS 10003053 AEDECOD   "
  )
  expect_identical(capture.output(print(cq)), expected_print_output)
})
