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
  input <- structure(letters, names = rev(letters), label = "Letters") # nolint: undesirable_function_linter
  input[c(1, 9, 23)] <- NA
  output <- convert_blanks_to_na(input)

  expect_identical(attr(output, "label"), "Letters")
  expect_identical(names(output), rev(letters))
})

## Test 3: blank strings are turned into `NA` inside data frames ----
test_that("convert_blanks_to_na Test 3: blank strings are turned into `NA` inside data frames", {
  input <- tibble::tribble(
    ~a, ~b, ~c, ~d,
    "a", 1, TRUE, "",
    "b", NA, FALSE, "",
    "", 21, TRUE, "s",
    "c", 9, TRUE, "q"
  )
  attr(input$a, "label") <- "A"
  attr(input$b, "label") <- "B"
  attr(input$c, "label") <- "C"
  attr(input$d, "label") <- "D"

  expected_output <- tibble::tribble(
    ~a, ~b, ~c, ~d,
    "a", 1, TRUE, NA,
    "b", NA, FALSE, NA,
    NA, 21, TRUE, "s",
    "c", 9, TRUE, "q"
  )
  attr(expected_output$a, "label") <- "A"
  attr(expected_output$b, "label") <- "B"
  attr(expected_output$c, "label") <- "C"
  attr(expected_output$d, "label") <- "D"

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

# convert_na_to_blanks ----
## Test 5: `NA` strings are turned into blank  ----
test_that("convert_na_to_blanks Test 5: `NA` strings are turned into blank ", {
  expect_identical(
    convert_na_to_blanks(c("a", NA, "b")),
    c("a", "", "b")
  )
})

## Test 6: attributes are preserved when converting `NA` to blanks ----
test_that("convert_na_to_blanks Test 6: attributes are preserved when converting `NA` to blanks", {
  input <- structure(letters, names = rev(letters), label = "Letters") # nolint: undesirable_function_linter
  input[c(1, 9, 23)] <- NA_character_
  output <- convert_na_to_blanks(input)

  expect_identical(attr(output, "label"), "Letters")
  expect_identical(names(output), rev(letters))
})

# convert_na_to_blanks.data.frame ----
## Test 7: `NA` are turned into blank strings inside data frames ----
test_that("convert_na_to_blanks.data.frame Test 7: `NA` are turned into blank strings inside data frames", { # nolint
  input <- tibble::tribble(
    ~a, ~b, ~c, ~d,
    "a", 1, TRUE, NA,
    "b", NA, FALSE, NA,
    NA, 21, TRUE, "s",
    "c", 9, TRUE, "q"
  )
  attr(input$a, "label") <- "A"
  attr(input$b, "label") <- "B"
  attr(input$c, "label") <- "C"
  attr(input$d, "label") <- "D"

  expected_output <- tibble::tribble(
    ~a, ~b, ~c, ~d,
    "a", 1, TRUE, "",
    "b", NA, FALSE, "",
    "", 21, TRUE, "s",
    "c", 9, TRUE, "q"
  )
  attr(expected_output$a, "label") <- "A"
  attr(expected_output$b, "label") <- "B"
  attr(expected_output$c, "label") <- "C"
  attr(expected_output$d, "label") <- "D"

  expect_equal(convert_na_to_blanks.data.frame(input), expected_output, ignore_attr = TRUE)
})

# convert_na_to_blanks.list ----
## Test 8: `convert_na_to_blanks.list` produces a lists ----
test_that("convert_na_to_blanks.list Test 8: `convert_na_to_blanks.list` produces a lists", {
  x <- c(NA_character_, NA_character_, NA_character_)
  expected_output <- lapply(x, convert_na_to_blanks)
  actual_output <- convert_na_to_blanks.list(x)

  expect_equal(expected_output, actual_output)
})

# negate_vars ----
## Test 9: negate_vars returns list of negated variables ----
test_that("negate_vars Test 9: negate_vars returns list of negated variables", {
  expect_identical(negate_vars(exprs(var1, var2)), rlang::exprs(-var1, -var2))
})

## Test 10: negate_vars returns NULL if input is NULL ----
test_that("negate_vars Test 10: negate_vars returns NULL if input is NULL", {
  expect_identical(negate_vars(NULL), NULL)
})

# get_one_to_many_dataset ----
## Test 11: returns a data frame after a previous error ----
test_that("get_one_to_many_dataset Test 11: returns a data frame after a previous error", {
  try(assert_one_to_one(admiral_adsl, exprs(STUDYID), exprs(SITEID)), silent = TRUE)

  expect_true(is.data.frame(get_one_to_many_dataset()))
})

# get_many_to_one_dataset ----
## Test 12: returns a data frame after a previous error ----
test_that("get_many_to_one_dataset Test 12: returns a data frame after a previous error", {
  try(assert_one_to_one(admiral_adsl, exprs(SITEID), exprs(STUDYID)), silent = TRUE)

  expect_true(is.data.frame(get_many_to_one_dataset()))
})

# print.source ----
## Test 13: `source` objects are printed as intended ----
test_that("print.source Test 13: `source` objects are printed as intended", {
  ttae <- event_source(
    dataset_name = "ae",
    date = AESTDTC,
    set_values_to = exprs(
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
    "  SRCSEQ: AESEQ",
    "order: NULL"
  )
  expect_identical(capture.output(print(ttae)), expected_print_output)
})

## Test 14: `source` objects containing `source` objects ----
test_that("print.source Test 14: `source` objects containing `source` objects", {
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

## Test 15: `source` objects containing `data.frame` ----
test_that("print.source Test 15: `source` objects containing `data.frame`", {
  cqterms <- tibble::tribble(
    ~TERMCHAR, ~TERMNUM,
    "APPLICATION SITE ERYTHEMA", 10003041L,
    "APPLICATION SITE PRURITUS", 10003053L
  ) %>%
    mutate(SRCVAR = "AEDECOD")

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
    "  TERMCHAR                   TERMNUM SRCVAR ",
    "  <chr>                        <int> <chr>  ",
    "1 APPLICATION SITE ERYTHEMA 10003041 AEDECOD",
    "2 APPLICATION SITE PRURITUS 10003053 AEDECOD"
  )
  # replace × with x due to differences between R versions and remove formatting
  expect_identical(
    str_replace_all(capture.output(print(cq)), "×", "x") %>%
      str_replace_all("\033\\[[\\d;]+m", ""),
    expected_print_output
  )
})

# print_named_list ----
## Test 16: named list ----
test_that("print_named_list Test 16: named list", {
  expect_identical(
    capture.output(print_named_list(list(a = 1, b = 2))),
    c(
      "a: 1",
      "b: 2"
    )
  )
})

## Test 17: unnamed list ----
test_that("print_named_list Test 17: unnamed list", {
  expect_identical(
    capture.output(print_named_list(list(1, 2))),
    c(
      "1: 1",
      "2: 2"
    )
  )
})

## Test 18: named list with unamed list ----
test_that("print_named_list Test 18: named list with unamed list", {
  expect_snapshot(
    print_named_list(list(
      list_item = list("Hello World!", expr(universe), list(42)),
      another_one = ymd("2020-02-02")
    ))
  )
})
