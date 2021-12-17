test_that("IDVAR is missing, join by USUBJID", {
  input <- tibble::tribble(
    ~USUBJID, ~DOMAIN,
    1, "DM",
    2, "DM",
    3, "DM",
    4, "DM"
  )

  input_supp <- tibble::tribble(
    ~USUBJID, ~RDOMAIN, ~IDVAR, ~IDVARVAL, ~QNAM, ~QLABEL, ~QVAL,
    1, "DM", "", "", "TEST1", "Test QNAM1", "Response 1",
    1, "DM", "", "", "TEST2", "Test QNAM2", "Response 2",
    2, "DM", "", "", "TEST1", "Test QNAM1", "Response 1",
    3, "DM", "", "", "TEST1", "Test QNAM1", "Response 1"
  )

  expected_output <- tibble::tribble(
    ~USUBJID, ~DOMAIN, ~TEST1, ~TEST2,
    1, "DM", "Response 1", "Response 2",
    2, "DM", "Response 1", NA,
    3, "DM", "Response 1", NA,
    4, "DM", NA, NA
  )
  attr(expected_output$TEST1, "label") <- "Test QNAM1"
  attr(expected_output$TEST2, "label") <- "Test QNAM2"

  actual_output <- derive_vars_suppqual(
    input,
    input_supp
  )

  expect_dfs_equal(
    base = expected_output,
    compare = actual_output,
    keys = "USUBJID"
  )
})

test_that("Multiple IDVARs, differing types", {
  input <- tibble::tribble(
    ~USUBJID, ~DOMAIN, ~NIDVAR, ~CIDVAR,
    1, "LB", 1, "A",
    1, "LB", 2, "B",
    2, "LB", 1, "A",
    3, "LB", 1, "A",
    3, "LB", 2, "A",
    4, "LB", 2, "B"
  )

  input_supp <- tibble::tribble(
    ~USUBJID, ~RDOMAIN, ~IDVAR, ~IDVARVAL, ~QNAM, ~QLABEL, ~QVAL,
    1, "LB", "NIDVAR", "1", "TEST1", "Test QNAM1", "Response 1.1",
    1, "LB", "NIDVAR", "2", "TEST1", "Test QNAM1", "Response 1.2",
    1, "LB", "CIDVAR", "B", "TEST2", "Test QNAM2", "Response 2B",
    2, "LB", "NIDVAR", "1", "TEST1", "Test QNAM1", "Response 1.1",
    3, "LB", "NIDVAR", "1", "TEST1", "Test QNAM1", "Response 1.1",
    3, "LB", "NIDVAR", "2", "TEST1", "Test QNAM1", "Response 2.2"
  )

  expected_output <- tibble::tribble(
    ~USUBJID, ~DOMAIN, ~NIDVAR, ~CIDVAR, ~TEST1, ~TEST2,
    1, "LB", 1, "A", "Response 1.1", NA_character_,
    1, "LB", 2, "B", "Response 1.2", "Response 2B",
    2, "LB", 1, "A", "Response 1.1", NA_character_,
    3, "LB", 1, "A", "Response 1.1", NA_character_,
    3, "LB", 2, "A", "Response 2.2", NA_character_,
    4, "LB", 2, "B", NA_character_, NA_character_,
  )
  attr(expected_output$TEST1, "label") <- "Test QNAM1"
  attr(expected_output$TEST2, "label") <- "Test QNAM2"

  actual_output <- derive_vars_suppqual(
    input,
    input_supp
  )

  expect_dfs_equal(
    base = expected_output,
    compare = actual_output,
    keys = c("USUBJID", "NIDVAR", "CIDVAR")
  )
})


test_that("Multiple Records for each IDVAR", {
  input <- tibble::tribble(
    ~USUBJID, ~DOMAIN, ~NIDVAR, ~CIDVAR,
    1, "LB", 1, "A",
    1, "LB", 2, "B",
    1, "LB", 3, "B",
    2, "LB", 1, "A"
  )

  input_supp <- tibble::tribble(
    ~USUBJID, ~RDOMAIN, ~IDVAR, ~IDVARVAL, ~QNAM, ~QLABEL, ~QVAL,
    1, "LB", "NIDVAR", "1", "TEST1", "Test QNAM1", "Response 1.1",
    1, "LB", "CIDVAR", "B", "TEST1", "Test QNAM1", "Response 1.2"
  )

  expected_output <- tibble::tribble(
    ~USUBJID, ~DOMAIN, ~NIDVAR, ~CIDVAR, ~TEST1,
    1, "LB", 1, "A", "Response 1.1",
    1, "LB", 2, "B", "Response 1.2",
    1, "LB", 3, "B", "Response 1.2",
    2, "LB", 1, "A", NA_character_
  )
  attr(expected_output$TEST1, "label") <- "Test QNAM1"

  actual_output <- derive_vars_suppqual(
    input,
    input_supp
  )

  expect_dfs_equal(
    base = expected_output,
    compare = actual_output,
    keys = c("USUBJID", "NIDVAR", "CIDVAR")
  )
})


test_that("Test domain paramter", {
  input <- tibble::tribble(
    ~USUBJID, ~DOMAIN,
    1, "DM",
    2, "DM",
    3, "DM",
    4, "DM"
  )

  input_supp <- tibble::tribble(
    ~USUBJID, ~RDOMAIN, ~IDVAR, ~IDVARVAL, ~QNAM, ~QLABEL, ~QVAL,
    1, "DM", "", "", "TEST1", "Test QNAM1", "Response 1",
    1, "AE", "", "", "TEST2", "Test QNAM2", "Response 2",
    2, "CM", "", "", "TEST1", "Test QNAM1", "Response 1",
    3, "DM", "", "", "TEST1", "Test QNAM1", "Response 1"
  )

  expected_output <- tibble::tribble(
    ~USUBJID, ~DOMAIN, ~TEST1,
    1, "DM", "Response 1",
    2, "DM", NA,
    3, "DM", "Response 1",
    4, "DM", NA
  )
  attr(expected_output$TEST1, "label") <- "Test QNAM1"

  actual_output <- derive_vars_suppqual(
    input,
    input_supp,
    domain = "DM"
  )

  expect_dfs_equal(
    base = expected_output,
    compare = actual_output,
    keys = "USUBJID"
  )
})

test_that("Errors", {
  input <- tibble::tribble(
    ~USUBJID, ~DOMAIN,
    1, "DM",
  )

  input_supp <- tibble::tribble(
    ~USUBJID, ~RDOMAIN, ~IDVAR, ~IDVARVAL, ~QNAM, ~QLABEL, ~QVAL,
    1, "DM", "", "", "TEST1", "Test QNAM1", "Response 1",
    1, "AE", "", "", "TEST2", "Test QNAM2", "Response 2",
    2, "CM", "", "", "TEST1", "Test QNAM1", "Response 1",
    3, "DM", "", "", "TEST1", "Test QNAM1", "Response 1"
  )

  expect_error(
    derive_vars_suppqual(
      input,
      input_supp,
      domain = "DS"
      ),
    "Can't find the domain `DS` in `dataset_suppqual`."
  )

  input <- tibble::tribble(
    ~USUBJID, ~DOMAIN,
    1, "DS",
  )

  input_supp <- tibble::tribble(
    ~USUBJID, ~RDOMAIN, ~IDVAR, ~IDVARVAL, ~QNAM, ~QLABEL, ~QVAL,
    1, "DM", "", "", "TEST1", "Test QNAM1", "Response 1"
  )

  expect_error(
    derive_vars_suppqual(
      input,
      input_supp
    ),
    "DOMAIN of `dataset` and RDOMAIN of `dataset_suppqual` do not match."
  )
})
