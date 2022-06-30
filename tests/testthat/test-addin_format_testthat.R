# ---- addin_format_testthat, test 1: works as expected ----
test_that("addin_format_testthat, test 1: works as expected", {

  # test: file exists
  expect_error(
    prepare_test_that_file("file_does_not_exist"),
    "Invalid file path, the file does not exist."
  )

  # test: correct file is used
  tf <- tempfile()
  writeLines("Hello", tf)
  expect_error(
    prepare_test_that_file(tf),
    "This Addin works only on unit test files that follow a testthat structure"
  )

  # test: test_that not found
  dir.create(file.path(tempdir(), "tests", "testthat"), recursive = TRUE)
  tf <- tempfile(pattern = "tests/testthat/test-")
  writeLines("Hello", tf)
  expect_null(prepare_test_that_file(tf))

  # test: headers and descriptions works as expected on dummy file
  tf <- tempfile(pattern = "tests/testthat/test-")
  writeLines(
    c(
      "# some stuff in comment",
      'test_that("my description", {',
      "  expect_true(TRUE)",
      "}"
    ),
    tf
  )
  expected <- list(
    file_content = c(
      "# some stuff in comment",
      paste0(
        "# ---- ", sub("test-", "", basename(tf)), ", ",
        "test 1: my description ----"
      ),
      paste0(
        'test_that("', sub("test-", "", basename(tf)), ", ",
        'test 1: my description", {'
      ),
      "  expect_true(TRUE)",
      "}"
    )
  )
  expect_identical(prepare_test_that_file(tf), expected)
})
