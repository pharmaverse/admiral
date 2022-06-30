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

  # test: headers and descriptions works as expected on dummy file
  dir.create(file.path(tempdir(), "tests", "testthat"), recursive = TRUE)
  tf <- tempfile(pattern = "tests/testthat/test-")#, tmpdir = "tests/testthat")
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
    descriptions = list(
      location = list(c(2, 1, 2, Inf)),
      text = paste0(
        'test_that("', sub("test-", "", basename(tf)), ", ",
        'test 1: my description", {'
      )
    ),
    headers_cond = TRUE,
    headers = list(
      location = list(c(2, 1)),
      text = paste0(
        "# ---- ", sub("test-", "", basename(tf)), ", ",
        "test 1: my description ----\n"
      )
    )
  )
  expect_identical(prepare_test_that_file(tf), expected)

})
