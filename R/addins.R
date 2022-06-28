
# Function for the RStudio addin, see inst/rstudio/addins.dcf.
# Updates current test_that test file by adding a function name, a test number, and a section
format_test_that_file <- function() {

  # get file info
  file_info <- rstudioapi::getActiveDocumentContext()

  # check that testthat is used
  uses_testhat <- devtools::uses_testthat()
  if (!uses_testhat) {
    stop("This addin requires that testthat package is used.")
  }

  # check that testing file is opened
  uses_test_file <- grepl("tests/testthat/test-", file_info$path, fixed = T)
  if (!uses_test_file) {
    stop("This addin works only on unit test files that follow a testthat structure.")
  }

  # parse the name of the testing function
  testing_fun <- sub("^test-", "", sub(".R$", "", basename(file_info$path)))

  # get file content
  file_content <- file_info$contents

  # get locations of tests - match 'test_that("' strings
  test_that_loc <- grep('^test_that\\("', file_content)

  # define test title: add function name according to the filename and numbering
  test_title <- paste0(
    testing_fun, ", ",
    "test ", seq_along(test_that_loc), ": "
  )

  # get test descriptions
  test_that_lines <- file_info$contents[test_that_loc]
  test_that_desc_parsed <- stringr::str_extract(
    string = test_that_lines,
    pattern = paste0(
      '(?<=test_that\\(")', # positive look-ahead - search matching expression after test_that("
      '.*',    # matching expression - match everything
      '(?=")'  # positive look-behind - search matching expression before "
    )
  )

  # formulate new test descriptions (update only those that don't include test_title)
  new_descriptions <- ifelse(
    stringr::str_detect(test_that_desc_parsed, test_title),
    test_that_desc_parsed,
    paste0(test_title, test_that_desc_parsed)
  )

  # insert new test descriptions into test_that lines
  test_that_lines_updated <- stringr::str_replace(
    string = test_that_lines,
    pattern = '(?<=test_that\\(").*"',
    replacement = paste0(new_descriptions, '"')
  )

  # formulate headers according to RStudio editor functionality
  headers <- paste0("# ---- ", new_descriptions, " ----")

  # modify the file content - replace the new test descriptions
  rstudioapi::insertText(
    location = Map(c, Map(c, test_that_loc, 1), Map(c, test_that_loc, Inf)),
    text = test_that_lines_updated,
    id = file_info$id
  )

  # modify the file content - add headers if not present
  idx_valid <- !headers %in% file_content
  if (any(idx_valid)) {
    rstudioapi::insertText(
      location = Map(c, test_that_loc, 1)[idx_valid],
      text = paste0(headers[idx_valid], "\n"),
      id = file_info$id
    )
  }

  # save document
  rstudioapi::documentSave(id = file_info$id)

  return(invisible(NULL))
}
