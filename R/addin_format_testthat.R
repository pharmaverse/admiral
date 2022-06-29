
# Returns the call for updating a given test_that test file
# by adding a function name, a test number, and a section.
# Call the function either by using wrapper "format_test_that_file" or
# programmatically in a for loop on list.files("tests/testthat", pattern = "^test-")
prepare_test_that_file <- function(file_info) {

  assert_s3_class(file_info, "document_context")

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
      ".*",    # matching expression - match everything
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

  # arguments to modify the file content - replace the new test descriptions
  l_desc <- list(
    location = Map(c, Map(c, test_that_loc, 1), Map(c, test_that_loc, Inf)),
    text = test_that_lines_updated,
    id = file_info$id
  )

  # arguments to modify the file content - add headers if not present
  idx_valid <- !headers %in% file_content
  l_header <- list(
    location = Map(c, test_that_loc, 1)[idx_valid],
    text = paste0(headers[idx_valid], "\n"),
    id = file_info$id
  )

  return(
    list(
      descriptions = l_desc,
      headers_cond = any(idx_valid),
      headers = l_header,
      file_info = file_info
    )
  )
}

# Function for the RStudio Addin, see inst/rstudio/addins.dcf.
# Wrapper of prepare_test_that_file.
format_test_that_file <- function() {
  result <- prepare_test_that_file(file_info = rstudioapi::getActiveDocumentContext())
  do.call(rstudioapi::insertText, result$descriptions)
  `if`(result$headers_cond, do.call(rstudioapi::insertText, result$headers))
  rstudioapi::documentSave(id = result$file_info$id)
}
