
# Returns the call for updating a given test_that test file
# by adding a function name, a test number, and a section.
# Call the function either by using RStudio Addin "format_test_that_file" or
# programmatically in a for loop on the test files and running
# rstudioapi::navigateToFile and format_test_that_file
prepare_test_that_file <- function(path) {

  assert_character_scalar(path)

  # check that file exists
  if (!file.exists(path)) {
    stop("Invalid file path, the file does not exist.")
  }

  # check that testthat is used and testing file is opened
  uses_test_file <- grepl("tests/testthat/test-", path, fixed = T)
  if (!uses_test_file) {
    stop("This Addin works only on unit test files that follow a testthat structure.")
  }

  # parse the name of the testing function
  testing_fun <- sub("^test-", "", sub(".R$", "", basename(path)))

  # get file content
  file_content <- readLines(path)

  # get locations of tests - match 'test_that("' strings
  test_that_loc <- grep('^test_that\\("', file_content)

  # define test title: add function name according to the filename and numbering
  test_title <- paste0(
    testing_fun, ", ",
    "test ", seq_along(test_that_loc), ": "
  )

  # get test descriptions
  test_that_lines <- file_content[test_that_loc]
  test_that_desc_parsed <- stringr::str_extract(
    string = test_that_lines,
    pattern = paste0(
      '(?<=test_that\\(")', # positive look-ahead - search matching expression after test_that("
      ".*",    # matching expression - match everything
      '(?=")'  # positive look-behind - search matching expression before "
    )
  )

  # formulate new test descriptions (update only those that don't include test_title)
  new_desc <- ifelse(
    stringr::str_detect(test_that_desc_parsed, test_title),
    test_that_desc_parsed,
    paste0(test_title, test_that_desc_parsed)
  )

  # insert new test descriptions into test_that lines
  test_that_lines_updated <- stringr::str_replace(
    string = test_that_lines,
    pattern = '(?<=test_that\\(").*"',
    replacement = paste0(new_desc, '"')
  )

  # formulate headers according to RStudio editor functionality
  headers <- paste0("# ---- ", new_desc, " ----")

  # arguments to modify the file content - replace the new test descriptions
  l_desc <- list(
    location = Map(c, Map(c, test_that_loc, 1), Map(c, test_that_loc, Inf)),
    text = test_that_lines_updated
  )

  # arguments to modify the file content - add headers if not present
  idx_valid <- !headers %in% file_content
  l_header <- list(
    location = Map(c, test_that_loc, 1)[idx_valid],
    text = paste0(headers[idx_valid], "\n")
  )

  list(
    descriptions = l_desc,
    headers_cond = any(idx_valid),
    headers = l_header
  )
}

# Function for the RStudio Addin, see inst/rstudio/addins.dcf.
# Wrapper of prepare_test_that_file.
format_test_that_file <- function() {
  file_info <- rstudioapi::getActiveDocumentContext()
  result <- prepare_test_that_file(path = file_info$path)
  do.call(rstudioapi::insertText, append(result$descriptions, list(id = file_info$id)))
  if (result$headers_cond)
    do.call(rstudioapi::insertText, append(result$headers, list(id = file_info$id)))
  rstudioapi::documentSave(id = file_info$id)
}
