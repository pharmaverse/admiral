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
  uses_test_file <- grepl("tests/testthat/test-", path, fixed = TRUE)
  if (!uses_test_file) {
    stop("This Addin works only on unit test files that follow a testthat structure.")
  }

  # parse the name of the testing function
  testing_file <- sub("^test-", "", sub(".R$", "", basename(path)))

  # get file content
  file_content <- readLines(path)

  # get locations of tests - match 'test_that("' strings
  test_that_loc <- grep('^test_that\\("', file_content)

  if (length(test_that_loc) == 0) {
    return(invisible(NULL))
  }

  ####
  ## HANDLE test_that DESCRIPTIONS
  ####

  # get and parse test descriptions
  test_that_lines <- file_content[test_that_loc]
  test_that_desc_parsed <- stringr::str_extract(
    string = test_that_lines,
    pattern = paste0(
      '(?<=test_that\\(")', # positive look-ahead - search matching expression after test_that("
      ".*", # matching expression - match everything
      '(?=")' # positive look-behind - search matching expression before "
    )
  )
  test_that_desc_cleaned <- stringr::str_remove(
    string = test_that_desc_parsed,
    pattern = paste0("([\\w\\.]+,? )?[Tt]est \\d{1,} ?: ")
  )

  # determine name of function which is tested
  # the function name can be specified by # function_name ---- comments
  function_name <- str_match(file_content, "# ([\\w\\.]+) ----")[, 2]
  if (is.na(function_name[1])) {
    function_name[1] <- testing_file
  }
  function_name <- tidyr::fill(data.frame(name = function_name), name)$name
  function_name <- function_name[test_that_loc]

  # formulate new test descriptions (update only those that don't include test_title)
  new_desc <- paste0(
    "Test ", seq_along(test_that_loc), ": ",
    test_that_desc_cleaned
  )

  # insert new test descriptions into test_that lines
  test_that_lines_updated <- stringr::str_replace(
    string = test_that_lines,
    pattern = '(?<=test_that\\(").*"',
    replacement = paste0(function_name, " ", new_desc, '"')
  )

  # modify the file content
  file_content[test_that_loc] <- test_that_lines_updated

  ####
  ## HANDLE HEADERS
  ####

  # formulate headers according to RStudio editor functionality
  headers <- paste0("## ", new_desc, " ----")

  # get locations of headers created by this function
  header_loc_lgl <- grepl(paste0("^##?( ----)?( \\w+)?,? [tT]est \\d{1,} ?: "), file_content)

  # remove those headers
  file_content <- file_content[!header_loc_lgl]

  # add new headers just before test_that calls
  header_loc <- grep('^test_that\\("', file_content) + seq_along(headers) - 1
  file_content_new <- vector(mode = "character", length = length(file_content) + length(headers))
  file_content_new[header_loc] <- headers
  file_content_new[-header_loc] <- file_content

  list(file_content = file_content_new)
}

# Function for the RStudio Addin, see inst/rstudio/addins.dcf.
# Wrapper of prepare_test_that_file.
format_test_that_file <- function() {
  file_info <- rstudioapi::getActiveDocumentContext()
  rstudioapi::documentSave(id = file_info$id)
  result <- prepare_test_that_file(path = file_info$path)
  rstudioapi::setDocumentContents(paste0(result$file_content, collapse = "\n"), id = file_info$id)
  rstudioapi::documentSave(id = file_info$id)
}
