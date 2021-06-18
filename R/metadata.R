#' Read a DAP M3
#'
#' Read all sheets of a DAP M3 into R
#'
#' @param file `character`. Path to a DAP M3 excel file
#'
#' @author Thomas Neitmann
#'
#' @return
#' A named `list` where each elements corresponds to one sheet in the original
#' excel workbook
#'
#' @export
#'
#' @examples
#' \dontrun{
#' read_dap_m3("DAP_M3.xlsx")
#' }
read_dap_m3 <- function(file) {
  assert_that(rlang::is_scalar_character(file), is_xlsx_file(file))

  is_entimice_path <- grepl("^root/", file)
  if (is_entimice_path) {
    assert_that(is_installed("rice"))
    bytes <- rice::rice_read(path, raw = TRUE)
    path <- tempfile(fileext = ".xlsx")
    # Using `c()` here is mandatory to strip off attributes from the raw
    # vector. Not doing this will lead to an error when calling `writeBin()`.
    writeBin(c(bytes), path)
  } else {
    assert_that(file.exists(file))
    path <- file
  }

  all_sheets <- readxl::excel_sheets(path)
  sheets_to_read <- setdiff(all_sheets, c("Title", "Template Instructions"))
  read <- function(sheet) {
    suppressMessages(readxl::read_excel(path, sheet, progress = FALSE))
  }
  structure(
    map(sheets_to_read, read),
    names = gsub(" ", "_", sheets_to_read),
    class = "DAP_M3"
  )
}
