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
  assert_that(is.character(file), file.exists(file))

  all_sheets <- readxl::excel_sheets(file)
  sheets_to_read <- setdiff(all_sheets, c("Title", "Template Instructions"))
  read <- function(sheet) {
    suppressMessages(readxl::read_excel(file, sheet, progress = FALSE))
  }
  structure(
    map(sheets_to_read, read),
    names = gsub(" ", "_", sheets_to_read),
    class = "DAP_M3"
  )
}
