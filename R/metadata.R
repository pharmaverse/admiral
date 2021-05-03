read_dap_m3 <- function(file) {
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
