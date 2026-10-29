read_input_workbook <- function(project_root, log_message) {
  input_path <- file.path(project_root, "data", "raw", "每月3日动态 - 副本.xlsx")
  stopifnot(requireNamespace("readxl", quietly = TRUE))
  sheets <- readxl::excel_sheets(input_path)
  if (!all(c("Sheet1", "土壤") %in% sheets)) {
    stop("输入工作簿缺少 Sheet1 或 土壤 工作表")
  }
  sheet1 <- readxl::read_excel(input_path, sheet = "Sheet1", .name_repair = "check_unique")
  soil <- readxl::read_excel(input_path, sheet = "土壤", .name_repair = "check_unique")
  log_message(sprintf("READ sheet=Sheet1 rows=%d cols=%d", nrow(sheet1), ncol(sheet1)))
  log_message(sprintf("READ sheet=土壤 rows=%d cols=%d", nrow(soil), ncol(soil)))
  list(path = input_path, sheets = sheets, sheet1 = sheet1, soil = soil)
}

validate_workbook_structure <- function(workbook) {
  expected <- expected_columns()
  actual <- list(Sheet1 = names(workbook$sheet1), 土壤 = names(workbook$soil))
  data.frame(
    sheet = names(expected),
    expected = vapply(expected, paste, "", collapse = "|"),
    actual = vapply(actual, paste, "", collapse = "|"),
    pass = mapply(identical, expected, actual),
    stringsAsFactors = FALSE
  )
}

expected_columns <- function() {
  list(
    Sheet1 = c("grazing", "month", "rep", "plot", "W", "T", "SR", "日期", "AGB", "BGB", "coverage", "LB", "H", "E", "D", "R"),
    土壤 = c("grazing", "plot", "rep", "month", "BD", "pH", "C", "SOC", "P", "NH4_N", "NO3_N", "AK", "AP", "N")
  )
}
