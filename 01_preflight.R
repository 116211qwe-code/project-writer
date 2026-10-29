preflight <- function(project_root, log_message) {
  required <- c("readxl")
  installed <- vapply(required, requireNamespace, logical(1), quietly = TRUE)
  input_path <- file.path(project_root, "data", "raw", "每月3日动态 - 副本.xlsx")
  requirements_path <- file.path(project_root, "docs", "三份文档整合版.docx")
  requirements_text_path <- file.path(project_root, "docs", "requirements_from_docx.md")
  manifest_path <- file.path(project_root, "config", "path_manifest.csv")
  manifest_ok <- file.exists(manifest_path)
  manifest_rows <- if (manifest_ok) utils::read.csv(manifest_path, stringsAsFactors = FALSE, check.names = FALSE) else data.frame()
  manifest_required_ok <- manifest_ok && all(manifest_rows$required != "yes" | file.exists(manifest_rows$absolute_path))
  checks <- data.frame(
    check = c("project_root", "input_workbook", "requirements_docx", "requirements_text", "path_manifest", "path_manifest_required_paths", paste0("package:", required)),
    actual = c(project_root, input_path, requirements_path, requirements_text_path, manifest_path, if (manifest_ok) "validated" else "missing", names(installed)),
    pass = c(dir.exists(project_root), file.exists(input_path), file.exists(requirements_path), file.exists(requirements_text_path), manifest_ok, manifest_required_ok, unname(installed)),
    stringsAsFactors = FALSE
  )
  for (i in seq_len(nrow(checks))) {
    log_message(sprintf("PREFLIGHT %s pass=%s actual=%s", checks$check[i], checks$pass[i], checks$actual[i]))
  }
  checks
}
