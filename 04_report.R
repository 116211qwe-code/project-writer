write_audit_outputs <- function(project_root, preflight, workbook, audit, log_lines) {
  report_dir <- file.path(project_root, "reports")
  log_dir <- file.path(project_root, "logs")
  dir.create(report_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(preflight, file.path(report_dir, "preflight.csv"), row.names = FALSE, fileEncoding = "UTF-8")
  utils::write.csv(audit$column_checks, file.path(report_dir, "audit_columns.csv"), row.names = FALSE, fileEncoding = "UTF-8")
  utils::write.csv(audit$missing_checks, file.path(report_dir, "audit_missing.csv"), row.names = FALSE, fileEncoding = "UTF-8")
  utils::write.csv(audit$strata, file.path(report_dir, "audit_strata.csv"), row.names = FALSE, fileEncoding = "UTF-8")
  utils::write.csv(audit$soil_strata, file.path(report_dir, "audit_soil_strata.csv"), row.names = FALSE, fileEncoding = "UTF-8")
  utils::write.csv(audit$order_signatures, file.path(report_dir, "audit_order_signatures.csv"), row.names = FALSE, fileEncoding = "UTF-8")
  utils::write.csv(data.frame(metric = names(audit$keys), value = unlist(audit$keys, use.names = FALSE)), file.path(report_dir, "audit_keys.csv"), row.names = FALSE, fileEncoding = "UTF-8")
  utils::write.csv(data.frame(
    sheet = c("Sheet1", "土壤"),
    condition_number = c(audit$matrix_audit$condition_number, audit$soil_matrix_audit$condition_number),
    positive_definite = c(audit$matrix_audit$positive_definite, audit$soil_matrix_audit$positive_definite),
    missing = c(audit$matrix_audit$missing, audit$soil_matrix_audit$missing),
    singular = c(audit$matrix_audit$singular, audit$soil_matrix_audit$singular),
    singular_objects = c(paste(audit$matrix_audit$singular_objects, collapse = ","), paste(audit$soil_matrix_audit$singular_objects, collapse = ",")),
    error = c(audit$matrix_audit$error, audit$soil_matrix_audit$error),
    pass = c(audit$matrix_audit$pass, audit$soil_matrix_audit$pass),
    stringsAsFactors = FALSE
  ), file.path(report_dir, "audit_matrices.csv"), row.names = FALSE, fileEncoding = "UTF-8")
  utils::write.csv(audit$semantic_mapping, file.path(report_dir, "semantic_mapping.csv"), row.names = FALSE, fileEncoding = "UTF-8")
  saveRDS(list(preflight = preflight, workbook = workbook[c("path", "sheets")], audit = audit), file.path(report_dir, "audit_state.rds"))
  report <- c(
    "# Clean-room R 工程骨架审计报告",
    "",
    "本报告以工程内 三份文档整合版.docx（完整需求）及其转换文本 requirements_from_docx.md、data/raw/每月3日动态 - 副本.xlsx 为依据。旧 requirements.md 已标记为非权威补充文件，不与完整需求静默合并。各模块的代码实现和证据状态由主流程继续记录；统计条件未通过时不伪造 PASS。",
    "",
    sprintf("- RUN_OK: %s", all(preflight$pass)),
    sprintf("- 输入: %s", workbook$path),
    sprintf("- 工作表: %s", paste(workbook$sheets, collapse = ", ")),
    sprintf("- 审计 PASS: %s", audit$pass),
    sprintf("- 土壤表重复连接键数: %d", audit$keys$soil_duplicate_keys),
    sprintf("- 未匹配连接键数: %d", audit$keys$unmatched_join_keys),
    sprintf("- 矩阵条件数: %s", audit$matrix_audit$condition_number),
    sprintf("- 矩阵正定性: %s", audit$matrix_audit$positive_definite),
    sprintf("- 土壤矩阵条件数: %s", audit$soil_matrix_audit$condition_number),
    sprintf("- 土壤矩阵正定性: %s", audit$soil_matrix_audit$positive_definite),
    "",
    "## 运行状态",
    "",
    "- 依赖检查、路径检查、工作簿读取和基础审计：已实现。",
    "- 数据构造、预检、PLS、搜索、统计验收、图件和最终 CSV/RDS：由集成主流程执行并保留证据；未通过条件仍明确标记 PASS=FALSE。",
    "- PASS=FALSE 不会阻断报告写出。",
    "",
    "## 关键结构冲突",
    "",
    "- 完整需求文档对数据对象的描述需与输入工作簿逐项核对；实际土壤指标位于独立的 土壤 工作表，而非 Sheet1。",
    "- PlantDiv 要求 D1=1-D；原始表没有 D1，骨架只登记派生规则，不修改原始数据。",
    "",
    "## 输出文件",
    "",
    "- preflight.csv",
    "- audit_columns.csv",
    "- audit_missing.csv",
    "- audit_strata.csv",
    "- audit_soil_strata.csv",
    "- audit_order_signatures.csv",
    "- audit_keys.csv",
    "- audit_matrices.csv",
    "- audit_state.rds",
    "- semantic_mapping.csv",
    "- audit_state.rds",
    ""
  )
  writeLines(report, file.path(report_dir, "audit_report.md"), useBytes = TRUE)
  writeLines(log_lines, file.path(log_dir, "run.log"), useBytes = TRUE)
  invisible(report_dir)
}
