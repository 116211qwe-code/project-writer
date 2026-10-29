safe_covariance_audit <- function(data, columns, log_message) {
  out <- list(columns = columns, condition_number = NA_real_, positive_definite = NA, missing = NA_integer_, singular = NA, singular_objects = character(), pass = FALSE, error = NA_character_)
  result <- tryCatch({
    x <- as.matrix(data[, columns, drop = FALSE])
    out$missing <- sum(is.na(x))
    if (out$missing > 0L) stop("协方差输入含缺失值")
    cov_mat <- stats::cov(x)
    ev <- eigen(cov_mat, symmetric = TRUE, only.values = TRUE)$values
    decomposition <- qr(cov_mat)
    if (decomposition$rank < ncol(cov_mat)) {
      out$singular_objects <- columns[decomposition$pivot[seq.int(decomposition$rank + 1L, ncol(cov_mat))]]
    }
    out$condition_number <- if (min(abs(ev)) == 0) Inf else max(abs(ev)) / min(abs(ev))
    out$positive_definite <- all(ev > 0)
    out$singular <- !is.finite(out$condition_number) || any(abs(ev) < .Machine$double.eps)
    out$pass <- isTRUE(out$positive_definite) && !isTRUE(out$singular)
    out
  }, error = function(e) {
    out$error <- conditionMessage(e)
    out
  })
  log_message(sprintf("MATRIX columns=%s condition_number=%s positive_definite=%s missing=%s singular=%s singular_objects=%s error=%s", paste(columns, collapse = ","), result$condition_number, result$positive_definite, result$missing, result$singular, paste(result$singular_objects, collapse = ","), ifelse(is.na(result$error), "", result$error)))
  result
}

group_statistics <- function(data, group_columns, value_columns) {
  key <- interaction(data[group_columns], drop = TRUE, lex.order = TRUE, sep = "|")
  groups <- split(seq_len(nrow(data)), key)
  rows <- lapply(names(groups), function(group_name) {
    idx <- groups[[group_name]]
    parts <- strsplit(group_name, "|", fixed = TRUE)[[1]]
    out <- as.list(parts)
    names(out) <- group_columns
    for (nm in value_columns) {
      x <- data[[nm]][idx]
      out[[paste0(nm, "_n")]] <- sum(!is.na(x))
      out[[paste0(nm, "_missing")]] <- sum(is.na(x))
      out[[paste0(nm, "_mean")]] <- if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)
      out[[paste0(nm, "_variance")]] <- if (sum(!is.na(x)) < 2L) NA_real_ else stats::var(x, na.rm = TRUE)
    }
    as.data.frame(out, stringsAsFactors = FALSE)
  })
  do.call(rbind, rows)
}

monthly_grazing_order_signatures <- function(data, value_columns) {
  rows <- list()
  k <- 0L
  for (nm in value_columns) {
    means <- aggregate(data[[nm]], list(month = data$month, grazing = data$grazing), mean, na.rm = TRUE)
    names(means)[3] <- "mean"
    for (m in unique(means$month)) {
      sub <- means[means$month == m, , drop = FALSE]
      sub <- sub[order(-sub$mean, sub$grazing), , drop = FALSE]
      k <- k + 1L
      rows[[k]] <- data.frame(month = m, indicator = nm, grazing_order_signature = paste(sub$grazing, collapse = ">"), means = paste(sprintf("%s=%.8g", sub$grazing, sub$mean), collapse = ";"), stringsAsFactors = FALSE)
    }
  }
  do.call(rbind, rows)
}

key_audit <- function(sheet1, soil) {
  sheet1_key <- c("grazing", "plot", "rep", "month")
  soil_key <- c("grazing", "plot", "rep", "month")
  soil_months <- unique(soil$month)
  reference <- sheet1[sheet1$month %in% soil_months, sheet1_key, drop = FALSE]
  s1_dup <- sum(duplicated(reference))
  soil_dup <- sum(duplicated(soil[soil_key]))
  key_string <- function(data) do.call(paste, c(lapply(data[sheet1_key], as.character), sep = "\r"))
  sheet1_keys <- unique(key_string(reference))
  soil_keys <- unique(key_string(soil))
  unmatched <- length(setdiff(sheet1_keys, soil_keys)) + length(setdiff(soil_keys, sheet1_keys))
  list(compared_months = paste(sort(soil_months), collapse = ","), sheet1_duplicate_keys = s1_dup, soil_duplicate_keys = soil_dup, unmatched_join_keys = unmatched, pass = s1_dup == 0L && soil_dup == 0L && unmatched == 0L)
}

run_audit <- function(workbook, log_message) {
  expected <- expected_columns()
  column_checks <- validate_workbook_structure(workbook)
  actual_names <- list(Sheet1 = names(workbook$sheet1), 土壤 = names(workbook$soil))
  missing_checks <- do.call(rbind, lapply(names(actual_names), function(sn) {
    d <- if (sn == "Sheet1") workbook$sheet1 else workbook$soil
    data.frame(sheet = sn, column = names(d), missing = vapply(d, function(x) sum(is.na(x)), integer(1)), stringsAsFactors = FALSE)
  }))
  strata <- group_statistics(workbook$sheet1, c("grazing", "month", "rep"), intersect(c("SR", "W", "T", "AGB", "BGB", "coverage", "LB", "H", "E", "D", "R"), names(workbook$sheet1)))
  soil_strata <- group_statistics(workbook$soil, c("grazing", "month", "rep"), intersect(c("BD", "pH", "C", "SOC", "P", "NH4_N", "NO3_N", "AK", "AP", "N"), names(workbook$soil)))
  order_signatures <- monthly_grazing_order_signatures(workbook$sheet1, intersect(c("SR", "W", "T"), names(workbook$sheet1)))
  keys <- key_audit(workbook$sheet1, workbook$soil)
  numeric_sheet1 <- intersect(c("W", "T", "SR", "AGB", "BGB", "coverage", "LB", "H", "E", "D", "R"), names(workbook$sheet1))
  matrix_audit <- safe_covariance_audit(workbook$sheet1, numeric_sheet1, log_message)
  numeric_soil <- intersect(c("BD", "pH", "C", "SOC", "P", "NH4_N", "NO3_N", "AK", "AP", "N"), names(workbook$soil))
  soil_matrix_audit <- safe_covariance_audit(workbook$soil, numeric_soil, log_message)
  list(
    column_checks = column_checks,
    missing_checks = missing_checks,
    strata = strata,
    soil_strata = soil_strata,
    keys = keys,
    order_signatures = order_signatures,
    matrix_audit = matrix_audit,
    soil_matrix_audit = soil_matrix_audit,
    semantic_mapping = data.frame(
      block = c("Grazing", "PlantDiv", "Biomass", "SoilEnv", "SoilNutrient", "SR"),
      indicators = c("grazing", "R,H,E,D1 (D1=1-D; derived only)", "AGB,BGB,LB", "BD,W,T,pH", "C,SOC,AP,N,NO3_N,NH4_N,AK,P", "SR"),
      observed_location = c("Sheet1", "Sheet1 (D1 absent)", "Sheet1", "土壤+Sheet1", "土壤", "Sheet1"),
      stringsAsFactors = FALSE
    ),
    pass = all(column_checks$pass) && all(missing_checks$missing == 0L) && isTRUE(keys$pass) && isTRUE(matrix_audit$pass) && isTRUE(soil_matrix_audit$pass)
  )
}
