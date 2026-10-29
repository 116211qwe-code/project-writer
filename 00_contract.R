load_acceptance_contract <- function(project_root) {
  path <- file.path(project_root, "config", "acceptance_contract.csv")
  if (!file.exists(path)) stop("acceptance contract is missing")
  tab <- utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  required <- c("requirement_id", "data_object", "statistical_method", "target",
                "boundary_inclusive", "hardness", "approved_version", "threshold",
                "comparator", "unit", "scope")
  if (!all(required %in% names(tab)) || anyDuplicated(tab$requirement_id))
    stop("acceptance contract schema or IDs are invalid")
  if (any(!nzchar(trimws(as.character(tab$requirement_id)))) ||
      any(!nzchar(trimws(as.character(tab$data_object)))) ||
      any(!nzchar(trimws(as.character(tab$statistical_method)))) ||
      any(!nzchar(trimws(as.character(tab$target)))) ||
      any(!nzchar(trimws(as.character(tab$approved_version)))) ||
      any(!nzchar(trimws(as.character(tab$unit)))) ||
      any(!nzchar(trimws(as.character(tab$scope)))) ||
      any(!tab$hardness %in% c("HARD", "SOFT")) ||
      any(!tab$comparator %in% c("eq", "ge", "le", "gt", "lt")) ||
      any(!is.finite(as.numeric(tab$threshold))))
    stop("acceptance contract contains empty, unknown, or invalid rule values")
  tab$threshold <- as.numeric(tab$threshold)
  if (!requireNamespace("digest", quietly = TRUE)) stop("package 'digest' is required")
  contract <- list(path = normalizePath(path, mustWork = TRUE), rules = tab,
       sha256 = digest::digest(file = path, algo = "sha256", serialize = FALSE),
       loaded_at_utc = format(Sys.time(), tz = "UTC", format = "%Y-%m-%dT%H:%M:%OS3Z"))
  mapping_path <- file.path(project_root, "config", "hard_rule_mapping.csv")
  if (!file.exists(mapping_path)) stop("hard rule mapping is missing")
  mapping <- utils::read.csv(mapping_path, stringsAsFactors = FALSE,
                             check.names = FALSE)
  validate_hard_rule_mapping(contract, mapping)
  contract
}

contract_rule <- function(contract, id) {
  hit <- contract$rules[contract$rules$requirement_id == id, , drop = FALSE]
  if (nrow(hit) != 1L) stop("contract rule missing or duplicated: ", id)
  hit[1L, , drop = FALSE]
}

validate_hard_rule_mapping <- function(contract, mapping) {
  required <- c("requirement_id", "data_object", "actual_expression", "comparator",
                "threshold", "lower_bound", "upper_bound", "hardness",
                "implementation_function", "corresponding_test")
  if (!is.data.frame(mapping) || !all(required %in% names(mapping)))
    stop("hard rule mapping is missing required columns")
  if (anyDuplicated(mapping$requirement_id)) stop("hard rule mapping contains duplicate IDs")
  hard_ids <- sort(contract$rules$requirement_id[contract$rules$hardness == "HARD"])
  mapping_ids <- sort(mapping$requirement_id)
  if (!identical(hard_ids, mapping_ids)) {
    stop("hard rule mapping IDs differ from contract HARD IDs; missing=",
         paste(setdiff(hard_ids, mapping_ids), collapse = ","), "; extra=",
         paste(setdiff(mapping_ids, hard_ids), collapse = ","))
  }
  if (any(mapping$hardness != "HARD") || any(!nzchar(trimws(mapping$actual_expression))) ||
      any(!nzchar(trimws(mapping$implementation_function))) ||
      any(!nzchar(trimws(mapping$corresponding_test))))
    stop("hard rule mapping contains incomplete entries")
  rule_rows <- match(mapping$requirement_id, contract$rules$requirement_id)
  if (anyNA(rule_rows) || any(mapping$comparator != contract$rules$comparator[rule_rows]) ||
      any(as.numeric(mapping$threshold) != contract$rules$threshold[rule_rows]))
    stop("hard rule mapping comparator/threshold differs from contract")
  invisible(TRUE)
}

contract_pass <- function(contract, id, actual) {
  rule <- contract_rule(contract, id)
  if (!is.numeric(actual) || length(actual) != 1L || !is.finite(actual))
    stop("actual value for contract rule must be one finite numeric scalar: ", id)
  cmp <- as.character(rule$comparator)
  target <- rule$threshold
  if (!is.finite(target)) stop("contract threshold is not finite: ", id)
  switch(cmp, eq = actual == target, ge = actual >= target,
         le = actual <= target, gt = actual > target, lt = actual < target,
         FALSE)
}

contract_gate <- function(id, actual, contract = getOption("acceptance_contract", NULL)) {
  if (is.null(contract)) stop("acceptance contract is required; gate blocked: ", id)
  contract_pass(contract, id, actual)
}

contract_boolean <- function(contract, id, condition) {
  if (length(condition) != 1L || is.na(condition)) return(FALSE)
  contract_pass(contract, id, as.numeric(isTRUE(condition)))
}

write_acceptance_contract_evidence <- function(contract, report_dir, code_version = NA_character_) {
  out <- contract$rules
  out$config_path <- contract$path
  out$config_sha256 <- contract$sha256
  out$code_version <- code_version
  out$loaded_at_utc <- contract$loaded_at_utc
  utils::write.csv(out, file.path(report_dir, "acceptance_contract_used.csv"),
                   row.names = FALSE, fileEncoding = "UTF-8")
  invisible(out)
}
