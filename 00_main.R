args <- commandArgs(trailingOnly = TRUE)
root_arg <- grep("^--project-root=", args, value = TRUE)
file_arg <- grep("^--file=", commandArgs(), value = TRUE)
file_candidate <- if (length(file_arg)) sub("^--file=", "", file_arg[[1L]]) else ""
script_path <- if (nzchar(file_candidate) && file.exists(file_candidate)) file_candidate else ""
project_root <- if (length(root_arg)) sub("^--project-root=", "", root_arg[[1L]]) else
  if (nzchar(script_path)) file.path(dirname(normalizePath(script_path, mustWork = TRUE)), "..") else getwd()
project_root <- normalizePath(project_root, mustWork = TRUE)

pipeline_stage_statuses <- c("EXECUTED_PASS", "EXECUTED_FAIL", "SKIPPED_BY_GATE",
                             "BLOCKED", "NOT_IMPLEMENTED", "NOT_APPLICABLE")

validate_pipeline_stage_status <- function(stages, path = NULL) {
  required <- c("stage", "status", "executed", "reason", "evidence_path")
  if (!is.data.frame(stages) || !all(required %in% names(stages)))
    stop("stage status table is missing required columns", call. = FALSE)
  executed_states <- stages$status %in% c("EXECUTED_PASS", "EXECUTED_FAIL")
  invalid <- is.na(stages$status) | !stages$status %in% pipeline_stage_statuses |
    is.na(stages$executed) | stages$executed != executed_states
  if (any(invalid)) {
    stages$reason[invalid] <- paste0("invalid_stage_status:", stages$status[invalid])
    if (!is.null(path)) utils::write.csv(stages, path, row.names = FALSE,
                                         fileEncoding = "UTF-8")
    stop("Illegal stage status or execution flag: ",
         paste(stages$stage[invalid], collapse = ", "), call. = FALSE)
  }
  invisible(stages)
}

candidate_generation_gate <- function(feasibility) {
  if (!is.null(feasibility) && identical(feasibility$generator_status, "BLOCKED"))
    return(list(status = "BLOCKED", executed = FALSE, generator_called = FALSE))
  list(status = "READY", executed = NA, generator_called = NA)
}

record_proxy_search_gate <- function(record_stage, calibration) {
  if (!is.null(calibration) &&
      identical(get("proxy_search_gate_status", envir = parent.frame())(calibration$has_pass), "BLOCKED")) {
    record_stage("search_diagnostics", "BLOCKED", "PROXY_NOT_CALIBRATED")
    return(invisible("BLOCKED"))
  }
  invisible(NULL)
}

finalize_pls_artifacts <- function(export, report_dir, full_pls_executed,
                                   proxy_available) {
  manifest <- export$manifest
  is_path_plot <- grepl("fig_(最终PLS路径图|路径图)\\.(emf|png|pdf)$",
                        basename(manifest$path))
  if (!full_pls_executed && any(is_path_plot)) {
    suffix <- if (proxy_available) "_PROXY_ONLY" else "_完整PLS未执行"
    for (i in which(is_path_plot)) {
      old_path <- manifest$path[i]
      extension <- tools::file_ext(old_path)
      new_path <- file.path(dirname(old_path),
                            paste0(tools::file_path_sans_ext(basename(old_path)),
                                   suffix, ".", extension))
      if (file.exists(new_path)) unlink(new_path)
      if (!file.rename(old_path, new_path))
        stop("Could not label non-final PLS diagram: ", old_path, call. = FALSE)
      manifest$path[i] <- new_path
      manifest$kind[i] <- if (proxy_available) "proxy_diagnostic" else
        "unavailable_diagnostic"
    }
  }
  utils::write.csv(manifest, file.path(report_dir, "export_manifest.csv"),
                   row.names = FALSE, fileEncoding = "UTF-8")
  artifact_type <- manifest$kind
  if (full_pls_executed) artifact_type[is_path_plot] <- "final_pls"
  if (!full_pls_executed && proxy_available) artifact_type[is_path_plot] <- "proxy_diagnostic"
  artifacts <- data.frame(path = manifest$path, artifact_type = artifact_type,
                          kind = manifest$kind, exists = file.exists(manifest$path),
                          required = TRUE, stringsAsFactors = FALSE)
  utils::write.csv(artifacts, file.path(report_dir, "artifacts_manifest.csv"),
                   row.names = FALSE, fileEncoding = "UTF-8")
  report_path <- file.path(report_dir, "export_report.md")
  pls_note <- if (full_pls_executed) "PLS status: complete PLS executed." else
    if (proxy_available) "PLS status: 完整 PLS 未执行；代理图仅供诊断（PROXY_ONLY）。" else
      "PLS status: 完整 PLS 未执行；没有可用代理图。"
  writeLines(c(readLines(report_path, encoding = "UTF-8", warn = FALSE), "", pls_note),
             report_path, useBytes = TRUE)
  export$manifest <- manifest
  export$paths <- manifest$path
  export
}

run_pipeline <- function(project_root, seed = 42L, bootstrap_reps = 200L,
                         search_iterations = 100L, search_population = 20L,
                         continue_exploration = FALSE,
                         precheck_runner = precheck_feasibility,
                         run_output_root = file.path(project_root, "reports", "runs")) {
  project_root <- normalizePath(project_root, mustWork = TRUE)
  dir.create(run_output_root, recursive = TRUE, showWarnings = FALSE)
  run_id <- sprintf("run-%s-pid-%d-seed-%d",
                    format(Sys.time(), "%Y%m%dT%H%M%OS3"), Sys.getpid(), as.integer(seed))
  run_root <- file.path(run_output_root, run_id)
  report_dir <- run_root
  log_dir <- file.path(run_root, "logs")
  dir.create(report_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)
  log_lines <- character()
  log_message <- function(message) {
    line <- sprintf("[%s] %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), message)
    log_lines <<- c(log_lines, line)
    cat(line, "\n")
  }
  stages <- data.frame(stage = character(), status = character(),
                       executed = logical(), reason = character(),
                       evidence_path = character(), generator_called = logical(),
                       stringsAsFactors = FALSE)
  stage_evidence <- function(stage) {
    file.path(run_root, paste0(gsub("[^A-Za-z0-9_-]", "_", stage), ".csv"))
  }
  record_stage <- function(stage, status, reason = "", evidence_path = NULL,
                           generator_called = NA) {
    if (is.null(evidence_path)) evidence_path <- stage_evidence(stage)
    if (stage %in% stages$stage) {
      i <- match(stage, stages$stage)
      stages$status[i] <<- status
      stages$executed[i] <<- status %in% c("EXECUTED_PASS", "EXECUTED_FAIL")
      stages$reason[i] <<- reason
      stages$evidence_path[i] <<- evidence_path
      stages$generator_called[i] <<- generator_called
    } else {
      stages <<- rbind(stages, data.frame(
        stage = stage, status = status,
        executed = status %in% c("EXECUTED_PASS", "EXECUTED_FAIL"),
        reason = reason, evidence_path = evidence_path,
        generator_called = generator_called, stringsAsFactors = FALSE))
    }
  }
  safe_stage <- function(stage, expr) {
    error <- NULL
    value <- tryCatch(force(expr), error = function(e) {
      error <<- e
      NULL
    })
    if (!is.null(error)) {
      reason <- conditionMessage(error)
      log_message(sprintf("STAGE %s ERROR %s", stage, reason))
      record_stage(stage, "EXECUTED_FAIL", reason)
      return(NULL)
    }
    passed <- !is.list(value) || is.null(value$has_pass) || isTRUE(value$has_pass)
    record_stage(stage, if (passed) "EXECUTED_PASS" else "EXECUTED_FAIL",
                 if (passed) "" else "stage_gate_not_passed")
    value
  }
  source_module <- function(name) source(file.path(project_root, "R", name), local = parent.frame())

  source_module("00_contract.R")
  contract <- load_acceptance_contract(project_root)
  hard_rule_mapping <- utils::read.csv(
    file.path(project_root, "config", "hard_rule_mapping.csv"),
    stringsAsFactors = FALSE, check.names = FALSE)
  validate_hard_rule_mapping(contract, hard_rule_mapping)
  options(acceptance_contract = contract)
  write_acceptance_contract_evidence(contract, report_dir, code_version = "working-tree")

  for (module in c("01_preflight.R", "02_read_data.R", "03_audit.R", "04_report.R",
                   "05_feasibility.R", "06_candidate_generation.R")) source_module(module)
  candidate_count_order_violations <- count_order_violations
  source_module("07_order_lock.R")
  order_count_order_violations <- count_order_violations
  # The generator (module 3) and order-lock checks (module 4) expose different
  # signatures under the same name; dispatch explicitly at the integration boundary.
  count_order_violations <- function(original, candidate, groups = NULL,
                                     columns = c("SR", "W", "T")) {
    if (is.data.frame(original) && is.data.frame(groups)) {
      order_count_order_violations(original, candidate,
                                   c("grazing", "month", "rep"), columns)
    } else {
      candidate_count_order_violations(original, candidate, groups,
                                       c("grazing", "month", "rep"), columns)
    }
  }
  check_order <- function(original_df, candidate_df,
                          groups = c("grazing", "month", "rep"),
                          columns = c("SR", "W", "T")) {
    order_count_order_violations(original_df, candidate_df, groups, columns)
  }
  for (module in c("08_pls_proxy.R", "08_search.R", "09_pls_full.R",
                   "10_statistical_acceptance.R", "09_export.R")) source_module(module)

  preflight <- safe_stage("preflight", preflight(project_root, log_message))
  preflight_pass <- !is.null(preflight) && all(preflight$pass)
  record_stage("preflight", if (preflight_pass) "EXECUTED_PASS" else "EXECUTED_FAIL",
               if (preflight_pass) "" else "preflight_not_passed")
  if (!preflight_pass) log_message("PREFLIGHT pass=FALSE; continuing for evidence and exports")

  workbook <- safe_stage("read_data_audit", {
    wb <- read_input_workbook(project_root, log_message)
    aud <- run_audit(wb, log_message)
    write_audit_outputs(project_root, preflight, wb, aud, log_lines)
    list(workbook = wb, audit = aud)
  })
  raw <- if (is.null(workbook)) NULL else workbook$workbook
  audit <- if (is.null(workbook)) NULL else workbook$audit

  feasibility <- NULL
  if (!is.null(raw)) {
    feasibility <- safe_stage("feasibility_precheck", precheck_runner(raw,
      list(conflict_core_path = file.path(report_dir, "conflict_core.csv")), log_message))
  } else record_stage("feasibility_precheck", "SKIPPED_BY_GATE", "read_data_audit_unavailable")
  if (!is.null(feasibility)) {
    feasibility_certified <- identical(feasibility$status, "PASS")
    record_stage("feasibility_precheck",
                 if (feasibility_certified) "EXECUTED_PASS" else "EXECUTED_FAIL",
                 if (feasibility_certified) "" else as.character(feasibility$status))
    for (nm in c("boundary_evidence", "distortion", "direction_evidence", "swap_pools",
                 "lock_evidence", "group_means", "residual_covariance",
                 "conflict_matrix", "triple_evidence", "sampled_ranges",
                 "constraints_evidence", "order_evidence")) {
      utils::write.csv(feasibility[[nm]], file.path(report_dir, paste0("precheck_", nm, ".csv")),
                       row.names = FALSE, fileEncoding = "UTF-8")
    }
    utils::write.csv(feasibility$boundary_evidence, file.path(report_dir, "boundary_evidence.csv"), row.names = FALSE, fileEncoding = "UTF-8")
    utils::write.csv(feasibility$constraints_evidence, file.path(report_dir, "precheck_constraint_evidence.csv"), row.names = FALSE, fileEncoding = "UTF-8")
    solver_rows <- data.frame(
      constraint = c("C1_order", "C2_bounds", "C3_statistics_pls"),
      method = c("row-wise interval arithmetic and pairwise order audit", "row-wise finite interval intersection", "no real joint solver or candidate search executed"),
      status = c(if (all(feasibility$order_evidence$feasible)) "CERTIFIED_BY_INTERVAL" else "INFEASIBLE", if (all(feasibility$boundary_evidence$lower <= feasibility$boundary_evidence$upper)) "CERTIFIED_BY_INTERVAL" else "INFEASIBLE", "NOT_ESTABLISHED"),
      samples_or_solves = c(0L, 0L, 0L),
      limitation = c("Does not establish statistical/PLS reachability", "Physical bounds must be supplied for zero values", "Requires actual candidates and complete bootstrap"), stringsAsFactors = FALSE)
    utils::write.csv(solver_rows, file.path(report_dir, "precheck_solver_evidence.csv"), row.names = FALSE)
    utils::write.csv(data.frame(alternative = c("保留原始数据并仅报告不可行证据",
                                                "仅放宽明确标注为软约束的诊断项",
                                                "重新采集或补充满足硬约束的数据")),
                     file.path(report_dir, "precheck_alternatives.csv"), row.names = FALSE,
                     fileEncoding = "UTF-8")
    log_message(sprintf("FEASIBILITY status=%s; generation_gate=%s", feasibility$status,
                        identical(feasibility$status, "PASS")))
  }

  candidate <- NULL
  feasibility_certified <- !is.null(feasibility) && identical(feasibility$status, "PASS")
  generation_allowed <- !is.null(feasibility) &&
    (feasibility_certified || isTRUE(continue_exploration))
  generation_gate <- candidate_generation_gate(feasibility)
  if (!is.null(raw) && identical(generation_gate$status, "BLOCKED")) {
    record_stage("candidate_generation", generation_gate$status, "precheck_infeasible",
                 generator_called = generation_gate$generator_called)
  } else if (!is.null(raw) && generation_allowed) {
    candidate <- safe_stage("candidate_generation", generate_parameterized_candidate(raw, feasibility, seed = seed))
    record_stage("candidate_generation", stages$status[stages$stage == "candidate_generation"],
                 if (!feasibility_certified) "user_elected_continue_exploration" else "",
                 generator_called = TRUE)
  } else record_stage("candidate_generation", "BLOCKED", "precheck_not_certified",
                      generator_called = FALSE)

  order_lock <- NULL
  if (!is.null(candidate)) {
    order_lock <- safe_stage("order_lock_check", {
      order_cols <- intersect(c("SR", "W", "T"), names(raw$sheet1))
      violations <- check_order(raw$sheet1, candidate$sheet1,
                                c("grazing", "month", "rep"), order_cols)
      tuple <- validate_tuple_lock(raw$sheet1, candidate$sheet1)
      out <- list(order_violations = violations, locked_tuples = tuple,
                  has_pass = violations == 0L && isTRUE(tuple$pass))
      utils::write.csv(data.frame(indicator = order_cols,
        group = paste(c("grazing", "month", "rep"), collapse = "*"),
        original_order = "row-level pairwise", candidate_order = "row-level pairwise",
        violations = violations), file.path(report_dir, "order_lock_evidence.csv"), row.names = FALSE)
      log_message(sprintf("ORDER_LOCK violations=%d tuple_pass=%s", violations, tuple$pass))
      out
    })
  } else record_stage("order_lock_check", "SKIPPED_BY_GATE", "candidate_unavailable")
  if (!file.exists(file.path(report_dir, "order_lock_evidence.csv")))
    utils::write.csv(data.frame(indicator = NA_character_, group = NA_character_,
      original_order = NA_character_, candidate_order = NA_character_, violations = NA_integer_,
      status = "NOT_RUN", reason = "feasibility_NOT_ESTABLISHED"),
      file.path(report_dir, "order_lock_evidence.csv"), row.names = FALSE)

  blocks <- list(Grazing = "grazing", PlantDiv = c("R", "H", "E", "D1"),
                 Biomass = c("AGB", "BGB", "LB"), SoilEnv = c("BD", "W", "T", "pH"),
                 SoilNutrient = c("C", "SOC", "AP", "N", "NO3_N", "NH4_N", "AK", "P"),
                 SR = "SR")
  modes <- c(Grazing = "A", PlantDiv = "A", Biomass = "B", SoilEnv = "B",
             SoilNutrient = "B", SR = "A")
  paths <- data.frame(from = c("Grazing", "Grazing", "Grazing", "Grazing", "Grazing",
                               "PlantDiv", "Biomass", "SoilEnv", "SoilEnv", "Biomass",
                               "SoilNutrient"),
                      to = c("PlantDiv", "Biomass", "SoilEnv", "SoilNutrient", "SR",
                              "SR", "SR", "SR", "SoilNutrient", "SoilNutrient", "SR"),
                      stringsAsFactors = FALSE)
  july_data <- function(value) {
    main <- value$sheet1[as.character(value$sheet1$month) == "7", , drop = FALSE]
    soil <- value$soil[as.character(value$soil$month) == "7", , drop = FALSE]
    keys <- c("grazing", "plot", "rep", "month")
    merged <- merge(main, soil, by = keys, all = FALSE, suffixes = c("", ".soil"), sort = FALSE)
    merged$D1 <- 1 - merged$D
    merged
  }
  analysis_data <- if (!is.null(candidate)) july_data(candidate) else NULL
  fast <- NULL
  if (!is.null(analysis_data)) {
    fast <- safe_stage("pls_fast_proxy", {
      proxy <- pls_fast_proxy(analysis_data, blocks, paths)
      gate <- evaluate_fast(proxy)
      list(proxy = proxy, gate = gate, has_pass = gate$has_pass)
    })
  } else record_stage("pls_fast_proxy", "SKIPPED_BY_GATE", "candidate_unavailable")

  calibration <- NULL
  if (generation_allowed && !is.null(candidate) && !is.null(fast)) {
    calibration <- safe_stage("pls_proxy_calibration", {
      calibration_candidates <- lapply(seq_len(100L), function(i)
        generate_parameterized_candidate(raw, feasibility, seed = seed + 10000L + i))
      proxy_fn <- function(x) pls_fast_proxy(july_data(x), blocks, paths)
      full_fn <- function(x) fit_pls_full(july_data(x), blocks, paths,
                                          modes = modes, bootstrap_reps = 100L,
                                          seed = seed + 20000L)
      out <- calibrate_pls_proxy_set(calibration_candidates, proxy_fn, full_fn,
        output_path = file.path(report_dir, "proxy_calibration.csv"))
      list(n = length(calibration_candidates), paired = out$paired,
           calibration = out$calibration,
           has_pass = isTRUE(out$calibration$calibration_pass))
    })
  } else record_stage("pls_proxy_calibration", "SKIPPED_BY_GATE",
                      if (!feasibility_certified) "precheck_not_certified" else
                        "candidate_or_proxy_unavailable")

  full <- NULL
  full_attempted <- !is.null(fast) && isTRUE(fast$has_pass) &&
    !is.null(calibration) && isTRUE(calibration$has_pass)
  if (full_attempted) {
    full <- safe_stage("pls_full_bootstrap", {
      fit <- fit_pls_full(analysis_data, blocks, paths, modes = modes,
                          bootstrap_reps = bootstrap_reps, seed = seed)
      list(result = fit, evaluation = evaluate_pls_full(fit), has_pass = fit$has_pass)
    })
  } else record_stage("pls_full_bootstrap", "SKIPPED_BY_GATE",
                      if (is.null(fast) || !isTRUE(fast$has_pass)) "proxy_not_passed" else
                        "proxy_not_calibrated")

  acceptance <- NULL
  statistical_inputs <- NULL
  if (!is.null(raw)) {
    statistics_workbook <- if (is.null(candidate)) raw else candidate
    july_statistics_data <- july_data(statistics_workbook)
    acceptance <- safe_stage("statistical_acceptance", {
      spec <- statistical_acceptance_columns()
      agg75 <- aggregate_pearson_data(statistics_workbook$sheet1, spec$pearson)
      anova_entry <- accept_anova(july_statistics_data)
      pearson_entry <- accept_pearson(agg75)
      soil_july <- statistics_workbook$soil[
        as.character(statistics_workbook$soil$month) == "7", , drop = FALSE]
      mantel_entry <- accept_mantel(statistics_workbook$sheet1, soil_july,
                                    seed = seed)
      soil_entry <- accept_soil_internal(july_statistics_data)
      assemble_statistical_acceptance(
        anova_entry, agg75, pearson_entry, mantel_entry, soil_entry,
        ecology_source = merge(statistics_workbook$sheet1,
          statistics_workbook$soil,
          by = c("grazing", "plot", "rep", "month"),
          all = FALSE, sort = FALSE),
        permutations = 100L, seed = seed)
    })
    statistical_inputs <- list(pearson = if (is.null(acceptance)) NULL else acceptance$pearson_data,
                               anova = july_statistics_data,
                               mantel = july_statistics_data,
                               soil_correlation = july_statistics_data)
    if (!is.null(acceptance)) {
      utils::write.csv(acceptance$checks, file.path(report_dir, "statistical_acceptance.csv"), row.names = FALSE, fileEncoding = "UTF-8")
      utils::write.csv(acceptance$anova, file.path(report_dir, "statistical_anova.csv"), row.names = FALSE, fileEncoding = "UTF-8")
      utils::write.csv(acceptance$lsd, file.path(report_dir, "statistical_lsd.csv"), row.names = FALSE, fileEncoding = "UTF-8")
      utils::write.csv(acceptance$pearson, file.path(report_dir, "statistical_pearson.csv"), row.names = FALSE, fileEncoding = "UTF-8")
      utils::write.csv(acceptance$mantel, file.path(report_dir, "statistical_mantel.csv"), row.names = FALSE, fileEncoding = "UTF-8")
      utils::write.csv(acceptance$sr_soil_pearson, file.path(report_dir, "statistical_sr_soil_pearson.csv"), row.names = FALSE, fileEncoding = "UTF-8")
      utils::write.csv(acceptance$ecology_full, file.path(report_dir, "ecology_signs_full.csv"), row.names = FALSE, fileEncoding = "UTF-8")
      utils::write.csv(acceptance$ecology_july, file.path(report_dir, "ecology_signs_july.csv"), row.names = FALSE, fileEncoding = "UTF-8")
    }
  } else record_stage("statistical_acceptance", "SKIPPED_BY_GATE",
                      "analysis_data_unavailable")

  search <- NULL
  search_allowed <- generation_allowed && !is.null(calibration) &&
    isTRUE(calibration$has_pass)
  if (!is.null(raw) && search_allowed && !is.null(candidate)) {
    search <- safe_stage("search_diagnostics", run_global_search(
      candidate_generator = function(n, distribution, iteration_seed, restart_id) {
        lapply(seq_len(n), function(i) generate_parameterized_candidate(
          raw, feasibility, seed = iteration_seed + i + restart_id * 1000L,
          scale = distribution$scale, mean_shift = distribution$mean_shift))
      },
      evaluator = function(x) {
        d <- july_data(x)
        p <- pls_fast_proxy(d, blocks, paths)
        e <- evaluate_fast(p)
        ord <- check_order(raw$sheet1, x$sheet1,
                           c("grazing", "month", "rep"),
                           intersect(c("SR", "W", "T"), names(x$sheet1)))
        lock <- validate_tuple_lock(raw$sheet1, x$sheet1)
        gates <- c(has_lock = isTRUE(lock$pass) && isTRUE(x$metadata$bounds_preserved),
                   has_order = ord == 0L,
                   has_ecology = p$sign_violations == 0L,
                   has_pls_skeleton = isTRUE(e$has_pass))
        residuals <- c(e$residual, lock_violation = as.numeric(!gates[["has_lock"]]),
                       order_violations = ord, ecology_sign_violations = p$sign_violations)
        score <- 0 - sum(residuals)
        list(score = max(-1e9, score), residuals = residuals, gates = gates,
             has_pass = all(gates), sign_violations = p$sign_violations,
             max_vif = p$max_vif, r2_proxy = max(p$r2, na.rm = TRUE),
             pls_proxy_diff = e$max_residual)
      },
      config = list(iterations = search_iterations, population_size = search_population,
                    elite_size = min(2L, search_population), seed = seed,
                    sign_penalty = 100, temperature_start = 1, temperature_decay = .95),
      initial_distribution = list(covariance_id = "parameterized_candidate")))
    if (!is.null(search)) utils::write.csv(search$records, file.path(report_dir, "search_records.csv"), row.names = FALSE)
  } else if (identical(record_proxy_search_gate(record_stage, calibration), "BLOCKED")) {
    invisible(NULL)
  } else record_stage("search_diagnostics", "SKIPPED_BY_GATE",
                      if (!feasibility_certified && !isTRUE(continue_exploration))
                        "precheck_not_certified" else "proxy_not_calibrated")
  if (is.null(search)) {
    utils::write.csv(data.frame(candidate_id=NA_character_, seed=NA_integer_, parent_id=NA_character_,
      restart_id=NA_integer_, hard_gate=NA_character_, residual_vector=NA_character_,
      status="SKIPPED_BY_GATE", reason="feasibility_NOT_ESTABLISHED", attempts=0L,
      stringsAsFactors=FALSE),
      file.path(report_dir, "search_records.csv"), row.names=FALSE)
    utils::write.csv(data.frame(restart_id=NA_integer_, seed=NA_integer_, mean_shift=NA_real_,
      scale=NA_real_, parameter_hash=NA_character_, status="SKIPPED_BY_GATE",
      reason="feasibility_NOT_ESTABLISHED", stringsAsFactors=FALSE),
      file.path(report_dir, "search_restart_evidence.csv"), row.names=FALSE)
  }
  if (!file.exists(file.path(report_dir, "generation_trace.csv")))
    utils::write.csv(data.frame(seed=NA_integer_, scale=NA_real_, mean_shift=NA_real_,
      stage="generation", before_target=NA_real_, after_target=NA_real_, violations=NA_integer_,
      status="SKIPPED_BY_GATE", reason="feasibility_NOT_ESTABLISHED", attempts=0L,
      stringsAsFactors=FALSE),
      file.path(report_dir, "generation_trace.csv"), row.names=FALSE)

  # The selected search candidate is the object carried into final validation
  # and export; never report statistics from a different candidate.
  if (!is.null(search$best$candidate)) {
    candidate <- search$best$candidate
    analysis_data <- july_data(candidate)
    ord <- check_order(raw$sheet1, candidate$sheet1,
                       c("grazing", "month", "rep"),
                       intersect(c("SR", "W", "T"), names(raw$sheet1)))
    tuple <- validate_tuple_lock(raw$sheet1, candidate$sheet1)
    order_lock <- list(order_violations=ord, locked_tuples=tuple,
                       has_pass=ord == 0L && isTRUE(tuple$pass) && isTRUE(candidate$metadata$bounds_preserved))
    proxy <- pls_fast_proxy(analysis_data, blocks, paths)
    gate <- evaluate_fast(proxy)
    fast <- list(proxy=proxy, gate=gate, has_pass=gate$has_pass)
    selected_full_attempted <- isTRUE(fast$has_pass) && !is.null(calibration) &&
      isTRUE(calibration$has_pass)
    full <- if (selected_full_attempted) tryCatch({
      fit <- fit_pls_full(analysis_data, blocks, paths, modes=modes,
                          bootstrap_reps=bootstrap_reps, seed=seed + 30000L)
      list(result=fit, evaluation=evaluate_pls_full(fit), has_pass=fit$has_pass)
    }, error=function(e) {
      log_message(sprintf("FINAL_SELECTED_PLS ERROR %s", conditionMessage(e)))
      NULL
    }) else NULL
    july_statistics_data <- analysis_data
    acceptance <- tryCatch({
      spec <- statistical_acceptance_columns()
      agg75 <- aggregate_pearson_data(candidate$sheet1, spec$pearson)
      anova_entry <- accept_anova(july_statistics_data)
      pearson_entry <- accept_pearson(agg75)
      soil_july <- candidate$soil[
        as.character(candidate$soil$month) == "7", , drop = FALSE]
      mantel_entry <- accept_mantel(candidate$sheet1, soil_july,
                                    seed = seed + 40000L)
      soil_entry <- accept_soil_internal(july_statistics_data)
      assemble_statistical_acceptance(
        anova_entry, agg75, pearson_entry, mantel_entry, soil_entry,
        ecology_source = merge(candidate$sheet1, candidate$soil,
          by = c("grazing", "plot", "rep", "month"),
          all = FALSE, sort = FALSE),
        permutations = 100L, seed = seed + 40000L)
    },
      error=function(e) {
        log_message(sprintf("FINAL_SELECTED_STATISTICS ERROR %s", conditionMessage(e)))
        NULL
      })
    statistical_inputs <- list(
      pearson=if (is.null(acceptance)) NULL else acceptance$pearson_data,
      anova=july_statistics_data, mantel=july_statistics_data,
      soil_correlation=july_statistics_data)
    if (!is.null(acceptance)) {
      utils::write.csv(acceptance$checks, file.path(report_dir, "statistical_acceptance.csv"), row.names=FALSE, fileEncoding="UTF-8")
      utils::write.csv(acceptance$anova, file.path(report_dir, "statistical_anova.csv"), row.names=FALSE, fileEncoding="UTF-8")
      utils::write.csv(acceptance$lsd, file.path(report_dir, "statistical_lsd.csv"), row.names=FALSE, fileEncoding="UTF-8")
      utils::write.csv(acceptance$pearson, file.path(report_dir, "statistical_pearson.csv"), row.names=FALSE, fileEncoding="UTF-8")
      utils::write.csv(acceptance$mantel, file.path(report_dir, "statistical_mantel.csv"), row.names=FALSE, fileEncoding="UTF-8")
      utils::write.csv(acceptance$sr_soil_pearson, file.path(report_dir, "statistical_sr_soil_pearson.csv"), row.names=FALSE, fileEncoding="UTF-8")
      utils::write.csv(acceptance$ecology_full, file.path(report_dir, "ecology_signs_full.csv"), row.names=FALSE, fileEncoding="UTF-8")
      utils::write.csv(acceptance$ecology_july, file.path(report_dir, "ecology_signs_july.csv"), row.names=FALSE, fileEncoding="UTF-8")
    }
    record_stage("order_lock_check",
                 if (isTRUE(order_lock$has_pass)) "EXECUTED_PASS" else "EXECUTED_FAIL",
                 if (isTRUE(order_lock$has_pass)) "" else "order_or_lock_check_failed")
    record_stage("pls_fast_proxy",
                 if (isTRUE(fast$has_pass)) "EXECUTED_PASS" else "EXECUTED_FAIL",
                 if (isTRUE(fast$has_pass)) "" else "proxy_gate_not_passed")
    if (is.null(full) && selected_full_attempted) {
      record_stage("pls_full_bootstrap", "EXECUTED_FAIL",
                   "selected_candidate_full_pls_failed")
    } else if (is.null(full)) {
      record_stage("pls_full_bootstrap", "SKIPPED_BY_GATE",
                   "selected_candidate_full_pls_not_executed")
    } else {
      record_stage("pls_full_bootstrap",
                   if (isTRUE(full$has_pass)) "EXECUTED_PASS" else "EXECUTED_FAIL",
                   if (isTRUE(full$has_pass)) "" else "full_pls_gate_not_passed")
    }
    if (is.null(acceptance)) {
      record_stage("statistical_acceptance", "EXECUTED_FAIL",
                   "selected_candidate_statistical_validation_failed")
    } else {
      record_stage("statistical_acceptance",
                   if (isTRUE(acceptance$has_pass)) "EXECUTED_PASS" else "EXECUTED_FAIL",
                   if (isTRUE(acceptance$has_pass)) "" else "statistical_acceptance_failed")
    }
  }

  pass <- isTRUE(preflight_pass) && !is.null(audit) && isTRUE(audit$pass) &&
    !is.null(feasibility) && identical(feasibility$status, "PASS") &&
    !is.null(candidate) && isTRUE(candidate$metadata$has_pass) &&
    !is.null(order_lock) && isTRUE(order_lock$has_pass) &&
    !is.null(search) && isTRUE(search$best$has_pass) &&
    !is.null(calibration) && isTRUE(calibration$has_pass) &&
    !is.null(fast) && isTRUE(fast$has_pass) && !is.null(full) && isTRUE(full$has_pass) &&
    !is.null(acceptance) && isTRUE(acceptance$has_pass)
  results <- list(run_ok = TRUE, pass = pass, stages = stages, preflight = preflight,
                  audit = audit, feasibility = feasibility, candidate = candidate,
                  order_lock = order_lock, fast = fast, calibration = calibration, full = full,
                  acceptance = acceptance, statistical_inputs = statistical_inputs,
                  search = search)
  record_stage("export", "EXECUTED_PASS", "export_in_progress")
  validate_pipeline_stage_status(stages, file.path(report_dir, "stage_status.csv"))
  export <- safe_stage("export", finalize_pls_artifacts(
    export_artifacts(
      results, run_root, if (is.null(candidate)) raw$sheet1 else candidate$sheet1,
      tables = list(stage_status = stages),
      metadata = list(RUN_OK = TRUE, PASS = pass, project_root = project_root),
      export_schema = list(sheet1 = raw$sheet1, soil = raw$soil)),
    report_dir, full_pls_executed = !is.null(full),
    proxy_available = !is.null(fast) && !is.null(fast$proxy$paths)))
  hash_files <- function(paths) {
    paths <- sort(normalizePath(paths[file.exists(paths)], mustWork = TRUE))
    if (!length(paths)) return(NA_character_)
    hashes <- vapply(paths, function(path) digest::digest(
      file = path, algo = "sha256", serialize = FALSE), character(1))
    digest::digest(hashes, algo = "sha256", serialize = TRUE)
  }
  input_files <- list.files(file.path(project_root, "data"),
    pattern = "\\.(xlsx|xls|csv)$", recursive = TRUE, full.names = TRUE,
    ignore.case = TRUE)
  code_files <- list.files(file.path(project_root, "R"), pattern = "\\.R$",
                           full.names = TRUE)
  utils::write.csv(data.frame(input="raw workbook + acceptance contract",
    input_sha256=hash_files(input_files), code_sha256=hash_files(code_files),
    config_sha256=contract$sha256, seed=seed,
    attempts=if (is.null(search)) 0L else nrow(search$records),
    full_pls_executions=as.integer(!is.null(full)),
    final_status=if (isTRUE(pass)) "PASS" else "FAIL",
    output_role=if (isTRUE(pass)) "FINAL" else "CANDIDATE_DIAGNOSTIC",
    run_id=run_id, stringsAsFactors=FALSE), file.path(report_dir, "run_manifest.csv"), row.names=FALSE)
  required_outputs <- c("stage_status.csv", "export_manifest.csv",
                        "artifacts_manifest.csv", "export_report.md", "run_manifest.csv")
  run_ok <- !is.null(export) && run_outputs_pass(
    pass, report_dir, export$manifest, required_outputs, run_id)
  output_hash <- if (!is.null(export) && is.data.frame(export$manifest) && requireNamespace("digest", quietly = TRUE))
    paste(vapply(export$manifest$path, function(p) if (file.exists(p)) digest::digest(file=p, algo="sha256", serialize=FALSE) else NA_character_, character(1)), collapse=";") else NA_character_
  run_manifest <- utils::read.csv(file.path(report_dir, "run_manifest.csv"), stringsAsFactors=FALSE)
  run_manifest$output_sha256 <- output_hash
  utils::write.csv(run_manifest, file.path(report_dir, "run_manifest.csv"), row.names=FALSE)
  record_stage("export", if (run_ok) "EXECUTED_PASS" else "EXECUTED_FAIL",
               if (run_ok) "" else "required_export_missing")
  validate_pipeline_stage_status(stages, file.path(report_dir, "stage_status.csv"))
  utils::write.csv(stages, file.path(report_dir, "stage_status.csv"),
                   row.names = FALSE, fileEncoding = "UTF-8")
  module_names <- c("路径和数据审计","可行域预检和冲突证据","参数化候选生成",
    "排序、值域和 H/E/D/R 锁定","PLS 快速代理及误差校准",
    "完整 PLS 与 Bootstrap","Pearson、ANOVA、Mantel 等验收",
    "搜索和重启","图、表、数据和日志导出")
  module_stages <- c("read_data_audit","feasibility_precheck","candidate_generation",
    "order_lock_check","pls_proxy_calibration","pls_full_bootstrap",
    "statistical_acceptance","search_diagnostics","export")
  module_status <- data.frame(module=module_names, stage=module_stages,
    status=vapply(module_stages, function(s) stages$status[match(s, stages$stage)], character(1)),
    executed=vapply(module_stages, function(s) any(stages$stage == s & stages$executed), logical(1)),
    stringsAsFactors=FALSE)
  utils::write.csv(module_status, file.path(report_dir, "module_status.csv"), row.names=FALSE, fileEncoding="UTF-8")
  # Persist reproducibility and checklist evidence from this exact run. These
  # files distinguish executed code from data-level PASS decisions.
  if (!is.null(search) && is.data.frame(search$records)) {
    records <- search$records
    checkpoints <- records[records$iteration %% 100L == 0L, , drop = FALSE]
    if (!nrow(checkpoints) && nrow(records)) checkpoints <- records[nrow(records), , drop = FALSE]
    utils::write.csv(checkpoints, file.path(report_dir, "diagnostic_checkpoints.csv"), row.names=FALSE, fileEncoding="UTF-8")
    seed_cols <- intersect(c("iteration", "candidate_id", "iteration_seed", "parent_seed", "restart_id", "parent_id", "temperature", "annealing_accepted"), names(records))
    utils::write.csv(records[seed_cols], file.path(report_dir, "seed_chain.csv"), row.names=FALSE, fileEncoding="UTF-8")
  }
  writeLines(c(
    "Rscript --vanilla tests/run_tests.R",
    "Rscript --vanilla -e \"source('R/00_main.R', local=TRUE); run_pipeline(getwd(), seed=7L, bootstrap_reps=2L, search_iterations=1L, search_population=1L)\"",
    sprintf("RUN_OK=%s", run_ok), sprintf("PASS=%s", pass),
    sprintf("timestamp=%s", format(Sys.time(), "%Y-%m-%d %H:%M:%S%z"))
  ), file.path(report_dir, "verification_commands.txt"), useBytes=TRUE)
  capture.output(sessionInfo(), file=file.path(report_dir, "sessionInfo.txt"))
  if (!is.null(acceptance) && is.data.frame(acceptance$checks)) {
    md <- c("# Statistical acceptance", "", sprintf("RUN_OK: %s; PASS: %s", run_ok, pass), "",
      "| Requirement | Target | Actual | Pass |", "|---|---|---|---|")
    md <- c(md, vapply(seq_len(nrow(acceptance$checks)), function(i)
      sprintf("| %s | %s | %s | %s |", acceptance$checks$requirement[i],
              acceptance$checks$target[i], acceptance$checks$actual[i], acceptance$checks$pass[i]), character(1)))
    writeLines(md, file.path(report_dir, "statistical_acceptance.md"), useBytes=TRUE)
  }
  results$module_status <- module_status
  results$run_ok <- run_ok
  results$pass <- pass
  results$stages <- stages
  results$run_id <- run_id
  results$run_dir <- normalizePath(run_root, mustWork = TRUE)
  saveRDS(results, file.path(report_dir, "pipeline_state.rds"))
  report_lines <- c("# Integrated pipeline report", "", paste0("RUN_OK: ", run_ok),
                    paste0("PASS: ", pass), "", "## Stage status", capture.output(print(stages)),
                    "", if (is.null(full)) "完整 PLS 未执行；PLS 图仅为代理诊断或不可用占位图。" else "完整 PLS 已执行。",
                    "", "## Evidence", "Precheck evidence and alternatives are retained even when feasibility is not certified.")
  writeLines(report_lines, file.path(report_dir, "pipeline_report.md"), useBytes = TRUE)
  writeLines(c(log_lines, sprintf("RUN_OK=%s PASS=%s", run_ok, pass)), file.path(log_dir, "run.log"), useBytes = TRUE)
  cat(sprintf("RUN_OK=%s PASS=%s\n", run_ok, pass))
  invisible(results)
}

if (sys.nframe() == 0L) run_pipeline(project_root)
