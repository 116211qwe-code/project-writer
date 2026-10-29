# Module 6: complete PLS-PM fit and auditable bootstrap diagnostics.
# Teaching simulation only; not for real scientific inference.

pls_default_signs <- c(
  "Grazing->PlantDiv" = 0, "Grazing->Biomass" = -1,
  "Grazing->SoilEnv" = 0, "Grazing->SoilNutrient" = -1,
  "Grazing->SR" = -1, "PlantDiv->SR" = 0, "Biomass->SR" = 1,
  "SoilNutrient->SR" = 1,
  "SoilEnv->SR" = 0, "SoilEnv->SoilNutrient" = 0,
  "Biomass->SoilNutrient" = 1
)

check_r042_path_signs <- function(path_coefs) {
  required <- c("Grazing", "Biomass", "SR", "SoilNutrient")
  if (!is.matrix(path_coefs) || !all(required %in% rownames(path_coefs)) ||
      !all(required %in% colnames(path_coefs)))
    stop("required structural path coefficients are missing", call. = FALSE)
  grazing_biomass <- path_coefs["Biomass", "Grazing"]
  nutrient_sr <- path_coefs["SR", "SoilNutrient"]
  if (is.finite(grazing_biomass) && grazing_biomass > 0)
    stop("生态逻辑冲突：放牧未抑制生物量", call. = FALSE)
  if (is.finite(nutrient_sr) && nutrient_sr < 0)
    stop("生态逻辑冲突：土壤养分对 SR 为负", call. = FALSE)
  invisible(TRUE)
}

validate_pls_inputs <- function(data, blocks, paths, modes = NULL) {
  if (!requireNamespace("plspm", quietly = TRUE))
    stop("package 'plspm' is required for complete PLS-PM", call. = FALSE)
  if (!is.data.frame(data) || nrow(data) < 10L)
    stop("data must be a data.frame with at least 10 observations", call. = FALSE)
  if (!is.list(blocks) || !length(blocks) || is.null(names(blocks)) || any(!nzchar(names(blocks))))
    stop("blocks must be a named non-empty list", call. = FALSE)
  indicators <- unique(unlist(blocks, use.names = FALSE))
  missing <- setdiff(indicators, names(data))
  if (length(missing)) stop("missing block indicators: ", paste(missing, collapse = ", "), call. = FALSE)
  if (!all(vapply(data[indicators], is.numeric, logical(1))) ||
      any(!is.finite(as.matrix(data[indicators]))))
    stop("block indicators must be finite numeric columns", call. = FALSE)
  if (anyDuplicated(indicators)) stop("an indicator cannot occur in multiple blocks", call. = FALSE)
  if (!is.data.frame(paths) || !all(c("from", "to") %in% names(paths)))
    stop("paths must contain from and to columns", call. = FALSE)
  if (!nrow(paths)) stop("at least one structural path is required", call. = FALSE)
  if (any(!as.character(paths$from) %in% names(blocks)) || any(!as.character(paths$to) %in% names(blocks)))
    stop("paths refer to unknown latent blocks", call. = FALSE)
  if (any(as.character(paths$from) == as.character(paths$to))) stop("self paths are not allowed", call. = FALSE)
  if (is.null(modes)) modes <- setNames(rep("A", length(blocks)), names(blocks))
  if (is.null(names(modes)) || !all(names(blocks) %in% names(modes)))
    stop("modes must be named for every block", call. = FALSE)
  modes <- modes[names(blocks)]
  if (any(!modes %in% c("A", "B"))) stop("modes must contain only A or B", call. = FALSE)
  list(indicators = indicators, modes = modes)
}

make_path_matrix <- function(block_names, paths) {
  pm <- matrix(0, length(block_names), length(block_names),
               dimnames = list(block_names, block_names))
  for (i in seq_len(nrow(paths))) {
    from <- as.character(paths$from[i]); to <- as.character(paths$to[i])
    if (pm[to, from] == 1) stop("duplicate structural path: ", from, " -> ", to, call. = FALSE)
    pm[to, from] <- 1
  }
  pm
}

fit_pls_once <- function(data, blocks, paths, modes) {
  pm <- make_path_matrix(names(blocks), paths)
  fit <- tryCatch(
    plspm::plspm(data, pm, blocks, modes = unname(modes), boot.val = FALSE),
    error = function(e) stop("PLS-PM fit failed: ", conditionMessage(e), call. = FALSE)
  )
  required_paths <- c("Grazing->Biomass", "SoilNutrient->SR")
  path_keys <- paste(paths$from, paths$to, sep = "->")
  if (!all(required_paths %in% path_keys))
    stop("required R-042 structural path is absent", call. = FALSE)
  for (key in required_paths) {
    edge <- strsplit(key, "->", fixed = TRUE)[[1L]]
    coefficient <- fit$path_coefs[edge[2L], edge[1L]]
    if (length(coefficient) != 1L || !is.finite(coefficient))
      stop("required path coefficient is missing or non-finite: ", key, call. = FALSE)
  }
  check_r042_path_signs(fit$path_coefs)
  if (!is.finite(fit$gof)) fit$gof <- NA_real_
  fit
}

path_table_from_fit <- function(fit, paths, sign_rules) {
  out <- transform(paths,
    from = as.character(from), to = as.character(to),
    coefficient = mapply(function(f, t) fit$path_coefs[t, f], from, to),
    key = paste(from, to, sep = "->"))
  out$p_value <- mapply(function(f, t) {
    tab <- fit$inner_model[[t]]
    if (is.null(tab) || !f %in% rownames(tab)) return(NA_real_)
    as.numeric(tab[f, "Pr(>|t|)"])
  }, out$from, out$to)
  out$sign_expected <- unname(sign_rules[out$key]); out$sign_expected[is.na(out$sign_expected)] <- 0
  out$sign_pass <- out$sign_expected == 0 | sign(out$coefficient) == out$sign_expected
  out
}

compute_block_diagnostics <- function(fit, blocks, modes) {
  outer <- fit$outer_model
  rows <- lapply(names(blocks), function(block) {
    x <- outer[as.character(outer$block) == block, , drop = FALSE]
    loading <- as.numeric(x$loading)
    if (modes[[block]] == "A") {
      ave <- mean(loading^2)
      cr <- sum(loading)^2 / (sum(loading)^2 + sum(1 - loading^2))
      single <- block == "SR" && length(blocks[[block]]) == 1L
      data.frame(block = block, mode = "A", min_loading = min(loading),
                 AVE = if (single) NA_real_ else ave,
                 CR = if (single) NA_real_ else cr,
                 max_vif = NA_real_, redundancy = mean(x$redundancy),
                 weight_significance = NA_real_, direction_consistency = NA_real_,
                 applicability = if (single) "not applicable (single-indicator block)" else "reflective",
                 stringsAsFactors = FALSE)
    } else {
      vals <- as.data.frame(scale(fit$data[blocks[[block]]]), check.names = FALSE)
      vifs <- vapply(seq_along(vals), function(j) {
        if (ncol(vals) == 1L) return(1)
        r2 <- summary(stats::lm(vals[[j]] ~ ., data = vals[-j]))$r.squared
        if (r2 >= 1) Inf else 1 / (1 - r2)
      }, numeric(1))
      weights <- if ("weight" %in% names(x)) as.numeric(x$weight) else rep(NA_real_, nrow(x))
      direction_consistency <- if (all(is.na(weights))) NA_real_ else mean(sign(weights) == sign(loading))
      data.frame(block = block, mode = "B", min_loading = min(loading), AVE = NA_real_,
                 CR = NA_real_, max_vif = max(vifs), redundancy = mean(x$redundancy),
                 weight_significance = NA_real_,
                 direction_consistency = direction_consistency,
                 applicability = "formative; AVE/CR not applicable",
                 stringsAsFactors = FALSE)
    }
  })
  do.call(rbind, rows)
}

compute_srmr <- function(fit, paths) {
  scores <- as.data.frame(fit$scores)
  residuals <- list()
  for (to in unique(as.character(paths$to))) {
    incoming <- as.character(paths$from[as.character(paths$to) == to])
    if (!length(incoming)) next
    predicted <- as.matrix(scores[, incoming, drop = FALSE]) %*%
      as.numeric(fit$path_coefs[to, incoming, drop = TRUE])
    residuals[[to]] <- scores[[to]] - as.numeric(predicted)
  }
  if (!length(residuals)) return(NA_real_)
  sqrt(mean(unlist(lapply(residuals, function(x) x^2))))
}

compute_inner_vif <- function(fit, paths) {
  scores <- as.data.frame(fit$scores)
  endogenous <- unique(as.character(paths$to))
  do.call(rbind, lapply(endogenous, function(to) {
    predictors <- unique(as.character(paths$from[as.character(paths$to) == to]))
    if (length(predictors) == 1L) return(data.frame(block=to, predictor=predictors, VIF=1))
    values <- scores[predictors]
    do.call(rbind, lapply(predictors, function(predictor) {
      others <- setdiff(predictors, predictor)
      r2 <- summary(stats::lm(values[[predictor]] ~ ., data=values[others]))$r.squared
      data.frame(block=to, predictor=predictor, VIF=if (r2 >= 1) Inf else 1/(1-r2))
    }))
  }))
}

passes_r047_sr_r2 <- function(r2) {
  if (!is.numeric(r2) || length(r2) != 1L || !is.finite(r2)) return(FALSE)
  contract <- getOption("acceptance_contract", NULL)
  if (is.null(contract)) stop("acceptance contract required for R-PLS-r2", call. = FALSE)
  contract_pass(contract, "R-PLS-r2-lower", r2) && contract_pass(contract, "R-PLS-r2-upper", r2)
}

passes_r048_inner_vif <- function(vif) {
  if (!is.numeric(vif) || length(vif) == 0L || any(!is.finite(vif))) return(FALSE)
  contract <- getOption("acceptance_contract", NULL)
  if (is.null(contract)) stop("acceptance contract required for R-PLS-vif", call. = FALSE)
  all(vapply(vif, function(x) contract_pass(contract, "R-PLS-vif", x), logical(1)))
}

passes_r050_total_effects <- function(effect_values) {
  required <- c("Grazing", "SoilNutrient", "SoilEnv", "Biomass", "PlantDiv")
  if (!is.numeric(effect_values) || is.null(names(effect_values)) ||
      !all(required %in% names(effect_values)) ||
      any(!is.finite(effect_values[required]))) return(FALSE)
  contract <- getOption("acceptance_contract", NULL)
  if (is.null(contract)) stop("acceptance contract required for total effects", call. = FALSE)
  checks <- c(
    contract_pass(contract, "R-PLS-te-grazing-negative", as.numeric(effect_values[["Grazing"]] < 0)),
    contract_pass(contract, "R-PLS-te-nutrient-positive", as.numeric(effect_values[["SoilNutrient"]] > 0)),
    contract_pass(contract, "R-PLS-te-grazing-dominance-margin",
                  abs(effect_values[["Grazing"]]) - abs(effect_values[["SoilNutrient"]])),
    contract_pass(contract, "R-PLS-te-downstream-positive",
                  sum(effect_values[c("SoilEnv", "Biomass", "PlantDiv")] > 0))
  )
  all(checks)
}

evaluate_pls_acceptance_gates <- function(metrics,
    contract = getOption("acceptance_contract", NULL)) {
  if (is.null(contract)) stop("acceptance contract is required for PLS evaluation")
  required <- c("significant_paths", "grazing_significant_paths", "gof", "sr_r2",
                "inner_vif", "critical_bootstrap_p", "sign_violations",
                "complete_significant_paths", "effect_values", "formative_vif",
                "formative_direction", "plantdiv_loading", "plantdiv_ave", "plantdiv_cr")
  missing <- setdiff(required, names(metrics))
  if (length(missing)) stop("PLS metrics missing: ", paste(missing, collapse = ", "))
  pass <- function(id, value) {
    if (!is.numeric(value) || length(value) != 1L || !is.finite(value)) return(FALSE)
    contract_pass(contract, id, value)
  }
  total <- function(id, value) {
    if (!is.numeric(value) || length(value) != 1L || !is.finite(value)) return(FALSE)
    contract_pass(contract, id, value)
  }
  effects <- metrics$effect_values
  needed_effects <- c("Grazing", "SoilNutrient", "SoilEnv", "Biomass", "PlantDiv")
  if (!is.numeric(effects) || is.null(names(effects)) ||
      !all(needed_effects %in% names(effects)) || any(!is.finite(effects[needed_effects])))
    effects <- stats::setNames(rep(NA_real_, 5L), needed_effects)
  gates <- c(
    significant_paths = pass("R-PLS-paths", metrics$significant_paths),
    grazing_significant_paths = pass("R-PLS-grazing", metrics$grazing_significant_paths),
    gof = pass("R-PLS-gof", metrics$gof),
    sr_r2 = pass("R-PLS-r2-lower", metrics$sr_r2) && pass("R-PLS-r2-upper", metrics$sr_r2),
    inner_vif = is.numeric(metrics$inner_vif) && length(metrics$inner_vif) > 0L &&
      all(vapply(metrics$inner_vif, function(x) pass("R-PLS-vif", x), logical(1))),
    critical_bootstrap = is.numeric(metrics$critical_bootstrap_p) &&
      length(metrics$critical_bootstrap_p) > 0L &&
      all(vapply(metrics$critical_bootstrap_p, function(x) pass("R-PLS-bootstrap", x), logical(1))),
    path_sign = total("R-PLS-path-sign", as.numeric(metrics$sign_violations == 0L)),
    complete_paths = pass("R-PLS-complete-paths", metrics$complete_significant_paths),
    te_grazing_negative = total("R-PLS-te-grazing-negative", as.numeric(effects[["Grazing"]] < 0)),
    te_nutrient_positive = total("R-PLS-te-nutrient-positive", as.numeric(effects[["SoilNutrient"]] > 0)),
    te_dominance = total("R-PLS-te-grazing-dominance-margin",
      abs(effects[["Grazing"]]) - abs(effects[["SoilNutrient"]])),
    te_downstream = total("R-PLS-te-downstream-positive",
      sum(effects[c("SoilEnv", "Biomass", "PlantDiv")] > 0)),
    formative_vif = pass("R-PLS-formative-vif", metrics$formative_vif),
    formative_direction = pass("R-PLS-formative-direction", metrics$formative_direction),
    plantdiv_loading = pass("R-PLS-plantdiv-loading", metrics$plantdiv_loading),
    plantdiv_ave = pass("R-PLS-plantdiv-ave", metrics$plantdiv_ave),
    plantdiv_cr = pass("R-PLS-plantdiv-cr", metrics$plantdiv_cr)
  )
  c(gates, C7_pls = contract_boolean(contract, "C7_pls", all(gates)))
}

evaluate_complete_bootstrap_paths <- function(bootstrap_summary) {
  path_edges <- list(
    "Grazing->PlantDiv->SR" = c("Grazing->PlantDiv", "PlantDiv->SR"),
    "Grazing->Biomass->SR" = c("Grazing->Biomass", "Biomass->SR"),
    "Grazing->SoilEnv->SR" = c("Grazing->SoilEnv", "SoilEnv->SR"),
    "Grazing->SoilNutrient->SR" = c("Grazing->SoilNutrient", "SoilNutrient->SR"),
    "Grazing->SoilEnv->SoilNutrient->SR" = c("Grazing->SoilEnv", "SoilEnv->SoilNutrient", "SoilNutrient->SR"),
    "Grazing->Biomass->SoilNutrient->SR" = c("Grazing->Biomass", "Biomass->SoilNutrient", "SoilNutrient->SR"),
    "Grazing->SR" = "Grazing->SR"
  )
  if (!is.data.frame(bootstrap_summary) || !all(c("key", "p_value") %in% names(bootstrap_summary)))
    stop("bootstrap summary must contain key and p_value", call.=FALSE)
  path_results <- lapply(names(path_edges), function(path) {
    edges <- path_edges[[path]]
    rows <- match(edges, bootstrap_summary$key)
    p <- rep(NA_real_, length(edges))
    present <- !is.na(rows)
    p[present] <- bootstrap_summary$p_value[rows[present]]
    data.frame(path=path, edges=paste(edges, collapse=";"),
      edge_p_values=paste(ifelse(is.finite(p), format(p, digits=6), "NA"), collapse=";"),
      all_edges_significant=all(vapply(p, function(x) is.finite(x) && contract_gate("R-PLS-bootstrap", x), logical(1))), stringsAsFactors=FALSE)
  })
  do.call(rbind, path_results)
}

evaluate_critical_bootstrap_paths <- function(bootstrap_summary) {
  critical_keys <- c("Grazing->SoilNutrient", "SoilNutrient->SR",
                     "Grazing->Biomass", "Biomass->SR")
  if (!is.data.frame(bootstrap_summary) || !all(c("key", "p_value") %in% names(bootstrap_summary)))
    return(FALSE)
  rows <- match(critical_keys, bootstrap_summary$key)
  if (anyNA(rows)) return(FALSE)
  p <- bootstrap_summary$p_value[rows]
  all(vapply(p, function(x) is.finite(x) && contract_gate("R-PLS-bootstrap", x), logical(1)))
}

bootstrap_pls <- function(data, blocks, paths, modes = NULL, reps = 200L,
                          seed = NULL, sign_rules = pls_default_signs) {
  had_seed <- exists(".Random.seed", .GlobalEnv, inherits=FALSE)
  old_seed <- if (had_seed) get(".Random.seed", .GlobalEnv) else NULL
  on.exit(if (had_seed) assign(".Random.seed", old_seed, .GlobalEnv) else if (exists(".Random.seed", .GlobalEnv, inherits=FALSE)) rm(list=".Random.seed", envir=.GlobalEnv), add=TRUE)
  checked <- validate_pls_inputs(data, blocks, paths, modes)
  modes <- checked$modes
  if (length(reps) != 1L || !is.finite(reps) || reps < 2 || reps != as.integer(reps))
    stop("reps must be an integer >= 2", call. = FALSE)
  if (!is.null(seed)) set.seed(seed)
  keys <- paste(paths$from, paths$to, sep = "->")
  estimates <- matrix(NA_real_, reps, nrow(paths), dimnames = list(NULL, keys))
  original_fit <- fit_pls_once(data, blocks, paths, modes)
  weight_keys <- paste(original_fit$outer_model$block, original_fit$outer_model$name, sep="::")
  weight_estimates <- matrix(NA_real_, reps, length(weight_keys), dimnames=list(NULL, weight_keys))
  failures <- data.frame(replicate = integer(), message = character(), stringsAsFactors = FALSE)
  for (b in seq_len(reps)) {
    idx <- sample.int(nrow(data), nrow(data), replace = TRUE)
    one <- tryCatch(fit_pls_once(data[idx, , drop = FALSE], blocks, paths, modes),
                    error = function(e) e)
    if (inherits(one, "error")) {
      failures <- rbind(failures, data.frame(replicate = b, message = conditionMessage(one)))
    } else {
      estimates[b, ] <- mapply(function(f, t) one$path_coefs[t, f], paths$from, paths$to)
      one_keys <- paste(one$outer_model$block, one$outer_model$name, sep="::")
      matched <- match(weight_keys, one_keys)
      weight_estimates[b, ] <- one$outer_model$weight[matched]
    }
  }
  valid <- stats::complete.cases(estimates)
  if (!any(valid)) stop("all bootstrap replicates failed; see failure diagnostics", call. = FALSE)
  summary <- do.call(rbind, lapply(seq_len(ncol(estimates)), function(j) {
    z <- estimates[valid, j]; original <- original_fit$path_coefs[paths$to[j], paths$from[j]]
    se <- stats::sd(z)
    data.frame(key = keys[j], Original = original, Mean.Boot = mean(z), Std.Error = se,
               perc.025 = unname(stats::quantile(z, .025)), perc.975 = unname(stats::quantile(z, .975)),
               p_value = if (is.finite(se) && se > 0) 2 * stats::pnorm(-abs(original / se)) else NA_real_,
               stringsAsFactors = FALSE)
  }))
  weight_valid <- stats::complete.cases(weight_estimates)
  weight_summary <- do.call(rbind, lapply(seq_along(weight_keys), function(j) {
    z <- weight_estimates[weight_valid, j]
    se <- stats::sd(z)
    data.frame(key=weight_keys[j], original=original_fit$outer_model$weight[j],
      bootstrap_mean=mean(z), std_error=se,
      p_value=if (is.finite(se) && se > 0) 2 * stats::pnorm(-abs(mean(z) / se)) else NA_real_,
      stringsAsFactors=FALSE)
  }))
  list(estimates = estimates, summary = summary, failures = failures,
       weight_estimates=weight_estimates, weight_summary=weight_summary,
       attempted = reps, successful = sum(valid), failed = sum(!valid))
}

fit_pls_full <- function(data, blocks, paths, modes = NULL, bootstrap_reps = 200L,
                         seed = NULL, sign_rules = pls_default_signs) {
  had_seed <- exists(".Random.seed", .GlobalEnv, inherits=FALSE)
  old_seed <- if (had_seed) get(".Random.seed", .GlobalEnv) else NULL
  on.exit(if (had_seed) assign(".Random.seed", old_seed, .GlobalEnv) else if (exists(".Random.seed", .GlobalEnv, inherits=FALSE)) rm(list=".Random.seed", envir=.GlobalEnv), add=TRUE)
  checked <- validate_pls_inputs(data, blocks, paths, modes); modes <- checked$modes
  if (!is.null(seed)) set.seed(seed)
  fit <- fit_pls_once(data, blocks, paths, modes)
  sign_rules <- c(pls_default_signs, sign_rules)
  paths_out <- path_table_from_fit(fit, paths, sign_rules)
  diagnostics <- compute_block_diagnostics(fit, blocks, modes)
  inner_vif <- compute_inner_vif(fit, paths)
  bootstrap <- bootstrap_pls(data, blocks, paths, modes, reps = bootstrap_reps,
                             seed = seed, sign_rules = sign_rules)
  diagnostics$weight_significant_share <- NA_real_
  for (i in which(diagnostics$mode == "B")) {
    prefix <- paste0(diagnostics$block[i], "::")
    rows <- startsWith(bootstrap$weight_summary$key, prefix)
    p <- bootstrap$weight_summary$p_value[rows]
    diagnostics$weight_significant_share[i] <- if (length(p)) mean(p < .05, na.rm=TRUE) else NA_real_
  }
  r2 <- stats::setNames(fit$inner_summary$R2, rownames(fit$inner_summary))
  critical_bootstrap_pass <- evaluate_critical_bootstrap_paths(bootstrap$summary)
  significant <- is.finite(paths_out$p_value) & paths_out$p_value < .05
  grazing_significant <- significant & paths_out$from == "Grazing"
  complete_path_results <- evaluate_complete_bootstrap_paths(bootstrap$summary)
  complete_path_pass <- complete_path_results$all_edges_significant
  effects <- fit$effects
  effect_sr <- effects[grepl(" -> SR$", effects$relationships), , drop = FALSE]
  effect_values <- stats::setNames(effect_sr$total, sub(" -> SR$", "", effect_sr$relationships))
  total_effect_order <- passes_r050_total_effects(effect_values)
  srmr <- compute_srmr(fit, paths)
  plantdiv <- diagnostics[diagnostics$block == "PlantDiv", , drop = FALSE]
  formative <- diagnostics[diagnostics$mode == "B", , drop = FALSE]
  critical_keys <- c("Grazing->SoilNutrient", "SoilNutrient->SR",
                     "Grazing->Biomass", "Biomass->SR")
  critical_rows <- match(critical_keys, bootstrap$summary$key)
  critical_p <- if (anyNA(critical_rows)) NA_real_ else bootstrap$summary$p_value[critical_rows]
  metrics <- list(
    significant_paths = sum(significant), grazing_significant_paths = sum(grazing_significant),
    gof = fit$gof, sr_r2 = if ("SR" %in% names(r2)) unname(r2[["SR"]]) else NA_real_,
    inner_vif = inner_vif$VIF, critical_bootstrap_p = critical_p,
    sign_violations = sum(!paths_out$sign_pass),
    complete_significant_paths = sum(complete_path_pass), effect_values = effect_values,
    formative_vif = if (nrow(formative)) max(formative$max_vif, na.rm = TRUE) else NA_real_,
    formative_direction = if (nrow(formative)) min(formative$direction_consistency, na.rm = TRUE) else NA_real_,
    plantdiv_loading = if (nrow(plantdiv)) plantdiv$min_loading[[1L]] else NA_real_,
    plantdiv_ave = if (nrow(plantdiv)) plantdiv$AVE[[1L]] else NA_real_,
    plantdiv_cr = if (nrow(plantdiv)) plantdiv$CR[[1L]] else NA_real_
  )
  gates <- evaluate_pls_acceptance_gates(metrics)
  list(fit = fit, paths = paths_out, diagnostics = diagnostics, inner_vif=inner_vif, bootstrap = bootstrap,
       r2 = r2, gof = fit$gof, srmr = srmr, effects = effects,
       significant_paths = sum(significant), grazing_significant_paths = sum(grazing_significant),
       complete_significant_paths = complete_path_pass,
       complete_path_results = complete_path_results, gates = gates,
       has_pass = isTRUE(all(gates)),
       sign_violations = sum(!paths_out$sign_pass))
}

evaluate_pls_full <- function(result) {
  if (!is.list(result) || is.null(result$paths) || is.null(result$bootstrap) || is.null(result$gates))
    stop("result must be returned by fit_pls_full", call. = FALSE)
  list(gates = result$gates, has_pass = result$has_pass,
       failed_bootstrap = result$bootstrap$failed,
       sign_violations = result$sign_violations,
       critical_paths = result$bootstrap$summary[result$bootstrap$summary$key %in%
                                                   c("Grazing->Biomass", "SoilNutrient->SR"), , drop = FALSE])
}
