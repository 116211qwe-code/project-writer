# Module 5: fast, transparent PLS-style proxy and error calibration.
# Teaching simulation only; this is a gate, not a substitute for full PLS-PM.

`%||%` <- function(x, y) if (is.null(x)) y else x

proxy_search_gate_status <- function(calibration_pass) {
  if (isTRUE(calibration_pass)) "READY" else "BLOCKED"
}

pls_fast_proxy <- function(data, blocks, paths,
                           sign_rules = c("Grazing->Biomass" = -1,
                                           "SoilNutrient->SR" = 1)) {
  if (!is.data.frame(data) || !length(blocks) || !is.data.frame(paths))
    stop("data, blocks, and paths must be supplied", call. = FALSE)
  if (!all(c("from", "to") %in% names(paths)))
    stop("paths must contain from and to columns", call. = FALSE)
  if (nrow(data) < 3L) stop("at least three observations are required", call. = FALSE)
  used <- unique(unlist(blocks, use.names = FALSE))
  missing <- setdiff(used, names(data))
  if (length(missing)) stop("missing block indicators: ", paste(missing, collapse = ", "), call. = FALSE)
  if (!all(vapply(data[used], is.numeric, logical(1))) || any(!is.finite(as.matrix(data[used]))))
    stop("block indicators must be finite numeric columns", call. = FALSE)
  scores <- lapply(blocks, function(cols) {
    x <- as.matrix(data[cols])
    if (any(apply(x, 2L, stats::sd) == 0)) stop("zero-variance block indicator", call. = FALSE)
    as.numeric(scale(rowMeans(x)))
  })
  scores <- as.data.frame(scores, check.names = FALSE)
  if (any(!is.finite(as.matrix(scores)))) stop("latent score construction failed", call. = FALSE)

  path_rows <- lapply(seq_len(nrow(paths)), function(i) {
    from <- as.character(paths$from[i]); to <- as.character(paths$to[i])
    if (!all(c(from, to) %in% names(scores))) stop("path refers to unknown block", call. = FALSE)
    incoming <- as.character(paths$from[as.character(paths$to) == to])
    if (!length(incoming)) stop("path target has no predictors", call. = FALSE)
    fit <- stats::lm(scores[[to]] ~ ., data = scores[incoming], model = FALSE)
    co <- stats::coef(fit)[from]
    r2 <- summary(fit)$r.squared
    pred <- scores[incoming]
    vifs <- vapply(seq_along(incoming), function(j) {
      if (length(incoming) == 1L) return(1)
      others <- pred[-j]
      if (!ncol(others)) return(1)
      1 / (1 - summary(stats::lm(pred[[j]] ~ ., data = others))$r.squared)
    }, numeric(1))
    data.frame(from = from, to = to, coefficient = unname(co), r2 = r2,
               max_vif = max(vifs), stringsAsFactors = FALSE)
  })
  path_table <- do.call(rbind, path_rows)
  loadings <- do.call(rbind, lapply(names(blocks), function(block) {
    z <- scores[[block]]
    data.frame(block = block, indicator = blocks[[block]],
               loading = vapply(blocks[[block]], function(nm) stats::cor(data[[nm]], z), numeric(1)),
               stringsAsFactors = FALSE)
  }))
  measurement <- do.call(rbind, lapply(split(loadings, loadings$block), function(x) {
    l <- x$loading; ave <- mean(l^2); cr <- sum(l)^2 / (sum(l)^2 + sum(1 - l^2))
    data.frame(block = x$block[1], min_loading = min(l), AVE = ave, CR = cr,
               stringsAsFactors = FALSE)
  }))
  path_table$path_abs <- abs(path_table$coefficient)
  path_table$sign_expected <- unname(sign_rules[paste(path_table$from, path_table$to, sep = "->")])
  path_table$sign_expected[is.na(path_table$sign_expected)] <- 0
  path_table$sign_pass <- path_table$sign_expected == 0 |
    sign(path_table$coefficient) == path_table$sign_expected
  list(scores = scores, paths = path_table, measurement = measurement,
       measurement_diagnostics = loadings,
       sign_violations = sum(!path_table$sign_pass),
       max_vif = max(path_table$max_vif),
       r2 = stats::setNames(path_table$r2, paste(path_table$from, path_table$to, sep = "->")))
}

evaluate_fast <- function(proxy, targets = list(min_loading = 0.4, AVE = 0.5, CR = 0.7,
                                                  r2_min = 0.5, r2_max = 0.8, max_vif = 5)) {
  if (!is.list(proxy) || is.null(proxy$paths) || is.null(proxy$measurement))
    stop("proxy must be the result of pls_fast_proxy", call. = FALSE)
  m <- proxy$measurement; p <- proxy$paths
  residual <- c(
    loading_deficit = sum(pmax(0, targets$min_loading - m$min_loading)),
    AVE_deficit = sum(pmax(0, targets$AVE - m$AVE)),
    CR_deficit = sum(pmax(0, targets$CR - m$CR)),
    r2_deficit = sum(pmax(0, targets$r2_min - p$r2) + pmax(0, p$r2 - targets$r2_max)),
    vif_excess = sum(pmax(0, p$max_vif - targets$max_vif)),
    sign_violations = proxy$sign_violations
  )
  gates <- c(has_measurement = all(m$min_loading >= targets$min_loading & m$AVE >= targets$AVE & m$CR >= targets$CR),
             has_r2 = all(p$r2 >= targets$r2_min & p$r2 <= targets$r2_max),
             has_vif = all(p$max_vif < targets$max_vif), has_sign = proxy$sign_violations == 0)
  # Keep the compact residual vector stable for existing callers, while
  # exposing every diagnostic required for calibration and audit reports.
  measurement_residual <- c(
    loading = max(pmax(0, targets$min_loading - m$min_loading)),
    AVE = max(pmax(0, targets$AVE - m$AVE)),
    CR = max(pmax(0, targets$CR - m$CR)))
  path_residual <- setNames(pmax(0, targets$r2_min - p$r2) + pmax(0, p$r2 - targets$r2_max),
                            paste(p$from, p$to, sep = "->"))
  list(gates = gates, residual = residual, max_residual = max(residual),
       measurement_residual = measurement_residual, path_residual = path_residual,
       path_coefficients = setNames(p$coefficient, paste(p$from, p$to, sep = "->")),
       has_pass = all(gates), score = -sum(residual))
}

calibrate_pls_proxy <- function(fast, full) {
  if (is.data.frame(fast) && is.data.frame(full) && (nrow(fast) != nrow(full) || nrow(fast) < 2L))
    stop("calibration requires at least two paired rows", call. = FALSE)
  required <- c("r2", "max_vif", "sign_violations", "loading", "AVE", "CR", "path_coefficient")
  if (!is.data.frame(fast) || !is.data.frame(full) || !all(required %in% names(fast)) || !all(required %in% names(full)))
    stop("fast and full must contain r2, max_vif, and sign_violations", call. = FALSE)
  if (any(!is.finite(as.matrix(fast[required]))) ||
      any(!is.finite(as.matrix(full[required]))))
    stop("calibration values must be finite", call. = FALSE)
  direction <- vapply(required, function(nm) {
    a <- fast[[nm]] - mean(fast[[nm]]); b <- full[[nm]] - mean(full[[nm]])
    if (all(a == 0) && all(b == 0)) return(all(fast[[nm]] == full[[nm]]))
    if (all(a == 0) || all(b == 0)) return(FALSE)
    stats::cor(a, b) >= 0
  }, logical(1))
  bias <- vapply(required, function(nm) mean(fast[[nm]] - full[[nm]]), numeric(1))
  spearman <- vapply(required, function(nm) {
    if (length(unique(fast[[nm]])) == 1L && length(unique(full[[nm]])) == 1L)
      return(if (all(fast[[nm]] == full[[nm]])) 1 else 0)
    value <- suppressWarnings(stats::cor(fast[[nm]], full[[nm]], method = "spearman"))
    if (is.finite(value)) value else 0
  }, numeric(1))
  max_abs_error <- vapply(required, function(nm) max(abs(fast[[nm]] - full[[nm]])), numeric(1))
  tolerance <- c(r2 = .05, max_vif = .5, sign_violations = 0, loading = .05,
                 AVE = .05, CR = .05, path_coefficient = .05)
  tolerance <- tolerance[required]
  fast_score <- if ("score" %in% names(fast)) fast$score else rowMeans(fast[required])
  full_score <- if ("score" %in% names(full)) full$score else rowMeans(full[required])
  rank_spearman <- suppressWarnings(stats::cor(fast_score, full_score, method = "spearman"))
  if (!is.finite(rank_spearman)) rank_spearman <- 0
  gate_cols <- intersect(c("gate_pass", "pass"), intersect(names(fast), names(full)))
  gate_value <- function(x) if (length(gate_cols)) as.logical(x[[gate_cols[1]]]) else
    x$r2 > .5 & x$r2 <= .8 & x$max_vif < 5 & x$sign_violations == 0 &
      x$loading > .4 & x$AVE > .5 & x$CR > .7
  gate_agreement <- mean(gate_value(fast) == gate_value(full))
  list(bias = bias, max_abs_error = max_abs_error, tolerance = tolerance,
       error_pass = all(max_abs_error <= tolerance), direction = direction,
       direction_pass = all(direction, na.rm = TRUE),
       spearman = spearman, spearman_min = min(spearman),
       rank_spearman = rank_spearman, gate_agreement = gate_agreement,
       n = nrow(fast), calibration_pass = nrow(fast) >= 100L &&
         all(direction, na.rm = TRUE) && all(spearman >= .8, na.rm = TRUE) &&
         rank_spearman >= .8 && gate_agreement >= .9 && all(max_abs_error <= tolerance),
       reason = if (nrow(fast) >= 100L && all(direction, na.rm = TRUE) && all(spearman >= .8, na.rm = TRUE) && rank_spearman >= .8 && gate_agreement >= .9 && all(max_abs_error <= tolerance)) "CALIBRATED" else "PROXY_NOT_CALIBRATED")
}

calibrate_pls_proxy_set <- function(candidates, proxy_fn, full_fn,
                                   output_path = "proxy_calibration.csv") {
  if (!is.list(candidates) || length(candidates) != 100L) stop("R-098 校准必须正好包含 100 个完整生成解")
  rows <- lapply(seq_along(candidates), function(i) {
    tryCatch({
    fast <- proxy_fn(candidates[[i]])
    full <- full_fn(candidates[[i]])
    data.frame(candidate_id=i, error="",
      fast_r2 = mean(as.numeric(fast$r2), na.rm = TRUE),
      full_r2 = mean(as.numeric(full$r2), na.rm = TRUE),
      fast_max_vif = fast$max_vif, full_max_vif = max(full$diagnostics$max_vif, na.rm = TRUE),
      fast_sign_violations = fast$sign_violations,
      full_sign_violations = full$sign_violations,
      fast_loading = min(fast$measurement$min_loading[fast$measurement$block == "PlantDiv"], na.rm = TRUE),
      full_loading = min(full$diagnostics$min_loading, na.rm = TRUE),
      fast_AVE = fast$measurement$AVE[fast$measurement$block == "PlantDiv"][1],
      full_AVE = full$diagnostics$AVE[full$diagnostics$block == "PlantDiv"][1],
      fast_CR = fast$measurement$CR[fast$measurement$block == "PlantDiv"][1],
      full_CR = full$diagnostics$CR[full$diagnostics$block == "PlantDiv"][1],
      fast_path_coefficient = unname(fast$path_coefficients["SoilNutrient->SR"] %||% 0),
      full_path_coefficient = unname(full$paths$coefficient[full$paths$key == "SoilNutrient->SR"][1]),
      stringsAsFactors = FALSE)
    }, error=function(e) data.frame(candidate_id=i, error=conditionMessage(e),
      fast_r2=NA_real_, full_r2=NA_real_, fast_max_vif=NA_real_, full_max_vif=NA_real_,
      fast_sign_violations=NA_real_, full_sign_violations=NA_real_, fast_loading=NA_real_, full_loading=NA_real_,
      fast_AVE=NA_real_, full_AVE=NA_real_, fast_CR=NA_real_, full_CR=NA_real_,
      fast_path_coefficient=NA_real_, full_path_coefficient=NA_real_, stringsAsFactors=FALSE))
  })
  paired <- do.call(rbind, rows)
  required_paired <- c("fast_r2", "full_r2", "fast_max_vif", "full_max_vif",
                       "fast_sign_violations", "full_sign_violations",
                       "fast_loading", "full_loading", "fast_AVE", "full_AVE",
                       "fast_CR", "full_CR", "fast_path_coefficient", "full_path_coefficient")
  complete <- !nzchar(paired$error) & stats::complete.cases(paired[required_paired]) &
    apply(paired[required_paired], 1, function(x) all(is.finite(as.numeric(x))))
  paired_complete <- paired[complete, , drop = FALSE]
  if (nrow(paired_complete) < 2L) {
    paired$sample_id <- paired$candidate_id
    paired$proxy_score <- paired$full_score <- paired$spearman <- paired$gate_agreement <- NA_real_
    paired$metric_errors <- vapply(seq_len(nrow(paired)), function(i) {
      if (nzchar(paired$error[i])) return(paired$error[i])
      missing <- required_paired[!vapply(paired[i, required_paired, drop = FALSE],
        function(x) is.finite(as.numeric(x)), logical(1))]
      paste0("insufficient finite paired PLS calibration rows; missing/nonfinite=",
        paste(missing, collapse = ","))
    }, character(1))
    paired$pass <- "FALSE"
    calibration_rows <- paired[c("sample_id", "proxy_score", "full_score", "spearman",
      "gate_agreement", "metric_errors", "pass")]
    utils::write.csv(calibration_rows, output_path, row.names = FALSE, na = "")
    return(list(paired = paired, calibration = list(n = nrow(paired_complete), attempted=length(candidates),
      successful=nrow(paired_complete), calibration_pass = FALSE, error_pass=FALSE,
      direction_pass=FALSE, reason = "insufficient finite paired PLS calibration rows",
      failures=paired[!complete, c("candidate_id", "error"), drop=FALSE])))
  }
  calibration <- calibrate_pls_proxy(
    setNames(paired_complete[c("fast_r2", "fast_max_vif", "fast_sign_violations", "fast_loading", "fast_AVE", "fast_CR", "fast_path_coefficient")],
             c("r2", "max_vif", "sign_violations", "loading", "AVE", "CR", "path_coefficient")),
    setNames(paired_complete[c("full_r2", "full_max_vif", "full_sign_violations", "full_loading", "full_AVE", "full_CR", "full_path_coefficient")],
             c("r2", "max_vif", "sign_violations", "loading", "AVE", "CR", "path_coefficient"))
  )
  calibration$finite_pairs <- nrow(paired_complete)
  calibration$attempted <- length(candidates)
  calibration$successful <- nrow(paired_complete)
  calibration$failures <- paired[!complete, c("candidate_id", "error"), drop=FALSE]
  calibration$calibration_pass <- isTRUE(calibration$calibration_pass) && nrow(paired_complete) == length(candidates)
  proxy_score <- -apply(paired_complete[c("fast_r2", "fast_max_vif", "fast_loading", "fast_AVE", "fast_CR")], 1, function(x) sum(abs(x)))
  full_score <- -apply(paired_complete[c("full_r2", "full_max_vif", "full_loading", "full_AVE", "full_CR")], 1, function(x) sum(abs(x)))
  proxy_gate <- paired_complete$fast_r2 > .5 & paired_complete$fast_r2 <= .8 & paired_complete$fast_max_vif < 5 & paired_complete$fast_sign_violations == 0
  full_gate <- paired_complete$full_r2 > .5 & paired_complete$full_r2 <= .8 & paired_complete$full_max_vif < 5 & paired_complete$full_sign_violations == 0
  gate_agreement <- mean(proxy_gate == full_gate)
  spearman <- suppressWarnings(stats::cor(proxy_score, full_score, method = "spearman"))
  if (!is.finite(spearman)) spearman <- 0
  calibration$rank_spearman <- spearman
  calibration$gate_agreement <- gate_agreement
  calibration$calibration_pass <- calibration$calibration_pass && spearman >= .8 && gate_agreement >= .9
  paired$proxy_score <- paired$full_score <- paired$spearman <- paired$gate_agreement <- NA_real_
  paired$metric_errors <- paired$pass <- NA_character_
  paired$proxy_score[complete] <- proxy_score
  paired$full_score[complete] <- full_score
  paired$spearman[complete] <- spearman
  paired$gate_agreement[complete] <- gate_agreement
  paired$metric_errors[complete] <- vapply(which(complete), function(i) paste(sprintf("%s=%.6g", c("r2","max_vif","loading","AVE","CR","path_coefficient"), abs(unlist(paired[i, c("fast_r2","fast_max_vif","fast_loading","fast_AVE","fast_CR","fast_path_coefficient")])-unlist(paired[i, c("full_r2","full_max_vif","full_loading","full_AVE","full_CR","full_path_coefficient")]))), collapse=";"), character(1))
  paired$pass[complete] <- as.character(calibration$calibration_pass)
  paired$sample_id <- paired$candidate_id
  paired$pass[!complete] <- "FALSE"
  calibration_rows <- paired[c("sample_id", "proxy_score", "full_score", "spearman",
    "gate_agreement", "metric_errors", "pass")]
  utils::write.csv(calibration_rows, output_path, row.names=FALSE, na="")
  calibration$reason <- if (calibration$calibration_pass) "CALIBRATED" else "PROXY_NOT_CALIBRATED"
  list(paired = paired, calibration = calibration)
}
