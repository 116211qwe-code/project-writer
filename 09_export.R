# Module 9: deterministic export and artifact validation.

open_plot_device <- function(path, width = 8, height = 6, type = c("emf", "png", "pdf")) {
  type <- match.arg(type)
  if (type == "emf") {
    if (!requireNamespace("devEMF", quietly = TRUE)) stop("devEMF is required for EMF export")
    devEMF::emf(file = path, width = width, height = height, emfPlus = TRUE)
  } else if (type == "png") {
    grDevices::png(path, width = width, height = height, units = "in", res = 300)
  } else grDevices::pdf(path, width = width, height = height, useDingbats = FALSE)
}

write_plot_bundle <- function(path_stem, plot_fun, width = 8, height = 6) {
  for (ext in c("emf", "png", "pdf")) {
    p <- paste0(path_stem, ".", ext)
    open_plot_device(p, width, height, ext)
    tryCatch(plot_fun(), finally = grDevices::dev.off())
    if (!file.exists(p) || file.info(p)$size <= 0) stop("plot export failed: ", p)
  }
}

plot_pls_paths <- function(results) {
  paths <- if (!is.null(results$full$result$paths)) results$full$result$paths else NULL
  proxy_only <- is.null(paths) && !is.null(results$fast$proxy$paths)
  if (proxy_only) paths <- results$fast$proxy$paths
  if (is.null(paths) || !nrow(paths)) { plot.new(); title("PLS path diagram"); text(.5, .5, "PLS not available"); return(invisible()) }
  nodes <- unique(c(paths$from, paths$to))
  layers <- list(c("Grazing"), c("PlantDiv", "Biomass", "SoilEnv", "SoilNutrient"), c("SR"))
  x <- setNames(rep(NA_real_, length(nodes)), nodes)
  y <- x
  for (i in seq_along(layers)) {
    present <- intersect(layers[[i]], nodes)
    x[present] <- c(.14, .50, .86)[i]
    y[present] <- if (length(present) == 1L) .5 else seq(.79, .21, length.out = length(present))
  }
  unplaced <- names(x)[!is.finite(x)]
  if (length(unplaced)) {
    x[unplaced] <- .5; y[unplaced] <- seq(.8, .2, length.out = length(unplaced))
  }
  plot.new(); plot.window(xlim = c(.02, .98), ylim = c(.08, .92), xaxs = "i", yaxs = "i")
  if (proxy_only) {
    title("PLS fast-proxy path diagram",
          sub = "Proxy diagnostic only; full PLS/Bootstrap skipped because the proxy gate failed",
          cex.sub = .68, col.sub = "firebrick")
  } else title("PLS path diagram")
  for (i in seq_len(nrow(paths))) {
    a <- paths$from[i]; b <- paths$to[i]; co <- paths$coefficient[i]
    p <- if ("p_value" %in% names(paths)) paths$p_value[i] else NA_real_
    col <- if (!is.finite(p) || p >= .05) "grey60" else if (co > 0) "red3" else "forestgreen"
    dx <- x[b] - x[a]; dy <- y[b] - y[a]
    x0 <- x[a] + sign(dx) * .075; y0 <- y[a] + dy * .10
    x1 <- x[b] - sign(dx) * .075; y1 <- y[b] - dy * .10
    arrows(x0, y0, x1, y1, col = col, lwd = if (is.finite(p) && p < .05) 1.6 else 1,
           lty = if (is.finite(p) && p < .05) 1 else 3, length = .065, angle = 20)
    stars <- if (!is.finite(p)) "" else if (p < .001) "***" else if (p < .01) "**" else if (p < .05) "*" else ""
    label_x <- (x0 + x1) / 2 + if (i %% 2) .012 else -.012
    label_y <- (y0 + y1) / 2 + if (i %% 2) .035 else -.035
    text(label_x, label_y, sprintf("%.2f%s", co, stars), col = col, cex = .68, font = 2)
  }
  for (node in nodes) {
    rect(x[node] - .075, y[node] - .035, x[node] + .075, y[node] + .035,
         col = "white", border = "grey25", lwd = 1.2)
    text(x[node], y[node], node, cex = .62)
  }
  legend("bottomleft", legend = c("Not significant", "Positive significant", "Negative significant"),
         col = c("grey60", "red3", "forestgreen"), lty = 1, bty = "n", cex = .75)
  attr(nodes, "plot_source") <- if (proxy_only) "fast_proxy" else "full_pls"
  invisible(nodes)
}

plot_mantel <- function(results, new_plot = TRUE) {
  tab <- results$acceptance$mantel
  if (is.null(tab) || !nrow(tab)) { plot.new(); title("Mantel: SR versus soil indicators"); text(.5, .5, "Mantel not available"); return(invisible()) }
  r <- tab$r %||% tab$correlation %||% rep(NA_real_, nrow(tab)); labs <- tab$indicator
  plot(r, seq_along(r), xlim = c(-1, 1), yaxt = "n", xlab = "Mantel r", ylab = "", pch = NA,
       col = ifelse(r >= 0, "red3", "royalblue3")); axis(2, at = seq_along(r), labels = labs, las = 2, cex.axis = .7)
  title("Mantel: SR versus soil indicators")
  abline(v = 0, lty = 2)
  for (i in seq_along(r)) {
    col <- if (r[i] >= 0) "red3" else "royalblue3"
    segments(0, i, r[i], i, col=col, lwd=3)
    points(r[i], i, pch=19, col=col)
    text(r[i], i, sprintf("%.2f", r[i]), pos=if (r[i] >= 0) 4 else 2, cex=.7)
  }
  legend("bottom", inset = -.24, xpd = TRUE, legend = c("Positive", "Negative", "0 reference"),
         col = c("red3", "royalblue3", "grey30"), pch = c(19, 19, NA), lty = c(NA, NA, 2),
         horiz = TRUE, bty = "n", cex = .8)
}

plot_pearson_heatmap <- function(results) {
  tab <- results$acceptance$pearson
  if (is.null(tab) || !nrow(tab)) { plot.new(); title("Vegetation and SR Pearson correlations"); text(.5, .5, "Pearson not available"); return(invisible()) }
  vars <- c("AGB", "BGB", "coverage", "LB", "H", "E", "D", "R", "SR")
  vars <- vars[vars %in% unique(c(tab$variable_1, tab$variable_2))]
  r <- p <- matrix(NA_real_, length(vars), length(vars), dimnames = list(vars, vars))
  diag(r) <- 1; diag(p) <- 0
  for (i in seq_len(nrow(tab))) {
    a <- tab$variable_1[i]; b <- tab$variable_2[i]
    if (a %in% vars && b %in% vars) { r[a,b] <- r[b,a] <- tab$r[i]; p[a,b] <- p[b,a] <- tab$p_value[i] }
  }
  # Keep the requested domain order, clustering within plant/diversity blocks.
  domains <- list(植物 = c("AGB", "BGB", "coverage", "LB"), 多样性 = c("H", "E", "D", "R", "SR"))
  order_vars <- unlist(lapply(domains, function(group) {
    group <- intersect(group, vars)
    if (length(group) < 3L) return(group)
    d <- as.dist(1 - abs(r[group, group, drop = FALSE])); group[hclust(d)$order]
  }), use.names = FALSE)
  order_vars <- unique(c(order_vars, setdiff(vars, order_vars)))
  r <- r[order_vars, order_vars, drop = FALSE]; p <- p[order_vars, order_vars, drop = FALSE]
  n <- length(order_vars)
  plot.new(); plot.window(xlim = c(.5, n + .5), ylim = c(.5, n + .5), xaxs = "i", yaxs = "i", asp = 1)
  title("Vegetation and SR Pearson correlations")
  palette <- grDevices::colorRampPalette(c("#62B8D1", "white", "#E56B6F"))(101)
  for (i in seq_len(n)) for (j in seq_len(n)) if (i > j) {
    value <- r[order_vars[i], order_vars[j]]; pv <- p[order_vars[i], order_vars[j]]
    rect(j - .45, n - i + .55, j + .45, n - i + 1.45, border = "grey75", col = "white")
    if (is.finite(value) && is.finite(pv) && pv < .05) {
      color <- palette[round((value + 1) * 50) + 1]
      symbols(j, n - i + 1, circles = .34 * sqrt(abs(value)), inches = FALSE,
              add = TRUE, bg = color, fg = "grey35")
      stars <- if (pv < .001) "***" else if (pv < .01) "**" else "*"
      text(j, n - i + 1, stars, cex = .62, font = 2)
    }
  }
  for (i in seq_len(n)) {
    text(i, n + .05, order_vars[i], srt = 45, adj = 0, xpd = TRUE, cex = .72)
    text(.35, n - i + 1, order_vars[i], adj = 1, xpd = TRUE, cex = .72)
  }
  legend("bottom", inset = -.16, xpd = TRUE, legend = c("positive", "negative", "blank = p >= 0.05"),
         pch = c(21,21,NA), pt.bg = c("#E56B6F", "#62B8D1", NA), col = "grey35",
         horiz = TRUE, bty = "n", cex = .72)
}

compact_lsd_letters <- function(indicators, groups, lsd) {
  out <- matrix("", nrow = length(indicators), ncol = length(groups),
                dimnames = list(indicators, as.character(groups)))
  for (v in indicators) {
    pairs <- lsd[lsd$indicator == v, , drop = FALSE]
    nonsig <- matrix(TRUE, length(groups), length(groups), dimnames = list(groups, groups))
    for (i in seq_len(nrow(pairs))) {
      a <- as.character(pairs$group_1[i]); b <- as.character(pairs$group_2[i])
      if (a %in% groups && b %in% groups) nonsig[a,b] <- nonsig[b,a] <- !pairs$significant[i]
    }
    subsets <- lapply(seq_len(2^length(groups)-1L), function(mask)
      which(as.logical(intToBits(mask)[seq_along(groups)])))
    subsets <- subsets[vapply(subsets, function(ix) all(nonsig[ix, ix, drop = FALSE]), logical(1))]
    maximal <- subsets[!vapply(subsets, function(ix) any(vapply(subsets, function(jx) length(jx) > length(ix) && all(ix %in% jx), logical(1))), logical(1))]
    alphabet <- base::letters[seq_along(maximal)]
    for (j in seq_along(groups)) {
      assigned <- alphabet[vapply(maximal, function(ix) j %in% ix, logical(1))]
      out[v, j] <- paste(assigned, collapse = "")
    }
  }
  out
}

plot_anova_bars <- function(results) {
  tab <- results$acceptance$anova
  data <- results$statistical_inputs$anova
  if (is.null(tab) || is.null(data) || !nrow(tab)) { plot.new(); title("July indicator means by grazing (mean +/- SE)"); text(.5, .5, "ANOVA data not available"); return(invisible()) }
  indicators <- tab$indicator
  groups <- sort(unique(data$grazing))
  means <- ses <- matrix(NA_real_, length(indicators), length(groups),
                         dimnames = list(indicators, groups))
  for (i in seq_along(indicators)) for (j in seq_along(groups)) {
    x <- data[[indicators[i]]][data$grazing == groups[j]]
    means[i,j] <- mean(x); ses[i,j] <- stats::sd(x)/sqrt(length(x))
  }
  cld <- compact_lsd_letters(indicators, groups, results$acceptance$lsd)
  plot.new(); plot.window(xlim = c(0, 3), ylim = c(0, 6), xaxs = "i", yaxs = "i")
  title("July ANOVA by grazing: mean +/- SE; LSD letters")
  cols <- c("#282828", "#888888", "#D2D2D2")
  for (k in seq_along(indicators)) {
    row <- (k - 1L) %/% 3L; col <- (k - 1L) %% 3L
    left <- col + .08; bottom <- 5 - row + .08; w <- .84; h <- .82
    rect(left, bottom, left + w, bottom + h, border = "grey65")
    v <- indicators[k]; center <- means[k, ]; err <- ses[k, ]
    lo <- min(0, center - err); hi <- max(center + err); span <- if (hi == lo) 1 else hi - lo
    yy <- function(z) bottom + .12 + .55 * (z - lo) / span
    barx <- left + c(.23, .43, .63)
    for (j in seq_along(groups)) {
      rect(barx[j] - .07, yy(0), barx[j] + .07, yy(center[j]), col = cols[j], border = "grey25")
      segments(barx[j], yy(center[j] - err[j]), barx[j], yy(center[j] + err[j]), lwd = .8)
      segments(barx[j] - .035, yy(center[j] - err[j]), barx[j] + .035, yy(center[j] - err[j]))
      segments(barx[j] - .035, yy(center[j] + err[j]), barx[j] + .035, yy(center[j] + err[j]))
      text(barx[j], min(bottom + h - .1, yy(center[j] + err[j]) + .08), cld[k,j], cex = .68)
      text(barx[j], bottom + .045, as.character(groups[j]), cex = .58)
    }
    text(left + .04, bottom + h - .07, v, adj = 0, cex = .68, font = 2)
    axis_y <- pretty(c(lo, hi), n = 2)
    for (tick in axis_y) text(left + .02, yy(tick), format(signif(tick, 2)), adj = 1, cex = .47)
    pval <- tab$p_value[k]
    text(left + w - .03, bottom + h - .07, if (is.finite(pval)) sprintf("p=%.2g", pval) else "p=NA", adj = 1, cex = .5)
  }
  legend("bottom", inset = -.12, xpd = TRUE, legend = paste("Grazing", groups),
         fill = cols[seq_along(groups)], horiz = TRUE, bty = "n", cex = .75)
}

plot_sr_soil_composite <- function(results) {
  plot.new(); plot.window(xlim = c(0, 13), ylim = c(0, 11), xaxs = "i", yaxs = "i")
  title("Mantel links: SR versus soil; soil-indicator Pearson correlations")
  pearson <- results$acceptance$soil_pearson
  if (!is.data.frame(pearson)) pearson <- data.frame(variable_1 = character(), variable_2 = character(),
                                                       r = numeric(), p_value = numeric(), significant = logical())
  soil <- c("BD", "pH", "C", "SOC", "AP", "N", "NO3_N", "NH4_N", "AK", "P")
  soil_tab <- results$acceptance$mantel
  n <- length(soil); base_x <- 2.2; base_y <- 9.5
  palette <- grDevices::colorRampPalette(c("#62B8D1", "white", "#E56B6F"))(101)
  for (i in seq_len(n)) for (j in seq_len(i)) {
    a <- soil[i]; b <- soil[j]
    hit <- pearson[(pearson$variable_1 == a & pearson$variable_2 == b) |
                   (pearson$variable_1 == b & pearson$variable_2 == a), , drop = FALSE]
    x <- base_x + j - 1; y <- base_y - i + 1
    value <- if (nrow(hit)) hit$r[1] else NA_real_
    fill <- if (is.finite(value) && nrow(hit) && isTRUE(hit$significant[1])) palette[round((value + 1) * 50) + 1] else "white"
    rect(x - .43, y - .43, x + .43, y + .43, col = fill, border = "grey65")
    if (is.finite(value) && nrow(hit) && isTRUE(hit$significant[1])) {
      stars <- if (hit$p_value[1] < .001) "***" else if (hit$p_value[1] < .01) "**" else if (hit$p_value[1] < .05) "*" else ""
      text(x, y, sprintf("%.2f%s", value, stars), cex = .58)
    }
  }
  for (i in seq_len(n)) {
    text(base_x + i - 1, base_y + .65, soil[i], srt = 45, adj = 0, xpd = TRUE, cex = .66)
    text(base_x + n - .42, base_y - i + 1, soil[i], adj = 0, xpd = TRUE, cex = .66)
  }
  sr_x <- .65; sr_y <- 5.2
  rect(.08, sr_y - .28, 1.2, sr_y + .28, col = "white", border = "grey25")
  text(.64, sr_y, "Soil respiration (SR)", cex = .62)
  mantel_r <- if (is.data.frame(soil_tab) && nrow(soil_tab)) soil_tab$r[match(soil, soil_tab$indicator)] else rep(NA_real_, n)
  mantel_p <- if (is.data.frame(soil_tab) && nrow(soil_tab)) soil_tab$p_value[match(soil, soil_tab$indicator)] else rep(NA_real_, n)
  for (i in seq_len(n)) {
    y <- base_y - i + 1
    if (is.finite(mantel_r[i])) {
      col <- if (mantel_r[i] >= 0) "#E66A72" else "#4FA5D0"
      segments(sr_x + .55, sr_y, base_x - .48, y, col = col,
               lwd = 1 + 3 * abs(mantel_r[i]), lty = if (mantel_p[i] < .05) 1 else 3)
      points(base_x - .48, y, pch = 21, bg = col, col = "blue", cex = .75)
      text(base_x - .56, y + .18, sprintf("%.2f", mantel_r[i]), adj = 1, cex = .52, col = col)
    }
  }
  legend("right", inset = -.01, xpd = TRUE, legend = c("Mantel positive", "Mantel negative", "p < 0.05", "p >= 0.05"),
         col = c("#E66A72", "#4FA5D0", "grey30", "grey30"), lty = c(1,1,1,3), lwd = c(2,2,1,1),
         bty = "n", cex = .68)
  legend("bottomright", inset = -.04, xpd = TRUE, legend = c("Pearson positive", "Pearson negative", "* p < .05"),
         fill = c("#E56B6F", "#62B8D1", "white"), border = "grey50", bty = "n", cex = .68)
}

plot_diagnostic <- function(results) {
  rec <- results$search$records
  if (is.null(rec) || !nrow(rec)) { plot.new(); title("Search diagnostics"); text(.5, .5, "Search diagnostics not available"); return(invisible()) }
  plot(rec$iteration, rec$score, pch = 19, col = "steelblue", xlab = "Iteration", ylab = "Score")
  if ("pls_proxy_diff" %in% names(rec)) lines(rec$iteration, -rec$pls_proxy_diff, col = "firebrick")
  legend("bottomleft", c("score", "-PLS proxy residual"), col = c("steelblue", "firebrick"), lty = c(NA, 1), pch = c(19, NA), bty = "n")
}

plot_conflict <- function(results) {
  tab <- results$feasibility$conflict_matrix
  if (is.null(tab) || !nrow(tab)) {
    plot.new(); title("Constraint conflict evidence")
    text(.5, .5, "Conflict matrix not available"); return(invisible())
  }
  rows <- unique(tab$row_requirement); cols <- unique(tab$column_requirement)
  status_code <- ifelse(grepl("EMPIRICALLY_FEASIBLE", tab$status), 1,
                        ifelse(grepl("NO_SAMPLE_FOUND", tab$status), -1, 0))
  matrix <- matrix(0, length(rows), length(cols), dimnames=list(rows, cols))
  for (i in seq_len(nrow(tab))) matrix[tab$row_requirement[i], tab$column_requirement[i]] <- status_code[i]
  oldpar <- par(mar=c(8,8,3,2)); on.exit(par(oldpar), add=TRUE)
  # image() creates its own plotting page; do not pre-open one on single-page devices such as EMF.
  image(seq_along(cols), seq_along(rows), t(matrix[nrow(matrix):1,,drop=FALSE]),
        zlim=c(-1,1), col=c("#C85A54", "#D8DADD", "#4C8B72"), axes=FALSE,
        xlab="Constraint", ylab="Constraint")
  axis(1, seq_along(cols), cols, las=2, cex.axis=.75)
  axis(2, seq_along(rows), rev(rows), las=2, cex.axis=.75)
  title("Pairwise empirical constraint matrix")
  for (i in seq_along(rows)) for (j in seq_along(cols))
    text(j, length(rows)-i+1, if (matrix[i,j] > 0) "OK" else if (matrix[i,j] < 0) "?" else "-", cex=.7)
}

`%||%` <- function(x, y) if (is.null(x)) y else x

run_outputs_pass <- function(pass, report_dir, export_manifest, required_outputs,
                             run_id) {
  if (!isTRUE(pass) || !dir.exists(report_dir) || !is.data.frame(export_manifest) ||
      !all(c("path", "role", "exists") %in% names(export_manifest))) return(FALSE)
  root <- normalizePath(report_dir, mustWork = TRUE)
  required <- file.path(root, required_outputs)
  if (!all(file.exists(required))) return(FALSE)
  paths <- normalizePath(export_manifest$path, mustWork = FALSE)
  prefix <- paste0(root, .Platform$file.sep)
  current <- startsWith(paths, prefix) & file.exists(paths) & export_manifest$exists
  final <- export_manifest$role %in% c("FINAL_DATA", "FINAL_RESULTS")
  if (!any(final) || any(final & !current)) return(FALSE)
  manifest_path <- file.path(root, "run_manifest.csv")
  manifest <- tryCatch(utils::read.csv(manifest_path, stringsAsFactors = FALSE),
                       error = function(e) NULL)
  is.data.frame(manifest) && nrow(manifest) == 1L &&
    all(c("run_id", "final_status", "output_role") %in% names(manifest)) &&
    identical(as.character(manifest$run_id[[1L]]), as.character(run_id)) &&
    identical(as.character(manifest$final_status[[1L]]), "PASS") &&
    identical(as.character(manifest$output_role[[1L]]), "FINAL")
}

export_required_plots <- function(results, report_dir) {
  specs <- list(
    fig_植被_SR_Pearson = plot_pearson_heatmap,
    fig_18指标_ANOVA柱状图 = plot_anova_bars,
    fig_SR_土壤_Pearson_Mantel综合图 = plot_sr_soil_composite,
    fig_最终PLS路径图 = plot_pls_paths,
    fig_Mantel连线图 = plot_mantel, fig_Pearson显著性热图 = plot_pearson_heatmap,
    fig_冲突矩阵热图 = plot_conflict, fig_诊断图 = plot_diagnostic,
    fig_路径图 = plot_pls_paths, fig_Pearson热图 = plot_pearson_heatmap)
  paths <- character()
  for (nm in names(specs)) {
    stem <- file.path(report_dir, nm)
    dimensions <- if (identical(nm, "fig_18指标_ANOVA柱状图")) c(12, 10) else
      if (identical(nm, "fig_SR_土壤_Pearson_Mantel综合图")) c(12, 7) else
        if (identical(nm, "fig_最终PLS路径图")) c(12, 8) else c(10, 8)
    write_plot_bundle(stem, function() specs[[nm]](results), width = dimensions[1], height = dimensions[2])
    paths <- c(paths, paste0(stem, ".", c("emf", "png", "pdf")))
  }
  paths
}

export_artifacts <- function(results, project_root, original_data,
                             tables = list(), metadata = list(),
                             export_schema = list()) {
  if (!is.list(results) || is.null(results$run_ok)) stop("results must contain run_ok")
  if (!is.data.frame(original_data) || !ncol(original_data)) stop("original_data must be a non-empty data.frame")
  if (is.null(names(original_data)) || anyNA(names(original_data)) || any(!nzchar(names(original_data))))
    stop("original_data must have non-empty column names")
  root <- normalizePath(project_root, mustWork = FALSE)
  report_dir <- if (file.exists(file.path(root, "R", "00_main.R")))
    file.path(root, "reports") else root
  dir.create(report_dir, recursive = TRUE, showWarnings = FALSE)
  is_final <- isTRUE(results$pass)
  stem <- if (is_final) "final" else "candidate_diagnostic"
  data_out <- file.path(report_dir, paste0(stem, "_data.csv"))
  soil_out <- file.path(report_dir, paste0(stem, "_soil.csv"))
  rds_out <- file.path(report_dir, paste0(stem, "_results.rds"))
  log_out <- file.path(report_dir, "export_manifest.csv")
  sheet_schema <- export_schema$sheet1 %||% original_data
  write_schema_csv(original_data, data_out, sheet_schema)
  if (!is.null(results$candidate$soil)) {
    soil_schema <- export_schema$soil %||% results$candidate$soil
    write_schema_csv(results$candidate$soil, soil_out, soil_schema)
    validate_export_artifacts(results$candidate$soil, soil_out, soil_schema)
  }
  saveRDS(results, rds_out)
  manifest <- data.frame(path = c(data_out, rds_out), kind = c("CSV", "RDS"),
                         role = if (is_final) c("FINAL_DATA", "FINAL_RESULTS") else c("CANDIDATE_DIAGNOSTIC", "CANDIDATE_DIAGNOSTIC"),
                         exists = file.exists(c(data_out, rds_out)), stringsAsFactors = FALSE)
  if (file.exists(soil_out)) manifest <- rbind(manifest, data.frame(path = soil_out, kind = "CSV", role = if (is_final) "FINAL_DATA" else "CANDIDATE_DIAGNOSTIC", exists = TRUE))
  for (nm in names(tables)) {
    tab <- tables[[nm]]
    if (!is.data.frame(tab)) stop(sprintf("table '%s' must be a data.frame", nm))
    p <- file.path(report_dir, paste0(nm, ".csv"))
    utils::write.csv(tab, p, row.names = FALSE, na = "")
    manifest <- rbind(manifest, data.frame(path = p, kind = "CSV", role = "EVIDENCE", exists = file.exists(p)))
  }
  plot_paths <- export_required_plots(results, report_dir)
  manifest <- rbind(manifest, data.frame(path = plot_paths,
                                         kind = toupper(tools::file_ext(plot_paths)),
                                         role = "EVIDENCE",
                                         exists = file.exists(plot_paths), stringsAsFactors = FALSE))
  summary_lines <- c("# Export report", "", "教学模拟数据：仅用于教学演示，不可用于真实科研。",
                     "", paste0("RUN_OK: ", isTRUE(results$run_ok)),
                     paste0("Rows exported: ", nrow(original_data)),
                     paste0("Columns exported: ", paste(names(original_data), collapse = ", ")))
  if (length(metadata)) summary_lines <- c(summary_lines, "", capture.output(dput(metadata)))
  report_path <- file.path(report_dir, "export_report.md")
  writeLines(summary_lines, report_path, useBytes = TRUE)
  manifest <- rbind(manifest, data.frame(path = report_path, kind = "Markdown", role = "EVIDENCE", exists = file.exists(report_path)))
  utils::write.csv(manifest, log_out, row.names = FALSE)
  validate_export_artifacts(original_data, data_out, sheet_schema)
  list(run_ok = TRUE, paths = manifest$path, manifest = manifest)
}

export_schema_frame <- function(data, schema) {
  if (!is.data.frame(data) || !is.data.frame(schema))
    stop("export data and schema must be data.frames")
  if (!identical(names(data), names(schema)))
    stop("exported columns do not match original structure")
  if (nrow(data) != nrow(schema))
    stop("exported row count does not match original structure")
  if (anyNA(data) || any(vapply(data, function(x)
      any(trimws(as.character(x)) == "", na.rm = TRUE), logical(1))))
    stop("export data contain NA or empty values")
  data
}

serialize_export_column <- function(value, prototype) {
  if (inherits(prototype, "POSIXct")) {
    tz <- attr(prototype, "tzone")
    tz <- if (length(tz) && nzchar(tz[[1L]])) tz[[1L]] else "UTC"
    return(format(as.POSIXct(value, origin = "1970-01-01", tz = tz),
                  format = "%Y-%m-%dT%H:%M:%OS6%z", tz = tz))
  }
  if (inherits(prototype, "Date"))
    return(format(as.Date(value, origin = "1970-01-01"), "%Y-%m-%d"))
  if (is.factor(prototype)) return(as.character(value))
  if (is.integer(prototype)) return(format(as.integer(value), scientific = FALSE,
                                            trim = TRUE))
  if (is.numeric(prototype))
    return(format(as.numeric(value), digits = 17L, scientific = FALSE,
                  trim = TRUE, decimal.mark = "."))
  if (is.character(prototype)) return(as.character(value))
  stop("unsupported export column class: ", paste(class(prototype), collapse = "/"))
}

write_schema_csv <- function(data, path, schema) {
  export_schema_frame(data, schema)
  encoded <- as.data.frame(lapply(names(schema), function(nm)
    serialize_export_column(data[[nm]], schema[[nm]])),
    check.names = FALSE, stringsAsFactors = FALSE)
  names(encoded) <- names(schema)
  old_options <- options("scipen")
  on.exit(options(old_options), add = TRUE)
  options(scipen = 999)
  utils::write.csv(encoded, path, row.names = FALSE, na = "", quote = TRUE)
  invisible(path)
}

read_schema_csv <- function(csv_path, schema) {
  encoded <- utils::read.csv(csv_path, check.names = FALSE,
                             stringsAsFactors = FALSE, colClasses = "character",
                             na.strings = character())
  if (!identical(names(encoded), names(schema)))
    stop("exported columns do not match original structure")
  decoded <- lapply(names(schema), function(nm) {
    value <- encoded[[nm]]
    prototype <- schema[[nm]]
    if (inherits(prototype, "POSIXct")) {
      tz <- attr(prototype, "tzone")
      tz <- if (length(tz) && nzchar(tz[[1L]])) tz[[1L]] else "UTC"
      return(as.POSIXct(value, format = "%Y-%m-%dT%H:%M:%OS%z", tz = tz))
    }
    if (inherits(prototype, "Date")) return(as.Date(value, format = "%Y-%m-%d"))
    if (is.factor(prototype))
      return(factor(value, levels = levels(prototype), ordered = is.ordered(prototype)))
    if (is.integer(prototype)) return(as.integer(value))
    if (is.numeric(prototype)) return(as.numeric(value))
    if (is.character(prototype)) return(as.character(value))
    stop("unsupported export schema class: ", paste(class(prototype), collapse = "/"))
  })
  names(decoded) <- names(schema)
  as.data.frame(decoded, check.names = FALSE, stringsAsFactors = FALSE)
}

validate_export_artifacts <- function(original_data, csv_path,
                                      schema = original_data) {
  if (!file.exists(csv_path)) stop("required CSV export is missing")
  export_schema_frame(original_data, schema)
  exported <- read_schema_csv(csv_path, schema)
  if (nrow(exported) != nrow(original_data)) stop("exported row count does not match original data")
  if (!identical(names(exported), names(schema)))
    stop("exported columns do not match original structure")
  class_mismatch <- names(schema)[!vapply(names(schema), function(nm)
    identical(class(exported[[nm]]), class(schema[[nm]])), logical(1))]
  if (length(class_mismatch))
    stop("exported column types do not match original structure: ",
         paste(class_mismatch, collapse = ", "))
  values_match <- vapply(names(schema), function(nm) {
    x <- original_data[[nm]]
    y <- exported[[nm]]
    if (is.numeric(x)) identical(as.numeric(x), as.numeric(y)) else
      identical(x, y)
  }, logical(1))
  if (!all(values_match))
    stop("exported values changed during CSV round trip: ",
         paste(names(values_match)[!values_match], collapse = ", "))
  if (anyNA(exported) || any(vapply(exported, function(x)
      any(trimws(as.character(x)) == "", na.rm = TRUE), logical(1))))
    stop("exported data contain NA or empty values")
  TRUE
}
