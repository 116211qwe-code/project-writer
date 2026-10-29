# Module 7: statistical acceptance for ANOVA/LSD, Pearson, and Mantel.
# Teaching simulation only; not for real scientific inference.

statistical_acceptance_columns <- function() {
  list(
    anova = c("AGB", "BGB", "coverage", "LB", "H", "E", "D", "R",
              "BD", "pH", "C", "SOC", "AP", "N", "NO3_N", "NH4_N", "AK", "P"),
    pearson = c("AGB", "BGB", "coverage", "LB", "H", "E", "D", "R", "SR"),
    soil = c("BD", "pH", "C", "SOC", "AP", "N", "NO3_N", "NH4_N", "AK", "P")
  )
}

validate_statistical_data <- function(data, columns, groups = character(), min_rows = 3L) {
  if (!is.data.frame(data)) stop("data must be a data.frame", call. = FALSE)
  required <- unique(c(groups, columns))
  missing <- setdiff(required, names(data))
  if (length(missing)) stop("missing statistical columns: ", paste(missing, collapse = ", "), call. = FALSE)
  if (nrow(data) < min_rows) stop("insufficient rows for statistical acceptance", call. = FALSE)
  if (anyNA(data[required])) stop("statistical columns contain missing values", call. = FALSE)
  bad_numeric <- columns[!vapply(data[columns], is.numeric, logical(1))]
  if (length(bad_numeric)) stop("statistical indicators must be numeric: ", paste(bad_numeric, collapse = ", "), call. = FALSE)
  non_finite <- columns[vapply(data[columns], function(x) any(!is.finite(x)), logical(1))]
  if (length(non_finite)) stop("statistical indicators must be finite: ", paste(non_finite, collapse = ", "), call. = FALSE)
  invisible(data)
}

statistical_check <- function(requirement, target, actual, display, margin, repair_suggestion) {
  contract <- getOption("acceptance_contract", NULL)
  if (is.null(contract)) stop("acceptance contract is required for statistical evaluation")
  contract_rule(contract, requirement)
  if (!is.numeric(actual) || length(actual) != 1L)
    stop("statistical actual must be one typed numeric scalar: ", requirement)
  status <- if (is.finite(actual)) "EVALUATED" else "NOT_EVALUATED"
  passed <- is.finite(actual) && contract_pass(contract, requirement, actual)
  data.frame(requirement = requirement, target = target, actual = actual,
             actual_display = display, pass = passed, status = status,
             critical_margin = margin, repair_suggestion = repair_suggestion,
             stringsAsFactors = FALSE)
}

statistical_object_signature <- function(data, source_object, key_columns) {
  if (!requireNamespace("digest", quietly = TRUE))
    stop("package 'digest' is required for statistical provenance", call. = FALSE)
  list(source_object = source_object, nrow = nrow(data), ncol = ncol(data),
       key_columns = paste(key_columns, collapse = "+"),
       derivation_hash = digest::digest(
         list(source_object = source_object, key_columns = key_columns, data = data),
         algo = "sha256", serialize = TRUE))
}

stamp_statistical_table <- function(table, signature) {
  table$source_object <- signature$source_object
  table$nrow <- signature$nrow
  table$ncol <- signature$ncol
  table$key_columns <- signature$key_columns
  table$derivation_hash <- signature$derivation_hash
  attr(table, "signature") <- signature
  table
}

assert_july_45 <- function(data, columns, label) {
  validate_statistical_data(data, columns,
                            groups = c("month", "plot", "rep", "grazing"), min_rows = 45L)
  if (nrow(data) != 45L)
    stop(label, " requires exactly 45 July rows; got ", nrow(data), call. = FALSE)
  if (!all(as.character(data$month) == "7"))
    stop(label, " requires July-only data", call. = FALSE)
  invisible(data)
}

accept_anova <- function(july45) {
  spec <- statistical_acceptance_columns()
  assert_july_45(july45, spec$anova, "accept_anova(july45)")
  signature <- statistical_object_signature(
    july45, "july45", c("month", "plot", "rep", "grazing"))
  fit <- run_july_anova_lsd(july45, spec$anova)
  list(anova = stamp_statistical_table(fit$anova, signature),
       lsd = stamp_statistical_table(fit$lsd, signature),
       july = july45, signature = signature)
}

accept_pearson <- function(agg75) {
  columns <- statistical_acceptance_columns()$pearson
  validate_statistical_data(agg75, columns,
                            groups = c("month", "plot", "grazing"), min_rows = 1L)
  if (nrow(agg75) != 75L)
    stop("accept_pearson(agg75) requires exactly 75 aggregated rows; got ",
         nrow(agg75), call. = FALSE)
  if (anyDuplicated(agg75[c("month", "plot")]))
    stop("accept_pearson(agg75) requires unique month+plot keys", call. = FALSE)
  signature <- statistical_object_signature(
    agg75, "agg75", c("month", "plot", "grazing"))
  pairs <- utils::combn(columns, 2L, simplify = FALSE)
  pearson_rows <- lapply(pairs, function(pair) {
    test <- stats::cor.test(agg75[[pair[[1L]]]], agg75[[pair[[2L]]]],
                            method = "pearson")
    data.frame(variable_1 = pair[[1L]], variable_2 = pair[[2L]],
               r = unname(test$estimate), p_value = test$p.value,
               significant = test$p.value < .05, stringsAsFactors = FALSE)
  })
  stamp_statistical_table(do.call(rbind, pearson_rows), signature)
}

accept_mantel <- function(sheet1_full, soil_july, seed = 1L,
                          permutations = 100L) {
  columns <- statistical_acceptance_columns()$soil
  key_columns <- c("grazing", "plot", "rep")
  validate_statistical_data(sheet1_full, "SR",
                            groups = c(key_columns, "month"), min_rows = 1L)
  validate_statistical_data(soil_july, columns,
                            groups = c(key_columns, "month"), min_rows = 1L)
  if (length(permutations) != 1L || !is.finite(permutations) ||
      permutations < 1L || permutations != as.integer(permutations))
    stop("permutations must be a positive integer", call. = FALSE)
  if (as.integer(permutations) != 100L)
    stop("accept_mantel() requires the specified 100 permutations",
         call. = FALSE)
  permutations <- as.integer(permutations)
  if (nrow(sheet1_full) != 225L)
    stop("accept_mantel() requires all 225 Sheet1 rows; got ",
         nrow(sheet1_full),
         call. = FALSE)
  if (!setequal(as.character(sheet1_full$month), as.character(6:10)))
    stop("accept_mantel() requires months 6, 7, 8, 9, and 10",
         call. = FALSE)
  if (nrow(soil_july) != 45L ||
      !all(as.character(soil_july$month) == "7"))
    stop("accept_mantel() requires exactly 45 July soil rows",
         call. = FALSE)
  make_key <- function(data, columns) {
    do.call(paste, c(lapply(data[columns], as.character), sep = "\r"))
  }
  full_key <- make_key(sheet1_full, key_columns)
  full_month <- as.character(sheet1_full$month)
  key_order <- unique(full_key)
  if (length(key_order) != 45L || anyDuplicated(data.frame(
      key = full_key, month = full_month)))
    stop("Sheet1 must contain one row per plot/rep key and month",
         call. = FALSE)
  month_matrix <- vapply(as.character(6:10), function(month) {
    rows <- which(full_month == month)
    positions <- match(key_order, full_key[rows])
    if (length(rows) != 45L || anyNA(positions))
      stop("each Sheet1 month must contain the same 45 plot/rep keys",
           call. = FALSE)
    sheet1_full$SR[rows[positions]]
  }, numeric(45L))
  colnames(month_matrix) <- as.character(6:10)
  rownames(month_matrix) <- key_order
  if (length(unique(as.vector(month_matrix))) < 2L)
    stop("Mantel requires variation in the five-month SR matrix",
         call. = FALSE)

  soil_key <- make_key(soil_july, key_columns)
  if (anyDuplicated(soil_key) || !setequal(key_order, soil_key))
    stop("July soil keys must exactly match the 45 Sheet1 plot/rep keys",
         call. = FALSE)
  soil_match <- match(key_order, soil_key)
  if (any(vapply(soil_july[columns], function(x) stats::sd(x) == 0,
                 logical(1))))
    stop("Mantel requires variation in every July soil indicator",
         call. = FALSE)

  signature <- list(
    source_object = "sheet1_full+soil_july", nrow = 45L, ncol = 15L,
    key_columns = paste(key_columns, collapse = "+"),
    derivation_hash = digest::digest(
      list(month_matrix = month_matrix,
           soil_july = soil_july[soil_match, c(key_columns, "month", columns),
                                 drop = FALSE]),
      algo = "sha256", serialize = TRUE))

  had_seed <- exists(".Random.seed", .GlobalEnv, inherits = FALSE)
  old_seed <- if (had_seed) get(".Random.seed", .GlobalEnv) else NULL
  on.exit(if (had_seed) assign(".Random.seed", old_seed, .GlobalEnv) else
    if (exists(".Random.seed", .GlobalEnv, inherits = FALSE))
      rm(list = ".Random.seed", envir = .GlobalEnv), add = TRUE)

  sr_distance <- as.matrix(stats::dist(month_matrix, method = "euclidean"))
  lower <- lower.tri(sr_distance)
  mantel_rows <- lapply(seq_along(columns), function(i) {
    if (!is.null(seed)) set.seed(seed + i - 1L)
    soil_distance <- as.matrix(stats::dist(
      soil_july[[columns[[i]]]][soil_match], method = "euclidean"))
    observed <- stats::cor(sr_distance[lower], soil_distance[lower])
    null <- numeric(permutations)
    for (b in seq_len(permutations)) {
      permutation <- sample.int(nrow(soil_distance))
      permuted <- soil_distance[permutation, permutation, drop = FALSE]
      null[[b]] <- stats::cor(sr_distance[lower], permuted[lower])
    }
    p_value <- (1 + sum(abs(null) >= abs(observed))) / (permutations + 1L)
    data.frame(indicator = columns[[i]], r = observed, p_value = p_value,
               significant = p_value < .05, stringsAsFactors = FALSE)
  })
  stamp_statistical_table(do.call(rbind, mantel_rows), signature)
}

accept_soil_internal <- function(july) {
  columns <- statistical_acceptance_columns()$soil
  validate_statistical_data(july, columns,
                            groups = c("month", "grazing"), min_rows = 1L)
  if (nrow(july) != 45L)
    stop("accept_soil_internal(july) requires exactly 45 July rows; got ",
         nrow(july), call. = FALSE)
  if (!all(as.character(july$month) == "7"))
    stop("accept_soil_internal(july) requires July-only data", call. = FALSE)
  key_columns <- intersect(c("month", "plot", "rep", "grazing"), names(july))
  signature <- statistical_object_signature(
    july, "july", key_columns)
  pairs <- utils::combn(columns, 2L, simplify = FALSE)
  pearson_rows <- lapply(pairs, function(pair) {
    test <- stats::cor.test(july[[pair[[1L]]]], july[[pair[[2L]]]],
                            method = "pearson")
    data.frame(variable_1 = pair[[1L]], variable_2 = pair[[2L]],
               r = unname(test$estimate), p_value = test$p.value,
               significant = test$p.value < .05, stringsAsFactors = FALSE)
  })
  stamp_statistical_table(do.call(rbind, pearson_rows), signature)
}

attach_acceptance_provenance <- function(checks, signatures) {
  source_for <- function(requirement) {
    if (requirement %in% c("R-060", "R-061", "R-062")) return(signatures$anova)
    if (requirement %in% sprintf("R-%03d", 64:69)) return(signatures$pearson)
    if (requirement %in% sprintf("R-%03d", 70:73)) return(signatures$mantel)
    if (requirement %in% c("R-074", "R-075")) return(signatures$soil_internal)
    if (requirement %in% c("R-027", "R-028")) {
      combined <- c(signatures$pearson, signatures$mantel,
                    signatures$soil_internal)
      combined$source_object <- "agg75+july"
      combined$nrow <- "75+45"
      combined$ncol <- paste(signatures$pearson$ncol,
                             signatures$mantel$ncol, sep = "+")
      combined$key_columns <- paste(signatures$pearson$key_columns,
                                    signatures$mantel$key_columns, sep = ";")
      combined$derivation_hash <- digest::digest(
        lapply(signatures[c("pearson", "mantel", "soil_internal")],
               `[[`, "derivation_hash"), algo = "sha256")
      return(combined)
    }
    list(source_object = "combined_acceptance",
         nrow = "45+75+45", ncol = paste(signatures$anova$ncol,
                                          signatures$pearson$ncol,
                                          signatures$mantel$ncol, sep = "+"),
         key_columns = paste(signatures$anova$key_columns,
                             signatures$pearson$key_columns,
                             signatures$mantel$key_columns, sep = ";"),
         derivation_hash = digest::digest(
           lapply(signatures, `[[`, "derivation_hash"), algo = "sha256"))
  }
  provenance <- lapply(checks$requirement, source_for)
  for (field in c("source_object", "nrow", "ncol", "key_columns", "derivation_hash"))
    checks[[field]] <- vapply(provenance, function(x) as.character(x[[field]]),
                              character(1))
  checks
}

evaluate_ecology_signs <- function(data) {
  required <- c("grazing", "SR", "AGB", "BGB", "BD", "C", "SOC", "AP", "N", "NO3_N", "NH4_N", "AK", "P", "E", "D", "H")
  # Merges of the Sheet1 and soil tables may add suffixes to overlapping
  # columns. Resolve a required semantic column only when the source is
  # unambiguous; never silently choose between conflicting duplicates.
  resolve_column <- function(nm) {
    if (nm %in% names(data)) return(nm)
    candidates <- grep(paste0("^", nm, "\\.(soil|x|y)$"), names(data), value = TRUE)
    if (length(candidates) == 1L) return(candidates[[1L]])
    if (length(candidates) > 1L) stop("ambiguous ecology sign columns for ", nm, ": ", paste(candidates, collapse=","), call. = FALSE)
    NA_character_
  }
  resolved <- setNames(vapply(required, resolve_column, character(1)), required)
  missing <- names(resolved)[is.na(resolved)]
  if (length(missing)) stop("missing ecology sign columns: ", paste(missing, collapse=","), call. = FALSE)
  pairs <- data.frame(
    object=c("Grazing->SR","Grazing->AGB","Grazing->BGB","Grazing->BD",
             "SR->E","E->D","AGB->BGB","SR->C","SR->SOC","SR->N","SR->P","SR->BD",
             "AGB->C","AGB->SOC","AGB->N","AGB->P","NH4_N->NO3_N","NH4_N->SR","NO3_N->SR"),
    x=c("grazing","grazing","grazing","grazing","SR","E","AGB","SR","SR","SR","SR","SR","AGB","AGB","AGB","AGB","NH4_N","NH4_N","NO3_N"),
    y=c("SR","AGB","BGB","BD","E","D","BGB","C","SOC","N","P","BD","C","SOC","N","P","NO3_N","SR","SR"),
    expected=c(-1,-1,-1,1,1,-1,1,1,1,1,1,-1,1,1,1,1,-1,1,-1), stringsAsFactors=FALSE)
  pairs$r <- mapply(function(x,y) suppressWarnings(stats::cor(data[[resolved[[x]]]], data[[resolved[[y]]]])), pairs$x, pairs$y)
  pairs$pass <- sign(pairs$r) == pairs$expected
  pairs$method <- "Pearson correlation; full input or July subset"
  pairs
}

run_july_anova_lsd <- function(data, columns) {
  july <- data[as.character(data$month) == "7", , drop = FALSE]
  if (!nrow(july)) stop("ANOVA requires July observations", call. = FALSE)
  if (length(unique(as.character(july$grazing))) < 2L || any(table(july$grazing) < 2L))
    stop("July ANOVA requires at least two grazing groups with two observations each", call. = FALSE)
  anova_rows <- vector("list", length(columns))
  lsd_rows <- list()
  for (i in seq_along(columns)) {
    variable <- columns[[i]]
    fit <- stats::aov(july[[variable]] ~ factor(july$grazing))
    tab <- summary(fit)[[1L]]
    p <- as.numeric(tab[1L, "Pr(>F)"])
    anova_rows[[i]] <- data.frame(indicator = variable, F = as.numeric(tab[1L, "F value"]),
                                  p_value = p, significant = is.finite(p) && p < .05,
                                  stringsAsFactors = FALSE)
    means <- tapply(july[[variable]], as.character(july$grazing), mean)
    groups <- names(means)
    pairs <- utils::combn(groups, 2L, simplify = FALSE)
    mse <- as.numeric(tail(tab[, "Mean Sq"], 1L))
    df_error <- as.numeric(tail(tab[, "Df"], 1L))
    lsd_rows[[i]] <- do.call(rbind, lapply(pairs, function(pair) {
      ni <- sum(as.character(july$grazing) == pair[[1L]])
      nj <- sum(as.character(july$grazing) == pair[[2L]])
      se <- sqrt(mse * (1 / ni + 1 / nj))
      diff <- means[[pair[[1L]]]] - means[[pair[[2L]]]]
      t_value <- diff / se
      p_lsd <- 2 * stats::pt(-abs(t_value), df = df_error)
      data.frame(indicator = variable, group_1 = pair[[1L]], group_2 = pair[[2L]],
                 mean_difference = diff, p_value = p_lsd,
                 significant = is.finite(p_lsd) && p_lsd < .05, stringsAsFactors = FALSE)
    }))
  }
  list(anova = do.call(rbind, anova_rows), lsd = do.call(rbind, lsd_rows), july = july)
}

aggregate_pearson_data <- function(data, columns) {
  key <- interaction(data$month, data$plot, drop = TRUE, lex.order = TRUE)
  groups <- split(seq_len(nrow(data)), key)
  rows <- lapply(groups, function(idx) {
    if (length(unique(as.character(data$grazing[idx]))) != 1L)
      stop("month+plot aggregation has conflicting grazing labels", call. = FALSE)
    out <- data.frame(month = data$month[idx[[1L]]], plot = data$plot[idx[[1L]]],
                      grazing = data$grazing[idx[[1L]]], stringsAsFactors = FALSE)
    for (column in columns) out[[column]] <- mean(data[[column]][idx])
    out
  })
  aggregated <- do.call(rbind, rows)
  rownames(aggregated) <- NULL
  if (any(vapply(aggregated[columns], stats::sd, numeric(1)) == 0))
    stop("Pearson aggregation contains a zero-variance indicator", call. = FALSE)
  aggregated
}

correlation_table <- function(data, columns) {
  pairs <- utils::combn(columns, 2L, simplify = FALSE)
  do.call(rbind, lapply(pairs, function(pair) {
    test <- stats::cor.test(data[[pair[[1L]]]], data[[pair[[2L]]]], method = "pearson")
    data.frame(variable_1 = pair[[1L]], variable_2 = pair[[2L]], r = unname(test$estimate),
               p_value = test$p.value, significant = test$p.value < .05,
               stringsAsFactors = FALSE)
  }))
}

mantel_permutation <- function(x, y, permutations, seed = NULL) {
  had_seed <- exists(".Random.seed", .GlobalEnv, inherits=FALSE)
  old_seed <- if (had_seed) get(".Random.seed", .GlobalEnv) else NULL
  on.exit(if (had_seed) assign(".Random.seed", old_seed, .GlobalEnv) else if (exists(".Random.seed", .GlobalEnv, inherits=FALSE)) rm(list=".Random.seed", envir=.GlobalEnv), add=TRUE)
  if (!is.null(seed)) set.seed(seed)
  dx <- as.matrix(stats::dist(x, method = "euclidean"))
  dy <- as.matrix(stats::dist(y, method = "euclidean"))
  lower <- lower.tri(dx)
  observed <- stats::cor(dx[lower], dy[lower])
  null <- numeric(permutations)
  for (b in seq_len(permutations)) {
    idx <- sample.int(nrow(dy))
    permuted <- dy[idx, idx, drop = FALSE]
    null[[b]] <- stats::cor(dx[lower], permuted[lower])
  }
  p <- (1 + sum(abs(null) >= abs(observed))) / (permutations + 1)
  c(r = observed, p_value = p)
}

run_july_mantel <- function(july, soil_columns, permutations, seed = NULL) {
  if (length(unique(as.character(july$SR))) < 2L)
    stop("Mantel requires variation in July SR", call. = FALSE)
  if (any(vapply(july[soil_columns], function(x) stats::sd(x) == 0, logical(1))))
    stop("Mantel requires variation in every July soil indicator", call. = FALSE)
  do.call(rbind, lapply(seq_along(soil_columns), function(i) {
    column <- soil_columns[[i]]
    result <- mantel_permutation(matrix(july$SR, ncol = 1L),
                                 matrix(july[[column]], ncol = 1L),
                                 permutations = permutations,
                                 seed = if (is.null(seed)) NULL else seed + i - 1L)
    data.frame(indicator = column, r = unname(result[["r"]]),
               p_value = unname(result[["p_value"]]), significant = result[["p_value"]] < .05,
               stringsAsFactors = FALSE)
  }))
}

evaluate_statistical_acceptance <- function(july_data, pearson_source,
                                            permutations = 100L, seed = 1L,
                                            ecology_source = pearson_source) {
  warning("evaluate_statistical_acceptance() is deprecated; use the four accept_*() entries",
          call. = FALSE, immediate. = TRUE)
  spec <- statistical_acceptance_columns()
  validate_statistical_data(july_data, unique(c(spec$anova, spec$soil)),
                            groups = c("month", "plot", "grazing"), min_rows = 45L)
  if (nrow(july_data) != 45L || !all(as.character(july_data$month) == "7"))
    stop("deprecated wrapper requires a 45-row July object", call. = FALSE)
  if (length(permutations) != 1L || !is.finite(permutations) || permutations < 1L ||
      permutations != as.integer(permutations))
    stop("permutations must be a positive integer", call. = FALSE)
  if (as.integer(permutations) != 100L)
    stop("the July Mantel acceptance entry requires the specified 100 permutations",
         call. = FALSE)
  july_data$month <- as.character(july_data$month)
  pearson_source$month <- as.character(pearson_source$month)
  anova <- accept_anova(july_data)
  pearson_data <- if (nrow(pearson_source) == 75L) pearson_source else
    aggregate_pearson_data(pearson_source, spec$pearson)
  pearson <- accept_pearson(pearson_data)
  soil_july <- anova$july
  mantel <- accept_mantel(pearson_source, soil_july, seed = seed)
  soil_pearson <- accept_soil_internal(soil_july)
  assemble_statistical_acceptance(anova, pearson_data, pearson, mantel,
                                  soil_pearson, ecology_source,
                                  as.integer(permutations), seed)
}

assemble_statistical_acceptance <- function(anova, pearson_data, pearson,
                                            mantel, soil_pearson,
                                            ecology_source, permutations,
                                            seed) {
  spec <- statistical_acceptance_columns()
  soil_july <- anova$july
  ecology_full <- evaluate_ecology_signs(ecology_source)
  ecology_july <- evaluate_ecology_signs(soil_july)
  sr_soil_pearson <- do.call(rbind, lapply(spec$soil, function(column) {
    test <- stats::cor.test(soil_july$SR, soil_july[[column]], method = "pearson")
    data.frame(variable_1 = "SR", variable_2 = column, r = unname(test$estimate),
               p_value = test$p.value, significant = test$p.value < .05,
               stringsAsFactors = FALSE)
  }))
  sr_soil_pearson <- stamp_statistical_table(
    sr_soil_pearson, attr(soil_pearson, "signature"))

  sr_cor <- pearson[pearson$variable_2 == "SR", , drop = FALSE]
  significant_sr <- sum(sr_cor$significant)
  required_sr <- c("AGB", "LB", "D")
  required_sr_rows <- match(required_sr, sr_cor$variable_1)
  required_sr_pass <- !anyNA(required_sr_rows) && all(sr_cor$significant[required_sr_rows])
  agb_h <- pearson[pearson$variable_1 == "AGB" & pearson$variable_2 == "H", , drop = FALSE]
  agb_lb <- pearson[pearson$variable_1 == "AGB" & pearson$variable_2 == "LB", , drop = FALSE]
  eligible_cap <- !(pearson$variable_1 %in% c("H", "E", "D", "R") |
                    pearson$variable_2 %in% c("H", "E", "D", "R"))
  capped <- abs(pearson$r[eligible_cap])
  mantel_bd_ph <- mantel[mantel$indicator %in% c("BD", "pH"), , drop = FALSE]
  nh4 <- mantel[mantel$indicator == "NH4_N", , drop = FALSE]
  no3 <- mantel[mantel$indicator == "NO3_N", , drop = FALSE]
  bd_pearson <- sr_soil_pearson[sr_soil_pearson$variable_2 == "BD", , drop=FALSE]
  ph_pearson <- sr_soil_pearson[sr_soil_pearson$variable_2 == "pH", , drop=FALSE]
  nh4_pearson <- sr_soil_pearson[sr_soil_pearson$variable_2 == "NH4_N", , drop=FALSE]
  no3_pearson <- sr_soil_pearson[sr_soil_pearson$variable_2 == "NO3_N", , drop=FALSE]
  nh4_no3 <- soil_pearson[((soil_pearson$variable_1 == "NH4_N" & soil_pearson$variable_2 == "NO3_N") |
                            (soil_pearson$variable_1 == "NO3_N" & soil_pearson$variable_2 == "NH4_N")), , drop=FALSE]
  soil_sig_count <- sum(soil_pearson$significant)
  soil_high_count <- sum(abs(soil_pearson$r) > .7)
  soil_over_08 <- any(abs(soil_pearson$r) > .8)
  sr_monthly <- do.call(rbind, lapply(6:9, function(m) {
    d <- ecology_source[as.character(ecology_source$month) == as.character(m), , drop = FALSE]
    fit <- tryCatch(stats::aov(SR ~ factor(grazing), data = d), error = function(e) NULL)
    p <- if (is.null(fit)) NA_real_ else as.numeric(summary(fit)[[1L]][1L, "Pr(>F)"])
    data.frame(month = m, p_value = p, significant = is.finite(p) && p < .05, n = nrow(d))
  }))
  grazing3_highest <- sum(vapply(setdiff(spec$anova, c("BD", "pH")), function(v) {
    means <- tapply(soil_july[[v]], as.character(soil_july$grazing), mean)
    "3" %in% names(means) && which.max(means) == match("3", names(means))
  }, logical(1)))
  checks <- rbind(
    statistical_check("R-060", "18 indicators and LSD comparisons", as.numeric(nrow(anova$anova) == 18L && nrow(anova$lsd) > 0L), sprintf("ANOVA=%d; LSD=%d", nrow(anova$anova), nrow(anova$lsd)), NA_real_, "Check July ANOVA and LSD output"),
    statistical_check("R-064", "75 aggregated rows", as.numeric(nrow(pearson_data)), sprintf("%d rows", nrow(pearson_data)), nrow(pearson_data) - 75L, "Check month+plot aggregation"),
    statistical_check("R-061", "At least 12 significant ANOVAs", sum(anova$anova$significant, na.rm=TRUE), sprintf("%d / 18", sum(anova$anova$significant, na.rm=TRUE)), NA_real_, "Review July group differences"),
    statistical_check("R-062", "At most 2 grazing=3 maxima", grazing3_highest, as.character(grazing3_highest), NA_real_, "Review grazing means"),
    statistical_check("R-063", "SR significant in months 6-9", sum(sr_monthly$significant, na.rm=TRUE), paste(sr_monthly$month[sr_monthly$significant], collapse=","), NA_real_, "Review monthly SR ANOVA"),
    statistical_check("R-065", "At least 5 significant SR links", significant_sr, as.character(significant_sr), NA_real_, "Review aggregated Pearson results"),
    statistical_check("R-066", "AGB, LB, D significant with SR", sum(sr_cor$significant[required_sr_rows], na.rm=TRUE), paste(required_sr, as.integer(sr_cor$significant[required_sr_rows]), collapse=";"), NA_real_, "Review required SR correlations"),
    statistical_check("R-067", "AGB-H and AGB-LB both positive", as.numeric(nrow(agb_h)==1L && nrow(agb_lb)==1L && agb_h$r > 0 && agb_lb$r > 0) * 2, sprintf("AGB-H=%g; AGB-LB=%g", agb_h$r, agb_lb$r), NA_real_, "Review required correlation signs"),
    statistical_check("R-068-lower", "At least 10 significant Pearson pairs", sum(pearson$significant, na.rm=TRUE), as.character(sum(pearson$significant, na.rm=TRUE)), NA_real_, "Review Pearson network density"),
    statistical_check("R-068-upper", "At most 25 significant Pearson pairs", sum(pearson$significant, na.rm=TRUE), as.character(sum(pearson$significant, na.rm=TRUE)), NA_real_, "Review Pearson network density"),
    statistical_check("R-069", "Non-exempt |r| <= 0.7", if(length(capped)) max(capped) else 0, if(length(capped)) sprintf("max=%g",max(capped)) else "no pairs", NA_real_, "Review correlation cap"),
    statistical_check("R-070", "10 Mantel indicators, 100 permutations", as.numeric(nrow(mantel)==10L && nrow(soil_july)==45L && permutations==100L), sprintf("n=%d; permutations=%d",nrow(mantel),permutations), NA_real_, "Review Mantel inputs"),
    statistical_check("R-071", "At least 7 significant Mantel relations", sum(mantel$significant, na.rm=TRUE), as.character(sum(mantel$significant, na.rm=TRUE)), NA_real_, "Review Mantel results"),
    statistical_check("R-072", "No positive-significant Mantel BD/pH", sum(mantel_bd_ph$r > 0 & mantel_bd_ph$significant, na.rm=TRUE), as.character(sum(mantel_bd_ph$r > 0 & mantel_bd_ph$significant, na.rm=TRUE)), NA_real_, "Review BD/pH Mantel direction"),
    statistical_check("R-073", "NH4_N and NO3_N Mantel signs oppose", as.numeric(nrow(nh4)==1L && nrow(no3)==1L && sign(nh4$r)!=sign(no3$r)), sprintf("NH4=%g; NO3=%g",nh4$r,no3$r), NA_real_, "Review nitrogen Mantel signs"),
    statistical_check("R-027", "No positive-significant SR-BD/pH Pearson", sum(c(bd_pearson$r,ph_pearson$r)>0 & c(bd_pearson$significant,ph_pearson$significant), na.rm=TRUE), sprintf("BD=%g; pH=%g",bd_pearson$r,ph_pearson$r), NA_real_, "Review SR-soil directions"),
    statistical_check("R-028", "Nitrogen SR signs oppose and internal pair negative", as.numeric(nrow(nh4_pearson)==1L && nrow(no3_pearson)==1L && nrow(nh4_no3)==1L && sign(nh4_pearson$r)!=sign(no3_pearson$r) && nh4_no3$r<0), sprintf("NH4-SR=%g; NO3-SR=%g; pair=%g",nh4_pearson$r,no3_pearson$r,nh4_no3$r), NA_real_, "Review nitrogen Pearson signs"),
    statistical_check("R-074-lower", "At least 10 significant soil pairs", soil_sig_count, as.character(soil_sig_count), NA_real_, "Review soil Pearson network"),
    statistical_check("R-074-upper", "At most 25 significant soil pairs", soil_sig_count, as.character(soil_sig_count), NA_real_, "Review soil Pearson network"),
    statistical_check("R-075-count", "At most 4 soil pairs with |r| > 0.7", soil_high_count, as.character(soil_high_count), NA_real_, "Review soil correlation magnitudes"),
    statistical_check("R-075-max", "Every soil pair has |r| < 0.8", if(nrow(soil_pearson)) max(abs(soil_pearson$r),na.rm=TRUE) else NA_real_, as.character(if(nrow(soil_pearson)) max(abs(soil_pearson$r),na.rm=TRUE) else NA_real_), NA_real_, "Review soil correlation magnitudes")
  )
  check_shape <- all(c("requirement", "target", "actual", "actual_display", "pass", "critical_margin", "repair_suggestion") %in% names(checks))
  checks <- rbind(checks, statistical_check("R-076", "Required acceptance columns", as.numeric(check_shape), as.character(check_shape), NA_real_, "Restore acceptance schema"))
  checks <- rbind(checks, statistical_check("C6_statistics", "All statistical hard rules pass", as.numeric(all(checks$pass)), as.character(all(checks$pass)), NA_real_, "Resolve every failed statistical contract rule"))
  checks <- attach_acceptance_provenance(checks, list(
    anova = anova$signature, pearson = attr(pearson, "signature"),
    mantel = attr(mantel, "signature"),
    soil_internal = attr(soil_pearson, "signature")))
  pending <- sr_monthly
  list(checks = checks, ecology_full=ecology_full, ecology_july=ecology_july,
       anova = anova$anova, lsd = anova$lsd,
       pearson_data = pearson_data, pearson = pearson,
       mantel = mantel, soil_pearson = soil_pearson,
       sr_soil_pearson = sr_soil_pearson,
       pending = pending,
       settings = list(july_month = "7", pearson_rows = nrow(pearson_data),
                       mantel_permutations = as.integer(permutations), seed = seed),
       has_pass = all(checks$pass))
}
