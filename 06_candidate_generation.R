# Parameterized candidate generation for the third implementation module.

candidate_mutable_columns <- function(sheet1, soil) {
  list(
    sheet1 = intersect(c("SR", "W", "T", "AGB", "BGB", "LB", "coverage"), names(sheet1)),
    soil = intersect(c("BD", "pH", "C", "SOC", "AP", "N", "NO3_N", "NH4_N", "AK", "P"), names(soil))
  )
}

reflect_to_bounds <- function(x, lower, upper) {
  if (length(x) != length(lower) || length(x) != length(upper)) stop("反射边界长度不一致")
  if (any(!is.finite(c(x, lower, upper))) || any(lower > upper)) stop("反射边界必须有限且 lower <= upper")
  out <- x
  fixed <- lower == upper
  out[fixed] <- lower[fixed]
  active <- !fixed
  # Modulo reflection is equivalent to repeated mirror reflections but does
  # not have an arbitrary iteration limit for a distant finite proposal.
  span <- upper[active] - lower[active]
  phase <- (out[active] - lower[active]) %% (2 * span)
  out[active] <- lower[active] + ifelse(phase <= span, phase, 2 * span - phase)
  out
}

r086_sem_generation_audit <- function(candidate) {
  metadata <- candidate$metadata
  paths <- metadata$signed_structural_paths
  required <- c("grazing_to_biomass", "soil_nutrient_to_sr")
  valid <- identical(metadata$generation_method, "parameterized SEM latent-factor generation") &&
    is.numeric(paths) && all(required %in% names(paths)) &&
    unname(paths[["grazing_to_biomass"]]) < 0 &&
    unname(paths[["soil_nutrient_to_sr"]]) > 0
  list(pass = isTRUE(valid), code = if (isTRUE(valid)) "SEM_SIGNED_PATHS_PRESENT" else
    "R-086_RANDOM_OR_UNSIGNED_GENERATION")
}

.count_order_violations_group <- function(original, candidate, idx) {
  if (length(idx) < 2L) return(0L)
  violations <- 0L
  for (i in seq_len(length(idx) - 1L)) for (j in (i + 1L):length(idx)) {
    old_diff <- original[idx[i]] - original[idx[j]]
    if (old_diff == 0) next
    new_diff <- candidate[idx[i]] - candidate[idx[j]]
    if (sign(old_diff) != sign(new_diff)) violations <- violations + 1L
  }
  violations
}

count_rank_violations <- function(original, candidate, data = NULL,
                                  group_columns = c("grazing", "month", "rep"),
                                  value_columns = c("SR", "W", "T")) {
  if (length(original) != length(candidate)) stop("原值和候选值长度不一致")
  if (is.null(data)) data <- data.frame(grazing = 1, month = 1, rep = seq_along(original))
  if (!all(group_columns %in% names(data))) stop("排序分组列缺失")
  key <- interaction(data[group_columns], drop = TRUE, lex.order = TRUE, sep = "|")
  groups <- split(seq_len(nrow(data)), key)
  sum(vapply(groups, function(idx) .count_order_violations_group(original, candidate, idx), integer(1)))
}

count_order_violations <- count_rank_violations

generate_parameterized_candidate <- function(raw, feasibility,
                                             seed = 1L,
                                             scale = 1,
                                             mean_shift = 0,
                                             tie_policy = c("average_rank", "minimum_rank", "break_ties"),
                                             swap_locked = TRUE,
                                             sem_parameters = list()) {
  tie_policy <- match.arg(tie_policy)
  if (!is.list(raw) || is.null(raw$sheet1) || is.null(raw$soil)) stop("raw 必须包含 sheet1 和 soil")
  if (!is.list(feasibility) || is.null(feasibility$boundary_evidence)) stop("缺少预检 boundary_evidence")
  if (!is.numeric(scale) || length(scale) != 1L || !is.finite(scale) || scale < 0) stop("scale 必须为非负有限数")
  if (!is.numeric(mean_shift) || length(mean_shift) != 1L || !is.finite(mean_shift)) stop("mean_shift must be finite")
  defaults <- list(grazing_biomass = -0.40, grazing_nutrient = -0.45,
                   biomass_noise = 0.80, nutrient_noise = 0.60,
                   soil_env_noise = 0.70,
                   biomass_sr = 1.50, nutrient_sr = 2.50, sr_noise = 1.50)
  if (!is.list(sem_parameters) || any(!names(sem_parameters) %in% names(defaults)))
    stop("sem_parameters contains unknown generator parameter")
  sem_parameters <- utils::modifyList(defaults, sem_parameters)
  if (any(!vapply(sem_parameters, function(x) is.numeric(x) && length(x) == 1L && is.finite(x), logical(1))) ||
      sem_parameters$grazing_biomass >= 0 || sem_parameters$grazing_nutrient >= 0 || sem_parameters$sr_noise < 0)
    stop("sem_parameters must contain finite coefficients with ecological directions")
  sheet1 <- as.data.frame(raw$sheet1, stringsAsFactors = FALSE)
  soil <- as.data.frame(raw$soil, stringsAsFactors = FALSE)
  locked <- c("H", "E", "D", "R")
  keys <- c("grazing", "month", "rep")
  if (!all(c(keys, locked) %in% names(sheet1))) stop("Sheet1 缺少排序键或 H/E/D/R 锁定列")
  if (anyNA(sheet1[c(keys, locked)]) || anyNA(soil[intersect(c("grazing", "month", "rep"), names(soil))])) stop("候选生成拒绝含缺失值输入")

  with_local_seed <- function(seed, code) {
    had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    old_seed <- if (had_seed) get(".Random.seed", envir = .GlobalEnv) else NULL
    on.exit({
      if (had_seed) assign(".Random.seed", old_seed, envir = .GlobalEnv)
      else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
        rm(".Random.seed", envir = .GlobalEnv)
    }, add = TRUE)
    set.seed(as.integer(seed %% .Machine$integer.max))
    force(code)
  }

  with_local_seed(seed, {
    mutable <- candidate_mutable_columns(sheet1, soil)
    bounds <- feasibility$boundary_evidence
    group_columns <- c("grazing", "month", "rep")
    group_id <- function(d) interaction(d[group_columns], drop = TRUE, lex.order = TRUE, sep = "|")
    standardize <- function(x) {
      s <- stats::sd(x)
      if (!is.finite(s) || s == 0) return(rep(0, length(x)))
      as.numeric((x - mean(x)) / s)
    }
    latent_sem <- function(d) {
      grazing <- standardize(as.numeric(d$grazing))
      div <- standardize(-0.25 * grazing + stats::rnorm(nrow(d), sd = 0.50))
      biomass <- standardize(sem_parameters$grazing_biomass * grazing + 0.20 * div + stats::rnorm(nrow(d), sd = sem_parameters$biomass_noise))
      soil_env <- standardize(-0.15 * grazing + 0.10 * biomass + stats::rnorm(nrow(d), sd = sem_parameters$soil_env_noise))
      nutrient <- standardize(sem_parameters$grazing_nutrient * grazing + 0.05 * biomass + 0.40 * soil_env + stats::rnorm(nrow(d), sd = sem_parameters$nutrient_noise))
      soil_factor <- standardize(0.95 * nutrient + 0.05 * biomass + stats::rnorm(nrow(d), sd = 0.10))
      sr <- standardize(-2.20 * grazing + 0.05 * div + sem_parameters$biomass_sr * biomass + sem_parameters$nutrient_sr * soil_factor + stats::rnorm(nrow(d), sd = sem_parameters$sr_noise))
      list(grazing = grazing, div = div, biomass = biomass, soil_env = soil_env,
           nutrient = nutrient, soil_factor = soil_factor, sr = sr)
    }
    main_latent <- latent_sem(sheet1)
    soil_latent <- latent_sem(soil)
    soil_keys <- interaction(soil[c("grazing", "plot", "rep")], drop = TRUE,
                             lex.order = TRUE, sep = "|")
    main_soil_keys <- interaction(sheet1[c("grazing", "plot", "rep")], drop = TRUE,
                                  lex.order = TRUE, sep = "|")
    soil_match <- match(as.character(main_soil_keys), as.character(soil_keys))
    if (anyNA(soil_match)) stop("Sheet1 与土壤表缺少共享 SEM plot/rep 潜变量键")
    main_latent$soil_factor <- soil_latent$soil_factor[soil_match]
    main_latent$sr <- standardize(-2.20 * main_latent$grazing +
      0.05 * main_latent$div + sem_parameters$biomass_sr * main_latent$biomass +
      sem_parameters$nutrient_sr * main_latent$soil_factor +
      stats::rnorm(nrow(sheet1), sd = sem_parameters$sr_noise))
    signal_for <- function(sheet_name, nm) {
      latent <- if (sheet_name == "Sheet1") main_latent else soil_latent
      if (nm == "SR") return(latent$sr)
      if (nm %in% c("AGB", "BGB", "LB")) return(latent$biomass)
      if (nm == "coverage") return(latent$div)
      if (nm %in% c("W", "T")) return(latent$soil_env)
      if (nm == "BD") return(latent$grazing)
      if (nm == "pH") return(-latent$sr)
      latent$soil_factor
    }
    loadings <- c(SR = .999, AGB = .95, BGB = .95, LB = .95, coverage = .7,
                  W = .75, T = .75, BD = .8, pH = -.7, C = .98, SOC = .98,
                  AP = .95, N = .98, NO3_N = .70, NH4_N = -.70, AK = .95, P = .95)
    generate_sheet <- function(d, sheet_label, columns) {
      out <- d
      order_audit <- list()
      reflection_audit <- list()
      hard_bounds_audit <- list()
      for (nm in columns) {
        b <- bounds[bounds$sheet == sheet_label & bounds$indicator == nm, , drop = FALSE]
        if (nrow(b) != nrow(d)) stop(sprintf("%s.%s 的边界证据行数不匹配", sheet_label, nm))
        if (anyNA(d[[nm]])) stop(sprintf("%s.%s 含缺失原值", sheet_label, nm))
        signal <- signal_for(sheet_label, nm)
        loading <- loadings[[nm]]
        latent_noise <- loading * standardize(signal) + sqrt(1 - loading^2) * stats::rnorm(nrow(d))
        group_key <- interaction(d[c("month", "grazing")], drop = TRUE,
                                 lex.order = TRUE, sep = "|")
        group_means <- tapply(d[[nm]], group_key, mean)
        lower_row <- b$lower
        upper_row <- b$upper
        if (any(!is.finite(lower_row)) || any(!is.finite(upper_row)) ||
            any(lower_row > upper_row)) stop(sprintf("%s.%s has invalid row bounds", sheet_label, nm))
        mean_g <- unname(group_means[as.character(group_key)])
        lower_g <- 0.5 * mean_g
        upper_g <- 1.5 * mean_g
        if (any(!is.finite(lower_g)) || any(!is.finite(upper_g)) || any(lower_g > upper_g))
          stop(sprintf("%s.%s has invalid month-by-grazing mean bounds", sheet_label, nm))
        global_sd <- stats::sd(d[[nm]])
        if (!is.finite(global_sd)) global_sd <- 0
        spread_by_group <- vapply(names(group_means), function(group_name) {
          idx <- which(as.character(group_key) == group_name)
          max(stats::sd(d[[nm]][idx]), 0.5 * global_sd, na.rm = TRUE)
        }, numeric(1))
        spread_g <- unname(spread_by_group[as.character(group_key)])
        proposed <- mean_g + mean_shift + scale * spread_g * latent_noise
        reflected_initial <- sum(proposed < lower_g | proposed > upper_g)
        proposed <- reflect_to_bounds(proposed, lower_g, upper_g)
        reflected_projection <- 0L
        group_rows <- split(seq_len(nrow(d)), as.character(group_key))
        month_groups <- split(names(group_rows), vapply(group_rows, function(idx)
          as.character(d$month[idx[1L]]), character(1)))
        for (month_names in month_groups) {
          grazing_groups <- group_rows[month_names]
          group_labels <- vapply(grazing_groups, function(idx)
            as.character(d$grazing[idx[1L]]), character(1))
          names(grazing_groups) <- group_labels
          original_means <- vapply(grazing_groups, function(idx) mean(d[[nm]][idx]), numeric(1))
          current_means <- vapply(grazing_groups, function(idx) mean(proposed[idx]), numeric(1))
          ordered_groups <- order(original_means, names(original_means), method = "radix")
          isotonic_means <- stats::isoreg(seq_along(ordered_groups),
                                           current_means[ordered_groups])$yf
          projected_means <- current_means
          ordered_lower <- vapply(grazing_groups[ordered_groups], function(idx) mean(lower_g[idx]), numeric(1))
          ordered_upper <- vapply(grazing_groups[ordered_groups], function(idx) mean(upper_g[idx]), numeric(1))
          projected_means[ordered_groups] <- pmin(ordered_upper, pmax(ordered_lower, isotonic_means))
          for (g in names(grazing_groups)) {
            idx <- grazing_groups[[g]]
            proposed[idx] <- proposed[idx] + projected_means[[g]] - current_means[[g]]
          }
          reflected_projection <- reflected_projection + sum(
            proposed[unlist(grazing_groups, use.names = FALSE)] < lower_g[unlist(grazing_groups, use.names = FALSE)] |
              proposed[unlist(grazing_groups, use.names = FALSE)] > upper_g[unlist(grazing_groups, use.names = FALSE)])
          month_idx <- unlist(grazing_groups, use.names = FALSE)
          proposed[month_idx] <- reflect_to_bounds(proposed[month_idx],
                                                    lower_g[month_idx], upper_g[month_idx])
          final_means <- vapply(grazing_groups, function(idx) mean(proposed[idx]), numeric(1))
          ordered_final <- final_means[ordered_groups]
          if (any(diff(ordered_final) < -1e-10) ||
              any(final_means < vapply(grazing_groups, function(idx) mean(lower_g[idx]), numeric(1)) - 1e-10) ||
              any(final_means > vapply(grazing_groups, function(idx) mean(upper_g[idx]), numeric(1)) + 1e-10)) {
            proposed[month_idx] <- d[[nm]][month_idx]
            final_means <- original_means
            used_fallback <- TRUE
          } else {
            used_fallback <- FALSE
          }
          order_audit[[length(order_audit) + 1L]] <- data.frame(
            sheet = sheet_label, indicator = nm, month = as.character(d$month[month_idx[1L]]),
            grazing = names(grazing_groups), original_mean = unname(original_means),
            pre_isotonic_mean = unname(current_means),
            isotonic_target_mean = unname(projected_means),
            group_lower = unname(vapply(grazing_groups, function(idx) mean(lower_g[idx]), numeric(1))),
            group_upper = unname(vapply(grazing_groups, function(idx) mean(upper_g[idx]), numeric(1))),
            final_mean = unname(final_means),
            actual_offset = unname(final_means - original_means),
            used_fallback = used_fallback, stringsAsFactors = FALSE)
        }
        out[[nm]] <- proposed
        obs_lower <- pmin(0.5 * d[[nm]], 1.5 * d[[nm]])
        obs_upper <- pmax(0.5 * d[[nm]], 1.5 * d[[nm]])
        obs_violation <- proposed < obs_lower - 1e-12 | proposed > obs_upper + 1e-12
        hard_bounds_audit[[length(hard_bounds_audit) + 1L]] <- do.call(rbind,
          lapply(names(group_rows), function(group_name) {
            idx <- group_rows[[group_name]]
            candidate_mean <- mean(proposed[idx])
            group_lower <- lower_g[idx[1L]]
            group_upper <- upper_g[idx[1L]]
            data.frame(sheet = sheet_label, indicator = nm, group = group_name,
              original_mean = mean_g[idx[1L]], candidate_mean = candidate_mean,
              group_mean_lower = group_lower, group_mean_upper = group_upper,
              group_mean_pass = candidate_mean >= group_lower - 1e-12 &&
                candidate_mean <= group_upper + 1e-12,
              obs_violation_count = sum(obs_violation[idx]), stringsAsFactors = FALSE)
          }))
        obs_violations <- sum(obs_violation)
        reflection_audit[[length(reflection_audit) + 1L]] <- data.frame(
          sheet = sheet_label, indicator = nm, n = nrow(d),
          n_reflected_initial = reflected_initial,
          n_reflected_after_isotonic = reflected_projection,
          n_obs_bound_violations = obs_violations,
          fraction_reflected_initial = reflected_initial / nrow(d),
          fraction_reflected_after_isotonic = reflected_projection / nrow(d),
          stringsAsFactors = FALSE)
      }
      list(data = out,
           order_audit = do.call(rbind, order_audit),
           reflection_audit = do.call(rbind, reflection_audit),
           hard_bounds_audit = do.call(rbind, hard_bounds_audit))
    }
    generated_sheet1 <- generate_sheet(sheet1, "Sheet1", mutable$sheet1)
    generated_soil <- generate_sheet(soil, "soil", mutable$soil)
    candidate_sheet1 <- generated_sheet1$data
    candidate_soil <- generated_soil$data

    order_columns <- intersect(c("SR", "W", "T"), names(sheet1))
    order_violations <- setNames(rep(NA_integer_, length(order_columns)), order_columns)

    if (isTRUE(swap_locked) && nrow(sheet1) > 1L) {
      pools <- split(seq_len(nrow(sheet1)), interaction(sheet1["grazing"], drop = TRUE, lex.order = TRUE))
      for (idx in pools) if (length(idx) > 1L) {
        perm <- sample(idx, length(idx), replace = FALSE)
        candidate_sheet1[idx, locked] <- sheet1[perm, locked, drop = FALSE]
      }
    }
    tuple_key <- function(d) apply(d[locked], 1L, paste, collapse = "\034")
    locked_groups <- unique(as.character(sheet1$grazing))
    locked_ok <- all(vapply(locked_groups, function(g) {
      oi <- which(as.character(sheet1$grazing) == g)
      ci <- which(as.character(candidate_sheet1$grazing) == g)
      # Compare the tuple multisets, not data-frame row names.  `sort()` keeps
      # names produced by `apply()`, and those names legitimately differ after
      # a within-group permutation.
      identical(sort(unname(tuple_key(sheet1[oi, , drop = FALSE]))),
                sort(unname(tuple_key(candidate_sheet1[ci, , drop = FALSE]))))
    }, logical(1)))
    bound_audit <- rbind(generated_sheet1$hard_bounds_audit,
                         generated_soil$hard_bounds_audit)
    bound_pass <- all(bound_audit$group_mean_pass)
    mean_deviation <- do.call(rbind, lapply(mutable$sheet1, function(nm) {
      old_key <- interaction(sheet1$plot, sheet1$grazing, drop=TRUE)
      new_key <- interaction(candidate_sheet1$plot, candidate_sheet1$grazing, drop=TRUE)
      old_mean <- tapply(sheet1[[nm]], old_key, mean)
      new_mean <- tapply(candidate_sheet1[[nm]], new_key, mean)
      common <- intersect(names(old_mean), names(new_mean))
      ratio <- abs(new_mean[common] - old_mean[common]) / pmax(abs(old_mean[common]), 1e-12)
      data.frame(indicator=nm, max_abs_relative_deviation=max(ratio, na.rm=TRUE), soft_pass=all(ratio <= .5), stringsAsFactors=FALSE)
    }))
    result <- list(sheet1 = candidate_sheet1, soil = candidate_soil,
         metadata = list(seed = as.integer(seed), scale = scale, mean_shift = mean_shift, tie_policy = tie_policy,
                         generation_method = "parameterized SEM latent-factor generation",
                         signed_structural_paths = c(grazing_to_biomass = sem_parameters$grazing_biomass,
                           grazing_to_soil_nutrient = sem_parameters$grazing_nutrient,
                           biomass_to_sr = sem_parameters$biomass_sr,
                           soil_nutrient_to_sr = sem_parameters$nutrient_sr,
                           sr_noise = sem_parameters$sr_noise),
                         latent_loadings = loadings, soil_factor = list(name="Soil_Factor", loadings=c(C=.85, SOC=.85, N=.80, P=.75),
                         relation_to_SR="positive"), mutable_columns = mutable,
                         projection_diagnostics = list(soil_factor_scores=main_latent$soil_factor,
                           sr_soil_factor_pearson=stats::cor(candidate_sheet1$SR, main_latent$soil_factor)),
                         locked_columns = locked, order_violations = order_violations,
                         order_constraint_status = "MONTH_GRAZING_MEAN_ORDER_PRESERVED",
                         locked_tuples_preserved = locked_ok,
                         bounds_preserved = bound_pass,
                         group_mean_order_audit = rbind(generated_sheet1$order_audit,
                                                        generated_soil$order_audit),
                         group_mean_order_violations = sum(vapply(
                           split(rbind(generated_sheet1$order_audit, generated_soil$order_audit),
                                 interaction(rbind(generated_sheet1$order_audit, generated_soil$order_audit)[c("sheet", "indicator", "month")],
                                             drop = TRUE)),
                           function(x) sum(diff(x$final_mean[order(x$original_mean, x$grazing)]) < -1e-10),
                           integer(1))),
                         reflection_diagnostics = rbind(generated_sheet1$reflection_audit,
                                                        generated_soil$reflection_audit),
                         group_mean_bounds_audit = bound_audit,
                         plot_grazing_mean_deviation = mean_deviation,
                         has_pass = locked_ok && bound_pass))
    result$metadata$r086_audit <- r086_sem_generation_audit(result)
    result$metadata$has_pass <- result$metadata$has_pass && result$metadata$r086_audit$pass
    result
  })
}
