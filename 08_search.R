# Module 8: population-based global search and distribution-changing restarts.
# This module is independent of the earlier clean-room modules and accepts
# candidate/evaluation functions as explicit dependencies.

local_rng <- function(seed, expr) {
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  old_seed <- if (had_seed) get(".Random.seed", envir = .GlobalEnv) else NULL
  on.exit({
    if (had_seed) assign(".Random.seed", old_seed, envir = .GlobalEnv)
    else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
      rm(".Random.seed", envir = .GlobalEnv)
  }, add = TRUE)
  set.seed(as.integer(seed %% .Machine$integer.max))
  force(expr)
}

restart_distribution <- function(distribution, restart_id) {
  stopifnot(is.list(distribution), length(distribution) > 0)
  out <- distribution
  # Every restart changes at least one generation characteristic.
  out$mean_shift <- (distribution$mean_shift %||% 0) + restart_id * 0.15
  out$scale <- (distribution$scale %||% 1) * (1 + 0.20 * restart_id)
  out$covariance_id <- paste0(distribution$covariance_id %||% "base", "-r", restart_id)
  out$restart_id <- restart_id
  out$changed_fields <- c("mean_shift", "scale", "covariance_id")
  out
}

`%||%` <- function(x, y) if (is.null(x)) y else x

validate_evaluation <- function(x) {
  if (!is.list(x) || is.null(x$score) || !is.finite(x$score))
    stop("evaluator must return a finite score")
  if (is.null(x$residuals)) x$residuals <- numeric()
  if (is.null(x$gates)) x$gates <- list()
  x$sign_violations <- as.integer(x$sign_violations %||% 0L)
  x
}

dominates <- function(a, b) {
  ra <- unlist(a$residuals, use.names = FALSE)
  rb <- unlist(b$residuals, use.names = FALSE)
  if (!length(ra) || !length(rb) || length(ra) != length(rb)) return(FALSE)
  all(ra <= rb) && any(ra < rb) && a$score >= b$score
}

run_global_search <- function(candidate_generator, evaluator,
                              config = list(), initial_distribution = list()) {
  stopifnot(is.function(candidate_generator), is.function(evaluator))
  cfg <- modifyList(list(iterations = 100L, population_size = 20L,
                         elite_size = 5L, stagnation_limit = 50L,
                         seed = 1L, base_score = 0, sign_penalty = 100,
                         temperature_start = 1, temperature_decay = .95,
                         score_floor = -1e9), config)
  stopifnot(cfg$iterations > 0, cfg$population_size > 0,
            cfg$elite_size > 0, cfg$stagnation_limit > 0)
  dist <- modifyList(list(mean_shift = 0, scale = 1, covariance_id = "base"),
                     initial_distribution)
  records <- vector("list", cfg$iterations * cfg$population_size)
  best <- NULL; current <- NULL; elites <- list(); stagnant <- 0L; restart_id <- 0L; k <- 0L
  max_residual_name <- function(e) {
    r <- unlist(e$residuals, use.names = TRUE)
    if (!length(r)) return(NA_character_)
    names(r)[which.max(r)]
  }
  for (iteration in seq_len(cfg$iterations)) {
    if (stagnant >= cfg$stagnation_limit) {
      restart_id <- restart_id + 1L
      dist <- restart_distribution(dist, restart_id)
      stagnant <- 0L
    }
    iteration_seed <- as.integer((cfg$seed + iteration * 104729L + restart_id * 1009L) %% .Machine$integer.max)
    candidates <- local_rng(iteration_seed,
      candidate_generator(cfg$population_size, dist, iteration_seed,
                          restart_id))
    if (!is.list(candidates) || length(candidates) != cfg$population_size)
      stop("candidate_generator must return a list of population_size candidates")
    evaluated <- lapply(seq_along(candidates), function(candidate_index) {
      candidate <- candidates[[candidate_index]]
      e <- validate_evaluation(evaluator(candidate))
      e$candidate <- candidate
      e$candidate_id <- paste0("i", iteration, "-c", candidate_index)
      e$parent_id <- if (is.null(current)) NA_character_ else current$candidate_id
      e$raw_score <- e$score
      e$score <- max(cfg$score_floor, cfg$base_score + e$score - cfg$sign_penalty * e$sign_violations)
      e
    })
    scores <- vapply(evaluated, `[[`, numeric(1), "score")
    order_idx <- order(scores, decreasing = TRUE)
    pool <- evaluated[order_idx[seq_len(min(cfg$elite_size, length(order_idx)))]]
    # Keep non-dominated residual profiles in addition to score elites.
    for (e in evaluated) if (!any(vapply(elites, dominates, logical(1), b = e))) elites <- c(elites, list(e))
    if (length(elites) > cfg$elite_size * 4L) {
      elites <- elites[order(vapply(elites, `[[`, numeric(1), "score"), decreasing = TRUE)[seq_len(cfg$elite_size * 4L)]]
    }
    iteration_best <- pool[[1L]]
    temperature <- cfg$temperature_start * cfg$temperature_decay^(iteration - 1L)
    accept <- is.null(current) || iteration_best$score >= current$score
    if (!accept && is.finite(temperature) && temperature > 0) {
      u <- local_rng(iteration_seed + 65537L, stats::runif(1L))
      accept <- u < exp((iteration_best$score - current$score) / temperature)
    }
    if (accept) {
      current <- iteration_best
      dist$mean_shift <- current$candidate$metadata$mean_shift %||% dist$mean_shift
      dist$scale <- current$candidate$metadata$scale %||% dist$scale
    }
    improved <- is.null(best) || iteration_best$score > best$score
    if (improved) { best <- iteration_best; stagnant <- 0L } else stagnant <- stagnant + 1L
    for (j in seq_along(evaluated)) {
      e <- evaluated[[j]]; k <- k + 1L
      gate_values <- unlist(e$gates, use.names = TRUE)
      residual_values <- unlist(e$residuals, use.names = TRUE)
      gate_or <- function(primary, fallback) {
        if (primary %in% names(gate_values)) isTRUE(gate_values[[primary]])
        else if (fallback %in% names(gate_values)) isTRUE(gate_values[[fallback]])
        else FALSE
      }
      records[[k]] <- data.frame(iteration = iteration, candidate_id = j,
        iteration_seed = iteration_seed, restart_id = restart_id,
        elite_pool_id = paste0("r", restart_id, "-i", iteration),
        candidate_source = "stratified_population", accept_flag = accept && j == order_idx[1L],
        score = e$score, sign_violations = e$sign_violations,
        max_residual_item = max_residual_name(e),
        residual_vector = paste(sprintf("%s=%.8g", names(residual_values), residual_values), collapse = ";"),
        has_pass = isTRUE(e$has_pass),
        pass1 = gate_or("has_lock", "has_measurement"),
        pass2 = gate_or("has_order", "has_r2"),
        pass3 = gate_or("has_ecology", "has_sign"),
        max_vif = as.numeric(e$max_vif %||% NA_real_),
        r2_proxy = as.numeric(e$r2_proxy %||% NA_real_),
        pls_proxy_diff = as.numeric(e$pls_proxy_diff %||% NA_real_),
        parent_seed = iteration_seed, parent_id = e$parent_id,
        temperature=temperature, annealing_accepted=accept && j == order_idx[1L],
        stringsAsFactors = FALSE)
    }
  }
  list(best = best, elites = elites, records = do.call(rbind, records[seq_len(k)]),
       config = cfg, restarts = restart_id, run_ok = TRUE)
}
