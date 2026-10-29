# Teaching simulation only; not for real scientific inference.

order_lock_bounds <- function(data, columns, physical_bounds = list()) {
  if (!is.data.frame(data)) stop("data must be a data.frame", call. = FALSE)
  missing_columns <- setdiff(columns, names(data))
  if (length(missing_columns)) stop("missing columns: ", paste(missing_columns, collapse = ", "), call. = FALSE)
  result <- lapply(columns, function(column) {
    old <- data[[column]]
    if (!is.numeric(old)) stop("column must be numeric: ", column, call. = FALSE)
    if (anyNA(old) || any(!is.finite(old))) stop("missing or non-finite values in: ", column, call. = FALSE)
    lower <- upper <- old
    positive <- old > 0
    negative <- old < 0
    lower[positive] <- old[positive] * 0.5
    upper[positive] <- old[positive] * 1.5
    endpoints <- cbind(old[negative] * 0.5, old[negative] * 1.5)
    lower[negative] <- apply(endpoints, 1L, min)
    upper[negative] <- apply(endpoints, 1L, max)
    zero <- which(old == 0)
    if (length(zero)) {
      physical <- physical_bounds[[column]]
      if (is.null(physical) || length(physical) != 2L || any(!is.finite(physical)) || physical[1L] < 0 || physical[1L] > physical[2L]) {
        stop("nonnegative physical bounds required for zero-valued column: ", column, call. = FALSE)
      }
      lower[zero] <- physical[1L]
      upper[zero] <- physical[2L]
    }
    data.frame(row = seq_along(old), column = column, old = old, lower = lower, upper = upper)
  })
  do.call(rbind, result)
}

reflect_bounds <- function(x, lower, upper) {
  if (length(x) != length(lower) || length(x) != length(upper) || anyNA(c(x, lower, upper)) || any(!is.finite(c(x, lower, upper)))) {
    stop("x, lower, and upper must be equal-length finite vectors", call. = FALSE)
  }
  if (any(lower > upper)) stop("lower bound exceeds upper bound", call. = FALSE)
  width <- upper - lower
  fixed <- width == 0
  result <- x
  result[fixed] <- lower[fixed]
  movable <- !fixed
  offset <- (x[movable] - lower[movable]) %% (2 * width[movable])
  result[movable] <- lower[movable] + ifelse(offset <= width[movable], offset, 2 * width[movable] - offset)
  result
}

validate_tuple_lock <- function(original, candidate, grazing = "grazing", tuple = c("H", "E", "D", "R")) {
  required <- c(grazing, tuple)
  if (!is.data.frame(original) || !is.data.frame(candidate)) stop("original and candidate must be data.frames", call. = FALSE)
  if (!all(required %in% names(original)) || !all(required %in% names(candidate))) stop("missing grazing or tuple columns", call. = FALSE)
  if (nrow(original) != nrow(candidate)) stop("row count changed", call. = FALSE)
  if (anyNA(original[required]) || anyNA(candidate[required])) stop("missing tuple values", call. = FALSE)
  key <- function(x) do.call(paste, c(lapply(x[tuple], as.character), sep = "\034"))
  groups <- unique(as.character(original[[grazing]]))
  violations <- 0L
  for (group in groups) {
    oi <- which(as.character(original[[grazing]]) == group)
    ci <- which(as.character(candidate[[grazing]]) == group)
    if (!identical(sort(key(original[oi, , drop = FALSE])), sort(key(candidate[ci, , drop = FALSE])))) violations <- violations + 1L
  }
  list(pass = violations == 0L, violations = violations)
}

count_order_violations <- function(original, candidate, groups, columns = c("SR", "W", "T")) {
  required <- c(groups, columns)
  if (!is.data.frame(original) || !is.data.frame(candidate) || !all(required %in% names(original)) || !all(required %in% names(candidate))) stop("missing data or required columns", call. = FALSE)
  if (nrow(original) != nrow(candidate)) stop("row count changed", call. = FALSE)
  if (anyNA(original[required]) || anyNA(candidate[required])) stop("missing order data", call. = FALSE)
  group_key <- function(data) do.call(paste, c(lapply(data[groups], as.character), sep = "\034"))
  key <- group_key(original)
  if (!identical(key, group_key(candidate))) stop("group membership/order changed", call. = FALSE)
  violations <- 0L
  for (column in columns) {
    if (!is.numeric(original[[column]]) || !is.numeric(candidate[[column]])) stop("order columns must be numeric", call. = FALSE)
    for (idx in split(seq_len(nrow(original)), key)) {
      if (length(idx) < 2L) next
      pairs <- combn(idx, 2L)
      before <- sign(original[[column]][pairs[1L, ]] - original[[column]][pairs[2L, ]])
      after <- sign(candidate[[column]][pairs[1L, ]] - candidate[[column]][pairs[2L, ]])
      violations <- violations + sum(before != 0 & before != after)
    }
  }
  as.integer(violations)
}

rank_map_candidates <- function(original, candidate_values, groups, tie_policy = c("average", "minimum", "break")) {
  tie_policy <- match.arg(tie_policy)
  if (length(original) != length(candidate_values) || anyNA(c(original, candidate_values)) || any(!is.finite(c(original, candidate_values)))) stop("rank mapping requires equal-length finite vectors", call. = FALSE)
  if (length(groups) != length(original) || anyNA(groups)) stop("groups must align with observations and contain no missing values", call. = FALSE)
  out <- candidate_values
  for (idx in split(seq_along(original), as.character(groups))) {
    if (length(idx) < 2L) next
    method <- if (tie_policy == "average") "average" else if (tie_policy == "minimum") "min" else "first"
    ranks <- rank(original[idx], ties.method = method)
    order_values <- sort(candidate_values[idx])
    if (tie_policy == "break") {
      mapped <- numeric(length(idx))
      mapped[order(ranks)] <- order_values
    } else {
      mapped <- approx(seq_along(order_values), order_values, xout = ranks, ties = "ordered", rule = 2)$y
    }
    out[idx] <- mapped
  }
  out
}
