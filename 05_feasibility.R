precheck_feasibility <- function(raw, constraints = list(), log_message = function(...) invisible(NULL)) {
  `%||%` <- function(x, y) if (is.null(x)) y else x
  if (is.null(raw$sheet1) || is.null(raw$soil)) stop("raw 必须包含 sheet1 和 soil")
  main <- as.data.frame(raw$sheet1)
  soil <- as.data.frame(raw$soil)
  key <- c("grazing", "month", "rep")
  locked <- c("H", "E", "D", "R")
  mutable <- c("SR", "W", "T", "AGB", "BGB", "LB", "coverage", "BD", "pH", "C", "SOC", "AP", "N", "NO3_N", "NH4_N", "AK", "P")
  required <- c(key, "plot", locked, intersect(mutable, names(main)), intersect(mutable, names(soil)))
  missing <- setdiff(required, c(names(main), names(soil)))
  if (length(missing)) stop(sprintf("预检缺少必需列: %s", paste(missing, collapse = ",")))
  if (anyNA(main[c(key, locked, intersect(mutable, names(main)))]) || anyNA(soil[c(key, intersect(mutable, names(soil)))])) {
    stop("预检输入含缺失值；拒绝静默填补")
  }
  if (!all(vapply(main[c("grazing", "month", "rep", locked)], is.numeric, logical(1))) ||
      !all(vapply(soil[c("grazing", "month", "rep")], is.numeric, logical(1)))) stop("分层键及锁定四元组必须为数值列")

  values <- list(Sheet1 = main, soil = soil)
  bounds <- list()
  for (sheet in names(values)) {
    d <- values[[sheet]]
    cols <- intersect(mutable, names(d))
    for (nm in cols) {
      x <- d[[nm]]
      if (!is.numeric(x)) stop(sprintf("可改指标 %s.%s 必须为数值列", sheet, nm))
      physical <- constraints$physical_bounds[[nm]]
      if (any(x == 0)) {
        if (is.null(physical) || length(physical) != 2L || any(!is.finite(physical)) || physical[1] > physical[2] || physical[1] < 0) {
          stop(sprintf("指标 %s 含零值，必须提供有效非负物理边界", nm))
        }
        lo <- pmin(0.5 * x, 1.5 * x)
        hi <- pmax(0.5 * x, 1.5 * x)
        lo[x == 0] <- physical[1]; hi[x == 0] <- physical[2]
        source <- ifelse(x == 0, "constraints$physical_bounds", "原值端点排序")
      } else {
        lo <- pmin(0.5 * x, 1.5 * x)
        hi <- pmax(0.5 * x, 1.5 * x)
        source <- ifelse(x < 0, "原值端点排序", "原值±50%")
      }
      if (!is.null(physical)) {
        if (length(physical) != 2L || any(!is.finite(physical)) || physical[1] > physical[2]) stop(sprintf("指标 %s 的物理边界无效", nm))
        lo <- pmax(lo, physical[1]); hi <- pmin(hi, physical[2])
        if (any(lo > hi)) { lo[lo > hi] <- NA_real_; hi[lo > hi] <- NA_real_ }
        source <- paste0(source, "; 与物理范围取交集")
      }
      bounds[[paste(sheet, nm, sep = ".")]] <- data.frame(sheet = sheet, row = seq_along(x), indicator = nm, old = x, lower = lo, upper = hi, source = source, stringsAsFactors = FALSE)
    }
  }
  boundary_evidence <- do.call(rbind, bounds)

  group_key <- function(d, cols) interaction(d[cols], drop = TRUE, lex.order = TRUE, sep = "|")
  distortion <- list(); k <- 0L
  for (sheet in names(values)) for (nm in intersect(mutable, names(values[[sheet]]))) {
    d <- values[[sheet]]; groups <- split(seq_len(nrow(d)), group_key(d, c(key)))
    for (gn in names(groups)) {
      x <- d[[nm]][groups[[gn]]]; n <- length(x); variance <- if (n > 1) stats::var(x) else NA_real_
      diffs <- diff(sort(x)); diff_ratio <- if (is.finite(variance) && variance > 0) if (length(diffs) > 1) stats::var(diffs) / variance else 0 else NA_real_
      skew <- if (n >= 3 && stats::sd(x) > 0) mean((x-mean(x))^3) / stats::sd(x)^3 else NA_real_
      kurt <- if (n >= 4 && stats::sd(x) > 0) mean((x-mean(x))^4) / stats::sd(x)^4 - 3 else NA_real_
      k <- k + 1L
      distortion[[k]] <- data.frame(sheet=sheet, group=gn, indicator=nm, unique_n=length(unique(x)), diff_variance_ratio=diff_ratio, skewness=skew, excess_kurtosis=kurt,
        unique_flag=length(unique(x)) < 3, diff_flag=is.finite(diff_ratio) && diff_ratio < .01, skew_flag=is.finite(skew) && abs(skew)>3, kurtosis_flag=is.finite(kurt) && abs(kurt)>10, status="SOFT_DIAGNOSTIC", stringsAsFactors=FALSE)
    }
  }
  distortion <- do.call(rbind, distortion)

  # A swap pool preserves each observed H/E/D/R tuple and is restricted to grazing.
  tuple <- do.call(paste, c(lapply(main[locked], as.character), sep="\r"))
  pool_key <- interaction(main$grazing, drop=TRUE, lex.order=TRUE)
  pools <- split(seq_len(nrow(main)), pool_key)
  swap_pools <- do.call(rbind, lapply(names(pools), function(g) {
    idx <- pools[[g]]
    data.frame(grazing=main$grazing[idx[1]], n=nrow(data.frame(i=idx)), distinct_tuples=length(unique(tuple[idx])),
      swap_degrees=max(0L, length(idx)-1L), tuple_preserved=TRUE, stringsAsFactors=FALSE)
  }))
  strata <- split(seq_len(nrow(main)), group_key(main, key))
  low_mobility <- names(strata)[vapply(strata, function(idx) length(unique(tuple[idx])) <= 1L, logical(1))]
  lock_evidence <- data.frame(pool=names(pools), grazing=vapply(pools, function(i) as.character(main$grazing[i[1]]), character(1)),
    rows=lengths(pools), distinct_tuples=vapply(pools, function(i) length(unique(tuple[i])), integer(1)),
    movable_rows=vapply(pools, function(i) if(length(unique(tuple[i]))>1) length(i) else 0L, integer(1)), stringsAsFactors=FALSE)

  direction_rows <- lapply(names(values), function(sheet_name) {
    d <- values[[sheet_name]]
    vars <- intersect(c("SR","AGB","BGB","BD","C","SOC","AP","N","NO3_N","NH4_N","AK","P"), names(d))
    if (!length(vars)) return(NULL)
    do.call(rbind, lapply(vars, function(nm) {
      means <- aggregate(d[[nm]], list(grazing=d$grazing), mean); names(means)[2] <- "mean"
      r <- suppressWarnings(stats::cor(d$grazing, d[[nm]]))
      data.frame(sheet=sheet_name, indicator=nm,
        grazing_mean_order=paste(means$grazing[order(means$mean)], collapse="<"),
        correlation=r,
        inferred_direction=if(is.finite(r) && r>0) "positive" else if(is.finite(r) && r<0) "negative" else "undetermined",
        stringsAsFactors=FALSE)
    }))
  })
  direction <- do.call(rbind, direction_rows)

  constraints_evidence <- data.frame(requirement=c("R-087","R-088","R-089","R-090","R-091","R-092","R-093"),
    target=c("低维摘要与换位/耦合审计","锁定池自由度","区间/线性约束可行性","新排序/边界/失真约束","PLS骨架联合可行性","PLS及测量模型预检","两两及三元最小冲突核"),
    status=c("EVALUATED","EVALUATED","NOT_EVALUATED","NOT_EVALUATED","NOT_EVALUATED","NOT_EVALUATED","NOT_EVALUATED"),
    evidence=c("组均值、组内残差协方差及换位池摘要由输入数据给出","已列出放牧内四元组自由度；零自由度层: ",
      "未配置真实联合求解器；结论 NOT_ESTABLISHED","仅输出输入边界；候选排序尚未执行",
      "未运行代理抽样或真实求解；NOT_ESTABLISHED",
      "未运行完整 PLS/Bootstrap；NOT_ESTABLISHED",
      "未运行真实两两/三元联合求解；NOT_EVALUATED"), stringsAsFactors=FALSE)
  constraints_evidence$evidence[2] <- paste0(constraints_evidence$evidence[2], paste(low_mobility, collapse=","))
  hard_ids <- c("R-088", "R-089", "R-090", "R-091", "R-092")
  conflict_matrix <- expand.grid(row_requirement=hard_ids, column_requirement=hard_ids,
                                 stringsAsFactors=FALSE)
  conflict_matrix$relation <- ifelse(conflict_matrix$row_requirement == conflict_matrix$column_requirement,
                                     "单项检查", "待联合抽样判定")
  conflict_matrix$evidence <- ""
  conflict_matrix$hardness <- "HARD"
  conflict_matrix$status <- "EVALUATED"

  # Constraint evidence is derived from observed inputs. A proxy sample may
  # not manufacture R2/VIF/loading values or certify formal feasibility.
  sample_n <- 0L
  sampled_ranges <- data.frame(metric=c("r2","vif","loading","AVE","CR","path_sign"),
    minimum=NA_real_, maximum=NA_real_, samples=0L, seed=NA_integer_,
    target_reachable=NA, status="NOT_ESTABLISHED",
    evidence="No real joint solver executed; proxy sampling is not run", stringsAsFactors=FALSE)
  # This deterministic audit proves only order feasibility for supplied
  # intervals; it does not imply statistical or PLS reachability.
  order_evidence <- do.call(rbind, lapply(c("SR", "W", "T"), function(nm) {
    if (!nm %in% names(main)) return(NULL)
    b <- boundary_evidence[boundary_evidence$sheet == "Sheet1" & boundary_evidence$indicator == nm, , drop = FALSE]
    rows <- split(seq_len(nrow(main)), group_key(main, key))
    do.call(rbind, lapply(names(rows), function(g) {
      idx <- rows[[g]][order(main[[nm]][rows[[g]]], rows[[g]])]
      lower <- b$lower[idx]; upper <- b$upper[idx]
      data.frame(indicator = nm, group = g, rows = length(idx),
                 feasible = all(cummax(lower) <= upper + 1e-12),
                 max_lower_upper_gap = max(cummax(lower) - upper),
                 stringsAsFactors = FALSE)
    }))
  }))
  conflict_matrix$relation <- "NOT_ESTABLISHED"
  conflict_matrix$evidence <- "No real joint solver executed; pair feasibility is NOT_EVALUATED"
  conflict_matrix$status <- "NOT_EVALUATED"
  conflict_matrix$joint_samples <- NA_integer_
  triples <- utils::combn(hard_ids, 3L, simplify=FALSE)
  triple_evidence <- do.call(rbind, lapply(triples, function(x) {
    data.frame(constraints=paste(x, collapse=" + "), individually_feasible=NA,
               pairwise_feasible=NA, triple_feasible=NA, joint_samples=NA_integer_, sample_n=0L,
               status="NOT_EVALUATED",
               evidence="No real joint solver executed; triple feasibility is NOT_ESTABLISHED", stringsAsFactors=FALSE)
  }))
  numeric_cols <- intersect(mutable, names(main))
  group_means <- aggregate(main[numeric_cols], main[key], mean)
  residual_covariance <- lapply(split(seq_len(nrow(main)), group_key(main,key)), function(idx) {
    x <- main[idx, intersect(numeric_cols,names(main)), drop=FALSE]
    if(nrow(x)<2 || ncol(x)<1) return(data.frame(group=character(), indicator_1=character(), indicator_2=character(), covariance=numeric()))
    cm <- stats::cov(x)
    ij <- which(upper.tri(cm, diag=TRUE), arr.ind=TRUE)
    data.frame(group=as.character(group_key(main,key)[idx[1]]), indicator_1=rownames(cm)[ij[,1]], indicator_2=colnames(cm)[ij[,2]], covariance=cm[ij], stringsAsFactors=FALSE)
  })
  residual_covariance <- do.call(rbind, residual_covariance)
  order_hard_conflict <- any(!order_evidence$feasible, na.rm=TRUE)
  pls_skeleton <- list(
    blocks_present = all(c("R", "H", "E", "D", "AGB", "BGB", "LB", "BD", "W", "T",
                           "pH", "C", "SOC", "AP", "N", "NO3_N", "NH4_N", "AK", "P", "SR") %in% names(main)),
    critical_paths_declared = TRUE,
    status = if (order_hard_conflict) "BLOCKED_BY_ORDER_INTERVAL" else "STRUCTURALLY_REPRESENTABLE",
    note = "Full loading/R2/VIF/bootstrap are evaluated in the PLS proxy/full stages."
  )
  boundary_conflict <- any(!is.finite(boundary_evidence$lower)) || any(!is.finite(boundary_evidence$upper)) || any(boundary_evidence$lower > boundary_evidence$upper, na.rm=TRUE)
  lock_conflict <- any(lock_evidence$rows < 2L & lock_evidence$movable_rows > 0L)
  hard_conflict <- order_hard_conflict || lock_conflict || boundary_conflict
  precheck_status <- if (hard_conflict) "INFEASIBLE" else "NOT_ESTABLISHED"
  core_rows <- list()
  add_core <- function(pair, violation, hardness="HARD", relaxation) {
    core_rows[[length(core_rows)+1L]] <<- data.frame(constraint_pair=pair, violation=violation,
      hardness=hardness, relaxation_order=relaxation, stringsAsFactors=FALSE)
  }
  if(lock_conflict) add_core("R-088 locked mobility x R-089 ordering", sum(lock_evidence$rows < 2L), relaxation="补充同 grazing 可动四元组观测")
  if(order_hard_conflict) add_core("R-089 SR/W/T ordering x R-090 bounds", max(order_evidence$max_lower_upper_gap, na.rm=TRUE), relaxation="复核排序目标后再考虑调整排序要求")
  if(boundary_conflict) add_core("R-090 bounds x physical bounds", sum(!is.finite(boundary_evidence$lower) | !is.finite(boundary_evidence$upper)), relaxation="仅凭数据/领域依据调整物理范围")
  conflict_core <- if(length(core_rows)) do.call(rbind, core_rows) else data.frame(constraint_pair=character(), violation=numeric(), hardness=character(), relaxation_order=character())
  conflict_core$relax_object <- if(nrow(conflict_core)) c("锁定池自由度", "排序/区间约束", "物理边界") [seq_len(nrow(conflict_core))] else character()
  conflict_core$relax_amount <- if(nrow(conflict_core)) c("增加至 n>=2", "移除/调整冲突排序", "按证据扩展区间")[seq_len(nrow(conflict_core))] else character()
  conflict_core$expected_impact <- if(nrow(conflict_core)) c("增加可交换自由度", "改变候选排序关系", "扩大可取值域")[seq_len(nrow(conflict_core))] else character()
  conflict_core$infeasibility_evidence <- if(nrow(conflict_core)) paste0(precheck_status, "; see precheck order/lock/boundary evidence") else character()
  conflict_core$risk <- if(nrow(conflict_core)) "会改变数据约束；完整验收前不得接受" else character()
  alternatives <- data.frame(option=1:3, relax_object=c("排序约束", "锁定池自由度", "物理边界"),
    relax_amount=c("调整一个冲突排序关系", "同 grazing 池补至至少 2 行", "仅按独立证据扩展区间"),
    expected_impact=c("改变 SR/W/T 的目标关系", "增加可置换候选", "扩大候选值空间"),
    infeasibility_evidence=if(hard_conflict) paste(conflict_core$constraint_pair, collapse="; ") else "无硬冲突",
    risk="改变约束语义；必须完整验收", stringsAsFactors=FALSE)
  alternative_rows <- data.frame(constraint_pair=paste0("ALTERNATIVE_", alternatives$option),
    violation=0, hardness="SOFT_PROPOSAL", relaxation_order=as.character(alternatives$option),
    relax_object=alternatives$relax_object, relax_amount=alternatives$relax_amount,
    expected_impact=alternatives$expected_impact,
    infeasibility_evidence=alternatives$infeasibility_evidence, risk=alternatives$risk,
    stringsAsFactors=FALSE)
  conflict_core <- rbind(conflict_core, alternative_rows)
  log_message(sprintf("FEASIBILITY status=%s boundaries=%d distortion_flags=%d order_conflict=%s pls_skeleton=%s",
    precheck_status, nrow(boundary_evidence),
    sum(distortion$unique_flag|distortion$diff_flag|distortion$skew_flag|distortion$kurtosis_flag),
    order_hard_conflict, pls_skeleton$status))
  result <- list(status=precheck_status, group_means=group_means, residual_covariance=residual_covariance, swap_pools=swap_pools,
    lock_evidence=lock_evidence, low_mobility_strata=low_mobility, boundary_evidence=boundary_evidence, distortion=distortion,
    order_evidence=order_evidence, pls_skeleton=pls_skeleton,
    direction_evidence=direction, conflict_matrix=conflict_matrix, triple_evidence=triple_evidence,
    sampled_ranges=sampled_ranges, constraints_evidence=constraints_evidence,
    conflict_core=conflict_core, alternatives=alternatives,
    generator_status=if(hard_conflict) "BLOCKED" else "READY")
  core_path <- constraints$conflict_core_path %||% "conflict_core.csv"
  utils::write.csv(result$conflict_core, core_path, row.names=FALSE, na="")
  utils::write.csv(result$alternatives, file.path(dirname(core_path), "precheck_alternatives.csv"), row.names=FALSE, na="")
  grDevices::png(file.path(dirname(core_path), "conflict_heatmap.png"), width=900, height=700)
  on.exit(grDevices::dev.off(), add=TRUE)
  image(matrix(as.numeric(conflict_matrix$status == "NO_SAMPLE_FOUND"), nrow=length(hard_ids)), axes=FALSE,
        col=c("white", "firebrick"), main="Precheck pairwise conflict diagnostics")
  axis(1, at=seq(0,1,length.out=length(hard_ids)), labels=hard_ids, las=2)
  axis(2, at=seq(0,1,length.out=length(hard_ids)), labels=hard_ids, las=2)
  result
}
