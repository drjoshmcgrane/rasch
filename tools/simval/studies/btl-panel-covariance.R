# Numerical conformance of the BTL-EFRM conditional panel covariance.
# This is not a null-calibration or coverage study. It holds the established
# block-precision point estimator fixed and independently stacks the fitted
# stage-one judge scores to verify its shared-judge sandwich covariance.

suppressWarnings(pkgload::load_all(".", quiet = TRUE))
source("tools/simval/harness.R")

STUDY <- "btl-panel-covariance"
NREP <- as.integer(Sys.getenv("SV_REPS", "10"))
if (length(NREP) != 1L || is.na(NREP) || NREP < 1L)
  stop("SV_REPS must be a positive integer")

remap_judges <- function(d, truth, allocation) {
  sa <- names(truth$object_sets)[truth$set_of[d$object_a]]
  sb <- names(truth$object_sets)[truth$set_of[d$object_b]]
  within <- sa == sb
  original <- d$judge
  if (allocation == "shared") {
    mapped <- original
  } else if (allocation == "disjoint") {
    mapped <- paste(original, ifelse(within, sa, "cross"), sep = ":")
  } else {
    panel_members <- split(names(truth$judge_panel), truth$judge_panel)
    shared <- unlist(lapply(panel_members, function(z)
      sort(z)[seq_len(floor(length(z) / 2L))]), use.names = FALSE)
    mapped <- ifelse(original %in% shared, original,
                     paste(original, ifelse(within, sa, "cross"), sep = ":"))
  }
  d$judge <- mapped
  list(data = d, set_a = sa, set_b = sb, within = within)
}

plant_outcomes <- function(layout, truth, dependence) {
  d <- layout$data
  u <- numeric(nrow(d))
  within_rows <- which(layout$within)
  keys <- unique(data.frame(judge = d$judge[within_rows],
                            set = layout$set_a[within_rows],
                            stringsAsFactors = FALSE))
  ue <- numeric(nrow(keys))
  for (judge in unique(keys$judge)) {
    take <- which(keys$judge == judge)
    z <- rnorm(1L, 0, 0.85)
    if (length(take) == 2L && dependence != 0) {
      take <- take[order(keys$set[take])]
      ue[take] <- c(z, dependence * z)
    } else {
      ue[take] <- rnorm(length(take), 0, 0.85)
    }
  }
  row_key <- paste(d$judge[within_rows], layout$set_a[within_rows], sep = "\r")
  key <- paste(keys$judge, keys$set, sep = "\r")
  u[within_rows] <- ue[match(row_key, key)]
  lp <- truth$v[d$object_a] - truth$v[d$object_b]
  lp[within_rows] <- exp(u[within_rows]) *
    (truth$beta[d$object_a[within_rows]] -
       truth$beta[d$object_b[within_rows]])
  y <- rbinom(nrow(d), 1L, plogis(lp))
  d$winner <- ifelse(y == 1L, d$object_a, d$object_b)
  d
}

# Rebuild the stage-one score contributions rather than consuming the
# implementation's returned influence_lrho field.
stage_block <- function(d, objects) {
  rows <- d$object_a %in% objects & d$object_b %in% objects
  ia <- match(d$object_a[rows], objects)
  ib <- match(d$object_b[rows], objects)
  y <- as.integer(d$winner[rows] == d$object_a[rows])
  panel <- as.character(d$panel[rows])
  judge <- as.character(d$judge[rows])
  K <- length(objects)
  f <- rasch:::.btlef_stage1(ia, ib, y, panel, judge, K, 60, 1e-8)

  B <- rbind(diag(K - 1L), rep(-1, K - 1L))
  Bd <- B[ia, , drop = FALSE] - B[ib, , drop = FALSE]
  rho <- unname(f$rho[panel])
  delta <- f$beta[ia] - f$beta[ib]
  Jr <- matrix(0, length(y), length(f$free))
  if (length(f$free)) for (g in seq_along(f$free)) {
    use <- panel == f$free[g]
    Jr[use, g] <- rho[use] * delta[use]
  }
  J <- cbind(Bd * rho, Jr)
  p <- plogis(rho * delta)
  bread <- solve(crossprod(J, J * (p * (1 - p))))
  score <- rowsum(J * (y - p), judge)
  nc <- nrow(score)
  influence <- score %*% bread * sqrt(nc / (nc - 1))
  li <- (K - 1L) + seq_along(f$free)
  influence <- influence[, li, drop = FALSE]
  dimnames(influence) <- list(rownames(score), f$free)
  list(ref = f$ref, free = f$free, lrho = log(f$rho[f$free]),
       influence = influence, cov = crossprod(influence),
       converged = f$converged, rank_ok = f$rank_ok)
}

independent_reconciliation <- function(panels, blocks) {
  y <- numeric(0)
  pan <- ref <- character(0)
  positions <- vector("list", length(blocks))
  Vblock <- matrix(0, 0, 0)
  for (i in seq_along(blocks)) {
    z <- blocks[[i]]
    pos <- length(y) + seq_along(z$free)
    positions[[i]] <- pos
    y <- c(y, z$lrho[z$free])
    pan <- c(pan, z$free)
    ref <- c(ref, rep(z$ref, length(z$free)))
    Vnew <- matrix(0, nrow(Vblock) + length(pos),
                   ncol(Vblock) + length(pos))
    if (nrow(Vblock))
      Vnew[seq_len(nrow(Vblock)), seq_len(ncol(Vblock))] <- Vblock
    Vnew[pos, pos] <- z$cov[z$free, z$free, drop = FALSE]
    Vblock <- Vnew
  }
  anchor <- panels[1L]
  cols <- setdiff(panels, anchor)
  X <- matrix(0, length(y), length(cols))
  for (r in seq_along(y)) {
    cp <- match(pan[r], cols)
    cr <- match(ref[r], cols)
    if (!is.na(cp)) X[r, cp] <- X[r, cp] + 1
    if (!is.na(cr)) X[r, cr] <- X[r, cr] - 1
  }
  W <- solve(Vblock)
  L <- solve(crossprod(X, W %*% X), crossprod(X, W))
  bred <- drop(L %*% y)

  judges <- unique(unlist(lapply(blocks, function(z)
    rownames(z$influence)), use.names = FALSE))
  stacked <- matrix(0, length(judges), length(y),
                    dimnames = list(judges, NULL))
  for (i in seq_along(blocks)) {
    z <- blocks[[i]]
    rr <- match(rownames(z$influence), judges)
    stacked[rr, positions[[i]]] <- z$influence
  }
  Vjoint <- crossprod(stacked)
  Vred <- L %*% Vjoint %*% t(L)
  Vold_red <- L %*% Vblock %*% t(L)

  G <- length(panels)
  full <- old_full <- matrix(0, G, G)
  full[match(cols, panels), match(cols, panels)] <- Vred
  old_full[match(cols, panels), match(cols, panels)] <- Vold_red
  lphi <- setNames(numeric(G), panels)
  lphi[cols] <- bred
  centre <- diag(G) - matrix(1 / G, G, G)
  cov <- centre %*% full %*% t(centre)
  old_cov <- centre %*% old_full %*% t(centre)
  list(lphi = setNames(drop(centre %*% lphi), panels),
       cov = cov, old_cov = old_cov,
       cross_cov = if (length(blocks) == 2L)
         Vjoint[positions[[1L]], positions[[2L]]] else NA_real_)
}

run_one <- function(seed, allocation, dependence) {
  set.seed(seed)
  made <- tryCatch({
    d <- simulate_btl_efrm(
      n_objects_per_set = 5, n_sets = 2, n_judges_per_panel = 12,
      n_panels = 2, reps_within = 24, reps_cross = 2, seed = seed)
    truth <- attr(d, "truth")
    layout <- remap_judges(d, truth, allocation)
    d <- plant_outcomes(layout, truth, dependence)
    list(data = d, truth = truth)
  }, error = function(e) e)
  if (inherits(made, "error"))
    return(list(status = "error", reason = conditionMessage(made)))

  fit <- tryCatch(btl_efrm(
    made$data, "object_a", "object_b", winner = "winner", judge = "judge",
    panels = "panel", object_sets = made$truth$object_sets,
    se_method = "conditional"), error = function(e) e)
  if (inherits(fit, "error"))
    return(list(status = "refused", reason = conditionMessage(fit)))
  if (!isTRUE(fit$converged))
    return(list(status = "nonconv", reason = "public fit did not converge"))

  checked <- tryCatch({
    blocks <- lapply(made$truth$object_sets, function(objects)
      stage_block(made$data, objects))
    if (!all(vapply(blocks, function(z)
      isTRUE(z$converged) && isTRUE(z$rank_ok), logical(1))))
      return(list(status = "nonconv",
                  reason = "independent stage-one reconstruction did not converge"))
    expected <- independent_reconciliation(fit$panels, blocks)
    actual_lphi <- setNames(log(fit$phi_table$phi), fit$phi_table$panel)
    actual_se <- setNames(fit$phi_table$se_log_phi, fit$phi_table$panel)
    expected_se <- setNames(sqrt(pmax(diag(expected$cov), 0)), fit$panels)
    old_se <- setNames(sqrt(pmax(diag(expected$old_cov), 0)), fit$panels)
    list(status = "ok",
         point_error = max(abs(actual_lphi[fit$panels] - expected$lphi)),
         covariance_error = max(abs(actual_se[fit$panels]^2 -
                                      diag(expected$cov))),
         reported_to_stacked = mean(actual_se[fit$panels] / expected_se),
         old_to_corrected = mean(old_se / expected_se),
         cross_cov = unname(expected$cross_cov[1L]))
  }, error = function(e) list(status = "error",
                              reason = conditionMessage(e)))
  checked
}

scenarios <- list(
  list(label = "shared judges; positive dependence",
       allocation = "shared", dependence = 1),
  list(label = "shared judges; negative dependence",
       allocation = "shared", dependence = -1),
  list(label = "partially shared judges; positive dependence",
       allocation = "partial", dependence = 1),
  list(label = "partially shared judges; negative dependence",
       allocation = "partial", dependence = -1),
  list(label = "disjoint judges; no cross-set dependence",
       allocation = "disjoint", dependence = 0)
)

rows <- lapply(seq_along(scenarios), function(i) {
  sc <- scenarios[[i]]
  ans <- lapply(seq_len(NREP), function(r)
    run_one(741000L + 1000L * i + r, sc$allocation, sc$dependence))
  status <- vapply(ans, function(z) z[["status"]], "")
  ok <- ans[status == "ok"]
  take <- function(name)
    if (length(ok)) vapply(ok, function(z) z[[name]], 0) else numeric(0)
  point_error <- take("point_error")
  covariance_error <- take("covariance_error")
  reported_ratio <- take("reported_to_stacked")
  old_ratio <- take("old_to_corrected")
  cross_cov <- take("cross_cov")
  reasons <- unique(vapply(ans[status != "ok"], function(z)
    z$reason %||% "", ""))

  row <- sv_row(
    STUDY, sc$label,
    "same-estimator conditional panel covariance agreement",
    n_reps = length(ok), effect = sc$dependence,
    n_attempted = NREP, n_refused = sum(status == "refused"),
    n_nonconv = sum(status == "nonconv"), n_error = sum(status == "error"),
    n_withheld = 0L, n_metric_unavailable = 0L,
    notes = paste0(
      NREP, " fixed seeds in this run; actual btl_efrm conditional panel variances ",
      "versus independently rebuilt CR1 stage-one judge scores and their ",
      "same-estimator stacked sandwich. Old/corrected SE ratios are descriptive ",
      "comparisons with the former block-diagonal covariance, not a null-",
      "calibration or coverage claim",
      if (length(reasons)) paste0("; first failure: ", reasons[1L]) else ""))
  row$max_absolute_log_phi_error <- if (length(point_error))
    max(point_error) else NA_real_
  row$max_absolute_covariance_error <- if (length(covariance_error))
    max(covariance_error) else NA_real_
  row$mean_reported_to_stacked_se_ratio <- if (length(reported_ratio))
    mean(reported_ratio) else NA_real_
  row$max_absolute_reported_to_stacked_se_ratio_minus_one <-
    if (length(reported_ratio)) max(abs(reported_ratio - 1)) else NA_real_
  row$mean_old_to_corrected_se_ratio <- if (length(old_ratio))
    mean(old_ratio) else NA_real_
  row$min_old_to_corrected_se_ratio <- if (length(old_ratio))
    min(old_ratio) else NA_real_
  row$max_old_to_corrected_se_ratio <- if (length(old_ratio))
    max(old_ratio) else NA_real_
  row$mean_stacked_cross_set_covariance <- if (length(cross_cov))
    mean(cross_cov) else NA_real_
  row
})

rows <- do.call(rbind, rows)
positive <- grepl("positive dependence", rows$scenario, fixed = TRUE)
negative <- grepl("negative dependence", rows$scenario, fixed = TRUE)
disjoint <- grepl("disjoint judges", rows$scenario, fixed = TRUE)
stopifnot(rows$n_reps + rows$n_refused + rows$n_nonconv + rows$n_error ==
            rows$n_attempted,
          all(rows$max_absolute_log_phi_error[rows$n_reps > 0] < 1e-8),
          all(rows$max_absolute_covariance_error[rows$n_reps > 0] < 1e-8),
          all(rows$mean_stacked_cross_set_covariance[positive] > 0),
          all(rows$mean_stacked_cross_set_covariance[negative] < 0),
          all(abs(rows$mean_stacked_cross_set_covariance[disjoint]) < 1e-12))
sv_write(rows, STUDY)
