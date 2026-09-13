suppressWarnings(pkgload::load_all(".", quiet=TRUE))
t0 <- Sys.time()

configs <- list(
  null        = list(set_unit_ratio = 1.0, group_unit_ratio = 1.0, mcar = 0),
  set1.3      = list(set_unit_ratio = 1.3, group_unit_ratio = 1.0, mcar = 0),
  group1.3    = list(set_unit_ratio = 1.0, group_unit_ratio = 1.3, mcar = 0),
  both1.5     = list(set_unit_ratio = 1.5, group_unit_ratio = 1.5, mcar = 0),
  set1.3_mcar = list(set_unit_ratio = 1.3, group_unit_ratio = 1.0, mcar = 0.25)
)

n_reps <- 6
n_per_group <- 300
items_per_set <- 8
n_sets <- 3
n_groups <- 3

apply_mcar <- function(d, p, seed) {
  set.seed(seed + 99999L)
  item_cols <- setdiff(names(d), c("id", "group"))
  X <- as.matrix(d[item_cols])
  mask <- matrix(runif(length(X)) < p, nrow(X), ncol(X))
  X[mask] <- NA
  d[item_cols] <- as.data.frame(X)
  d
}

# custom recovery of item location (common-unit) and person theta, since
# sim_recovery()'s efrm branch only covers alpha (set units) and phi (group
# units) -- item/person recovery done by hand here, mean-centred as an
# origin-identified location parameter (matches sim_recovery's own convention
# for the "rasch" layout)
item_person_recovery <- function(fit, tr) {
  ei <- setNames(fit$item_arbitrary$location, fit$item_arbitrary$item)
  cm <- intersect(names(tr$difficulty), names(ei))
  ti <- tr$difficulty[cm] - mean(tr$difficulty[cm])
  ei <- ei[cm] - mean(ei[cm])
  item_r <- c(cor = cor(ti, ei), rmse = sqrt(mean((ei - ti)^2)))

  th <- fit$person$theta
  keep <- is.finite(th)
  tt <- tr$theta[keep] - mean(tr$theta[keep])
  ee <- th[keep] - mean(th[keep])
  pers_r <- c(cor = cor(tt, ee), rmse = sqrt(mean((ee - tt)^2)))
  list(item = item_r, person = pers_r)
}

results <- list()
for (cfg_name in names(configs)) {
  cfg <- configs[[cfg_name]]
  rows <- list()
  for (r in seq_len(n_reps)) {
    seed <- 1000L * which(names(configs) == cfg_name) + r
    d <- simulate_efrm(n_per_group, items_per_set, n_sets = n_sets, n_groups = n_groups,
                        set_unit_ratio = cfg$set_unit_ratio,
                        group_unit_ratio = cfg$group_unit_ratio,
                        seed = seed)
    tr <- attr(d, "truth")
    if (cfg$mcar > 0) d <- apply_mcar(d, cfg$mcar, seed)
    fit <- tryCatch(rasch_efrm(d, id = "id", item_sets = tr$item_sets, groups = "group",
                                se_method = "hybrid"),
                     error = function(e) NULL)
    if (is.null(fit)) next
    rec <- sim_recovery(fit, d)
    ip <- item_person_recovery(fit, tr)
    rows[[length(rows) + 1L]] <- list(
      alpha_cor = rec$summary$correlation[rec$summary$parameter == "set unit (log)"],
      alpha_rmse = rec$summary$rmse[rec$summary$parameter == "set unit (log)"],
      phi_cor = rec$summary$correlation[rec$summary$parameter == "group unit (log)"],
      phi_rmse = rec$summary$rmse[rec$summary$parameter == "group unit (log)"],
      item_cor = ip$item["cor"], item_rmse = ip$item["rmse"],
      pers_cor = ip$person["cor"], pers_rmse = ip$person["rmse"])
  }
  df <- do.call(rbind, lapply(rows, as.data.frame))
  results[[cfg_name]] <- df
  se <- function(x) sd(x, na.rm=TRUE) / sqrt(sum(is.finite(x)))
  cat(sprintf(
    "[%s] n_ok=%d\n  alpha_cor=%.3f(mcse %.3f) alpha_rmse=%.3f\n  phi_cor=%.3f(mcse %.3f) phi_rmse=%.3f\n  item_cor=%.3f item_rmse=%.3f\n  pers_cor=%.3f pers_rmse=%.3f\n",
    cfg_name, nrow(df),
    mean(df$alpha_cor, na.rm=TRUE), se(df$alpha_cor), mean(df$alpha_rmse, na.rm=TRUE),
    mean(df$phi_cor, na.rm=TRUE), se(df$phi_cor), mean(df$phi_rmse, na.rm=TRUE),
    mean(df$item_cor, na.rm=TRUE), mean(df$item_rmse, na.rm=TRUE),
    mean(df$pers_cor, na.rm=TRUE), mean(df$pers_rmse, na.rm=TRUE)))
}

saveRDS(results, "tools/simval/round1/efrm/check1_results.rds")
cat("elapsed:", as.numeric(Sys.time() - t0, units = "secs"), "s\n")
