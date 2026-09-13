suppressWarnings(pkgload::load_all(".", quiet=TRUE))
t0 <- Sys.time()

apply_mcar <- function(d, p, seed) {
  set.seed(seed + 99999L)
  item_cols <- setdiff(names(d), c("id", "group"))
  X <- as.matrix(d[item_cols])
  mask <- matrix(runif(length(X)) < p, nrow(X), ncol(X))
  X[mask] <- NA
  d[item_cols] <- as.data.frame(X)
  d
}

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

n_reps <- 8
n_per_group <- 400
items_per_set <- 10

rows <- list(); n_refused <- 0L
for (r in seq_len(n_reps)) {
  seed <- 90000L + r
  d <- simulate_efrm(n_per_group, items_per_set, n_sets = 2, n_groups = 2,
                      set_unit_ratio = 1.3, group_unit_ratio = 1.2, seed = seed)
  tr <- attr(d, "truth")
  d <- apply_mcar(d, 0.25, seed)
  fit <- tryCatch(rasch_efrm(d, id = "id", item_sets = tr$item_sets, groups = "group", se_method = "hybrid"),
                   error = function(e) { cat("  refused (rep", r, "):", conditionMessage(e), "\n"); NULL })
  if (is.null(fit)) { n_refused <- n_refused + 1L; next }
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
cat(sprintf("n_ok=%d  n_refused=%d\n", nrow(df), n_refused))
print(round(colMeans(df, na.rm = TRUE), 4))
saveRDS(list(df = df, n_refused = n_refused),
        "tools/simval/round1/efrm/check1b_results.rds")
cat("elapsed:", as.numeric(Sys.time() - t0, units = "secs"), "s\n")
