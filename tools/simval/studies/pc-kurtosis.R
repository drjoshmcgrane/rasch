# Published-polynomial checks and fixed-truth PC kurtosis recovery.
# Run from the package root with Rscript tools/simval/studies/pc-kurtosis.R.
options(rasch.max_workers = 1L, device = function(...) grDevices::pdf(NULL))
pkgload::load_all(".", quiet = TRUE)
source("tools/simval/harness.R")

# Independent adjacent-threshold expression: Linacre, Andrich & Luo (2003),
# RMT 17:3, p. 944. Production code differences cumulative polynomials.
published <- function(m) {
  z <- seq_len(m) - (m + 1) / 2
  cbind(1, 2 * z, 2 * (3 * z^2 - (m^2 - 1) / 4),
        5 * (4 * z^3 - z * (3 * m^2 - 7) / 5))
}
reps <- 100L
rows <- list()
for (m in 4:6) {
  G <- published(m)
  stopifnot(max(abs(.pc_gcoefs(m) - G)) < 1e-10, qr(G)$rank == 4L)
  coefficients <- rbind(seq(-1, 1, length.out = 6), rep(.35, 6),
                         rep(c(-.01, .01), 3), rep(c(-.012, .006, .018), 2))
  tau <- G %*% coefficients
  estimate <- se <- matrix(NA_real_, reps, 6)
  errors <- nonconv <- 0L
  max_tau_error <- max_cov_error <- max_loglik_error <- 0
  for (r in seq_len(reps)) {
    set.seed(96200 + 1000 * m + r)
    theta <- rnorm(1200)
    X <- vapply(seq_len(6), function(i) vapply(theta, function(th)
      sample(0:m, 1L, prob = exp((0:m) * th - c(0, cumsum(tau[, i])))), 0L),
      integer(length(theta)))
    colnames(X) <- paste0("I", 1:6)
    fit <- tryCatch(pcml_pc(X, n_components = 4), error = identity)
    if (inherits(fit, "error")) { errors <- errors + 1L; next }
    if (!isTRUE(fit$converged)) { nonconv <- nonconv + 1L; next }
    estimate[r, ] <- fit$components$kurtosis
    se[r, ] <- fit$components$kurtosis_se
    if (m == 4L) {
      free <- pcml(X)
      stopifnot(isTRUE(free$converged), fit$n_parameters == free$n_parameters)
      max_tau_error <- max(max_tau_error, max(abs(fit$thr$tau - free$thr$tau)))
      max_cov_error <- max(max_cov_error, max(abs(fit$cov_tau - free$cov_tau)))
      max_loglik_error <- max(max_loglik_error, abs(fit$loglik - free$loglik))
    }
    if (r %% 25 == 0) cat("thresholds", m, "replicate", r, "\n")
  }
  ok <- apply(is.finite(estimate) & is.finite(se) & se > 0, 1, all)
  err <- sweep(estimate[ok, , drop = FALSE], 2, coefficients[4, ], "-")
  cover <- rowMeans(abs(err) <= qnorm(.975) * se[ok, , drop = FALSE])
  # Bias/SD/mean SE describe I3, whose truth is fixed at .018. Coverage is
  # averaged over items with Monte Carlo uncertainty at the replicate level.
  row <- sv_row("pc-kurtosis", paste("1200 persons, 6 items,", m, "thresholds"),
    "I3 kurtosis recovery; mean item coverage", sum(ok),
    bias = mean(err[, 3]), emp_sd = sd(estimate[ok, 3]), mean_se = mean(se[ok, 3]),
    coverage95 = mean(cover), mc_override = list(coverage95 = sd(cover) / sqrt(sum(ok))),
    effect = coefficients[4, 3], n_attempted = reps, n_refused = 0L,
    n_nonconv = nonconv, n_error = errors, n_withheld = 0L,
    n_metric_unavailable = sum(!ok) - nonconv - errors,
    notes = paste("Fixed truth generated using independent published adjacent polynomials;",
      "100 replicates per condition, not a principal null-size study.",
      if (m == 4L) paste("Maximum free-PCM discrepancies: thresholds", max_tau_error,
        "covariance", max_cov_error, "log likelihood", max_loglik_error) else
        "Five/six thresholds use the correctly restricted four-component model."))
  rows[[length(rows) + 1L]] <- row
  print(row[, c("scenario", "bias", "se_ratio", "coverage95", "mc_se_coverage")])
}
sv_write(do.call(rbind, rows), "pc-kurtosis")
