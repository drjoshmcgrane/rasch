simP <- function(theta, tau) { x <- 0:length(tau); p <- exp(x * theta - c(0, cumsum(tau))); p / sum(p) }

test_that("pc reparameterisation matches pcml exactly through four thresholds", {
  set.seed(99); Np <- 1500
  mvec <- rep(c(1, 2, 3, 4), length.out = 12)
  tt <- lapply(mvec, function(m) sort(rnorm(m, 0, 0.9)) + rnorm(1, 0, 1.1))
  tt <- lapply(tt, function(t) t - mean(sapply(tt, mean)))
  th <- rnorm(Np, 0, 1.5)
  X <- sapply(seq_along(mvec), function(i)
    sapply(th, function(t) sample(0:mvec[i], 1, prob = simP(t, tt[[i]]))))
  colnames(X) <- sprintf("P%02d", seq_along(mvec))

  free <- pcml(X, model = "PCM")
  pc <- pcml_pc(X, n_components = 4)

  expect_equal(pc$loglik, free$loglik, tolerance = 1e-6)
  expect_equal(pc$thr$tau, free$thr$tau, tolerance = 1e-6)
  expect_equal(pc$cov_tau, free$cov_tau, tolerance = 1e-6)
  expect_equal(pc$n_parameters, free$n_parameters)
  expect_true(all(is.na(pc$components$spread[mvec == 1])))
  expect_true(all(!is.na(pc$components$spread[mvec >= 2])))
  expect_true(all(is.na(pc$components$skewness[mvec < 3])))
})

test_that("kurtosis is identified at four thresholds and longer scales are reduced", {
  set.seed(3); Np <- 3000

  mk <- function(mvec) {
    tt <- lapply(mvec, function(m) sort(rnorm(m, 0, 0.7)) + rnorm(1, 0, 1.1))
    tt <- lapply(tt, function(t) t - mean(sapply(tt, mean)))
    th <- rnorm(Np, 0, 1.5)
    X <- sapply(seq_along(mvec), function(i)
      sapply(th, function(t) sample(0:mvec[i], 1, prob = simP(t, tt[[i]]))))
    colnames(X) <- sprintf("P%02d", seq_along(mvec))
    X
  }

  # Four thresholds have four independent polynomial components.
  X4 <- mk(rep(4, 6))
  pc4 <- pcml_pc(X4, n_components = 4)
  free4 <- pcml(X4, model = "PCM")
  expect_true(all(is.finite(pc4$components$kurtosis)))
  expect_true(all(is.finite(pc4$components$kurtosis_se)))
  expect_equal(pc4$n_parameters, free4$n_parameters)
  expect_equal(pc4$loglik, free4$loglik, tolerance = 1e-6)
  expect_equal(pc4$thr$tau, free4$thr$tau, tolerance = 1e-6)

  # 5 thresholds: kurtosis now identified, but the family caps at 4
  # components while 5 thresholds need 5 for an exact free-threshold match
  X5 <- mk(rep(5, 6))
  pc5 <- pcml_pc(X5, n_components = 4)
  free5 <- pcml(X5, model = "PCM")
  expect_true(all(!is.na(pc5$components$kurtosis)))
  expect_lt(pc5$n_parameters, free5$n_parameters)
  expect_lt(pc5$loglik, free5$loglik + 1e-8)
})

test_that("PC coefficients match the published threshold polynomials", {
  # Linacre, Andrich & Luo (2003), RMT 17:3, p. 944: an independent
  # adjacent-threshold expression, not a copy of .pc_gcoefs' cumulative one.
  for (m in 1:14) {
    z <- seq_len(m) - (m + 1) / 2
    published <- cbind(1, 2 * z,
      2 * (3 * z^2 - (m^2 - 1) / 4),
      5 * (4 * z^3 - z * (3 * m^2 - 7) / 5))
    G <- .pc_gcoefs(m)
    nc <- min(m, 4L)
    expect_equal(G[, seq_len(nc), drop = FALSE],
                 published[, seq_len(nc), drop = FALSE], ignore_attr = TRUE)
    expect_equal(qr(G[, seq_len(nc), drop = FALSE])$rank, nc)
  }
  expect_equal(.pc_gcoefs(4)[, 4], c(-6, 18, -18, 6))
})

test_that("location-only pc model recovers item locations and is a valid reduction", {
  set.seed(11); Np <- 2000
  mvec <- rep(3, 8)
  tt <- lapply(mvec, function(m) sort(rnorm(m, 0, 0.3)) + rnorm(1, 0, 1.2))
  tt <- lapply(tt, function(t) t - mean(sapply(tt, mean)))
  loc_true <- vapply(tt, mean, 0)
  th <- rnorm(Np, 0, 1.5)
  X <- sapply(seq_along(mvec), function(i)
    sapply(th, function(t) sample(0:mvec[i], 1, prob = simP(t, tt[[i]]))))
  colnames(X) <- sprintf("Q%02d", seq_along(mvec))

  full <- pcml_pc(X, n_components = 4)
  loc_only <- pcml_pc(X, n_components = 1)

  expect_gt(cor(loc_only$components$location, loc_true), 0.98)
  expect_lt(loc_only$loglik, full$loglik + 1e-8)
  expect_true(all(diff(loc_only$thr$tau[loc_only$thr$item == 1]) == 0))
  expect_equal(loc_only$n_parameters, length(mvec) - 1L)
})

test_that("pcml_pc handles missing data and a purely dichotomous set", {
  set.seed(5); Np <- 800; L <- 8
  dtrue <- scale(seq(-2, 2, length.out = L), scale = FALSE)[, 1]
  th <- rnorm(Np, 0, 1.3)
  X <- matrix(rbinom(Np * L, 1, plogis(outer(th, dtrue, "-"))), Np, L)
  X[sample(length(X), floor(0.1 * length(X)))] <- NA
  colnames(X) <- sprintf("I%02d", 1:L)

  est <- pcml_pc(X)
  expect_true(est$converged)
  expect_gt(cor(est$thr$tau, dtrue), 0.95)
  expect_true(all(is.na(est$components$spread)))
})
