simP <- function(theta, tau) { x <- 0:length(tau); p <- exp(x * theta - c(0, cumsum(tau))); p / sum(p) }

test_that("pc reparameterisation matches pcml exactly for items with at most 3 thresholds", {
  set.seed(99); Np <- 1500
  mvec <- rep(c(1, 2, 3), length.out = 12)
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
  expect_equal(pc$n_parameters, free$n_parameters)
  expect_true(all(is.na(pc$components$spread[mvec == 1])))
  expect_true(all(!is.na(pc$components$spread[mvec >= 2])))
  expect_true(all(is.na(pc$components$skewness[mvec < 3])))
})

test_that("kurtosis is identified at 5 thresholds but not at exactly 4, and either way pcml_pc is a genuine reduction once an item exceeds 3 thresholds", {
  # Andrich/Pedler's family stops at the quartic (kurtosis) component, so full
  # rank equivalence with the free-threshold pcml() only holds while every
  # item has at most 3 thresholds (location + spread + skewness span it
  # exactly). From 4 thresholds on, pcml_pc is necessarily a smoothed,
  # reduced-rank reparameterisation, however large n_components is set.
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

  # exactly 4 thresholds: quartic collinear with cubic, only 3 components ever
  X4 <- mk(rep(4, 6))
  pc4 <- pcml_pc(X4, n_components = 4)
  free4 <- pcml(X4, model = "PCM")
  expect_true(all(is.na(pc4$components$kurtosis)))
  expect_equal(pc4$n_parameters, free4$n_parameters - length(unique(pc4$thr$item)))
  expect_lt(pc4$loglik, free4$loglik + 1e-8)

  # 5 thresholds: kurtosis now identified, but the family caps at 4
  # components while 5 thresholds need 5 for an exact free-threshold match
  X5 <- mk(rep(5, 6))
  pc5 <- pcml_pc(X5, n_components = 4)
  free5 <- pcml(X5, model = "PCM")
  expect_true(all(!is.na(pc5$components$kurtosis)))
  expect_lt(pc5$n_parameters, free5$n_parameters)
  expect_lt(pc5$loglik, free5$loglik + 1e-8)
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
