simP <- function(theta, tau) { x <- 0:length(tau); p <- exp(x * theta - c(0, cumsum(tau))); p / sum(p) }

test_that("dichotomous item difficulties recover by pairwise conditional ML", {
  set.seed(7); Np <- 1500; L <- 20
  dtrue <- scale(seq(-2.5, 2.5, length.out = L), scale = FALSE)[, 1]
  th <- rnorm(Np, 0, 1.4)
  X <- matrix(rbinom(Np * L, 1, plogis(outer(th, dtrue, "-"))), Np, L)
  colnames(X) <- sprintf("I%02d", 1:L)

  est <- pcml(X)
  expect_true(est$converged)
  expect_gt(cor(est$thr$tau, dtrue), 0.99)
  expect_lt(sqrt(mean((est$thr$tau - dtrue)^2)), 0.12)
  expect_true(all(est$thr$se > 0 & est$thr$se < 0.5))
})

test_that("sandwich standard errors track the sampling variability", {
  set.seed(42); L <- 10; Np <- 500
  dtrue <- scale(seq(-2, 2, length.out = L), scale = FALSE)[, 1]
  reps <- 12
  est <- ses <- matrix(NA, reps, L)
  for (r in 1:reps) {
    th <- rnorm(Np, 0, 1.3)
    X <- matrix(rbinom(Np * L, 1, plogis(outer(th, dtrue, "-"))), Np, L)
    f <- pcml(X)
    est[r, ] <- f$thr$tau; ses[r, ] <- f$thr$se
  }
  ratio <- mean(colMeans(ses)) / mean(apply(est, 2, sd))
  expect_gt(ratio, 0.6); expect_lt(ratio, 1.6)
})

test_that("polytomous PCM thresholds recover", {
  set.seed(99); Np <- 2500
  mvec <- rep(c(2, 3), length.out = 15)
  tt <- lapply(mvec, function(m) sort(rnorm(m, 0, 0.9)) + rnorm(1, 0, 1.1))
  tt <- lapply(tt, function(t) t - mean(sapply(tt, mean)))
  th <- rnorm(Np, 0, 1.5)
  X <- sapply(seq_along(mvec), function(i)
    sapply(th, function(t) sample(0:mvec[i], 1, prob = simP(t, tt[[i]]))))
  colnames(X) <- sprintf("P%02d", seq_along(mvec))

  fit <- rasch(X, model = "PCM")
  expect_true(fit$est$converged)
  expect_gt(cor(fit$thresholds$tau, unlist(tt)), 0.99)
})

test_that("constrained RSM recovers locations and steps", {
  set.seed(123); Np <- 2500; L <- 12; m <- 4
  step_true <- c(-1.3, -0.4, 0.5, 1.2); step_true <- step_true - mean(step_true)
  loc_true <- scale(seq(-1.8, 1.8, length.out = L), scale = FALSE)[, 1]
  th <- rnorm(Np, 0, 1.4)
  X <- sapply(seq_len(L), function(i)
    sapply(th, function(t) sample(0:m, 1, prob = simP(t, loc_true[i] + step_true))))
  colnames(X) <- sprintf("R%02d", seq_len(L))

  fit <- rasch(X, model = "RSM")
  loc_est <- fit$items$location
  step_est <- vapply(1:m, function(k)
    mean(fit$thresholds$tau[fit$thresholds$k == k] -
         loc_est[fit$thresholds$item[fit$thresholds$k == k]]), 0)
  expect_gt(cor(loc_est, loc_true), 0.99)
  expect_lt(max(abs(step_est - step_true)), 0.12)
})

test_that("Warm WLE is finite at extreme scores and symmetric", {
  tau_list <- as.list(scale(seq(-2, 2, length.out = 20), scale = FALSE)[, 1])
  pe <- person_wle(tau_list)
  expect_true(is.finite(pe$theta["0"]))
  expect_true(is.finite(pe$theta["20"]))
  expect_equal(unname(pe$theta["0"]), -unname(pe$theta["20"]), tolerance = 1e-3)
  expect_true(all(is.finite(pe$se[c("0", "1", "19", "20")])))
})

test_that("missing data: persons estimated on their observed items", {
  set.seed(5); Np <- 600; L <- 10
  dtrue <- scale(seq(-2, 2, length.out = L), scale = FALSE)[, 1]
  th <- rnorm(Np, 0, 1.3)
  X <- matrix(rbinom(Np * L, 1, plogis(outer(th, dtrue, "-"))), Np, L)
  colnames(X) <- sprintf("I%02d", 1:L)
  X[sample(length(X), 600)] <- NA

  fit <- rasch(X)
  ok <- fit$person$n_items > 0
  expect_true(all(!is.na(fit$person$theta[ok])))
  expect_true(all(fit$person$max_raw[ok] <= L))
  expect_gt(cor(fit$person$theta, th, use = "complete.obs"), 0.6)
})

test_that("score table is monotone in the raw score", {
  set.seed(2); Np <- 800; L <- 12
  dtrue <- scale(seq(-2, 2, length.out = L), scale = FALSE)[, 1]
  X <- matrix(rbinom(Np * L, 1, plogis(outer(rnorm(Np, 0, 1.3), dtrue, "-"))), Np, L)
  colnames(X) <- sprintf("I%02d", 1:L)
  sc <- score_table(rasch(X))
  expect_equal(nrow(sc), L + 1)
  expect_true(all(diff(sc$theta) > 0))
  expect_true(all(sc$se > 0))
})
