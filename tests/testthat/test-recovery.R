test_that("dichotomous item difficulties recover (eigenvector and pairwise-LS)", {
  set.seed(7); Np <- 1500; L <- 20
  dtrue <- scale(seq(-2.5, 2.5, length.out = L), scale = FALSE)[, 1]
  th <- rnorm(Np, 0, 1.4)
  X <- matrix(rbinom(Np * L, 1, plogis(outer(th, dtrue, "-"))), Np, L)
  colnames(X) <- sprintf("I%02d", 1:L)

  de <- est_eigen_dich(X)
  fit <- rasch(X, model = "PCM", solver = "LS")
  dl <- fit$items$location

  expect_gt(cor(de, dtrue), 0.99)
  expect_gt(cor(dl, dtrue), 0.99)
  expect_lt(sqrt(mean((dl - dtrue)^2)), 0.15)
  expect_lt(max(abs(de - dl)), 0.10)          # the two solvers agree
})

test_that("polytomous PCM thresholds recover and LS == RA", {
  set.seed(99); Np <- 2500
  simP <- function(theta, tau) { x <- 0:length(tau); p <- exp(x * theta - c(0, cumsum(tau))); p / sum(p) }
  mvec <- rep(c(2, 3), length.out = 15)
  tt <- lapply(mvec, function(m) sort(rnorm(m, 0, 0.9)) + rnorm(1, 0, 1.1))
  tt <- lapply(tt, function(t) t - mean(unlist(tt)))
  th <- rnorm(Np, 0, 1.5)
  X <- sapply(seq_along(mvec), function(i) sapply(th, function(t) sample(0:mvec[i], 1, prob = simP(t, tt[[i]]))))
  colnames(X) <- sprintf("P%02d", seq_along(mvec))

  fLS <- rasch(X, model = "PCM", solver = "LS")
  fRA <- rasch(X, model = "PCM", solver = "RA")
  expect_gt(cor(fLS$thresholds$tau, unlist(tt)), 0.99)
  expect_lt(max(abs(fLS$thresholds$tau - fRA$thresholds$tau)), 1e-3)
})

test_that("Warm WLE is finite at extreme scores and symmetric", {
  tau_list <- as.list(scale(seq(-2, 2, length.out = 20), scale = FALSE)[, 1])
  pe <- person_wle(tau_list)
  expect_true(is.finite(pe$theta["0"]))
  expect_true(is.finite(pe$theta["20"]))
  expect_equal(unname(pe$theta["0"]), -unname(pe$theta["20"]), tolerance = 1e-3)
  expect_true(all(is.finite(pe$se[c("0", "1", "19", "20")])))
})

test_that("constrained RSM recovers locations and steps without shrinkage", {
  set.seed(123); Np <- 3000; L <- 12; m <- 4
  simP <- function(theta, tau) { x <- 0:length(tau); p <- exp(x * theta - c(0, cumsum(tau))); p / sum(p) }
  step_true <- c(-1.3, -0.4, 0.5, 1.2); step_true <- step_true - mean(step_true)
  loc_true <- scale(seq(-1.8, 1.8, length.out = L), scale = FALSE)[, 1]
  tau2 <- lapply(loc_true, function(b) b + step_true)
  th <- rnorm(Np, 0, 1.4)
  X <- sapply(seq_len(L), function(i) sapply(th, function(t) sample(0:m, 1, prob = simP(t, tau2[[i]]))))
  colnames(X) <- sprintf("R%02d", seq_len(L))

  rs <- rasch_rsm(X)
  expect_gt(cor(rs$location, loc_true), 0.99)
  expect_lt(max(abs(rs$step - step_true)), 0.10)   # shrinkage controlled
})
