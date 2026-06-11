simP <- function(theta, tau) { x <- 0:length(tau); p <- exp(x * theta - c(0, cumsum(tau))); p / sum(p) }

test_that("dimensionality test separates 1D from 2D data", {
  set.seed(1); Np <- 1500; L <- 20
  d <- scale(seq(-2, 2, length.out = L), scale = FALSE)[, 1]; th <- rnorm(Np, 0, 1.4)
  X1 <- matrix(rbinom(Np * L, 1, plogis(outer(th, d, "-"))), Np, L); colnames(X1) <- sprintf("U%02d", 1:L)
  dt1 <- dimensionality_test(rasch(X1, model = "PCM"))
  expect_false(dt1$multidimensional)

  set.seed(2)
  thA <- rnorm(Np, 0, 1.4); thB <- 0.3 * thA + sqrt(1 - 0.3^2) * rnorm(Np, 0, 1.4)
  XA <- matrix(rbinom(Np * 10, 1, plogis(outer(thA, d[1:10], "-"))), Np, 10)
  XB <- matrix(rbinom(Np * 10, 1, plogis(outer(thB, d[11:20], "-"))), Np, 10)
  X2 <- cbind(XA, XB); colnames(X2) <- sprintf("D%02d", 1:20)
  dt2 <- dimensionality_test(rasch(X2, model = "PCM"))
  expect_true(dt2$multidimensional)
  expect_gt(dt2$prop_significant, dt1$prop_significant)
})

test_that("local dependence is flagged for a near-duplicated item", {
  set.seed(1); Np <- 1500; L <- 20
  d <- scale(seq(-2, 2, length.out = L), scale = FALSE)[, 1]; th <- rnorm(Np, 0, 1.4)
  X <- matrix(rbinom(Np * L, 1, plogis(outer(th, d, "-"))), Np, L); colnames(X) <- sprintf("U%02d", 1:L)
  set.seed(3); X[, 5] <- ifelse(runif(Np) < 0.85, X[, 4], X[, 5])
  fl <- residual_correlations(rasch(X, model = "PCM"), flag = 0.2)$flagged
  expect_true(any((fl$item_a == "U04" & fl$item_b == "U05") |
                  (fl$item_a == "U05" & fl$item_b == "U04")))
})

test_that("uniform DIF is detected on planted items only", {
  set.seed(4); Np <- 2000; L <- 15
  d <- scale(seq(-2, 2, length.out = L), scale = FALSE)[, 1]
  grp <- rep(c("ref", "foc"), each = Np / 2); th <- rnorm(Np, 0, 1.4)
  shift <- matrix(0, Np, L); shift[grp == "foc", 3] <- 1.0; shift[grp == "foc", 10] <- -1.0
  X <- matrix(rbinom(Np * L, 1, plogis(outer(th, d, "-") - shift)), Np, L); colnames(X) <- sprintf("G%02d", 1:L)
  dif <- dif_anova(rasch(X, model = "PCM"), group = grp, n_groups = 6)
  expect_true(dif$uniform_DIF[3]); expect_true(dif$uniform_DIF[10])
  expect_equal(sum(dif$uniform_DIF[-c(3, 10)]), 0)
})

test_that("threshold disordering is detected", {
  set.seed(5); Np <- 2000; th <- rnorm(Np, 0, 1.4)
  tau_dis <- c(1.5, -1.5, 0.8); tau_ok <- c(-1.0, 0.0, 1.0)
  X <- cbind(sapply(th, function(t) sample(0:3, 1, prob = simP(t, tau_dis))),
             sapply(th, function(t) sample(0:3, 1, prob = simP(t, tau_ok))),
             sapply(th, function(t) sample(0:3, 1, prob = simP(t, tau_ok + 0.4))))
  colnames(X) <- c("DIS", "ok1", "ok2")
  td <- rasch(X, model = "PCM", n_groups = 6)$thresholds_diag
  expect_false(td[["1"]]$ordered)
  expect_true(td[["2"]]$ordered)
})
