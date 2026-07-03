# Bradley-Terry-Luce paired comparisons: equivalence with the conditional
# Rasch model, parameter recovery, judge diagnostics, and guards.

sim_btl <- function(beta, n_per_pair, seed = 1, judges = NULL) {
  set.seed(seed)
  pr <- t(utils::combn(names(beta), 2))
  d <- data.frame(a = rep(pr[, 1], each = n_per_pair),
                  b = rep(pr[, 2], each = n_per_pair))
  p <- plogis(beta[d$a] - beta[d$b])
  d$win <- ifelse(runif(nrow(d)) < p, d$a, d$b)
  if (!is.null(judges)) d$judge <- sample(judges, nrow(d), replace = TRUE)
  d
}

test_that("BTL is the conditional form of the dichotomous Rasch model", {
  # extract from Rasch data every pair response with exactly one correct:
  # 'the correct item beat the incorrect one'. The pairwise conditional
  # likelihood pcml() maximises is exactly the BTL likelihood of those
  # comparisons, so the estimates must agree up to solver tolerance
  # (with beta = -delta: winning means being easier).
  set.seed(4)
  d0 <- seq(-1.5, 1.5, length.out = 6)
  N <- 350
  X <- matrix(rbinom(N * 6, 1, plogis(outer(rnorm(N), d0, "-"))), N, 6)
  colnames(X) <- paste0("I", 1:6)
  cmp <- list()
  for (i in 1:5) for (j in (i + 1):6) {
    sel <- which(X[, i] + X[, j] == 1)
    if (!length(sel)) next
    cmp[[length(cmp) + 1L]] <- data.frame(
      person = sel,
      a = colnames(X)[i], b = colnames(X)[j],
      win = ifelse(X[sel, i] == 1, colnames(X)[i], colnames(X)[j]))
  }
  cmp <- do.call(rbind, cmp)
  bt <- btl(cmp, "a", "b", "win", judge = "person")
  pc <- pcml(X)
  delta <- vapply(1:6, function(i) mean(pc$thr$tau[pc$thr$item == i]), 0)
  loc <- bt$objects$location[match(colnames(X), bt$objects$object)]
  expect_equal(unname(loc), unname(-delta), tolerance = 1e-4)
  # person-clustered sandwich SEs agree with pcml's to a similar order
  se_pc <- vapply(1:6, function(i) {
    rows <- pc$thr$id[pc$thr$item == i]
    sqrt(mean(pc$cov_tau[rows, rows]))
  }, 0)
  se_bt <- bt$objects$se[match(colnames(X), bt$objects$object)]
  expect_equal(unname(se_bt), unname(se_pc), tolerance = 0.25)
})

test_that("BTL recovers simulated locations with calibrated fit", {
  beta <- c(A = -1.2, B = -0.4, C = 0, D = 0.5, E = 1.1)
  d <- sim_btl(beta, 60, seed = 7)
  ft <- btl(d, "a", "b", "win")
  expect_true(ft$converged)
  loc <- ft$objects$location[match(names(beta), ft$objects$object)]
  expect_gt(cor(loc, beta), 0.98)
  expect_equal(sum(loc), 0, tolerance = 1e-8)          # sum-zero
  expect_lt(abs(max(loc - (beta - mean(beta)))), 0.35) # centred recovery
  # model-true data: pairwise chi-square unremarkable, fit residuals tame
  expect_gt(ft$total_p, 0.01)
  expect_lt(max(abs(ft$objects$fit_resid), na.rm = TRUE), 3)
  expect_true(ft$osi$PSI > 0.7)
})

test_that("an erratic judge is flagged by the judge fit residual", {
  beta <- c(A = -1.2, B = -0.4, C = 0.2, D = 0.6, E = 0.8)
  d <- sim_btl(beta, 80, seed = 11, judges = paste0("J", 1:8))
  # judge J1 answers at random
  sel <- d$judge == "J1"
  d$win[sel] <- ifelse(runif(sum(sel)) < 0.5, d$a[sel], d$b[sel])
  ft <- btl(d, "a", "b", "win", judge = "judge")
  expect_true(ft$clustered)
  jt <- ft$judges
  expect_equal(which.max(jt$fit_resid), match("J1", jt$judge))
  expect_gt(jt$fit_resid[jt$judge == "J1"], 2.5)
  expect_lt(max(jt$fit_resid[jt$judge != "J1"]), 2.5)
})

test_that("ties, extremes, counts, and disconnection are handled", {
  beta <- c(A = -0.8, B = 0, C = 0.8)
  d <- sim_btl(beta, 40, seed = 3)
  d$win[1:5] <- "tie"
  expect_error(btl(d, "a", "b", "win", ties = "error"), "tie")
  ft_d <- btl(d, "a", "b", "win", ties = "drop")
  expect_true(any(grepl("tie", ft_d$notes)))
  ft_h <- btl(d, "a", "b", "win", ties = "half")
  expect_gt(ft_h$n_comparisons, ft_d$n_comparisons)
  # counts replicate rows
  dc <- data.frame(a = c("A", "A", "B"), b = c("B", "C", "C"),
                   win = c("A", "C", "B"), k = c(30, 30, 30))
  dc2 <- rbind(data.frame(a = "A", b = "B", win = "B", k = 10),
               data.frame(a = "B", b = "C", win = "C", k = 10),
               data.frame(a = "A", b = "C", win = "A", k = 10), dc)
  ft_c <- btl(dc2, "a", "b", "win", count = "k")
  expect_equal(ft_c$n_comparisons, 120)
  # an undefeated object is removed with a note
  d2 <- sim_btl(beta, 30, seed = 9)
  d2$win[d2$a == "C" | d2$b == "C"] <- "C"
  expect_error(btl(d2, "a", "b", "win"), "three comparable")
  d3 <- sim_btl(c(beta, D = 0.2, E = -0.2), 30, seed = 9)
  d3$win[d3$a == "E" | d3$b == "E"] <- "E"
  ft3 <- btl(d3, "a", "b", "win")
  expect_true(any(grepl("no wins or no losses", ft3$notes)))
  expect_false("E" %in% ft3$objects$object)
  # disconnected comparison graphs are refused with the components listed
  dd <- data.frame(a = c("A", "A", "C", "C"), b = c("B", "B", "D", "D"),
                   win = c("A", "B", "C", "D"))
  expect_error(btl(dd, "a", "b", "win"), "disconnected")
})

test_that("plot_btl draws and print method runs", {
  beta <- c(A = -1, B = 0, C = 0.4, D = 0.6)
  ft <- btl(sim_btl(beta, 40, seed = 5), "a", "b", "win")
  pdf(NULL); on.exit(dev.off())
  expect_no_error(plot_btl(ft))
  expect_output(print(ft), "Bradley-Terry-Luce")
})
