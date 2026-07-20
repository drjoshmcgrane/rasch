# Statistical-validity regressions from the 2026-07 external review, each
# verified by simulation before fixing: the recentring covariance transform
# (mixed max scores), the Warm WLE discrimination cancellation, honest
# chi-square degrees of freedom, covariance-correct equating drift tests,
# the few-judges clustering guard, and the judge-level DIF ANOVA.

test_that("mixed-max-score item location SEs are calibrated (cov transform)", {
  skip_on_cran()   # 60 replicate fits
  set.seed(9)
  L <- 6; m <- rep(c(1L, 3L), each = 3)
  btrue <- seq(-1, 1, length.out = L)
  tau_l <- lapply(1:L, function(j) btrue[j] + seq(-0.6, 0.6, length.out = m[j]))
  gen <- function(N) {
    th <- rnorm(N)
    X <- matrix(0L, N, L, dimnames = list(NULL, sprintf("I%d", 1:L)))
    for (j in 1:L) for (i in 1:N) {
      d <- th[i] - tau_l[[j]]
      p <- c(1, exp(cumsum(d))); X[i, j] <- sample(0:m[j], 1, prob = p / sum(p))
    }
    X
  }
  R <- 60; locs <- matrix(NA, R, L); ses <- matrix(NA, R, L)
  for (r in 1:R) {
    f <- pcml(gen(500))
    locs[r, ] <- vapply(1:L, function(j) mean(f$thr$tau[f$thr$item == j]), 0)
    ses[r, ] <- vapply(1:L, function(j) {
      rows <- f$thr$id[f$thr$item == j]
      sqrt(mean(f$cov_tau[rows, rows]))
    }, 0)
  }
  ratio <- apply(locs, 2, sd) / colMeans(ses)
  # before the transform, dichotomous items sat ~0.90 and 3-threshold items
  # ~1.20 systematically; after it every item is within noise of 1
  expect_true(all(ratio > 0.8 & ratio < 1.25))
  expect_lt(abs(mean(ratio) - 1), 0.12)
})

test_that("Warm WLE is invariant to a common discrimination", {
  tau_list <- list(c(-1), c(0), c(1), c(-0.5, 0.5))
  for (a in c(0.5, 2)) {
    w <- person_wle(tau_list, disc = a)
    for (R in 1:4) {
      # exact WLE: root of the weighted score a(R-E) + a^3 mu3 / (2 a^2 V)
      obj <- function(th) {
        mo <- lapply(tau_list, item_moments, theta = th, disc = a)
        E <- sum(sapply(mo, `[[`, "E")); V <- sum(sapply(mo, `[[`, "V"))
        m3 <- sum(sapply(mo, `[[`, "mu3"))
        a * (R - E) + a^3 * m3 / (2 * a^2 * V)
      }
      exact <- uniroot(obj, c(-30, 30), tol = 1e-12)$root
      expect_equal(unname(w$theta[as.character(R)]), exact, tolerance = 1e-6)
    }
  }
})

test_that("an item with fewer than two class intervals gets NA, not df = 1", {
  # two items, all non-extreme persons share one raw score -> one interval
  X <- cbind(I1 = rep(c(0L, 1L), 60), I2 = rep(c(1L, 0L), 60))
  f <- rasch(X, n_groups = 2)
  expect_true(all(is.na(f$items$df)))
  expect_true(all(is.na(f$items$chisq)))
  expect_true(all(is.na(f$items$p)))
  expect_equal(f$total_df, 0L)
})

test_that("equating drift tests are calibrated under the null", {
  skip_on_cran()   # 80 replicate pairs of fits
  set.seed(42)
  L <- 8; btrue <- seq(-1.5, 1.5, length.out = L)
  mk <- function() {
    X <- matrix(rbinom(400 * L, 1, plogis(outer(rnorm(400), btrue, "-"))),
                400, L, dimnames = list(NULL, paste0("I", 1:L)))
    rasch(as.data.frame(X))
  }
  rej <- 0; tot <- 0
  for (r in 1:80) {
    eq <- equate_tests(mk(), mk())
    rej <- rej + sum(eq$table$p < 0.05, na.rm = TRUE)
    tot <- tot + sum(is.finite(eq$table$p))
  }
  # naive sqrt(v) denominators (no shift covariance) were mis-calibrated;
  # the projected covariance restores ~nominal rejection
  expect_gt(rej / tot, 0.02)
  expect_lt(rej / tot, 0.09)
})

test_that("clustered SEs refuse a single judge and note few judges", {
  d <- data.frame(object_a = rep(paste0("O", 1:4), 30),
                  object_b = rep(paste0("O", c(2:4, 1)), 30))
  set.seed(2); d$winner <- ifelse(runif(120) < .5, d$object_a, d$object_b)
  d$judge <- "J1"
  expect_error(btl(d, "object_a", "object_b", "winner", judge = "judge"),
               "at least 2 judges")
  d$judge <- rep(sprintf("J%d", 1:6), 20)
  f <- btl(d, "object_a", "object_b", "winner", judge = "judge")
  expect_true(any(grepl("judge clusters", f$notes)))
})

test_that("btl_dif does not flag under judge heterogeneity with null groups", {
  skip_on_cran()   # several fits
  simnull <- function(seed) {
    set.seed(seed)
    K <- 8; beta <- seq(-1.2, 1.2, length.out = K); nj <- 12; npj <- 60
    rows <- list()
    for (j in 1:nj) {
      bj <- beta + rnorm(K, 0, 0.6)
      ia <- sample(K, npj, TRUE); ib <- (ia + sample(K - 1, npj, TRUE) - 1L) %% K + 1L
      win <- rbinom(npj, 1, plogis(bj[ia] - bj[ib]))
      rows[[j]] <- data.frame(object_a = paste0("O", ia),
                              object_b = paste0("O", ib),
                              winner = paste0("O", ifelse(win == 1, ia, ib)),
                              judge = sprintf("J%02d", j))
    }
    d <- do.call(rbind, rows)
    bt <- btl(d, "object_a", "object_b", "winner", judge = "judge")
    grp <- setNames(rep(c("A", "B"), each = nj / 2), sprintf("J%02d", 1:nj))
    df <- btl_dif(bt, factors = list(g = grp))
    sum(df$summary$uniform_DIF %in% TRUE)
  }
  # comparison-level pseudoreplication flagged 6 of 10 such nulls
  flags <- vapply(1:5, simnull, 0)
  expect_lte(sum(flags > 0), 1)
})
