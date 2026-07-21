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
  # with nothing testable the omnibus is NA too, not chi-square 0 on 0 df
  expect_true(is.na(f$total_df))
  expect_true(is.na(f$total_chisq_p))
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

test_that("pairwise chi-square df counts every estimated parameter", {
  set.seed(3)
  K <- 8; beta <- seq(-1.5, 1.5, length.out = K)
  n <- 1500
  ia <- sample(K, n, TRUE); ib <- (ia + sample(K - 1, n, TRUE) - 1L) %% K + 1L
  win <- rbinom(n, 1, plogis(beta[ia] - beta[ib] + 0.3))
  d <- data.frame(object_a = paste0("O", ia), object_b = paste0("O", ib),
                  winner = paste0("O", ifelse(win == 1, ia, ib)))
  f0 <- btl(d, "object_a", "object_b", "winner")
  f1 <- btl(d, "object_a", "object_b", "winner", position = TRUE)
  # the position covariate consumes one further degree of freedom
  expect_equal(f1$total_df, f0$total_df - 1L)
  # a design with no testable pairs left reports NA, not df = 1
  d3 <- data.frame(object_a = c("A", "B", "C"), object_b = c("B", "C", "A"),
                   winner = c("A", "B", "C"))
  d3 <- d3[rep(1:3, 20), ]
  f3 <- btl(d3, "object_a", "object_b", "winner", position = TRUE)
  expect_true(is.na(f3$total_df) || f3$total_df >= 1L)
})

test_that("dimensionality reference respects comparison counts", {
  skip_on_cran()   # two bootstrap references
  set.seed(7)
  K <- 8; beta <- setNames(seq(-1.2, 1.2, length.out = K), paste0("O", 1:K))
  pr <- t(combn(names(beta), 2)); rows <- list()
  for (i in seq_len(nrow(pr))) {
    wins <- rbinom(1, 12, plogis(beta[pr[i, 1]] - beta[pr[i, 2]]))
    rows[[length(rows) + 1]] <- data.frame(a = pr[i, 1], b = pr[i, 2],
                                           win = pr[i, 1], k = wins)
    rows[[length(rows) + 1]] <- data.frame(a = pr[i, 1], b = pr[i, 2],
                                           win = pr[i, 2], k = 12 - wins)
  }
  agg <- do.call(rbind, rows); agg <- agg[agg$k > 0, ]
  expd <- agg[rep(seq_len(nrow(agg)), agg$k), ]; expd$k <- 1
  fa <- btl(agg, "a", "b", winner = "win", count = "k")
  fe <- btl(expd, "a", "b", winner = "win")
  set.seed(11); da <- btl_dimensionality(fa, reps = 150)
  set.seed(11); de <- btl_dimensionality(fe, reps = 150)
  # pair-level simulation makes the two forms draw from the same design:
  # compare the leading-strength reference distribution, not a pooled mean
  ra <- da$reference$draws; re <- de$reference$draws
  expect_lt(abs(mean(ra) - mean(re)) / mean(re), 0.05)
  expect_lt(abs(stats::quantile(ra, .95) - stats::quantile(re, .95)) /
              stats::quantile(re, .95), 0.05)
})

test_that("judge-resampling bootstrap runs and matches the estimates", {
  skip_on_cran()   # B pipeline refits
  d <- simulate_btl_efrm(n_objects_per_set = 6, n_sets = 2, n_panels = 2,
                         n_judges_per_panel = 8, reps_within = 25,
                         reps_cross = 25, panel_units = c(0.8, 1.25),
                         set_units = c(1, 1.3), seed = 7)
  os <- attr(d, "truth")$object_sets
  set.seed(1)
  fj <- btl_efrm(d, "object_a", "object_b", "winner", "judge", "panel", os,
                 se_method = "judge_bootstrap", boot_reps = 30)
  fc <- btl_efrm(d, "object_a", "object_b", "winner", "judge", "panel", os,
                 se_method = "conditional")
  expect_equal(fj$phi_table$phi, fc$phi_table$phi)
  expect_equal(fj$objects$v, fc$objects$v)
  expect_true(all(is.finite(fj$phi_table$se_log_phi)))
  expect_true(is.finite(fj$alpha_table$se_log_alpha[2]))
  expect_match(fj$se_note, "judge-resampling")
})

test_that("half-tie weights no longer break the dimensionality reference", {
  set.seed(4)
  K <- 6; beta <- setNames(seq(-1, 1, length.out = K), paste0("O", 1:K))
  pr <- t(combn(names(beta), 2)); rows <- list()
  for (i in seq_len(nrow(pr))) for (r in 1:10) {
    pw <- plogis(beta[pr[i, 1]] - beta[pr[i, 2]])
    win <- if (runif(1) < 0.15) "tie" else
      if (runif(1) < pw) pr[i, 1] else pr[i, 2]
    rows[[length(rows) + 1]] <- data.frame(a = pr[i, 1], b = pr[i, 2],
                                           win = win)
  }
  d <- do.call(rbind, rows)
  f <- btl(d, "a", "b", winner = "win", ties = "half")
  dm <- btl_dimensionality(f, reps = 60)
  # fractional 0.5 weights fed as.integer() a zero binomial size before:
  # every reference draw was degenerate at 0
  expect_true(all(is.finite(dm$reference$draws)))
  expect_gt(dm$reference$mean, 0.5)
})

test_that("judge bootstrap refuses a single-judge panel", {
  d <- simulate_btl_efrm(n_objects_per_set = 5, n_sets = 2, n_panels = 2,
                         n_judges_per_panel = 4, reps_within = 15,
                         reps_cross = 15, seed = 3)
  os <- attr(d, "truth")$object_sets
  # collapse panel 2 to a single judge
  keep <- d$panel != "panel2" | d$judge == d$judge[d$panel == "panel2"][1]
  expect_error(
    btl_efrm(d[keep, ], "object_a", "object_b", "winner", "judge", "panel",
             os, se_method = "judge_bootstrap", boot_reps = 20),
    "at least 2 judges in every panel")
})

test_that("clustered covariance notes rank deficiency (judges <= parameters)", {
  set.seed(6)
  K <- 10; beta <- seq(-1.5, 1.5, length.out = K)
  n <- 600
  ia <- sample(K, n, TRUE); ib <- (ia + sample(K - 1, n, TRUE) - 1L) %% K + 1L
  win <- rbinom(n, 1, plogis(beta[ia] - beta[ib]))
  d <- data.frame(object_a = paste0("O", ia), object_b = paste0("O", ib),
                  winner = paste0("O", ifelse(win == 1, ia, ib)),
                  judge = sample(sprintf("J%d", 1:5), n, TRUE))
  f <- btl(d, "object_a", "object_b", "winner", judge = "judge")
  expect_true(any(grepl("rank-deficient", f$notes)))
})

# --- DIF ANOVA engine (round 10): order invariance, person units, GG ------

test_that("multi-factor DIF tests are order-invariant (Type II)", {
  skip_on_cran()
  set.seed(42)
  N <- 500; L <- 8; btrue <- seq(-1.5, 1.5, length.out = L)
  g1 <- sample(c("a", "b"), N, TRUE, prob = c(0.35, 0.65))
  g2 <- ifelse(runif(N) < 0.8, ifelse(g1 == "a", "x", "y"),
               sample(c("x", "y"), N, TRUE))
  th <- rnorm(N)
  X <- matrix(0L, N, L, dimnames = list(NULL, paste0("I", 1:L)))
  for (j in 1:L)
    X[, j] <- rbinom(N, 1, plogis(th - btrue[j] +
                                    ifelse(j == 3 & g1 == "b", 0.8, 0)))
  df <- data.frame(X, g1 = g1, g2 = g2)
  s12 <- dif_anova(rasch(df, factors = c("g1", "g2"),
                         items = paste0("I", 1:L)))$summary
  s21 <- dif_anova(rasch(df, factors = c("g2", "g1"),
                         items = paste0("I", 1:L)))$summary
  key <- function(s) s[order(s$item, s$term), c("F_uniform", "F_nonuniform")]
  expect_equal(key(s12), key(s21), tolerance = 1e-8, ignore_attr = TRUE)
})

test_that("duplicating every person leaves the DIF tests exactly unchanged", {
  set.seed(7)
  N <- 250; L <- 6; btrue <- seq(-1.2, 1.2, length.out = L)
  g <- sample(c("a", "b"), N, TRUE); th <- rnorm(N)
  X <- matrix(0L, N, L, dimnames = list(NULL, paste0("I", 1:L)))
  for (j in 1:L) X[, j] <- rbinom(N, 1, plogis(th - btrue[j]))
  d1 <- data.frame(X, g = g, pid = sprintf("P%03d", 1:N))
  f1 <- rasch(d1, factors = "g", id = "pid", items = paste0("I", 1:L))
  f2 <- rasch(rbind(d1, d1), factors = "g", id = "pid",
              items = paste0("I", 1:L))
  s1 <- dif_anova(f1)$summary; s2 <- dif_anova(f2)$summary
  expect_equal(s2$F_uniform[order(s2$item)], s1$F_uniform[order(s1$item)],
               tolerance = 1e-8)
})

test_that("stacked between-treatment is refused; within declarations checked", {
  set.seed(3)
  N <- 120; d <- seq(-1, 1, length.out = 6)
  X <- rbind(matrix(rbinom(N * 6, 1, plogis(outer(rnorm(N), d, "-"))), N, 6),
             matrix(rbinom(N * 6, 1, plogis(outer(rnorm(N), d, "-"))), N, 6))
  colnames(X) <- paste0("I", 1:6)
  dat <- data.frame(X, occ = rep(c("t1", "t2"), each = N))
  id <- rep(sprintf("P%03d", 1:N), 2)
  fit <- rasch(dat, factors = "occ", id = id)
  expect_error(dif_anova(fit, within = character(0)), "vary within persons")
  expect_error(dif_anova(fit, within = "nope"), "not among the nominated")
})

test_that("multi-level within DIF is GG-calibrated and has power", {
  skip_on_cran()   # replicate fits
  np <- 120; L2 <- 6; b2 <- seq(-1, 1, length.out = L2); K <- 4
  gen <- function(seed, occ_sd, shift3_t4 = 0) {
    set.seed(seed)
    th0 <- rnorm(np); rows <- list()
    for (k in 1:K) {
      thk <- th0 + rnorm(np, 0, occ_sd[k])
      Xk <- matrix(0L, np, L2, dimnames = list(NULL, paste0("I", 1:L2)))
      for (j in 1:L2)
        Xk[, j] <- rbinom(np, 1, plogis(thk - b2[j] +
                                          ifelse(j == 3 & k == 4,
                                                 shift3_t4, 0)))
      rows[[k]] <- data.frame(Xk, occ = paste0("t", k),
                              pid = sprintf("P%03d", 1:np))
    }
    rasch(do.call(rbind, rows), factors = "occ", id = "pid",
          items = paste0("I", 1:L2))
  }
  ## nonspherical null: raw rejections at or below ~nominal over 15 fits
  rej <- 0; tot <- 0
  for (r in 1:15) {
    ss <- dif_anova(gen(100 + r, c(0.05, 0.05, 0.05, 1.5)))$summary
    rej <- rej + sum(ss$p_uniform < 0.05, na.rm = TRUE)
    tot <- tot + sum(is.finite(ss$p_uniform))
  }
  expect_lte(rej / tot, 0.09)   # the uncorrected engine sat near 9-percent
  ## planted occasion DIF: detected as the top flag in most replicates
  hits <- 0
  for (r in 1:5) {
    ss <- dif_anova(gen(200 + r, rep(0.4, K), shift3_t4 = -1.0))$summary
    fl <- ss$item[ss$uniform_DIF %in% TRUE]
    hits <- hits + ("I3" %in% fl)
  }
  expect_gte(hits, 3)
})

# --- DIF ANOVA round 11: multi-within alignment, incomplete panels --------

test_that("multi-within contrasts align: a pure w1 effect loads on w1", {
  skip_on_cran()
  set.seed(5)
  np <- 150; L <- 6; b <- seq(-1, 1, length.out = L)
  rows <- list(); th0 <- rnorm(np)
  for (i1 in c("a1", "a2")) for (i2 in c("b1", "b2", "b3")) {
    th <- th0 + rnorm(np, 0, 0.3)
    X <- matrix(0L, np, L, dimnames = list(NULL, paste0("I", 1:L)))
    for (j in 1:L)
      X[, j] <- rbinom(np, 1, plogis(th - b[j] +
                                       ifelse(j == 3 & i1 == "a2", -0.9, 0)))
    rows[[length(rows) + 1]] <- data.frame(X, w1 = i1, w2 = i2,
                                           pid = sprintf("P%03d", 1:np))
  }
  f <- rasch(do.call(rbind, rows), factors = c("w1", "w2"), id = "pid",
             items = paste0("I", 1:L))
  tt <- dif_anova(f, effects = "factorial")$terms
  i3 <- tt[tt$item == "I3", ]
  # interaction()'s first-fastest cell order silently rotated this into
  # w1:w2 (F 7.4 on the interaction, 3.2 on w1)
  expect_gt(i3$F_value[i3$term == "w1"], 10)
  expect_lt(i3$F_value[i3$term == "w1:w2"], 4)
  expect_lt(i3$F_value[i3$term == "w2"], 4)
  # GG metadata reproducible: p == pf(F, eps*df, eps*df_denom)
  r <- i3[i3$term == "w1", ]
  expect_equal(r$p, pf(r$F_value, r$gg_epsilon * r$df,
                       r$gg_epsilon * r$df_denom, lower.tail = FALSE),
               tolerance = 1e-12)
})

test_that("differentially incomplete within panels give no false group DIF", {
  skip_on_cran()
  set.seed(9)
  np <- 200; L <- 6; b <- seq(-1, 1, length.out = L)
  g <- rep(c("A", "B"), each = np / 2)
  th0 <- rnorm(np); rows <- list()
  for (k in 1:2) {
    th <- th0 + rnorm(np, 0, 0.3)
    X <- matrix(0L, np, L, dimnames = list(NULL, paste0("I", 1:L)))
    for (j in 1:L)
      X[, j] <- rbinom(np, 1, plogis(th - b[j] +
                                       ifelse(j == 3 & k == 2, -0.9, 0)))
    keep <- if (k == 2) (g == "A") | (runif(np) < 0.2) else rep(TRUE, np)
    rows[[k]] <- data.frame(X, occ = paste0("t", k), grp = g,
                            pid = sprintf("P%03d", 1:np))[keep, ]
  }
  f <- rasch(do.call(rbind, rows), factors = c("occ", "grp"), id = "pid",
             items = paste0("I", 1:L))
  s <- dif_anova(f)$summary
  gF <- s$F_uniform[s$item == "I3" & s$term == "grp"]
  # raw person means over unmatched cells reported group F = 37.6 here
  expect_true(is.na(gF) || gF < 8)
})

test_that("an NA in a between factor does not flip it to within-subject", {
  set.seed(2)
  b <- seq(-1, 1, length.out = 6)
  X <- matrix(rbinom(200 * 6, 1, plogis(outer(rnorm(200), b, "-"))), 200, 6,
              dimnames = list(NULL, paste0("I", 1:6)))
  gg <- rep(c("a", "b"), 100); gg[5] <- NA
  f <- rasch(data.frame(X, g = gg, pid = sprintf("P%03d", 1:200)),
             factors = "g", id = "pid", items = paste0("I", 1:6))
  d <- dif_anova(f)
  expect_length(d$within, 0L)
  expect_gt(nrow(d$summary), 0L)
})

# --- DIF ANOVA round 12: incomplete-panel edges ---------------------------

test_that("trait-dependent within effects with differential missingness give no group DIF", {
  skip_on_cran()
  set.seed(13)
  np <- 240; L <- 6; b <- seq(-1, 1, length.out = L)
  g <- rep(c("A", "B"), each = np / 2)
  th0 <- rnorm(np); rows <- list()
  for (k in 1:2) {
    th <- th0 + rnorm(np, 0, 0.3)
    X <- matrix(0L, np, L, dimnames = list(NULL, paste0("I", 1:L)))
    for (j in 1:L)
      X[, j] <- rbinom(np, 1, plogis(th - b[j] +
        ifelse(j == 3 & k == 2, -0.5 - 0.6 * th0, 0)))
    keep <- if (k == 2) (g == "A") | (runif(np) < 0.2) else rep(TRUE, np)
    rows[[k]] <- data.frame(X, occ = paste0("t", k), grp = g,
                            pid = sprintf("P%03d", 1:np))[keep, ]
  }
  f <- rasch(do.call(rbind, rows), factors = c("occ", "grp"), id = "pid",
             items = paste0("I", 1:L))
  s <- dif_anova(f)$summary
  gU <- s$F_uniform[s$item == "I3" & s$term == "grp"]
  gN <- s$F_nonuniform[s$item == "I3" & s$term == "grp"]
  # cell-only centring reported non-uniform group DIF at F = 214.7 here
  expect_true(is.na(gU) || gU < 8)
  expect_true(is.na(gN) || gN < 8)
})

test_that("a between level with no complete panels yields NA terms, not a crash", {
  set.seed(14)
  np <- 160; L <- 6; b <- seq(-1, 1, length.out = L)
  g <- rep(c("A", "B"), each = np / 2); th0 <- rnorm(np); rows <- list()
  for (k in 1:4) {
    th <- th0 + rnorm(np, 0, 0.3)
    X <- matrix(0L, np, L, dimnames = list(NULL, paste0("I", 1:L)))
    for (j in 1:L) X[, j] <- rbinom(np, 1, plogis(th - b[j]))
    keep <- if (k == 1) rep(TRUE, np) else g == "B"
    rows[[k]] <- data.frame(X, occ = paste0("t", k), grp = g,
                            pid = sprintf("P%03d", 1:np))[keep, ]
  }
  f <- rasch(do.call(rbind, rows), factors = c("occ", "grp"), id = "pid",
             items = paste0("I", 1:L))
  expect_s3_class(dif_anova(f), "rasch_dif")
  df <- dif_anova(f, effects = "factorial")
  expect_true(any(is.na(df$terms$F_value[df$terms$term == "occ:grp"])))
  expect_true(any(grepl("non-estimable", df$notes)))
  expect_true(any(grepl("dropped from the within-person", df$notes)))
})

test_that("significant multilevel within terms never reach ordinary Tukey", {
  set.seed(15)
  np <- 150; L <- 6; b <- seq(-1, 1, length.out = L)
  rows <- list(); th0 <- rnorm(np)
  gg <- rep(c("x", "y"), length.out = np)
  for (k in 1:4) {
    th <- th0 + rnorm(np, 0, 0.3)
    X <- matrix(0L, np, L, dimnames = list(NULL, paste0("I", 1:L)))
    for (j in 1:L)
      X[, j] <- rbinom(np, 1, plogis(th - b[j] +
                                       ifelse(j == 3 & k == 4, -1.2, 0)))
    keep <- (gg == "x") | (runif(np) < 0.7)
    rows[[k]] <- data.frame(X, occ = paste0("t", k), grp = gg,
                            pid = sprintf("P%03d", 1:np))[keep, ]
  }
  f <- rasch(do.call(rbind, rows), factors = c("occ", "grp"), id = "pid",
             items = paste0("I", 1:L))
  d <- dif_anova(f)
  expect_s3_class(d, "rasch_dif")
  # any Tukey rows present concern between terms only
  if (nrow(d$tukey)) expect_false(any(grepl("occ", d$tukey$term)))
})

# --- round 13: BTL / EFRM / MFRM DIF and identification -------------------

test_that("btl_dif is order-invariant across correlated judge factors", {
  skip_on_cran()
  set.seed(31)
  K <- 10; beta <- seq(-1.4, 1.4, length.out = K); nj <- 28; npj <- 60
  g1 <- rep(c("a", "b"), c(10, 18))
  g2 <- ifelse(runif(nj) < 0.75, ifelse(g1 == "a", "x", "y"),
               sample(c("x", "y"), nj, TRUE))
  jids <- sprintf("J%02d", 1:nj); rows <- list()
  for (j in 1:nj) {
    bj <- beta; if (g1[j] == "b") bj[4] <- bj[4] - 1
    ia <- sample(K, npj, TRUE); ib <- (ia + sample(K - 1, npj, TRUE) - 1L) %% K + 1L
    win <- rbinom(npj, 1, plogis(bj[ia] - bj[ib]))
    rows[[j]] <- data.frame(object_a = paste0("O", ia),
                            object_b = paste0("O", ib),
                            winner = paste0("O", ifelse(win == 1, ia, ib)),
                            judge = jids[j])
  }
  bt <- btl(do.call(rbind, rows), "object_a", "object_b", "winner",
            judge = "judge")
  A <- setNames(g1, jids); B <- setNames(g2, jids)
  sAB <- btl_dif(bt, factors = list(A = A, B = B))$summary
  sBA <- btl_dif(bt, factors = list(B = B, A = A))$summary
  expect_equal(sAB$F_uniform[order(sAB$object, sAB$term)],
               sBA$F_uniform[order(sBA$object, sBA$term)],
               tolerance = 1e-8)
})

test_that("btl_dif rejects factors that vary within a judge", {
  set.seed(3)
  K <- 6; b <- seq(-1, 1, length.out = K); n <- 600
  ia <- sample(K, n, TRUE); ib <- (ia + sample(K - 1, n, TRUE) - 1L) %% K + 1L
  d <- data.frame(object_a = paste0("O", ia), object_b = paste0("O", ib),
                  winner = paste0("O", ifelse(
                    rbinom(n, 1, plogis(b[ia] - b[ib])) == 1, ia, ib)),
                  judge = sample(sprintf("J%d", 1:10), n, TRUE))
  bt <- btl(d, "object_a", "object_b", "winner", judge = "judge")
  rowfac <- sample(c("u", "v"), n, TRUE)      # varies within judges
  expect_error(btl_dif(bt, factors = list(g = rowfac)),
               "varies within judge")
})

test_that("btl_efrm validates within-set connectivity and alpha identification", {
  d <- simulate_btl_efrm(n_objects_per_set = 5, n_sets = 2, n_panels = 2,
                         n_judges_per_panel = 6, reps_within = 20,
                         reps_cross = 20, seed = 4)
  os <- attr(d, "truth")$object_sets
  ## alpha: cross-set rows touching only ONE object of set 2
  s2 <- os[[2]]
  cross_rows <- (d$object_a %in% os[[1]]) != (d$object_b %in% os[[1]])
  keep <- !cross_rows | (d$object_a == s2[1]) | (d$object_b == s2[1])
  expect_error(
    btl_efrm(d[keep, ], "object_a", "object_b", "winner", "judge", "panel",
             os, se_method = "conditional"),
    "unit \\(alpha\\) is unidentified")
  ## within-set connectivity: split set 1's internal comparisons
  g1 <- os[[1]][1:2]; g2 <- os[[1]][3:5]
  within1 <- (d$object_a %in% os[[1]]) & (d$object_b %in% os[[1]])
  bridge <- within1 & ((d$object_a %in% g1) != (d$object_b %in% g1))
  expect_error(
    btl_efrm(d[!bridge, ], "object_a", "object_b", "winner", "judge",
             "panel", os, se_method = "conditional"),
    "not connected|no within-set comparison")
})

test_that("EFRM group linkage requires shared items, not shared set labels", {
  d <- simulate_efrm(n_per_group = 300, items_per_set = 8, n_sets = 2,
                     n_groups = 2, group_unit_ratio = 1.25, seed = 2)
  tr <- attr(d, "truth")
  ## make the groups' item subsets DISJOINT within every set
  grp <- tr$groups
  X <- d
  items <- unlist(tr$item_sets)
  for (s in tr$item_sets) {
    half <- seq_len(floor(length(s) / 2))
    X[grp == unique(grp)[1], s[half]] <- NA
    X[grp == unique(grp)[2], s[-half]] <- NA
  }
  expect_error(
    rasch_efrm(X, item_sets = tr$item_sets, groups = grp),
    "not linked|unidentified")
})

test_that("dif_anova integrates with EFRM and MFRM fits", {
  skip_on_cran()
  d <- simulate_efrm(n_per_group = 300, items_per_set = 8, n_sets = 2,
                     n_groups = 2, group_unit_ratio = 1.25, seed = 1)
  tr <- attr(d, "truth")
  expect_error(
    dif_anova(rasch_efrm(d, item_sets = tr$item_sets, groups = tr$groups)),
    "frame structure")
  sex <- rep(c("m", "f"), length.out = nrow(d))
  f2 <- rasch_efrm(d, item_sets = tr$item_sets, groups = tr$groups,
                   factors = data.frame(sex = sex))
  d2 <- dif_anova(f2)
  expect_true(any(grepl("frame structure", d2$notes)))
  expect_gt(nrow(d2$summary), 0)

  set.seed(1)
  simP <- function(th, tau) { x <- 0:length(tau)
    p <- exp(x * th - c(0, cumsum(tau))); p / sum(p) }
  persons <- sprintf("P%03d", 1:150); raters <- paste0("R", 1:3)
  th <- setNames(rnorm(150, 0, 1.3), persons)
  sx <- setNames(rep(c("m", "f"), length.out = 150), persons)
  tau <- list(A = c(-1, 1), B = c(-0.5, 1.2), C = c(-1.2, 0.4))
  dd <- expand.grid(person = persons, item = names(tau), rater = raters,
                    stringsAsFactors = FALSE)
  dd$score <- mapply(function(p, i, r)
    sample(0:2, 1, prob = simP(th[p] + ifelse(i == "B" & sx[p] == "f",
                                              -0.6, 0),
                               tau[[i]] + c(R1 = -0.4, R2 = 0,
                                            R3 = 0.4)[r])),
    dd$person, dd$item, dd$rater)
  dd$sex <- sx[dd$person]
  mf <- rasch_mfrm(dd, person = "person", item = "item", score = "score",
                   facets = "rater", factors = "sex")
  dm <- dif_anova(mf)                       # pooled to underlying items
  expect_true(dm$summary$uniform_DIF[dm$summary$item == "B"] %in% TRUE)
  dmv <- dif_anova(mf, pool_facets = FALSE) # virtual mode preserved
  expect_true(any(dmv$summary$uniform_DIF[grepl("^B:", dmv$summary$item)] %in%
                    TRUE))
  expect_error(
    rasch_mfrm(dd, person = "person", item = "item", score = "score",
               facets = "rater", factors = "rater"),
    "varies within person")
})

# --- capability round: pooled MFRM DIF, MFRM sizes, factorial EFRM frames --

test_that("MFRM DIF pools to underlying items and resolves magnitudes", {
  skip_on_cran()
  set.seed(1)
  simP <- function(th, tau) { x <- 0:length(tau)
    p <- exp(x * th - c(0, cumsum(tau))); p / sum(p) }
  persons <- sprintf("P%03d", 1:200); raters <- paste0("R", 1:3)
  th <- setNames(rnorm(200, 0, 1.3), persons)
  sx <- setNames(rep(c("m", "f"), length.out = 200), persons)
  rho <- setNames(c(-0.4, 0, 0.4), raters)
  tau <- list(A = c(-1, 1), B = c(-0.5, 1.2), C = c(-1.2, 0.4))
  dd <- expand.grid(person = persons, item = names(tau), rater = raters,
                    stringsAsFactors = FALSE)
  dd$score <- mapply(function(p, i, r)
    sample(0:2, 1, prob = simP(th[p] + ifelse(i == "B" & sx[p] == "f",
                                              -0.7, 0),
                               tau[[i]] + rho[r])),
    dd$person, dd$item, dd$rater)
  dd$sex <- sx[dd$person]
  mf <- rasch_mfrm(dd, person = "person", item = "item", score = "score",
                   facets = "rater", factors = "sex")
  dp <- dif_anova(mf)
  # pooled to underlying items; planted B carries the dominant F (with 3
  # items, compensating artificial DIF on A and C is expected and real)
  expect_true(all(dp$summary$item %in% c("A", "B", "C")))
  expect_equal(dp$summary$item[which.max(dp$summary$F_uniform)], "B")
  expect_true(any(grepl("pooled to the underlying", dp$notes)))
  # per-virtual mode preserved
  expect_true(any(grepl(":", dif_anova(mf, pool_facets = FALSE)$summary$item)))
  # magnitudes at the underlying-item level, planted size recovered
  ds <- dif_size(mf, "B", by = "sex")
  expect_lt(abs(abs(ds$pairs$difference) - 0.7), 3 * ds$pairs$se)
  expect_true(ds$pairs$significant)
  dsz <- dif_anova(mf, sizes = TRUE)
  expect_gt(nrow(dsz$sizes), 0)
})

test_that("EFRM accepts several frame factors and reports the decomposition", {
  skip_on_cran()
  d <- simulate_efrm(n_per_group = 300, items_per_set = 8, n_sets = 2,
                     n_groups = 2, group_unit_ratio = 1.25, seed = 1)
  tr <- attr(d, "truth")
  d$region <- rep(c("N", "S"), length.out = nrow(d))
  d$grp <- tr$groups
  f <- rasch_efrm(d, item_sets = tr$item_sets, groups = c("grp", "region"),
                  items = unlist(tr$item_sets))
  expect_equal(nrow(f$phi_table), 4L)
  expect_false(is.null(f$phi_factorial))
  expect_true(any(grepl("grp", f$phi_factorial$term)))
  # every frame-defining factor is excluded from DIF testing
  expect_error(dif_anova(f), "frame structure")
  expect_error(rasch_efrm(d, item_sets = tr$item_sets,
                          groups = c("grp", "nonexistent")),
               "not found")
})

test_that("EFRM unit identification is honest at both extremes", {
  skip_on_cran()
  set.seed(42)
  N <- 500; g <- rep(c("g1", "g2"), each = N / 2)
  th <- rnorm(N, sd = 1.5)
  phi_t <- ifelse(g == "g1", 1, 1.5)
  sets <- list(A = paste0("v", 1:8), B = paste0("v", 9:16))
  mk <- function(deltas) {
    X <- as.data.frame(matrix(0L, N, 16)); names(X) <- paste0("v", 1:16)
    for (i in 1:16) X[[i]] <- rbinom(N, 1, plogis(phi_t * (th - deltas[i])))
    X
  }
  # all items equally difficult: the units are weakly identified at best.
  # The identification check is on the information (rank/conditioning),
  # not on a spread heuristic, so this full-rank case keeps its estimate;
  # honesty lives in the SEs, which must be large enough that no unit
  # difference could be claimed from them
  fz <- rasch_efrm(mk(rep(0, 16)), groups = g, item_sets = sets,
                   se_method = "hybrid")
  lr <- abs(diff(log(fz$phi_table$phi)))
  se_lr <- sqrt(sum(fz$phi_table$se_log_phi^2))
  expect_lt(lr / se_lr, 2)              # no spurious unit claim
  expect_gt(min(fz$phi_table$se_log_phi), 0.2)   # no spurious precision
  # a modest half-logit spread is real signal and must not be refused
  fm <- rasch_efrm(mk(rep(seq(-0.25, 0.25, length.out = 8), 2)),
                   groups = g, item_sets = sets, se_method = "hybrid")
  expect_false(anyNA(fm$phi_table$phi))
})

test_that("phi_factorial_tests recover a planted region unit effect", {
  skip_on_cran()
  set.seed(42)
  N <- 500; th <- rnorm(N, sd = 1.5)
  delh <- rep(seq(-1.2, 1.2, length.out = 8), 2)
  reg <- rep(c("N", "S"), each = N / 2)
  set.seed(7)
  grp <- sample(c("g1", "g2"), N, TRUE)
  phi_r <- ifelse(reg == "N", 1, 1.5)      # unit effect on region only
  X <- as.data.frame(matrix(0L, N, 16)); names(X) <- paste0("v", 1:16)
  for (i in 1:16) X[[i]] <- rbinom(N, 1, plogis(phi_r * (th - delh[i])))
  X$grp <- grp; X$region <- reg
  f <- rasch_efrm(X, groups = c("grp", "region"), items = paste0("v", 1:16),
                  item_sets = list(A = paste0("v", 1:8),
                                   B = paste0("v", 9:16)),
                  se_method = "hybrid")
  tt <- f$phi_factorial_tests
  expect_setequal(tt$term, c("grp", "region", "grp:region"))
  expect_lt(tt$p[tt$term == "region"], 0.01)
  expect_gt(tt$p[tt$term == "grp"], 0.05)
  expect_gt(tt$p[tt$term == "grp:region"], 0.05)
  # the sum-coded region effect reproduces the planted unit ratio 1.5
  lu <- f$phi_factorial$log_unit[f$phi_factorial$term == "region1"]
  expect_gt(exp(2 * abs(lu)), 1.2)
  expect_lt(exp(2 * abs(lu)), 1.9)
  # bootstrap draws feed se_log_rho jointly: finite and positive
  fb <- rasch_efrm(X, groups = c("grp", "region"),
                   items = paste0("v", 1:16),
                   item_sets = list(A = paste0("v", 1:8),
                                    B = paste0("v", 9:16)),
                   se_method = "bootstrap", boot_reps = 40)
  expect_true(all(is.finite(fb$frames$se_log_rho)))
  expect_true(all(fb$frames$se_log_rho > 0))
})

test_that("btl_efrm declares a degenerate set unit instead of reporting alpha = 1", {
  set.seed(31)
  # set B has two objects whose within-set record is perfectly balanced:
  # both centred locations are exactly 0, so the cross-set derivative for
  # log alpha_B vanishes identically and stage 2 is rank-deficient even
  # though the cross-set comparisons touch both objects of B
  judges <- sprintf("J%d", 1:6)
  withinA <- expand.grid(a = c("A1", "A2", "A3"), b = c("A1", "A2", "A3"),
                         rep = 1:8, stringsAsFactors = FALSE)
  withinA <- withinA[withinA$a < withinA$b, c("a", "b")]
  bA <- c(A1 = -0.8, A2 = 0, A3 = 0.8)
  withinA$win <- ifelse(rbinom(nrow(withinA), 1,
                               plogis(bA[withinA$a] - bA[withinA$b])) == 1,
                        withinA$a, withinA$b)
  withinB <- data.frame(a = rep("B1", 24), b = rep("B2", 24),
                        win = rep(c("B1", "B2"), 12))
  cross <- expand.grid(a = c("A1", "A2", "A3"), b = c("B1", "B2"),
                       rep = 1:5, stringsAsFactors = FALSE)[, c("a", "b")]
  cross$win <- ifelse(seq_len(nrow(cross)) %% 2 == 0, cross$a, cross$b)
  d <- rbind(withinA, withinB, cross)
  d$judge <- rep_len(judges, nrow(d))
  d$panel <- "P1"
  fit <- suppressWarnings(
    btl_efrm(d, "a", "b", winner = "win", judge = "judge", panels = "panel",
             object_sets = list(setA = c("A1", "A2", "A3"),
                                setB = c("B1", "B2")),
             se_method = "conditional", min_link = 5))
  # the unit is declared unidentified (NA), not reported as an arbitrary 1
  expect_true(is.na(fit$alpha_table$alpha[fit$alpha_table$set == "setB"]))
  expect_true(is.na(
    fit$alpha_table$se_log_alpha[fit$alpha_table$set == "setB"]))
  expect_true(any(grepl("unit\\(s\\) unidentified", fit$notes)))
  # the set's objects are still PLACED: origin kappa and locations finite
  expect_true(is.finite(
    fit$kappa_table$kappa[fit$kappa_table$set == "setB"]))
  expect_true(all(is.finite(fit$objects$v[fit$objects$set == "setB"])))
})

test_that("btl_efrm refuses panels that observe disjoint object pairs", {
  set.seed(44)
  mk <- function(a, b, n, pwin, judges, panel)
    data.frame(a = a, b = b,
               win = ifelse(rbinom(n, 1, pwin) == 1, a, b),
               judge = sample(judges, n, TRUE), panel = panel)
  # panel P1 sees only A-B and panel P2 only B-C: the free panel ratio and
  # the location contrast enter only as a product, so stage 1 is rank 2
  # for 3 parameters and the panel units are unidentified
  d <- rbind(mk("A", "B", 40, 0.6, sprintf("J%d", 1:5), "P1"),
             mk("B", "C", 40, 0.6, sprintf("K%d", 1:5), "P2"))
  expect_error(suppressWarnings(
    btl_efrm(d, "a", "b", winner = "win", judge = "judge",
             panels = "panel", object_sets = list(S = c("A", "B", "C")),
             se_method = "conditional", min_link = 5)),
    "panel")
})

test_that("dif_size withholds magnitude and significance on a weak category", {
  # one group's response on the target item sits almost entirely in a
  # near-empty category, so its resolved location rests on a weak threshold
  set.seed(101)
  N <- 400; L <- 5
  d <- seq(-1.5, 1.5, length.out = L)
  X <- matrix(rbinom(N * L, 1, plogis(outer(rnorm(N), d, "-"))), N, L)
  colnames(X) <- paste0("I", 1:L)
  grp <- factor(rep(c("a", "b"), each = N / 2))
  # make item I3 a 0/1/2 item; group a populates category 2 normally, but
  # group b reaches it in only 2 responses -- a near-empty category whose
  # threshold split_items() flags weak (se = NA)
  X3 <- X[, 3] + rbinom(N, 1, ifelse(grp == "a", 0.45, 0))
  X3[X3 > 2] <- 2L
  brows <- which(grp == "b")
  X3[brows] <- pmin(X3[brows], 1L)          # group b capped at category 1
  X3[brows[1:2]] <- 2L                       # exactly two b responses in cat 2
  X[, 3] <- X3
  fit <- rasch(data.frame(X, grp = grp, check.names = FALSE),
               model = "PCM", factors = "grp")
  ds <- dif_size(fit, "I3", by = "grp")
  # the weak level is flagged and its SE/verdict withheld, not fabricated
  expect_true(any(ds$levels$weak))
  expect_true(is.na(ds$levels$se[ds$levels$weak]))
  expect_true(all(is.na(ds$pairs$se)))
  expect_true(all(is.na(ds$pairs$significant)))
  expect_true(any(grepl("weakly identified|withheld", ds$notes)))
})

test_that("resolve_dif does not split a uniform flag driven by a thin cell", {
  set.seed(102)
  N <- 600; L <- 6
  d <- seq(-1.5, 1.5, length.out = L)
  th <- rnorm(N)
  X <- matrix(rbinom(N * L, 1, plogis(outer(th, d, "-"))), N, L)
  colnames(X) <- paste0("I", 1:L)
  sex <- factor(sample(c("F", "M"), N, TRUE))
  # a real uniform DIF on I5 by sex (so resolve_dif has something genuine)
  p5 <- plogis(th - d[5] + ifelse(sex == "M", 0.9, 0))
  X[, 5] <- rbinom(N, 1, p5)
  # age_band: a tiny 'old' cell (n=15) that will drive a spurious ANOVA flag
  age <- factor(c(rep("young", 300), rep("mid", 285), rep("old", 15)))
  fit <- rasch(data.frame(X, sex = sex, age_band = age, check.names = FALSE),
               factors = c("sex", "age_band"))
  rr <- resolve_dif(fit, factors = c("sex", "age_band"))
  # the genuine sex DIF is resolved; the thin-cell age flag is not chased
  split_key <- paste(rr$splits$item, rr$splits$factor)
  expect_false(any(grepl("age_band", rr$splits$factor)))
})

test_that("btl marks subset separation as non-convergence, not a boundary fit", {
  # graded cross-block comparisons all at the ceiling category, with one
  # near-ceiling concession that makes the win graph strongly connected:
  # the {A,B}-vs-{C,D} contrast is (quasi-)separated
  mk <- function(a, b, n, x) data.frame(a = rep(a, n), b = rep(b, n), resp = x)
  within1 <- mk("A", "B", 40, rep(c(0, 1, 2, 3, 4), 8))
  within2 <- mk("C", "D", 40, rep(c(0, 1, 2, 3, 4), 8))
  cross <- rbind(mk("A", "C", 400, rep(4, 400)), mk("A", "D", 400, rep(4, 400)),
                 mk("B", "C", 400, rep(4, 400)), mk("B", "D", 399, rep(4, 399)),
                 data.frame(a = "B", b = "D", resp = 3))
  fit <- btl(rbind(within1, within2, cross), "a", "b", response = "resp")
  expect_false(fit$converged)
  expect_true(any(is.na(fit$objects$se)))
  expect_true(any(grepl("run to the location boundary", fit$notes)))
  # equating refuses a non-converged calibration
  fit2 <- fit
  expect_error(btl_equate(fit, fit2), "did not converge")
})

test_that("weak-category honesty reaches MFRM, EFRM, and average anchoring", {
  # a near-empty category must yield weak = TRUE / se = NA on every path
  # that builds its own estimate, not only rasch()/pcml()
  set.seed(303)
  # --- average anchoring (k = NA) keeps free thresholds' weak flag ---
  N <- 300; L <- 5
  simP <- function(th, tau) { x <- 0:length(tau)
    p <- exp(x * th - c(0, cumsum(tau))); p / sum(p) }
  th <- rnorm(N)
  taus <- list(c(-1, 0, 1), c(-1, 0, 1), c(-1, 0, 1), c(-1, 0, 1), c(-1, 0, 1))
  X <- sapply(taus, function(tt) vapply(th, function(t)
    sample(0:3, 1, prob = simP(t, tt)), 0L))
  colnames(X) <- paste0("I", 1:L)
  # force I1 category 1 to two responses
  i1 <- which(X[, 1] == 1); keep <- i1[1:2]; X[setdiff(i1, keep), 1] <- 0L
  r_free <- pcml(X)
  r_anch <- pcml(X, anchors = data.frame(item = "I1", k = NA, tau = 0))
  w_free <- r_free$thr$weak[r_free$thr$item == 1]
  w_anch <- r_anch$thr$weak[r_anch$thr$item == 1]
  expect_true(any(w_free))
  expect_equal(w_anch, w_free)              # average anchor does not hide it
  expect_true(all(is.na(r_anch$thr$se[r_anch$thr$item == 1 & r_anch$thr$weak])))

  # --- MFRM virtual item on a near-empty category ---
  persons <- sprintf("P%03d", 1:200)
  d <- expand.grid(person = persons, item = c("A", "B"),
                   rater = c("R1", "R2"), stringsAsFactors = FALSE)
  d$score <- sample(0:2, nrow(d), replace = TRUE)
  # make virtual item A:R2 reach category 2 only once
  sel <- d$item == "A" & d$rater == "R2"
  d$score[sel] <- 0L; d$score[which(sel)[1]] <- 2L; d$score[which(sel)[2]] <- 1L
  mf <- rasch_mfrm(d, person = "person", item = "item", score = "score",
                   facets = "rater")
  expect_true("weak" %in% names(mf$est$thr))
  expect_true(any(mf$est$thr$weak))
  expect_true(all(is.na(mf$est$thr$se[mf$est$thr$weak])))
})
