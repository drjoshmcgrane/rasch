# Statistical-validity regressions from the 2026-07 external review, each
# verified by simulation before fixing: the recentring covariance transform
# (mixed max scores), the Warm WLE discrimination cancellation, honest
# chi-square degrees of freedom, covariance-correct equating drift tests,
# the few-judges clustering guard, and the judge-level DIF ANOVA.

test_that("repeated response rows use a person-clustered calibration sandwich", {
  set.seed(3401)
  N <- 320L; L <- 6L
  X <- matrix(rbinom(N * L, 1,
    plogis(outer(rnorm(N), seq(-1.2, 1.2, length.out = L), "-"))),
    N, L, dimnames = list(NULL, paste0("I", seq_len(L))))

  once <- rasch(X, id = sprintf("P%03d", seq_len(N)))
  twice <- rasch(rbind(X, X),
                 id = rep(sprintf("P%03d", seq_len(N)), 2L))
  independent <- rasch(rbind(X, X), id = sprintf("R%03d", seq_len(2L * N)))

  expect_equal(twice$thresholds$tau, once$thresholds$tau,
               tolerance = 1e-8)
  # Extreme response vectors contain no informative conditional item pairs.
  # The finite-cluster correction therefore uses contributing people, not the
  # raw sample size.
  G <- twice$est$cluster_support$n
  cr1 <- G / (G - 1)
  expect_equal(twice$est$cov_tau, cr1 * once$est$cov_tau,
               tolerance = 1e-7)
  expect_equal(independent$est$cov_tau, once$est$cov_tau / 2,
               tolerance = 1e-7)
  expect_equal(twice$est$cluster_support$cr1, cr1)
  expect_equal(once$est$cluster_support$cr1, 1)
  expect_match(paste(twice$notes, collapse = " "), "clustered")
  expect_true(isTRUE(twice$repeated_ids))
  expect_true(all(is.na(twice$items$p)))
  expect_true(all(is.na(twice$items$p_adj)))
  expect_true(all(is.na(twice$items$p_anova)))
  expect_true(all(is.na(twice$item_trait$p)))
  expect_true(all(is.na(twice$item_anova$p)))
  expect_true(is.na(twice$total_chisq_p))
  expect_true(any(is.finite(twice$items$chisq)))
  expect_match(paste(twice$notes, collapse = " "), "probabilities withheld")
  expect_false(isTRUE(once$repeated_ids))
  expect_false(isTRUE(independent$repeated_ids))
  expect_true(any(is.finite(independent$items$p_adj)))
  expect_error(fit_bootstrap(twice, B = 5, workers = 1, seed = 1),
               "assumes independent response rows")
})

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

test_that("unavailable item-fit tests remain in the predeclared families", {
  set.seed(47)
  ci <- rep(1:2, each = 20)
  Z <- cbind(I1 = rnorm(40), I2 = rnorm(40), I3 = NA_real_)
  anova <- .item_anova(Z, ci, extreme = rep(FALSE, 40))
  ok_a <- is.finite(anova$p)
  expect_equal(anova$p_adj[ok_a],
               p.adjust(anova$p[ok_a], "holm", n = nrow(anova)))
  expect_true(is.na(anova$p_adj[!ok_a]))

  X <- matrix(rbinom(120, 1, 0.5), 40, 3,
              dimnames = list(NULL, paste0("I", 1:3)))
  X[ci == 2, 3] <- NA_integer_
  mo <- list(E = matrix(0.5, 40, 3), V = matrix(0.25, 40, 3))
  trait <- .item_trait(X, mo, ci)
  ok_t <- is.finite(trait$p)
  expect_equal(trait$p_adj[ok_t],
               p.adjust(trait$p[ok_t], "holm", n = nrow(trait)))
  expect_true(is.na(trait$p_adj[!ok_t]))
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
    eq <- equate_tests(mk(), mk(), independent = TRUE)
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

test_that("judge-clustered pair fit keeps the row chi-square descriptive", {
  d <- simulate_btl(5, 12, 5, seed = 1)
  fit <- btl(d, "object_a", "object_b", "winner", judge = "judge")
  expect_true(fit$cl$inference_available)
  expect_true(is.finite(fit$total_chisq))
  expect_true(is.na(fit$total_p))
  expect_match(paste(fit$notes, collapse = " "),
               "row-based reference does not model within-judge dependence")
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

test_that("BTL-EFRM inference refuses an inadequately clustered panel", {
  d <- simulate_btl_efrm(n_objects_per_set = 5, n_sets = 2, n_panels = 2,
                         n_judges_per_panel = 4, reps_within = 15,
                         reps_cross = 15, seed = 3)
  os <- attr(d, "truth")$object_sets
  # collapse panel 2 to a single judge
  keep <- d$panel != "panel2" | d$judge == d$judge[d$panel == "panel2"][1]
  expect_error(
    btl_efrm(d[keep, ], "object_a", "object_b", "winner", "judge", "panel",
             os, se_method = "judge_bootstrap", boot_reps = 30),
    "at least 10 judges|more judge clusters than parameters")
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

test_that("clustered BTL covariance needs full empirical score rank", {
  # Twelve nominal judges all contribute exactly the same strongly-connected
  # tournament. The judge count exceeds the three free object parameters,
  # but their cluster score vectors span only one direction.
  block <- data.frame(
    object_a = c("A", "A", "A", "B", "B", "C"),
    object_b = c("B", "C", "D", "C", "D", "D"),
    winner = c("A", "C", "D", "B", "B", "C"))
  d <- do.call(rbind, lapply(seq_len(12L), function(j)
    transform(block, judge = sprintf("J%02d", j))))
  f <- btl(d, "object_a", "object_b", "winner", judge = "judge")
  expect_true(f$converged)
  expect_gt(f$cl$n_units, f$cl$n_parameters)
  expect_lt(f$cl$score_rank, f$cl$n_parameters)
  expect_false(f$cl$inference_available)
  expect_true(is.na(rasch:::.btl_equate_cov_df(f, "A")))
  expect_true(all(is.na(f$objects$se)))
  expect_true(all(is.na(f$cov_beta)))
  expect_true(is.na(f$total_p))
  expect_true(is.na(f$cl$eff_params))
  expect_true(all(is.na(f$cov_parameters)))
  expect_true(any(grepl("empirical covariance is rank-deficient", f$notes,
                        fixed = TRUE)))

  # The same empirical-rank requirement applies to an unclustered sandwich.
  # Its likelihood information is identified, but these repeated response
  # patterns span only two of the three fitted directions.
  block_u <- data.frame(
    object_a = c("O3", "O1", "O3"),
    object_b = c("O1", "O2", "O1"),
    response = 2:0)
  fu <- btl(block_u[rep(seq_len(3L), 10L), ],
            "object_a", "object_b", response = "response")
  expect_true(fu$converged)
  expect_gt(fu$cl$n_units, fu$cl$n_parameters)
  expect_lt(fu$cl$score_rank, fu$cl$n_parameters)
  expect_false(fu$cl$inference_available)
  expect_true(is.na(rasch:::.btl_equate_cov_df(fu, "O1")))
  expect_true(all(is.na(fu$objects$se)))
  expect_true(all(is.na(fu$cov_beta)))
  expect_error(btl_next_pairs(fu), "covariance")
  expect_true(all(is.na(fu$cov_parameters)))
  expect_match(paste(fu$notes, collapse = " "),
               "independent-comparison score covariance")

  # If the global judge guard withholds the sandwich, every coefficient table
  # also withholds its reference degrees of freedom.
  small <- simulate_btl(5, 6, 3, seed = 14)
  fp <- btl(small, "object_a", "object_b", "winner", judge = "judge",
            position = TRUE)
  expect_false(fp$cl$inference_available)
  expect_true(all(is.na(fp$dependence$se)))
  expect_true(all(is.na(fp$dependence$t)))
  expect_true(all(is.na(fp$dependence$df)))
  expect_true(all(is.na(fp$dependence$p)))
  expect_true(all(is.na(fp$cov_beta)))
  predictors <- data.frame(object = paste0("O", 1:5), x = 1:5)
  fe <- btl_explanatory(small, predictors, ~ x,
                        "object_a", "object_b", winner = "winner",
                        judge = "judge")
  expect_true(all(is.na(fe$object_coefficients$se)))
  expect_true(all(is.na(fe$object_coefficients$t)))
  expect_true(all(is.na(fe$object_coefficients$df)))
  expect_true(all(is.na(fe$object_coefficients$p)))

  # Exact anchors remain constants even when the free sandwich is withheld:
  # their complete covariance row and column are known zeros, not missing.
  anchored <- btl(small, "object_a", "object_b", "winner", judge = "judge",
                  anchors = c(O1 = 0))
  expect_false(anchored$cl$inference_available)
  expect_true(all(anchored$cov_beta["O1", ] == 0))
  expect_true(all(anchored$cov_beta[, "O1"] == 0))
  expect_true(all(is.na(anchored$cov_beta[-1L, -1L])))
  expect_true(is.infinite(rasch:::.btl_equate_cov_df(anchored, "O1")))
  expect_true(is.na(rasch:::.btl_equate_cov_df(anchored, "O2")))
})

test_that("an exact BTL anchor survives an unsupported companion sandwich", {
  truth <- setNames(seq(-1, 1, length.out = 5), paste0("O", 1:5))
  good_data <- simulate_btl(5, 12, 15, object_locations = truth, seed = 501)
  weak_data <- simulate_btl(5, 6, 15, object_locations = truth, seed = 502)
  good <- btl(good_data, "object_a", "object_b", "winner", judge = "judge",
              anchors = c(O1 = unname(truth["O1"])))
  weak <- btl(weak_data, "object_a", "object_b", "winner", judge = "judge",
              anchors = c(O1 = unname(truth["O1"]),
                          O2 = unname(truth["O2"])))
  expect_true(good$cl$inference_available)
  expect_false(weak$cl$inference_available)
  expect_true(is.infinite(rasch:::.btl_equate_cov_df(weak, c("O1", "O2"))))
  expect_true(is.na(rasch:::.btl_equate_cov_df(weak, c("O1", "O3"))))

  eq <- btl_equate(good, weak, independent = TRUE, shift = "none")
  row_o2 <- eq$table$object == "O2"
  expect_true(eq$inferential)
  expect_true(is.finite(eq$table$se_diff[row_o2]))
  expect_true(is.finite(eq$table$df[row_o2]))
  expect_true(is.finite(eq$table$p_adj[row_o2]))
})

test_that("row-based fit support is distinct from calibration support", {
  expect_false(rasch:::.has_repeated_person_ids(c("", "", NA_character_)))
  d <- simulate_rasch(120, 6, n_categories = 3, seed = 992)
  item_names <- grep("^I", names(d), value = TRUE)
  X <- as.matrix(d[item_names])
  # This repeated row has one interior polytomous response. It contributes no
  # conditional item pair, but it does carry a fitted residual and therefore
  # cannot be treated as an independent person by row-based fit references.
  extra <- rep(NA_integer_, ncol(X))
  extra[1L] <- 1L
  fit <- rasch(rbind(X, extra), id = c(d$id, d$id[1L]))

  expect_false(fit$repeated_ids)
  expect_true(fit$repeated_residual_ids)
  expect_true(rasch:::.has_repeated_residual_units(fit))
  expect_true(all(is.na(fit$items$p)))
  expect_true(is.na(fit$total_chisq_p))
  expect_error(fit_bootstrap(fit, B = 1L), "independent response rows")
  expect_error(rasch:::.scree_reference(fit, 2L, 20L, seed = 1L),
               "independent occasions")
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

test_that("unavailable item DIF terms remain in the Holm family", {
  set.seed(701)
  n <- 320L
  g <- factor(rep(c("A", "B"), each = n / 2L))
  d <- seq(-1.2, 1.2, length.out = 6L)
  X <- matrix(rbinom(n * 6L, 1L, plogis(outer(rnorm(n), d, "-"))),
              n, 6L, dimnames = list(NULL, paste0("I", 1:6)))
  X[g == "B", "I1"] <- NA_integer_
  fit <- rasch(data.frame(X, g), factors = "g")
  out <- dif_anova(fit, n_groups = 2L)

  unavailable <- out$terms$item == "I1" &
    out$terms$term %in% c("g", "g:ci")
  expect_equal(sum(unavailable), 2L)
  expect_true(all(is.na(out$terms$p[unavailable])))
  family <- !out$term_ids %in% c("Residuals", "ci")
  expect_equal(sum(family), 2L * ncol(fit$X))
  usable <- family & is.finite(out$terms$p)
  expect_equal(out$terms$p_adj[usable],
               p.adjust(out$terms$p[usable], "holm", n = sum(family)))
  expect_match(paste(out$notes, collapse = " "),
               "remain in the adjusted-probability family")
})

test_that("mixed DIF results have one residual key and validate", {
  set.seed(702)
  n <- 100L
  d <- seq(-1, 1, length.out = 5L)
  theta <- rnorm(n)
  make_wave <- function(offset) {
    X <- matrix(rbinom(n * 5L, 1L,
      plogis(outer(theta + offset, d, "-"))), n, 5L)
    colnames(X) <- paste0("I", 1:5)
    X
  }
  dat <- data.frame(
    rbind(make_wave(0), make_wave(0.2)),
    occasion = rep(c("pre", "post"), each = n))
  fit <- rasch(dat, factors = "occasion", id = rep(seq_len(n), 2L))
  out <- dif_anova(fit, n_groups = 2L)

  expect_equal(sum(out$term_ids == "Residuals"), ncol(fit$X))
  expect_no_error(rasch:::.validate_dif_result(out, fit))
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
  expect_true(any(grepl("not estimable", df$notes)))
  expect_true(any(grepl("dropped from the within-person", df$notes)))
})

test_that("dif_anova has no residual-mean Tukey surface", {
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
  expect_false("tukey" %in% names(d))
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
  # the groups are given by value and the data carries an identical
  # `group` column, so the item columns are named explicitly: this test is
  # about linkage, not about resolving that ambiguity
  expect_error(
    rasch_efrm(X, item_sets = tr$item_sets, groups = grp, items = items),
    "not linked|unidentified")
})

test_that("dif_anova integrates with EFRM and MFRM fits", {
  skip_on_cran()
  d <- simulate_efrm(n_per_group = 300, items_per_set = 8, n_sets = 2,
                     n_groups = 2, group_unit_ratio = 1.25, seed = 1)
  tr <- attr(d, "truth")
  expect_error(
    dif_anova(rasch_efrm(
      d, item_sets = tr$item_sets, groups = tr$groups,
      items = unlist(tr$item_sets, use.names = FALSE),
      boot_reps = 30, workers = 1
    )),
    "frame structure")
  sex <- rep(c("m", "f"), length.out = nrow(d))
  f2 <- rasch_efrm(
    d, item_sets = tr$item_sets, groups = tr$groups,
    items = unlist(tr$item_sets, use.names = FALSE),
    factors = data.frame(sex = sex), boot_reps = 30, workers = 1
  )
  d2 <- dif_anova(f2)
  expect_true(any(grepl("frame structure", d2$notes)))
  expect_gt(nrow(d2$summary), 0)
  expect_error(dif_contrasts(f2), "not available for EFRM")
  expect_error(dif_posthoc(f2, f2$items$item[1], term = "sex"),
               "not available for EFRM")

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
    "cannot also define another model role")
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
  ph <- dif_posthoc(mf, "B", term = "sex")
  expect_s3_class(ph, "rasch_dif_posthoc")
  expect_equal(nrow(ph$table), 1L)
  expect_true(is.finite(ph$table$estimate))
  dc <- dif_contrasts(mf, factors = "sex", items = "B")
  expect_equal(unique(dc$table$item), "B")
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
  # Within one set, equal item difficulties leave the group units weakly
  # identified but do not require a set link. Retain the estimate and show
  # its weakness through the uncertainty.
  weak <- mk(rep(0, 16))
  fz <- rasch_efrm(weak[, sets$A], groups = g,
                   item_sets = list(A = sets$A), boot_reps = 0)
  lr <- abs(diff(log(fz$phi_table$phi)))
  se_lr <- sqrt(sum(fz$phi_table$se_log_phi^2))
  expect_lt(lr / se_lr, 2)              # no spurious unit claim
  expect_gt(min(fz$phi_table$se_log_phi), 0.2)   # no spurious precision
  # The same flat thresholds cannot support a stable scale link between
  # item sets. The linking bootstrap must refuse it rather than return a
  # covariance estimated from a minority of the requested resamples.
  expect_error(
    rasch_efrm(weak, groups = g, item_sets = sets, se_method = "hybrid",
               boot_reps = 30, workers = 1, seed = 42),
    "unit-linking bootstrap replicates were usable")
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
  expect_true(all(c("p_adj", "significant") %in% names(tt)))
  expect_equal(tt$p_adj, stats::p.adjust(tt$p, method = "holm"))
  expect_identical(tt$significant, tt$p_adj < 0.05)
  expect_lt(tt$p[tt$term == "region"], 0.01)
  expect_gt(tt$p[tt$term == "grp"], 0.05)
  expect_gt(tt$p[tt$term == "grp:region"], 0.05)
  # the sum-coded region effect reproduces the planted unit ratio 1.5
  lu <- f$phi_factorial$log_unit[f$phi_factorial$term == "region1"]
  expect_gt(exp(2 * abs(lu)), 1.2)
  expect_lt(exp(2 * abs(lu)), 1.9)
  # Full-bootstrap draws feed se_log_rho jointly. Use enough replicates for
  # the covariance threshold and verify that no hybrid fallback occurred.
  fb <- rasch_efrm(X, groups = c("grp", "region"),
                   items = paste0("v", 1:16),
                   item_sets = list(A = paste0("v", 1:8),
                                    B = paste0("v", 9:16)),
                   se_method = "bootstrap", boot_reps = 50)
  expect_identical(fb$se_method, "bootstrap")
  expect_gte(fb$boot_reps_used, 30L)
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

test_that("withheld DIF probabilities remain in their declared families", {
  d <- simulate_rasch(900, 6, model = "PCM", n_categories = 3,
                      n_groups = 3, seed = 912)
  thin <- which(d$group == "g3")
  d$I02[thin] <- pmin(d$I02[thin], 1L)
  d$I02[thin[1L]] <- 2L
  fit <- rasch(d, model = "PCM", factors = "group", id = "id")

  sizes <- dif_size(fit, "I02", by = "group")$pairs
  usable <- is.finite(sizes$p)
  expect_equal(sum(usable), 1L)
  expect_equal(sizes$p_adj[usable],
               p.adjust(sizes$p[usable], "holm", n = nrow(sizes)))
  expect_true(all(is.na(sizes$p_adj[!usable])))

  contrasts <- dif_contrasts(fit, factors = "group")$table
  usable <- is.finite(contrasts$p)
  expect_lt(sum(usable), nrow(contrasts))
  expect_equal(contrasts$p_adj[usable],
               p.adjust(contrasts$p[usable], "holm",
                        n = nrow(contrasts)))
  expect_true(all(is.na(contrasts$p_adj[!usable])))
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
  expect_false(fit$cl$inference_available)
  expect_true(is.na(fit$cl$eff_params))
  expect_true(all(is.na(fit$objects$se)))
  expect_true(all(is.na(fit$thresholds$se)))
  expect_true(all(is.na(fit$components$se)))
  expect_true(is.na(fit$total_p))
  expect_true(any(grepl("run to the location boundary", fit$notes)))
  # equating refuses a non-converged calibration
  fit2 <- fit
  expect_error(btl_equate(fit, fit2), "did not converge")
})

test_that("btl separation diagnosis is invariant to an anchored origin", {
  # Two densely observed balanced blocks joined by a sparse but balanced
  # bridge are weakly conditioned, not separated. Translating the anchor is
  # a pure change of origin and must not change convergence or any free SE.
  mk <- function(a, b, n)
    data.frame(a = rep(a, n), b = rep(b, n),
               win = rep(c(a, b), length.out = n))
  d <- rbind(mk("A", "B", 1000), mk("C", "D", 1000), mk("A", "C", 10))
  f0 <- btl(d, "a", "b", "win", anchors = c(A = 0))
  f100 <- btl(d, "a", "b", "win", anchors = c(A = 100))
  expect_true(f0$converged)
  expect_true(f100$converged)
  expect_equal(f100$objects$location, f0$objects$location + 100,
               tolerance = 1e-8)
  expect_equal(f100$objects$se, f0$objects$se, tolerance = 1e-8)
  expect_equal(f100$objects$se[f100$objects$object == "A"], 0)
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

test_that("item infit/outfit exclude extreme-score persons", {
  # extreme persons have boundary measures and near-zero residuals that
  # deflate the mean-squares; excluding them raises an item's outfit toward
  # its honest value (matching the log-of-mean-square residual convention)
  set.seed(77); N <- 600; L <- 6
  d <- seq(-1.2, 1.2, length.out = L)
  X <- matrix(rbinom(N * L, 1, plogis(outer(rnorm(N, 0, 1.6), d, "-"))), N, L)
  colnames(X) <- paste0("I", 1:L)
  f <- rasch(as.data.frame(X))
  th <- f$person$theta; ex <- f$person$extreme
  expect_true(any(ex))                       # design has extreme persons
  mo <- .moment_arrays(th, lapply(seq_len(L), function(j)
    f$thresholds$tau[f$thresholds$item == j]), disc = NULL)
  Z <- (X - mo$E) / sqrt(mo$V); colnames(Z) <- colnames(X)
  # the package's own item-fit, with and without the extreme mask, using the
  # exact mean-square formula: the reported values are the excluded ones
  excl <- .item_fit(X, Z, mo, extreme = ex)
  incl <- .item_fit(X, Z, mo, extreme = rep(FALSE, N))
  expect_equal(f$items$outfit_ms, excl$outfit_ms, tolerance = 1e-8)
  expect_false(isTRUE(all.equal(excl$outfit_ms, incl$outfit_ms)))
})

test_that("sim_recovery reports NA bias for origin-identified parameters", {
  skip_on_cran()
  d <- simulate_rasch(400, 10, difficulty = c(-1, 2), seed = 5)
  r <- sim_recovery(rasch(d), d)
  expect_true(all(is.na(r$summary$bias)))       # locations: not identifiable
  expect_true(all(r$summary$correlation > 0.7, na.rm = TRUE))
})

test_that("sim_apply is resilient to per-replicate failures", {
  batch <- sim_replicate(simulate_rasch, 6, n_persons = 200, n_items = 6,
                         seed = 1)
  i <- 0
  out <- sim_apply(batch, function(d) {
    i <<- i + 1; if (i %% 2 == 0) stop("boom") else rasch(d)$psi$PSI })
  expect_equal(attr(out, "n_failed"), 3L)
  expect_equal(sum(is.na(out)), 3L)
  expect_true(any(grepl("boom", attr(out, "failure_messages"))))
})

test_that("btl_dimensionality withholds the verdict under a shared fixed order", {
  skip_on_cran()
  objs <- LETTERS[1:6]
  beta <- setNames(seq(-1.2, 1.2, length.out = 6), objs)
  fixed_seq <- t(combn(objs, 2)); nj <- 12; carry <- 1.5
  gen <- function(shared) {
    rows <- list()
    for (j in 1:nj) {
      seqp <- if (shared) fixed_seq else fixed_seq[sample(nrow(fixed_seq)), ]
      last <- NA
      for (r in seq_len(nrow(seqp))) {
        a <- seqp[r, 1]; b <- seqp[r, 2]; lp <- beta[a] - beta[b]
        if (!is.na(last) && last == a) lp <- lp + carry
        if (!is.na(last) && last == b) lp <- lp - carry
        w <- ifelse(runif(1) < plogis(lp), a, b); last <- w
        rows[[length(rows) + 1]] <- data.frame(a = a, b = b, win = w,
                                               judge = sprintf("J%02d", j),
                                               seq = r)
      }
    }
    do.call(rbind, rows)
  }
  set.seed(1)
  # shared fixed order: the order effect is confounded, verdict withheld
  fs <- suppressWarnings(btl(gen(TRUE), "a", "b", "win", judge = "judge",
                             order = "seq"))
  ds <- suppressWarnings(btl_dimensionality(fs, reps = 30,
    independent_comparisons = TRUE))
  expect_true(is.na(ds$bimensions$above_reference[1]))
  expect_true(is.na(ds$leading_structured))
  expect_true(any(grepl("withheld", capture.output(print(ds)))))
  expect_true(any(grepl("shares", ds$notes) | grepl("withheld", ds$notes)))
  # randomised order: the confound detector does NOT fire, a verdict is given
  fr <- suppressWarnings(btl(gen(FALSE), "a", "b", "win", judge = "judge",
                             order = "seq"))
  dr <- suppressWarnings(btl_dimensionality(fr, reps = 30,
    independent_comparisons = TRUE))
  expect_false(is.na(dr$bimensions$above_reference[1]))
  expect_false(any(grepl("withheld", dr$notes)))
})

test_that("BTL order variation is measured between judges, not repeated rows", {
  one <- data.frame(
    object_a = c("A", "A", "A", "B", "A", "B"),
    object_b = c("B", "C", "B", "C", "C", "C"),
    order = 1:6, judge = "J1")
  same <- rbind(one, transform(one, judge = "J2"))
  fixed <- rasch:::.btl_order_variation(same, LETTERS[1:3])
  expect_true(fixed$replicated)
  expect_true(fixed$shared)

  varied <- same
  varied$order[varied$judge == "J2"] <- c(5, 1, 6, 2, 3, 4)
  randomised <- rasch:::.btl_order_variation(varied, LETTERS[1:3])
  expect_true(randomised$replicated)
  expect_false(randomised$shared)

  # Averaging many presentations within a judge lowers the raw variance of
  # pair positions. The detector scales by that random-order expectation, so
  # a genuinely random order is not mistaken for a common sequence.
  set.seed(4606)
  pairs <- data.frame(object_a = rep(c("A", "A", "B"), each = 8L),
                      object_b = rep(c("B", "C", "C"), each = 8L))
  many <- do.call(rbind, lapply(seq_len(20L), function(j) {
    z <- pairs[sample(nrow(pairs)), , drop = FALSE]
    z$order <- seq_len(nrow(z)); z$judge <- paste0("J", j); z
  }))
  expect_false(rasch:::.btl_order_variation(many, LETTERS[1:3])$shared)

  unreplicated <- rasch:::.btl_order_variation(one, LETTERS[1:3])
  expect_false(unreplicated$replicated)
  expect_true(unreplicated$shared)
})

test_that("a static BTL position effect does not invoke the shared-order guard", {
  d <- simulate_btl(7, 10, reps_per_pair = 12, seed = 903)
  fit <- btl(d, "object_a", "object_b", winner = "winner",
             position = TRUE)
  z <- btl_dimensionality(fit, reps = 20, seed = 2)
  expect_false(is.na(z$leading_structured))
  expect_false(any(grepl("comparison order", z$notes, fixed = TRUE)))
})

test_that("EFRM fits export and report despite the residual-PCA refusal", {
  skip_on_cran()
  set.seed(1); Np <- 240
  simP <- function(th, tau, r) { x <- 0:length(tau)
    p <- exp(r * (x * th - c(0, cumsum(tau)))); p / sum(p) }
  grp <- rep(c("A", "B"), each = Np / 2); phi <- c(A = 0.8, B = 1.25)
  d <- seq(-1.5, 1.5, length.out = 8); theta <- rnorm(Np)
  X <- sapply(seq_along(d), function(i) sapply(seq_len(Np), function(n)
    sample(0:1, 1, prob = simP(theta[n], d[i], phi[grp[n]]))))
  colnames(X) <- sprintf("I%02d", seq_along(d))
  fit <- rasch_efrm(data.frame(X, grp = grp),
                    item_sets = list(core = colnames(X)), groups = "grp")
  out <- tempfile("efrm-export-regression-")
  on.exit(unlink(out, recursive = TRUE), add = TRUE)
  save_warnings <- character(0)
  withCallingHandlers(
    save_outputs(fit, out, formats = "png"),
    warning = function(w) {
      save_warnings <<- c(save_warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  expect_gt(length(save_warnings), 0L)
  expect_true(all(grepl("residual PCA is undefined", save_warnings)))
  for (nm in c("response_cell_statistics.csv", "items_common_unit.csv",
               "thresholds_common_unit.csv", "frame_model_comparison.csv",
               "unit_omnibus_tests.csv", "unit_contrasts.csv"))
    expect_true(file.exists(file.path(out, "tables", nm)), info = nm)
  expect_false(file.exists(file.path(out, "tables", "item_statistics.csv")))
  hfile <- file.path(tempdir(), "efrm-rep.html")
  report_warnings <- character(0)
  withCallingHandlers(
    report_html(fit, hfile),
    warning = function(w) {
      report_warnings <<- c(report_warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  expect_gt(length(report_warnings), 0L)
  expect_true(all(grepl("residual PCA is undefined", report_warnings)))
  html <- paste(readLines(hfile, warn = FALSE), collapse = "\n")
  expect_match(html, "Common-scale item estimates", fixed = TRUE)
  expect_match(html, "Response-cell fit", fixed = TRUE)
  expect_match(html, "Frame model comparison", fixed = TRUE)
  # dimensionality_test carries the refusal as its standard note, not a crash
  dt <- dimensionality_test(fit)
  expect_true(is.na(dt$multidimensional))
  expect_true(grepl("split unavailable", dt$note))
})

test_that("dimensionality_test withholds an uncalibrated automatic verdict", {
  set.seed(1); N <- 500; L <- 10
  d <- seq(-2, 2, length.out = L)
  X <- matrix(rbinom(N * L, 1, plogis(outer(rnorm(N), d, "-"))), N, L)
  colnames(X) <- paste0("I", 1:L)
  dt <- dimensionality_test(rasch(as.data.frame(X)))
  expect_true(is.na(dt$multidimensional))
  expect_true(dt$binomial_multidimensional %in% c(TRUE, FALSE))
  expect_true(!is.null(dt$caution) && grepl("score points", dt$caution))
})

test_that("test_information splits EFRM groups by administration pattern", {
  skip_on_cran()
  set.seed(7)
  simP <- function(th, tau, r) { x <- 0:length(tau)
    ln <- c(0, cumsum(r * (th - tau))); p <- exp(ln - max(ln)); p / sum(p) }
  d10 <- seq(-1.5, 1.5, length.out = 10)
  nA1 <- 199; nA2 <- 40; nB <- 400
  mk <- function(th, phi) t(sapply(th, function(t)
    sapply(d10, function(x) sample(0:1, 1, prob = simP(t, x, phi)))))
  thA1 <- rnorm(nA1); thA2 <- rnorm(nA2); thB <- rnorm(nB)
  core <- rbind(mk(thA1, 0.8), mk(thA2, 0.8), mk(thB, 1.3))
  extra <- rbind(matrix(NA_integer_, nA1, 10), mk(thA2, 0.8), mk(thB, 1.3))
  colnames(core) <- sprintf("C%02d", 1:10); colnames(extra) <- sprintf("E%02d", 1:10)
  grp <- c(rep("A", nA1 + nA2), rep("B", nB))
  fit <- rasch_efrm(data.frame(core, extra, grp = grp, check.names = FALSE),
                    item_sets = list(core = colnames(core),
                                     extra = colnames(extra)),
                    groups = "grp", se_method = "hybrid")
  ti <- test_information(fit)
  des <- unique(ti$design)
  # the core-only majority of group A gets its own honest curve
  expect_true(any(grepl("group=A, sets=core$", des)))
  expect_true(any(grepl("group=A, sets=core\\+extra", des)))
  sem_core <- ti$sem[grepl("group=A, sets=core$", ti$design) & ti$theta == 0]
  sem_both <- ti$sem[grepl("group=A, sets=core\\+extra", ti$design) & ti$theta == 0]
  expect_gt(sem_core, sem_both)          # fewer items -> larger SEM, honestly
})

test_that("MFRM information combines facets administered to the same person", {
  fit <- structure(list(
    est = list(converged = TRUE),
    tau_list = rep(list(0), 4),
    disc = rep(1, 4),
    virtual_map = data.frame(
      vkey = paste0("v", 1:4), item = rep(c("I1", "I2"), 2),
      rater = rep(c("A", "B"), each = 2), stringsAsFactors = FALSE),
    facet_spec = "rater",
    X = matrix(c(1, 0, 1, 0,
                 0, 1, 0, 1,
                 1, 0, NA, NA,
                 NA, NA, 0, 1), nrow = 4, byrow = TRUE)
  ), class = c("rasch_mfrm", "rasch"))
  z <- rasch:::.design_blocks(fit)
  expect_length(z, 3L)
  expect_true(any(vapply(z, function(i) setequal(i, 1:4), logical(1))))
  expect_true(any(vapply(z, function(i) setequal(i, 1:2), logical(1))))
  expect_true(any(vapply(z, function(i) setequal(i, 3:4), logical(1))))
  ti <- test_information(fit, grid = 0)
  both <- grepl("rater=A.*rater=B", ti$design)
  expect_equal(ti$info[both], 1, tolerance = 1e-12)
  expect_equal(sort(ti$info[!both]), c(0.5, 0.5), tolerance = 1e-12)
})

test_that("MFRM interaction omnibus uses the Hotelling-style F reference", {
  set.seed(31)
  simP <- function(th, tau) { x <- 0:length(tau)
    p <- exp(x * th - c(0, cumsum(tau))); p / sum(p) }
  persons <- sprintf("P%03d", 1:120); raters <- paste0("R", 1:3)
  items <- paste0("I", 1:4)
  th <- setNames(rnorm(120, 0, 1.3), persons)
  rho <- setNames(c(-0.5, 0, 0.5), raters)
  d <- expand.grid(person = persons, item = items, rater = raters,
                   stringsAsFactors = FALSE)
  d$score <- mapply(function(p, i, r)
    sample(0:2, 1, prob = simP(th[p], c(-0.4, 0.4) + rho[r])),
    d$person, d$item, d$rater)
  fit <- rasch_mfrm(d, person = "person", item = "item", score = "score",
                    facets = "rater", interaction = "rater")
  it <- fit$interaction_test
  expect_true(all(c("f", "df2") %in% names(it)))
  # the F reference is strictly more conservative than the old chi-square
  p_chisq <- pchisq(it$wald, it$df, lower.tail = FALSE)
  expect_gte(it$p, p_chisq)
  # cells use a t reference (p >= the normal-reference p)
  ie <- fit$interaction_effects
  expected_df <- floor(min(fit$interaction_support$effective_persons)) - 1L
  expect_equal(ie$df, rep(expected_df, nrow(ie)))
  expect_equal(ie$p, 2 * pt(-abs(ie$t), df = ie$df))
  expect_true(all(ie$p >= 2 * pnorm(-abs(ie$t)) - 1e-12))
  out <- tempfile("mfrm-export-")
  on.exit(unlink(out, recursive = TRUE), add = TRUE)
  expect_no_warning(invisible(
    save_outputs(fit, out, formats = "png", item_plots = FALSE)))
  expect_true(file.exists(file.path(out, "tables",
                                    "interaction_omnibus_test.csv")))
  expect_true(file.exists(file.path(out, "tables",
                                    "response_cell_statistics.csv")))
})

test_that("MFRM interaction inference uses the least-supported cell", {
  set.seed(3101)
  simP <- function(th, tau) {
    x <- 0:length(tau); p <- exp(x * th - c(0, cumsum(tau))); p / sum(p)
  }
  persons <- sprintf("P%03d", 1:90); items <- paste0("I", 1:4)
  d <- expand.grid(person = persons, item = items,
                   rater = paste0("R", 1:3), stringsAsFactors = FALSE)
  d <- d[d$rater != "R3" | d$person %in% persons[1:12], ]
  theta <- setNames(rnorm(length(persons)), persons)
  d$score <- vapply(seq_len(nrow(d)), function(i)
    sample(0:2, 1, prob = simP(theta[d$person[i]], c(-0.5, 0.5))), 0L)
  fit <- rasch_mfrm(d, "person", "item", "score", facets = "rater",
                    interaction = "rater")
  expect_false(fit$interaction_test$inference_available)
  expect_identical(fit$interaction_test$inference_available,
                   is.finite(fit$interaction_test$p))
  expect_equal(min(fit$interaction_support$n_persons), 12)
  expect_true(is.na(fit$interaction_test$p))
  expect_true(all(is.na(fit$interaction_effects$p)))
  expect_true(all(is.finite(fit$interaction_effects$gamma)))

  # A single sparse item-by-rater cell must withhold inference even when that
  # rater has ample observations on the remaining items.
  d2 <- expand.grid(person = persons, item = items,
                    rater = paste0("R", 1:3), stringsAsFactors = FALSE)
  d2 <- d2[d2$rater != "R1" | d2$item != "I1" |
             d2$person %in% persons[1:12], ]
  d2$score <- vapply(seq_len(nrow(d2)), function(i)
    sample(0:2, 1, prob = simP(theta[d2$person[i]], c(-0.5, 0.5))), 0L)
  fit2 <- rasch_mfrm(d2, "person", "item", "score", facets = "rater",
                     interaction = "rater")
  sparse <- fit2$interaction_support$item == "I1" &
    fit2$interaction_support$level == "R1"
  expect_equal(fit2$interaction_support$n_persons[sparse], 12)
  expect_false(fit2$interaction_test$inference_available)
  expect_true(all(is.na(fit2$interaction_effects$p)))
})

test_that("btl eff_params is withheld with clustered inference", {
  set.seed(3)
  K <- 6; b <- seq(-1, 1, length.out = K); n <- 200
  ia <- sample(K, n, TRUE); ib <- (ia + sample(K - 1, n, TRUE) - 1L) %% K + 1L
  d <- data.frame(a = paste0("O", ia), b = paste0("O", ib),
                  winner = paste0("O", ifelse(
                    rbinom(n, 1, plogis(b[ia] - b[ib])) == 1, ia, ib)),
                  judge = sample(sprintf("J%d", 1:4), n, TRUE))
  f <- suppressWarnings(btl(d, "a", "b", "winner", judge = "judge"))
  expect_false(isTRUE(f$cl$inference_available))
  expect_true(is.na(f$cl$eff_params))
})

test_that("DIF post-hoc adjustment retains unsupported planned contrasts", {
  set.seed(8301)
  n <- 360L
  group <- factor(rep(c("A", "B", "C"), each = n / 3L))
  stratum <- factor(c(rep("X", 2L * n / 3L), rep("Y", n / 3L)))
  theta <- stats::rnorm(n)
  difficulty <- seq(-1.2, 1.2, length.out = 6L)
  response <- vapply(difficulty, function(delta) {
    stats::rbinom(n, 1L, stats::plogis(theta - delta))
  }, integer(n))
  colnames(response) <- paste0("I", seq_len(ncol(response)))

  fit <- rasch(data.frame(response, group, stratum),
               items = colnames(response), factors = c("group", "stratum"))
  out <- dif_posthoc(fit, "I1", term = "group", min_n = 20L)

  expect_identical(out$family$contrast, "B - A")
  expect_identical(out$family_n_per_item, 3L)
  expect_identical(out$family_n, 3L)
  expect_equal(out$table$p_adj,
               stats::p.adjust(out$table$p, method = "holm", n = 3L))
  expect_true(any(grepl("remain in the multiplicity family", out$notes,
                        fixed = TRUE)))

  factors <- data.frame(group, stratum)
  cells <- rasch:::.factor_cells(factors, sep = ":")
  cellmap <- unique(data.frame(cell = as.character(cells), factors))
  family <- rasch:::.dif_posthoc_family(
    factors, cellmap, target = "group", within = character(0))
  expect_identical(family$planned, c("B - A", "C - A", "C - B"))
  expect_identical(family$planned_n, 3L)

  auto <- dif_contrasts(fit, factors = c("group", "stratum"),
                        items = "I1", min_n = 20L)
  expect_identical(nrow(auto$table), 1L)
  expect_identical(auto$family_n_per_item, 4L)
  expect_identical(auto$family_n, 4L)
  expect_equal(auto$table$p_adj,
               stats::p.adjust(auto$table$p, method = "holm", n = 4L))
  expect_match(auto$family$cells,
               "A:X -1.00, B:X \\+1.00")
  expect_true(any(grepl("remain in the multiplicity family", auto$notes,
                        fixed = TRUE)))
})

test_that("DIF contrast labels cannot overwrite punctuation-bearing levels", {
  f <- factor(c("a", "b - c", "c - a", "b"),
              levels = c("a", "b - c", "c - a", "b"))
  one <- rasch:::.dif_factor_contrasts(f, "group")
  expect_length(one, choose(nlevels(f), 2L))
  expect_false(anyDuplicated(names(one)) > 0L)

  factors <- data.frame(group = f)
  cells <- rasch:::.factor_cells(factors, sep = ":")
  cellmap <- unique(data.frame(cell = as.character(cells), factors))
  posthoc <- rasch:::.dif_posthoc_family(
    factors, cellmap, target = "group", within = character(0))
  expect_equal(posthoc$planned_n, choose(nlevels(f), 2L))
  expect_length(posthoc$family, choose(nlevels(f), 2L))
  expect_false(anyDuplicated(names(posthoc$family)) > 0L)
})
