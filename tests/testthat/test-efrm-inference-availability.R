test_that("failed EFRM stage one withholds every exposed covariance", {
  set.seed(1)
  n <- 240L
  X <- matrix(rbinom(n * 4L, 1L, plogis(rnorm(n))), n, 4L)
  colnames(X) <- paste0("I", 1:4)
  d <- data.frame(X, group = rep(c("A", "B"), each = n / 2L))
  expect_warning(f <- rasch_efrm(
    d, list(all = colnames(X)), "group", boot_reps = 0, maxit = 1L),
    "did NOT converge")
  expect_false(f$est$stage1_converged)
  expect_false(f$est$converged)
  expect_false(f$unit_support$phi_inference)
  expect_false(f$unit_support$alpha_inference)
  expect_true(all(is.na(f$est$cov_tau)))
  expect_true(all(is.na(f$est$thr$se)))
  for (nm in setdiff(names(f$unit_cov), "method"))
    if (!is.null(f$unit_cov[[nm]]))
      expect_true(all(is.na(f$unit_cov[[nm]])), info = nm)
  expect_true(all(is.na(f$score_curves$sem)))
  expect_true(all(is.finite(f$score_curves$expected_score)))
  expect_match(paste(f$notes, collapse = " "),
               "probabilities are withheld because the conditional calibration did not converge")
})

test_that("converged weak EFRM group units remain descriptive", {
  # With almost no item spread, the relative unit can be estimated but is
  # much too uncertain for its normal/Wald reference. This ordinary fit
  # previously warned about weak identification while still reporting p.
  set.seed(16)
  n <- 2000L
  X <- matrix(rbinom(n * 2L, 1L, 0.5), n, 2L)
  colnames(X) <- c("I1", "I2")
  d <- data.frame(X, group = rep(c("A", "B"), each = n / 2L))
  expect_warning(f <- rasch_efrm(
    d, list(all = colnames(X)), "group", boot_reps = 0),
    "only weakly identified")
  expect_true(f$est$converged)
  expect_true(all(f$phi_table$se_log_phi > 5))
  expect_true(all(f$unit_support$group$n_persons >= 50))
  expect_false(f$unit_support$phi_inference)
  expect_true(all(is.finite(f$unit_cov$cov_log_phi)))
  # The one-set alpha is fixed at one, not estimated through a group link.
  expect_true(f$unit_support$alpha_inference)
  expect_true(all(is.finite(f$phi_table$phi)))
  for (part in c("unit_tests", "unit_omnibus")) {
    expect_true(all(is.na(f$efrm_vs_rasch[[part]]$p)))
    expect_true(all(is.na(f$efrm_vs_rasch[[part]]$p_adj)))
    expect_true(all(is.na(f$efrm_vs_rasch[[part]]$significant)))
  }
  interval_drawn <- FALSE
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  testthat::with_mocked_bindings(
    plot_frames(f),
    segments = function(...) interval_drawn <<- TRUE,
    .package = "rasch")
  expect_false(interval_drawn)
})

test_that("EFRM weak-unit decisions use the returned bootstrap SE", {
  d <- simulate_efrm(n_per_group = 150, items_per_set = 5,
                     n_sets = 1, n_groups = 2, seed = 4)
  original_solve <- rasch:::.efrm_solve
  # Only the intermediate reported SE changes; all estimates and covariance
  # matrices, including the actual full-bootstrap draws, remain unmodified.
  testthat::local_mocked_bindings(
    .efrm_solve = function(...) {
      z <- original_solve(...)
      z$se_log_phi[] <- 6
      z
    }, .package = "rasch")
  expect_no_warning(f <- rasch_efrm(
    d, attr(d, "truth")$item_sets, "group", se_method = "bootstrap",
    boot_reps = 30, workers = 1, seed = 41))
  expect_identical(f$se_method, "bootstrap")
  expect_true(f$est$converged)
  expect_true(all(f$phi_table$se_log_phi < 5))
  expect_true(f$unit_support$phi_inference)
  expect_true(all(is.finite(f$efrm_vs_rasch$unit_tests$p)))
  expect_false(any(grepl("weakly identified unit", f$notes)))
})

test_that("EFRM set inference follows the group units used in its link", {
  d <- simulate_efrm(n_per_group = 200, items_per_set = 4,
                     n_sets = 2, n_groups = 2, seed = 93)
  original_solve <- rasch:::.efrm_solve
  testthat::local_mocked_bindings(
    .efrm_solve = function(...) {
      z <- original_solve(...)
      z$se_log_phi[] <- 6
      z
    }, .package = "rasch")
  expect_warning(f <- rasch_efrm(
    d, attr(d, "truth")$item_sets, "group", boot_reps = 40,
    workers = 1, seed = 49), "only weakly identified")
  expect_true(f$est$converged)
  expect_true(all(f$unit_support$set$n_common_persons >= 50))
  expect_true(all(is.finite(f$alpha_table$se_log_alpha)))
  expect_true(all(is.finite(f$unit_cov$cov_log_alpha)))
  expect_false(f$unit_support$phi_inference)
  expect_false(f$unit_support$alpha_inference)
  expect_true(all(is.na(f$efrm_vs_rasch$unit_omnibus$p_adj)))
  expect_true(all(is.na(f$efrm_vs_rasch$unit_tests$p_adj)))
  expect_match(paste(f$notes, collapse = " "), "their links depend on")
})

test_that("EFRM failed set links retain only stage-one uncertainty", {
  d <- simulate_efrm(n_per_group = 120, items_per_set = 6,
                     n_sets = 2, n_groups = 1, seed = 1901)
  original_mass <- rasch:::efrm_fit_weights_cpp
  testthat::local_mocked_bindings(
    efrm_fit_weights_cpp = function(...) {
      z <- original_mass(...)
      z$converged <- FALSE
      z
    }, .package = "rasch")
  expect_warning(f <- rasch_efrm(
    d, attr(d, "truth")$item_sets, "group", boot_reps = 0, workers = 1),
    "nuisance masses")
  expect_true(f$est$stage1_converged)
  expect_false(f$est$converged)
  expect_true(all(is.finite(f$unit_cov$cov_dtilde)))
  expect_true(all(is.finite(f$unit_cov$cov_log_phi)))
  expect_true(all(is.na(f$unit_cov$cov_delta)))
  expect_true(all(is.na(f$score_curves$sem)))
  expect_true(all(is.finite(f$score_curves$expected_score)))
})

test_that("EFRM unit contrasts count a centred pair as one hypothesis", {
  d <- simulate_efrm(n_per_group = 250, items_per_set = 5, n_sets = 2,
                     n_groups = 2, set_unit_ratio = 1.6,
                     group_unit_ratio = 1.3, seed = 3)
  f <- rasch_efrm(d, attr(d, "truth")$item_sets, "group", id = "id",
                  boot_reps = 0, workers = 1)
  ut <- f$efrm_vs_rasch$unit_tests
  om <- f$efrm_vs_rasch$unit_omnibus
  expect_identical(nrow(ut), 4L)
  # The units are centred, so the two group rows state one hypothesis and
  # the two set rows another. The second row of a pair restates the first:
  # it keeps its own unadjusted probability, but the adjustment and the flag
  # are withheld, so one difference is not reported as two deviating units.
  expect_equal(ut$p[1], ut$p[2])
  expect_true(is.na(ut$p_adj[2]))
  expect_true(is.na(ut$significant[2]))
  expect_true(isTRUE(ut$significant[1]))
  expect_identical(sum(ut$significant, na.rm = TRUE), 1L)
  expect_match(paste(f$notes, collapse = " | "),
               "two centred unit rows are one hypothesis", fixed = TRUE)
  # The set link is unavailable at boot_reps = 0, but it stays a declared
  # member, so the follow-up family is the two the omnibus declares and the
  # two tables adjust the same hypothesis by the same rule.
  expect_true(all(is.na(ut$p[3:4])))
  expect_equal(ut$p_adj[1], min(1, 2 * ut$p[1]))
  expect_equal(ut$p_adj[1], om$p_adj[om$term == "group units (phi)"])
  expect_true(all(is.na(ut$p_adj[3:4])))
  # An unavailable probability has no reference distribution to report.
  expect_identical(ut$df, c(Inf, Inf, NA_real_, NA_real_))
})

test_that("three EFRM group units stay three members of the unit family", {
  d <- simulate_efrm(n_per_group = 150, items_per_set = 5, n_sets = 1,
                     n_groups = 3, group_unit_ratio = 1.3, seed = 11)
  f <- rasch_efrm(d, attr(d, "truth")$item_sets, "group", id = "id",
                  boot_reps = 0, workers = 1)
  ut <- f$efrm_vs_rasch$unit_tests
  expect_identical(nrow(ut), 3L)
  # Beyond two groups no two rows are the same hypothesis, so each reported
  # contrast is its own member.
  expect_equal(ut$p_adj, stats::p.adjust(ut$p, "holm"))
})

test_that("EFRM unit contrasts refer a bootstrap standard error to t", {
  d <- simulate_efrm(n_per_group = 150, items_per_set = 5, n_sets = 2,
                     n_groups = 2, set_unit_ratio = 1.5,
                     group_unit_ratio = 1.3, seed = 17)
  f <- rasch_efrm(d, attr(d, "truth")$item_sets, "group", id = "id",
                  se_method = "hybrid", boot_reps = 40, workers = 1, seed = 5)
  ut <- f$efrm_vs_rasch$unit_tests
  # The group units keep the analytic sandwich, so they keep the normal
  # reference; the set link is a person bootstrap, and the standard
  # deviation of its retained draws has B - 1 degrees of freedom.
  expect_identical(ut$df[1:2], c(Inf, Inf))
  expect_equal(unique(ut$df[3:4]),
               f$linking$boot_reps_used - 1)
  expect_equal(ut$p, 2 * stats::pt(-abs(ut$z), ut$df))
  expect_gt(ut$p[3], 2 * stats::pnorm(-abs(ut$z[3])))
  # Two available centred pairs are two hypotheses, not four: each pair is
  # adjusted once and the row that restates it is withheld.
  expect_equal(ut$p_adj[c(1, 3)], p.adjust(ut$p[c(1, 3)], "holm"))
  expect_true(all(is.na(ut$p_adj[c(2, 4)])))
  expect_true(all(is.na(ut$significant[c(2, 4)])))
  # The omnibus on the same estimated covariance is referred to F, so a
  # one-dimensional omnibus and its contrast report one probability rather
  # than two for the same hypothesis.
  om <- f$efrm_vs_rasch$unit_omnibus
  expect_identical(om$df, c(1L, 1L))
  expect_identical(om$df2[1], Inf)
  expect_equal(om$df2[2], f$linking$boot_reps_used - 1)
  expect_equal(om$p, ut$p[c(1, 3)])
})

test_that("an analytic EFRM omnibus covariance keeps the chi-square reference", {
  d <- simulate_efrm(n_per_group = 150, items_per_set = 5, n_sets = 1,
                     n_groups = 3, group_unit_ratio = 1.3, seed = 11)
  f <- rasch_efrm(d, attr(d, "truth")$item_sets, "group", id = "id",
                  boot_reps = 0, workers = 1)
  om <- f$efrm_vs_rasch$unit_omnibus
  expect_identical(om$df2, Inf)
  expect_equal(om$f, om$wald / om$df)
  expect_equal(om$p, pchisq(om$wald, om$df, lower.tail = FALSE))
})

test_that("the crossed unit decomposition shares the omnibus reference", {
  d <- simulate_efrm(n_per_group = 140, items_per_set = 6, n_sets = 1,
                     n_groups = 2, seed = 5120)
  d$band <- rep(c("N", "S"), length.out = nrow(d))
  sets <- attr(d, "truth")$item_sets
  # An analytic cell-unit covariance keeps the asymptotic chi-square.
  fa <- rasch_efrm(d, sets, c("group", "band"), id = "id", boot_reps = 0,
                   workers = 1)
  ta <- fa$phi_factorial_tests
  expect_identical(ta$df2, rep(Inf, nrow(ta)))
  expect_equal(ta$p, pchisq(ta$wald, ta$df, lower.tail = FALSE))
  expect_identical(fa$efrm_vs_rasch$unit_omnibus$df2, Inf)
  # A full bootstrap estimates that covariance instead, so the crossed
  # decomposition and the unit omnibus drawn from the same resamples are
  # both referred to F. One fit must not report two references for two
  # tests on one set of draws.
  fb <- rasch_efrm(d, sets, c("group", "band"), id = "id",
                   se_method = "bootstrap", boot_reps = 60, workers = 1,
                   seed = 9)
  expect_identical(fb$se_method, "bootstrap")
  tb <- fb$phi_factorial_tests
  B <- fb$boot_reps_used
  expect_equal(tb$df2, B - tb$df)
  expect_equal(tb$p, pf(tb$f, tb$df, tb$df2, lower.tail = FALSE))
  expect_true(all(tb$p > pchisq(tb$wald, tb$df, lower.tail = FALSE)))
  ob <- fb$efrm_vs_rasch$unit_omnibus
  expect_equal(ob$df2, B - ob$df)
  expect_true(all(is.finite(ob$f)))
})

test_that("a one-group EFRM still reports one adjusted set-unit probability", {
  # The shape the shipped case study fits: no group units, two centred set
  # units. The whole follow-up family is then the single set-unit hypothesis,
  # so the row that restates it is withheld and the reduction the case study
  # applies to that column must still find a probability to reduce.
  d <- simulate_efrm(n_per_group = 200, items_per_set = 5, n_sets = 2,
                     n_groups = 1, set_unit_ratio = 1.5, seed = 23)
  f <- rasch_efrm(d, attr(d, "truth")$item_sets, "group", id = "id",
                  se_method = "hybrid", boot_reps = 60, workers = 1, seed = 4)
  ut <- f$efrm_vs_rasch$unit_tests
  expect_identical(nrow(ut), 2L)
  expect_equal(ut$p[1], ut$p[2])
  # One declared member leaves Holm nothing to adjust against.
  expect_equal(ut$p_adj[1], ut$p[1])
  expect_true(is.na(ut$p_adj[2]))
  expect_true(is.finite(min(ut$p_adj, na.rm = TRUE)))
  expect_identical(unique(ut$df), f$linking$boot_reps_used - 1)
  expect_equal(f$efrm_vs_rasch$unit_omnibus$p, ut$p[1])
})
