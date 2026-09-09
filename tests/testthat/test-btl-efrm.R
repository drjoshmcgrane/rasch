# fast path for tests whose subject is not the standard errors: the
# conditional path avoids bootstrapping; the point estimates are identical
befit <- function(...) btl_efrm(..., se_method = "conditional")

test_that("panel-unit reconciliation uses scale-free precision weights", {
  block <- function(value, variance) list(
    ref = "A", free = "B", lrho = c(B = value),
    cov = matrix(variance, 1L, 1L, dimnames = list("B", "B")))
  z <- rasch:::.btlef_reconcile_phi(
    c("A", "B"), list(block(1, 1e-14), block(3, 4e-14)))
  # inverse-variance mean = (1 + 3 / 4) / (1 + 1 / 4) = 1.4;
  # centring changes the two origins but not their difference
  expect_equal(unname(z$lphi["B"] - z$lphi["A"]), 1.4,
               tolerance = 1e-10)
  expect_error(rasch:::.btlef_reconcile_phi(
    c("A", "B"), list(block(1, -1))), "positive definite")
  asym_block <- list(
    ref = "A", free = c("B", "C"), lrho = c(B = 1, C = 2),
    cov = matrix(c(1, 0.5, 0, 1), 2L,
                 dimnames = list(c("B", "C"), c("B", "C"))))
  expect_error(rasch:::.btlef_reconcile_phi(
    c("A", "B", "C"), list(asym_block)), "asymmetric")
  bad <- rasch:::.btlef_wald_unit(
    c(0.2, -0.2), diag(c(1, -0.5)), "unit")
  expect_true(is.na(bad$wald))
  expect_true(is.na(bad$p))
  asymmetric <- matrix(c(1, 0.5, 0, 1), 2L)
  bad <- rasch:::.btlef_wald_unit(c(0.2, -0.2), asymmetric, "unit")
  expect_true(is.na(bad$wald))
})

test_that("panel reconciliation carries cross-set judge covariance", {
  block <- function(value, influence) {
    influence <- matrix(influence, ncol = 1L,
                        dimnames = list(names(influence), "B"))
    list(ref = "A", free = "B", lrho = c(B = value),
         cov = crossprod(influence), influence_lrho = influence)
  }
  base <- c(J1 = 0.5, J2 = -0.5, J3 = 0.5, J4 = -0.5)

  # The established point estimate uses only the two within-set variances:
  # weights .8 and .2 give a log-ratio of 1.4 in every case below. With all
  # judges shared, the joint observation covariance is singular, but it is
  # propagated rather than inverted. Cross-covariance is 2, so the fixed
  # estimator has variance .8 + 2(.8)(.2)(2) = 1.44 before centring.
  shared <- rasch:::.btlef_reconcile_phi(c("A", "B"), list(
    block(1, base), block(3, 2 * base)))
  expect_equal(unname(shared$lphi["B"] - shared$lphi["A"]), 1.4,
               tolerance = 1e-12)
  expect_equal(unname(shared$se_log_phi), rep(sqrt(1.44) / 2, 2),
               tolerance = 1e-12)

  # Name alignment uses only J3/J4 here: their product is 1. The remaining
  # judges belong to only one set and therefore add no cross-set covariance.
  partial2 <- c(J3 = 1, J4 = -1, J5 = 1, J6 = -1)
  partial <- rasch:::.btlef_reconcile_phi(c("A", "B"), list(
    block(1, base), block(3, partial2)))
  expect_equal(unname(partial$se_log_phi), rep(sqrt(1.12) / 2, 2),
               tolerance = 1e-12)

  # Disjoint judge pools reproduce the old block-diagonal covariance.
  disjoint2 <- c(K1 = 1, K2 = -1, K3 = 1, K4 = -1)
  disjoint <- rasch:::.btlef_reconcile_phi(c("A", "B"), list(
    block(1, base), block(3, disjoint2)))
  expect_equal(unname(disjoint$se_log_phi), rep(sqrt(0.8) / 2, 2),
               tolerance = 1e-12)

  # Negative cross-set dependence must reduce, not inflate, the same
  # estimator's covariance. Reversing the partially shared contributions
  # gives cross-covariance -1 and variance .48 before centring.
  negative <- rasch:::.btlef_reconcile_phi(c("A", "B"), list(
    block(1, base), block(3, -partial2)))
  expect_equal(unname(negative$se_log_phi), rep(sqrt(0.48) / 2, 2),
               tolerance = 1e-12)
})

test_that("conditional BTL-EFRM SE retains shared-judge covariance", {
  set.seed(831)
  d <- simulate_btl_efrm(n_objects_per_set = 6, n_sets = 2,
                         n_judges_per_panel = 30, n_panels = 2,
                         reps_within = 30, reps_cross = 10, seed = 831)
  truth <- attr(d, "truth")
  # A persistent judge discrimination creates visible covariance between the
  # separately estimated panel ratios without changing the sampling unit.
  q <- setNames(exp(rnorm(length(unique(d$judge)), 0, 0.85)), unique(d$judge))
  eta <- q[d$judge] * (truth$v[d$object_a] - truth$v[d$object_b])
  ya <- rbinom(nrow(d), 1L, plogis(eta))
  d$winner <- ifelse(ya == 1L, d$object_a, d$object_b)

  fit <- btl_efrm(d, "object_a", "object_b", winner = "winner",
                  judge = "judge", panels = "panel",
                  object_sets = truth$object_sets,
                  se_method = "conditional")
  blocks <- lapply(names(truth$object_sets), function(s) {
    os <- truth$object_sets[[s]]
    rows <- d$object_a %in% os & d$object_b %in% os
    f <- rasch:::.btlef_stage1(
      match(d$object_a[rows], os), match(d$object_b[rows], os),
      as.integer(d$winner[rows] == d$object_a[rows]),
      as.character(d$panel[rows]), as.character(d$judge[rows]),
      length(os), 60, 1e-8)
    expect_equal(crossprod(f$influence_lrho), f$cov_lrho,
                 tolerance = 1e-12)
    list(ref = f$ref, free = f$free, lrho = log(f$rho[f$free]),
         cov = f$cov_lrho, influence_lrho = f$influence_lrho)
  })
  joint <- rasch:::.btlef_reconcile_phi(c("panel1", "panel2"), blocks)
  independent <- rasch:::.btlef_reconcile_phi(
    c("panel1", "panel2"), lapply(blocks, function(z) {
      z$influence_lrho <- NULL
      z
    }))

  expect_equal(log(fit$phi_table$phi), unname(joint$lphi), tolerance = 1e-10)
  expect_equal(fit$phi_table$se_log_phi, unname(joint$se_log_phi),
               tolerance = 1e-10)
  expect_equal(joint$lphi, independent$lphi, tolerance = 1e-12)
  expect_gt(min(joint$se_log_phi - independent$se_log_phi), 0.01)
  expect_equal(round(unname(independent$se_log_phi[1]), 4), 0.1011)
  expect_equal(round(unname(joint$se_log_phi[1]), 4), 0.1189)
})

# Extended frame of reference for paired comparisons: reduction to btl(),
# recovery of panel and set units, null calibration, and the guards.

test_that("G = 1, S = 1 reduces exactly to btl()", {
  d <- simulate_btl_efrm(n_objects_per_set = 7, n_sets = 1, n_panels = 1,
                         reps_within = 40, seed = 1)
  fit <- befit(d, "object_a", "object_b", winner = "winner", judge = "judge",
                  panels = "panel", object_sets = attr(d, "truth")$object_sets)
  expect_s3_class(fit, "rasch_btl_efrm")
  bt <- btl(d, "object_a", "object_b", winner = "winner")
  loc_ef <- fit$objects$beta_set[match(bt$objects$object, fit$objects$object)]
  expect_equal(loc_ef, bt$objects$location, tolerance = 1e-6)
  # a single frame: phi = 1, common-scale value equals the frame location, and
  # there is no set unit beyond the reference row
  expect_equal(fit$phi_table$phi, 1)
  expect_equal(fit$objects$v, fit$objects$beta_set)
  expect_equal(nrow(fit$alpha_table), 1L)
  expect_equal(fit$alpha_table$alpha, 1)
  expect_error(btl_dif(fit, factors = rep("g", nrow(fit$comparisons))),
               "not defined after a BTL-EFRM frame adjustment")
})

test_that("a non-converged BTL-EFRM fit retains no inferential fields", {
  d <- simulate_btl_efrm(n_objects_per_set = 8, n_sets = 1, n_panels = 1,
                         reps_within = 40, seed = 1)
  expect_warning(
    fit <- btl_efrm(d, "object_a", "object_b", winner = "winner",
                    judge = "judge", panels = "panel",
                    object_sets = attr(d, "truth")$object_sets,
                    se_method = "conditional", maxit = 2),
    "did NOT converge")
  expect_false(fit$converged)
  expect_true(all(is.na(fit$objects$se)))
  expect_true(all(is.na(fit$objects$se_beta)))
  expect_true(all(is.na(fit$phi_table$se_log_phi)))
  expect_true(all(is.na(fit$phi_table$p)))
  expect_true(all(is.na(fit$phi_table$p_adj)))
  expect_true(is.na(fit$total_p))
  expect_true(is.na(fit$osi$PSI))
  expect_null(fit$cov_beta)
  expect_true(is.na(fit$equal_unit$two_delta_ll))
  expect_match(paste(fit$notes, collapse = " "), "probabilities withheld")
})

test_that("a failed equal-unit refit cannot supply a likelihood difference", {
  d <- simulate_btl_efrm(n_objects_per_set = 6, n_sets = 2, n_panels = 2,
                         n_judges_per_panel = 10, reps_within = 20,
                         reps_cross = 20, seed = 381)
  original <- rasch:::.btlef_stage1
  testthat::local_mocked_bindings(
    .btlef_stage1 = function(ia, ib, y, panel, judge, K, maxit, tol,
                             rho_fixed = NULL) {
      z <- original(ia, ib, y, panel, judge, K, maxit, tol,
                    rho_fixed = rho_fixed)
      if (length(unique(panel)) == 1L && identical(unique(panel), "all"))
        z$converged <- FALSE
      z
    },
    .package = "rasch"
  )
  fit <- btl_efrm(d, "object_a", "object_b", winner = "winner",
                  judge = "judge", panels = "panel",
                  object_sets = attr(d, "truth")$object_sets,
                  se_method = "conditional")
  expect_true(fit$converged)
  expect_true(is.na(fit$equal_unit$loglik_single))
  expect_true(is.na(fit$equal_unit$difference))
  expect_true(is.na(fit$equal_unit$two_delta_ll))
  expect_match(fit$equal_unit$note, "unavailable")
})

test_that("conditional panel units are recovered but inference is withheld", {
  phi_true <- c(0.7, 1.0, 1.43); phi_true <- phi_true / exp(mean(log(phi_true)))
  d <- simulate_btl_efrm(n_objects_per_set = 8, n_sets = 1, n_panels = 3,
                         n_judges_per_panel = 20, reps_within = 120,
                         panel_units = phi_true, seed = 1)
  fit <- befit(d, "object_a", "object_b", winner = "winner", judge = "judge",
                  panels = "panel", object_sets = attr(d, "truth")$object_sets)
  expect_true(fit$converged)
  pt <- fit$phi_table
  lp_true <- log(phi_true)[match(pt$panel, sprintf("panel%d", 1:3))]
  z <- (log(pt$phi) - lp_true) / pt$se_log_phi
  expect_lt(max(abs(z)), 3)                              # each within 3 SE
  # geometric-mean-one normalisation
  expect_equal(exp(mean(log(pt$phi))), 1, tolerance = 1e-8)
  expect_true(all(is.na(pt$p)))
  expect_true(all(is.na(pt$p_adj)))
  # single set: no set units estimated
  expect_equal(nrow(fit$alpha_table), 1L)
})

test_that("planted set units (alpha) and origins (kappa) are recovered", {
  d <- simulate_btl_efrm(n_objects_per_set = 8, n_sets = 2, n_panels = 2,
                         n_judges_per_panel = 15, reps_within = 120,
                         reps_cross = 50, set_units = c(1, 1.4),
                         set_origins = c(0, 0.8), seed = 1)
  fit <- befit(d, "object_a", "object_b", winner = "winner", judge = "judge",
                  panels = "panel", object_sets = attr(d, "truth")$object_sets)
  expect_true(fit$converged)
  a2 <- fit$alpha_table[fit$alpha_table$set == "set2", ]
  za <- (log(a2$alpha) - log(1.4)) / a2$se_log_alpha
  expect_lt(abs(za), 3)                                 # alpha_2 within 3 SE
  k2 <- fit$kappa_table[fit$kappa_table$set == "set2", ]
  zk <- (k2$kappa - 0.8) / k2$se_kappa
  expect_lt(abs(zk), 3)                                 # kappa_2 within 3 SE
  # reference set carries alpha = 1, kappa = 0 with no standard error
  expect_equal(fit$alpha_table$alpha[fit$alpha_table$set == "set1"], 1)
  expect_true(is.na(fit$alpha_table$se_log_alpha[fit$alpha_table$set == "set1"]))
  # the common-scale values recover the truth across ALL objects
  tr <- attr(d, "truth"); vhat <- setNames(fit$objects$v, fit$objects$object)
  expect_gt(cor(vhat, tr$v[names(vhat)]), 0.99)
  # the frames model fits the differing-unit data better than one unit
  expect_gt(fit$equal_unit$difference, 0)
  expect_setequal(fit$unit_omnibus$term,
                  c("panel units (phi)", "set units (alpha)",
                    "set origins (kappa)"))
  expect_true(all(is.na(fit$unit_omnibus$p)))
  expect_true(all(is.na(fit$unit_omnibus$p_adj)))
  expect_true(all(c("p_adj", "significant") %in%
                    names(fit$alpha_table)))
})

test_that("conditional equal-unit fits remain numerically stable under the null", {
  base <- 300L; diffs <- numeric(12)
  for (k in seq_len(12)) {
    d <- simulate_btl_efrm(n_objects_per_set = 6, n_sets = 2, n_panels = 2,
                           n_judges_per_panel = 10, reps_within = 20,
                           reps_cross = 20, seed = base + k)
    fit <- befit(d, "object_a", "object_b", winner = "winner",
                    judge = "judge", panels = "panel",
                    object_sets = attr(d, "truth")$object_sets)
    expect_true(all(is.na(fit$phi_table$p)))
    expect_true(all(is.na(fit$alpha_table$p)))
    diffs[k] <- fit$equal_unit$difference
  }
  expect_lt(max(abs(diffs)), 25)                        # equal-unit gap small
})

test_that("guards fire with informative errors", {
  os <- list(set1 = sprintf("S1O%02d", 1:6), set2 = sprintf("S2O%02d", 1:6))
  d <- simulate_btl_efrm(n_objects_per_set = 6, n_sets = 2, n_panels = 2,
                         reps_within = 15, reps_cross = 15, seed = 2)

  # graded response is not supported in this first implementation
  d$grade <- 1L
  expect_error(
    befit(d, "object_a", "object_b", winner = "winner", judge = "judge",
             panels = "panel", object_sets = os, response = "grade"),
    "dichotomous")

  # an object present in the data but not assigned to any set
  os_missing <- list(set1 = sprintf("S1O%02d", 1:6),
                     set2 = sprintf("S2O%02d", 1:5))    # drops S2O06
  expect_error(
    befit(d, "object_a", "object_b", winner = "winner", judge = "judge",
             panels = "panel", object_sets = os_missing),
    "belong to exactly one set")

  # an object assigned to two sets
  os_dup <- list(set1 = c(sprintf("S1O%02d", 1:6), "S2O01"),
                 set2 = sprintf("S2O%02d", 1:6))
  expect_error(
    befit(d, "object_a", "object_b", winner = "winner", judge = "judge",
             panels = "panel", object_sets = os_dup),
    "more than one set")

  # External mappings use the same canonical labels as the comparison data.
  judges <- unique(d$judge)
  pmap_pad <- stats::setNames(
    d$panel[match(judges, d$judge)], paste0(" ", judges, " "))
  os_pad <- lapply(os, function(x) paste0(" ", x, " "))
  names(os_pad) <- paste0(" ", names(os_pad), " ")
  f_pad <- befit(d, "object_a", "object_b", winner = "winner",
                 judge = "judge", panels = pmap_pad, object_sets = os_pad)
  expect_setequal(f_pad$alpha_table$set, c("set1", "set2"))
  bad_set_names <- os
  names(bad_set_names) <- c("set", " set ")
  expect_error(
    befit(d, "object_a", "object_b", winner = "winner", judge = "judge",
          panels = "panel", object_sets = bad_set_names),
    "after trimming")

  shaped_panels <- matrix(d$panel[match(judges, d$judge)], nrow = 1L)
  names(shaped_panels) <- judges
  expect_error(
    befit(d, "object_a", "object_b", winner = "winner", judge = "judge",
          panels = shaped_panels, object_sets = os),
    "plain named")
  shaped_sets <- os
  shaped_sets[[1L]] <- matrix(shaped_sets[[1L]], nrow = 1L)
  expect_error(
    befit(d, "object_a", "object_b", winner = "winner", judge = "judge",
          panels = "panel", object_sets = shaped_sets),
    "named list")

  # insufficient cross-set links: no set pair reaches min_link
  expect_error(
    befit(d, "object_a", "object_b", winner = "winner", judge = "judge",
             panels = "panel", object_sets = os, min_link = 100000),
    "not reachable from the reference set")

  # disconnected panel-by-set graph: each set is judged by its own panel only,
  # so no set contains comparisons from more than one panel
  set.seed(9)
  o1 <- LETTERS[1:3]; o2 <- LETTERS[4:6]
  mk <- function(objs, judges) {
    pr <- t(utils::combn(objs, 2))
    dd <- data.frame(object_a = rep(pr[, 1], 12), object_b = rep(pr[, 2], 12),
                     stringsAsFactors = FALSE)
    dd$judge <- sample(judges, nrow(dd), replace = TRUE)
    dd$winner <- ifelse(runif(nrow(dd)) < 0.5, dd$object_a, dd$object_b)
    dd
  }
  dw <- rbind(mk(o1, c("j1", "j2")), mk(o2, c("j3", "j4")))
  pmap <- c(j1 = "P1", j2 = "P1", j3 = "P2", j4 = "P2")   # set1->P1, set2->P2
  expect_error(
    befit(dw, "object_a", "object_b", winner = "winner", judge = "judge",
             panels = pmap, object_sets = list(set1 = o1, set2 = o2)),
    "panel")

  # single set: the print states the panel-units model
  ds <- simulate_btl_efrm(n_objects_per_set = 6, n_sets = 1, n_panels = 2,
                          reps_within = 25, seed = 3)
  fs <- befit(ds, "object_a", "object_b", winner = "winner", judge = "judge",
                 panels = "panel", object_sets = attr(ds, "truth")$object_sets)
  expect_output(print(fs), "panel units only")
})

test_that("plot_btl_units draws without error", {
  d <- simulate_btl_efrm(n_objects_per_set = 6, n_sets = 2, n_panels = 2,
                         reps_within = 25, reps_cross = 25,
                         set_units = c(1, 1.3), seed = 4)
  fit <- befit(d, "object_a", "object_b", winner = "winner", judge = "judge",
                  panels = "panel", object_sets = attr(d, "truth")$object_sets)
  pdf(NULL); on.exit(dev.off())
  expect_silent(plot_btl_units(fit))

  failed <- fit
  failed$converged <- FALSE
  expect_error(plot_btl_units(failed), "did not converge")
})

test_that("plot_btl_units uses each unit's inferential reference", {
  d <- simulate_btl_efrm(n_objects_per_set = 6, n_sets = 2, n_panels = 2,
                         reps_within = 25, reps_cross = 25, seed = 404)
  fit <- befit(d, "object_a", "object_b", winner = "winner", judge = "judge",
               panels = "panel", object_sets = attr(d, "truth")$object_sets)
  fit$phi_table$phi <- 1
  fit$phi_table$se_log_phi <- 1
  fit$phi_table$df <- c(4, Inf)
  fit$alpha_table$alpha <- 1
  fit$alpha_table$se_log_alpha <- 1
  fit$alpha_table$df <- c(NA, 9)
  widths <- NULL
  testthat::local_mocked_bindings(
    segments = function(x0, y0, x1, y1, ...)
      widths <<- (x1 - x0) / 2,
    .package = "rasch")
  grDevices::pdf(NULL); on.exit(grDevices::dev.off(), add = TRUE)
  plot_btl_units(fit)
  expect_equal(sort(widths),
               sort(stats::qt(0.975, c(4, Inf, 9))), tolerance = 1e-12)

  # Saved objects from before the df field used the normal reference.
  fit$phi_table$df <- NULL
  fit$alpha_table$df <- NULL
  widths <- NULL
  plot_btl_units(fit)
  expect_length(widths, nrow(fit$phi_table) + nrow(fit$alpha_table))
  expect_equal(widths, rep(stats::qnorm(0.975), length(widths)),
               tolerance = 1e-12)
})

test_that("frame estimates propagate through paired-comparison diagnostics", {
  d <- simulate_btl_efrm(n_objects_per_set = 6, n_sets = 2, n_panels = 2,
                         n_judges_per_panel = 10, reps_within = 30,
                         reps_cross = 30, panel_units = c(0.8, 1.25),
                         set_units = c(1, 1.4), set_origins = c(0, 0.6),
                         seed = 44)
  fit <- befit(d, "object_a", "object_b", "winner", "judge", "panel",
               attr(d, "truth")$object_sets)

  expect_s3_class(fit, "rasch_btl")
  expect_equal(fit$objects$location, fit$objects$v)
  expect_equal(fit$objects$se, fit$objects$se_v)
  expect_true(all(c("expected", "frame_slope", "information") %in%
                    names(fit$comparisons)))

  cmp <- fit$comparisons
  loc <- setNames(fit$objects$location, fit$objects$object)
  expected <- plogis(cmp$frame_slope *
                       (loc[cmp$object_a] - loc[cmp$object_b]))
  expect_equal(cmp$expected, unname(expected), tolerance = 1e-8)

  info <- btl_information(fit)
  expect_equal(info$total, sum(cmp$information), tolerance = 1e-10)
  expect_equal(sum(info$objects$information), 2 * info$total,
               tolerance = 1e-10)
  expect_true(is.finite(fit$total_chisq))
  expect_true(is.finite(fit$total_df))
  expect_true(is.na(fit$total_p))
  expect_match(paste(fit$notes, collapse = " "),
               "row-based probability is withheld")
  expect_no_error(judge_surprise(fit, fit$judges$judge[1]))
  expect_no_error(btl_dimensionality(fit, reps = 20))
  expect_error(btl_next_pairs(fit), "panel and object set")

  tab <- fit_summary_table(fit)
  expect_equal(tab$value[tab$statistic == "Object sets"], "2")
  expect_identical(tab$value[tab$statistic == "Pairwise fit probability"],
                   "unavailable (judge clustering)")
  pdf(NULL); on.exit(dev.off())
  expect_silent(plot_btl_icc(fit, fit$objects$object[1]))
  expect_silent(plot_btl_targeting(fit))
})

test_that("conditional frame fits withhold equating drift inference", {
  d1 <- simulate_btl_efrm(n_objects_per_set = 6, n_sets = 2, n_panels = 2,
                          reps_within = 25, reps_cross = 25, seed = 71)
  d2 <- simulate_btl_efrm(n_objects_per_set = 6, n_sets = 2, n_panels = 2,
                          reps_within = 25, reps_cross = 25, seed = 72)
  f1 <- befit(d1, "object_a", "object_b", "winner", "judge", "panel",
              attr(d1, "truth")$object_sets)
  f2 <- befit(d2, "object_a", "object_b", "winner", "judge", "panel",
              attr(d2, "truth")$object_sets)
  eq <- btl_equate(f1, f2, independent = TRUE)
  expect_false(eq$inferential)
  expect_true(any(grepl("joint object-location covariance", eq$notes,
                        fixed = TRUE)))
})

test_that("bootstrap SEs propagate linking uncertainty (estimates unchanged)", {
  skip_on_cran()
  d <- simulate_btl_efrm(n_objects_per_set = 6, n_sets = 2, n_panels = 2,
                         set_units = c(1, 1.4), set_origins = c(0, 0.8),
                         reps_within = 60, reps_cross = 50, seed = 1)
  os <- attr(d, "truth")$object_sets
  set.seed(9)
  fb <- btl_efrm(d, "object_a", "object_b", "winner", "judge", "panel", os,
                 boot_reps = 40)
  fc <- btl_efrm(d, "object_a", "object_b", "winner", "judge", "panel", os,
                 se_method = "conditional")
  expect_equal(fb$alpha_table$alpha, fc$alpha_table$alpha)   # same estimator
  expect_equal(fb$objects$v, fc$objects$v)
  # The bootstrap refits stage one rather than conditioning on it. Its
  # realised SE need not exceed the conditional SE in every finite set of
  # resamples, so test the distinct, finite result rather than its direction.
  expect_true(is.finite(fb$alpha_table$se_log_alpha[2]))
  expect_gt(fb$alpha_table$se_log_alpha[2], 0)
  expect_false(isTRUE(all.equal(fb$alpha_table$se_log_alpha[2],
                                fc$alpha_table$se_log_alpha[2],
                                tolerance = 1e-8)))
  expect_identical(fb$boot_reps_requested, 40L)
  expect_identical(fb$boot_reps_used + fb$boot_reps_failed,
                   fb$boot_reps_requested)
  expect_true(all(is.finite(fb$unit_omnibus$df2)))
  expect_equal(fb$unit_omnibus$p_adj,
               p.adjust(fb$unit_omnibus$p, "holm"))
  follow_p <- c(fb$phi_table$p, fb$alpha_table$p, fb$kappa_table$p)
  follow_adj <- c(fb$phi_table$p_adj, fb$alpha_table$p_adj,
                  fb$kappa_table$p_adj)
  ok <- is.finite(follow_p)
  expect_equal(follow_adj[ok], p.adjust(follow_p[ok], "holm"))

  # The model-based bootstrap draws comparison outcomes independently. Its
  # reference is therefore normal/chi-square, not the finite-judge reference
  # used after resampling judges.
  set.seed(10)
  fp <- btl_efrm(d, "object_a", "object_b", "winner", "judge", "panel", os,
                 se_method = "bootstrap", boot_reps = 30)
  expect_identical(fp$boot_reps_requested, 30L)
  expect_identical(fp$boot_reps_used + fp$boot_reps_failed,
                   fp$boot_reps_requested)
  expect_true(all(is.infinite(fp$unit_omnibus$df2)))
  expect_true(all(is.infinite(fp$alpha_table$df)))
})

test_that("BTL-EFRM omnibus families are not truncated by unavailable covariance", {
  V <- diag(c(0.04, 0.09))
  full <- .btlef_wald_unit(c(0.2, -0.3), V, "set units")
  expect_equal(full$df, 2L)
  expect_true(is.finite(full$p))

  V[2, 2] <- NA_real_
  unavailable <- .btlef_wald_unit(c(0.2, -0.3), V, "set units")
  expect_identical(unavailable$term, "set units")
  expect_true(is.na(unavailable$df))
  expect_true(is.na(unavailable$wald))
  expect_true(is.na(unavailable$p))

  singular <- .btlef_wald_unit(c(0.2, 0.2),
                               matrix(c(1, 1, 1, 1), 2), "set units")
  expect_equal(singular$df, 1L)
  expect_true(is.finite(singular$p))
  outside <- .btlef_wald_unit(c(0.2, -0.3),
                              matrix(c(1, 1, 1, 1), 2), "set units")
  expect_true(is.na(outside$df))
  expect_true(is.na(outside$p))
})

test_that("judge-bootstrap unit tests respect panel-specific judge support", {
  skip_on_cran()
  d <- simulate_btl_efrm(n_objects_per_set = 6, n_sets = 1, n_panels = 2,
                         n_judges_per_panel = 10, reps_within = 35, seed = 81)
  by_judge <- unique(d[c("judge", "panel")])
  move <- by_judge$judge[by_judge$panel == "panel2"][1:8]
  d$panel[d$judge %in% move] <- "panel1"
  fit <- btl_efrm(d, "object_a", "object_b", "winner", "judge", "panel",
                  attr(d, "truth")$object_sets,
                  se_method = "judge_bootstrap", boot_reps = 30)
  expect_lt(min(fit$unit_support$panel$n_judges), 6)
  expect_true(all(is.na(fit$phi_table$p)))
  expect_true(all(is.na(fit$unit_omnibus$p)))
  expect_true(all(is.finite(fit$phi_table$phi)))
  expect_false(any(grepl("set-unit and set-origin inference is withheld",
                         fit$notes, fixed = TRUE)))
})

test_that("set-unit df follow the support path to the reference set", {
  skip_on_cran()
  d <- simulate_btl_efrm(
    n_objects_per_set = 4, n_sets = 3, n_judges_per_panel = 30,
    n_panels = 1, reps_within = 40, reps_cross = 40, seed = 181
  )
  os <- attr(d, "truth")$object_sets
  set_of <- setNames(rep(names(os), lengths(os)), unlist(os))
  sa <- unname(set_of[d$object_a])
  sb <- unname(set_of[d$object_b])
  lo <- pmin(sa, sb)
  hi <- pmax(sa, sb)
  same <- sa == sb
  keep <- same |
    (lo == "set1" & hi == "set2" &
       d$judge %in% sprintf("J%03d", 1:9)) |
    (lo == "set2" & hi == "set3" &
       d$judge %in% sprintf("J%03d", 1:20))

  fit <- btl_efrm(
    d[keep, ], "object_a", "object_b", "winner", "judge", "panel", os,
    se_method = "judge_bootstrap", boot_reps = 30, workers = 1, seed = 20
  )
  edge <- fit$unit_support$edge
  weak <- edge$effective_judges[
    edge$set_a == "set1" & edge$set_b == "set2"
  ]
  support <- setNames(
    fit$unit_support$set$effective_judges,
    fit$unit_support$set$set
  )
  dfs <- setNames(fit$alpha_table$df, fit$alpha_table$set)

  expect_equal(unname(support[c("set2", "set3")]), rep(weak, 2))
  expect_equal(unname(dfs[c("set2", "set3")]), rep(floor(weak) - 1, 2))
  # The sole panel unit is fixed at one and is not a hypothesis. It must not
  # enlarge the Holm family for the four free set-unit/origin questions.
  free_set <- fit$alpha_table$set != fit$reference_set
  follow_p <- c(fit$alpha_table$p[free_set],
                fit$kappa_table$p[free_set])
  follow_adj <- c(fit$alpha_table$p_adj[free_set],
                  fit$kappa_table$p_adj[free_set])
  expect_true(all(is.finite(follow_p)))
  expect_equal(follow_adj, p.adjust(follow_p, "holm"))
  expect_true(is.na(fit$phi_table$p_adj))
  set3_objects <- fit$objects$object[fit$objects$set == "set3"]
  expect_equal(.btl_equate_cov_df(fit, set3_objects), floor(weak) - 1)
})

test_that("set pairs below min_link do not enter the linking fit", {
  d <- simulate_btl_efrm(
    n_objects_per_set = 3, n_sets = 3, n_judges_per_panel = 8,
    n_panels = 2, reps_within = 20, reps_cross = 3, seed = 44
  )
  os <- attr(d, "truth")$object_sets
  set_of <- setNames(rep(names(os), lengths(os)), unlist(os))
  sa <- unname(set_of[d$object_a])
  sb <- unname(set_of[d$object_b])
  lo <- pmin(sa, sb)
  hi <- pmax(sa, sb)
  weak_edge <- lo == "set1" & hi == "set3"
  keep <- !weak_edge
  keep[which(weak_edge)[1L]] <- TRUE
  sparse <- d[keep, ]

  with_weak <- btl_efrm(
    sparse, "object_a", "object_b", "winner", "judge", "panel", os,
    min_link = 20, se_method = "conditional", boot_reps = 0
  )
  without_weak <- btl_efrm(
    sparse[!(unname(set_of[sparse$object_a]) %in% c("set1", "set3") &
               unname(set_of[sparse$object_b]) %in% c("set1", "set3") &
               unname(set_of[sparse$object_a]) !=
                 unname(set_of[sparse$object_b])), ],
    "object_a", "object_b", "winner", "judge", "panel", os,
    min_link = 20, se_method = "conditional", boot_reps = 0
  )

  expect_equal(with_weak$alpha_table$alpha, without_weak$alpha_table$alpha)
  expect_equal(with_weak$kappa_table$kappa, without_weak$kappa_table$kappa)
  expect_equal(with_weak$objects$location, without_weak$objects$location)
  expect_false(with_weak$n_cross$used[
    with_weak$n_cross$set_a == "set1" & with_weak$n_cross$set_b == "set3"
  ])
  expect_true(any(grepl("omitted from estimation", with_weak$notes,
                        fixed = TRUE)))
})

test_that("bootstrap SEs are calibrated on the chain-linked design", {
  skip_on_cran()   # ~6 x 41 pipeline fits; the full battery: tools/calibration.R
  la <- se <- numeric(6)
  for (s in 1:6) {
    d <- simulate_btl_efrm(n_objects_per_set = 7, n_sets = 3,
                           n_judges_per_panel = 8, n_panels = 2,
                           reps_within = 60, reps_cross = 60,
                           panel_units = c(0.8, 1.25),
                           set_units = c(1, 1.3, 0.75),
                           set_origins = c(0, 0.5, -0.4), seed = 100 + s)
    so <- attr(d, "truth")$set_of
    sa <- so[d$object_a]; sb <- so[d$object_b]
    d2 <- d[!((sa == 1 & sb == 3) | (sa == 3 & sb == 1)), ]   # sever A-C
    set.seed(500 + s)
    f <- btl_efrm(d2, "object_a", "object_b", "winner", "judge", "panel",
                  attr(d, "truth")$object_sets, boot_reps = 40)
    la[s] <- log(f$alpha_table$alpha[2]); se[s] <- f$alpha_table$se_log_alpha[2]
  }
  covered <- sum(abs(la - log(1.3)) <= 1.96 * se)
  expect_gte(covered, 4L)   # was 4/12 with conditional SEs; bootstrap restores
})

test_that("a set with no stable panel-ratio information is screened, not fatal", {
  # The 'weak' set has one within pair whose preference DIRECTION flips
  # between the panels: its panel ratio logit(p_yes)/logit(p_no) is negative,
  # which no positive-unit parameterisation can represent, and before the
  # screen the stage-1 solver diverged and poisoned the phi reconciliation
  # (found live on GermanParties2009, where the CDU/CSU:FDP pair is
  # near-even). The strong set identifies phi; the weak set is refit at the
  # reconciled units and noted.
  strong <- simulate_btl_efrm(n_objects_per_set = 5, n_sets = 1, n_panels = 2,
                              n_judges_per_panel = 8, reps_within = 40,
                              panel_units = c(0.8, 1.25), seed = 42)
  sobj <- sort(unique(c(strong$object_a, strong$object_b)))   # "S1O01".."S1O05"
  jd_of <- function(pnl, i) sprintf("%s_J%d", pnl, (i %% 8) + 1)
  flip <- do.call(rbind, lapply(seq_len(160), function(i) {
    pnl <- if (i <= 80) "panel1" else "panel2"
    win <- if (i <= 80) (i %% 10 < 8) else (i %% 10 >= 8)   # 80/20 vs 20/80
    data.frame(object_a = "w1", object_b = "w2",
               winner = if (win) "w1" else "w2",
               judge = jd_of(pnl, i), panel = pnl)
  }))
  cross <- do.call(rbind, lapply(seq_len(200), function(i) {
    pnl <- if (i <= 100) "panel1" else "panel2"
    data.frame(object_a = sobj[(i %% 5) + 1],
               object_b = paste0("w", (i %% 2) + 1),
               # Preserve positive cross-set support for the reconciled
               # weak-set ordering (w1 below w2). The former nearly equal
               # cross rates hit alpha = 0 and silently used unit one.
               winner = if (i %% 10 < if (i %% 2 == 0) 2 else 8)
                          paste0("w", (i %% 2) + 1)
                        else sobj[(i %% 5) + 1],
               judge = jd_of(pnl, i), panel = pnl)
  }))
  d <- rbind(strong[, names(flip)], flip, cross)
  os <- list(strong = sobj, weak = c("w1", "w2"))
  expect_no_error(fit <- befit(d, "object_a", "object_b", "winner", "judge",
                               "panel", os))
  expect_true(fit$converged)
  expect_true(any(grepl("weak", fit$notes) & grepl("panel-ratio", fit$notes)))
  # phi comes from the strong set alone and is finite and sane
  expect_true(all(is.finite(fit$phi_table$phi)))
  expect_true(all(fit$phi_table$phi > 0.2 & fit$phi_table$phi < 5))
  # the weak set's objects are still located on the common scale
  expect_true(all(is.finite(fit$objects$v[fit$objects$set == "weak"])))
})

test_that("estimates and convergence are invariant to duplicating the data", {
  skip_on_cran()
  # an absolute gradient threshold is scale-dependent: on k-fold duplicated
  # data a converged fit was flagged unconverged, which (with the stability
  # screen) silently rerouted a set's estimation and CHANGED the estimates;
  # the per-comparison criterion is invariant
  d <- simulate_btl_efrm(n_objects_per_set = 6, n_sets = 2, n_panels = 2,
                         n_judges_per_panel = 6, reps_within = 25,
                         reps_cross = 25, panel_units = c(0.8, 1.25),
                         set_units = c(1, 1.3), seed = 7)
  os <- attr(d, "truth")$object_sets
  f1 <- befit(d, "object_a", "object_b", "winner", "judge", "panel", os)
  d50 <- d[rep(seq_len(nrow(d)), 50), ]
  f2 <- befit(d50, "object_a", "object_b", "winner", "judge", "panel", os)
  expect_true(f1$converged && f2$converged)
  expect_equal(f2$phi_table$phi, f1$phi_table$phi, tolerance = 1e-6)
  expect_equal(f2$alpha_table$alpha, f1$alpha_table$alpha, tolerance = 1e-6)
  expect_equal(f2$objects$v, f1$objects$v, tolerance = 1e-6)
})

test_that("btl_efrm refuses (quasi-)complete cross-set separation", {
  set.seed(3)
  K <- 5
  o1 <- sprintf("S1O%02d", seq_len(K)); o2 <- sprintf("S2O%02d", seq_len(K))
  judges <- sprintf("J%03d", 1:8)
  gen_within <- function(objs, beta, reps) {
    pr <- t(combn(objs, 2)); aa <- bb <- character(0)
    for (r in seq_len(nrow(pr))) {
      aa <- c(aa, rep(pr[r, 1], reps)); bb <- c(bb, rep(pr[r, 2], reps)) }
    p <- plogis(beta[aa] - beta[bb])
    data.frame(object_a = aa, object_b = bb,
               winner = ifelse(runif(length(aa)) < p, aa, bb),
               stringsAsFactors = FALSE)
  }
  b1 <- setNames(as.numeric(scale(seq_len(K))), o1)
  b2 <- setNames(as.numeric(scale(seq_len(K))), o2)
  w1 <- gen_within(o1, b1, 20); w2 <- gen_within(o2, b2, 20)
  # every set2 object beats every set1 object, always: perfect separation
  grid <- expand.grid(oa = o1, ob = o2, stringsAsFactors = FALSE)
  c12 <- data.frame(object_a = rep(grid$oa, 25), object_b = rep(grid$ob, 25),
                    winner = rep(grid$ob, 25), stringsAsFactors = FALSE)
  d <- rbind(w1, w2, c12)
  d$judge <- sample(judges, nrow(d), TRUE); d$panel <- "panel1"
  os <- list(set1 = o1, set2 = o2)
  # both the conditional and the default bootstrap path must refuse, not
  # report a boundary alpha/kappa with a fabricated SE
  expect_error(suppressWarnings(
    btl_efrm(d, "object_a", "object_b", "winner", judge = "judge",
             panels = "panel", object_sets = os,
             se_method = "conditional", min_link = 20)),
    "separated")
  expect_error(suppressWarnings(
    btl_efrm(d, "object_a", "object_b", "winner", judge = "judge",
             panels = "panel", object_sets = os,
             se_method = "bootstrap", boot_reps = 40, min_link = 20)),
    "separated")
})

test_that("BTL-EFRM judge bootstrap reports progress and restores the RNG", {
  skip_on_cran()
  d <- simulate_btl_efrm(n_objects_per_set = 6, n_sets = 2, n_panels = 2,
                         n_judges_per_panel = 8, reps_within = 25,
                         reps_cross = 25, seed = 701)
  os <- attr(d, "truth")$object_sets
  seen <- list()
  set.seed(702)
  before <- .Random.seed
  f <- btl_efrm(d, "object_a", "object_b", "winner", "judge", "panel", os,
                boot_reps = 30, workers = 1, seed = 703,
                progress = function(stage, current, total)
                  seen[[length(seen) + 1L]] <<- c(stage, current, total))
  expect_identical(.Random.seed, before)
  expect_identical(f$workers, 1L)
  expect_identical(f$seed, 703L)
  stages <- vapply(seen, `[`, "", 1L)
  expect_true(all(c("two-stage fit", "judge bootstrap", "finalising") %in%
                    stages))
  js <- seen[stages == "judge bootstrap"]
  expect_identical(as.integer(tail(js, 1L)[[1L]][2:3]), c(30L, 30L))

  completed <- 0L
  expect_condition(
    btl_efrm(d, "object_a", "object_b", "winner", "judge", "panel", os,
      boot_reps = 30, workers = 1, seed = 704,
      progress = function(stage, current, total)
        if (stage == "judge bootstrap") completed <<- as.integer(current),
      cancel = function() completed >= 2L),
    class = "rasch_cancelled")
  expect_lte(completed, 2L)
})

test_that("parallel BTL-EFRM judge bootstraps are seed-identical", {
  skip_on_cran()
  skip_if_not(rasch:::.rasch_namespace_is_installed(),
              "parallel integration test needs an installed package namespace")
  expect_true(rasch:::.rasch_namespace_is_installed())
  old_workers <- options(rasch.max_workers = 2L)
  on.exit(options(old_workers), add = TRUE)
  probe <- try(parallel::makePSOCKcluster(2L), silent = TRUE)
  skip_if(inherits(probe, "try-error"), "local socket clusters unavailable")
  parallel::stopCluster(probe)

  d <- simulate_btl_efrm(n_objects_per_set = 6, n_sets = 2, n_panels = 2,
                         n_judges_per_panel = 8, reps_within = 25,
                         reps_cross = 25, seed = 705)
  os <- attr(d, "truth")$object_sets
  serial <- btl_efrm(d, "object_a", "object_b", "winner", "judge", "panel",
                     os, boot_reps = 30, workers = 1, seed = 706)
  para <- btl_efrm(d, "object_a", "object_b", "winner", "judge", "panel",
                   os, boot_reps = 30, workers = 2, seed = 706)
  expect_identical(para$phi_table, serial$phi_table)
  expect_identical(para$alpha_table, serial$alpha_table)
  expect_identical(para$kappa_table, serial$kappa_table)
  expect_identical(para$objects, serial$objects)
  expect_identical(para$unit_omnibus, serial$unit_omnibus)
})

test_that("BTL-EFRM defaults to four available workers", {
  expect_identical(formals(btl_efrm)$workers, 4L)
})

test_that("BTL-EFRM sizes bootstrap rank to its actual covariance blocks", {
  d <- simulate_btl_efrm(n_objects_per_set = 16, n_sets = 2,
                         n_panels = 2, n_judges_per_panel = 4,
                         reps_within = 1, reps_cross = 1, seed = 707)
  expect_error(
    btl_efrm(d, "object_a", "object_b", "winner", "judge", "panel",
             attr(d, "truth")$object_sets, boot_reps = 30, workers = 1),
    "at least 32 replicates")
})
