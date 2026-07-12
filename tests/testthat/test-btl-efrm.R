# Extended frame of reference for paired comparisons: reduction to btl(),
# recovery of panel and set units, null calibration, and the guards.

test_that("G = 1, S = 1 reduces exactly to btl()", {
  d <- simulate_btl_efrm(n_objects_per_set = 7, n_sets = 1, n_panels = 1,
                         reps_within = 40, seed = 1)
  fit <- btl_efrm(d, "object_a", "object_b", winner = "winner", judge = "judge",
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
})

test_that("planted panel units (phi) are recovered and Wald-flagged", {
  phi_true <- c(0.7, 1.0, 1.43); phi_true <- phi_true / exp(mean(log(phi_true)))
  d <- simulate_btl_efrm(n_objects_per_set = 8, n_sets = 1, n_panels = 3,
                         n_judges_per_panel = 20, reps_within = 120,
                         panel_units = phi_true, seed = 1)
  fit <- btl_efrm(d, "object_a", "object_b", winner = "winner", judge = "judge",
                  panels = "panel", object_sets = attr(d, "truth")$object_sets)
  expect_true(fit$converged)
  pt <- fit$phi_table
  lp_true <- log(phi_true)[match(pt$panel, sprintf("panel%d", 1:3))]
  z <- (log(pt$phi) - lp_true) / pt$se_log_phi
  expect_lt(max(abs(z)), 3)                              # each within 3 SE
  # geometric-mean-one normalisation
  expect_equal(exp(mean(log(pt$phi))), 1, tolerance = 1e-8)
  # the two panels with phi != 1 flag; the phi = 1 panel does not
  flagged <- pt$p < 0.05
  expect_true(flagged[pt$panel == "panel1"])            # 0.70
  expect_true(flagged[pt$panel == "panel3"])            # 1.43
  expect_false(flagged[pt$panel == "panel2"])           # 1.00
  # single set: no set units estimated
  expect_equal(nrow(fit$alpha_table), 1L)
})

test_that("planted set units (alpha) and origins (kappa) are recovered", {
  d <- simulate_btl_efrm(n_objects_per_set = 8, n_sets = 2, n_panels = 2,
                         n_judges_per_panel = 15, reps_within = 120,
                         reps_cross = 50, set_units = c(1, 1.4),
                         set_origins = c(0, 0.8), seed = 1)
  fit <- btl_efrm(d, "object_a", "object_b", winner = "winner", judge = "judge",
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
})

test_that("units are not spuriously flagged under the null", {
  base <- 300L; flags <- 0L; diffs <- numeric(12)
  for (k in seq_len(12)) {
    d <- simulate_btl_efrm(n_objects_per_set = 6, n_sets = 2, n_panels = 2,
                           n_judges_per_panel = 10, reps_within = 20,
                           reps_cross = 20, seed = base + k)
    fit <- btl_efrm(d, "object_a", "object_b", winner = "winner",
                    judge = "judge", panels = "panel",
                    object_sets = attr(d, "truth")$object_sets)
    flags <- flags + any(fit$phi_table$p < 0.05, na.rm = TRUE) +
      any(fit$alpha_table$p < 0.05, na.rm = TRUE)
    diffs[k] <- fit$equal_unit$difference
  }
  expect_lte(flags, 2L)                                 # loose binomial bound
  expect_lt(max(abs(diffs)), 25)                        # equal-unit gap small
})

test_that("guards fire with informative errors", {
  os <- list(set1 = sprintf("S1O%02d", 1:6), set2 = sprintf("S2O%02d", 1:6))
  d <- simulate_btl_efrm(n_objects_per_set = 6, n_sets = 2, n_panels = 2,
                         reps_within = 15, reps_cross = 15, seed = 2)

  # graded response is not supported in this first implementation
  d$grade <- 1L
  expect_error(
    btl_efrm(d, "object_a", "object_b", winner = "winner", judge = "judge",
             panels = "panel", object_sets = os, response = "grade"),
    "dichotomous")

  # an object present in the data but not assigned to any set
  os_missing <- list(set1 = sprintf("S1O%02d", 1:6),
                     set2 = sprintf("S2O%02d", 1:5))    # drops S2O06
  expect_error(
    btl_efrm(d, "object_a", "object_b", winner = "winner", judge = "judge",
             panels = "panel", object_sets = os_missing),
    "belong to exactly one set")

  # an object assigned to two sets
  os_dup <- list(set1 = c(sprintf("S1O%02d", 1:6), "S2O01"),
                 set2 = sprintf("S2O%02d", 1:6))
  expect_error(
    btl_efrm(d, "object_a", "object_b", winner = "winner", judge = "judge",
             panels = "panel", object_sets = os_dup),
    "more than one set")

  # insufficient cross-set links: no set pair reaches min_link
  expect_error(
    btl_efrm(d, "object_a", "object_b", winner = "winner", judge = "judge",
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
    btl_efrm(dw, "object_a", "object_b", winner = "winner", judge = "judge",
             panels = pmap, object_sets = list(set1 = o1, set2 = o2)),
    "panel")

  # single set: the print states the panel-units model
  ds <- simulate_btl_efrm(n_objects_per_set = 6, n_sets = 1, n_panels = 2,
                          reps_within = 25, seed = 3)
  fs <- btl_efrm(ds, "object_a", "object_b", winner = "winner", judge = "judge",
                 panels = "panel", object_sets = attr(ds, "truth")$object_sets)
  expect_output(print(fs), "panel units only")
})

test_that("plot_btl_units draws without error", {
  d <- simulate_btl_efrm(n_objects_per_set = 6, n_sets = 2, n_panels = 2,
                         reps_within = 25, reps_cross = 25,
                         set_units = c(1, 1.3), seed = 4)
  fit <- btl_efrm(d, "object_a", "object_b", winner = "winner", judge = "judge",
                  panels = "panel", object_sets = attr(d, "truth")$object_sets)
  pdf(NULL); on.exit(dev.off())
  expect_silent(plot_btl_units(fit))
})
