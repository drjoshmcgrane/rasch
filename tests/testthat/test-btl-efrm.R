# fast path for tests whose subject is not the standard errors: the
# conditional SEs are exact for beta/phi and the estimates are identical
befit <- function(...) btl_efrm(..., se_method = "conditional")

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
})

test_that("planted panel units (phi) are recovered and Wald-flagged", {
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
})

test_that("units are not spuriously flagged under the null", {
  base <- 300L; flags <- 0L; diffs <- numeric(12)
  for (k in seq_len(12)) {
    d <- simulate_btl_efrm(n_objects_per_set = 6, n_sets = 2, n_panels = 2,
                           n_judges_per_panel = 10, reps_within = 20,
                           reps_cross = 20, seed = base + k)
    fit <- befit(d, "object_a", "object_b", winner = "winner",
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
  # the bootstrap carries stage-one noise the conditional errors omit
  expect_gt(fb$alpha_table$se_log_alpha[2], fc$alpha_table$se_log_alpha[2])
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
               winner = if (i %% 3 == 0) paste0("w", (i %% 2) + 1)
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
