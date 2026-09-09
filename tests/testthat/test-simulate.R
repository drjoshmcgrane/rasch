# Data simulation: each planted departure must trip its matching diagnostic
# (the plant -> detect loop), and the truth is attached for recovery checks.

test_that("simulate_rasch plants misfit the Rasch diagnostics detect", {
  # discrimination: a central over-discriminating item overfits (low outfit)
  d <- simulate_rasch(600, 11, discrimination = c(rep(1, 5), 3, rep(1, 5)),
                      seed = 1)
  f <- rasch(d)
  # the over-discriminating item overfits: its outfit is the lowest and
  # clearly below expectation (extreme persons are excluded from item fit,
  # so the mean-square is not deflated by their boundary residuals)
  expect_lt(f$items$outfit_ms[6], 0.75)
  expect_equal(which.min(f$items$outfit_ms), 6L)
  expect_s3_class(d, "rasch_sim")

  # DIF flags the planted item and (essentially) nothing else
  d <- simulate_rasch(800, 10, dif = list(items = "I05", uniform = 1.2),
                      n_groups = 2, seed = 2)
  s <- dif_anova(rasch(d, factors = "group"))$summary
  expect_true(s$uniform_DIF[s$item == "I05"])
  expect_equal(sum(s$uniform_DIF[s$item != "I05"]), 0L)

  # local dependence flags the planted pair on Q3*
  d <- simulate_rasch(1000, 10,
                      dependence = list(pairs = list(c("I03", "I04")),
                                        strength = 2.5), seed = 4)
  h <- residual_correlations(rasch(d), flag = 0.2)$flagged
  expect_true(any((h$item_a == "I03" & h$item_b == "I04") |
                  (h$item_a == "I04" & h$item_b == "I03")))

  # careless responders inflate person outfit
  d <- simulate_rasch(600, 15, careless = 0.12, seed = 6)
  po <- rasch(d)$person$outfit_ms; ci <- attr(d, "truth")$careless_idx
  expect_gt(mean(po[ci], na.rm = TRUE), mean(po[-ci], na.rm = TRUE) + 0.5)

  expect_output(print(d), "careless")
})

test_that("response-style probabilities remain stable at large strengths", {
  extreme <- simulate_rasch(
    40, 5, model = "PCM", n_categories = 3,
    response_style = list(type = "extreme", prop = 0.5, strength = 1000),
    seed = 91)
  sx <- attr(extreme, "truth")$style_idx
  expect_true(all(as.matrix(extreme[sx, paste0("I0", 1:5)]) %in% c(0L, 2L)))

  middle <- simulate_rasch(
    40, 5, model = "PCM", n_categories = 3,
    response_style = list(type = "middle", prop = 0.5, strength = 1000),
    seed = 92)
  sm <- attr(middle, "truth")$style_idx
  expect_true(all(as.matrix(middle[sm, paste0("I0", 1:5)]) == 1L))
})

test_that("a requested threshold disorder is present in the generating truth", {
  # PCM threshold spans vary by item. The disorder must therefore be imposed
  # as an adjacent reversal, not as a subtraction that a wide random gap can
  # sometimes absorb.
  for (seed in 1:30) {
    d <- simulate_rasch(
      40, 4, model = "PCM", n_categories = 3,
      disordered = "I01", seed = seed)
    tr <- attr(d, "truth")
    expect_named(tr$thresholds, sprintf("I%02d", 1:4))
    expect_true(any(diff(tr$thresholds$I01) < 0))
    expect_equal(mean(tr$thresholds$I01), tr$difficulty[["I01"]],
                 tolerance = 1e-12)
  }
})

test_that("simulate_btl plants misfit the paired-comparison diagnostics detect", {
  # erratic judges carry large fit residuals and low consistency
  d <- simulate_btl(8, 12, erratic_judges = 0.17, seed = 1)
  bt <- btl(d, "object_a", "object_b", winner = "winner", judge = "judge")
  er <- attr(d, "truth")$erratic
  expect_gt(mean(bt$judges$fit_resid[bt$judges$judge %in% er]),
            mean(bt$judges$fit_resid[!bt$judges$judge %in% er]) + 1)

  # graded comparisons recover the object locations
  d <- simulate_btl(8, 12, model = "graded", n_categories = 4, seed = 4)
  bt <- btl(d, "object_a", "object_b", response = "response", judge = "judge")
  expect_gt(cor(bt$objects$location,
                attr(d, "truth")$location[bt$objects$object]), 0.95)

  # a planted carry-over dependence is recovered
  d <- simulate_btl(6, 10, reps_per_pair = 40,
                    dependence = list(exposure = 0, carry_over = 1.2), seed = 3)
  bt <- btl(d, "object_a", "object_b", winner = "winner", judge = "judge",
            order = "order")
  co <- bt$dependence$estimate[bt$dependence$effect == "carry_over"]
  expect_gt(co, 0.5)
})

test_that("simulate_btl accepts explanatory object locations", {
  loc <- c(O1 = -1.2, O2 = -0.1, O3 = 0.2, O4 = 1.8)
  d <- simulate_btl(4, 8, 20, object_locations = loc, seed = 81)
  expect_equal(attr(d, "truth")$location, loc - mean(loc))
  expect_error(simulate_btl(4, object_locations = c(1, 2, 3)),
               "length 4")
  expect_error(simulate_btl(4, object_locations = c(A = 1, B = 2, C = 3, D = 4)),
               "every generated object")
  expect_error(simulate_btl(4, object_locations = rep(1, 4),
                            second_attribute = list(rho = 0.2)),
               "positive spread")
})

test_that("explanatory simulation departures cannot be absorbed by the design", {
  p <- data.frame(exposure = seq(-1, 1, length.out = 5),
                  type = rep(c("A", "B"), length.out = 5))
  B <- stats::model.matrix(~ exposure * type, p)
  d <- .sim_explanatory_departure(B, 1.25)
  expect_equal(unname(drop(crossprod(B, d$values))), rep(0, ncol(B)),
               tolerance = 1e-10)
  expect_equal(d$values[d$index], 1.25, tolerance = 1e-12)

  saturated <- stats::model.matrix(
    ~ exposure * type,
    data.frame(exposure = seq(-1, 1, length.out = 4),
               type = rep(c("A", "B"), length.out = 4)))
  expect_error(.sim_explanatory_departure(saturated, 1), "saturated")
})

test_that("simulate_mfrm plants rater severity, misfit, and interaction", {
  d <- simulate_mfrm(120, 5, 6, rater_severity_sd = 0.7, erratic_raters = 0.17,
                     interaction = list(rater = "R3", item = "I2", bias = 1.8),
                     seed = 1)
  mf <- rasch_mfrm(d, person = "person", item = "item", score = "score",
                   facets = "rater", interaction = "rater")
  tr <- attr(d, "truth"); fe <- mf$facet_effects$rater
  rec <- fe$severity[match(names(tr$severity), fe$level)]
  # severities are recovered for the well-behaved raters (the erratic one's
  # true severity is meaningless once it rates at random)
  keep <- !(names(tr$severity) %in% tr$erratic)
  expect_gt(abs(cor(rec[keep], tr$severity[keep])), 0.9)
  er <- tr$erratic
  expect_gt(mean(fe$fit_resid[fe$level %in% er]),
            mean(fe$fit_resid[!fe$level %in% er]) + 1)      # erratic rater misfits
  ie <- mf$interaction_effects                             # interaction at R3xI2
  expect_equal(ie[which.max(abs(ie$gamma)), c("item", "level")],
               data.frame(item = "I2", level = "R3"), ignore_attr = TRUE)
})

test_that("simulation SD arguments are realised sample standard deviations", {
  rr <- simulate_rasch(20, 5, theta_mean = 0.4, theta_sd = 1.6, seed = 119)
  expect_equal(mean(attr(rr, "truth")$theta), 0.4, tolerance = 1e-12)
  expect_equal(stats::sd(attr(rr, "truth")$theta), 1.6,
               tolerance = 1e-12)

  d <- simulate_mfrm(20, 5, 3, theta_sd = 1.4, item_sd = 1.7,
                      rater_severity_sd = 0.8, seed = 120)
  expect_equal(stats::sd(attr(d, "truth")$theta), 1.4,
               tolerance = 1e-12)
  expect_equal(stats::sd(attr(d, "truth")$difficulty), 1.7,
               tolerance = 1e-12)
  expect_equal(stats::sd(attr(d, "truth")$severity), 0.8,
               tolerance = 1e-12)
  equal <- simulate_mfrm(20, 5, 3, item_sd = 0, seed = 120)
  expect_identical(unname(attr(equal, "truth")$difficulty), rep(0, 5))

  bt <- simulate_btl(5, 4, 2, object_sd = 1.3, seed = 121)
  expect_equal(stats::sd(attr(bt, "truth")$location), 1.3,
               tolerance = 1e-12)
  ef <- simulate_efrm(10, 3, theta_sd = 1.1, seed = 122)
  expect_equal(stats::sd(attr(ef, "truth")$theta), 1.1,
               tolerance = 1e-12)
  bf <- simulate_btl_efrm(4, 2, 2, 2, 2, 2,
                           object_sd = 0.9, seed = 123)
  tr <- attr(bf, "truth")
  for (s in unique(tr$set_of))
    expect_equal(stats::sd(tr$beta[tr$set_of == s]), 0.9,
                 tolerance = 1e-12)

  expect_error(simulate_rasch(
    20, 3, theta_mean = 1e300, theta_sd = 1, seed = 124
  ), "cannot be represented")
  expect_error(simulate_mfrm(20, 3, 3, item_sd = 1e-320, seed = 125),
               "item-difficulty standard deviation")
  expect_error(simulate_btl(4, 4, 2, object_sd = 1e-320, seed = 126),
               "object-location standard deviation")
  expect_error(simulate_btl_efrm(
    3, 2, 2, 2, 2, 2, object_sd = 1e-320, seed = 127
  ), "within-set object-location standard deviation")
})

test_that("simulate_efrm plants a frame-unit ratio rasch_efrm recovers", {
  d <- simulate_efrm(300, 8, set_unit_ratio = 1.35, seed = 2)
  tr <- attr(d, "truth")
  ef <- rasch_efrm(d, item_sets = tr$item_sets, groups = "group")
  ratio <- max(ef$alpha_table$alpha) / min(ef$alpha_table$alpha)
  expect_gt(ratio, 1.2); expect_lt(ratio, 1.55)          # ~1.35 recovered
  expect_output(print(d), "set-unit ratio")
  expect_output(print(d), "Generating values and departures")
})

test_that("simulation printing separates generating values from departures", {
  ordinary <- simulate_rasch(20, 3, seed = 201)
  expect_output(print(ordinary),
                "Default model-conforming settings \\(no departures planted\\)")

  # Unequal frame units are part of the fitted EFRM, not model misfit.
  framed <- simulate_efrm(20, 3, set_unit_ratio = 1.2, seed = 202)
  shown <- capture.output(print(framed))
  expect_true(any(grepl("Generating values and departures", shown,
                        fixed = TRUE)))
  expect_false(any(grepl("Planted departures", shown, fixed = TRUE)))
})

test_that("simulate_efrm generates partial credit items on request", {
  d <- simulate_efrm(250, 6, n_sets = 2, n_groups = 1, set_unit_ratio = 1.3,
                     n_categories = 4, seed = 5)
  tr <- attr(d, "truth")
  X <- as.matrix(d[, unlist(tr$item_sets)])
  expect_setequal(sort(unique(as.vector(X))), 0:3)
  expect_length(tr$thresholds, ncol(X))
  expect_true(all(vapply(tr$thresholds, length, 1L) == 3L))
  # thresholds centre on the item locations
  expect_equal(unname(vapply(tr$thresholds, mean, 0)),
               unname(tr$difficulty), tolerance = 1e-10)
  # the dichotomous draw stream is untouched by the generalisation
  d2 <- simulate_efrm(50, 4, n_sets = 2, n_groups = 2,
                      set_unit_ratio = 1.3, seed = 1)
  expect_identical(sum(as.matrix(d2[, 2:9])), 402L)
})

test_that("EFRM simulators plant frame-specific and judge departures", {
  clean <- simulate_efrm(600, 5, n_sets = 2, n_groups = 2,
                         set_unit_ratio = 1, seed = 710)
  drift <- simulate_efrm(600, 5, n_sets = 2, n_groups = 2,
    set_unit_ratio = 1,
    item_drift = list(items = "S1I03", group = "g2", shift = 1.2),
    careless = 0.05, missing = 0.03, seed = 710)
  tr <- attr(drift, "truth")
  expect_identical(tr$item_drift$group, "g2")
  expect_identical(tr$item_drift$items, "S1I03")
  expect_length(tr$careless_idx, 60L)
  expect_equal(length(tr$missing_cells), round(0.03 * 1200 * 10))
  # The same seed isolates the shifted item before the later contamination.
  g2 <- clean$group == "g2"
  expect_lt(mean(drift$S1I03[g2], na.rm = TRUE),
            mean(clean$S1I03[g2], na.rm = TRUE) - 0.08)
  expect_true(any(grepl("item drift", tr$planted)))

  b <- simulate_btl_efrm(5, 2, n_judges_per_panel = 5, n_panels = 2,
                         reps_within = 4, reps_cross = 4,
                         erratic_judges = 0.2, seed = 711)
  expect_length(attr(b, "truth")$erratic, 2L)
  expect_true(any(grepl("erratic judge", attr(b, "truth")$planted)))
  expect_setequal(unique(b$judge), sprintf("J%03d", 1:10))
  expect_setequal(unique(b$panel), c("panel1", "panel2"))
  expect_true(all(attr(b, "truth")$erratic %in% b$judge))
  expect_equal(length(unique(b$panel[b$judge %in%
                                    attr(b, "truth")$erratic])), 2L)
  expect_lte(diff(range(table(b$judge))), 1L)
  set_of <- attr(b, "truth")$set_of
  stratum <- paste(pmin(set_of[b$object_a], set_of[b$object_b]),
                   pmax(set_of[b$object_a], set_of[b$object_b]), sep = ":")
  for (ss in unique(stratum))
    expect_lte(diff(range(table(factor(b$panel[stratum == ss],
                                       levels = c("panel1", "panel2"))))), 1L)
})

test_that("BTL-EFRM simulator unit arguments are plain vectors", {
  expect_error(
    simulate_btl_efrm(n_panels = 2, panel_units = matrix(c(1, 1.2), 1)),
    "panel_units")
  expect_error(
    simulate_btl_efrm(n_sets = 2, set_units = matrix(c(1, 1.2), 1)),
    "set_units")
  expect_error(
    simulate_btl_efrm(n_sets = 2, set_origins = matrix(c(0, 0.2), 1)),
    "set_origins")
  expect_error(
    simulate_btl_efrm(n_sets = 2, set_units = c(1e-300, 1e300)),
    "numerically representable")
  expect_error(
    simulate_btl_efrm(n_sets = 2, set_origins = c(-1e308, 1e308)),
    "numerically representable")
})

test_that("EFRM simulator normalises very small ratios without underflow", {
  d <- simulate_efrm(n_per_group = 2, items_per_set = 2, n_sets = 2,
                     n_groups = 1, set_unit_ratio = 1e-300, seed = 713)
  expect_true(all(is.finite(attr(d, "truth")$alpha)))
  expect_true(all(attr(d, "truth")$alpha > 0))
})

test_that("planted-misfit selectors are plain vectors", {
  expect_error(
    simulate_rasch(40, 4, disordered = matrix("I01", 1),
                   model = "PCM"),
    "plain vectors")
  expect_error(
    simulate_mfrm(10, 3, 3,
      interaction = list(rater = matrix("R1", 1), item = "I1", bias = 1)),
    "each name one")
  expect_error(
    simulate_efrm(10, 3,
      item_drift = list(items = matrix("S1I01", 1), group = "g1", shift = 1)),
    "plain vector")
  expect_error(
    simulate_efrm(10, 3,
      item_drift = list(items = "S1I01", group = matrix("g1", 1), shift = 1)),
    "name one generated group")
})

test_that("simulation logits refuse finite inputs that overflow in combination", {
  expect_error(suppressWarnings(simulate_rasch(
    10, 2, model = "PCM", n_categories = 3,
    difficulty = rep(.Machine$double.xmax, 2), seed = 706
  )), "numerically representable")
})

test_that("a planted second dimension leaves indicators on both traits", {
  all_items <- sprintf("I%02d", 1:4)
  expect_error(
    simulate_rasch(
      40, 4, second_dim = list(items = all_items, rho = 0.2), seed = 714
    ),
    "at least one item on the primary trait"
  )
  expect_no_error(simulate_rasch(
    40, 4, second_dim = list(items = all_items[1:3], rho = 0.2), seed = 714
  ))
  expect_error(simulate_rasch(
    40, 4, second_dim = list(items = all_items[1:2], rho = 1), seed = 714
  ), "still use the primary trait")
})

test_that("planted item differences retain an invariant reference", {
  all_items <- sprintf("I%02d", 1:4)
  expect_error(simulate_rasch(
    40, 4, n_groups = 2,
    dif = list(items = all_items, uniform = 1), seed = 715
  ), "at least one invariant item")
  # A formally supplied zero effect is removed before this check because it
  # does not assert a DIF contrast.
  expect_no_error(simulate_rasch(
    40, 4, n_groups = 2,
    dif = list(items = all_items, uniform = 0), seed = 715
  ))

  set1 <- sprintf("S1I%02d", 1:3)
  expect_error(simulate_efrm(
    40, 3, n_sets = 2, n_groups = 2,
    item_drift = list(items = set1, group = "g2", shift = 1), seed = 716
  ), "invariant item in every affected set")
  expect_no_error(simulate_efrm(
    40, 3, n_sets = 2, n_groups = 2,
    item_drift = list(items = set1[1:2], group = "g2", shift = 1), seed = 716
  ))
  expect_error(simulate_efrm(
    20, 3, n_sets = 2, n_groups = 2,
    set_unit_ratio = .Machine$double.xmax,
    group_unit_ratio = .Machine$double.xmax, seed = 716
  ), "numerically representable")
})

test_that("Rasch recovery uses the generated logit unit", {
  d <- simulate_rasch(300, 6, discrimination = 2, seed = 717)
  tr <- attr(d, "truth")
  expect_equal(tr$calibration_truth$theta, 2 * tr$theta)
  expect_equal(tr$calibration_truth$difficulty, 2 * tr$difficulty)
  expect_equal(tr$calibration_truth$thresholds,
               lapply(tr$thresholds, function(z) 2 * z))
  r <- sim_recovery(rasch(d), d)
  expect_equal(r$pieces[["item difficulty"]]$true,
               as.numeric(2 * (tr$difficulty - mean(tr$difficulty))))
  expect_equal(r$pieces[["person ability"]]$true,
               as.numeric(2 * (tr$theta - mean(tr$theta))))
  expect_null(r$note)

  legacy <- d
  attr(legacy, "truth")$calibration_truth <- NULL
  r_legacy <- sim_recovery(rasch(legacy), legacy)
  expect_equal(r_legacy$pieces[["item difficulty"]]$true,
               r$pieces[["item difficulty"]]$true, tolerance = 1e-12)
  expect_equal(r_legacy$pieces[["person ability"]]$true,
               r$pieces[["person ability"]]$true, tolerance = 1e-12)

  misspecified <- simulate_rasch(
    300, 6, discrimination = c(0.6, 0.8, 1, 1.2, 1.4, 1.6), seed = 718
  )
  rm <- sim_recovery(rasch(misspecified), misspecified)
  expect_setequal(names(rm$pieces),
    c("item location (misspecified slopes)",
      "person ability (misspecified slopes)"))
  expect_match(rm$note, "not unbiased recovery targets")

  guessed <- simulate_rasch(300, 6, guessing = 0.15, seed = 719)
  rg <- sim_recovery(rasch(guessed), guessed)
  expect_match(rg$note, "guessing")
  expect_match(rg$note, "not unbiased recovery targets")

  speeded <- simulate_rasch(300, 6, speeded = 0.2, seed = 720)
  rs <- sim_recovery(rasch(speeded), speeded)
  # Independently assigned missing tails do not change the response model.
  expect_null(rs$note)

  explanatory <- simulate_rasch(200, 6, seed = 721)
  et <- attr(explanatory, "truth")
  et$departure_types <- "a fixed explanatory item departure"
  et$explanatory_departure <- list(target = "I03", level = "item")
  attr(explanatory, "truth") <- et
  free <- rasch(explanatory)
  expect_null(sim_recovery(free, explanatory)$note)
  restricted <- free
  restricted$explanatory <- list(relaxations = data.frame(
    item = character(0), component = character(0)))
  class(restricted) <- c("rasch_explanatory", "rasch")
  expect_match(sim_recovery(restricted, explanatory)$note,
               "fixed explanatory item departure")
  restricted$explanatory$relaxations <- data.frame(
    item = "I03", component = "Item location")
  expect_null(sim_recovery(restricted, explanatory)$note)
})

test_that("paired-comparison simulators retain every declared judge", {
  d <- simulate_btl(n_objects = 3, n_judges = 20, reps_per_pair = 7,
                    erratic_judges = 0.2, seed = 713)
  expect_setequal(unique(d$judge), sprintf("J%d", 1:20))
  expect_lte(diff(range(table(d$judge))), 1L)
  expect_true(all(attr(d, "truth")$erratic %in% d$judge))
  expect_error(simulate_btl(n_objects = 3, n_judges = 20,
                            reps_per_pair = 1, seed = 714),
               "fewer than the 20 requested judges")
})

test_that("the extra misfit types plant detectable signals", {
  # extreme response style: style persons over-use the end categories
  d <- simulate_rasch(600, 12, model = "PCM", n_categories = 4,
                      response_style = list(type = "extreme", prop = 0.3), seed = 1)
  si <- attr(d, "truth")$style_idx; cats <- as.matrix(d[, grep("^I", names(d))])
  expect_gt(mean(cats[si, ] %in% c(0, 3)), mean(cats[-si, ] %in% c(0, 3)) + 0.1)

  # speededness: a missing tail growing toward the last item
  d <- simulate_rasch(800, 15, speeded = 0.5, seed = 2)
  miss <- colMeans(is.na(as.matrix(d[, grep("^I", names(d))])))
  expect_lt(miss[8], 0.02)
  expect_gt(miss[15], 0.3)
  expect_true(miss[15] > miss[13] && miss[13] > miss[11])   # monotone gradient

  # MFRM halo: halo raters barely differentiate items -> large interaction
  d <- simulate_mfrm(140, 6, 6, rater_severity_sd = 0.5, item_sd = 1.3,
                     halo = 0.17, seed = 5)
  mf <- rasch_mfrm(d, person = "person", item = "item", score = "score",
                   facets = "rater", interaction = "rater")
  hr <- attr(d, "truth")$halo; ie <- mf$interaction_effects
  expect_gt(mean(abs(ie$gamma[ie$level %in% hr])),
            2 * mean(abs(ie$gamma[!ie$level %in% hr])))
})

test_that("sim_replicate and sim_recovery support Monte Carlo and recovery", {
  b <- sim_replicate(simulate_rasch, 6, n_persons = 300, n_items = 8, seed = 1)
  expect_s3_class(b, "rasch_sim_batch")
  expect_length(b, 6)
  expect_false(identical(b[[1]]$I01, b[[2]]$I01))          # different datasets

  # recovery: a clean fit gets its planted parameters back
  d <- simulate_rasch(600, 12, seed = 1)
  fit <- rasch(d)
  rec <- sim_recovery(fit, d)
  expect_s3_class(rec, "rasch_recovery")
  s <- rec$summary
  expect_gt(s$correlation[s$parameter == "item difficulty"], 0.95)
  # person ability is noisier (WLE precision from only 12 items limits it)
  expect_gt(s$correlation[s$parameter == "person ability"], 0.75)
  # bias is not identifiable for an origin-centred location parameter, so it
  # is reported NA rather than a structurally-zero value
  expect_true(is.na(s$bias[s$parameter == "item difficulty"]))
  pdf(NULL); on.exit(dev.off()); expect_no_error(plot_recovery(rec))

  bad <- fit
  bad$est$converged <- FALSE
  expect_error(sim_recovery(bad, d), "did not converge")

  # recovery across the other layouts
  d <- simulate_btl(8, 12, seed = 2)
  rb <- sim_recovery(btl(d, "object_a", "object_b", winner = "winner",
                         judge = "judge"), d)
  expect_gt(rb$summary$correlation[1], 0.9)

  # paired-comparison frames use their fitted identifying conventions: the
  # common object scale, geometric-mean-one panel units, and the first set as
  # the unit and origin reference
  d <- simulate_btl_efrm(5, 2, 5, 2, 8, 8,
    panel_units = c(1, 1.25), set_units = c(1, 1.3),
    set_origins = c(0, 0.4), seed = 71)
  tr <- attr(d, "truth")
  bf <- btl_efrm(d, "object_a", "object_b", winner = "winner",
                 judge = "judge", panels = "panel",
                 object_sets = tr$object_sets, se_method = "conditional",
                 boot_reps = 0)
  rf <- sim_recovery(bf, d)
  expect_setequal(rf$summary$parameter,
    c("object location", "panel unit (log)", "set unit (log)", "set origin"))
  expect_gt(rf$summary$correlation[
    rf$summary$parameter == "object location"], 0.9)

  # Set names are arbitrary labels: recovery follows the object partition,
  # not the simulator's conventional set1/set2 spelling.
  renamed <- bf
  set_names <- c(set1 = "Form A", set2 = "Form B")
  renamed$objects$set <- unname(set_names[as.character(renamed$objects$set)])
  renamed$alpha_table$set <- unname(set_names[as.character(
    renamed$alpha_table$set)])
  renamed$kappa_table$set <- unname(set_names[as.character(
    renamed$kappa_table$set)])
  rr <- sim_recovery(renamed, d)
  expect_setequal(rr$summary$parameter,
    c("object location", "panel unit (log)", "set unit (log)", "set origin"))

  # Panel names are presentation labels too.  The generating panel units are
  # matched through the judges who belong to each panel.
  panel_renamed <- bf
  panel_names <- c(panel1 = "Panel A", panel2 = "Panel B")
  panel_renamed$comparisons$panel <- unname(panel_names[
    panel_renamed$comparisons$panel])
  panel_renamed$phi_table$panel <- unname(panel_names[
    panel_renamed$phi_table$panel])
  rp <- sim_recovery(panel_renamed, d)
  expect_identical(rp$pieces[["panel unit (log)"]]$label,
                   unname(panel_names[names(tr$phi)]))
  unavailable <- bf
  unavailable$kappa_table$kappa[2L] <- NA_real_
  expect_error(sim_recovery(unavailable, d),
               "set origins cannot be matched to finite")
})

test_that("EFRM recovery follows set membership rather than set labels", {
  d <- simulate_efrm(20, 4, n_sets = 2, n_groups = 2, seed = 72)
  tr <- attr(d, "truth")
  set_of <- stats::setNames(
    rep(c("Form A", "Form B"), lengths(tr$item_sets)),
    unlist(tr$item_sets, use.names = FALSE))
  fit <- list(
    est = list(converged = TRUE),
    linking = list(alpha_edges = data.frame(converged = TRUE)),
    set_of = set_of,
    alpha_table = data.frame(set = c("Form A", "Form B"), alpha = tr$alpha),
    phi_table = data.frame(group = c("Group A", "Group B"), phi = tr$phi),
    person = data.frame(id = tr$person_id, theta = tr$theta),
    factors = data.frame(frame = factor(ifelse(
      as.character(tr$groups) == names(tr$phi)[1L], "Group A", "Group B"))),
    frame_group = "frame")
  class(fit) <- c("rasch_efrm", "rasch")
  r <- sim_recovery(fit, d)
  expect_setequal(r$summary$parameter,
                  c("set unit (log)", "group unit (log)"))
  expect_identical(r$pieces[["set unit (log)"]]$label,
                   c("Form A", "Form B"))
  expect_identical(r$pieces[["group unit (log)"]]$label,
                   c("Group A", "Group B"))

  repartitioned <- fit
  repartitioned$set_of[tr$item_sets[[1]][1]] <- "Form B"
  expect_error(sim_recovery(repartitioned, d),
               "partition does not match")
  unavailable <- fit
  unavailable$alpha_table$alpha[2L] <- NA_real_
  expect_error(sim_recovery(unavailable, d),
               "set units cannot be matched to finite")
})

test_that("sim_recovery verifies the fitted responses in every model family", {
  d1 <- simulate_rasch(120, 6, seed = 801)
  wrong_ids <- rasch(d1)
  wrong_ids$person$id[2L] <- wrong_ids$person$id[1L]
  expect_error(sim_recovery(wrong_ids, d1), "person allocation")

  d2 <- simulate_rasch(120, 6, seed = 802)
  expect_error(sim_recovery(rasch(d2), d1), "not fitted to the responses")
  subset_fit <- rasch(d1[, c("id", "I01", "I02", "I03", "I04")], id = "id")
  expect_error(sim_recovery(subset_fit, d1), "fitted items do not match")
  predictors <- data.frame(item = names(attr(d2, "truth")$difficulty),
                           position = seq_len(6L))
  fx <- rasch_explanatory(d2[predictors$item], predictors, ~ position)
  expect_error(sim_recovery(fx, d1), "not fitted to the responses")

  db <- simulate_btl(6, 8, reps_per_pair = 5, seed = 803)
  fb <- btl(db, "object_a", "object_b", winner = "winner", judge = "judge")
  db$winner[1L] <- if (db$winner[1L] == db$object_a[1L])
    db$object_b[1L] else db$object_a[1L]
  expect_error(sim_recovery(fb, db), "comparison, outcome or judge")

  dm <- simulate_mfrm(24, 4, 4, seed = 804)
  fm <- rasch_mfrm(dm, person = "person", item = "item", score = "score",
                   facets = "rater")
  dm$score[1L] <- (dm$score[1L] + 1L) %% 4L
  expect_error(sim_recovery(fm, dm), "rating or facet allocation")

  de <- simulate_efrm(100, 5, n_sets = 2, n_groups = 2, seed = 805)
  te <- attr(de, "truth")
  fe <- rasch_efrm(de, item_sets = te$item_sets, groups = "group", id = "id",
                   boot_reps = 0)
  ie <- te$item_sets[[1L]][1L]
  de[[ie]][1L] <- 1L - de[[ie]][1L]
  expect_error(sim_recovery(fe, de), "fitted response differs")

  dc <- simulate_btl_efrm(5, 2, 5, 2, 6, 6, seed = 806)
  tc <- attr(dc, "truth")
  fc <- btl_efrm(dc, "object_a", "object_b", winner = "winner",
                 judge = "judge", panels = "panel",
                 object_sets = tc$object_sets, se_method = "conditional",
                 boot_reps = 0)
  wrong_panel <- dc
  wrong_panel$panel[1L] <- setdiff(unique(dc$panel), dc$panel[1L])[1L]
  expect_error(sim_recovery(fc, wrong_panel),
               "judge/frame allocation")
  dc$winner[1L] <- if (dc$winner[1L] == dc$object_a[1L])
    dc$object_b[1L] else dc$object_a[1L]
  expect_error(sim_recovery(fc, dc), "comparison, outcome or judge")
})

test_that("sim_replicate preserves a seed at the integer boundary", {
  seen <- integer(0)
  generator <- function(seed) {
    seen <<- c(seen, seed)
    structure(data.frame(x = 1), truth = list(layout = "test"))
  }
  out <- sim_replicate(generator, 1, seed = .Machine$integer.max)
  expect_identical(seen, .Machine$integer.max)
  expect_length(out, 1L)
  expect_error(sim_replicate(generator, 2, seed = .Machine$integer.max),
               "exceeds the integer range")
})

test_that("person recovery matches generating truth by ID", {
  d <- simulate_rasch(500, 12, seed = 721)
  ord <- sample.int(nrow(d))
  shuffled <- d[ord, , drop = FALSE]
  fit <- rasch(shuffled, id = "id")
  rec <- sim_recovery(fit, shuffled)
  pp <- rec$pieces[["person ability"]]
  expect_identical(pp$label, as.character(fit$person$id)[
    is.finite(fit$person$theta)])
  expect_gt(rec$summary$correlation[
    rec$summary$parameter == "person ability"], 0.7)
})

test_that("recovery retains an externally anchored origin", {
  d <- simulate_rasch(600, 8, seed = 720)
  tr <- attr(d, "truth")
  anchors <- data.frame(item = names(tr$difficulty)[c(1, 8)], k = 1L,
                        tau = unname(tr$difficulty[c(1, 8)]))
  fit <- rasch(d, id = "id", anchors = anchors)
  rec <- sim_recovery(fit, d)
  ip <- rec$pieces[["item difficulty"]]
  pp <- rec$pieces[["person ability"]]
  expect_equal(ip$estimated,
               fit$items$location[match(ip$label, fit$items$item)])
  expect_equal(pp$true,
               unname(tr$theta[match(pp$label, tr$person_id)]))
  expect_equal(pp$estimated,
               fit$person$theta[match(pp$label, fit$person$id)])
  expect_true(all(is.finite(rec$summary$bias)))
  expect_equal(rec$summary$bias[rec$summary$parameter == "item difficulty"],
               mean(fit$items$location - unname(tr$difficulty[
                 match(fit$items$item, names(tr$difficulty))])))

  d <- simulate_btl(7, 20, 10, seed = 721)
  tr <- attr(d, "truth")
  fit <- btl(d, "object_a", "object_b", winner = "winner", judge = "judge",
             anchors = tr$location[c(1, 7)])
  rec <- sim_recovery(fit, d)
  expect_true(is.finite(rec$summary$bias[
    rec$summary$parameter == "object location"]))
})

test_that("the recovery plot refuses an unrelated object directly", {
  expect_error(plot_recovery(list()), "recovery result")
})

test_that("audit fixes hold: PCM structure, truth honesty, recovery centring", {
  # PCM and RSM now genuinely differ: per-item threshold patterns for PCM,
  # one common pattern for RSM; PCM thresholds stay ordered
  d1 <- simulate_rasch(200, 8, model = "PCM", n_categories = 4, seed = 42)
  d2 <- simulate_rasch(200, 8, model = "RSM", n_categories = 4, seed = 42)
  expect_false(identical(as.matrix(d1[, 2:9]), as.matrix(d2[, 2:9])))
  rel1 <- lapply(attr(d1, "truth")$thresholds, function(t) round(t - mean(t), 3))
  expect_gt(length(unique(rel1)), 1)                       # PCM varies
  rel2 <- lapply(attr(d2, "truth")$thresholds, function(t) round(t - mean(t), 3))
  expect_length(unique(rel2), 1)                           # RSM common
  expect_true(all(vapply(attr(d1, "truth")$thresholds,
                         function(t) !is.unsorted(t), TRUE)))
  expect_false("item-specific threshold structure" %in%
                 attr(d1, "truth")$departure_types)
  expect_null(sim_recovery(rasch(d1, model = "PCM"), d1)$note)
  expect_match(sim_recovery(rasch(d1, model = "RSM"), d1)$note,
               "item-specific threshold structure")

  # dif without groups is an error, not a silent no-op with a false truth
  expect_error(simulate_rasch(50, 6, dif = list(items = "I03", uniform = 1)),
               "n_groups")
  # polytomous guessing warns and is not recorded as planted
  expect_warning(d <- simulate_rasch(50, 4, model = "PCM", n_categories = 4,
                                     guessing = 0.3, seed = 1), "dichotomous")
  expect_false(any(grepl("guessing", attr(d, "truth")$planted)))
  # disordered thresholds work at 3 categories and preserve the location
  t2 <- rasch:::.sim_thresholds(0.5, 2, 1.2, disordered = TRUE)
  expect_true(is.unsorted(t2)); expect_equal(mean(t2), 0.5)
  rsm_all <- simulate_rasch(
    80, 4, model = "RSM", n_categories = 3,
    disordered = sprintf("I%02d", 1:4), seed = 2)
  expect_false("item-specific threshold structure" %in%
                 attr(rsm_all, "truth")$departure_types)
  rsm_one <- simulate_rasch(
    80, 4, model = "RSM", n_categories = 3,
    disordered = "I01", seed = 2)
  expect_true("item-specific threshold structure" %in%
                attr(rsm_one, "truth")$departure_types)
  # careless overwrite is not double-counted in the truth
  d <- simulate_rasch(300, 10, model = "PCM", n_categories = 4, careless = 0.5,
                      response_style = list(type = "extreme", prop = 0.5),
                      seed = 3)
  tr <- attr(d, "truth")
  expect_length(intersect(tr$style_idx, tr$careless_idx), 0)
  # requested proportions are either realised exactly or refused; one
  # mechanism cannot silently erase part of another
  expect_error(simulate_mfrm(30, 4, 5, erratic_raters = 0.4, halo = 0.8,
                             seed = 1), "cannot coexist")
  expect_error(simulate_mfrm(30, 4, 5, erratic_raters = 0.4, halo = 0.6,
                             seed = 1), "at least one ordinary rater")
  d <- simulate_mfrm(30, 4, 5, erratic_raters = 0.4, halo = 0.4, seed = 1)
  expect_length(attr(d, "truth")$erratic, 2L)
  expect_length(attr(d, "truth")$halo, 2L)

  # person ability is centred in recovery: an asymmetric difficulty range
  # must not masquerade as person-ability bias. Bias is not identifiable
  # up to the origin, so it is reported NA; the alignment shows instead as
  # a high correlation with no residual scale error
  d <- simulate_rasch(400, 10, difficulty = c(0, 3), seed = 2)
  r <- sim_recovery(rasch(d), d)
  expect_true(is.na(r$summary$bias[r$summary$parameter == "person ability"]))
  expect_gt(r$summary$correlation[r$summary$parameter == "person ability"], 0.75)
  # MFRM recovery reports item difficulties from the item margins
  d <- simulate_mfrm(60, 5, 5, seed = 1)
  mf <- rasch_mfrm(d, person = "person", item = "item", score = "score",
                   facets = "rater")
  r <- sim_recovery(mf, d)
  expect_true("item difficulty" %in% r$summary$parameter)
  expect_gt(r$summary$correlation[r$summary$parameter == "item difficulty"], 0.9)
})

test_that("MFRM recovery identifies the simulated rater facet", {
  d <- simulate_mfrm(12, 3, 4, seed = 73)
  tr <- attr(d, "truth")
  # A fitted object may contain additional facets and need not put the rater
  # facet first.  Its levels, rather than its position, identify the planted
  # severity parameter when the rater column has also been renamed.
  fit <- structure(list(
    est = list(converged = TRUE),
    facet_effects = list(
      occasion = data.frame(level = c("O1", "O2"), severity = c(-0.2, 0.2)),
      judge = data.frame(level = names(tr$severity),
                         severity = unname(tr$severity))
    ),
    item_effects = data.frame(item = names(tr$difficulty),
                              location = unname(tr$difficulty)),
    person = data.frame(id = tr$person_id, theta = tr$theta)
  ), class = c("rasch_mfrm", "rasch"))
  r <- sim_recovery(fit, d)
  rs <- r$pieces[["rater severity"]]
  expect_identical(rs$label, names(tr$severity))
  expect_equal(rs$estimated, as.numeric(tr$severity - mean(tr$severity)))

  ambiguous <- fit
  ambiguous$facet_effects$occasion <- ambiguous$facet_effects$judge
  expect_error(sim_recovery(ambiguous, d), "cannot be matched to one")

  misleading <- fit
  misleading$facet_effects$rater <- data.frame(
    level = c(names(tr$severity)[1], "unrelated"),
    severity = c(99, -99))
  r2 <- sim_recovery(misleading, d)
  expect_identical(r2$pieces[["rater severity"]]$label,
                   names(tr$severity))

  ordinary <- rasch(simulate_rasch(80, 4, seed = 74))
  expect_error(sim_recovery(ordinary, d), "does not match")
})

test_that("MFRM recovery verifies every fitted facet allocation", {
  d <- simulate_mfrm(20, 3, 3, seed = 731)
  p <- match(d$person, unique(d$person))
  i <- match(d$item, unique(d$item))
  r <- match(d$rater, unique(d$rater))
  d$occasion <- ifelse((p + i + r) %% 2L, "O1", "O2")
  fit <- rasch_mfrm(d, "person", "item", "score",
                    facets = c("occasion", "rater"))
  expect_s3_class(sim_recovery(fit, d), "rasch_recovery")
  changed <- d
  changed$occasion[1L] <- if (changed$occasion[1L] == "O1") "O2" else "O1"
  expect_error(sim_recovery(fit, changed), "facet allocation")

  renamed <- d
  renamed$judge <- renamed$rater
  renamed$rater <- NULL
  renamed_fit <- rasch_mfrm(renamed, "person", "item", "score",
                            facets = c("occasion", "judge"))
  expect_s3_class(sim_recovery(renamed_fit, renamed), "rasch_recovery")
})

test_that("recovery distinguishes fitted and unfitted generating departures", {
  bd <- simulate_btl(
    5, 8, 15, dependence = list(exposure = 0.3, carry_over = 0.2),
    seed = 75
  )
  b0 <- btl(bd, "object_a", "object_b", winner = "winner", judge = "judge")
  expect_match(sim_recovery(b0, bd)$note, "within-judge dependence")
  b1 <- btl(bd, "object_a", "object_b", winner = "winner", judge = "judge",
            order = "order")
  expect_null(sim_recovery(b1, bd)$note)

  md <- simulate_mfrm(
    60, 4, 4, halo = 0.25,
    interaction = list(rater = "R2", item = "I2", bias = 1.2), seed = 76
  )
  m0 <- rasch_mfrm(md, "person", "item", "score", facets = "rater")
  expect_match(sim_recovery(m0, md)$note, "halo ratings")
  m1 <- rasch_mfrm(md, "person", "item", "score", facets = "rater",
                   interaction = "rater")
  mr <- sim_recovery(m1, md)
  expect_null(mr$note)

  tr <- attr(md, "truth")
  surface <- outer(tr$difficulty, tr$severity, `+`)
  surface[, tr$halo] <- matrix(
    mean(tr$difficulty) + tr$severity[tr$halo], nrow(surface),
    length(tr$halo), byrow = TRUE)
  surface[tr$interaction$item, tr$interaction$rater] <-
    surface[tr$interaction$item, tr$interaction$rater] + tr$interaction$bias
  expected_items <- rowMeans(surface) - mean(surface)
  expected_raters <- colMeans(surface) - mean(surface)
  expect_equal(mr$pieces[["item difficulty"]]$true,
               unname(expected_items[mr$pieces[["item difficulty"]]$label]))
  expect_equal(mr$pieces[["rater severity"]]$true,
               unname(expected_raters[mr$pieces[["rater severity"]]$label]))
})

test_that("misfit layers compose: dependence and style respect DIF / 2nd dim", {
  # dependence regeneration keeps the target item's DIF
  d <- simulate_rasch(2000, 6, dif = list(items = "I02", uniform = 1.5),
                      n_groups = 2,
                      dependence = list(pairs = list(c("I01", "I02")),
                                        strength = 1), seed = 7)
  g <- attr(d, "truth")$groups
  expect_gt(mean(d$I02[g != "g2"]) - mean(d$I02[g == "g2"]), 0.12)
  # ...and the target item's second dimension
  d <- simulate_rasch(2000, 6, second_dim = list(items = "I02", rho = 0.2),
                      dependence = list(pairs = list(c("I01", "I02")),
                                        strength = 1), seed = 8)
  expect_lt(cor(d$I02, attr(d, "truth")$theta), 0.2)
  # response style keeps DIF for style-affected persons
  d <- simulate_rasch(3000, 6, model = "PCM", n_categories = 4,
                      dif = list(items = "I02", uniform = 2.5), n_groups = 2,
                      response_style = list(type = "extreme", prop = 0.5),
                      seed = 10)
  tr <- attr(d, "truth"); g <- tr$groups
  sty <- seq_len(nrow(d)) %in% tr$style_idx
  expect_gt(mean(d$I02[sty & g != "g2"]) - mean(d$I02[sty & g == "g2"]), 0.8)

  # Response style redraws the source response. The target must therefore be
  # regenerated from that styled source, not from the response it replaced.
  d <- simulate_rasch(
    4000, 6, model = "PCM", n_categories = 4,
    dependence = list(pairs = list(c("I01", "I02")), strength = 2),
    response_style = list(type = "extreme", prop = 0.5, strength = 2),
    seed = 771)
  tr <- attr(d, "truth")
  sty <- seq_len(nrow(d)) %in% tr$style_idx
  dep_coef <- function(rows) unname(stats::coef(stats::lm(
    d$I02[rows] ~ d$I01[rows] +
      stats::poly(tr$theta[rows], 5, raw = TRUE)))[2L])
  expect_gt(dep_coef(sty), 0.25)
  expect_gt(dep_coef(!sty), 0.20)
})

test_that("btl_dimensionality reference honours fitted dependence effects", {
  skip_on_cran()   # heavy simulation; verified locally and on CI
  # one-dimensional data whose only structure is within-judge order effects,
  # fitted WITH order: the dependence-aware reference must not read the
  # order structure as a second attribute
  verdicts <- vapply(1:4, function(s) {
    d <- simulate_btl(8, 12, reps_per_pair = 30,
                      dependence = list(exposure = 1, carry_over = 0.8),
                      seed = 300 + s)
    bt <- btl(d, "object_a", "object_b", winner = "winner", judge = "judge",
              order = "order")
    btl_dimensionality(bt, reps = 50,
                       independent_comparisons = TRUE)$leading_structured
  }, logical(1))
  # The simulator randomises order separately by judge, so this is an
  # identified test rather than four withheld results that happen to count
  # as FALSE through isTRUE().
  expect_false(anyNA(verdicts))
  expect_lte(sum(verdicts), 1L)   # was ~36% false-positive before the fix
})

test_that("simulated within-judge order is not shared across judges", {
  d <- simulate_btl(8, 12, reps_per_pair = 8,
                    dependence = list(exposure = 0.4, carry_over = 0.3),
                    seed = 411)
  variation <- rasch:::.btl_order_variation(
    d, sort(unique(c(d$object_a, d$object_b))))
  expect_true(variation$replicated)
  expect_false(variation$shared)
})

test_that("a seeded simulator call leaves the caller's RNG stream alone", {
  set.seed(99); before <- runif(3)
  set.seed(99); invisible(simulate_rasch(60, 5, seed = 7)); after <- runif(3)
  expect_equal(before, after)                     # stream not commandeered
  expect_identical(simulate_rasch(60, 5, seed = 7)$I01,
                   simulate_rasch(60, 5, seed = 7)$I01)   # still reproducible
  # and an absent stream is left absent rather than seeded behind the caller
  if (exists(".Random.seed", envir = globalenv(), inherits = FALSE))
    rm(".Random.seed", envir = globalenv())
  invisible(simulate_btl(n_objects = 4, n_judges = 3, reps_per_pair = 2,
                         seed = 3))
  expect_false(exists(".Random.seed", envir = globalenv(), inherits = FALSE))
})

test_that("every positive planted proportion affects at least one observation", {
  styled <- simulate_rasch(
    10, 4, model = "PCM", n_categories = 3,
    response_style = list(type = "extreme", prop = 0.01, strength = 1),
    seed = 41)
  expect_length(attr(styled, "truth")$style_idx, 1L)

  careless <- simulate_rasch(10, 4, careless = 0.01, seed = 42)
  expect_length(attr(careless, "truth")$careless_idx, 1L)

  incomplete <- simulate_rasch(10, 4, speeded = 0.01,
                               missing = 0.001, seed = 43)
  truth <- attr(incomplete, "truth")
  expect_length(truth$speeded_idx, 1L)
  expect_length(truth$missing_cells, 1L)
  expect_true(all(is.na(as.matrix(incomplete[, paste0("I0", 1:4)]))[
    truth$missing_cells]))
  expect_true(any(grepl("1 response missing", truth$planted, fixed = TRUE)))

  speeded_only <- simulate_rasch(20, 5, speeded = 1, seed = 430)
  both <- simulate_rasch(20, 5, speeded = 1, missing = 0.2, seed = 430)
  speeded_cells <- which(is.na(as.matrix(speeded_only[, paste0("I0", 1:5)])))
  both_truth <- attr(both, "truth")
  expect_length(both_truth$missing_cells, 20L)
  expect_length(intersect(speeded_cells, both_truth$missing_cells), 0L)
  expect_error(simulate_rasch(10, 3, speeded = 1, missing = 1, seed = 431),
               "finite value")

  erratic <- simulate_mfrm(10, 3, 2, erratic_raters = 0.01, seed = 44)
  expect_length(attr(erratic, "truth")$erratic, 1L)
  halo <- simulate_mfrm(10, 3, 2, halo = 0.01, seed = 45)
  expect_length(attr(halo, "truth")$halo, 1L)
  expect_error(simulate_mfrm(10, 3, 2, erratic_raters = 1,
                             halo = 0.01, seed = 45),
               "at least one rater")

  expect_error(simulate_rasch(
    10, 4, model = "PCM", n_categories = 3, careless = 1,
    response_style = list(type = "extreme", prop = 0.01, strength = 1),
    seed = 46), "at least one person")

  framed <- simulate_efrm(n_per_group = 2, items_per_set = 2,
                          n_sets = 2, n_groups = 2,
                          careless = 0.01, missing = 0.001, seed = 47)
  truth <- attr(framed, "truth")
  expect_length(truth$careless_idx, 1L)
  expect_length(truth$missing_cells, 1L)

  no_drift <- simulate_efrm(n_per_group = 4, items_per_set = 2,
                            n_sets = 2, n_groups = 2, seed = 48)
  zero_drift <- simulate_efrm(
    n_per_group = 4, items_per_set = 2, n_sets = 2, n_groups = 2,
    item_drift = list(items = "S1I01", group = "g2", shift = 0), seed = 48)
  expect_identical(as.data.frame(zero_drift), as.data.frame(no_drift))
  expect_null(attr(zero_drift, "truth")$item_drift)
  expect_error(simulate_efrm(
    n_per_group = 10, items_per_set = 3, n_sets = 2, n_groups = 1,
    item_drift = list(items = "S1I01", group = "g1", shift = 0.5), seed = 49),
    "requires n_groups >= 2")
})

test_that("contamination cannot erase every model-based observation unit", {
  expect_error(simulate_rasch(2, 3, careless = 0.8, seed = 51),
               "at least one person")
  expect_error(simulate_efrm(2, 2, careless = 1, seed = 52),
               "at least one person")
  expect_error(simulate_btl(3, 2, 3, erratic_judges = 0.8, seed = 53),
               "at least one judge")
  expect_error(simulate_btl_efrm(
    3, 2, 2, 2, 2, 2, erratic_judges = 1, seed = 54
  ), "at least one judge")
  expect_error(simulate_mfrm(10, 3, 2, erratic_raters = 0.8, seed = 55),
               "at least one rater")
  expect_error(simulate_mfrm(10, 3, 2, halo = 0.8, seed = 56),
               "at least one ordinary rater")
  expect_error(simulate_mfrm(10, 3, 3, item_sd = 0, halo = 0.3, seed = 56),
               "halo requires item_sd > 0")
  expect_error(simulate_rasch(
    10, 3, model = "PCM", n_categories = 3,
    response_style = list(type = "extreme", prop = 1), seed = 56
  ), "absorbed by the fitted thresholds")

  # At rho = 1 the response data are exactly those from the primary object
  # locations, so labelling the same vector as a second attribute is refused.
  expect_error(simulate_btl(
    4, 4, 3, second_attribute = list(rho = 1), seed = 57
  ), "finite value in")
  reversed <- simulate_btl(
    4, 4, 3, second_attribute = list(rho = -1), seed = 57
  )
  expect_true(any(grepl("reversal stress", attr(reversed, "truth")$planted,
                        fixed = TRUE)))
})

test_that("curve-shape departures require realised person variation", {
  expect_error(simulate_rasch(
    20, 3, theta_sd = 0, discrimination = c(0.7, 1, 1.3), seed = 62
  ), "heterogeneous discrimination requires theta_sd > 0")
  expect_error(simulate_rasch(
    20, 3, theta_sd = 0, guessing = 0.2, seed = 62
  ), "guessing requires theta_sd > 0")
  expect_error(simulate_rasch(
    20, 3, theta_sd = 0, n_groups = 2,
    dif = list(items = "I02", nonuniform = 0.5), seed = 62
  ), "non-uniform DIF requires theta_sd > 0")
  expect_no_error(simulate_rasch(
    20, 3, theta_sd = 0, discrimination = 2, n_groups = 2,
    dif = list(items = "I02", uniform = 0.5), seed = 62
  ))
})

test_that("contamination retains every group, panel and attribute camp", {
  expect_error(simulate_rasch(
    6, 3, n_groups = 2, careless = 0.83, seed = 58
  ), "in every generated group")
  grouped <- simulate_rasch(6, 3, n_groups = 2, careless = 0.67, seed = 58)
  g <- attr(grouped, "truth")$groups
  bad <- attr(grouped, "truth")$careless_idx
  expect_true(all(vapply(levels(g), function(z)
    any(!which(g == z) %in% bad), logical(1))))

  expect_error(simulate_efrm(
    2, 2, n_sets = 2, n_groups = 2, careless = 0.75, seed = 59
  ), "in every generated group")
  one_group <- simulate_efrm(
    2, 2, n_sets = 1, n_groups = 1, careless = 0.5, seed = 59
  )
  expect_length(attr(one_group, "truth")$careless_idx, 1L)
  expect_true(attr(one_group, "truth")$careless_idx %in% 1:2)
  expect_error(simulate_btl(
    4, 4, 2, second_attribute = list(rho = 0),
    erratic_judges = 0.75, seed = 60
  ), "each second-attribute camp")
  expect_error(simulate_btl_efrm(
    3, 2, 2, 2, 2, 2, erratic_judges = 0.75, seed = 61
  ), "in every panel")
  framed <- simulate_btl_efrm(
    3, 2, 2, 2, 2, 2, erratic_judges = 0.5, seed = 61
  )
  tr <- attr(framed, "truth")
  expect_true(all(vapply(unique(tr$judge_panel), function(z)
    any(!names(tr$judge_panel)[tr$judge_panel == z] %in% tr$erratic),
    logical(1))))
})

test_that("Rasch simulation preserves numbered group order", {
  x <- simulate_rasch(100, 4, n_groups = 10,
    dif = list(items = 1, uniform = 1), seed = 930)
  expect_identical(levels(x$group), paste0("g", 1:10))
  expect_true(any(grepl("DIF (group g10)", attr(x, "truth")$planted,
                       fixed = TRUE)))
  speeded <- simulate_rasch(100, 4, speeded = 0.2, seed = 930)
  expect_length(attr(speeded, "truth")$departure_types, 0L)
  expect_length(attr(speeded, "truth")$speeded_idx, 20L)
})

test_that("response-data simulators cannot remove every observation", {
  expect_error(simulate_rasch(20, 5, missing = 1), "finite value")
  expect_error(simulate_efrm(10, 3, missing = 1), "finite value")
  expect_error(simulate_rasch(2, 2, missing = 0.99), "every remaining response")
  expect_error(simulate_efrm(2, 2, missing = 0.99), "every response")
})

test_that("EFRM simulation refuses an unidentified cross-set scale", {
  expect_error(
    simulate_efrm(
      n_per_group = 20, items_per_set = 3, n_sets = 2,
      n_groups = 2, theta_sd = 0
    ),
    "relative set units require person variation"
  )
  expect_no_error(simulate_efrm(
    n_per_group = 20, items_per_set = 3, n_sets = 1,
    n_groups = 2, theta_sd = 0
  ))
})
