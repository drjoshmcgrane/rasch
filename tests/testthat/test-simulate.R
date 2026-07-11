# Data simulation: each planted departure must trip its matching diagnostic
# (the plant -> detect loop), and the truth is attached for recovery checks.

test_that("simulate_rasch plants misfit the Rasch diagnostics detect", {
  # discrimination: a central over-discriminating item overfits (low outfit)
  d <- simulate_rasch(600, 11, discrimination = c(rep(1, 5), 3, rep(1, 5)),
                      seed = 1)
  f <- rasch(d)
  expect_lt(f$items$outfit_ms[6], 0.7)
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
  h <- residual_correlations(rasch(d))$flagged
  expect_true(any((h$item_a == "I03" & h$item_b == "I04") |
                  (h$item_a == "I04" & h$item_b == "I03")))

  # careless responders inflate person outfit
  d <- simulate_rasch(600, 15, careless = 0.12, seed = 6)
  po <- rasch(d)$person$outfit_ms; ci <- attr(d, "truth")$careless_idx
  expect_gt(mean(po[ci], na.rm = TRUE), mean(po[-ci], na.rm = TRUE) + 0.5)

  expect_output(print(d), "careless")
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

test_that("simulate_efrm plants a frame-unit ratio rasch_efrm recovers", {
  d <- simulate_efrm(300, 8, set_unit_ratio = 1.35, seed = 2)
  tr <- attr(d, "truth")
  ef <- rasch_efrm(d, item_sets = tr$item_sets, groups = "group")
  ratio <- max(ef$alpha_table$alpha) / min(ef$alpha_table$alpha)
  expect_gt(ratio, 1.2); expect_lt(ratio, 1.55)          # ~1.35 recovered
  expect_output(print(d), "set-unit ratio")
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
  rec <- sim_recovery(rasch(d), d)
  expect_s3_class(rec, "rasch_recovery")
  s <- rec$summary
  expect_gt(s$correlation[s$parameter == "item difficulty"], 0.95)
  # person ability is noisier (WLE precision from only 12 items limits it)
  expect_gt(s$correlation[s$parameter == "person ability"], 0.75)
  expect_lt(abs(s$bias[s$parameter == "item difficulty"]), 0.1)
  pdf(NULL); on.exit(dev.off()); expect_no_error(plot_recovery(rec))

  # recovery across the other layouts
  d <- simulate_btl(8, 12, seed = 2)
  rb <- sim_recovery(btl(d, "object_a", "object_b", winner = "winner",
                         judge = "judge"), d)
  expect_gt(rb$summary$correlation[1], 0.9)
})
