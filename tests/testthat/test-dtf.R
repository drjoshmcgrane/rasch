# Bundle DIF through dif_anova(bundles = ) and test- or bundle-level
# differential functioning through dtf() on a resolved calibration.

dtf_sim <- function(seed, n = 800, L = 10, shift = NULL) {
  set.seed(seed)
  d <- seq(-1.5, 1.5, length.out = L)
  g <- rep(c("a", "b"), each = n / 2)
  sh <- matrix(0, n, L)
  for (nm in names(shift)) sh[g == "b", as.integer(nm)] <- shift[[nm]]
  X <- matrix(rbinom(n * L, 1, plogis(outer(rnorm(n), d, "-") - sh)), n, L)
  colnames(X) <- paste0("I", seq_len(L))
  list(fit = rasch(data.frame(X, grp = g), factors = "grp"), g = g)
}

test_that("a bundle shift too small to flag item by item flags the bundle", {
  s <- dtf_sim(3, L = 12, shift = list("2" = 0.35, "3" = 0.35, "4" = 0.35))
  da <- dif_anova(s$fit, bundles = list(passage = c("I2", "I3", "I4"),
                                        other = c("I9", "I10", "I11")))
  sm <- da$summary
  expect_false(any(sm$uniform_DIF[!sm$item %in% c("passage", "other")]))
  expect_true(sm$uniform_DIF[sm$item == "passage"])
  expect_false(sm$uniform_DIF[sm$item == "other"])
  expect_identical(names(da$bundles), c("passage", "other"))
  expect_identical(da$bootstrap_design$bundles, da$bundles)
  expect_match(paste(da$notes, collapse = " "), "passage = {I2, I3, I4}",
               fixed = TRUE)
  out <- paste(capture.output(print(da)), collapse = "\n")
  expect_match(out, "Bundles: passage = {I2, I3, I4}; other = {I9, I10, I11}",
               fixed = TRUE)
  # the bundle joins the adjustment family with the items
  expect_equal(nrow(sm), 14L)
  # a flagged bundle has no post-hoc follow-up but says so
  ds <- dif_anova(s$fit, bundles = list(passage = c("I2", "I3", "I4")),
                  sizes = TRUE)
  expect_false("passage" %in% ds$posthoc$item)
  expect_match(paste(ds$notes, collapse = " "), "no post-hoc comparison")
  # a plain analysis carries no bundles
  expect_null(dif_anova(s$fit)$bundles)
})

test_that("bundles are validated against the items", {
  s <- dtf_sim(4, n = 300, L = 6)
  f <- s$fit
  expect_error(dif_anova(f, bundles = list(c("I2", "I3"))), "needs a name")
  expect_error(dif_anova(f, bundles = list(I1 = c("I2", "I3"))),
               "also item names: I1")
  expect_error(dif_anova(f, bundles = list(one = "I2")), "at least two")
  expect_error(dif_anova(f, bundles = list(all = paste0("I", 1:6))),
               "whole test")
  expect_error(dif_anova(f, bundles = list(x = c("I2", "Z"))),
               "not in the analysis: Z")
  expect_error(dif_anova(f, bundles = list(a = c("I2", "I3"),
                                           a = c("I4", "I5"))),
               "more than once: a")
  expect_error(dif_anova(f, bundles = "I2"), "named list")
})

test_that("bundle terms replay through the DIF bootstrap", {
  s <- dtf_sim(5, n = 200, L = 6)
  da <- dif_anova(s$fit, bundles = list(pair = c("I2", "I3")))
  db <- suppressWarnings(dif_bootstrap(s$fit, da, B = 2, workers = 1,
                                       seed = 9))
  expect_equal(db$B_used, 2L)
  expect_true("pair" %in% db$summary$item)
  expect_true(all(db$summary$n_boot_uniform[db$summary$item == "pair"] == 2L))
})

test_that("dtf recovers a planted test-level shift with delta-method errors", {
  s <- dtf_sim(11, n = 1000, shift = list("2" = 0.8, "3" = 0.8))
  d <- dtf(s$fit, by = "grp", items = c("I2", "I3"),
           bundles = list(pair = c("I2", "I3"), none = c("I5", "I6", "I7")))
  expect_s3_class(d, "rasch_dtf")
  expect_identical(d$anchors, paste0("I", c(1, 4:10)))
  expect_identical(d$groups, "b")
  expect_identical(d$reference, "a")
  it <- d$items
  expect_identical(it$split, it$item %in% c("I2", "I3"))
  expect_true(all(it$shift[!it$split] == 0))
  expect_true(all(it$se[!it$split] == 0))
  expect_true(all(is.na(it$p[!it$split])))
  sp <- it[it$split, ]
  expect_true(all(abs(sp$shift - 0.8) < 3 * sp$se))
  expect_true(all(sp$significant))
  # the mean shift is the planted shifts spread over the test
  tt <- d$test
  expect_true(abs(tt$shift_mean - 0.16) < 3 * tt$se)
  expect_equal(tt$shift_mean, mean(it$shift))
  expect_true(tt$p < 0.001)
  # harder for b everywhere, so the signed and unsigned summaries agree
  expect_true(all(d$scores$shift > 0))
  expect_equal(tt$sDTF_logit, tt$uDTF_logit)
  expect_equal(tt$sDTF_score, tt$uDTF_score)
  expect_equal(tt$sDTF_pct, 100 * tt$sDTF_score / 10)
  expect_true(tt$max_score >= tt$uDTF_score)
  cu <- d$curves
  expect_true(all(cu$shift_score > 0) && all(cu$shift_logit > 0))
  expect_equal(cu$shift_score, cu$expected_reference - cu$expected_group)
  # the score-to-measure difference at a score is the logit-shift curve at
  # the group's measure for that score
  r <- 5L
  sc <- d$scores[d$scores$score == r, ]
  f2 <- split_items(s$fit, c("I2", "I3"), by = "grp")
  idx_b <- match(c("I1", "I2 (b)", "I3 (b)", paste0("I", 4:10)),
                 f2$items$item)
  idx_a <- match(c("I1", "I2 (a)", "I3 (a)", paste0("I", 4:10)),
                 f2$items$item)
  expected <- function(th, idx) sum(vapply(f2$tau_list[idx], function(tt)
    item_moments(th, tt)$E, 0))
  expect_equal(expected(sc$theta_group, idx_b), r, tolerance = 1e-6)
  expect_equal(expected(sc$theta_reference, idx_a), r, tolerance = 1e-6)
  # the delta-method standard error agrees with finite differences
  measure <- function(f, idx) uniroot(function(t)
    expected_f(f, t, idx) - r, c(-10, 10), tol = 1e-10)$root
  expected_f <- function(f, th, idx) sum(vapply(f$tau_list[idx],
    function(tt) item_moments(th, tt)$E, 0))
  base <- measure(f2, idx_b) - measure(f2, idx_a)
  expect_equal(base, sc$shift, tolerance = 1e-6)
  grad <- vapply(seq_len(nrow(f2$thresholds)), function(j) {
    fj <- f2; ii <- f2$thresholds$item[j]; kk <- f2$thresholds$k[j]
    fj$tau_list[[ii]][kk] <- fj$tau_list[[ii]][kk] + 1e-5
    (measure(fj, idx_b) - measure(fj, idx_a) - base) / 1e-5
  }, 0)
  se_num <- sqrt(drop(t(grad) %*% f2$est$cov_tau %*% grad))
  expect_equal(sc$se, se_num, tolerance = 1e-3)
  # bundles: the pair shifts homogeneously; the anchors' bundle is zero
  b <- d$bundles
  pair <- b[b$bundle == "pair", ]
  expect_true(abs(pair$shift_mean - 0.8) < 3 * pair$se)
  expect_equal(pair$df_hom, 1)
  expect_true(pair$p_hom > 0.05)
  expect_true(pair$significant)
  expect_equal(pair$sDBF_pct, 100 * pair$sDBF_score / 2)
  none <- b[b$bundle == "none", ]
  expect_equal(none$shift_mean, 0)
  expect_true(is.na(none$chisq_hom))
  expect_equal(none$uDBF_score, 0)
  out <- paste(capture.output(print(d)), collapse = "\n")
  expect_match(out, "8 of 10 items anchor the groups", fixed = TRUE)
  expect_match(out, "Positive values: harder for the group", fixed = TRUE)
  expect_match(out, "Bundles:", fixed = TRUE)
  expect_match(out, "Split items:", fixed = TRUE)
})

test_that("shifts that cancel leave a signed difference near zero", {
  s <- dtf_sim(12, n = 1000, shift = list("2" = 0.9, "8" = -0.9))
  d <- dtf(s$fit, by = "grp", items = c("I2", "I8"),
           bundles = list(both = c("I2", "I8")))
  tt <- d$test
  expect_true(abs(tt$sDTF_logit) < 3 * tt$se_sDTF_logit)
  expect_true(tt$uDTF_logit > 3 * tt$se_uDTF_logit)
  expect_true(abs(tt$sDTF_score) < tt$uDTF_score)
  expect_true(abs(tt$shift_mean) < 3 * tt$se)
  sp <- d$items[d$items$split, ]
  expect_true(sp$shift[sp$item == "I2"] > 0 && sp$shift[sp$item == "I8"] < 0)
  both <- d$bundles
  expect_true(both$p_hom < 0.001)
  expect_true(abs(both$shift_mean) < 3 * both$se)
})

test_that("dtf takes its grouping from a resolution and measures several groups", {
  s <- dtf_sim(13, n = 900, shift = list("2" = 1))
  res <- resolve_dif(s$fit)
  expect_true("I2" %in% res$splits$item)
  d <- dtf(res)
  expect_identical(d$by, "grp")
  expect_true(all(d$items$split == (d$items$item == "I2")))
  expect_true(d$test$shift_mean > 0)
  expect_error(dtf(res, items = "I3"), "already split")
  # nothing split: every difference is zero by construction
  s0 <- dtf_sim(14, n = 400, L = 8)
  res0 <- resolve_dif(s0$fit)
  expect_equal(res0$n_splits, 0L)
  d0 <- dtf(res0)
  expect_equal(d0$test$shift_mean, 0)
  expect_equal(d0$test$uDTF_score, 0)
  expect_true(is.na(d0$test$p))
  expect_true(all(d0$scores$shift == 0))
  expect_match(paste(d0$notes, collapse = " "), "zero by construction")
  # three groups, a grouping vector, and a chosen reference
  set.seed(15); n <- 900; L <- 8
  g3 <- rep(c("x", "y", "z"), each = n / 3)
  dd <- seq(-1.5, 1.5, length.out = L)
  sh <- matrix(0, n, L); sh[g3 == "z", 1:2] <- 0.7
  X <- matrix(rbinom(n * L, 1, plogis(outer(rnorm(n), dd, "-") - sh)), n, L)
  colnames(X) <- paste0("I", 1:L)
  f3 <- rasch(data.frame(X, g = g3), factors = "g")
  f3 <- split_items(f3, c("I1", "I2"), by = "g")
  d3 <- dtf(f3, by = g3, reference = "x")
  expect_identical(d3$groups, c("y", "z"))
  expect_identical(d3$by, "group")
  tz <- d3$test[d3$test$group == "z", ]
  ty <- d3$test[d3$test$group == "y", ]
  expect_true(tz$p < 0.01)
  expect_true(abs(ty$shift_mean) < 3 * ty$se)
  expect_equal(nrow(d3$items), 2L * L)
  expect_equal(sum(is.finite(d3$items$p_adj)), 4L)
})

test_that("dtf refuses what it cannot compare", {
  s <- dtf_sim(16, n = 300, L = 6)
  f <- s$fit
  expect_error(dtf(f), "no item is split")
  expect_error(dtf(f, by = "grp"), "no item is split")
  f2 <- split_items(f, "I2", by = "grp")
  expect_error(dtf(f2, by = "grp", items = "I3"), "already carries split")
  expect_error(dtf(f2, by = "grp", reference = "q"), "one level")
  expect_error(dtf(f2, by = "nope"), "not a person factor")
  expect_error(dtf(f2, by = rep(c("u", "v"), length.out = 300)),
               "several copies answered by group 'u'")
  expect_error(dtf(f2, by = "grp", bundles = list(one = "I2")),
               "at least two")
  expect_error(dtf(f2, by = "grp", bundles = list(x = c("I2", "Q"))),
               "not calibrated for every group: Q")
  expect_error(dtf(f2, by = "grp", grid = 2), "grid")
  expect_error(dtf(list()), "ordinary rasch fit")
})

test_that("an item without a copy for every group leaves the comparison", {
  set.seed(22); n <- 900; th <- rnorm(n)
  g <- rep(c("a", "b", "c"), each = n / 3)
  tau <- lapply(1:6, function(i) c(-0.8, 0.8) + seq(-1, 1, length.out = 6)[i])
  X <- sapply(1:6, function(i) vapply(seq_len(n), function(p) {
    lp <- c(0, cumsum(th[p] - tau[[i]] - (g[p] == "b" && i == 2) * 0.6))
    pr <- exp(lp - max(lp)); sample(0:2, 1, prob = pr / sum(pr)) }, 0))
  colnames(X) <- paste0("P", 1:6)
  X[g == "c" & X[, 1] == 2, 1] <- 1
  fp <- rasch(data.frame(X, g = g), factors = "g")
  f3 <- split_items(fp, c("P1", "P2"), by = "g")
  expect_false("P1 (c)" %in% f3$items$item)
  d <- dtf(f3, by = "g")
  expect_identical(d$dropped, "P1")
  expect_equal(d$n_items, 5L)
  expect_equal(d$max_score, 10L)
  expect_match(paste(d$notes, collapse = " "), "left out of the comparison: P1")
  tb <- d$test[d$test$group == "b", ]
  expect_true(tb$p < 0.01)
  expect_true(all(d$items$item != "P1"))
})

test_that("plot_dtf draws a group and refuses an unknown one", {
  s <- dtf_sim(17, n = 300, L = 6, shift = list("2" = 0.8))
  d <- dtf(s$fit, by = "grp", items = "I2")
  pdf(NULL)
  on.exit(dev.off())
  expect_identical(plot_dtf(d), d)
  expect_error(plot_dtf(d, group = "zz"), "must be one of: b")
})
