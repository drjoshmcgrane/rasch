# dif_wald(): the conditional Wald test of DIF on the resolved calibration,
# and resolve_dif(criterion = "wald") built on it.

wald_sim <- function(seed, n = 600, L = 8, shift = NULL) {
  set.seed(seed)
  d <- seq(-2, 2, length.out = L)
  g <- rep(c("a", "b"), each = n / 2)
  sh <- matrix(0, n, L)
  if (!is.null(shift))
    for (j in names(shift)) sh[g == "b", as.integer(j)] <- shift[[j]]
  th <- rnorm(n)
  X <- matrix(rbinom(n * L, 1, plogis(outer(th, d, "-") - sh)), n, L)
  colnames(X) <- paste0("I", 1:L)
  list(fit = rasch(data.frame(X, grp = g), factors = "grp"), g = g, X = X)
}

test_that("dif_wald recovers a planted shift and agrees with dif_size", {
  s <- wald_sim(4, shift = list("3" = 0.8))
  w <- dif_wald(s$fit)
  expect_s3_class(w, "rasch_dif_wald")
  sm <- w$summary
  expect_identical(names(sm), c("item", "factor", "n_levels", "shift", "se",
                                "wald", "df", "ref_df", "p", "p_adj",
                                "significant"))
  expect_identical(sm$item, paste0("I", 1:8))
  expect_true(all(sm$factor == "grp"))
  expect_true(all(sm$n_levels == 2L))
  expect_true(all(sm$df == 1))
  expect_true(all(is.infinite(sm$ref_df)))
  expect_true(all(is.finite(sm$shift)) && all(is.finite(sm$se)))
  r3 <- sm[sm$item == "I3", ]
  expect_true(r3$significant)
  expect_true(abs(r3$shift - 0.8) < 3 * r3$se)
  expect_true(r3$shift > 0)
  expect_false(any(sm$significant[sm$item != "I3"]))
  # the statistic is the squared ratio and the reference is chi-square
  expect_equal(sm$wald, (sm$shift / sm$se)^2)
  expect_equal(sm$p, stats::pchisq(sm$wald, 1, lower.tail = FALSE))
  expect_equal(sm$p_adj, stats::p.adjust(sm$p, "holm"))
  # the per-level table carries the resolved locations behind the shift
  lv <- w$levels
  expect_equal(nrow(lv), 16L)
  expect_true(all(lv$n == 300L))
  l3 <- lv[lv$item == "I3", ]
  expect_equal(diff(l3$location), r3$shift)
  expect_true(all(is.finite(lv$se)))
  # dif_size resolves the same item on the same anchor
  ds <- dif_size(s$fit, "I3", by = "grp")
  expect_equal(ds$pairs$difference, -r3$shift)
  expect_equal(ds$pairs$se, r3$se)
  expect_equal(ds$pairs$t^2, r3$wald)
  expect_length(w$notes, 0L)
  out <- capture.output(print(w))
  expect_match(paste(out, collapse = "\n"), "Conditional Wald test of DIF by grp")
  expect_match(paste(out, collapse = "\n"), "chi-square reference")
  # a subset of items restricts the table and the adjustment family
  w1 <- dif_wald(s$fit, items = c("I3", "I5"))
  expect_identical(w1$summary$item, c("I3", "I5"))
  expect_equal(w1$summary$p_adj, stats::p.adjust(w1$summary$p, "holm"))
})

test_that("dif_wald leaves out split copies and withholds an anchored item", {
  s <- wald_sim(4, shift = list("3" = 0.8))
  f3 <- split_items(s$fit, "I3", by = "grp")
  w <- dif_wald(f3)
  expect_false(any(grepl("I3", w$summary$item)))
  expect_equal(nrow(w$summary), 7L)
  expect_match(paste(w$notes, collapse = " "),
               "grp: not tested because persons in one level only answered: I3 \\(a\\), I3 \\(b\\)")
  expect_false(any(w$summary$significant))
  anc <- data.frame(item = "I4", k = 1L, tau = 0)
  fa <- rasch(data.frame(s$X, grp = s$g), factors = "grp", anchors = anc)
  wa <- dif_wald(fa, items = c("I4", "I5"))
  r4 <- wa$summary[wa$summary$item == "I4", ]
  expect_true(is.na(r4$shift) && is.na(r4$wald) && is.na(r4$p) &&
                is.na(r4$significant))
  expect_true(is.finite(wa$summary$p[wa$summary$item == "I5"]))
  expect_match(paste(wa$notes, collapse = " "), "an anchored item cannot be split")
})

test_that("dif_wald handles three levels of a polytomous factor", {
  set.seed(22); n <- 900; th <- rnorm(n)
  g <- rep(c("a", "b", "c"), each = n / 3)
  tau <- lapply(1:6, function(i) c(-0.8, 0.8) + seq(-1, 1, length.out = 6)[i])
  X <- sapply(1:6, function(i) vapply(seq_len(n), function(p) {
    lp <- c(0, cumsum(th[p] - tau[[i]] - (g[p] == "b" && i == 2) * 0.6))
    pr <- exp(lp - max(lp)); sample(0:2, 1, prob = pr / sum(pr)) }, 0))
  colnames(X) <- paste0("P", 1:6)
  X[g == "c" & X[, 1] == 2, 1] <- 1
  fp <- rasch(data.frame(X, g = g), factors = "g")
  w <- dif_wald(fp)
  sm <- w$summary
  expect_true(all(sm$n_levels == 3L))
  expect_true(all(sm$df == 2))
  expect_true(all(is.na(sm$se)))
  r2 <- sm[sm$item == "P2", ]
  expect_true(r2$significant)
  l2 <- w$levels[w$levels$item == "P2", ]
  expect_identical(l2$level, c("a", "b", "c"))
  expect_equal(r2$shift, diff(range(l2$location)))
  expect_true(which.max(l2$location) == 2L)
  expect_equal(r2$p, stats::pf(r2$wald / 2, 2, Inf, lower.tail = FALSE))
  # a group that lost a score category has no common structure to compare
  r1 <- sm[sm$item == "P1", ]
  expect_true(is.na(r1$wald) && is.na(r1$p))
  expect_match(paste(w$notes, collapse = " "),
               "P1: resolved contrasts withheld because groups have different observed response-category structures")
  expect_match(paste(capture.output(print(w)), collapse = "\n"),
               "notes:")
})

test_that("dif_wald refuses bad input", {
  s <- wald_sim(1)
  expect_error(dif_wald(s$fit, items = "Z"), "item\\(s\\) not in the fit: Z")
  expect_error(dif_wald(s$fit, items = character()), "at least one item")
  expect_error(dif_wald(s$fit, factors = "sex"), "not stored in the fit")
  expect_error(dif_wald(list()), "ordinary rasch fit")
  expect_error(dif_wald(s$fit, alpha = 2), "alpha")
})

test_that("resolve_dif with the Wald criterion splits both planted items", {
  s <- wald_sim(1, shift = list("2" = 1, "6" = -1))
  r <- resolve_dif(s$fit, criterion = "wald")
  expect_s3_class(r, "rasch_resolve_dif")
  expect_identical(r$criterion, "wald")
  expect_identical(r$algorithm, "factor-design-resolution-3")
  expect_setequal(r$splits$item, c("I2", "I6"))
  expect_true(all(is.na(r$splits$eta2)))
  expect_true(all(is.finite(r$splits$magnitude)))
  expect_true(all(r$splits$magnitude > 0.5))
  expect_equal(r$n_remaining_dif, 0L)
  expect_true(is.na(r$n_nonuniform))
  expect_equal(r$n_untested, 0L)
  expect_null(r$dif)
  expect_match(r$stopped, "no significant DIF remains")
  expect_true(all(c("I2 (a)", "I2 (b)", "I6 (a)", "I6 (b)") %in% r$fit$items$item))
  out <- paste(capture.output(print(r)), collapse = "\n")
  expect_match(out, "Non-uniform DIF is not tested by the Wald criterion")
  expect_match(out, "2 split\\(s\\)")
  # the criterion tests each factor on its own
  expect_error(resolve_dif(s$fit, criterion = "wald", effects = "factorial"),
               "effects = \"main\"")
  # no DIF planted: nothing to split, and the default criterion is unchanged
  s0 <- wald_sim(4)
  r0 <- resolve_dif(s0$fit, criterion = "wald")
  expect_equal(r0$n_splits, 0L)
  expect_equal(r0$n_remaining_dif, 0L)
  ra <- resolve_dif(s0$fit)
  expect_identical(ra$criterion, "anova")
  expect_equal(ra$n_splits, 0L)
})
