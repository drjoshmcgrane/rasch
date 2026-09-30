# Regressions for the round-five DIF fix: a model that fits every residual
# mean exactly leaves rounding noise, not zero, as its residual sum of
# squares, and every test built on that residual is withheld with the reason
# rather than computed from the noise.

# Two groups by four class intervals, every person in a cell sharing that
# cell's residual mean, so the factorial design fits the response exactly.
sim_exact_cells <- function() {
  d <- data.frame(pid = sprintf("p%03d", 1:40),
                  ci = factor(rep(1:4, each = 10)),
                  f1 = factor(rep(c("a", "b"), 20)))
  cell <- c(a1 = 0.1, b1 = -0.1, a2 = 0.37, b2 = -0.37,
            a3 = 0.23, b3 = -0.23, a4 = 0.5, b4 = -0.5)
  d$z <- unname(cell[paste0(d$f1, d$ci)])
  d
}

test_that("zero residual variation is judged at the scale of the response", {
  zero <- rasch:::.residual_variation_is_zero
  y <- c(0.1, -0.1, 0.37, -0.37, 0.23, -0.23)
  expect_true(zero(9e-31, y))
  expect_true(zero(0, y))
  expect_false(zero(0.05, y))
  # the judgement follows the response's scale in both directions
  expect_true(zero(9e-31 * 1e-12, y * 1e-6))
  expect_false(zero(0.05 * 1e-12, y * 1e-6))
  expect_true(zero(9e-31 * 1e12, y * 1e6))
  expect_false(zero(0.05 * 1e12, y * 1e6))
  # weights enter the scale as they enter the residual
  expect_true(zero(1e-30, y, w = rep(0.5, 6)))
  # a response of zeros is fitted exactly by any model
  expect_true(zero(0, rep(0, 6)))
  # inputs the check cannot judge are left to the caller
  expect_false(zero(NA_real_, y))
  expect_false(zero(1e-30, c(y, NA)))
})

test_that("an exactly fitted response withholds every test with the reason", {
  d <- sim_exact_cells()
  terms <- c("f1", "ci", "f1:ci")
  for (v in c("classical", "hc3", "cr3")) {
    got <- rasch:::.dif_type2(d, terms, variance = v,
      robust_terms = if (v == "hc3") "f1" else NULL,
      cluster = if (v == "cr3") factor(d$pid) else NULL)
    tested <- got$term != "Residuals"
    expect_identical(got$term[tested], terms, info = v)
    expect_true(all(is.na(got$F_value[tested])), info = v)
    expect_true(all(is.na(got$p[tested])), info = v)
    expect_true(all(grepl("no residual variation to test against",
                          got$unavailable_reason[tested], fixed = TRUE)),
                info = v)
    # the sums of squares of the terms remain reportable; the residual is
    # reported as the zero it was judged to be
    expect_equal(got$sum_sq[got$term == "f1"], 3.6, info = v)
    expect_true(all(is.na(got$resid_ss[tested])), info = v)
    expect_identical(got$sum_sq[!tested], 0, info = v)
    expect_identical(got$mean_sq[!tested], 0, info = v)
    expect_equal(got$df[!tested], 32, info = v)
  }

  # residual variation within the cells restores every test
  set.seed(51)
  d$z <- d$z + rnorm(nrow(d), sd = 0.3)
  for (v in c("classical", "hc3", "cr3")) {
    got <- rasch:::.dif_type2(d, terms, variance = v,
      robust_terms = if (v == "hc3") "f1" else NULL,
      cluster = if (v == "cr3") factor(d$pid) else NULL)
    tested <- got$term != "Residuals"
    expect_true(all(is.finite(got$F_value[tested])), info = v)
    expect_true(all(is.finite(got$p[tested])), info = v)
    expect_true(all(is.na(got$unavailable_reason[tested])), info = v)
    expect_true(got$sum_sq[!tested] > 0, info = v)
  }
})

test_that("the within-person tests apply the same judgement", {
  # every person of a group has that group's pair of cell means, so the
  # between design fits the within contrast scores exactly
  f1 <- factor(rep(c("a", "b"), each = 12))
  y <- cbind(t1 = ifelse(f1 == "a", 0.2, -0.3),
             t2 = ifelse(f1 == "a", -0.2, 0.3))
  rownames(y) <- paste0("p", seq_along(f1))
  pd <- data.frame(pid = rownames(y), f1 = f1)
  got <- rasch:::.dif_within_tests(y, pd, "occ", list(occ = 2L),
                                   c("occ", "f1:occ"), "f1")
  tested <- got$term != "Residuals"
  expect_identical(got$term[tested], c("occ", "f1:occ"))
  expect_true(all(is.na(got$F_value[tested])))
  expect_true(all(grepl("no residual variation to test against",
                        got$unavailable_reason[tested], fixed = TRUE)))

  set.seed(52)
  y <- y + matrix(rnorm(length(y), sd = 0.3), nrow(y))
  got <- rasch:::.dif_within_tests(y, pd, "occ", list(occ = 2L),
                                   c("occ", "f1:occ"), "f1")
  tested <- got$term != "Residuals"
  expect_true(all(is.finite(got$F_value[tested])))
  expect_true(all(is.na(got$unavailable_reason[tested])))
})

# Item six is scored from the other five: a group-a person passes it with
# two or more of them right, a group-b person with four or more. With five
# class intervals, every person of a group in an interval then has the same
# raw score and the same response, so the same residual mean, and the
# factorial design fits the item's residual means exactly.
sim_two_cut_item <- function() {
  set.seed(1)
  n <- 300L
  delta <- seq(-1.5, 1.5, length.out = 5)
  theta <- rnorm(n)
  x <- matrix(rbinom(n * 5, 1, plogis(outer(theta, delta, "-"))), n, 5)
  g <- rep(c("a", "b"), length.out = n)
  s <- rowSums(x)
  x6 <- ifelse(g == "a", as.integer(s >= 2), as.integer(s >= 4))
  x <- cbind(x, x6)
  colnames(x) <- paste0("I", 1:6)
  rasch(data.frame(x, g = g), factors = "g")
}

test_that("dif_anova withholds the tests of an exactly fitted item", {
  skip_on_cran()
  fit <- sim_two_cut_item()
  da <- dif_anova(fit, n_groups = 5)
  six <- da$terms[da$terms$item == "I6", ]
  expect_setequal(six$term, c("g", "ci", "g:ci", "Residuals"))
  expect_true(all(is.na(six$F_value)))
  expect_true(all(is.na(six$p)))
  expect_true(all(is.na(six$p_adj)))
  expect_false(any(six$significant))
  expect_identical(six$sum_sq[six$term == "Residuals"], 0)
  expect_match(da$notes, "I6 \\[g\\]: unavailable because the retained design fits every residual mean exactly", all = FALSE)
  expect_match(da$notes, "I6 \\[g:ci\\]: unavailable because the retained design fits every residual mean exactly", all = FALSE)
  expect_match(da$notes, "I6 \\[ci\\]: unavailable because the retained design fits every residual mean exactly.*not a requested DIF test", all = FALSE)
  # the withheld tests stay in the family: the other items' adjustments
  # count them
  others <- da$terms[da$terms$item != "I6" & da$terms$term != "Residuals" &
                     da$terms$term != "ci", ]
  expect_true(all(is.finite(others$F_value)))
  expect_true(max(others$F_value) < 100)
  expect_equal(others$p_adj, pmin(p.adjust(others$p, "holm", n = 12), 1),
               tolerance = 1e-12)
  s6 <- da$summary[da$summary$item == "I6", ]
  expect_true(is.na(s6$F_uniform))
  expect_true(is.na(s6$F_nonuniform))
  expect_false(isTRUE(s6$uniform_DIF))
  expect_false(isTRUE(s6$nonuniform_DIF))
})
