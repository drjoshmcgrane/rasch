# dif_anova() on a split fit: class intervals formed on a statistic common
# to the groups, so a split does not separate the groups into different
# intervals and drain the power of the tests on the remaining items.

interval_sim <- function(seed, n = 800, L = 10, shift = NULL) {
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

test_that("a split fit forms class intervals on the first copy's mapping", {
  s <- interval_sim(12, shift = list("2" = 0.9, "8" = -0.9))
  # an unsplit fit keeps its own intervals and says nothing about them
  it0 <- rasch:::.dif_interval_theta(s$fit)
  expect_false(it0$common)
  expect_identical(it0$theta, s$fit$person$theta)
  da0 <- dif_anova(s$fit, "grp")
  expect_false(any(grepl("class intervals formed on", da0$notes)))
  f2 <- split_items(s$fit, "I2", by = "grp")
  it <- rasch:::.dif_interval_theta(f2)
  expect_true(it$common)
  # the first copy is group a's, so group a keeps its own locations and
  # group b is mapped through the same copy
  expect_equal(it$theta[s$g == "a"], f2$person$theta[s$g == "a"])
  expect_false(isTRUE(all.equal(it$theta[s$g == "b"],
                                f2$person$theta[s$g == "b"])))
  # with complete data the intervals are a partition of the merged raw score
  ci <- rasch:::.dif_class_intervals(f2, f2$n_groups)
  rs <- rowSums(s$X)
  shared <- tapply(ci, rs, function(z) length(unique(z)))
  expect_true(all(shared == 1L))
  expect_true(all(diff(tapply(rs, ci, max)) > 0))
  # every interval holds both groups
  expect_true(all(table(ci, s$g) > 0))
  da <- dif_anova(f2, "grp")
  expect_match(paste(da$notes, collapse = " "),
               "class intervals formed on the locations under the first copy of each split item \\(I2\\)")
  r8 <- da$summary[da$summary$item == "I8", ]
  expect_true(r8$uniform_DIF)
  expect_true(r8$p_uniform_adj < 0.01)
  # the remaining planted item is now resolved as well
  r <- resolve_dif(s$fit)
  expect_setequal(r$splits$item, c("I2", "I8"))
  expect_equal(r$n_remaining_dif, 0L)
})
