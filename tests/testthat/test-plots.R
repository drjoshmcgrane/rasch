# Targeting displays: the person-item threshold distribution and the
# conventional Wright map (Wright & Stone 1979), including the class-interval
# and scale-range controls the app exposes.

test_that("targeting plots render for dichotomous and polytomous fits", {
  set.seed(4)
  d <- seq(-2, 2, length.out = 8)
  X <- matrix(rbinom(300 * 8, 1, plogis(outer(rnorm(300), d, "-"))), 300, 8)
  colnames(X) <- sprintf("I%02d", 1:8)
  f <- rasch(X)
  simP <- function(th, t) {
    x <- 0:length(t); p <- exp(x * th - c(0, cumsum(t))); p / sum(p)
  }
  th <- rnorm(250)
  Xp <- sapply(1:5, function(i) vapply(th, function(t)
    sample(0:3, 1, prob = simP(t, c(-1, 0, 1))), 0L))
  colnames(Xp) <- sprintf("P%02d", 1:5)
  fp <- rasch(Xp)

  pdf(NULL); on.exit(dev.off())
  expect_no_error(plot_pimap(f))
  expect_no_error(plot_pimap(fp, bins = 15, xlim = c(-2, 2)))
  expect_no_error(plot_wright(f))
  expect_no_error(plot_wright(fp, bins = 20, xlim = c(-3, 3)))
  # a scale range excluding some persons and thresholds still renders
  expect_no_error(plot_wright(f, xlim = c(-1, 1)))
  expect_no_error(plot_pimap(f, xlim = c(-1, 1)))
  # class-interval and grid-range controls on the expected value curve
  expect_no_error(plot_icc(fp, "P03", n_groups = 8, grid = seq(-3, 3, 0.05)))
})
