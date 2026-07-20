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

test_that("the kidmap and batch savers work, and Q3 pairs are complete", {
  set.seed(11)
  d <- seq(-2, 2, length.out = 10)
  X <- matrix(rbinom(200 * 10, 1, plogis(outer(rnorm(200), d, "-"))), 200, 10)
  colnames(X) <- sprintf("I%02d", 1:10)
  f <- rasch(X)

  pdf(NULL); on.exit(dev.off())
  expect_no_error(plot_kidmap(f, person = 1))
  expect_no_error(plot_kidmap(f, person = 2, level = 0.9, xlim = c(-3, 3)))
  expect_error(plot_kidmap(f, person = "no-such-id"), "not found")

  # Yen's Q3: every off-diagonal pair, star = excess over the average
  rc <- residual_correlations(f)
  expect_equal(nrow(rc$pairs), choose(10, 2))
  expect_equal(rc$pairs$q3_star, rc$pairs$q3 - rc$average)
  expect_true(all(rc$flagged$q3_star > 0.2))
  expect_identical(names(rc$pairs),
                   c("item_a", "item_b", "q3", "q3_star", "flagged"))

  # batch savers: multi-page PDF and ZIP-of-PNGs by extension
  pdf_path <- file.path(tempdir(), "icc_all.pdf")
  zip_path <- file.path(tempdir(), "kidmaps.zip")
  expect_equal(save_item_plots(f, "icc", pdf_path), pdf_path)
  expect_true(file.exists(pdf_path))
  save_person_plots(f, zip_path, persons = 1:4)
  expect_equal(length(utils::unzip(zip_path, list = TRUE)$Name), 4L)
  expect_error(save_item_plots(f, "icc", "bad.txt"), "pdf or")
  unlink(c(pdf_path, zip_path))
})

test_that("the DIF overlay accepts one or several nominated factor names", {
  set.seed(2); n <- 300
  g <- rep(c("ref", "foc"), n / 2); s <- sample(c("m", "f"), n, TRUE)
  d <- seq(-1, 1, length.out = 5)
  X <- matrix(rbinom(n * 5, 1, plogis(outer(rnorm(n), d, "-"))), n, 5)
  colnames(X) <- paste0("I", 1:5)
  f <- rasch(data.frame(X, group = g, sex = s), factors = c("group", "sex"))
  pdf(NULL); on.exit(dev.off())
  expect_no_error(plot_icc(f, "I2", group = "sex"))
  # several names draw the factor-combination cells (the factorial display)
  expect_no_error(plot_icc(f, "I2", group = c("group", "sex")))
  # a raw person vector still works
  expect_no_error(plot_icc(f, "I2", group = s))
})

test_that("residual components beyond the first can be inspected and tested", {
  set.seed(4)
  d <- seq(-2, 2, length.out = 10)
  X <- matrix(rbinom(500 * 10, 1, plogis(outer(rnorm(500), d, "-"))), 500, 10)
  colnames(X) <- sprintf("I%02d", 1:10)
  f <- rasch(X)
  pdf(NULL); on.exit(dev.off())
  expect_no_error(plot_pca(f, component = 1))
  expect_no_error(plot_pca(f, component = 2))
  expect_no_error(plot_pca(f, component = 3))
  expect_error(plot_pca(f, component = 99), "not available")
  # the t-test default split follows the chosen component
  expect_match(dimensionality_test(f, component = 1)$split, "component 1")
  expect_match(dimensionality_test(f, component = 2)$split, "component 2")
})

test_that("residual dependence displays generalise to MFRM and EFRM fits", {
  # MFRM and EFRM inherit from "rasch" and carry residuals over their virtual
  # items, so the residual-PCA / Q3 / biplot suite must run on them unchanged
  set.seed(11); Np <- 300; L <- 6
  d <- seq(-1.5, 1.5, length.out = L)
  X <- matrix(rbinom(Np * L, 1, plogis(outer(rnorm(Np), d, "-"))), Np, L)
  colnames(X) <- sprintf("I%02d", 1:L)
  # each person appears under BOTH raters: a nested one-rater-per-person
  # design leaves severity confounded with the person blocks and is now
  # (correctly) refused by the connectivity check
  X2 <- matrix(rbinom(Np * L, 1, plogis(outer(rnorm(Np), d, "-"))), Np, L)
  colnames(X2) <- colnames(X)
  mf <- rasch_mfrm(data.frame(person = rep(seq_len(Np), 2), rbind(X, X2),
                              rater = rep(c("A", "B"), each = Np),
                              check.names = FALSE),
                   person = "person", facets = "rater", items = colnames(X))
  expect_false(is.null(mf$residuals))
  expect_s3_class(residual_pca(mf)$loadings_matrix, "data.frame")
  expect_true("star_matrix" %in% names(residual_correlations(mf)))
  pdf(NULL); on.exit(dev.off())
  expect_no_error(plot_pca_biplot(mf))
  expect_no_error(plot_resid_cor(mf, stat = "q3"))
  expect_no_error(plot_resid_cor(mf, stat = "q3star"))

  # EFRM: one item set, two groups differing in discrimination so the sets link
  set.seed(12); per_g <- 300; glev <- c("G1", "G2")
  phi <- c(0.7, 1.3); grp <- rep(glev, each = per_g); Np2 <- length(grp)
  th <- rnorm(Np2, 0, 1.3); dd <- scale(seq(-2, 2, length.out = 10), scale = FALSE)[, 1]
  XE <- sapply(1:10, function(i)
    rbinom(Np2, 1, plogis(phi[match(grp, glev)] * (th - dd[i]))))
  colnames(XE) <- sprintf("E%02d", 1:10)
  ef <- rasch_efrm(data.frame(XE, g = grp),
                   item_sets = list(core = colnames(XE)), groups = "g")
  expect_false(is.null(ef$residuals))
  expect_no_error(plot_pca_biplot(ef))
  expect_no_error(plot_resid_cor(ef))
})
