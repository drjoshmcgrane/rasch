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
  expect_no_error(plot_pimap(f, information = TRUE))
  expect_no_error(plot_wright(f))
  expect_no_error(plot_wright(fp, bins = 20, xlim = c(-3, 3)))
  # a scale range excluding some persons and thresholds still renders
  expect_no_error(plot_wright(f, xlim = c(-1, 1)))
  expect_no_error(plot_pimap(f, xlim = c(-1, 1)))
  # class-interval and grid-range controls on the expected value curve
  expect_no_error(plot_icc(fp, "P03", n_groups = 8, grid = seq(-3, 3, 0.05)))
  # multi-item overlays use a common proportional score scale and may omit
  # the observed class-interval points
  expect_no_error(plot_icc(fp, c("P01", "P03", "P05"), observed = TRUE))
  expect_no_error(plot_icc(fp, c(1, 2), observed = FALSE))
  expect_error(plot_icc(fp, c("P01", "P01")), "same item more than once")
})

test_that("default model-curve grids follow the calibration origin", {
  fit <- rasch(simulate_rasch(n_persons = 300, n_items = 6, seed = 411))
  anchors <- data.frame(
    item = fit$items$item[1:3], k = NA_real_,
    tau = fit$items$location[1:3] + 100, average = TRUE)
  shifted <- rasch(fit$X, anchors = anchors)
  grid0 <- rasch:::.model_grid(fit)
  grid1 <- rasch:::.model_grid(shifted)
  expect_equal(grid1, grid0 + 100, tolerance = 1e-7)

  failed <- fit
  failed$est$converged <- FALSE
  expect_error(plot_icc(failed, failed$items$item[1]), "did not converge")
  expect_error(plot_tcc(failed), "did not converge")
  expect_error(plot_pimap(failed), "did not converge")
  expect_error(plot_wright(failed), "did not converge")
  expect_error(plot_threshold_map(failed), "did not converge")
  expect_error(plot_item_map(failed), "did not converge")
  expect_error(plot_person_fit(failed), "did not converge")
  expect_error(plot_kidmap(failed, 1), "did not converge")
  expect_error(plot_resid_dist(failed), "did not converge")

  failed_mfrm <- failed
  class(failed_mfrm) <- c("rasch_mfrm", "rasch")
  expect_error(plot_facets(failed_mfrm), "did not converge")
})

test_that("targeting plots require a converged EFRM link", {
  d <- simulate_efrm(n_per_group = 55, items_per_set = 4, n_sets = 2,
                     n_groups = 2, n_categories = 2, seed = 73)
  truth <- attr(d, "truth")
  fit <- rasch_efrm(d, item_sets = truth$item_sets, groups = "group",
                    id = "id", boot_reps = 0)
  fit$linking$alpha_edges$converged[1] <- FALSE

  expect_error(plot_pimap(fit), "set-unit link did not converge")
  expect_error(plot_wright(fit), "set-unit link did not converge")
  expect_error(plot_frames(fit), "set-unit link did not converge")
})

test_that("the default person-item scale labels beyond every estimate", {
  s <- rasch:::.pimap_scale(c(-4.2, 3.3))
  expect_lte(s$range[1], -4.2)
  expect_gte(s$range[2], 3.3)
  expect_equal(range(s$ticks), s$range)
  expect_gt(max(s$ticks), 3.3)

  fixed <- rasch:::.pimap_scale(c(-4.2, 3.3), c(-3.5, 2.5))
  expect_equal(fixed$range, c(-3.5, 2.5))
  expect_true(all(fixed$ticks >= -3.5 & fixed$ticks <= 2.5))
})

test_that("the kidmap and batch savers work, and Q3 pairs are complete", {
  set.seed(11)
  d <- seq(-2, 2, length.out = 10)
  X <- matrix(rbinom(200 * 10, 1, plogis(outer(rnorm(200), d, "-"))), 200, 10)
  colnames(X) <- sprintf("I%02d", 1:10)
  f <- rasch(X)

  pdf(NULL); on.exit(dev.off())
  expect_error(plot_icc(f, f$items$item[1:9]), "At most eight")
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
  expect_error(save_person_plots(f, zip_path, persons = 1.5),
               "whole row numbers")
  expect_error(save_person_plots(f, zip_path, persons = matrix(1L)),
               "ordinary vector")
  expect_error(save_item_plots(f, "icc", pdf_path, items = c(1L, 1L)),
               "same item")
  expect_error(save_item_plots(f, "icc", pdf_path, items = matrix(1L)),
               "ordinary vector")
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
  expect_error(plot_icc(f, c("I1", "I2"), group = s),
               "single item")
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
  expect_match(dimensionality_test(f, component = 1, min_score_points = 2)$split,
               "component 1")
  expect_match(dimensionality_test(f, component = 2, min_score_points = 2)$split,
               "component 2")
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
  mi <- test_information(mf, grid = c(-1, 0, 1))
  expect_true("design" %in% names(mi))
  expect_equal(length(unique(mi$design)), 1L)
  expect_match(unique(mi$design), "rater=A.*rater=B")
  expect_no_error(plot_tif(mf, grid = c(-1, 0, 1)))
  expect_no_error(plot_tcc(mf, grid = c(-1, 0, 1)))
  # A factor is still an ordinary categorical item selector; do not route it
  # to the virtual-cell lookup merely because it is not a character vector.
  expect_no_error(plot_icc(
    mf, factor(unique(as.character(mf$virtual_map$item))[1L])))

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
  expect_error(plot_pca_biplot(ef), "no respondents in common")
  expect_no_error(plot_resid_cor(ef))
  ei <- test_information(ef, grid = c(-1, 0, 1))
  expect_equal(length(unique(ei$design)), 2L)
  expect_equal(nrow(ei), 6L)
  line_types <- numeric(0)
  expect_no_error(testthat::with_mocked_bindings(
    plot_tif(ef, grid = c(-1, 0, 1)),
    lines = function(x, y = NULL, ..., lty = 1) {
      line_types <<- c(line_types, lty)
      graphics::lines(x, y, ..., lty = lty)
    },
    .package = "rasch"))
  expect_gte(sum(line_types == 5), 2L)
  expect_no_error(plot_tcc(ef, grid = c(-1, 0, 1)))
})

test_that("fit plots display standardised infit and outfit with the count note", {
  set.seed(4)
  d <- seq(-1, 1, length.out = 6)
  X <- matrix(rbinom(200 * 6, 1, plogis(outer(rnorm(200), d, "-"))), 200, 6)
  colnames(X) <- paste0("I", 1:6)
  f <- suppressWarnings(rasch(X))
  expect_true(all(c("infit_z", "outfit_z") %in% names(f$person)))
  expect_true(any(is.finite(f$person$infit_z)))
  png(tf <- tempfile(fileext = ".png"))
  on.exit({dev.off(); unlink(tf)}, add = TRUE)
  expect_no_error(plot_person_fit(f, statistic = "infit"))
  expect_no_error(plot_person_fit(f, statistic = "outfit"))
  expect_no_error(plot_item_map(f, statistic = "infit"))
  expect_no_error(plot_item_map(f, statistic = "outfit"))
  boundary <- f
  boundary$items$fit_resid[1] <- -Inf
  boundary$person$fit_resid[1] <- -Inf
  expect_no_error(plot_item_map(boundary))
  expect_no_error(plot_person_fit(boundary))
  expect_no_error(plot_resid_dist(boundary, what = "items"))
  expect_no_error(plot_resid_dist(boundary, what = "persons"))
  boundary$items$fit_resid[] <- -Inf
  boundary$person$fit_resid[] <- -Inf
  expect_error(plot_item_map(boundary), "does not carry")
  expect_error(plot_person_fit(boundary), "does not carry")
  expect_error(plot_resid_dist(boundary, what = "items"), "fewer than 3")
  expect_error(plot_resid_dist(boundary, what = "persons"), "fewer than 3")
})

test_that("the PCA biplot handles identically zero displayed loadings", {
  pc <- list(
    loadings_matrix = data.frame(
      item = c("I1", "I2"), PC1 = c(0, 0), PC2 = c(0, 0)),
    eigen_table = data.frame(proportion = c(0, 0)))
  pdf(NULL); on.exit(dev.off())
  expect_no_error(testthat::with_mocked_bindings(
    plot_pca_biplot(structure(list(), class = "rasch")),
    residual_pca = function(...) pc,
    .package = "rasch"))
})

test_that("the single-component PCA plot handles zero displayed loadings", {
  pc <- list(
    loadings_matrix = data.frame(item = c("I1", "I2"), PC1 = c(0, 0)),
    eigen_table = data.frame(eigenvalue = 0, proportion = 0))
  fit <- structure(list(items = data.frame(
    item = c("I1", "I2"), location = c(-0.5, 0.5))), class = "rasch")
  pdf(NULL); on.exit(dev.off())
  expect_no_error(testthat::with_mocked_bindings(
    plot_pca(fit),
    residual_pca = function(...) pc,
    .package = "rasch"))
})

test_that("a restricted person-item map carries its own information", {
  d <- simulate_rasch(300, 10, seed = 31)
  f <- rasch(d, id = "id")
  grid <- seq(-2, 2, by = 0.5)
  full <- test_information(f, grid)
  sub <- test_information(f, grid, items = c("I01", "I02", "I03"))
  expect_true(all(sub$info < full$info))
  # the subset curve is exactly the sum of its own item variances
  by_hand <- vapply(grid, function(t)
    sum(vapply(1:3, function(i) item_moments(t, f$tau_list[[i]])$V, 0)), 0)
  expect_equal(sub$info, by_hand, tolerance = 1e-10)
  expect_error(test_information(f, grid, items = "NOPE"), "not in the fit")
})

test_that("an ambiguous group level is refused, and the qualified form works", {
  d <- simulate_rasch(200, 6, seed = 32)
  d$grp <- rep(c("a", "b"), 100)
  d$site <- rep(c("a", "z"), each = 100)   # level "a" in both factors
  f <- rasch(d, id = "id", factors = c("grp", "site"))
  expect_error(.pimap_persons(f, "a"), "several factors")
  expect_equal(sum(.pimap_persons(f, "grp: a")), 100)
  expect_equal(sum(.pimap_persons(f, "site: a")), 100)
  expect_false(identical(.pimap_persons(f, "grp: a"),
                         .pimap_persons(f, "site: a")))
  expect_error(.pimap_persons(f, "grp: q"), "not a level")
  # an unambiguous bare level still works
  expect_equal(sum(.pimap_persons(f, "z")), 100)
})

test_that("test_information refuses fractional indices and finds EFRM items", {
  d <- simulate_rasch(200, 6, seed = 33)
  f <- rasch(d, id = "id")
  expect_error(test_information(f, items = 1.9), "whole numbers")
  de <- simulate_efrm(n_per_group = 80, items_per_set = 5, n_sets = 2,
                      n_groups = 2, seed = 34)
  fe <- rasch_efrm(de, item_sets = attr(de, "truth")$item_sets,
                   groups = "group", boot_reps = 0)
  base <- as.character(fe$virtual_map$item[1])
  ti <- test_information(fe, seq(-1, 1, by = 0.5), items = base)
  expect_true(nrow(ti) > 0)
  # an EFRM map restricted to one group carries only that group's designs
  grid <- seq(-1, 1, by = 0.5)
  g1 <- levels(factor(fe$factors[[fe$frame_group[1]]]))[1]
  gcols <- which(as.character(fe$virtual_map$group) == g1)
  ti_g <- test_information(fe, grid, items = gcols)
  expect_true(all(grepl(paste0("group=", g1), unique(ti_g$design), fixed = TRUE)))
})

test_that("qualified group names split on the factor name, not any colon", {
  d <- simulate_rasch(200, 6, seed = 35)
  d$group <- rep(c("a", "b"), 100)
  d$cohort <- rep(c("2024: spring", "2024: autumn"), 100)   # colons in levels
  f <- rasch(d, id = "id", factors = c("group", "cohort"))
  k <- .pimap_persons(f, "cohort: 2024: spring")
  expect_equal(sum(k), 100)
  expect_error(.pimap_persons(f, "cohort: 2024: winter"), "not a level")
  # a bare qualified-looking string that is actually a plain level
  expect_error(.pimap_persons(f, "nofactor: x"), "does not match")
})

test_that("a crossed-frame EFRM map keeps only the selected group designs", {
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  de <- simulate_efrm(n_per_group = 80, items_per_set = 5, n_sets = 2,
                      n_groups = 2, seed = 36)
  tr <- attr(de, "truth")
  de$grp <- tr$groups
  de$region <- rep(c("N", "S"), length.out = nrow(de))
  fe <- rasch_efrm(de, item_sets = tr$item_sets,
                   groups = c("grp", "region"), boot_reps = 0)
  keep_p <- .pimap_persons(fe, "grp: g1")
  frames_of <- as.character(fe$factors[[fe$frame_group[1]]])
  glev <- unique(frames_of[keep_p])
  gcols <- which(as.character(fe$virtual_map$group) %in% glev)
  ti <- test_information(fe, seq(-1, 1, by = 0.5), items = gcols)
  expect_setequal(glev, c("g1:N", "g1:S"))
  expect_true(all(startsWith(unique(ti$design), "group=g1:")))
  expect_equal(length(unique(ti$design)), 2L)
  expect_equal(length(unique(test_information(
    fe, seq(-1, 1, by = 0.5))$design)), 4L)
  # The crossed factor and its level both contain colons; the qualified
  # selector still addresses that one cell exactly.
  cell <- sprintf("%s: %s", fe$frame_group[1], glev[1])
  expect_equal(sum(.pimap_persons(fe, cell)), sum(frames_of == glev[1]))

  exact_g1 <- as.character(fe$virtual_map$vkey[
    which(as.character(fe$virtual_map$group) == "g1:N")[1L]])
  expect_error(
    plot_pimap(fe, information = TRUE, group = "grp: g2", items = exact_g1),
    "select no common EFRM response cells")
  expect_no_error(plot_pimap(
    fe, information = TRUE, group = "grp: g2",
    items = as.character(fe$virtual_map$item[1L])))
})

test_that("category-frequency plots validate the fitted model first", {
  expect_error(
    plot_catfreq(1, "I1"),
    "fitted response-data Rasch model"
  )

  set.seed(117)
  X <- matrix(
    rbinom(300 * 5, 1, 0.5), 300, 5,
    dimnames = list(NULL, paste0("I", 1:5))
  )
  failed <- rasch(X)
  failed$est$converged <- FALSE
  expect_error(plot_catfreq(failed, "I1"), "did not converge")
})

test_that("ICC overlays retain fitted intervals and honour observed = FALSE", {
  set.seed(118)
  X <- matrix(
    rbinom(360 * 6, 1, 0.5), 360, 6,
    dimnames = list(NULL, paste0("I", 1:6))
  )
  # Different missingness makes the fitted allocation item-specific.
  X[seq(1, 360, by = 3), "I1"] <- NA
  X[seq(2, 360, by = 4), "I2"] <- NA
  fit <- rasch(X)
  used <- integer(0)
  original <- rasch:::.fit_class_intervals
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  testthat::with_mocked_bindings(
    plot_icc(fit, c("I1", "I2"), observed = TRUE),
    .fit_class_intervals = function(fit, i, ...) {
      used <<- c(used, i)
      original(fit, i, ...)
    },
    .package = "rasch"
  )
  expect_identical(used, c(1L, 2L))

  n_points <- 0L
  testthat::with_mocked_bindings(
    plot_icc(fit, "I1", group = rep(c("A", "B"), each = 180),
             observed = FALSE),
    points = function(...) n_points <<- n_points + 1L,
    .package = "rasch"
  )
  expect_identical(n_points, 0L)

  legend_cols <- NULL
  many_groups <- rep(paste0("g", 1:9), each = nrow(fit$X) / 9)
  testthat::with_mocked_bindings(
    plot_icc(fit, "I1", group = many_groups, observed = TRUE,
             n_groups = 2),
    .rr_legend = function(pos, ...) legend_cols <<- list(...)$col,
    .package = "rasch"
  )
  expect_length(legend_cols, 9L)
  expect_false(anyNA(legend_cols))
  expect_identical(legend_cols[9L], legend_cols[1L])
})

test_that("threshold and PCA plot selectors do not truncate their displays", {
  fit <- structure(list(
    est = list(converged = TRUE), model = "PCM",
    tau_list = list(seq(-4.5, 4.5, length.out = 9)),
    items = data.frame(item = "P1", location = 0),
    n_groups = 5L
  ), class = "rasch")
  legend_col <- NULL
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  testthat::with_mocked_bindings(
    plot_threshold_prob(fit, "P1"),
    .rr_legend = function(pos, ...) legend_col <<- list(...)$col,
    .package = "rasch"
  )
  expect_length(legend_col, 9L)
  expect_false(anyNA(legend_col))
  expect_identical(legend_col[9L], legend_col[1L])

  requested <- NULL
  pc <- list(
    loadings_matrix = data.frame(
      item = paste0("I", 1:12), PC12 = seq(-0.5, 0.5, length.out = 12)),
    eigen_table = data.frame(eigenvalue = rep(1, 12),
                             proportion = rep(1 / 12, 12))
  )
  pfit <- structure(list(items = data.frame(
    item = paste0("I", 1:12), location = seq(-1, 1, length.out = 12))),
    class = "rasch")
  testthat::with_mocked_bindings(
    plot_pca(pfit, component = 12),
    residual_pca = function(fit, n_components) {
      requested <<- n_components
      pc
    },
    .package = "rasch"
  )
  expect_identical(requested, 12L)
})

test_that("kidmap does not turn missing uncertainty into a zero-width band", {
  fit <- rasch(simulate_rasch(250, 8, seed = 119), id = "id")
  n <- which(is.finite(fit$person$theta))[1L]
  fit$person$se[n] <- NA_real_
  expect_error(plot_kidmap(fit, n), "standard error.*unavailable")
})

test_that("person-item map item selections cannot overstate their size", {
  fit <- rasch(simulate_rasch(250, 8, seed = 120), id = "id")
  expect_error(.pimap_items(fit, c("I01", "I01")), "same item")
  expect_error(.pimap_items(fit, " "), "must name")
})

test_that("test-level displays say so when designs outrun the palette", {
  # A linking design can carry far more administrations than the palette can
  # name. Cycling colours there would imply distinctions no reader can follow
  # back to a curve, so the legend counts them instead. The statistics are
  # untouched: every design still has its own curve.
  set.seed(21); Np <- 160L; L <- 5L
  d <- seq(-1.2, 1.2, length.out = L)
  raters <- c("A", "B", "C", "D")
  th <- rnorm(Np)
  blocks <- lapply(seq_along(raters), function(k) {
    X <- matrix(rbinom(Np * L, 1, plogis(outer(th, d, "-"))), Np, L)
    colnames(X) <- sprintf("I%02d", seq_len(L))
    X
  })
  dd <- data.frame(person = rep(seq_len(Np), length(raters)),
                   do.call(rbind, blocks),
                   rater = rep(raters, each = Np), check.names = FALSE)
  # each person is rated by a random subset of at least two of the raters
  keep <- runif(nrow(dd)) >= 0.4
  for (p in seq_len(Np)) {
    idx <- which(dd$person == p)
    if (sum(keep[idx]) < 2L) keep[idx[1:2]] <- TRUE
  }
  mf <- rasch_mfrm(dd[keep, ], person = "person", facets = "rater",
                   items = sprintf("I%02d", seq_len(L)))
  n_des <- length(unique(test_information(mf, grid = c(-1, 0, 1))$design))
  expect_gt(n_des, length(.rr$pal))

  legends <- list()
  spy <- function(...) {
    legends[[length(legends) + 1L]] <<- list(...)
    invisible(NULL)
  }
  grDevices::png(tf <- tempfile(fileext = ".png"))
  on.exit({grDevices::dev.off(); unlink(tf)}, add = TRUE)
  testthat::with_mocked_bindings(
    {
      plot_tif(mf, grid = seq(-2, 2, by = 0.5))
      plot_tcc(mf, grid = seq(-2, 2, by = 0.5))
      plot_pimap(mf, information = TRUE)
    },
    .rr_legend = spy, .package = "rasch")
  counted <- vapply(legends, function(a)
    length(a[[2L]]) == 1L && grepl("more than the palette can name", a[[2L]]),
    NA)
  expect_equal(sum(counted), 3L)
  for (a in legends[counted]) {
    expect_match(a[[2L]], sprintf("^%d administrations;", n_des))
    expect_true(all(a$col == .rr$soft))
  }

  # below the limit each design is still named in its own colour
  de <- simulate_efrm(n_per_group = 80, items_per_set = 5, n_sets = 1,
                      n_groups = 2, seed = 36)
  tr <- attr(de, "truth")
  de$grp <- tr$groups
  ef <- rasch_efrm(de, item_sets = tr$item_sets, groups = "grp",
                   boot_reps = 0)
  legends <- list()
  testthat::with_mocked_bindings(
    {
      plot_tif(ef, grid = c(-1, 0, 1))
      plot_tcc(ef, grid = c(-1, 0, 1))
      plot_pimap(ef, information = TRUE)
    },
    .rr_legend = spy, .package = "rasch")
  named <- vapply(legends, function(a)
    length(a[[2L]]) == 2L && all(grepl("^group=", a[[2L]])), NA)
  expect_equal(sum(named), 3L)
  for (a in legends[named]) expect_equal(length(unique(a$col)), 2L)
})
