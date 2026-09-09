# Output tables: chi-square class-interval detail, score-table
# estimators and extreme-score extrapolation, the LR test block, and the
# summary distribution statistics.

simd2 <- function(N, d, seed = 1) {
  set.seed(seed)
  th <- rnorm(N)
  X <- matrix(rbinom(N * length(d), 1, plogis(outer(th, d, "-"))),
              N, length(d))
  colnames(X) <- paste0("I", seq_along(d))
  X
}

test_that("chisq_detail reproduces the item-trait chi-square from its components", {
  fit <- rasch(simd2(400, seq(-1.5, 1.5, length.out = 8), seed = 2))
  cd <- chisq_detail(fit, "I4")
  expect_equal(sum(cd$intervals$chisq[cd$intervals$used]),
               fit$item_trait$chisq[4], tolerance = 1e-10)
  expect_equal(cd$df, fit$item_trait$df[4])
  # ES is the residual with the interval-size factor removed
  used <- cd$intervals$used
  expect_equal(cd$intervals$residual[used] / sqrt(cd$intervals$n[used]),
               cd$intervals$es[used], tolerance = 1e-10)
  # per-interval category proportions sum to 1; dichotomous OBS.T = OBS.P(1)
  cats <- cd$categories
  for (g in unique(cats$interval)) {
    rows <- cats$interval == g
    if (any(is.na(cats$obs_p[rows]))) next
    expect_equal(sum(cats$obs_p[rows]), 1, tolerance = 1e-10)
    expect_equal(sum(cats$est_p[rows]), 1, tolerance = 1e-8)
    expect_equal(cats$obs_t[rows][2], cats$obs_p[rows][2], tolerance = 1e-10)
  }
  # interval means increase along the trait
  expect_true(all(diff(na.omit(cd$intervals$theta_mean)) > 0))
})

test_that("score_table extrapolates extremes by the geometric rule", {
  fit <- rasch(simd2(500, seq(-2, 2, length.out = 10), seed = 6))
  # MLE table: infinite extremes are NA, interior solves the score equation
  mle <- score_table(fit, method = "mle")
  expect_true(is.na(mle$theta[1]) && is.na(mle$theta[11]))
  r5 <- mle$theta[mle$score == 5]
  expect_equal(sum(vapply(fit$tau_list, function(tt)
    item_moments(r5, tt)$E, 0)), 5, tolerance = 1e-6)
  # WLE is finite everywhere and flatter than MLE (shrunk towards centre)
  wle <- score_table(fit)
  expect_true(all(is.finite(wle$theta)))
  expect_lt(wle$theta[10], mle$theta[10])
  # geometric extrapolation: d_next = b^2 / a on the preceding differences
  ex <- score_table(fit, method = "mle", extremes = "extrapolated")
  th <- mle$theta
  b <- th[10] - th[9]; a <- th[9] - th[8]
  expect_equal(ex$theta[11], th[10] + b^2 / a, tolerance = 1e-10)
  b0 <- th[3] - th[2]; a0 <- th[4] - th[3]
  expect_equal(ex$theta[1], th[2] - b0^2 / a0, tolerance = 1e-10)
  expect_true(all(ex$extrapolated[c(1, 11)]))
  # extrapolated SE exceeds the SE one score in
  expect_gt(ex$se[11], ex$se[10])
  # frequencies count complete responders
  expect_equal(sum(ex$freq), sum(stats::complete.cases(fit$X)))
  expect_equal(ex$cum_pct[11], 100)

  failed <- fit
  failed$est$converged <- FALSE
  expect_error(score_table(failed), "did not converge")
  expect_error(person_extrapolated(failed), "did not converge")
  expect_error(chisq_detail(failed, "I1"), "did not converge")
})

test_that("person score roots follow an externally translated origin", {
  tau <- list(c(-1, 0, 1), c(-0.5, 0.5), c(0, 1))
  base <- person_wle(tau)
  shifted <- person_wle(lapply(tau, `+`, 100))
  expect_true(all(is.finite(shifted$theta)))
  expect_equal(unname(shifted$theta), unname(base$theta) + 100,
               tolerance = 1e-8)
  expect_equal(unname(shifted$se), unname(base$se), tolerance = 1e-8)

  fit <- rasch(simd2(500, seq(-2, 2, length.out = 8), seed = 61))
  anchors <- data.frame(
    item = fit$items$item[1:3], k = NA_real_,
    tau = fit$items$location[1:3] + 100, average = TRUE)
  translated <- rasch(fit$X, anchors = anchors)
  expect_identical(
    lapply(translated$threshold_diagnostics, `[[`,
           "never_modal_categories"),
    lapply(fit$threshold_diagnostics, `[[`,
           "never_modal_categories"))
  mle_base <- score_table(fit, method = "mle")
  mle_translated <- score_table(translated, method = "mle")
  interior <- is.finite(mle_base$theta)
  expect_equal(mle_translated$theta[interior],
               mle_base$theta[interior] + 100, tolerance = 1e-7)
  expect_equal(mle_translated$se[interior], mle_base$se[interior],
               tolerance = 1e-7)

  info_base <- test_information(fit)
  info_translated <- test_information(translated)
  expect_equal(info_translated$theta, info_base$theta + 100,
               tolerance = 1e-7)
  expect_equal(info_translated$info, info_base$info, tolerance = 1e-7)

  weights <- stats::setNames(rep(1, nrow(fit$items)), fit$items$item)
  weighted_base <- weighted_person_estimates(fit, weights)
  weighted_translated <- weighted_person_estimates(translated, weights)
  expect_true(all(is.finite(weighted_translated$theta)))
  expect_equal(weighted_translated$theta, weighted_base$theta + 100,
               tolerance = 1e-7)
  expect_equal(weighted_translated$se, weighted_base$se, tolerance = 1e-7)

  unequal_score <- getFromNamespace(".efrm_person_estimates", "rasch")
  X <- fit$X[, 1:3, drop = FALSE]
  units <- c(0.7, 1.1, 1.5)
  e0 <- unequal_score(X, fit$tau_list[1:3], units)
  e1 <- unequal_score(X, lapply(fit$tau_list[1:3], `+`, 100), units)
  expect_true(all(is.finite(e1$theta)))
  expect_equal(e1$theta, e0$theta + 100, tolerance = 1e-7)
  expect_equal(e1$se, e0$se, tolerance = 1e-7)
})

test_that("extreme-score extrapolation refuses collapsed interior spacings", {
  geometric <- getFromNamespace(".geometric_score_extremes", "rasch")
  expect_equal(unname(geometric(c(NA, 1, 2, 4, NA))), c(0.5, 8))
  expect_null(geometric(c(NA, 1, 1, 2, NA)))
  expect_null(geometric(c(NA, 1, 2, 2, NA)))
  expect_null(geometric(c(NA, 1, NA, 2, NA)))

  fit <- rasch(simd2(500, seq(-2, 2, length.out = 10), seed = 61))
  fit$score_table$theta[3] <- fit$score_table$theta[2]
  expect_error(score_table(fit, extremes = "extrapolated"),
               "stable increasing spacings")
})

test_that("modal-category diagnostics do not depend on a theta grid", {
  modal <- getFromNamespace(".modal_score_categories", "rasch")
  expect_identical(modal(c(-1, 0, 1)), 0:3)
  expect_identical(modal(c(1, -1)), c(0L, 2L))
  # Category 2 is modal only in the one-millionth-logit interval between the
  # two central thresholds. A regular plotting grid does not resolve it.
  expect_true(2L %in% modal(c(-1, 0.123, 0.123001, 2)))
  # Equal thresholds leave only the extreme categories modal. Their crossing
  # points differ by an ulp once cumulated, which must not become an interval.
  expect_identical(modal(c(0.1, 0.1, 0.1)), c(0L, 3L))
  tied <- 0.40785616827135612
  expect_identical(modal(c(tied, tied + 5e-17, tied)), c(0L, 3L))
  expect_identical(modal(0.5), 0:1)
  expect_identical(modal(c(2, 1, 0)), c(0L, 3L))
})

test_that("lr_test prefers PCM when thresholds differ and RSM when they do not", {
  set.seed(9)
  simpoly <- function(taus) {
    X <- sapply(taus, function(tt) vapply(rnorm(400), function(b)
      sample(0:2, 1, prob = item_moments(b, tt)$P), 0L))
    colnames(X) <- paste0("Q", seq_along(taus)); X
  }
  d <- seq(-1, 1, length.out = 6)
  # common threshold structure: the adjusted test retains RSM; the raw
  # composite statistic is anticonservative (that is why it is adjusted)
  same <- lapply(d, function(dd) c(-0.8, 0.8) + dd)
  lr1 <- lr_test(rasch(simpoly(same)))
  expect_gt(lr1$p_adj, 0.01)
  expect_equal(lr1$df, 5)  # 6 items x 2 thresholds - 1 vs (6 - 1) + 1
  expect_lt(lr1$chisq_adj, lr1$chisq)
  # the Godambe eigenvalues carry the inflation of the composite statistic
  expect_gt(mean(lr1$lambda), 1)
  # item-specific spreads: PCM needed, decisively
  diff_ <- lapply(seq_along(d), function(i) c(-0.2, 0.2) * i + d[i])
  lr2 <- lr_test(rasch(simpoly(diff_)))
  expect_lt(lr2$p_adj, 0.001)
  expect_gt(lr2$chisq, 0)
  # guards
  expect_error(lr_test(rasch(simd2(150, c(-1, 0, 1)))), "dichotomous")
})

test_that("summary blocks carry the distribution statistics", {
  X <- simd2(400, seq(-1.5, 1.5, length.out = 8), seed = 13)
  fit <- rasch(X)
  for (blk in list(fit$item_fit_summary, fit$person_fit_summary,
                   fit$summary_stats$item_location,
                   fit$summary_stats$person_location_noext)) {
    expect_true(all(c("mean", "sd", "skewness", "kurtosis") %in% names(blk)))
    expect_true(is.finite(blk$sd))
  }
  expect_equal(fit$summary_stats$item_location$mean, 0, tolerance = 1e-8)
  expect_true(is.finite(fit$summary_stats$cor_item_fit_location))
  expect_true(fit$summary_stats$df_factor < 1 && fit$summary_stats$df_factor > 0.8)
  # item separation index behaves like a reliability
  expect_true(fit$isi$PSI >= 0 && fit$isi$PSI <= 1)
  # spread-out items with a large sample are well separated
  expect_gt(fit$isi$PSI, 0.9)
  # alpha applicability flag
  expect_true(fit$alpha$applicable)
  Xm <- X; Xm[1, 1] <- NA
  expect_false(rasch(Xm)$alpha$applicable)
})

test_that("report_html writes a complete self-contained report", {
  set.seed(3)
  d <- seq(-1.5, 1.5, length.out = 6)
  X <- matrix(rbinom(250 * 6, 1, plogis(outer(rnorm(250), d, "-"))), 250, 6)
  colnames(X) <- paste0("I", 1:6)
  fit <- rasch(data.frame(X, g = rep(c("a", "b"), each = 125)), factors = "g")
  out <- file.path(tempdir(), "rasch_report_test.html")
  on.exit(unlink(out), add = TRUE)
  report_html(fit, out, title = "Test report")
  expect_true(file.exists(out))
  html <- paste(readLines(out, warn = FALSE), collapse = "\n")
  # self-contained: plots embedded, no external references
  expect_gt(lengths(regmatches(html, gregexpr("data:image/png;base64", html))), 5)
  expect_false(grepl("src=\"http", html))
  for (sec in c("Summary", "Item statistics", "Thresholds", "Score to measure",
                "Dimensionality", "Local dependence", "Classical companions",
                "Differential item functioning", "Person estimates"))
    expect_true(grepl(paste0("<h2>", sec), html), label = sec)
  # the base-R base64 encoder matches the RFC 4648 test vectors
  enc <- function(txt) {
    p <- tempfile(); on.exit(unlink(p), add = TRUE)
    writeBin(charToRaw(txt), p)
    .b64(p)
  }
  expect_equal(enc("Man"), "TWFu")
  expect_equal(enc("Ma"), "TWE=")
  expect_equal(enc("M"), "TQ==")
  expect_equal(enc("foobar"), "Zm9vYmFy")
})

test_that("export entry points refuse invalid fits before writing", {
  out <- tempfile("invalid-export-")
  expect_error(save_outputs(1, out, formats = "png", item_plots = FALSE),
               "fitted Rasch or paired-comparison model")
  expect_false(dir.exists(out))
  expect_error(report_html(1, tempfile(fileext = ".html")),
               "person-by-item Rasch model")
  expect_error(save_item_plots(1, "icc", tempfile(fileext = ".pdf")),
               "person-by-item Rasch model")
  expect_error(save_person_plots(1, tempfile(fileext = ".pdf")),
               "person-by-item Rasch model")

  fit <- rasch(simd2(120, seq(-1, 1, length.out = 5), seed = 231))
  fit$est$converged <- FALSE
  out_failed <- tempfile("failed-export-")
  html_failed <- tempfile(fileext = ".html")
  expect_error(save_outputs(fit, out_failed, formats = "png",
                            item_plots = FALSE), "did not converge")
  expect_false(dir.exists(out_failed))
  expect_error(report_html(fit, html_failed), "did not converge")
  expect_false(file.exists(html_failed))
  expect_error(report_document(fit, html_failed), "did not converge")

  bd <- simulate_btl(5, 10, reps_per_pair = 4, seed = 232)
  bt <- btl(bd, "object_a", "object_b", winner = "winner", judge = "judge")
  bt$converged <- FALSE
  out_btl <- tempfile("failed-btl-export-")
  expect_error(save_outputs(bt, out_btl, formats = "png",
                            item_plots = FALSE), "did not converge")
  expect_false(dir.exists(out_btl))
  expect_error(report_document(bt, html_failed), "did not converge")
})

test_that("exports accept only a compatible fit-bootstrap result", {
  fit <- rasch(simd2(160, seq(-1, 1, length.out = 5), seed = 24))
  other <- rasch(simd2(160, seq(-1, 1, length.out = 5), seed = 25))
  bs <- suppressWarnings(fit_bootstrap(fit, B = 5, workers = 1, seed = 2))
  saved <- tempfile(fileext = ".rds")
  on.exit(unlink(saved), add = TRUE)
  saveRDS(list(fit = fit, bootstrap = bs), saved)
  restored <- readRDS(saved)
  expect_no_error(.validate_fit_bootstrap(restored$bootstrap, restored$fit))

  bad_dir <- tempfile("bad-bootstrap-")
  expect_error(save_outputs(fit, bad_dir, formats = "png",
                            item_plots = FALSE, bootstrap = list()),
               "current fit_bootstrap")
  expect_false(dir.exists(bad_dir))

  wrong_dir <- tempfile("wrong-bootstrap-")
  expect_error(save_outputs(other, wrong_dir, formats = "png",
                            item_plots = FALSE, bootstrap = bs),
               "different fitted model")
  expect_false(dir.exists(wrong_dir))

  reordered <- fit
  reordered$X <- reordered$X[nrow(reordered$X):1L, , drop = FALSE]
  expect_error(.validate_fit_bootstrap(bs, reordered),
               "different fitted model")

  bd <- simulate_btl(5, 12, reps_per_pair = 3, seed = 26)
  bt <- btl(bd, "object_a", "object_b", winner = "winner", judge = "judge")
  expect_error(save_outputs(bt, tempfile("wrong-family-"),
                            formats = "png", item_plots = FALSE,
                            bootstrap = bs), "different model family")

  html <- tempfile(fileext = ".html")
  expect_error(report_html(fit, html, bootstrap = list()),
               "current fit_bootstrap")
  expect_false(file.exists(html))
  expect_error(report_document(fit, html, bootstrap = list()),
               "current fit_bootstrap")
})

test_that("exports accept dimensionality only from the fitted model being reported", {
  fit <- rasch(simd2(160, seq(-1, 1, length.out = 5), seed = 241))
  other <- rasch(simd2(160, seq(-1, 1, length.out = 5), seed = 242))
  dimensionality <- .scree_analysis(fit, n_components = 3,
                                    parallel = FALSE)
  expect_no_error(.validate_scree_result(dimensionality, fit))
  expect_error(.validate_scree_result(dimensionality, other),
               "this fitted model")

  out <- tempfile("wrong-dimensionality-")
  expect_error(save_outputs(other, out, formats = "png", item_plots = FALSE,
                            dimensionality = dimensionality),
               "this fitted model")
  expect_false(dir.exists(out))

  changed <- dimensionality
  changed$eigenvalue[1] <- changed$eigenvalue[1] + 1
  expect_error(report_html(fit, tempfile(fileext = ".html"),
                           dimensionality = changed),
               "plot_scree")

  subtest <- dimensionality_test(
    fit, items_positive = colnames(fit$X)[1:2],
    items_negative = colnames(fit$X)[3:5], min_score_points = 2)
  expect_no_error(.validate_dimensionality_test(subtest, fit))
  out <- tempfile("subtest-export-")
  on.exit(unlink(out, recursive = TRUE), add = TRUE)
  save_outputs(fit, out, formats = "png", dpi = 72,
               item_plots = FALSE, subtest = subtest)
  exported <- read.csv(file.path(
    out, "tables", "unidimensionality_t_test.csv"),
    stringsAsFactors = FALSE)
  expect_identical(exported$split, "manual")
  expect_equal(exported$prop_significant, subtest$prop_significant)
  expect_equal(exported$items_positive,
               paste(subtest$items_positive, collapse = ";"))
  expect_error(.validate_dimensionality_test(subtest, other),
               "this fitted model")
  expect_error(save_outputs(other, tempfile("wrong-subtest-"),
                            formats = "png", item_plots = FALSE,
                            subtest = subtest), "this fitted model")
  altered <- subtest
  altered$prop_significant <- altered$prop_significant + .01
  expect_error(report_html(fit, tempfile(fileext = ".html"),
                           subtest = altered), "dimensionality_test")
})

test_that("exports turn an unavailable person-subset test into a note", {
  d <- simulate_rasch(100, 6, seed = 243)
  repeated <- d[rep(seq_len(nrow(d)), each = 2L), ]
  fit <- rasch(repeated, id = "id", items = sprintf("I%02d", 1:6))
  z <- .dimensionality_or_note(fit)
  expect_match(z$note, "requires one response row per person")
  expect_true(is.na(z$multidimensional))
})

test_that("exports accept DIF only from the fitted model being reported", {
  X <- simd2(160, seq(-1, 1, length.out = 5), seed = 27)
  g1 <- rep(c("A", "B"), each = 80)
  g2 <- rep(c("A", "B"), times = 80)
  fit1 <- rasch(data.frame(X, group = g1), factors = "group")
  fit2 <- rasch(data.frame(X, group = g2), factors = "group")
  da <- dif_anova(fit1)
  expect_no_error(.validate_dif_result(da, fit1))
  expect_error(.validate_dif_result(da, fit2), "different fitted model")

  out <- tempfile("wrong-dif-")
  expect_error(save_outputs(fit2, out, formats = "png",
                            item_plots = FALSE, dif = da),
               "different fitted model")
  expect_false(dir.exists(out))
  html <- tempfile(fileext = ".html")
  expect_error(report_html(fit2, html, dif = da), "different fitted model")
  expect_false(file.exists(html))
  expect_error(report_document(fit2, html, dif = da),
               "different fitted model")

  # A report carries the substantive follow-up estimates, not only the
  # omnibus DIF flags. Use a small valid result-shaped table so the output
  # path is tested independently of whether this random fixture flags DIF.
  da$sizes <- data.frame(
    item = "I1", term = "group", level_a = "A", level_b = "B",
    difference = .5, se = .2, z = 2.5, p_adj = .02, practical = TRUE)
  da$followup_algorithm <- "normalized-design-1"
  # This fixture extends a valid result solely to exercise the report block.
  # Re-seal it as an internal constructor would; an unsigned edited result is
  # deliberately rejected by the public export path.
  unsigned <- unclass(da)
  unsigned$result_signature <- NULL
  da$result_signature <- .fit_boot_md5(unsigned)
  magnitude_html <- tempfile(fileext = ".html")
  on.exit(unlink(magnitude_html), add = TRUE)
  expect_no_error(report_html(fit1, magnitude_html, dif = da))
  magnitude_text <- paste(readLines(magnitude_html, warn = FALSE),
                          collapse = "\n")
  expect_match(magnitude_text, "DIF magnitude", fixed = TRUE)
  expect_match(magnitude_text, "<th>z</th>", fixed = TRUE)
  template <- testthat::test_path("..", "..", "inst", "rmarkdown",
                                 "rasch-report.Rmd")
  if (!file.exists(template))
    template <- system.file("rmarkdown", "rasch-report.Rmd", package = "rasch")
  template_text <- paste(readLines(template, warn = FALSE), collapse = "\n")
  expect_match(template_text,
               'cat("\\n## DIF magnitude\\n\\n")', fixed = TRUE)
  expect_match(template_text, '"t", "z", "df"', fixed = TRUE)
  expect_match(template_text,
               'c("item", "level", "gamma", "se", "t", "df", "z", "p_adj"',
               fixed = TRUE)

  # Person identity is part of the statistical design: the same response
  # rows and factors are a repeated-person analysis under one identifier and
  # a between-person analysis under another.
  Xr <- simd2(160, seq(-1, 1, length.out = 5), seed = 28)
  occasion <- rep(c("pre", "post"), each = 80)
  repeated <- rasch(data.frame(id = rep(sprintf("P%03d", 1:80), 2),
                               occasion, Xr),
                    id = "id", factors = "occasion")
  independent <- rasch(data.frame(id = sprintf("R%03d", 1:160),
                                  occasion, Xr),
                       id = "id", factors = "occasion")
  repeated_dif <- dif_anova(repeated)
  expect_identical(repeated_dif$within, "occasion")
  expect_error(.validate_dif_result(repeated_dif, independent),
               "different fitted model")
})

test_that("report_document writes a self-contained HTML report", {
  skip_if_not_installed("rmarkdown")
  skip_if_not(rmarkdown::pandoc_available())
  set.seed(31)
  X <- matrix(rbinom(120 * 5, 1, .5), 120, 5,
              dimnames = list(NULL, paste0("I", 1:5)))
  fit <- rasch(X)
  out <- tempfile(fileext = ".html")
  on.exit(unlink(out), add = TRUE)

  expect_identical(report_document(fit, out, title = "Test analysis"),
                   normalizePath(out))
  expect_gt(file.info(out)$size, 10000)
  html <- paste(readLines(out, warn = FALSE), collapse = "\n")
  expect_match(html, "Test analysis", fixed = TRUE)
  expect_match(html, "Analysis summary", fixed = TRUE)
  expect_match(html, "Reproducibility", fixed = TRUE)
  expect_error(report_document(fit, sub("html$", "pdf", out),
                               format = "html"),
               "extension")
})

test_that("report_document escapes fitted names used as headings", {
  skip_if_not_installed("rmarkdown")
  skip_if_not(rmarkdown::pandoc_available())
  d <- simulate_mfrm(25, 4, 3, seed = 311)
  facet <- "<em>Injected facet</em>\r\n---\r\nForged block"
  names(d)[names(d) == "rater"] <- facet
  fit <- rasch_mfrm(d, "person", "item", "score", facets = facet)
  out <- tempfile(fileext = ".html")
  on.exit(unlink(out), add = TRUE)

  report_document(fit, out, title = "T <script>alert(1)</script>")
  html <- paste(readLines(out, warn = FALSE), collapse = "\n")
  expect_false(grepl("<h2><em>Injected facet</em></h2>", html,
                     fixed = TRUE))
  expect_match(html, "&lt;em&gt;Injected facet&lt;/em&gt;", fixed = TRUE)
  expect_false(grepl("<hr[^>]*>\\s*<p>Forged block", html, perl = TRUE))
  expect_match(html,
               "<h2>&lt;em&gt;Injected facet&lt;/em&gt;[^<]*Forged block</h2>",
               perl = TRUE)
  expect_false(grepl("<script>alert(1)</script>", html, fixed = TRUE))
  expect_match(html, "&lt;script&gt;alert", fixed = TRUE)
})

test_that("the fit and targeting summaries are complete tidy tables", {
  expect_error(fit_summary_table(list()), "fitted model")
  expect_error(targeting_table(list()), "response-data fit")
  set.seed(1)
  d <- seq(-2, 2, length.out = 6)
  X <- matrix(rbinom(300 * 6, 1, plogis(outer(rnorm(300), d, "-"))), 300, 6)
  colnames(X) <- paste0("I", 1:6)
  f <- rasch(X)
  ft <- fit_summary_table(f)
  tt <- targeting_table(f)
  expect_identical(names(ft), c("statistic", "value"))
  expect_identical(names(tt), c("statistic", "value"))
  expect_false(anyNA(ft$value))
  # spot-check against the fit object
  expect_equal(ft$value[ft$statistic == "Model"], f$model)
  expect_lt(abs(as.numeric(ft$value[
    ft$statistic == "Approximate asymptotic total item-trait chi-square"]) -
                f$total_chisq), 5e-4)
  expect_lt(abs(as.numeric(tt$value[tt$statistic == "PSI"]) - f$psi$PSI), 5e-4)
  expect_lt(abs(as.numeric(tt$value[tt$statistic == "Coefficient alpha"]) -
                f$alpha$alpha), 5e-4)
  failed <- f
  failed$est$converged <- FALSE
  expect_error(targeting_table(failed), "did not converge")

  # EFRM records convergence per set-link edge, not in a scalar `link` field.
  # A failed real edge must stop every output that depends on the linked scale.
  bad_link <- f
  class(bad_link) <- c("rasch_efrm", class(bad_link))
  bad_link$alpha_table <- data.frame(set = c("set1", "set2"))
  bad_link$linking <- list(alpha_edges = data.frame(converged = FALSE))
  expect_false(rasch:::.efrm_link_converged(bad_link))
  expect_identical(fit_summary_table(bad_link)$value[
    fit_summary_table(bad_link)$statistic == "Converged"], "NO")
  expect_error(targeting_table(bad_link), "set-unit link did not converge")
  expect_error(score_table(bad_link), "set-unit link did not converge")
  expect_error(person_extrapolated(bad_link), "set-unit link did not converge")
  expect_error(chisq_detail(bad_link, "I1"), "set-unit link did not converge")
  expect_error(test_information(bad_link), "set-unit link did not converge")
  expect_error(guttman_table(bad_link), "set-unit link did not converge")
  expect_error(rasch:::.wright_map_data(bad_link),
               "set-unit link did not converge")
  bad_link_dir <- tempfile("failed-link-export-")
  bad_link_html <- tempfile(fileext = ".html")
  expect_error(save_outputs(bad_link, bad_link_dir, formats = "png",
                            item_plots = FALSE),
               "set-unit link did not converge")
  expect_false(dir.exists(bad_link_dir))
  expect_error(report_html(bad_link, bad_link_html),
               "set-unit link did not converge")
  expect_false(file.exists(bad_link_html))
  expect_error(weighted_person_estimates(
    bad_link, stats::setNames(rep(1, ncol(f$X)), colnames(f$X))),
    "set-unit link did not converge")
  unknown_link <- bad_link
  unknown_link$linking$alpha_edges$converged <- NA
  expect_false(rasch:::.efrm_link_converged(unknown_link))

  one_set <- bad_link
  one_set$alpha_table <- data.frame(set = "set1")
  one_set$linking$alpha_edges <- data.frame()
  expect_true(rasch:::.efrm_link_converged(one_set))
  missing_edges <- bad_link
  missing_edges$linking <- NULL
  expect_false(rasch:::.efrm_link_converged(missing_edges))
  invalid_detail <- f
  invalid_detail$moments$V[, 1] <- 0
  invalid_detail$item_trait$chisq[1] <- NA_real_
  invalid_detail$item_trait$df[1] <- NA_integer_
  invalid_detail$item_trait$p[1] <- NA_real_
  detail <- chisq_detail(invalid_detail, 1)
  expect_true(all(is.na(detail$intervals$residual)))
  expect_true(all(is.na(detail$intervals$chisq)))
  expect_false(any(detail$intervals$used))
  # robust when alpha is not applicable (missing data)
  Xm <- X; Xm[1:150, 1] <- NA; Xm[151:300, 6] <- NA
  fm <- rasch(Xm)
  ttm <- targeting_table(fm)
  expect_true("NA" %in% ttm$value[ttm$statistic == "Coefficient alpha"] ||
              is.finite(as.numeric(ttm$value[ttm$statistic == "Coefficient alpha"])))
  expect_no_error(fit_summary_table(fm))

  unavailable <- f
  unavailable$items$p_adj <- NA_real_
  fu <- fit_summary_table(unavailable)
  expect_identical(fu$value[grepl("Holm p", fu$statistic)], "unavailable")
  expect_match(paste(capture.output(summary(unavailable)), collapse = "\n"),
               "Holm p < 0.05: unavailable", fixed = TRUE)
})

test_that("structural summaries name response cells and withhold alpha", {
  d <- simulate_mfrm(40, 4, 3, seed = 18)
  f <- rasch_mfrm(d, person = "person", item = "item", score = "score",
                  facets = "rater")
  fs <- fit_summary_table(f)
  ts <- targeting_table(f)
  expect_identical(fs$value[fs$statistic == "Estimation"],
                   "pairwise conditional ML over response cells")
  expect_true(any(fs$statistic ==
                    "Approximate asymptotic total response-cell-trait chi-square"))
  expect_true(any(fs$statistic ==
                    "Response cells with approximate asymptotic Holm p < .05"))
  expect_true(any(ts$statistic == "Response-cell location SD"))
  expect_true(any(ts$statistic == "Calibration threshold minimum"))
  expect_identical(ts$value[ts$statistic == "Coefficient alpha"],
                   "not applicable")
  expect_true(is.na(f$alpha$alpha))
  expect_false(f$alpha$design_applicable)
  expect_null(score_table(f))
  expect_error(person_extrapolated(f), "common raw-score conversion")
  saved_before_flag <- f
  saved_before_flag$alpha$design_applicable <- NULL
  expect_error(ctt_table(saved_before_flag), "several frame or facet response cells")
  expect_error(guttman_table(saved_before_flag), "several frame or facet response cells")

  e <- simulate_efrm(n_per_group = 80, items_per_set = 4, n_sets = 1,
                     n_groups = 2, seed = 19)
  tr <- attr(e, "truth")
  ef <- rasch_efrm(e, item_sets = tr$item_sets, groups = "group", id = "id",
                   boot_reps = 0)
  efs <- fit_summary_table(ef)
  expect_identical(efs$value[efs$statistic == "Estimation"],
                   "pairwise conditional ML over response cells")
  expect_true(is.na(ef$alpha$alpha))
  expect_false(ef$alpha$design_applicable)
  expect_null(score_table(ef))
  expect_error(person_extrapolated(ef), "common raw-score conversion")
})

test_that("lr_test refuses PCM fits that are already constrained", {
  set.seed(46)
  tau <- list(c(-1, 0, 1), c(-.7, .1, 1.1), c(-1.2, -.1, .8),
              c(-.8, .2, 1.2), c(-1.1, 0, .9))
  X <- sapply(tau, function(tt) vapply(rnorm(350), function(th)
    sample(0:3, 1, prob = item_moments(th, tt)$P), 0L))
  colnames(X) <- paste0("Q", seq_len(ncol(X)))
  anchored <- rasch(X, anchors = data.frame(item = "Q1", k = 1, tau = -1))
  expect_error(lr_test(anchored), "unrestricted PCM")
  expect_error(lr_test(rasch(X, pc_components = 2)), "unrestricted PCM")
})
