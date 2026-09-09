simP <- function(theta, tau) { x <- 0:length(tau); p <- exp(x * theta - c(0, cumsum(tau))); p / sum(p) }

test_that("location anchoring fixes item means with free thresholds", {
  set.seed(9); Np <- 1200
  loc_true <- seq(-1.6, 1.6, length.out = 8) + 0.5
  tt <- lapply(seq_len(8), function(i) loc_true[i] + c(-0.9, 0, 0.9))
  th <- rnorm(Np, 0.5, 1.3)
  X <- sapply(seq_len(8), function(i)
    sapply(th, function(t) sample(0:3, 1, prob = simP(t, tt[[i]]))))
  colnames(X) <- sprintf("P%d", 1:8)

  fit <- rasch(X, anchors = data.frame(item = c("P1", "P8"), k = NA,
                                       tau = loc_true[c(1, 8)]))
  expect_equal(fit$items$location[c(1, 8)], loc_true[c(1, 8)],
               tolerance = 1e-9, ignore_attr = TRUE)
  expect_equal(fit$items$se[c(1, 8)], c(0, 0), tolerance = 1e-6)
  expect_true(all(fit$thresholds$se[fit$thresholds$item == 1] > 0))
  expect_lt(sqrt(mean((fit$items$location[2:7] - loc_true[2:7])^2)), 0.15)

  # mixed: mean anchor and threshold anchor together
  fit2 <- rasch(X, anchors = data.frame(item = c("P1", "P8"), k = c(NA, 2),
                                        tau = c(loc_true[1], tt[[8]][2])))
  expect_equal(mean(fit2$tau_list[[1]]), loc_true[1], tolerance = 1e-9)
  expect_identical(fit2$tau_list[[8]][2], tt[[8]][2])
  expect_error(rasch(X, anchors = data.frame(item = "P1", k = c(NA, 1),
                                             tau = c(0, 0))),
               "both a location anchor and threshold anchors")

  # Translating all location anchors changes only the origin. Thresholds and
  # persons move together; likelihood, residuals and covariance do not.
  shifted <- data.frame(item = c("P1", "P8"), k = NA_real_,
                        tau = loc_true[c(1, 8)] + 3)
  fit_shifted <- rasch(X, anchors = shifted)
  expect_equal(fit_shifted$thresholds$tau, fit$thresholds$tau + 3,
               tolerance = 1e-8)
  expect_equal(fit_shifted$person$theta, fit$person$theta + 3,
               tolerance = 1e-7)
  expect_equal(fit_shifted$est$loglik, fit$est$loglik, tolerance = 1e-8)
  expect_equal(fit_shifted$est$cov_tau, fit$est$cov_tau, tolerance = 1e-10)
  expect_equal(fit_shifted$residuals, fit$residuals, tolerance = 1e-7)
  expect_equal(fit$est$n_parameters, sum(fit$m) - 2L)
})

test_that("split_items resolves planted uniform DIF", {
  set.seed(4); Np <- 1500; L <- 10
  d <- scale(seq(-2, 2, length.out = L), scale = FALSE)[, 1]
  grp <- rep(c("ref", "foc"), each = Np / 2); th <- rnorm(Np, 0, 1.4)
  shift <- matrix(0, Np, L); shift[grp == "foc", 3] <- 1.0
  X <- matrix(rbinom(Np * L, 1, plogis(outer(th, d, "-") - shift)), Np, L)
  colnames(X) <- sprintf("G%02d", 1:L)
  fit <- rasch(data.frame(X, grp = grp), factors = "grp", n_groups = 6)
  expect_true(dif_anova(fit)$summary$uniform_DIF[3])

  fit2 <- split_items(fit, "G03", by = "grp")
  locs <- setNames(fit2$items$location, fit2$items$item)
  expect_true(all(c("G03 (foc)", "G03 (ref)") %in% names(locs)))
  expect_equal(unname(locs["G03 (foc)"] - locs["G03 (ref)"]), 1.0, tolerance = 0.3)
  expect_false(any(dif_anova(fit2)$summary$uniform_DIF, na.rm = TRUE))
  expect_true(any(grepl("split", fit2$notes)))
  expect_error(split_items(fit, "G03", by = "nope"), "not a person factor")
})

test_that("equate_tests flags drifted common items only", {
  set.seed(5); L <- 12; d <- seq(-2, 2, length.out = L)
  mk <- function(drift = 0) {
    dd <- d; dd[4] <- dd[4] + drift
    X <- matrix(rbinom(900 * L, 1, plogis(outer(rnorm(900), dd, "-"))), 900, L)
    colnames(X) <- sprintf("I%02d", 1:L); rasch(X)
  }
  f1 <- mk()
  eq_dep <- equate_tests(f1, mk())
  expect_false(eq_dep$inferential)
  expect_true(all(is.na(eq_dep$table$p)))
  expect_match(eq_dep$note, "independence")
  eq0 <- equate_tests(f1, mk(), independent = TRUE)
  expect_equal(sum(eq0$table$drift), 0)
  expect_gt(eq0$correlation, 0.99)

  bad_covariance <- f1
  bad_covariance$est$cov_tau[1, 1] <- -1e6
  guarded <- equate_tests(bad_covariance, mk(), independent = TRUE)
  expect_false(guarded$inferential)
  expect_true(all(is.na(guarded$table$p_adj)))
  expect_match(guarded$note, "not positive semidefinite")

  eq1 <- equate_tests(f1, mk(drift = 0.8), independent = TRUE)
  expect_equal(eq1$table$p_adj,
               p.adjust(eq1$table$p, method = "holm"))
  expect_identical(eq1$table$item[eq1$table$drift], "I04")

  # reference table path (item bank style)
  bank <- data.frame(item = sprintf("I%02d", 1:L), location = d - mean(d), se = 0.05)
  eqb <- equate_tests(f1, bank)
  expect_equal(eqb$n, L)
  expect_false(eqb$inferential)
  expect_match(eqb$note, "joint item-location covariance")
  expect_error(equate_tests(f1, data.frame(item = "ZZ", location = 0)),
               "at least two common items")
})

test_that("equating refuses fitted models from the other response family", {
  rf <- rasch(simulate_rasch(100, 5, seed = 740))
  d <- simulate_btl(5, 12, 6, seed = 741)
  bf <- btl(d, "object_a", "object_b", winner = "winner", judge = "judge")

  expect_error(equate_tests(rf, bf), "different response scale")
  expect_error(btl_equate(bf, rf), "different response scale")
})

test_that("split-item grouping must be an ordinary vector", {
  d <- simulate_rasch(80, 5, seed = 742)
  g <- rep(c("a", "b"), each = 40)
  fit <- rasch(data.frame(d[, -1], group = g), factors = "group")
  expect_error(split_items(fit, "I01", by = matrix(g, ncol = 1L)),
               "ordinary grouping vector")
})

test_that("interactive facet mode recovers a planted item-by-rater effect", {
  set.seed(6); Np <- 500
  persons <- sprintf("P%03d", 1:Np); raters <- paste0("R", 1:4)
  th <- setNames(rnorm(Np, 0, 1.3), persons)
  tau <- list(A = c(-1, 1), B = c(-0.5, 1.2), C = c(-1.2, 0.4), D = c(0, 0.8))
  d <- expand.grid(person = persons, item = names(tau), rater = raters,
                   stringsAsFactors = FALSE)
  d$score <- mapply(function(p, i, r) {
    extra <- if (i == "B" && r == "R2") 0.9 else 0
    sample(0:2, 1, prob = simP(th[p], tau[[i]] + extra))
  }, d$person, d$item, d$rater)

  fit <- rasch_mfrm(d, "person", "item", "score", facets = "rater",
                    interaction = "rater")
  expect_true(fit$est$converged)
  ie <- fit$interaction_effects
  expect_equal(nrow(ie), 16)
  # double sum-to-zero margins
  expect_equal(max(abs(tapply(ie$gamma, ie$item, sum))), 0, tolerance = 1e-8)
  expect_equal(max(abs(tapply(ie$gamma, ie$level, sum))), 0, tolerance = 1e-8)
  top <- ie[which.max(abs(ie$gamma)), ]
  expect_identical(top$item, "B"); expect_identical(top$level, "R2")
  expect_equal(top$gamma, 0.9 * (1 - 1/4 - 1/4 + 1/16), tolerance = 0.2)
  expect_equal(fit$interaction_test$df, (4 - 1) * (4 - 1))
  expect_true(is.finite(fit$interaction_test$p))
  expect_true(all(c("t", "df", "p_adj", "significant") %in%
                    names(fit$interaction_effects)))
  expect_false("z" %in% names(fit$interaction_effects))
  expect_error(rasch_mfrm(d, "person", "item", "score", facets = "rater",
                          interaction = "nope"), "must name one of the facets")
})

test_that("MFRM retains independent-person sandwich support metadata", {
  d <- simulate_mfrm(8, 3, 2, n_categories = 2, seed = 1)
  fit <- suppressWarnings(
    rasch_mfrm(d, "person", "item", "score", facets = "rater"))
  expect_true(fit$est$converged)
  expect_false(fit$est$cluster_support$repeated)
  expect_false(fit$est$cluster_inference)
  expect_true(all(is.na(fit$est$cov_beta)))
  expect_true(all(is.na(fit$item_effects$se)))
  expect_match(paste(fit$notes, collapse = " "),
               "fewer than 10 independent person units")
})

test_that("MFRM interaction inference also needs calibration support", {
  d <- simulate_mfrm(35, 6, 6, n_categories = 2, seed = 913)
  fit <- suppressWarnings(rasch_mfrm(
    d, "person", "item", "score", facets = "rater",
    interaction = "rater"))
  expect_true(all(fit$interaction_support$n_persons >=
                    fit$interaction_support$minimum_required))
  expect_false(fit$est$cluster_inference)
  expect_true(all(is.na(fit$interaction_effects$df)))
  expect_true(all(is.na(fit$interaction_effects$p)))
  expect_false(fit$interaction_test$inference_available)
  expect_true(is.na(fit$interaction_test$p))
  printed <- capture.output(print(fit))
  expect_true(any(grepl("omnibus test: inference unavailable", printed,
                        fixed = TRUE)))
  expect_true(any(grepl("exploratory cell probabilities unavailable", printed,
                        fixed = TRUE)))
  expect_false(any(grepl("p = $", printed)))
  expect_match(paste(fit$notes, collapse = " "),
               "underlying item-parameter sandwich covariance")
})

test_that("factorial DIF: full table, logit follow-ups, main-effects mode", {
  set.seed(1); n <- 1500
  d <- seq(-1.5, 1.5, length.out = 8)
  g1 <- rep(c("a", "b"), each = n / 2)                      # 2 levels, DIF on I3
  g2 <- sample(c("x", "y", "z"), n, replace = TRUE)         # 3 levels, DIF on I6
  th <- rnorm(n)
  sh <- matrix(0, n, 8)
  sh[g1 == "b", 3] <- 1
  sh[g2 == "z", 6] <- 1
  X <- sapply(seq_along(d), function(i)
    rbinom(n, 1, plogis(th - d[i] - sh[, i])))
  colnames(X) <- paste0("I", 1:8)
  fit <- rasch(data.frame(X, g1 = g1, g2 = g2), factors = c("g1", "g2"))

  df <- dif_anova(fit, effects = "factorial")
  # the complete ANOVA table: SS/MS columns and the residual row per item
  expect_true(all(c("sum_sq", "mean_sq") %in% names(df$terms)))
  expect_true(all(table(df$terms$item[df$terms$term == "Residuals"]) == 1))
  expect_true(all(is.na(df$terms$p_adj[df$terms$term == "Residuals"])))

  t3 <- df$terms[df$terms$item == "I3", ]
  expect_true(t3$significant[t3$term == "g1"])
  expect_false(t3$superseded[t3$term == "g1"])
  expect_false("tukey" %in% names(df))
  # A three-level main effect is followed by covariance-aware logit contrasts.
  t6 <- df$terms[df$terms$item == "I6", ]
  expect_true(t6$significant[t6$term == "g2"])
  ph6 <- dif_posthoc(fit, "I6", term = "g2")$table
  expect_equal(nrow(ph6), 3)
  expect_lt(min(ph6$p_adj), 0.01)

  # main-effects mode drops the factor-by-factor terms
  dm <- dif_anova(fit, effects = "main")
  expect_false(any(grepl("g1:g2", dm$terms$term)))
  expect_true(dm$terms$significant[dm$terms$item == "I3" & dm$terms$term == "g1"])

  # the joint main-effects model flags the planted g1 DIF on I3
  su <- dif_anova(fit)$summary
  expect_true(su$uniform_DIF[su$term == "g1" & su$item == "I3"])
  expect_true("p_uniform_adj" %in% names(su))
})

test_that("a significant interaction supersedes its main effects", {
  set.seed(2); n <- 1600
  d <- seq(-1, 1, length.out = 6)
  g1 <- rep(c("a", "b"), each = n / 2)
  g2 <- rep(c("x", "y"), times = n / 2)
  th <- rnorm(n)
  # DIF only in the (b, y) cell: pure g1:g2 interaction structure
  sh <- ifelse(g1 == "b" & g2 == "y", 1.2, 0)
  X <- sapply(seq_along(d), function(i)
    rbinom(n, 1, plogis(th - d[i] - if (i == 2) sh else 0)))
  colnames(X) <- paste0("I", 1:6)
  fit <- rasch(data.frame(X, g1 = g1, g2 = g2), factors = c("g1", "g2"))
  df <- dif_anova(fit, effects = "factorial")
  t2 <- df$terms[df$terms$item == "I2", ]
  expect_true(t2$significant[t2$term == "g1:g2"])
  # the cell-shift induces main effects too; they must be marked superseded
  for (tt in c("g1", "g2"))
    if (t2$significant[t2$term == tt]) expect_true(t2$superseded[t2$term == tt])
  expect_false(t2$superseded[t2$term == "g1:g2"])
  # The interaction follow-up is the logit difference-in-differences.
  phi <- dif_posthoc(fit, "I2", term = c("g1", "g2"))$table
  expect_equal(nrow(phi), 1)
  expect_true(phi$practical)
})

test_that("multiple-choice scoring and miskey detection work", {
  set.seed(5); Np <- 800
  th <- rnorm(Np)
  d <- seq(-1.2, 1.2, length.out = 8)
  keyv <- setNames(rep("A", 8), sprintf("M%02d", 1:8))
  keyv["M04"] <- "B"                                  # miskey: true correct is C
  raw <- sapply(seq_along(d), function(i) {
    correct <- if (i == 4) "C" else "A"
    ok <- rbinom(Np, 1, plogis(th - d[i]))
    ifelse(ok == 1, correct,
           sample(setdiff(c("A", "B", "C", "D"), correct), Np, replace = TRUE))
  })
  colnames(raw) <- sprintf("M%02d", 1:8)
  raw[sample(length(raw), 60)] <- ""                  # blanks become missing

  fit <- rasch(raw, key = keyv)
  expect_false(is.null(fit$mc))
  expect_true(all(fit$m == 1))
  expect_true(any(grepl("scored 0/1 against the key", fit$notes)))
  # scoring matches a manual comparison (blanks NA)
  manual <- ifelse(raw[, "M01"] == "", NA_integer_,
                   as.integer(raw[, "M01"] == "A"))
  expect_identical(unname(fit$X[, "M01"]), manual)

  da <- distractor_analysis(fit)
  m4 <- da[da$item == "M04", ]
  expect_true(m4$flag[m4$option == "C"])              # the real correct option
  expect_false(any(m4$flag[m4$option != "C"]))
  expect_equal(sum(da$flag[da$item != "M04"]), 0)     # clean items stay clean
  # the keyed option carries the top point-biserial on clean items
  clean <- da[da$item == "M01", ]
  expect_identical(clean$option[which.max(clean$point_biserial)], "A")

  # key as a data frame, case-insensitive matching
  fit2 <- rasch(raw, key = data.frame(item = names(keyv), key = tolower(keyv)))
  expect_equal(fit2$items$location, fit$items$location)
  expect_error(distractor_analysis(rasch(fit$X)), "no key")
})

test_that("dimensionality: 10-component PCA, scree and manual subsets", {
  set.seed(2); Np <- 1200; L <- 16
  d <- scale(seq(-2, 2, length.out = L), scale = FALSE)[, 1]
  thA <- rnorm(Np, 0, 1.4); thB <- 0.3 * thA + sqrt(1 - 0.3^2) * rnorm(Np, 0, 1.4)
  XA <- matrix(rbinom(Np * 8, 1, plogis(outer(thA, d[1:8], "-"))), Np, 8)
  XB <- matrix(rbinom(Np * 8, 1, plogis(outer(thB, d[9:16], "-"))), Np, 8)
  X <- cbind(XA, XB); colnames(X) <- sprintf("D%02d", 1:16)
  fit <- rasch(X)

  pc <- residual_pca(fit)
  expect_equal(ncol(pc$loadings_matrix), 11)        # item + PC1..PC10
  expect_equal(nrow(pc$eigen_table), 10)
  expect_true(all(diff(pc$eigen_table$eigenvalue) <= 1e-10))
  expect_equal(pc$eigen_table$cumulative[10],
               sum(pc$eigenvalues[1:10]) / sum(pc$eigenvalues))
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  et <- plot_scree(fit)
  expect_equal(nrow(et), 10)

  # default split detects the planted second dimension; binomial CI fields present
  dt <- dimensionality_test(fit, min_score_points = 2)
  expect_true(is.na(dt$multidimensional))
  expect_true(dt$binomial_multidimensional)
  expect_identical(dt$split, "residual component 1")
  expect_true(dt$ci[1] >= 0 && dt$ci[2] <= 1 && dt$ci[1] < dt$ci[2])
  expect_true(dt$n + dt$n_excluded_extreme >= dt$n)

  # manual subsets matching the true structure also detect it
  dtm <- dimensionality_test(fit, items_positive = sprintf("D%02d", 1:8),
                             items_negative = sprintf("D%02d", 9:16),
                             min_score_points = 2)
  expect_true(dtm$multidimensional)
  expect_identical(dtm$split, "manual")
  expect_gt(dtm$prop_significant, 0.05)

  expect_error(dimensionality_test(fit, items_positive = sprintf("D%02d", 1:8)),
               "both item subsets")
  expect_error(dimensionality_test(fit, items_positive = c("D01", "D02"),
                                   items_negative = c("D02", "D03")),
               "disjoint")
  expect_error(dimensionality_test(fit, items_positive = c("D01", "ZZ"),
                                   items_negative = c("D03", "D04")),
               "not in the fit")
  expect_error(dimensionality_test(
    fit, items_positive = matrix(c("D01", "D02"), nrow = 1L),
    items_negative = c("D03", "D04")),
    "ordinary vector")
  expect_error(dimensionality_test(
    fit, items_positive = c("D01", NA_character_),
    items_negative = c("D03", "D04")),
    "ordinary vector")
  expect_error(combine_items(
    fit, list(matrix(c("D01", "D02"), nrow = 1L))),
    "ordinary vector")
  expect_error(dimensionality_magnitude(
    fit, list(matrix(sprintf("D%02d", 1:8), nrow = 1L),
              sprintf("D%02d", 9:16))),
    "ordinary vector")
})

test_that("compare_fits contrasts nested models on the same data", {
  set.seed(3); Np <- 700
  simP <- function(th, tau) { x <- 0:length(tau); p <- exp(x * th - c(0, cumsum(tau))); p / sum(p) }
  th <- rnorm(Np)
  loc <- seq(-1, 1, length.out = 6); step <- c(-0.8, 0, 0.8)
  X <- sapply(loc, function(b) sapply(th, function(t)
    sample(0:3, 1, prob = simP(t, b + step))))
  colnames(X) <- paste0("R", 1:6)
  pcm <- rasch(X, model = "PCM"); rsm <- rasch(X, model = "RSM")

  cmp <- compare_fits(PCM = pcm, RSM = rsm)
  expect_s3_class(cmp, "rasch_compare")
  expect_identical(attr(cmp, "reference"), "PCM")
  expect_true(all(cmp$same_data))
  # RSM is nested in PCM: fewer parameters, lower (or equal) loglik
  expect_lt(cmp$parameters[2], cmp$parameters[1])
  expect_lte(cmp$two_delta_ll[2], 1e-8)
  expect_equal(cmp$delta_parameters[2],
               cmp$parameters[2] - cmp$parameters[1])
  # different data -> no loglik comparison, descriptive columns still there
  X2 <- X[, 1:5]
  cmp2 <- compare_fits(full = pcm, short = rasch(X2))
  expect_false(cmp2$same_data[2])
  expect_true(is.na(cmp2$two_delta_ll[2]))
  expect_true(all(is.finite(cmp2$chisq_per_df)))
  expect_error(compare_fits(pcm), "at least two")

  no_df <- pcm
  no_df$total_df <- 0L
  no_df$total_chisq <- 0
  cmp3 <- compare_fits(valid = pcm, unavailable = no_df)
  expect_true(is.na(cmp3$chisq_per_df[2L]))
})

test_that("maxit and tol are honoured by the estimators", {
  set.seed(4); Np <- 400; L <- 8
  d <- seq(-1.5, 1.5, length.out = L)
  X <- matrix(rbinom(Np * L, 1, plogis(outer(rnorm(Np), d, "-"))), Np, L)
  colnames(X) <- sprintf("I%02d", 1:L)
  # a very loose tolerance stops after the first step; the deliberate
  # non-convergence now warns loudly (that warning is itself under test)
  expect_warning(f_loose <- rasch(X, maxit = 60, tol = 10),
                 "did NOT converge")
  expect_lte(f_loose$est$iterations, 2)
  # tight settings converge as usual and agree with defaults
  f_tight <- rasch(X, maxit = 200, tol = 1e-10)
  f_def <- rasch(X)
  expect_equal(f_tight$items$location, f_def$items$location, tolerance = 1e-6)
})

test_that("MFRM virtual cells accept colon-bearing labels without collision", {
  set.seed(48)
  d <- expand.grid(person = sprintf("P%03d", 1:120),
                   item = c("A:B", "A", "D"),
                   rater = c("C", "B:C"), stringsAsFactors = FALSE)
  d$score <- rbinom(nrow(d), 2, .5)
  f <- rasch_mfrm(d, person = "person", item = "item", score = "score",
                  facets = "rater")
  expect_false(anyDuplicated(f$virtual_map$vkey) > 0L)
  expect_equal(nrow(f$virtual_map), 6L)
})

test_that("a numeric anchor index survives a dropped constant item", {
  set.seed(7)
  X <- matrix(rbinom(300 * 4, 1, plogis(outer(rnorm(300),
    c(-1, 0, 0.5, 1), "-"))), 300, 4)
  colnames(X) <- paste0("I", 1:4)
  anch <- data.frame(item = 3, k = 1, tau = -2.5)
  f_ok <- suppressWarnings(rasch(X, anchors = anch))
  X[, 2] <- 1L
  f_drop <- suppressWarnings(rasch(X, anchors = anch))
  expect_identical(f_ok$refit_spec$anchors$item, "I3")
  expect_identical(f_drop$refit_spec$anchors$item, "I3")
})

test_that("class-interval detail refuses an item with no classified persons", {
  set.seed(7)
  N <- 200; L <- 6
  X <- matrix(rbinom(N * L, 1, plogis(outer(rnorm(N),
    seq(-1.5, 1.5, length.out = L), "-"))), N, L)
  colnames(X) <- paste0("I", 1:L)
  extra <- rbind(matrix(1L, 10, L + 1), matrix(0L, 10, L + 1))
  X7 <- cbind(X, I7 = NA_integer_)
  X7 <- rbind(X7, `colnames<-`(extra, colnames(X7)))
  anch <- data.frame(item = c("I1", "I7"), k = c(1L, 1L), tau = c(-1.5, 0))
  f <- suppressWarnings(rasch(X7, anchors = anch))
  expect_error(chisq_detail(f, "I7"), "no persons in any class interval")
  td <- tempfile("gap-export-")
  on.exit(unlink(td, recursive = TRUE), add = TRUE)
  expect_no_error(suppressWarnings(save_outputs(f, td, formats = "png",
                                               item_plots = FALSE)))
})

test_that("dependence magnitude withholds inference on weak resolved thresholds", {
  set.seed(3)
  N <- 700
  X <- matrix(rbinom(N * 8, 1, plogis(outer(rnorm(N),
    seq(-1.5, 1.5, length.out = 8), "-"))), N, 8)
  X[, 5] <- ifelse(runif(N) < 0.985, X[, 4], 1 - X[, 4])
  colnames(X) <- paste0("I", 1:8)
  dm <- dependence_magnitude(suppressWarnings(rasch(X)),
                             dependent = "I5", independent = "I4")
  expect_true(is.finite(dm$d))
  expect_true(is.na(dm$se) && is.na(dm$p))
  expect_match(dm$note, "weakly identified")
})

test_that("wide-format MFRM scores factor columns by label, not level code", {
  set.seed(31)
  base <- data.frame(person = rep(sprintf("P%02d", 1:60), 3),
    item = rep(c("A", "B", "C"), each = 60),
    score = sample(0:2, 180, TRUE),
    rater = sample(c("r1", "r2"), 180, TRUE))
  wide <- reshape(base, idvar = c("person", "rater"), timevar = "item",
                  direction = "wide")
  names(wide) <- sub("score\\.", "", names(wide))
  w2 <- wide
  for (cn in c("A", "B", "C")) w2[[cn]] <- factor(w2[[cn]], levels = c(0, 2, 1))
  f1 <- rasch_mfrm(wide, person = "person", facets = "rater",
                   items = c("A", "B", "C"))
  f2 <- rasch_mfrm(w2, person = "person", facets = "rater",
                   items = c("A", "B", "C"))
  expect_equal(f1$items$location, f2$items$location)
})

test_that("item and anchor indices are validated, not truncated", {
  set.seed(32)
  X <- matrix(rbinom(300 * 6, 1, plogis(outer(rnorm(300),
    seq(-1, 1, length.out = 6), "-"))), 300, 6)
  colnames(X) <- paste0("I", 1:6)
  f <- suppressWarnings(rasch(X))
  expect_error(chisq_detail(f, "NOPE"), "no such item")
  expect_error(chisq_detail(f, 1.9), "whole numbers")
  expect_error(dependence_magnitude(f, 2.9, 1.9), "whole numbers")
  expect_error(rasch(X, anchors = data.frame(item = 1.9, k = 1, tau = 0)),
               "whole numbers")
  expect_error(rasch(X, n_groups = 1.9), "whole number")
  expect_error(rasch(X, tol = 0), "positive")
})

test_that("anchor table columns cannot hide matrices", {
  set.seed(731)
  X <- matrix(rbinom(120 * 4, 1, 0.5), 120, 4,
              dimnames = list(NULL, paste0("I", 1:4)))
  shaped <- function(column, value) {
    a <- data.frame(item = "I1", k = 1, tau = 0)
    a[[column]] <- I(matrix(value, nrow = 1L))
    a
  }
  for (fun in list(pcml, rasch)) {
    expect_error(fun(X, anchors = shaped("item", "I1")), "anchor `item`")
    expect_error(fun(X, anchors = shaped("k", 1)), "anchor `k`")
    expect_error(fun(X, anchors = shaped("tau", 0)), "anchor `tau`")
    expect_error(fun(X, anchors = shaped("average", TRUE)), "average")
  }
})

test_that("BTL anchors are a plain named vector", {
  d <- data.frame(a = c("A", "A", "B", "B", "C", "C"),
                  b = c("B", "C", "A", "C", "A", "B"),
                  win = c("A", "C", "B", "B", "A", "C"))
  shaped <- matrix(0, 1L, 1L)
  names(shaped) <- "A"
  expect_error(btl(d, "a", "b", "win", anchors = shaped),
               "named numeric vector")
  d$margin <- I(matrix(rep(1, nrow(d)), ncol = 1L))
  expect_error(btl(d, "a", "b", "win", margin = "margin"),
               "margin")
})

test_that("duplicate named mappings are refused", {
  set.seed(34)
  d <- simulate_btl(n_objects = 6, n_judges = 20, reps_per_pair = 4)
  expect_error(btl(d, "object_a", "object_b", winner = "winner",
                   judge = "judge", anchors = c(O2 = -1, O2 = 3)),
               "duplicate anchor")
})

test_that("weak-item explanatory diagnostics withhold their probabilities", {
  set.seed(35)
  Xw <- vapply(1:8, function(i) {
    if (i == 3) sample(c(0L, 1L, 2L), 500, TRUE, prob = c(0.985, 0.01, 0.005))
    else sample(0:2, 500, TRUE)
  }, integer(500))
  colnames(Xw) <- paste0("I", 1:8)
  qq <- data.frame(item = colnames(Xw), z = rep(c(0, 1), 4))
  few <- suppressWarnings(rasch_explanatory(Xw, predictors = qq, formula = ~ z))
  skip_if_not(any(few$thresholds$weak[few$thresholds$item == 3]))
  dg <- explanatory_diagnostics(few)
  expect_true(all(is.na(dg$p[dg$item == "I3"])))
  expect_true(any(is.finite(dg$p_adj[dg$item != "I3"])))
  use <- is.finite(dg$p)
  expect_equal(dg$p_adj[use], p.adjust(dg$p[use], "holm", n = nrow(dg)))
})

test_that("judge diagnostics report count-weighted comparisons", {
  set.seed(36)
  objs <- paste0("O", 1:6)
  beta <- setNames(seq(-1.2, 1.2, length.out = 6), objs)
  pr <- t(combn(objs, 2))
  d <- do.call(rbind, lapply(sprintf("J%d", 1:4), function(j)
    data.frame(judge = j, a = pr[, 1], b = pr[, 2])))
  d$winner <- ifelse(runif(nrow(d)) < plogis(beta[d$a] - beta[d$b]), d$a, d$b)
  d$count <- 5L
  f <- btl(d, "a", "b", winner = "winner", judge = "judge", count = "count")
  tr <- btl_transitivity(f)
  expect_true(all(tr$judges$n_comparisons == 75))
  js <- judge_surprise(f, as.character(tr$judges$judge[1]))
  expect_identical(js$n_comparisons, 75)
})

test_that("pooled MFRM items flow through DIF follow-ups", {
  set.seed(23)
  dm <- simulate_mfrm(n_persons = 200, n_items = 4, n_raters = 3)
  mf <- rasch_mfrm(dm, person = "person", item = "item", score = "score",
                   facets = "rater",
                   factors = data.frame(g = rep(c("x", "y"), length.out = 200)))
  old_split <- rasch:::split_items
  bad_covariance <- testthat::with_mocked_bindings(
    dif_size(mf, "I2", by = "g"),
    split_items = function(...) {
      z <- old_split(...)
      z$est$cov_tau[1L, 1L] <- -1e6
      z
    },
    .package = "rasch")
  expect_true(all(is.na(bad_covariance$levels$location)))
  expect_true(all(is.na(bad_covariance$pairs$difference)))
  expect_match(paste(bad_covariance$notes, collapse = " "),
               "pooling weights")
  dc <- dif_contrasts(mf, items = "I2")
  expect_identical(unique(dc$table$item), "I2")
  dc_all <- dif_contrasts(mf)
  expect_setequal(unique(dc_all$table$item), paste0("I", 1:4))
})

test_that("matrix items selection validates and subsets", {
  set.seed(41)
  X <- matrix(rbinom(300 * 5, 1, plogis(outer(rnorm(300),
    seq(-1, 1, length.out = 5), "-"))), 300, 5)
  colnames(X) <- paste0("I", 1:5)
  expect_error(rasch(X, items = c(1.9, 2.9, 3.9)), "whole numbers")
  expect_identical(suppressWarnings(rasch(X, items = c(2, 3, 4)))$items$item,
                   c("I2", "I3", "I4"))
  expect_error(pcml(X, anchors = data.frame(item = 1.9, k = 1, tau = 0)),
               "whole numbers")
  expect_error(pcml(X, anchors = data.frame(item = 1, k = 1.5, tau = 0)),
               "positive whole threshold")
  expect_error(pcml(X, anchors = data.frame(item = "I1", k = NaN, tau = 0)),
               "positive whole threshold")
  expect_error(rasch(X, anchors = data.frame(item = "I1", k = NaN, tau = 0)),
               "positive whole threshold")
  expect_error(rasch(X, anchors = data.frame(item = "I1", k = Inf, tau = 0)),
               "positive whole threshold")
  expect_error(pcml(X, anchors = data.frame(item = "I1", k = 1,
                                             tau = "zero")),
               "finite numeric")
  expect_error(pcml(X, maxit = 5.5), "iteration cap")
  expect_error(pcml_pc(X, tol = 0), "tolerance")
  expect_error(pcml(X, maxit = .Machine$integer.max + 1), "iteration cap")
  expect_error(residual_pca(suppressWarnings(rasch(X)),
                            n_components = .Machine$integer.max + 1),
               "integer range")
})

test_that("Wald ratios withhold zero without depending on coefficient scale", {
  z <- .wald_ratio(c(1, 0, 2, NA), c(0, 0, 0.5, 1))
  expect_true(all(is.na(z[1:2])))
  expect_equal(z[3], 4)
  expect_true(is.na(z[4]))
  expect_equal(.wald_ratio(1e9, 1), 1e9)
  expect_equal(.wald_ratio(1, 1e-12), 1e12)
  expect_equal(.wald_ratio(1e9, 1) / .wald_ratio(1, 1e-9), 1)
  expect_error(.wald_ratio(1:2, 1:3), "incompatible lengths")
})

test_that("inverse-variance link weights have no absolute numerical floor", {
  v <- c(1, 4, 9) * 1e-20
  expect_equal(.inverse_variance_weights(v), c(1, 1 / 4, 1 / 9))
  expect_equal(.inverse_variance_weights(v * 1e18),
               .inverse_variance_weights(v))
  expect_equal(.inverse_variance_weights(c(0, 0, 1)), c(1, 1, 0))
})

test_that("equating refuses zero estimated variance as an exact anchor", {
  set.seed(221)
  X <- matrix(rbinom(800 * 6, 1, 0.5), 800, 6,
              dimnames = list(NULL, paste0("I", 1:6)))
  f <- rasch(X)
  f$est$cov_tau[,] <- 0
  f$items$se <- 0
  eq <- equate_tests(f, f, independent = TRUE)
  expect_false(eq$inferential)
  expect_equal(eq$n, 0L)
  expect_equal(eq$n_common, 6L)
  expect_true(all(is.na(eq$table$t)))
  expect_true(all(is.na(eq$table$p_adj)))
  expect_match(eq$note, "zero pooled estimated variance but without exact anchors")

  d <- simulate_btl(n_objects = 6, n_judges = 30,
                    reps_per_pair = 20, seed = 222)
  bt <- btl(d, "object_a", "object_b", "winner", judge = "judge")
  bt$cov_beta[,] <- 0
  bt$objects$se <- 0
  be <- btl_equate(bt, bt, independent = TRUE)
  expect_false(be$inferential)
  expect_equal(be$n_inference, 0L)
  expect_equal(be$n_common, 6L)
  expect_true(all(is.na(be$table$t)))
  expect_true(all(is.na(be$table$p_adj)))
  expect_match(paste(be$notes, collapse = " "),
               "zero pooled estimated variance but without exact anchors")
})

test_that("DIF identifiers keep every unknown row distinct", {
  expect_identical(.role_text_values(c(1, NaN, 2)), c("1", NA, "2"))
  expect_true(.same_role_values(c(1, NaN), c("1", NA)))
  expect_false(.same_role_values(c(1, NaN), c("1", "NaN")))
  expect_identical(.factor_keys(data.frame(x = NaN)), "N;")
  zn <- .dif_ids(c(1, 1, NA_real_, NA_real_, NaN))
  expect_identical(zn[1:2], c("1", "1"))
  expect_length(unique(zn[3:5]), 3L)
  expect_false(any(zn[3:5] %in% zn[1:2]))
  zc <- .dif_ids(c("P1", "P1", NA, "", "  "))
  expect_identical(zc[1:2], c("1", "1"))
  expect_length(unique(zc[3:5]), 3L)
  expect_false(any(zc[3:5] %in% zc[1:2]))
})

test_that("numeric NaN judges are missing rather than a shared judge level", {
  d <- simulate_btl(n_objects = 5, n_judges = 12,
                    reps_per_pair = 12, seed = 223)
  d$judge_number <- as.numeric(sub("^J", "", d$judge))
  d$judge_number[1:2] <- NaN
  f <- btl(d, "object_a", "object_b", "winner", judge = "judge_number")
  expect_false(any(f$comparisons$judge == "NaN", na.rm = TRUE))
  expect_equal(f$n_comparisons, nrow(d) - 2L)
})

test_that("every estimator validates its controls and coercions", {
  set.seed(42)
  d6 <- simulate_btl(n_objects = 5, n_judges = 16, reps_per_pair = 3)
  expect_error(btl(d6, "object_a", "object_b", winner = "winner", maxit = Inf),
               "iteration cap")
  expect_error(btl(d6, "object_a", "object_b", winner = "winner",
                   anchors = c(O2 = NaN)), "finite")
  dm <- simulate_mfrm(n_persons = 120, n_items = 4, n_raters = 2)
  expect_error(rasch_mfrm(dm, person = "person", item = "item",
                          score = "score", facets = "rater", tol = -1),
               "tolerance")
  w <- data.frame(person = sprintf("P%02d", 1:40),
    A = factor(sample(c("0", "1", "poor"), 40, TRUE)),
    B = sample(0:1, 40, TRUE), rater = sample(c("r1", "r2"), 40, TRUE))
  expect_error(rasch_mfrm(w, person = "person", facets = "rater",
                          items = c("A", "B")), "non-numeric")
  f6 <- btl(d6, "object_a", "object_b", winner = "winner", judge = "judge")
  gmap <- setNames(rep(c("x", "y"), 9), rep(sprintf("J%d", 1:9), each = 2))
  png(tf <- tempfile(fileext = ".png"))
  on.exit({dev.off(); unlink(tf)}, add = TRUE)
  expect_error(plot_btl_icc(f6, "O2", group = gmap, min_n = 1),
               "duplicate judge")
})

test_that("MFRM controls and anchor values are validated", {
  set.seed(42)
  dm <- simulate_mfrm(n_persons = 120, n_items = 4, n_raters = 2)
  expect_error(rasch_mfrm(dm, person = "person", item = "item",
                          score = "score", facets = "rater", n_groups = 1.9),
               "whole number")
  X <- matrix(rbinom(300 * 5, 1, 0.5), 300, 5)
  colnames(X) <- paste0("I", 1:5)
  expect_error(pcml(X, anchors = data.frame(item = 2, k = 1, tau = Inf)),
               "finite")
})

test_that("EFRM data selection refuses silent misfits", {
  set.seed(51)
  de <- simulate_efrm(n_per_group = 100, items_per_set = 5, n_sets = 2,
                      n_groups = 2, seed = 51)
  ts <- attr(de, "truth")$item_sets
  expect_error(rasch_efrm(de, item_sets = ts, groups = "group",
                          id = "not_a_column", boot_reps = 0), "not found")
  expect_error(rasch_efrm(de, item_sets = ts, groups = "group",
                          id = 1:10, boot_reps = 0), "entries")
  expect_error(rasch_efrm(de, item_sets = ts, groups = "group", id = "id",
                          factors = data.frame(g = 1:5), boot_reps = 0),
               "rows")
  expect_error(rasch_efrm(de, item_sets = ts, groups = "group", id = "id",
                          items = c("S1I01", "NOPE"), boot_reps = 0),
               "not found")
  expect_error(rasch_efrm(de, item_sets = ts, groups = "group", id = "id",
                          min_link_persons = 2.5, boot_reps = 0),
               "whole number")
})

test_that("selection boundaries refuse silent alterations", {
  set.seed(52)
  db <- simulate_btl(n_objects = 5, n_judges = 16, reps_per_pair = 3)
  db$cnt <- 2L; db$cnt[1] <- 1.5
  expect_error(btl(db, "object_a", "object_b", winner = "winner",
                   count = "cnt"), "replication counts")
  resp <- matrix(sample(c("A", "B", "C"), 300 * 3, TRUE), 300, 3)
  colnames(resp) <- paste0("Q", 1:3)
  expect_error(rasch(resp, key = data.frame(item = "Q1", option = "A",
                                            score = 1.9)),
               "non-negative integers")
  expect_error(rasch(resp, key = data.frame(item = c("Q1", "Q1", "Q2", "Q3"),
                                            key = c("A", "B", "A", "C"))),
               "duplicate key")
  Xp <- matrix(sample(0:2, 300 * 4, TRUE), 300, 4)
  colnames(Xp) <- paste0("P", 1:4)
  prd <- expand.grid(item = colnames(Xp), threshold = c(1, 1.9))
  prd$w <- rep(c(0, 1), each = 4)
  expect_error(rasch_explanatory(Xp, predictors = prd, formula = ~ w,
                                 level = "threshold"), "positive integers")
  Xd <- matrix(rbinom(200 * 4, 1, 0.5), 200, 4)
  colnames(Xd) <- c("A", "B", "A", "C")
  expect_error(pcml(Xd), "unique")
  expect_error(pcml_pc(Xp, n_components = 1.9), "whole number")
  fd <- suppressWarnings(rasch(cbind(Xd[, c(1, 2, 4)],
    D = rbinom(200, 1, 0.5), E = rbinom(200, 1, 0.5))))
  expect_error(dif_size(fd, "D", by = rep(c("x", "y"), 100), alpha = 2),
               "probability")
  expect_error(tailored_analysis(fd, anchor_items = c("A", "ZZZ")),
               "not in the fit")
})

test_that("EFRM matrix input honours id and factors and validates selection", {
  set.seed(61)
  de <- simulate_efrm(n_per_group = 100, items_per_set = 5, n_sets = 2,
                      n_groups = 2, seed = 61)
  ts <- attr(de, "truth")$item_sets
  Xm <- as.matrix(de[, grep("^S", names(de), value = TRUE)])
  f <- rasch_efrm(Xm, item_sets = ts, groups = de$group, id = de$id,
                  factors = data.frame(g2 = rep(c("u", "v"), 100)),
                  boot_reps = 0)
  expect_identical(as.character(f$person$id[1:5]), as.character(de$id[1:5]))
  expect_true("g2" %in% names(f$person))
  expect_error(rasch_efrm(Xm, item_sets = ts, groups = de$group, id = 1:7,
                          boot_reps = 0), "entries")
  expect_error(rasch_efrm(de, item_sets = ts, groups = "group", id = "id",
                          items = c("S1I01", "S1I01", "S1I02"),
                          boot_reps = 0), "more than once")
  ts_dup <- ts; names(ts_dup) <- c("s1", "s1")
  expect_error(rasch_efrm(de, item_sets = ts_dup, groups = "group",
                          id = "id", boot_reps = 0), "duplicate set name")
})

test_that("remaining selector and count boundaries hold", {
  set.seed(62)
  X <- matrix(rbinom(300 * 6, 1, plogis(outer(rnorm(300),
    seq(-1, 1, length.out = 6), "-"))), 300, 6)
  colnames(X) <- paste0("I", 1:6)
  dd <- data.frame(id = rep(sprintf("P%03d", 1:150), 2), X,
                   occ = rep(c("t1", "t2"), each = 150))
  fd <- rasch(dd, id = "id", factors = "occ")
  expect_error(dif_contrasts(fd, within = "occasion"), "not among")
  Xd <- matrix(rbinom(200 * 4, 1, 0.5), 200, 4)
  colnames(Xd) <- c("A", "B", "A", "C")
  expect_error(pcml_pc(Xd), "unique")
  expect_error(threshold_index(c(1, 1.9)), "whole non-negative")
  set.seed(63)
  db <- simulate_btl(n_objects = 5, n_judges = 16, reps_per_pair = 3)
  db$cnt <- 2L; db$cnt[1:3] <- 0L
  f0 <- btl(db, "object_a", "object_b", winner = "winner", count = "cnt")
  expect_true(any(grepl("zero-count", f0$notes)))
  db$cnt[1] <- -1L
  expect_error(btl(db, "object_a", "object_b", winner = "winner",
                   count = "cnt"), "non-negative")
})

test_that("EFRM accepts a factor vector and refuses matrix duplicates", {
  set.seed(61)
  de <- simulate_efrm(n_per_group = 100, items_per_set = 5, n_sets = 2,
                      n_groups = 2, seed = 61)
  ts <- attr(de, "truth")$item_sets
  Xm <- as.matrix(de[, grep("^S", names(de), value = TRUE)])
  expect_error(rasch_efrm(Xm, item_sets = ts, groups = de$group, id = de$id,
                          items = c("S1I01", "S1I01", "S1I02"),
                          boot_reps = 0), "more than once")
  gv <- rep(c("u", "v"), 100)
  f_df <- rasch_efrm(de, item_sets = ts, groups = "group", id = "id",
                     factors = gv, boot_reps = 0)
  expect_true("gv" %in% names(f_df$person))
  f_mx <- rasch_efrm(Xm, item_sets = ts, groups = de$group, id = de$id,
                     factors = gv, boot_reps = 0)
  expect_true("gv" %in% names(f_mx$person))
})

test_that("EFRM refuses ambiguous value-matched role columns", {
  set.seed(71)
  base <- simulate_efrm(n_per_group = 100, items_per_set = 5, n_sets = 2,
                        n_groups = 2, seed = 71)
  ts <- attr(base, "truth")$item_sets
  item_cols <- unlist(ts)
  d1 <- base; d1$cohort <- rep(c(1, 2), 100)
  expect_error(rasch_efrm(d1, item_sets = ts, groups = "group", id = "id",
                          factors = d1$cohort, boot_reps = 0),
               "identical to a supplied role vector")
  f1 <- rasch_efrm(d1, item_sets = ts, groups = "group", id = "id",
                   factors = d1$cohort, items = item_cols, boot_reps = 0)
  expect_false(any(grepl("cohort", f1$item_arbitrary$item)))
  d2 <- base; d2$gnum <- as.integer(factor(d2$group))
  expect_error(rasch_efrm(d2, item_sets = ts, groups = d2$gnum, id = "id",
                          boot_reps = 0),
               "identical to a supplied role vector")
  f2 <- rasch_efrm(d2, item_sets = ts, groups = d2$gnum, id = "id",
                   items = item_cols, boot_reps = 0)
  expect_false(any(grepl("gnum", f2$item_arbitrary$item)))
})

test_that("a character groups vector by value works with explicit items", {
  set.seed(71)
  base <- simulate_efrm(n_per_group = 100, items_per_set = 5, n_sets = 2,
                        n_groups = 2, seed = 71)
  ts <- attr(base, "truth")$item_sets
  d <- base
  gchr <- as.character(as.integer(factor(d$group)))
  d$group <- gchr
  expect_error(rasch_efrm(d, item_sets = ts, groups = gchr, id = "id",
                          boot_reps = 0), "identical to a supplied role")
  f <- rasch_efrm(d, item_sets = ts, groups = gchr, id = "id",
                  items = unlist(ts), boot_reps = 0)
  expect_false(any(grepl("group", f$item_arbitrary$item)))
  # labels colliding with column names stay a by-value grouping
  d2 <- base
  d2$lab <- ifelse(d2$group == levels(factor(d2$group))[1], "id", "group")
  f2 <- rasch_efrm(d2, item_sets = ts, groups = d2$lab, id = "id",
                   items = unlist(ts), boot_reps = 0)
  gcol <- grep("group|frame", names(f2$person))[1]
  expect_true(nlevels(factor(f2$person[[gcol]])) >= 1)
})

test_that("the person characteristic curve uses the model expectation", {
  set.seed(81)
  th0 <- rnorm(200, 0, 1.2)
  X <- vapply(1:8, function(i) { tt <- c(-1, 0, 1) + (i - 4.5) * 0.45
    vapply(th0, function(t) sample(0:3, 1, prob = item_moments(t, tt)$P), 0L)
  }, integer(200))
  colnames(X) <- paste0("P", 1:8)
  f <- suppressWarnings(rasch(X))
  png(tf <- tempfile(fileext = ".png"))
  on.exit({dev.off(); unlink(tf)}, add = TRUE)
  expect_no_error(plot_pcc(f, 5))
  Xd <- matrix(rbinom(200 * 6, 1, 0.5), 200, 6)
  colnames(Xd) <- paste0("D", 1:6)
  expect_no_error(plot_pcc(suppressWarnings(rasch(Xd)), 3))
})

test_that("utility boundaries from the eighth review round hold", {
  set.seed(92)
  X <- matrix(rbinom(300 * 6, 1, plogis(outer(rnorm(300),
    seq(-1, 1, length.out = 6), "-"))), 300, 6)
  colnames(X) <- paste0("I", 1:6)
  fd <- suppressWarnings(rasch(X))
  expect_error(combine_items(fd, groups = list()), "non-empty")
  expect_error(compare_fits(fd, fd, reference = 1.9), "whole number")
  expect_error(compare_fits(a = fd, a = fd), "fit names must be unique")
  expect_error(compare_fits(fd, fit1 = fd), "fit names must be unique")
  expect_error(item_moments(c(0, 1), c(-1, 1)), "one finite location")
  expect_error(item_moments(0, c(-1, 1), disc = 0), "positive finite")
  dd <- data.frame(pid = rep(1:100, 2), t = rep(1:2, each = 100),
                   Q1 = rbinom(200, 1, 0.5), Q2 = rbinom(200, 1, 0.5))
  expect_error(stack_data(dd, "pid", "t", c("Q1", "Q1")), "more than once")
  expect_error(stack_data(dd, "pid", "t", c("pid", "Q1")), "cannot also")
  set.seed(93)
  db <- simulate_btl(n_objects = 5, n_judges = 16, reps_per_pair = 3)
  f6 <- btl(db, "object_a", "object_b", winner = "winner", judge = "judge")
  expect_error(judge_surprise(f6, c("J1", "J2")), "one judge")
  expect_error(judge_surprise(f6, data.frame(judge = c("J1", "J2"))),
               "one judge")
  expect_error(judge_pair_surprise(
    f6, data.frame(judge = c("J1", "J2"))), "one judge")
  expect_error(btl_next_pairs(f6, weight_se = NA), "TRUE or FALSE")
  expect_error(threshold_index(numeric(0)), "at least one")
  png(tf <- tempfile(fileext = ".png"))
  on.exit({dev.off(); unlink(tf)}, add = TRUE)
  expect_error(plot_kidmap(fd, 1, level = 2), "probability")
})

test_that("the ninth review round's regressions and boundaries hold", {
  set.seed(92)
  X <- matrix(rbinom(300 * 6, 1, plogis(outer(rnorm(300),
    seq(-1, 1, length.out = 6), "-"))), 300, 6)
  colnames(X) <- paste0("I", 1:6)
  fd <- suppressWarnings(rasch(X))
  fb <- suppressWarnings(rasch(X[, 1:4]))
  expect_no_error(compare_fits(full = fd, sub = fb, reference = "sub"))
  expect_error(compare_fits(fd, fb, reference = 1.9), "whole number")
  expect_error(person_wle(list(c(-1, 1), 0), disc = c(1, 2)),
               "common discrimination")
  expect_no_error(combine_items(fd, groups = c("I2", "I3")))
  expect_identical(nrow(threshold_index(c(0, 0))), 0L)
  dd <- data.frame(pid = rep(1:50, 2), t = rep(1:2, each = 50),
                   Q1 = rbinom(100, 1, 0.5))
  expect_error(rack_data(dd, "pid", "t", character(0)), "at least one")
  expect_error(stack_data(dd, "pid", "t", character(0)), "at least one")
  # the interval PCC follows the observed locations, not the grid
  set.seed(82)
  th0 <- rnorm(300, 0, 1.5)
  Xp <- vapply(1:8, function(i) { tt <- c(-1, 0, 1) + (i - 4.5) * 1.8
    vapply(th0, function(t) sample(0:3, 1, prob = item_moments(t, tt)$P), 0L)
  }, integer(300))
  colnames(Xp) <- paste0("W", 1:8)
  fp <- suppressWarnings(rasch(Xp))
  png(tf <- tempfile(fileext = ".png"))
  on.exit({dev.off(); unlink(tf)}, add = TRUE)
  expect_no_error(plot_pcc(fp, 4))
})

test_that("compare_fits validates character references helpfully", {
  set.seed(92)
  X <- matrix(rbinom(200 * 5, 1, 0.5), 200, 5)
  colnames(X) <- paste0("I", 1:5)
  fd <- suppressWarnings(rasch(X))
  fb <- suppressWarnings(rasch(X[, 1:4]))
  expect_error(compare_fits(a = fd, b = fb, reference = character(0)),
               "one fit name")
  expect_error(compare_fits(a = fd, b = fb, reference = c("a", "b")),
               "one fit name")
  expect_error(compare_fits(a = fd, b = fb, reference = matrix("a")),
               "one fit name")
})

test_that("the tenth review round's boundaries hold", {
  set.seed(101)
  dm <- simulate_mfrm(n_persons = 150, n_items = 4, n_raters = 3)
  expect_error(rasch_mfrm(dm, person = "person", item = "item",
                          score = "score", facets = "rater",
                          interaction = c("rater", "rater")),
               "exactly one facet")
  resp <- matrix(sample(c("A", "B", "C", "D"), 400 * 3, TRUE), 400, 3)
  colnames(resp) <- paste0("Q", 1:3)
  fmc <- suppressWarnings(rasch(resp, key = c(Q1 = "A", Q2 = "B", Q3 = "C")))
  expect_error(distractor_rescore(fmc, z = -1), "positive finite")
  expect_error(distractor_analysis(fmc, items = character(0)), "at least one")
  X <- matrix(rbinom(300 * 6, 1, plogis(outer(rnorm(300),
    seq(-1, 1, length.out = 6), "-"))), 300, 6)
  colnames(X) <- paste0("I", 1:6)
  fd <- suppressWarnings(rasch(X))
  png(tf <- tempfile(fileext = ".png"))
  on.exit({dev.off(); unlink(tf)}, add = TRUE)
  expect_error(plot_ccc(fd, c("I1", "I2")), "exactly one item")
  expect_error(chisq_detail(fd, c("I1", "I2")), "exactly one item")
  qq <- data.frame(item = colnames(X), z = rep(c(0, 1), 3))
  fe <- rasch_explanatory(X, predictors = qq, formula = ~ z)
  expect_error(relax_explanatory(fe, c("I1", "I2")), "exactly one item")
  expect_error(explanatory_diagnostics(fe, p_adjust = c("holm", "BH")),
               "one method")
})

test_that("equating displays use paired-finite rows and n_common", {
  set.seed(111)
  X <- matrix(rbinom(300 * 6, 1, plogis(outer(rnorm(300),
    seq(-1, 1, length.out = 6), "-"))), 300, 6)
  colnames(X) <- paste0("I", 1:6)
  f1 <- suppressWarnings(rasch(X))
  f2 <- suppressWarnings(rasch(X[1:200, ]))
  png(tf <- tempfile(fileext = ".png"))
  on.exit({dev.off(); unlink(tf)}, add = TRUE)
  expect_no_error(plot_equate(f1, f2, independent = TRUE))
  resp <- matrix(sample(c("A", "B", "C"), 300 * 3, TRUE), 300, 3)
  colnames(resp) <- paste0("Q", 1:3)
  fmc <- suppressWarnings(rasch(resp, key = c(Q1 = "A", Q2 = "B", Q3 = "C")))
  expect_error(plot_distractors(fmc, c("Q1", "Q2")), "exactly one item")
  expect_error(plot_distractors(fmc, matrix("Q1")), "exactly one item")
  expect_error(distractor_analysis(fmc, matrix("Q1")), "ordinary vector")
  expect_error(distractor_analysis(fmc, c("Q1", " Q1 ")),
               "same keyed item")
  # the paired-missing boundary: printing survives a table with no finite
  # location pairs, reporting the correlation and RMSD as unavailable
  set.seed(121)
  db1 <- simulate_btl(n_objects = 6, n_judges = 20, reps_per_pair = 3)
  db2 <- simulate_btl(n_objects = 6, n_judges = 20, reps_per_pair = 3)
  fb1 <- btl(db1, "object_a", "object_b", winner = "winner", judge = "judge")
  fb2 <- btl(db2, "object_a", "object_b", winner = "winner", judge = "judge")
  eq <- btl_equate(fb1, fb2, independent = TRUE)
  eq$table$location_2[] <- NA_real_
  out <- capture.output(print(eq))
  expect_match(out[1], "unavailable")
})

test_that("duplicate and ambiguous selectors cannot alter a test family", {
  set.seed(131)
  X <- matrix(rbinom(300 * 6, 1, plogis(outer(rnorm(300),
    seq(-1, 1, length.out = 6), "-"))), 300, 6)
  colnames(X) <- paste0("I", 1:6)
  dd <- data.frame(id = sprintf("P%03d", 1:300), X,
                   grp = rep(c("a", "b"), 150))
  fd <- rasch(dd, id = "id", factors = "grp")
  expect_error(dif_anova(fd, factors = c("grp", "grpp")), "not stored")
  expect_error(dif_anova(fd, factors = data.frame(g = 1:5)), "300 persons")
  expect_error(dif_contrasts(fd, items = c("I2", "I2")), "more than once")
  expect_error(dif_contrasts(fd, contrasts = list(
    first = c(a = 1, b = -1), first = c(a = -1, b = 1))),
    "duplicate contrast")
  expect_error(dimensionality_test(fd, items_positive = c("I1", "I1", "I2"),
                                   items_negative = c("I3", "I4")),
               "more than once")
  expect_error(rasch(dd, id = "id", factors = c("grp", "grp")),
               "more than once")
  expect_error(dif_size(fd, c("I1", "I2"), by = "grp"), "exactly one item")
  set.seed(132)
  db <- simulate_btl(n_objects = 5, n_judges = 24, reps_per_pair = 20)
  fb <- btl(db, "object_a", "object_b", winner = "winner", judge = "judge")
  jl <- sort(unique(fb$comparisons$judge))
  gmap <- setNames(rep(c("x", "y"), length.out = length(jl)), jl)
  expect_error(btl_dif(fb, factors = gmap, objects = c("O2", "GHOST")),
               "not in the fit")
  expect_error(btl_dif(fb, factors = gmap, objects = c("O2", "O2")),
               "more than once")
  expect_error(btl_dif(fb, factors = list(g = gmap, g = gmap)),
               "duplicate factor")
  expect_error(btl_dif(fb, factors = list(g = gmap, " " = gmap)),
               "non-empty name")
  expect_setequal(unique(btl_dif(
    fb, factors = gmap, objects = " O2 ", min_n = 1)$summary$object), "O2")
  png(tf <- tempfile(fileext = ".png"))
  on.exit({dev.off(); unlink(tf)}, add = TRUE)
  expect_error(plot_item_map(fd, band = -2.5), "positive finite")
  expect_error(plot_pimap(fd, information = NA), "TRUE or FALSE")
  expect_error(plot_scree(fd, parallel = NA), "TRUE or FALSE")
  expect_equal(
    test_information(fd, grid = c(-1, 0, 1), items = " I2 "),
    test_information(fd, grid = c(-1, 0, 1), items = "I2"))
  expect_error(test_information(
    fd, grid = c(-1, 0, 1), items = data.frame(item = "I2")),
    "vector of item names")
  expect_error(test_information(
    fd, grid = c(-1, 0, 1), items = matrix(2L)),
    "vector of item names")
  expect_error(test_information(
    fd, grid = c(-1, 0, 1), items = c(2L, 2L)),
    "same item")
  expect_error(test_information(
    fd, grid = c(-1, 0, 1), items = c("I2", " I2 ")),
    "more than once")
  unconverged <- fd
  unconverged$est$converged <- FALSE
  expect_error(test_information(unconverged), "did not converge")
})

test_that("the fourteenth round's residual boundaries hold", {
  set.seed(141)
  X <- matrix(rbinom(300 * 6, 1, plogis(outer(rnorm(300),
    seq(-1, 1, length.out = 6), "-"))), 300, 6)
  colnames(X) <- paste0("I", 1:6)
  dd <- data.frame(id = sprintf("P%03d", 1:300), X,
                   grp = rep(c("a", "b"), 150))
  fd <- rasch(dd, id = "id", factors = "grp")
  expect_error(dif_contrasts(fd, contrasts = list(
    h = c(a = 1, a = -1, b = 1))), "more than once")
  badf <- data.frame(g = rep(c("a", "b"), 150), g2 = rep(c("x", "y"), 150))
  names(badf) <- c("g", "g")
  expect_error(rasch(dd[, 1:7], id = "id", factors = badf),
               "duplicate factor column")
  de <- simulate_efrm(n_per_group = 100, items_per_set = 5, n_sets = 2,
                      n_groups = 2, seed = 142)
  fe <- rasch_efrm(de, item_sets = attr(de, "truth")$item_sets,
                   groups = "group", id = "id", boot_reps = 0)
  it1 <- fe$virtual_map$item[1]
  png(tf <- tempfile(fileext = ".png"))
  on.exit({dev.off(); unlink(tf)}, add = TRUE)
  expect_error(plot_icc_frames(fe, it1, n_groups = 0), "whole number")
  expect_error(plot_icc_frames(fe, it1, grid = NA), "finite locations")
})

test_that("degenerate reliability and targeting summaries remain available", {
  z <- rasch:::.psi(rep(0, 4), rep(1, 4))
  expect_true(is.na(z$PSI))
  expect_true(is.na(z$separation))
  expect_equal(z$n, 4L)
  expect_true(is.na(rasch:::.safe_cor(c(1, 2, Inf), c(1, 2, 3))))
  expect_true(is.na(rasch:::.dist_stats(rep(2, 4))$skewness))

  tg <- rasch:::.targeting(
    data.frame(theta = c(-1, 1, Inf), extreme = c(FALSE, FALSE, FALSE)),
    data.frame(item = c(1L, 1L, 2L), tau = c(10, 12, 20)))
  expect_equal(tg$item_mean, mean(c(11, 20)))
  expect_equal(tg$person_mean, 0)

  X <- matrix(c(0, 1, 0, 1), ncol = 1,
              dimnames = list(NULL, "I1"))
  it <- rasch:::.item_trait(
    X, list(E = matrix(0.5, 4, 1), V = matrix(c(0, 0, 0.25, 0.25), 4, 1)),
    ci = c(1L, 1L, 2L, 2L))
  expect_true(is.na(it$chisq))
  expect_true(is.na(it$df))
})

test_that("factor structures are unambiguous in every input branch", {
  set.seed(151)
  X <- matrix(rbinom(300 * 6, 1, plogis(outer(rnorm(300),
    seq(-1, 1, length.out = 6), "-"))), 300, 6)
  colnames(X) <- paste0("I", 1:6)
  dd <- data.frame(id = sprintf("P%03d", 1:300), X,
                   grp = rep(c("a", "b"), 150))
  fd <- rasch(dd, id = "id", factors = "grp")
  expect_error(dif_contrasts(fd, contrasts = list(h = c(a = Inf, b = -1))),
               "finite numbers")
  badf <- data.frame(g = rep(c("a", "b"), 150), h = rep(c("x", "y"), 150))
  names(badf) <- c("g", "g")
  expect_error(rasch(X, factors = badf), "duplicate factor column")
  emptyf <- data.frame(g = rep(c("a", "b"), 150))
  names(emptyf) <- ""
  expect_error(rasch(X, factors = emptyf), "non-empty name")
  dm <- simulate_mfrm(n_persons = 100, n_items = 4, n_raters = 2)
  badm <- data.frame(a = rep("u", 100), b = rep("v", 100))
  names(badm) <- c("g", "g")
  expect_error(rasch_mfrm(dm, person = "person", item = "item",
                          score = "score", facets = "rater", factors = badm),
               "duplicate factor column")
})

test_that("the sixteenth round's boundaries hold", {
  set.seed(161)
  db <- simulate_btl(n_objects = 5, n_judges = 16, reps_per_pair = 3)
  dbad <- db; dbad$object_a[3] <- "   "
  expect_error(btl(dbad, "object_a", "object_b", winner = "winner"),
               "blank object")
  dbj <- db; dbj$judge[2] <- ""
  expect_error(btl(dbj, "object_a", "object_b", winner = "winner",
                   judge = "judge"), "blank judge")
  expect_error(btl(db, "object_a", "object_b", winner = "winner",
                   position = NA), "TRUE or FALSE")
  set.seed(162)
  dm <- simulate_mfrm(n_persons = 150, n_items = 4, n_raters = 4)
  fm <- rasch_mfrm(dm, person = "person", item = "item", score = "score",
                   facets = "rater",
                   factors = data.frame(g = rep(c("x", "y"), 75)))
  expect_error(dif_anova(fm, pool_facets = NA), "TRUE or FALSE")
  expect_error(dif_anova(fm, sizes = c(TRUE, FALSE)), "TRUE or FALSE")
  X <- matrix(rbinom(300 * 6, 1, plogis(outer(rnorm(300),
    seq(-1, 1, length.out = 6), "-"))), 300, 6)
  colnames(X) <- paste0("I", 1:6)
  dd <- data.frame(id = sprintf("P%03d", 1:300), X,
                   grp = rep(c("a", "b"), 150))
  fd <- rasch(dd, id = "id", factors = "grp")
  png(tf <- tempfile(fileext = ".png"))
  on.exit({dev.off(); unlink(tf)}, add = TRUE)
  expect_error(plot_icc(fd, "I2", group = "typo"), "no fitted person factor")
  # a short character vector names factors; a short value vector is a
  # length error
  expect_error(plot_icc(fd, "I2", group = c("a", "b")),
               "no fitted person factor")
  expect_error(plot_icc(fd, "I2", group = c(1, 2)), "one value per")
  expect_error(plot_icc(fd, "I2", observed = c(TRUE, FALSE)), "TRUE or FALSE")
  dd2 <- data.frame(pid = c(NA, 1:99, NA, 2:100), t = rep(1:2, each = 100),
                    Q1 = rbinom(200, 1, 0.5))
  expect_error(stack_data(dd2, "pid", "t", "Q1"), "missing or blank person")
  expect_warning(rasch:::.rr_save_plot(function() stop("boom"), "broken",
                                       tempdir(), "png", 7, 5, 96),
                 "could not be drawn")
})

test_that("the seventeenth round's boundaries hold", {
  set.seed(171)
  X <- matrix(rbinom(300 * 6, 1, plogis(outer(rnorm(300),
    seq(-1, 1, length.out = 6), "-"))), 300, 6)
  colnames(X) <- paste0("I", 1:6)
  dd <- data.frame(id = sprintf("P%03d", 1:300), X,
                   group = rep(c("a", "b"), 150),
                   sex = rep(c("m", "f"), each = 150))
  fd <- rasch(dd, id = "id", factors = c("group", "sex"))
  png(tf <- tempfile(fileext = ".png"))
  on.exit({dev.off(); unlink(tf)}, add = TRUE)
  expect_no_error(plot_icc(fd, "I2", group = c("group", "sex")))
  expect_error(plot_icc(fd, "I2", grid = numeric(0)), "finite locations")
  expect_error(plot_icc(fd, "I2", grid = c(0, -1, 1)),
               "strictly increasing")
  expect_error(plot_icc(fd, "I2", n_groups = NA), "whole number")
  expect_error(plot_ccc(fd, "I2", observed = c(TRUE, FALSE)), "TRUE or FALSE")
  expect_error(plot_pimap(fd, xlim = 2), "ascending limits")
  expect_error(plot_threshold_map(fd, order_by_location = NA), "TRUE or FALSE")
  dd2 <- data.frame(pid = c(" ", 1:99, "x", 2:100),
                    t = rep(1:2, each = 100), Q1 = rbinom(200, 1, 0.5))
  expect_error(rack_data(dd2, "pid", "t", "Q1"), "blank person")
  expect_warning(save_item_plots(fd, what = "icc",
                                 file = file.path(tempdir(), "b2.pdf"),
                                 items = c("I1", "GHOST")),
                 "omitted: GHOST")
  dm <- simulate_mfrm(n_persons = 100, n_items = 4, n_raters = 2)
  fm <- rasch_mfrm(dm, person = "person", item = "item", score = "score",
                   facets = "rater")
  expect_error(plot_facets(fm, facet = c("rater", "rater")),
               "one non-empty facet name")
  expect_error(plot_facets(fm, facet = factor("rater")),
               "one non-empty facet name")
  expect_error(plot_facets(fm, facet = matrix("rater", ncol = 1L)),
               "one non-empty facet name")
})

test_that("the eighteenth round's display controls hold", {
  set.seed(181)
  X <- matrix(rbinom(300 * 6, 1, plogis(outer(rnorm(300),
    seq(-1, 1, length.out = 6), "-"))), 300, 6)
  colnames(X) <- paste0("I", 1:6)
  fd <- rasch(X)
  Xp <- matrix(sample(0:2, 300 * 5, TRUE), 300, 5)
  colnames(Xp) <- paste0("P", 1:5)
  fp <- rasch(Xp)
  png(tf <- tempfile(fileext = ".png"))
  on.exit({dev.off(); unlink(tf)}, add = TRUE)
  expect_error(plot_ccc(fp, "P1", observed = TRUE, n_groups = NA),
               "whole number")
  expect_error(plot_threshold_prob(fp, "P1", observed = c(TRUE, FALSE)),
               "TRUE or FALSE")
  expect_error(plot_threshold_prob(fp, "P1", n_groups = 2.5), "whole number")
  expect_error(plot_wright(fd, xlim = 2), "ascending limits")
  expect_error(plot_wright(fd, cex_labels = 0), "positive finite size")
  # limits that admit no thresholds are an empty display, not a drawing
  expect_error(plot_wright(fd, xlim = c(20, 24)), "no thresholds")
  expect_error(plot_kidmap(fd, 1, xlim = c(20, 24)), "no thresholds")
  expect_error(plot_pimap(fd, xlim = c(20, 24)), "no thresholds")
  expect_error(plot_resid_cor(fd, cap = 0), "positive finite correlation")
  expect_error(plot_resid_cor(fd, cap = -0.5), "positive finite correlation")
  expect_no_error(plot_wright(fd, xlim = c(-3, 3)))
  expect_no_error(plot_resid_cor(fd))
})

test_that("the nineteenth round's drawing grids are checked", {
  set.seed(191)
  X <- matrix(rbinom(300 * 6, 1, plogis(outer(rnorm(300),
    seq(-1, 1, length.out = 6), "-"))), 300, 6)
  colnames(X) <- paste0("I", 1:6)
  fd <- rasch(X)
  db <- simulate_btl(n_objects = 5, n_judges = 16, reps_per_pair = 3)
  fb <- btl(db, "object_a", "object_b", winner = "winner", judge = "judge")
  png(tf <- tempfile(fileext = ".png"))
  on.exit({dev.off(); unlink(tf)}, add = TRUE)
  expect_error(plot_tcc(fd, grid = numeric(0)), "finite locations")
  expect_error(plot_tcc(fd, grid = c(0, NA)), "finite locations")
  expect_error(plot_tif(fd, grid = c(1000, 1001)),
               "no positive test information")
  expect_error(plot_btl_targeting(fb, grid = numeric(0)), "finite locations")
  expect_error(plot_btl_categories(fb, grid = numeric(0)), "finite locations")
  expect_error(plot_btl_categories(list(m = 2L)), "paired-comparison")
  expect_error(plot_tif(fd, grid = 0), "finite locations")
  expect_error(plot_btl_icc(fb, fb$objects$object[1], grid = c(0, NA)),
               "finite locations")
  expect_no_error(plot_tcc(fd))
  expect_no_error(plot_btl_targeting(fb))
})

test_that("plot selectors use canonical labels and reject shaped scalars", {
  d <- simulate_rasch(180, 6, model = "PCM", n_categories = 3, seed = 192)
  f <- rasch(d, id = "id")
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  expect_no_error(plot_ccc(f, " I01 "))
  expect_no_error(plot_threshold_prob(f, " I01 "))
  expect_no_error(plot_catfreq(f, " I01 "))
  expect_error(plot_catfreq(f, matrix(1L)), "ordinary vector")
  expect_error(chisq_detail(f, matrix(1L)), "ordinary vector")
  expect_error(plot_pcc(f, data.frame(person = c(1, 2))),
               "exactly one person")
})

test_that("the twentieth round's roles and identifiers hold", {
  # a role vector equal to a real item's responses is ambiguous: refusing it
  # keeps a genuine item from vanishing from the fit, and items= resolves it
  set.seed(201)
  X <- matrix(rbinom(200 * 6, 1, 0.5), 200, 6)
  colnames(X) <- paste0("I", 1:6)
  dd <- data.frame(id = sprintf("P%03d", 1:200), X, stringsAsFactors = FALSE)
  expect_error(rasch(dd, id = "id", factors = dd$I1),
               "identical to a supplied role vector")
  f <- rasch(dd, id = "id", factors = dd$I1, items = paste0("I", 1:6))
  expect_equal(nrow(f$items), 6L)
  expect_true("I1" %in% f$items$item)

  # a recycled identifier would understate the sampling units
  set.seed(202)
  X2 <- matrix(rbinom(400 * 6, 1, plogis(outer(rnorm(400),
    seq(-1, 1, length.out = 6), "-"))), 400, 6)
  colnames(X2) <- paste0("I", 1:6)
  d2 <- data.frame(id = sprintf("P%03d", 1:400), X2,
                   grp = rep(c("a", "b"), 200))
  f2 <- rasch(d2, id = "id", factors = "grp")
  expect_error(dif_anova(f2, id = sprintf("S%03d", 1:200)),
               "one identifier per row")
  expect_no_error(dif_anova(f2, id = sprintf("S%03d", 1:400)))

  # a person with no frame group would be dropped from every set
  de <- simulate_efrm(n_per_group = 100, items_per_set = 5, n_sets = 2,
                      n_groups = 2, seed = 203)
  ts <- attr(de, "truth")$item_sets
  dna <- de; dna$group <- as.character(dna$group); dna$group[1:7] <- NA
  expect_error(rasch_efrm(dna, item_sets = ts, groups = "group", id = "id",
                          boot_reps = 0), "missing or blank frame group")
  dbl <- de; dbl$group <- as.character(dbl$group); dbl$group[1:3] <- "   "
  expect_error(rasch_efrm(dbl, item_sets = ts, groups = "group", id = "id",
                          boot_reps = 0), "missing or blank frame group")

  # blank structural identifiers are not levels
  dm <- simulate_mfrm(n_persons = 80, n_items = 4, n_raters = 3)
  for (col in c("person", "item", "rater")) {
    bad <- dm
    bad[[col]] <- as.character(bad[[col]])
    bad[[col]][bad[[col]] == unique(bad[[col]])[1]] <- " "
    expect_error(rasch_mfrm(bad, person = "person", item = "item",
                            score = "score", facets = "rater"),
                 "blank identifier")
  }
  dbe <- simulate_btl_efrm(n_objects_per_set = 4, n_sets = 2,
                           n_judges_per_panel = 8, n_panels = 2,
                           reps_within = 6, reps_cross = 6, seed = 204)
  os <- attr(dbe, "truth")$object_sets
  bj <- dbe; bj$judge <- as.character(bj$judge); bj$judge[1:3] <- " "
  expect_error(btl_efrm(bj, "object_a", "object_b", winner = "winner",
                        judge = "judge", panels = "panel", object_sets = os,
                        se_method = "conditional"), "blank judge identifier")
  bp <- dbe; bp$panel <- as.character(bp$panel); bp$panel[1:3] <- "  "
  expect_error(btl_efrm(bp, "object_a", "object_b", winner = "winner",
                        judge = "judge", panels = "panel", object_sets = os,
                        se_method = "conditional"), "blank panel identifier")

  fd <- rasch(X2)
  png(tf <- tempfile(fileext = ".png"))
  on.exit({dev.off(); unlink(tf)}, add = TRUE)
  expect_error(plot_kidmap(fd, c(1, 2)), "exactly one person")
  expect_error(plot_resid_dist(fd, bins = 10.5), "whole number")
  expect_error(plot_resid_dist(fd, bins = 0), "whole number")
  d10 <- data.frame(pid = rep(1:100, 2), t = rep(1:2, each = 100),
                    Q1 = rbinom(200, 1, 0.5))
  expect_error(rack_data(d10, character(0), "t", "Q1"),
               "exactly one column")
  expect_error(stack_data(d10, c("pid", "t"), "t", "Q1"),
               "exactly one column")
  expect_error(rack_data(d10, "pid", "nope", "Q1"), "column not found")
})

test_that("the twenty-first round's identifiers hold", {
  # the mapped judge -> panel route validates as the column route does
  dbe <- simulate_btl_efrm(n_objects_per_set = 4, n_sets = 2,
                           n_judges_per_panel = 8, n_panels = 2,
                           reps_within = 6, reps_cross = 6, seed = 204)
  os <- attr(dbe, "truth")$object_sets
  jn <- unique(as.character(dbe$judge))
  pmap <- stats::setNames(rep(c("p1", "p2"), length.out = length(jn)), jn)
  bad <- pmap; bad[1:2] <- " "
  expect_error(btl_efrm(dbe, "object_a", "object_b", winner = "winner",
                        judge = "judge", panels = bad, object_sets = os,
                        se_method = "conditional"),
               "blank panel identifier")
  expect_no_error(btl_efrm(dbe, "object_a", "object_b", winner = "winner",
                           judge = "judge", panels = pmap, object_sets = os,
                           se_method = "conditional"))

  # a blank component survives inside a crossed label, so each source
  # grouping is checked before the cells are crossed
  de <- simulate_efrm(n_per_group = 150, items_per_set = 5, n_sets = 2,
                      n_groups = 2, seed = 211)
  ts <- attr(de, "truth")$item_sets
  de$sex <- rep(c("m", "f"), length.out = nrow(de))
  db <- de; db$sex[1:5] <- "  "
  expect_error(rasch_efrm(db, item_sets = ts, groups = c("group", "sex"),
                          id = "id", boot_reps = 0),
               "blank value\\(s\\) in frame group column")
  dn <- de; dn$sex[1:5] <- NA
  expect_error(rasch_efrm(dn, item_sets = ts, groups = c("group", "sex"),
                          id = "id", boot_reps = 0),
               "missing or blank frame group")
  f <- rasch_efrm(de, item_sets = ts, groups = c("group", "sex"),
                  id = "id", boot_reps = 0)
  expect_equal(levels(factor(f$person[["group:sex"]])),
               c("g1:f", "g1:m", "g2:f", "g2:m"))
})

test_that("restricted models are not silently relaxed by a refit", {
  set.seed(221)
  Xp <- matrix(sample(0:2, 400 * 6, TRUE), 400, 6)
  colnames(Xp) <- paste0("Q", 1:6)
  prd <- data.frame(item = colnames(Xp), w = c(0, 0, 1, 1, 2, 2))
  fe <- rasch_explanatory(Xp, predictors = prd, formula = ~ w,
                          level = "item")
  # the rating refit would drop the design, so the models are not nested
  expect_error(lr_test(fe), "not nested")
  expect_error(lr_test(1), "fitted Rasch-family object")
  expect_no_error(lr_test(rasch(Xp, model = "PCM")))

  Xd <- matrix(rbinom(400 * 6, 1, 0.55), 400, 6)
  colnames(Xd) <- paste0("I", 1:6)
  prd2 <- data.frame(item = colnames(Xd), w = c(0, 0, 1, 1, 2, 2))
  fed <- rasch_explanatory(Xd, predictors = prd2, formula = ~ w,
                           level = "item")
  expect_error(tailored_analysis(fed, chance = 0.25),
               "restricts the item locations")

  # the parallel reference must be analysed under the model that was fitted
  png(tf <- tempfile(fileext = ".png"))
  on.exit({dev.off(); unlink(tf)}, add = TRUE)
  cnt <- 0L
  invisible(utils::capture.output(
    trace(".explanatory_refit_modified", where = asNamespace("rasch"),
          tracer = function() cnt <<- cnt + 1L, print = FALSE)))
  on.exit(try(invisible(utils::capture.output(
    untrace(".explanatory_refit_modified", where = asNamespace("rasch")))),
    silent = TRUE), add = TRUE)
  expect_no_error(plot_scree(fe, parallel = TRUE, reps = 20))
  expect_gt(cnt, 0L)

  db <- simulate_btl(n_objects = 6, n_judges = 20, reps_per_pair = 2,
                     seed = 1)
  obj <- sort(unique(c(as.character(db$object_a), as.character(db$object_b))))
  bp <- data.frame(object = obj, w = seq_along(obj) %% 2)
  fbe <- btl_explanatory(db, predictors = bp, formula = ~ w,
                         object_a = "object_a", object_b = "object_b",
                         winner = "winner", judge = "judge")
  gm <- rep(c("x", "y"), length.out = nrow(fbe$comparisons))
  expect_error(btl_dif(fbe, factors = gm), "explanatory comparison fit")
})

test_that("planted features are planted, or refused", {
  # an item named in several dependence pairs carries every one of them
  sd2 <- simulate_rasch(600, 8, dependence = list(
    pairs = list(c(1, 3), c(2, 3)), strength = 2.5), seed = 222)
  rc <- residual_correlations(rasch(sd2))$matrix
  expect_gt(rc[1, 3], 0.08)
  expect_gt(rc[2, 3], 0.08)

  # a dichotomous item has one threshold and cannot be disordered
  expect_warning(sd3 <- simulate_rasch(300, 6, disordered = 3, seed = 223),
                 "polytomous items only")
  expect_false(any(grepl("disordered",
                         attr(sd3, "truth")$planted %||% character(0))))

  # a unit ratio is a comparison between frames
  expect_error(simulate_efrm(100, 5, n_sets = 1, set_unit_ratio = 1.3),
               "must be 1 when `n_sets` is 1")
  expect_error(simulate_efrm(100, 5, n_groups = 1, group_unit_ratio = 1.2),
               "must be 1 when `n_groups` is 1")
  expect_no_error(simulate_efrm(100, 5, n_sets = 1, set_unit_ratio = 1))
})

test_that("a failed export keeps the caller's device and its own file name", {
  png(tf <- tempfile(fileext = ".png"))
  on.exit({try(dev.off(), silent = TRUE); unlink(tf)}, add = TRUE)
  d0 <- grDevices::dev.cur()
  suppressWarnings(.rr_save_plot(function() plot(1), "x",
                                 file.path(tempdir(), "no_such_dir_zz"),
                                 "png", 7, 5, 96))
  expect_identical(grDevices::dev.cur(), d0)

  # two item names that sanitise alike must not share one file
  expect_equal(length(unique(.rr_safe_stem(c("A/B", "A B", "A_B")))), 3L)
  set.seed(225)
  Xc <- matrix(rbinom(300 * 3, 1, 0.5), 300, 3)
  colnames(Xc) <- c("A/B", "A B", "C")
  od <- file.path(tempdir(), "outprobe")
  unlink(od, recursive = TRUE)
  invisible(save_outputs(rasch(Xc), od, formats = "png", item_plots = TRUE))
  expect_equal(length(list.files(od, pattern = "_icc\\.png$",
                                 recursive = TRUE)), 3L)
})

test_that("every export name path keeps its own file", {
  # objects, like items, can sanitise to the same stem
  set.seed(233)
  db <- simulate_btl(n_objects = 6, n_judges = 24, reps_per_pair = 4)
  ren <- c(O1 = "A/B", O2 = "A B", O3 = "C", O4 = "D", O5 = "E", O6 = "F")
  for (cc in c("object_a", "object_b", "winner"))
    db[[cc]] <- unname(ren[as.character(db[[cc]])])
  fb <- btl(db, "object_a", "object_b", winner = "winner", judge = "judge")
  od <- file.path(tempdir(), "gapsbtl"); unlink(od, recursive = TRUE)
  invisible(suppressWarnings(save_outputs(fb, od, formats = "png",
                                          item_plots = TRUE)))
  expect_equal(length(list.files(od, pattern = "_icc[.]png$",
                                 recursive = TRUE)), 6L)

  # facet tables and facet plots take their names the same way
  set.seed(235)
  np <- 40; items <- paste0("It", 1:3); raters <- paste0("R", 1:3)
  occ <- c("t1", "t2")
  g <- expand.grid(person = sprintf("P%02d", 1:np), item = items,
                   rater = raters, occasion = occ, stringsAsFactors = FALSE)
  th <- stats::setNames(rnorm(np), sprintf("P%02d", 1:np))
  di <- stats::setNames(seq(-1, 1, length.out = 3), items)
  sv <- stats::setNames(c(-0.4, 0, 0.4), raters)
  ov <- stats::setNames(c(-0.2, 0.2), occ)
  g$score <- rbinom(nrow(g), 2, plogis((th[g$person] - di[g$item] -
                                        sv[g$rater] - ov[g$occasion]) / 2))
  names(g)[names(g) == "rater"] <- "r/x"
  names(g)[names(g) == "occasion"] <- "r x"
  fm <- rasch_mfrm(g, person = "person", item = "item", score = "score",
                   facets = c("r/x", "r x"))
  md <- file.path(tempdir(), "gapsmfrm"); unlink(md, recursive = TRUE)
  invisible(suppressWarnings(save_outputs(fm, md, formats = "png")))
  expect_equal(length(grep("facet", basename(list.files(md, pattern = "csv$",
                                                        recursive = TRUE)),
                           value = TRUE)), 2L)
  expect_equal(length(list.files(md, pattern = "facet_severities.*png$",
                                 recursive = TRUE)), 2L)

  # the report closes only devices it opened
  X <- matrix(rbinom(300 * 5, 1, 0.5), 300, 5)
  colnames(X) <- paste0("I", 1:5)
  png(tf <- tempfile(fileext = ".png"))
  on.exit({try(dev.off(), silent = TRUE); unlink(tf)}, add = TRUE)
  d0 <- grDevices::dev.cur()
  rp <- file.path(tempdir(), "gapsrep.html"); unlink(rp)
  invisible(suppressWarnings(report_html(rasch(X), rp)))
  expect_identical(grDevices::dev.cur(), d0)
  expect_true(file.exists(rp))
})

test_that("the explanatory scree reference aligns its person rows", {
  set.seed(241)
  N <- 200
  Xp <- matrix(sample(0:2, N * 6, TRUE), N, 6)
  colnames(Xp) <- paste0("Q", 1:6)
  Xp[7, ] <- NA                     # a respondent with no person estimate
  d <- data.frame(id = sprintf("P%03d", 1:N), Xp,
                  sex = rep(c("m", "f"), length.out = N))
  prd <- data.frame(item = colnames(Xp), w = c(0, 0, 1, 1, 2, 2))
  fe <- rasch_explanatory(d, predictors = prd, formula = ~ w, level = "item",
                          id = "id", factors = "sex", items = colnames(Xp))
  expect_gt(sum(is.na(fe$person$theta)), 0L)
  ref <- .scree_reference(fe, k = 2, reps = 20)
  expect_true(all(is.finite(unlist(ref))))

  # the refit carries the retained rows' identifiers and factors, stays an
  # explanatory fit, and keeps an approved departure in force
  keep <- !is.na(fe$person$theta)
  fr <- relax_explanatory(fe, "Q3", "location")
  rf <- .explanatory_refit_modified(fr, fr$X[keep, , drop = FALSE],
                                    person_rows = which(keep))
  expect_s3_class(rf, "rasch_explanatory")
  expect_equal(nrow(rf$person), sum(keep))
  expect_equal(nrow(rf$factors), sum(keep))
  expect_identical(as.character(rf$person$id),
                   as.character(fe$person$id[keep]))
  expect_equal(ncol(rf$est$B), ncol(fr$est$B))
  expect_gt(ncol(rf$est$B), ncol(fe$est$B))
  expect_equal(nrow(rf$explanatory$relaxations),
               nrow(fr$explanatory$relaxations))
  # a source of the wrong height is refused rather than misaligned
  expect_error(.explanatory_refit_modified(fe, fe$X[keep, , drop = FALSE],
                                           person_rows = seq_len(N)),
               "one fitted row per row")
})

test_that("a batch export that drew nothing does not report success", {
  set.seed(242)
  X <- matrix(rbinom(300 * 4, 1, 0.5), 300, 4)
  colnames(X) <- paste0("I", 1:4)
  f <- rasch(X)
  png(tfd <- tempfile(fileext = ".png"))
  on.exit({try(dev.off(), silent = TRUE); unlink(tfd)}, add = TRUE)
  d0 <- grDevices::dev.cur()

  zp <- file.path(tempdir(), "gapsbatch.zip"); unlink(zp)
  expect_error(save_item_plots(f, what = "icc", file = zp, width = -1),
               "positive finite value")
  expect_false(file.exists(zp))

  zp2 <- file.path(tempdir(), "gapsbatch2.zip"); unlink(zp2)
  expect_error(save_item_plots(f, what = "icc", file = zp2,
                               items = c("GHOST1", "GHOST2")),
               "no archive was created")
  expect_false(file.exists(zp2))
  expect_identical(grDevices::dev.cur(), d0)

  # partial success still produces a real archive, and says what is missing
  zp3 <- file.path(tempdir(), "gapsbatch3.zip"); unlink(zp3)
  expect_warning(save_item_plots(f, what = "icc", file = zp3,
                                 items = c("I1", "GHOST")),
                 "omitted: GHOST")
  expect_true(file.exists(zp3))
  expect_equal(nrow(utils::unzip(zp3, list = TRUE)), 1L)

  zp4 <- file.path(tempdir(), "gapsbatch4.zip"); unlink(zp4)
  invisible(save_item_plots(f, what = "icc", file = zp4))
  cont <- utils::unzip(zp4, list = TRUE)$Name
  expect_equal(length(cont), 4L)
  expect_equal(length(unique(cont)), 4L)

  pdfp <- file.path(tempdir(), "gapsbatch.pdf"); unlink(pdfp)
  expect_error(save_item_plots(f, what = "icc", file = pdfp,
                               items = c("GHOST1", "GHOST2")),
               "no file was created")
  expect_false(file.exists(pdfp))
  invisible(save_item_plots(f, what = "icc", file = pdfp))
  expect_true(file.exists(pdfp))
})

test_that("derived fits keep the controls they were built from", {
  set.seed(251)
  N <- 700; tau <- c(-0.7, 0.7); th <- rnorm(N)
  X <- sapply(seq(-1, 1, length.out = 6), function(dd)
    vapply(th, function(b) sample(0:2, 1, prob = item_moments(b, tau + dd)$P),
           0L))
  colnames(X) <- paste0("Q", 1:6)
  d <- data.frame(id = sprintf("P%03d", 1:N), X,
                  grp = rep(c("a", "b"), N / 2))
  f <- rasch(d, id = "id", factors = "grp", model = "PCM")
  lr <- lr_test(f)
  # the refit must carry the person factors the follow-up analyses need
  expect_false(is.null(lr$fit_rsm$factors))
  expect_no_error(dif_anova(lr$fit_rsm, factors = "grp"))
  ratio <- sum(lr$fit_rsm$item_trait$chisq, na.rm = TRUE) /
    sum(f$item_trait$chisq, na.rm = TRUE)
  expect_gt(ratio, 0.5)

  set.seed(252)
  P <- 0.25 + 0.75 * plogis(outer(rnorm(800), seq(-2, 2, length.out = 10), "-"))
  Xd <- matrix(rbinom(800 * 10, 1, P), 800, 10)
  colnames(Xd) <- paste0("I", 1:10)
  fd <- rasch(Xd)
  ta <- tailored_analysis(fd, chance = 0.25)
  expect_gt(sum(ta$anchored$item_trait$chisq, na.rm = TRUE) /
              sum(fd$item_trait$chisq, na.rm = TRUE), 0.5)

  # a subtest total spans a wider range than any member item
  set.seed(253)
  Xp <- matrix(sample(0:2, 400 * 8, TRUE), 400, 8)
  Xp[1, 1:5] <- 0L; Xp[1, 6:8] <- 2L
  Xp[2, 1:5] <- 2L; Xp[2, 6:8] <- 0L
  colnames(Xp) <- paste0("I", 1:8)
  fp <- rasch(Xp, model = "PCM", na_codes = 9)
  cb <- combine_items(fp, list(c("I1", "I2", "I3", "I4", "I5")))
  sn <- "I1+I2+I3+I4+I5"
  expect_equal(sum(is.na(cb$X[, sn])), 0L)
  expect_true(any(grepl("missing-data code", cb$notes)))
  fp_text <- rasch(Xp, model = "PCM", na_codes = "09")
  cb_text <- combine_items(fp_text,
                           list(c("I1", "I2", "I3", "I4", "I5")))
  expect_equal(cb_text$X[, sn], cb$X[, sn])
  expect_true(any(grepl("missing-data code", cb_text$notes)))

  # the class-interval detail reconciles with the item table it quotes
  cd <- chisq_detail(fd, "I4")
  expect_equal(sum(cd$intervals$chisq[cd$intervals$used]), cd$chisq,
               tolerance = 1e-8)
})

test_that("frame invariance uses the covariance of the centred differences", {
  set.seed(21)
  K <- 6; npg <- 400
  inm <- sprintf("I%02d", 1:K)
  mvec <- c(1, 1, 1, 1, 1, 4)
  delta <- seq(-1, 1, length.out = K)
  tl <- lapply(seq_len(K), function(i)
    if (mvec[i] == 1) delta[i] else delta[i] + seq(-1.2, 1.2, length.out = 4))
  phi <- c(1 / sqrt(1.5), sqrt(1.5))
  grp <- factor(rep(c("g1", "g2"), each = npg)); N <- length(grp)
  theta <- rnorm(N, 0, 1.4)
  X <- matrix(NA_integer_, N, K, dimnames = list(NULL, inm))
  for (col in seq_len(K)) {
    rho <- phi[as.integer(grp)]; m <- mvec[col]; ct <- cumsum(tl[[col]])
    E <- cbind(0, sapply(seq_len(m), function(k) rho * (k * theta - ct[k])))
    P <- exp(E - apply(E, 1, max)); P <- P / rowSums(P)
    cum <- P %*% upper.tri(diag(m + 1L), diag = TRUE)
    X[, col] <- as.integer(rowSums(runif(N) > cum))
  }
  dd <- data.frame(id = sprintf("P%04d", seq_len(N)), X, group = grp,
                   check.names = FALSE)
  fit <- rasch_efrm(dd, item_sets = list(set1 = inm), groups = "group",
                    id = "id", boot_reps = 0)
  fi <- frame_invariance(fit)
  expect_true(all(is.finite(fi$locations$se)))
  # C = I - 1w' is not symmetric when the items differ in maximum score, so
  # Var(Cd) = C S C'. Under the mistaken C S C the five dichotomous items'
  # standard errors come out ~6% too small and the five-category item's 22%
  # too large: c(0.1751, 0.1657, 0.1553, 0.1579, 0.1611, 0.1120).
  expect_equal(round(fi$locations$se, 4),
               c(0.1857, 0.1768, 0.1644, 0.1676, 0.1703, 0.0914),
               tolerance = 1e-3)
  expect_equal(fi$summary$rmse, sqrt(mean(fi$locations$se^2)),
               tolerance = 1e-8)
})

test_that("displays and refusals match the analysis they report", {
  # an EFRM fitted without the bootstrap has no unit standard errors
  de <- simulate_efrm(n_per_group = 120, items_per_set = 5, n_sets = 2,
                      n_groups = 2, seed = 261)
  fe <- rasch_efrm(de, item_sets = attr(de, "truth")$item_sets,
                   groups = "group", id = "id", boot_reps = 0)
  png(tf <- tempfile(fileext = ".png"))
  on.exit({dev.off(); unlink(tf)}, add = TRUE)
  expect_true(all(is.na(fe$frames$se_log_rho)))
  expect_no_error(plot_frames(fe))

  # the observed points follow the fit's own per-item class intervals
  set.seed(262)
  N <- 600; th <- rnorm(N)
  X <- matrix(rbinom(N * 8, 1,
    plogis(outer(th, seq(-1.5, 1.5, length.out = 8), "-"))), N, 8)
  colnames(X) <- paste0("I", 1:8)
  X[1:250, 3] <- NA
  f <- rasch(X)
  expect_false(is.null(f$ci_item))
  i3 <- .item_idx(f, "I3")
  ex <- if (!is.null(f$person$extreme)) f$person$extreme else rep(FALSE, N)
  used <- .fit_class_intervals(f, i3, f$X[, i3], f$person$theta, ex,
                               f$n_groups, FALSE, NULL)
  expect_equal(max(used, na.rm = TRUE),
               nrow(chisq_detail(f, "I3")$intervals))
  # an explicit request is still honoured
  expect_equal(max(.fit_class_intervals(f, i3, f$X[, i3], f$person$theta, ex,
                                        4, TRUE, NULL), na.rm = TRUE), 4)
  expect_no_error(plot_icc(f, "I3", observed = TRUE))
  expect_no_error(plot_ccc(f, "I3", observed = TRUE))
  expect_no_error(plot_threshold_prob(f, "I3", observed = TRUE))
})

test_that("reports, batches and simulations refuse malformed requests", {
  set.seed(271)
  X <- matrix(rbinom(300 * 5, 1, 0.5), 300, 5)
  colnames(X) <- paste0("I", 1:5)
  f <- rasch(X)
  rp <- file.path(tempdir(), "gapsr1.html"); unlink(rp)
  expect_error(report_html(f, rp, dpi = -100), "positive finite value")
  expect_error(report_html(f, rp, dpi = c(100, 200)), "positive finite value")
  expect_error(report_html(f, rp, title = c("A", "B")), "one non-missing title")
  expect_error(report_html(f, rp, title = matrix("A")), "one non-missing title")
  expect_error(report_html(f, c(rp, rp)), "one non-missing path")
  expect_false(file.exists(rp))
  expect_error(report_document(f, file.path(tempdir(), "g.html"),
                               title = c("A", "B")), "one non-missing title")
  expect_error(report_document(f, file.path(tempdir(), "g.html"),
                               title = matrix("A")), "one non-missing title")

  # save_outputs checks before it creates anything
  od <- file.path(tempdir(), "gapsso"); unlink(od, recursive = TRUE)
  expect_error(save_outputs(f, od, width = -1), "positive finite value")
  expect_false(dir.exists(od))
  expect_error(save_outputs(f, od, item_plots = NA), "TRUE or FALSE")
  expect_error(save_outputs(f, c(od, od)), "one non-missing path")

  # sim_apply takes one atomic scalar per replicate
  b <- sim_replicate(simulate_rasch, 2, n_persons = 120, n_items = 5,
                     seed = 271)
  r1 <- sim_apply(b, function(d) data.frame(a = 1))
  expect_equal(attr(r1, "n_failed"), 2L)
  r2 <- sim_apply(b, function(d) list(c(1, 2, 3)))
  expect_equal(attr(r2, "n_failed"), 2L)
  expect_equal(length(r2), 2L)
  r3 <- sim_apply(b, function(d) mean(rasch(d)$items$location))
  expect_equal(attr(r3, "n_failed"), 0L)

  # a sequence-dependent design is not the same data under a different order
  set.seed(272)
  db <- simulate_btl(n_objects = 5, n_judges = 12, reps_per_pair = 3)
  db$ord <- seq_len(nrow(db))
  db2 <- db; db2$ord <- rev(db2$ord)
  p1 <- btl(db, "object_a", "object_b", winner = "winner", judge = "judge",
            order = "ord", position = TRUE)
  p2 <- btl(db2, "object_a", "object_b", winner = "winner", judge = "judge",
            order = "ord", position = TRUE)
  expect_false(compare_fits(seqA = p1, seqB = p2)$same_data[2])
  # orientation stays arbitrary when nothing sequence-dependent is fitted
  sw <- db; a <- as.character(sw$object_a)
  sw$object_a <- as.character(sw$object_b); sw$object_b <- a
  q1 <- btl(db, "object_a", "object_b", winner = "winner", judge = "judge")
  q2 <- btl(sw, "object_a", "object_b", winner = "winner", judge = "judge")
  expect_true(compare_fits(plainA = q1, plainB = q2)$same_data[2])

  # every judge needs a panel, and a lone frame carries no unit or origin
  dbe <- simulate_btl_efrm(n_objects_per_set = 4, n_sets = 2,
                           n_judges_per_panel = 8, n_panels = 2,
                           reps_within = 6, reps_cross = 6, seed = 273)
  os <- attr(dbe, "truth")$object_sets
  jn <- unique(as.character(dbe$judge))
  full <- stats::setNames(rep(c("p1", "p2"), length.out = length(jn)), jn)
  expect_error(btl_efrm(dbe, "object_a", "object_b", winner = "winner",
                        judge = "judge", panels = full[1:(length(full) - 3)],
                        object_sets = os, se_method = "conditional"),
               "missing from the panels map")
  expect_error(simulate_btl_efrm(4, 2, 6, n_panels = 1, panel_units = 1.4,
                                 seed = 1), "must be 1 when `n_panels` is 1")
  expect_error(simulate_btl_efrm(4, n_sets = 1, set_units = 1.3, seed = 1),
               "must be 1 when `n_sets` is 1")
  expect_error(simulate_btl_efrm(4, n_sets = 1, set_origins = 0.5, seed = 1),
               "must be 0 when `n_sets` is 1")
  expect_no_error(simulate_btl_efrm(4, n_sets = 1, set_units = 1,
                                    set_origins = 0, seed = 1))
})

test_that("design-restricted fits are refused where locations are not items", {
  set.seed(281)
  q <- data.frame(item = paste0("I", 1:8), operation = rep(0:1, each = 4),
                  format = rep(c("A", "B"), 4))
  d0 <- -1 + 0.7 * q$operation + 0.4 * (q$format == "B")
  mk <- function(dv, n = 800) { th <- rnorm(n)
    Y <- matrix(rbinom(n * 8, 1, plogis(outer(th, dv, "-"))), n, 8)
    colnames(Y) <- q$item; Y }
  X1 <- mk(d0); X2 <- mk(d0)
  expl <- rasch_explanatory(X1, predictors = q, formula = ~ operation + format)
  ref <- rasch(X2)
  expect_error(equate_tests(expl, ref, independent = TRUE),
               "not defined for an explanatory calibration")
  expect_no_error(equate_tests(rasch(X1), ref, independent = TRUE))
  subs <- list(paste0("I", 1:4), paste0("I", 5:8))
  expect_error(dimensionality_magnitude(expl, subs),
               "not defined for an explanatory calibration")
  expect_no_error(dimensionality_magnitude(rasch(X1), subs))

  # the object design centres on the objects actually calibrated
  set.seed(282)
  db <- simulate_btl(n_objects = 7, n_judges = 20, reps_per_pair = 3)
  ob <- sort(unique(c(as.character(db$object_a), as.character(db$object_b))))
  bp <- data.frame(object = ob, w = seq_along(ob) %% 2)
  fbe <- btl_explanatory(db, predictors = bp, formula = ~ w,
                         object_a = "object_a", object_b = "object_b",
                         winner = "winner", judge = "judge")
  lo <- fbe$objects$location[!(fbe$objects$extreme %in% TRUE)]
  expect_equal(mean(lo), 0, tolerance = 1e-6)
  expect_error(btl_equate(fbe, fbe), "explanatory comparison fit")
})

test_that("simulation plants what it records, and reports what it cannot", {
  # speededness needs a tail
  expect_warning(s2 <- simulate_rasch(200, 2, speeded = 0.3, seed = 284),
                 "at least three items")
  expect_false(any(grepl("peeded", attr(s2, "truth")$planted %||%
                           character(0))))
  # a chain is fine; a source regenerated as a later target is not
  expect_no_error(simulate_rasch(400, 6, dependence = list(
    pairs = list(c(1, 3), c(3, 5))), seed = 285))
  expect_error(simulate_rasch(400, 6, dependence = list(
    pairs = list(c(3, 5), c(1, 3))), seed = 285),
    "source of an earlier dependence")
  # group units follow the group labels, not the sorted level order
  dg <- simulate_efrm(n_per_group = 30, items_per_set = 4, n_sets = 2,
                      n_groups = 12, group_unit_ratio = 2, seed = 286)
  ph <- attr(dg, "truth")$phi
  expect_equal(names(ph), sprintf("g%d", 1:12))
  expect_false(is.unsorted(ph))
  # a bias planted on a rater who answers at random is not planted at all
  expect_error(simulate_mfrm(80, 4, 4, erratic_raters = 0.25,
                             interaction = list(rater = "R1", item = "I1",
                                                bias = 2), seed = 287),
               "erratic raters")
  expect_no_error(simulate_mfrm(80, 4, 4, erratic_raters = 0.25,
                                interaction = list(rater = "R4", item = "I1",
                                                   bias = 2), seed = 287))
})

test_that("a batch archive and a data signature say what they are", {
  set.seed(288)
  X <- matrix(rbinom(300 * 4, 1, 0.5), 300, 4)
  colnames(X) <- paste0("I", 1:4)
  f <- rasch(X)
  zp <- file.path(tempdir(), "gapsfresh.zip"); unlink(zp)
  invisible(save_item_plots(f, what = "icc", file = zp))
  expect_equal(nrow(utils::unzip(zp, list = TRUE)), 4L)
  # zip() appends: a second call must not carry the first call's plots
  invisible(save_item_plots(f, what = "icc", file = zp,
                            items = c("I1", "I2")))
  expect_equal(nrow(utils::unzip(zp, list = TRUE)), 2L)
  expect_equal(.estimation_label(f), "pairwise conditional ML")

  # position is a fitted covariate, presentation is data: a plain fit and a
  # position fit of the same comparisons are still the same data
  set.seed(291)
  db <- simulate_btl(n_objects = 5, n_judges = 12, reps_per_pair = 3)
  db$ord <- seq_len(nrow(db))
  plain <- btl(db, "object_a", "object_b", winner = "winner", judge = "judge")
  posn <- btl(db, "object_a", "object_b", winner = "winner", judge = "judge",
              position = TRUE, order = "ord")
  expect_true(compare_fits(plain = plain, position = posn)$same_data[2])
  sw <- db; k <- seq(1, nrow(sw), by = 2)
  a <- as.character(sw$object_a); b <- as.character(sw$object_b)
  sw$object_a[k] <- b[k]; sw$object_b[k] <- a[k]
  posn2 <- btl(sw, "object_a", "object_b", winner = "winner",
               judge = "judge", position = TRUE, order = "ord")
  expect_false(compare_fits(pA = posn, pB = posn2)$same_data[2])
})

test_that("curves, styles and selectors follow the design they describe", {
  # a score curve belongs to an administration, not only to a group
  de <- simulate_efrm(n_per_group = 120, items_per_set = 5, n_sets = 2,
                      n_groups = 2, seed = 301)
  ts <- attr(de, "truth")$item_sets
  g1 <- which(de$group == "g1")
  de[g1[1:60], ts[[2]]] <- NA          # half of g1 sat set 1 only
  fe <- rasch_efrm(de, item_sets = ts, groups = "group", id = "id",
                   boot_reps = 0)
  sc <- fe$score_curves
  expect_true("design" %in% names(sc))
  expect_setequal(unique(paste0(sc$group, ":", sc$design)),
                  c("g1:set1", "g1:set1 + set2", "g2:set1 + set2"))
  mx <- tapply(sc$expected_score, paste0(sc$group, ":", sc$design), max)
  expect_lt(mx[["g1:set1"]], mx[["g1:set1 + set2"]])

  # A partial form within one set is its exact scored-item pattern, not the
  # whole set.  These two forms each contribute three items to their curve.
  dp <- simulate_efrm(n_per_group = 120, items_per_set = 5, n_sets = 1,
                      n_groups = 2, seed = 311)
  ip <- attr(dp, "truth")$item_sets[[1]]
  form_a <- rep(c(TRUE, FALSE), length.out = nrow(dp))
  dp[form_a, ip[4:5]] <- NA
  dp[!form_a, ip[1:2]] <- NA
  fp <- rasch_efrm(dp, item_sets = list(core = ip), groups = "group",
                   id = "id", boot_reps = 0)
  sp <- fp$score_curves
  expect_equal(length(unique(paste(sp$group, sp$design))), 4L)
  expect_true(all(grepl("^core \\[", unique(sp$design))))
  expect_equal(unique(sp$n_persons), 60)
  hi <- stats::aggregate(expected_score ~ group + design, sp, max)
  ph <- stats::setNames(fp$phi_table$phi, fp$phi_table$group)
  expect_true(all(hi$expected_score <= 3 * ph[hi$group] + 1e-8))

  # a response style redraws from the probabilities the response was drawn
  # under, so planted local dependence survives it
  d1 <- simulate_rasch(800, 8, model = "PCM", n_categories = 3,
                       dependence = list(pairs = list(c(1, 2)), strength = 3),
                       response_style = list(type = "extreme", prop = 0.5,
                                             strength = 1.6), seed = 302)
  rc <- residual_correlations(rasch(d1, model = "PCM"))$matrix
  expect_gt(rc[1, 2], 0.05)
  d0 <- simulate_rasch(300, 6, model = "PCM", n_categories = 3,
                       response_style = list(type = "extreme", prop = 0.5,
                                             strength = 0), seed = 303)
  expect_false(any(grepl("style", attr(d0, "truth")$planted %||%
                           character(0))))

  # an object set aside at a boundary keeps its predictor row
  set.seed(304)
  db <- simulate_btl(n_objects = 6, n_judges = 20, reps_per_pair = 3)
  ob <- sort(unique(c(as.character(db$object_a), as.character(db$object_b))))
  win <- ob[1]
  db$winner[db$object_a == win | db$object_b == win] <- win
  bp <- data.frame(object = ob, w = seq_along(ob) %% 2)
  expect_no_error(btl_explanatory(db, predictors = bp, formula = ~ w,
                                  object_a = "object_a",
                                  object_b = "object_b", winner = "winner",
                                  judge = "judge"))

  # a repeated identifier does not select a row
  dd <- data.frame(pid = rep(sprintf("P%02d", 1:60), 2),
                   t = rep(1:2, each = 60))
  set.seed(305)
  for (j in 1:6) dd[[paste0("Q", j)]] <- rbinom(120, 1, 0.5)
  st <- stack_data(dd, "pid", "t", paste0("Q", 1:6))
  fs <- rasch(st, id = "id", factors = "time", items = paste0("Q", 1:6))
  png(tf <- tempfile(fileext = ".png"))
  on.exit({dev.off(); unlink(tf)}, add = TRUE)
  expect_error(plot_pcc(fs, "P01"), "appears in 2 rows")
  expect_error(plot_kidmap(fs, "P01"), "appears in 2 rows")
  expect_no_error(plot_pcc(fs, 1))
  expect_no_error(plot_kidmap(fs, 1))

  # the post-hoc identifier is checked where the analysis is described
  set.seed(306)
  X <- matrix(rbinom(400 * 6, 1, plogis(outer(rnorm(400),
    seq(-1, 1, length.out = 6), "-"))), 400, 6)
  colnames(X) <- paste0("I", 1:6)
  d2 <- data.frame(id = sprintf("P%03d", 1:400), X,
                   grp = rep(c("a", "b"), 200))
  f2 <- rasch(d2, id = "id", factors = "grp")
  expect_error(dif_posthoc(f2, "I2", term = "grp",
                           id = sprintf("S%03d", 1:200)),
               "one identifier per row")
  expect_no_error(dif_posthoc(f2, "I2", term = "grp"))
})

test_that("a chained dependence carries no unplanted shift", {
  # 1 -> 3 -> 5: item 3's own carry-over must enter its expectation, or the
  # residual it passes on has a systematic mean and shifts item 5
  d <- simulate_rasch(1500, 6, dependence = list(
    pairs = list(c(1, 3), c(3, 5)), strength = 3), seed = 312)
  rc <- residual_correlations(rasch(d))$matrix
  ctl <- rc[2, 4]
  expect_gt(rc[1, 3] - ctl, 0.08)     # both planted pairs survive the chain
  expect_gt(rc[3, 5] - ctl, 0.08)

  bias5 <- function(pairs, reps = 8) {
    mean(replicate(reps, {
      dd <- simulate_rasch(700, 6, dependence = list(pairs = pairs,
                                                     strength = 3))
      tr <- attr(dd, "truth"); f <- rasch(dd)
      (f$items$location - mean(f$items$location))[5] -
        (tr$difficulty - mean(tr$difficulty))[5]
    }))
  }
  set.seed(321)
  chain <- bias5(list(c(1, 3), c(3, 5)))
  single <- bias5(list(c(3, 5)))
  # the chain must not displace item 5 beyond what its own pair does
  expect_lt(abs(chain - single), 0.3)
})

test_that("names are carried as values, not parsed out of labels", {
  # a set name containing the label separator keeps its items
  de <- simulate_efrm(n_per_group = 120, items_per_set = 5, n_sets = 2,
                      n_groups = 2, seed = 321)
  ts <- attr(de, "truth")$item_sets
  names(ts) <- c("read+write", "maths")
  fe <- rasch_efrm(de, item_sets = ts, groups = "group", id = "id",
                   boot_reps = 0)
  sc <- fe$score_curves
  expect_true(all(grepl("read\\+write", unique(sc$design))))
  # both sets contribute: the curve tops out near the full weighted score
  expect_gt(max(sc$expected_score), 8)

  # a level no item or object carries is not a design column
  set.seed(322)
  X <- matrix(rbinom(500 * 6, 1, plogis(outer(rnorm(500),
    seq(-1, 1, length.out = 6), "-"))), 500, 6)
  colnames(X) <- paste0("I", 1:6)
  prd <- data.frame(item = colnames(X),
                    fmt = factor(rep(c("mc", "open"), 3),
                                 levels = c("mc", "open", "essay")))
  expect_no_error(rasch_explanatory(X, predictors = prd, formula = ~ fmt,
                                    level = "item"))
  db <- simulate_btl(n_objects = 6, n_judges = 20, reps_per_pair = 3)
  ob <- sort(unique(c(as.character(db$object_a), as.character(db$object_b))))
  bp <- data.frame(object = ob,
                   kind = factor(rep(c("a", "b"), 3),
                                 levels = c("a", "b", "never_used")))
  expect_no_error(btl_explanatory(db, predictors = bp, formula = ~ kind,
                                  object_a = "object_a",
                                  object_b = "object_b", winner = "winner",
                                  judge = "judge"))
})

test_that("a resolved comparison name cannot be an existing object", {
  ren <- c(O1 = "A", O2 = "A (g1)", O3 = "C", O4 = "D", O5 = "E")
  set.seed(41)
  d <- simulate_btl(n_objects = 5, n_judges = 60, reps_per_pair = 12,
                    seed = 41)
  for (cc in c("object_a", "object_b", "winner"))
    d[[cc]] <- unname(ren[as.character(d[[cc]])])
  jn <- unique(as.character(d$judge))
  jm <- stats::setNames(rep(c("g1", "g2"), length.out = length(jn)), jn)
  d$grp <- unname(jm[as.character(d$judge)])
  hit <- (d$object_a == "A" | d$object_b == "A") & d$grp == "g1"
  d$winner[hit] <- "A"
  f <- btl(d, "object_a", "object_b", winner = "winner", judge = "judge")
  g <- d$grp[match(paste(f$comparisons$object_a, f$comparisons$object_b,
                         f$comparisons$judge),
                   paste(d$object_a, d$object_b, d$judge))]
  dif <- btl_dif(f, factors = g, min_n = 5)
  expect_true(any(dif$summary$uniform_DIF %in% TRUE))
  # merging the resolved copy into the existing object would bias it
  expect_true(any(grepl("already name object", dif$notes)))
  expect_equal(nrow(f$objects), 5L)
})

test_that("defaults, reserved names and generated names behave", {
  # a default ratio describes the ordinary design, not an explicit request
  d1 <- simulate_efrm(100, 5, n_sets = 1)
  expect_true(any(grepl("ratio 1.00", attr(d1, "truth")$planted)))
  expect_error(simulate_efrm(100, 5, n_sets = 1, set_unit_ratio = 1.3),
               "must be 1 when `n_sets` is 1")
  expect_no_error(simulate_efrm(100, 5, n_groups = 1))
  expect_true(any(grepl("ratio 1.30",
                        attr(simulate_efrm(100, 5, n_sets = 2),
                             "truth")$planted)))

  # a person factor cannot claim a generated column's name
  set.seed(341)
  X <- matrix(rbinom(300 * 5, 1, 0.5), 300, 5)
  colnames(X) <- paste0("I", 1:5)
  expect_error(rasch(X, factors = data.frame(
    class_interval = rep(c("a", "b"), 150))), "reserved")
  expect_error(rasch(X, factors = data.frame(
    theta = rep(c("a", "b"), 150))), "reserved")
  expect_no_error(rasch(X, factors = data.frame(grp = rep(c("a", "b"), 150))))

  # the detail reconciles with the statistic when an item has missing data
  set.seed(342)
  N <- 800; th <- rnorm(N)
  Xm <- matrix(rbinom(N * 8, 1,
    plogis(outer(th, seq(-1.5, 1.5, length.out = 8), "-"))), N, 8)
  colnames(Xm) <- paste0("I", 1:8)
  Xm[1:300, 3] <- NA
  fm <- rasch(Xm)
  cd <- chisq_detail(fm, "I3")
  expect_equal(sum(cd$intervals$chisq[cd$intervals$used]), cd$chisq,
               tolerance = 1e-8)

  # generated and repeated names in the frame design
  de <- simulate_efrm(n_per_group = 100, items_per_set = 5, n_sets = 2,
                      n_groups = 2, seed = 343)
  ts <- attr(de, "truth")$item_sets
  expect_error(rasch_efrm(de, item_sets = list(`(rest)` = ts[[1]]),
                          groups = "group", id = "id", boot_reps = 0),
               "already contains a set named")
  de$sex <- rep(c("m", "f"), length.out = nrow(de))
  expect_error(rasch_efrm(de, item_sets = ts, groups = c("group", "group"),
                          id = "id", boot_reps = 0), "named more than once")
  fe <- rasch_efrm(de, item_sets = ts, groups = "group", id = "id",
                   boot_reps = 0)

  # the report names the estimator the fit used
  rp <- file.path(tempdir(), "gapslbl.html"); unlink(rp)
  invisible(suppressWarnings(report_html(fe, rp)))
  h <- paste(readLines(rp, warn = FALSE), collapse = " ")
  expect_false(grepl("Pairwise conditional estimation", h))
  expect_true(grepl("semiparametric set linking", h, fixed = TRUE))

  # a batch target is one path
  f <- rasch(X)
  expect_error(save_item_plots(f, what = "icc", file = character(0)),
               "one non-missing path")
  expect_error(save_item_plots(f, what = "icc",
                               file = c(file.path(tempdir(), "a.zip"),
                                        file.path(tempdir(), "b.zip"))),
               "one non-missing path")
})

test_that("empty definitions are refused and one-person designs hold", {
  de <- simulate_efrm(n_per_group = 100, items_per_set = 5, n_sets = 2,
                      n_groups = 2, seed = 351)
  ts <- attr(de, "truth")$item_sets
  expect_error(rasch_efrm(de, item_sets = c(ts, list(ghost = character(0))),
                          groups = "group", id = "id", boot_reps = 0),
               "no items")
  expect_error(rasch_efrm(de, item_sets = c(ts, list(ghost = "  ")),
                          groups = "group", id = "id", boot_reps = 0),
               "no items")
  dbe <- simulate_btl_efrm(n_objects_per_set = 4, n_sets = 2,
                           n_judges_per_panel = 8, n_panels = 2,
                           reps_within = 6, reps_cross = 6, seed = 352)
  os <- attr(dbe, "truth")$object_sets
  expect_error(btl_efrm(dbe, "object_a", "object_b", winner = "winner",
                        judge = "judge", panels = "panel",
                        object_sets = c(os, list(ghost = character(0))),
                        se_method = "conditional"), "no objects")

  db <- simulate_btl(n_objects = 5, n_judges = 20, reps_per_pair = 3)
  fb <- btl(db, "object_a", "object_b", winner = "winner", judge = "judge")
  expect_error(btl_dif(fb, factors = list()), "at least one judge factor")

  png(tf <- tempfile(fileext = ".png"))
  on.exit({dev.off(); unlink(tf)}, add = TRUE)
  # one person: the administration pattern is persons x cells, not its
  # transpose, so the design must cover every cell that person answered
  dm <- simulate_mfrm(n_persons = 40, n_items = 4, n_raters = 3, seed = 354)
  one <- dm[dm$person == unique(dm$person)[1], , drop = FALSE]
  fm1 <- rasch_mfrm(one, person = "person", item = "item", score = "score",
                    facets = "rater")
  expect_equal(nrow(fm1$X), 1L)
  blocks <- .design_blocks(fm1)
  expect_equal(length(blocks), 1L)
  expect_equal(length(blocks[[1]]), ncol(fm1$X))
  expect_no_error(test_information(fm1))
  expect_no_error(plot_tif(fm1))

  # one person in a group: their sets must be the sets they answered
  de2 <- simulate_efrm(n_per_group = 150, items_per_set = 4, n_sets = 2,
                       n_groups = 2, seed = 361)
  ts2 <- attr(de2, "truth")$item_sets
  keep <- c(which(de2$group == "g1"), which(de2$group == "g2")[1])
  fe <- rasch_efrm(de2[keep, , drop = FALSE], item_sets = ts2,
                   groups = "group", id = "id", boot_reps = 0)
  ti <- test_information(fe)
  expect_true(all(grepl("set1\\+set2", unique(ti$design))))
  expect_no_error(plot_tif(fe))
})

test_that("a dropped row cannot define the analysis it is dropped from", {
  set.seed(371)
  db <- simulate_btl(n_objects = 5, n_judges = 16, reps_per_pair = 3,
                     model = "polytomous", n_categories = 3)
  db$cnt <- 1L
  base <- btl(db, "object_a", "object_b", response = "response",
              judge = "judge", count = "cnt")
  ghost <- db[1, , drop = FALSE]
  ghost$response <- 5L; ghost$cnt <- 0L
  gone <- btl(rbind(db, ghost), "object_a", "object_b", response = "response",
              judge = "judge", count = "cnt")
  expect_equal(gone$m, base$m)
  expect_true(any(grepl("dropped", gone$notes)))

  # an external person factor is checked exactly as a named column is
  set.seed(375)
  dm <- simulate_mfrm(n_persons = 60, n_items = 4, n_raters = 3, seed = 372)
  v <- data.frame(g = sample(c("a", "b"), nrow(dm), TRUE))
  expect_error(rasch_mfrm(dm, person = "person", item = "item",
                          score = "score", facets = "rater", factors = v),
               "varies within person")
  pu <- unique(as.character(dm$person))
  cst <- data.frame(g = ifelse(match(as.character(dm$person), pu) %% 2 == 0,
                               "a", "b"))
  expect_no_error(rasch_mfrm(dm, person = "person", item = "item",
                             score = "score", facets = "rater",
                             factors = cst))
  expect_no_error(rasch_mfrm(dm, person = "person", item = "item",
                             score = "score", facets = "rater",
                             factors = data.frame(
                               g = rep(c("a", "b"), length.out = length(pu)))))

  # by-value group labels are values, not column names
  de <- simulate_efrm(n_per_group = 100, items_per_set = 5, n_sets = 2,
                      n_groups = 2, seed = 373)
  ts <- attr(de, "truth")$item_sets
  de$cohort <- rep(c("g1", "g2"), length.out = nrow(de))
  fe <- rasch_efrm(de, item_sets = ts, groups = de$group, id = "id",
                   factors = "cohort", items = unlist(ts), boot_reps = 0)
  expect_true("cohort" %in% names(fe$person))

  # every judge needs a stated panel, and every set a real name
  dbe <- simulate_btl_efrm(n_objects_per_set = 4, n_sets = 2,
                           n_judges_per_panel = 8, n_panels = 2,
                           reps_within = 6, reps_cross = 6, seed = 374)
  os <- attr(dbe, "truth")$object_sets
  jn <- unique(as.character(dbe$judge))
  pm <- stats::setNames(rep(c("p1", "p2"), length.out = length(jn)), jn)
  pm[2] <- NA
  expect_error(btl_efrm(dbe, "object_a", "object_b", winner = "winner",
                        judge = "judge", panels = pm, object_sets = os,
                        se_method = "conditional"),
               "missing or blank panel identifier")
  expect_error(rasch_efrm(de, item_sets = stats::setNames(ts, c("   ", "s2")),
                          groups = "group", id = "id", boot_reps = 0),
               "NAMED list")
  mapv <- stats::setNames(rep(c("s1", "s2"), each = 5), unlist(ts))
  mapv[1] <- "  "
  expect_error(rasch_efrm(de, item_sets = mapv, groups = "group", id = "id",
                          boot_reps = 0), "blank set name")
  expect_error(btl_efrm(dbe, "object_a", "object_b", winner = "winner",
                        judge = "judge", panels = "panel",
                        object_sets = stats::setNames(os, c(" ", "s2")),
                        se_method = "conditional"), "named list")
})

test_that("margin categories come from the rows the analysis keeps", {
  set.seed(382)
  db <- simulate_btl(n_objects = 5, n_judges = 20, reps_per_pair = 4)
  db$margin <- factor(sample(c("small", "clear"), nrow(db), TRUE),
                      levels = c("small", "clear"), ordered = TRUE)
  db$cnt <- 1L
  base <- btl(db, "object_a", "object_b", winner = "winner",
              margin = "margin", judge = "judge", count = "cnt")

  # a zero-count row carrying a margin nothing was judged in
  db2 <- db
  db2$margin <- factor(as.character(db2$margin),
                       levels = c("small", "clear", "huge"), ordered = TRUE)
  ghost <- db2[1, , drop = FALSE]
  ghost$margin <- factor("huge", levels = levels(db2$margin), ordered = TRUE)
  ghost$cnt <- 0L
  with_ghost <- btl(rbind(db2, ghost), "object_a", "object_b",
                    winner = "winner", margin = "margin", judge = "judge",
                    count = "cnt")
  expect_equal(with_ghost$m, base$m)
  expect_true(any(grepl("dropped", with_ghost$notes)))

  # a tie only in a dropped row is not a tie in the data
  gtie <- db[1, , drop = FALSE]; gtie$winner <- "tie"; gtie$cnt <- 0L
  with_tie <- btl(rbind(db, gtie), "object_a", "object_b", winner = "winner",
                  margin = "margin", judge = "judge", count = "cnt")
  expect_equal(with_tie$m, base$m)

  # a tie that was actually judged still adds its category
  rtie <- db[1, , drop = FALSE]; rtie$winner <- "tie"; rtie$cnt <- 1L
  with_real <- btl(rbind(db, rtie), "object_a", "object_b", winner = "winner",
                   margin = "margin", judge = "judge", count = "cnt")
  expect_gt(with_real$m, base$m)
})

test_that("ties, keys, identifiers and project files are read faithfully", {
  # a tie carries no margin, so its margin value opens no win/loss category
  set.seed(391)
  db <- simulate_btl(n_objects = 5, n_judges = 20, reps_per_pair = 4)
  db$margin <- factor(sample(c("small", "clear"), nrow(db), TRUE),
                      levels = c("small", "clear", "huge"), ordered = TRUE)
  db$cnt <- 1L
  base <- btl(db, "object_a", "object_b", winner = "winner",
              margin = "margin", judge = "judge", count = "cnt")
  tied <- db
  tied$winner[1] <- "tie"
  tied$margin[1] <- factor("huge", levels = levels(db$margin), ordered = TRUE)
  f2 <- btl(tied, "object_a", "object_b", winner = "winner",
            margin = "margin", judge = "judge", count = "cnt")
  expect_true(isTRUE(f2$converged))
  expect_equal(f2$m, base$m + 1L)          # the tie category, and no more

  # an option nobody can have chosen would score its item zero throughout
  set.seed(393)
  resp <- matrix(sample(c("A", "B", "C"), 400 * 3, TRUE), 400, 3)
  colnames(resp) <- paste0("Q", 1:3)
  expect_error(rasch(resp, key = data.frame(
    item = c("Q1", "Q1", "Q2", "Q3"), option = c("A", "", "A", "C"),
    score = c(1, 1, 1, 1))), "missing or blank option")
  expect_error(rasch(resp, key = data.frame(
    item = c("Q1", "Q2", "Q3"), option = c(NA, "A", "C"),
    score = c(1, 1, 1))), "missing or blank option")
  expect_error(rasch(resp, key = c(Q1 = "A", Q2 = "  ", Q3 = "C")),
               "blank key value")
  expect_no_error(rasch(resp, key = c(Q1 = "A", Q2 = "A", Q3 = "C")))

  # an unknown identifier is not a shared one
  set.seed(392)
  X <- matrix(rbinom(400 * 6, 1, plogis(outer(rnorm(400),
    seq(-1, 1, length.out = 6), "-"))), 400, 6)
  colnames(X) <- paste0("I", 1:6)
  ids <- sprintf("P%03d", 1:400); ids[c(5, 50, 120)] <- NA
  d <- data.frame(id = ids, X, grp = rep(c("a", "b"), 200))
  f <- rasch(d, id = "id", factors = "grp")
  expect_length(dif_anova(f)$within, 0L)
  # the bootstrap treats each unknown identifier as its own person
  bs <- .tailored_boot_rows(c("A", "A", NA, NA, "", "  ", "B"))
  expect_equal(length(unique(bs$id)), 6L)

  db2 <- simulate_btl(n_objects = 5, n_judges = 12, reps_per_pair = 3)
  db2$resp <- NA_integer_
  expect_error(btl(db2, "object_a", "object_b", response = "resp"),
               "no usable comparisons")

  mk <- function(mt) .seal_app_project(list(
    format = "rasch-shiny-project", schema = 2L,
    data = as.data.frame(X), base_fit = f, model_type = mt))
  expect_error(.validate_app_project(mk(c("rasch", "efrm"))),
               "unsupported model type")
  expect_error(.validate_app_project(mk(NA_character_)),
               "unsupported model type")
  expect_no_error(.validate_app_project(mk("rasch")))
})

test_that("unknown identifiers and blank key items do not stand in for data", {
  # a polytomous key row with no item names no item
  set.seed(401)
  resp <- matrix(sample(c("A", "B", "C"), 400 * 3, TRUE), 400, 3)
  colnames(resp) <- paste0("Q", 1:3)
  expect_error(rasch(resp, key = data.frame(
    item = c("Q1", "", "Q3"), option = c("A", "A", "C"),
    score = c(1, 1, 1))), "missing or blank item name")
  expect_error(rasch(resp, key = data.frame(
    item = c("Q1", NA, "Q3"), option = c("A", "A", "C"),
    score = c(1, 1, 1))), "missing or blank item name")
  expect_error(rasch(resp, key = data.frame(
    item = c(1, NaN, 3), option = c("A", "A", "C"),
    score = c(1, 1, 1))), "missing or blank item name")
  expect_error(rasch(resp, key = data.frame(
    item = c("Q1", "Q2", "Q3"), option = c(1, NaN, 3),
    score = c(1, 1, 1))), "missing or blank option")
  expect_equal(nrow(rasch(resp, key = data.frame(
    item = c("Q1", "Q2", "Q3"), option = c("A", "A", "C"),
    score = c(1, 1, 1)))$items), 3L)

  # several unknown identifiers are not one repeated person
  set.seed(402)
  N <- 400
  X <- matrix(rbinom(N * 6, 1, plogis(outer(rnorm(N),
    seq(-1, 1, length.out = 6), "-"))), N, 6)
  colnames(X) <- paste0("I", 1:6)
  ids <- sprintf("P%03d", 1:N); ids[c(5, 50, 120, 300)] <- NA
  d <- data.frame(id = ids, X, grp = rep(c("a", "b"), N / 2))
  ds <- dif_size(rasch(d, id = "id", factors = "grp"), "I2", by = "grp")
  expect_false(any(grepl("repeat across response rows",
                         ds$notes %||% character(0))))
  expect_true(all(is.finite(ds$pairs$se)))
  expect_true(all(is.infinite(ds$pairs$df)))
  # a genuinely repeated design uses person-clustered Wald inference
  dd <- data.frame(pid = rep(sprintf("P%03d", 1:200), 2),
                   t = rep(1:2, each = 200))
  set.seed(403)
  for (j in 1:6) dd[[paste0("Q", j)]] <- rbinom(400, 1, 0.5)
  st <- stack_data(dd, "pid", "t", paste0("Q", 1:6))
  fs <- rasch(st, id = "id", factors = "time", items = paste0("Q", 1:6))
  ds_rep <- dif_size(fs, "Q2", by = "time")
  expect_true(all(is.finite(ds_rep$levels$se)))
  expect_true(all(is.finite(ds_rep$pairs$se)))
  cluster_df <- length(unique(fs$person$id)) - 1L
  expect_equal(ds_rep$pairs$df, cluster_df)
  expect_equal(ds_rep$pairs$p,
               2 * stats::pt(-abs(ds_rep$pairs$t), df = cluster_df))
  expect_equal(ds_rep$pairs$lower,
               ds_rep$pairs$difference -
                 stats::qt(0.975, cluster_df) * ds_rep$pairs$se)
  expect_equal(ds_rep$pairs$p_beyond_A,
               .ets_p_beyond(ds_rep$pairs$difference,
                             ds_rep$pairs$se, cluster_df))
  # Repetition is defined by the rows that contribute conditional pair
  # information, not by duplicate identifiers on unusable rows.
  nonrepeated_support <- fs
  nonrepeated_support$est$cluster_support$repeated <- FALSE
  expect_identical(.dif_refit_df(nonrepeated_support), Inf)
  nonrepeated_support$est$cluster_inference <- FALSE
  expect_true(is.na(.dif_refit_df(nonrepeated_support)))
  legacy <- fs
  legacy$est$cluster_support <- NULL
  expect_true(is.na(.dif_refit_df(legacy)))

  # Duplicate IDs and apparent factor changes confined to empty rows do not
  # create a repeated-person DIF design. This holds for the fitted identifier
  # and for an explicitly supplied analysis identifier.
  set.seed(404)
  Ni <- 180L; Li <- 6L
  Xi <- matrix(rbinom(Ni * Li, 1L, 0.5), Ni, Li,
               dimnames = list(NULL, paste0("D", seq_len(Li))))
  gi <- rep(c("A", "B"), each = Ni / 2L)
  Xi <- rbind(Xi, matrix(NA_integer_, 2L, Li,
                         dimnames = list(NULL, colnames(Xi))))
  idi <- c(sprintf("I%03d", seq_len(Ni)), "empty", "empty")
  fi <- rasch(data.frame(Xi, group = c(gi, "A", "B")),
              id = idi, factors = "group", items = colnames(Xi))
  expect_false(fi$est$cluster_support$repeated)
  dai <- dif_anova(fi, factors = "group", n_groups = 3L)
  expect_length(dai$within, 0L)
  expect_error(dif_anova(fi, factors = "group", n_groups = 3L,
                         within = "group"),
               "within-subject factors need repeated person ids")
  expect_false(dif_contrasts(fi, factors = "group",
                             items = "D2")$paired)
  expect_false(dif_contrasts(fi, factors = "group", items = "D2",
                             id = idi)$paired)
  expect_false(dif_posthoc(fi, "D2", term = "group",
                           factors = "group")$paired)

  # a declared model type is a character scalar, not a factor code
  f0 <- rasch(X)
  mk <- function(mt) .seal_app_project(list(
    format = "rasch-shiny-project", schema = 2L,
    data = as.data.frame(X), base_fit = f0, model_type = mt))
  expect_error(.validate_app_project(mk(factor("rasch"))),
               "unsupported model type")
  expect_no_error(.validate_app_project(mk("rasch")))
})

test_that("one outcome, one column per role, one item name per key entry", {
  set.seed(411)
  db <- simulate_btl(n_objects = 5, n_judges = 16, reps_per_pair = 3,
                     model = "polytomous", n_categories = 3)
  # the polytomous branch never looks at `winner`
  expect_error(btl(db, "object_a", "object_b", winner = "object_a",
                   response = "response", judge = "judge"), "not both")
  expect_no_error(btl(db, "object_a", "object_b", response = "response",
                      judge = "judge"))
  dbw <- simulate_btl(n_objects = 5, n_judges = 16, reps_per_pair = 3)
  expect_no_error(btl(dbw, "object_a", "object_b", winner = "winner",
                      judge = "judge"))
  expect_error(btl(dbw, "object_a", "object_b", winner = "object_a"),
               "role columns must be distinct")
  expect_error(btl(dbw, "object_a", "object_b", winner = "winner",
                   judge = "object_a"), "role columns must be distinct")
  ob <- sort(unique(c(as.character(db$object_a), as.character(db$object_b))))
  bp <- data.frame(object = ob, w = seq_along(ob) %% 2)
  expect_error(btl_explanatory(db, predictors = bp, formula = ~ w,
                               object_a = "object_a", object_b = "object_b",
                               winner = "object_a", response = "response",
                               judge = "judge"), "not both")

  # a role names exactly one column
  expect_error(btl(dbw, c("object_a", "object_b"), "object_b",
                   winner = "winner"), "exactly one column")
  expect_error(btl(dbw, character(0), "object_b", winner = "winner"),
               "exactly one column")
  dm <- simulate_mfrm(n_persons = 40, n_items = 4, n_raters = 3, seed = 412)
  expect_error(rasch_mfrm(dm, person = c("person", "item"), item = "item",
                          score = "score", facets = "rater"),
               "exactly one column")
  expect_error(rasch_mfrm(dm, person = "person", item = "item",
                          score = character(0), facets = "rater"),
               "exactly one column")
  expect_no_error(rasch_mfrm(dm, person = "person", item = "item",
                             score = "score", facets = "rater"))

  dbe <- simulate_btl_efrm(4, 2, 3, 2, reps_within = 4, reps_cross = 4,
                           seed = 414)
  ose <- attr(dbe, "truth")$object_sets
  expect_error(btl_efrm(dbe, "object_a", "object_b", winner = "winner",
                        judge = "object_a", panels = "panel",
                        object_sets = ose, se_method = "conditional"),
               "role columns must be distinct")
  expect_error(btl_efrm(dbe, "object_a", "object_b", winner = "winner",
                        judge = "judge", panels = "judge",
                        object_sets = ose, se_method = "conditional"),
               "role columns must be distinct")

  # every key form needs an item name
  set.seed(413)
  resp <- matrix(sample(c("A", "B", "C"), 400 * 3, TRUE), 400, 3)
  colnames(resp) <- paste0("Q", 1:3)
  expect_error(rasch(resp, key = data.frame(item = c("Q1", "", "Q3"),
                                            key = c("A", "A", "C"))),
               "missing or blank item name")
  expect_error(rasch(resp, key = stats::setNames(c("A", "A", "C"),
                                                 c("Q1", "", "Q3"))),
               "missing or blank item name")
  expect_error(rasch(resp, key = stats::setNames(c("A", "A", "C"),
                                                 c("Q1", NA, "Q3"))),
               "missing or blank item name")
  expect_no_error(rasch(resp, key = c(Q1 = "A", Q2 = "A", Q3 = "C")))

  # the stored schema is read as stored
  X <- matrix(rbinom(300 * 5, 1, 0.5), 300, 5)
  colnames(X) <- paste0("I", 1:5)
  f0 <- rasch(X)
  mk <- function(sc) list(format = "rasch-shiny-project", schema = sc,
                          data = as.data.frame(X), base_fit = f0,
                          model_type = "rasch")
  for (bad in list(factor("future"), 1.5, TRUE, "2", 1L, 3L))
    expect_error(.validate_app_project(mk(bad)), "schema")
  good <- .seal_app_project(mk(2L))
  expect_no_error(.validate_app_project(good))
  expect_no_error(.validate_app_project(.seal_app_project(mk(2))))
})

test_that("banks, sequences and margins are read as the values they hold", {
  # a factor bank column read as level codes would rewrite every location
  set.seed(421)
  X <- matrix(rbinom(400 * 6, 1, plogis(outer(rnorm(400),
    seq(-2, 2, length.out = 6), "-"))), 400, 6)
  colnames(X) <- paste0("I", 1:6)
  f <- rasch(X)
  bank <- data.frame(item = colnames(X),
                     location = c(-2.5, -1.2, 0.25, 0.9, 2.1, 4.0), se = 0)
  bank_f <- bank; bank_f$location <- factor(bank_f$location)
  e1 <- equate_tests(f, bank)
  e2 <- equate_tests(f, bank_f)
  expect_equal(e1$shift, e2$shift)
  expect_equal(e1$table$location_2, e2$table$location_2)

  # the judging sequence has to order the comparisons it describes
  set.seed(423)
  db <- simulate_btl(n_objects = 5, n_judges = 12, reps_per_pair = 3,
                     dependence = list(exposure = 0.3, carry_over = 0.3))
  db$ord <- ave(seq_len(nrow(db)), db$judge, FUN = seq_along)
  fo <- btl(db, "object_a", "object_b", winner = "winner", judge = "judge",
            order = "ord")
  expect_false(is.null(fo$dependence))
  j <- names(which(table(db$judge) >= 2))[1]
  rows <- which(db$judge == j)[1:2]
  tied <- db; tied$ord[rows[2]] <- tied$ord[rows[1]]
  expect_error(btl(tied, "object_a", "object_b", winner = "winner",
                   judge = "judge", order = "ord"), "sequence repeats")
  inf <- db; inf$ord[3] <- Inf
  expect_error(btl(inf, "object_a", "object_b", winner = "winner",
                   judge = "judge", order = "ord"), "finite sequence")
  # the same data in another row order gives the same fit
  sh <- db[sample(nrow(db)), ]
  f2 <- btl(sh, "object_a", "object_b", winner = "winner", judge = "judge",
            order = "ord")
  expect_equal(sort(fo$objects$location), sort(f2$objects$location))

  # a margin is an ordered factor or a positive magnitude on retained wins
  mk <- function(mg) {
    d <- simulate_btl(n_objects = 5, n_judges = 16, reps_per_pair = 4)
    d$margin <- mg[seq_len(nrow(d))]; d
  }
  expect_error(btl(mk(rep(c(TRUE, FALSE), 200)), "object_a", "object_b",
                   winner = "winner", margin = "margin"), "not a margin")
  expect_error(btl(mk(as.Date("2020-01-01") + 0:399), "object_a", "object_b",
                   winner = "winner", margin = "margin"), "not a margin")
  expect_error(btl(mk(rep(c(-1, 2), 200)), "object_a", "object_b",
                   winner = "winner", margin = "margin"),
               "finite and positive")
  expect_error(btl(mk(rep(c(Inf, 2), 200)), "object_a", "object_b",
                   winner = "winner", margin = "margin"),
               "finite and positive")
  zero <- mk(rep(c(1, 2), 200)); zero$margin[1] <- 0
  expect_error(btl(zero, "object_a", "object_b", winner = "winner",
                   margin = "margin"), "zero is a tie")
  # a tie and a zero-count row do not define the analysed margin scale
  ignored_tie <- mk(rep(c(1, 2), 200))
  ignored_tie$winner[1] <- "tie"; ignored_tie$margin[1] <- -1
  expect_no_error(btl(ignored_tie, "object_a", "object_b", winner = "winner",
                      margin = "margin"))
  ignored_count <- mk(rep(c(1, 2), 200))
  ignored_count$n <- 1; ignored_count$n[1] <- 0
  ignored_count$margin[1] <- -1
  expect_no_error(btl(ignored_count, "object_a", "object_b", winner = "winner",
                      margin = "margin", count = "n"))
  expect_no_error(btl(mk(rep(c(1, 2), 200)), "object_a", "object_b",
                      winner = "winner", margin = "margin"))

  # a replication count is real
  dc <- simulate_btl(n_objects = 5, n_judges = 12, reps_per_pair = 3)
  dc$cnt <- complex(real = 2, imaginary = 1)
  expect_error(btl(dc, "object_a", "object_b", winner = "winner",
                   count = "cnt"), "complex")

  # the wide many-facet form checks its selectors too
  dm <- simulate_mfrm(n_persons = 40, n_items = 4, n_raters = 3, seed = 422)
  w <- reshape(dm, idvar = c("person", "rater"), timevar = "item",
               direction = "wide")
  names(w) <- sub("^score[.]", "", names(w))
  its <- setdiff(names(w), c("person", "rater"))
  expect_no_error(rasch_mfrm(w, person = "person", items = its,
                             facets = "rater"))
  expect_error(rasch_mfrm(w, person = c("person", "rater"), items = its,
                          facets = "rater"), "exactly one column")
  expect_error(rasch_mfrm(w, person = character(0), items = its,
                          facets = "rater"), "exactly one column")
  expect_error(rasch_mfrm(w, person = "person", items = character(0),
                          facets = "rater"), "at least one item column")
  expect_error(rasch_mfrm(w, person = "person", items = c(its, "person"),
                          facets = "rater"), "cannot also be the person")
  expect_error(rasch_mfrm(dm, person = "person", item = "item",
                          score = "person", facets = "rater"), "must be distinct")
  expect_error(rasch_mfrm(dm, person = "person", item = "item",
                          score = "score", facets = "rater", factors = "item"),
               "cannot also define another model role")
})

test_that("external calibration and design tables are interpreted strictly", {
  set.seed(431)
  X <- matrix(rbinom(500 * 6, 1, plogis(outer(rnorm(500),
    seq(-1.5, 1.5, length.out = 6), "-"))), 500, 6)
  colnames(X) <- paste0("I", 1:6)
  f <- rasch(X)
  bank <- data.frame(item = f$items$item, location = f$items$location,
                     se = f$items$se, max = f$items$max)

  bad <- bank; bad$location <- as.Date("2020-01-01") + seq_len(nrow(bad))
  expect_error(equate_tests(f, bad), "not a calibration field")
  bad <- bank; bad$location <- complex(real = bad$location, imaginary = 1)
  expect_error(equate_tests(f, bad), "not a calibration field")
  bad <- bank
  bad$location <- I(matrix(rep(bad$location, 2L), nrow = nrow(bad)))
  expect_error(equate_tests(f, bad), "not a calibration field")
  bad <- bank; bad$se <- as.character(bad$se); bad$se[1] <- "unknown"
  expect_error(equate_tests(f, bad), "non-numeric value")
  duplicate_bank <- bank
  names(duplicate_bank)[2] <- "item"
  expect_error(equate_tests(f, duplicate_bank), "column names must be unique")
  duplicate_list <- as.list(bank)
  names(duplicate_list)[2] <- "item"
  expect_error(equate_tests(f, duplicate_list), "column names must be unique")

  # An attached covariance completes individual missing marginal SEs.
  bank_cov <- bank; bank_cov$se <- NA_real_
  attr(bank_cov, "cov_location") <- .equate_loc_cov(f, bank_cov$item)
  eq <- equate_tests(f, bank_cov)
  expect_true(eq$inferential)
  expect_equal(eq$n, nrow(bank_cov))

  # Agreement between stated SEs and a bank covariance is relative to their
  # scale. A 100% discrepancy must not disappear below an absolute floor.
  tiny <- bank
  tiny$se <- rep(1e-8, nrow(tiny))
  attr(tiny, "cov_location") <- diag(rep((2e-8)^2, nrow(tiny)))
  expect_error(equate_tests(f, tiny), "standard errors must agree")
  attr(tiny, "cov_location") <- diag(rep((1e-8)^2, nrow(tiny)))
  expect_no_error(equate_tests(f, tiny))

  indefinite <- bank
  indefinite$se <- rep(1e-10, nrow(indefinite))
  Cbad <- diag(rep(1e-20, nrow(indefinite)))
  Cbad[1, 2] <- Cbad[2, 1] <- 2e-20
  attr(indefinite, "cov_location") <- Cbad
  expect_error(equate_tests(f, indefinite), "positive semidefinite")

  Casym <- diag(rep(1e-20, nrow(indefinite)))
  Casym[1, 2] <- 5e-21
  attr(indefinite, "cov_location") <- Casym
  expect_error(equate_tests(f, indefinite), "symmetric")

  item_predictors <- data.frame(item = colnames(X), x = seq_len(ncol(X)))
  names(item_predictors)[2] <- "item"
  expect_error(.explanatory_metadata(item_predictors, ~ item, X),
               "column names must be unique")
  object_predictors <- data.frame(object = LETTERS[1:5], x = seq_len(5))
  names(object_predictors)[2] <- "object"
  expect_error(.btl_explanatory_design(object_predictors, ~ object,
                                       LETTERS[1:5]),
               "column names must be unique")

  key <- data.frame(item = "I1", option = "A", score = 1)
  names(key)[3] <- "option"
  expect_error(.resolve_key(key), "column names must be unique")
  anchors <- data.frame(item = "I1", k = 1, tau = 0)
  names(anchors)[3] <- "k"
  expect_error(rasch(X, anchors = anchors), "column names must be unique")
})

test_that("paired-comparison banks complete SEs and use two-point precision links", {
  objects <- LETTERS[1:5]
  pairs <- t(utils::combn(objects, 2))
  d <- data.frame(a = rep(pairs[, 1], each = 20),
                  b = rep(pairs[, 2], each = 20))
  d$winner <- ifelse(rep(seq_len(20), nrow(pairs)) <= 10, d$a, d$b)
  fit <- btl(d, "a", "b", winner = "winner")

  bank_cov <- data.frame(object = fit$objects$object,
                         location = fit$objects$location, se = NA_real_)
  attr(bank_cov, "cov_location") <- fit$cov_beta
  eq_cov <- btl_equate(fit, bank_cov)
  expect_true(eq_cov$inferential)
  expect_equal(eq_cov$n_inference, nrow(bank_cov))

  delta <- c(0, 8, 80, 80, 80)
  bank_two <- data.frame(object = fit$objects$object,
                         location = fit$objects$location - delta,
                         se = c(0.05, 1, NA, NA, NA))
  eq_two <- btl_equate(fit, bank_two)
  usable <- is.finite(bank_two$se)
  variance <- fit$objects$se[usable]^2 + bank_two$se[usable]^2
  expected <- weighted.mean(delta[usable], 1 / pmax(variance, 1e-10))
  expect_equal(eq_two$shift, expected)
  expect_false(eq_two$inferential)

  duplicate_bank <- bank_two
  names(duplicate_bank)[2] <- "object"
  expect_error(btl_equate(fit, duplicate_bank),
               "column names must be unique")
  duplicate_list <- as.list(bank_two)
  names(duplicate_list)[2] <- "object"
  expect_error(btl_equate(fit, duplicate_list),
               "column names must be unique")
})

test_that("maps and public controls cannot silently select another analysis", {
  set.seed(441)
  X <- matrix(rbinom(240 * 6, 1, 0.5), 240, 6)
  colnames(X) <- paste0("I", 1:6)
  fac <- rep(c("a", "b"), each = 120)
  fit <- rasch(data.frame(X, fac), factors = "fac")

  blank <- data.frame(fac, check.names = FALSE)
  names(blank) <- ""
  expect_error(dif_anova(fit, factors = blank), "non-empty name")
  expect_error(residual_correlations(fit, flag = factor(1)), "flag")
  expect_error(residual_pca(fit, n_components = factor(1)), "n_components")
  expect_error(dimensionality_test(fit, alpha = factor(0.05),
                                   min_score_points = 2), "alpha")
  expect_error(dimensionality_test(fit, component = factor(1),
                                   min_score_points = 2), "component")
  expect_error(plot_guttman(fit, max_persons = factor(80)), "max_persons")
  expect_error(resolve_dif(fit, min_anchors = factor(3)), "min_anchors")
  expect_error(split_items(fit, c("I1", "I1"), by = "fac"),
               "more than once")
  expect_error(drop_items(fit, c("I1", "I1")), "more than once")
  expect_error(drop_items(fit, matrix("I1")), "at least one non-missing")
  expect_error(dif_posthoc(fit, "I1", term = c("fac", "fac")),
               "more than once")

  de <- simulate_efrm(25, 3, n_sets = 2, n_groups = 2, seed = 442)
  tr <- attr(de, "truth")
  item_map <- setNames(rep(names(tr$item_sets), lengths(tr$item_sets)),
                       unlist(tr$item_sets, use.names = FALSE))
  expect_error(rasch_efrm(de, item_sets = c(item_map, TYPO = "set1"),
                          groups = "group", boot_reps = 0),
               "not in the data")

  db <- simulate_btl_efrm(4, 2, 3, 2, reps_within = 4, reps_cross = 4,
                          seed = 443)
  os <- attr(db, "truth")$object_sets
  panel_map <- setNames(db$panel, db$judge)
  panel_map <- panel_map[!duplicated(names(panel_map))]
  expect_error(btl_efrm(db, "object_a", "object_b", "winner", "judge",
                        c(panel_map, TYPO = "panel1"), os,
                        se_method = "conditional"), "not present")
  dup_sets <- os; names(dup_sets) <- rep("set", length(dup_sets))
  expect_error(btl_efrm(db, "object_a", "object_b", "winner", "judge",
                        "panel", dup_sets, se_method = "conditional"),
               "duplicate set name")
  dup_member <- os; dup_member[[1]] <- c(dup_member[[1]], dup_member[[1]][1])
  expect_error(btl_efrm(db, "object_a", "object_b", "winner", "judge",
                        "panel", dup_member, se_method = "conditional"),
               "repeated within")

  dc <- simulate_btl(5, 16, 6, seed = 444)
  bf <- btl(dc, "object_a", "object_b", "winner", judge = "judge")
  judges <- unique(bf$comparisons$judge)
  group_map <- setNames(rep(c("a", "b"), length.out = length(judges)), judges)
  expect_error(btl_dif(bf, group_map[-1]), "missing from a named factor")
  # a map built from the source data names judges the fit set aside, so a
  # name the fit does not carry is ignored and reported rather than refused.
  # Every fitted judge must still have an entry, and that is what catches a
  # mistyped judge name: it leaves the real judge unmapped.
  extra_ok <- btl_dif(bf, c(group_map, TYPO = "a"))
  expect_true(any(grepl("not in the fitted comparisons", extra_ok$notes)))
  grDevices::pdf(NULL); on.exit(grDevices::dev.off(), add = TRUE)
  expect_error(plot_btl_icc(bf, "O3", group = group_map[-1]),
               "missing from the group map")

  # Numeric NaN is missing metadata, not a judge group literally named NaN.
  numeric_group <- setNames(rep(c(1, 2), length.out = length(judges)), judges)
  numeric_group[1] <- NaN
  nan_dif <- btl_dif(bf, numeric_group)
  used_group <- nan_dif$bootstrap_design$factors[[1L]]
  expect_true(anyNA(used_group))
  expect_false(any(used_group == "NaN", na.rm = TRUE))
  expect_no_error(plot_btl_icc(bf, "O3", group = numeric_group))

  blank_group <- setNames(rep(c("A", "B"), length.out = length(judges)), judges)
  blank_group[1] <- "   "
  blank_dif <- btl_dif(bf, blank_group)
  used_blank <- blank_dif$bootstrap_design$factors[[1L]]
  expect_true(anyNA(used_blank))
  expect_false(any(used_blank == "", na.rm = TRUE))
})

test_that("the person-item map can be restricted to a group or an item set", {
  set.seed(61)
  X <- matrix(rbinom(400 * 8, 1, plogis(outer(rnorm(400),
    seq(-1.5, 1.5, length.out = 8), "-"))), 400, 8)
  colnames(X) <- paste0("I", 1:8)
  d <- data.frame(id = sprintf("P%03d", 1:400), X,
                  grp = rep(c("a", "b"), 200),
                  sex = rep(c("m", "f"), each = 200))
  f <- rasch(d, id = "id", factors = c("grp", "sex"))
  png(tf <- tempfile(fileext = ".png"))
  on.exit({dev.off(); unlink(tf)}, add = TRUE)
  expect_no_error(plot_pimap(f))
  expect_no_error(plot_pimap(f, group = "a"))
  expect_no_error(plot_pimap(f, group = "m"))       # a second factor's level
  expect_no_error(plot_pimap(f, items = c("I1", "I2", "I3")))
  expect_no_error(plot_pimap(f, group = "a", items = c("I1", "I2", "I3")))
  # a selection that names nothing is refused, not silently ignored
  expect_error(plot_pimap(f, group = "zzz"), "not a level of any fitted")
  expect_error(plot_pimap(f, items = c("I1", "NOPE")), "not in the fit")
  expect_error(plot_pimap(f, group = c("a", "b")), "exactly one person-group")
  expect_error(plot_pimap(f, group = matrix("a")),
               "exactly one person-group")

  # an extended-frame fit calibrates item-by-group cells, so a set name must
  # match through the underlying items rather than the virtual keys
  de <- simulate_efrm(n_per_group = 150, items_per_set = 6, n_sets = 2,
                      n_groups = 2, seed = 25)
  fe <- rasch_efrm(de, item_sets = attr(de, "truth")$item_sets,
                   groups = "group", id = "id", boot_reps = 0)
  expect_true(all(c("set1", "set2") %in% unique(as.character(fe$set_of))))
  expect_no_error(plot_pimap(fe, items = "set1"))
  expect_no_error(plot_pimap(fe, group = "g1"))
  expect_no_error(plot_pimap(fe, group = "g1", items = "set2"))
})
test_that("covariance shape checks fail closed", {
  sym <- getFromNamespace(".covariance_is_symmetric", "rasch")
  wald <- getFromNamespace(".covariance_supports_wald", "rasch")

  expect_true(sym(diag(2)))
  expect_true(sym(matrix(numeric(0), 0, 0)))
  expect_false(sym(1:2))
  expect_false(sym(matrix(1:6, 2, 3)))
  expect_false(sym(matrix(c(1, NA, NA, 1), 2)))
  expect_false(sym(diag(2) + 0i))
  expect_false(wald(1:2))
  expect_false(wald(diag(2) + 0i))
})

test_that("shared scalar controls reject shaped and classed values", {
  expect_error(.check_whole(matrix(2), "reps"), "whole number")
  expect_error(.check_whole(structure(2, class = "units"), "reps"),
               "whole number")
  expect_error(.check_prob(matrix(0.05), "alpha"), "probability")
  expect_error(.check_controls(matrix(20), 1e-8), "iteration cap")
  expect_error(.check_controls(20, matrix(1e-8)), "tolerance")
  expect_error(threshold_index(matrix(c(1, 2))), "maximum score")
  expect_error(.check_grid(matrix(seq(-1, 1, length.out = 5))),
               "strictly increasing vector")
  expect_error(.check_band(matrix(2.5)), "half-width")
  expect_error(.sim_count(matrix(10), "n"), "whole number")
  expect_error(.sim_scalar(matrix(1), "effect"), "finite value")
  expect_error(.sim_vector(matrix(1:2), "effects", 2L),
               "plain finite numeric values")
  expect_error(.p_adjust_family(c(0.1, 0.2), n = matrix(2)),
               "whole number")
  expect_false(.app_scalar_text(matrix("rasch")))
})

test_that("DIF design inputs cannot acquire dimensions silently", {
  set.seed(2711)
  X <- matrix(rbinom(240 * 5, 1, 0.5), 240, 5,
              dimnames = list(NULL, paste0("I", 1:5)))
  g <- rep(c("a", "b"), each = 120)
  f <- rasch(data.frame(X, group = g), factors = "group")
  expect_error(dif_anova(f, factors = matrix(g, ncol = 1L)),
               "ordinary vectors")
  expect_error(rasch(data.frame(X, group = g),
                     factors = matrix(g, ncol = 1L)),
               "plain vector")
  expect_error(rasch(data.frame(X), items = matrix(c("I1", "I2"), nrow = 1L)),
               "plain vector")
  expect_error(rasch(data.frame(X), id = matrix(seq_len(nrow(X)), ncol = 1L)),
               "plain vector")
  expect_error(dif_anova(f, within = matrix("group")),
               "ordinary vector")
  expect_error(dif_anova(f, id = matrix(seq_len(240), ncol = 1L)),
               "ordinary vector")
  expect_error(dif_size(f, matrix("I1"), by = "group"),
               "exactly one item")
  expect_error(dif_contrasts(f, items = matrix("I1")),
               "ordinary vector")
  expect_error(dif_posthoc(f, matrix("I1"), term = "group"),
               "one item")
})

test_that("EFRM item roles cannot acquire dimensions silently", {
  d <- simulate_efrm(30, 3, seed = 733)
  tr <- attr(d, "truth")
  expect_error(
    rasch_efrm(d, item_sets = tr$item_sets, groups = "group",
               items = matrix(unlist(tr$item_sets), nrow = 1L)),
    "plain vector")
  expect_error(
    rasch_efrm(d, item_sets = tr$item_sets, groups = "group",
               id = matrix(seq_len(nrow(d)), ncol = 1L)),
    "plain vector")
  shaped_set <- tr$item_sets
  shaped_set[[1L]] <- matrix(shaped_set[[1L]], nrow = 1L)
  expect_error(rasch_efrm(d, item_sets = shaped_set, groups = "group"),
               "plain item-name vectors")
  map <- matrix(rep(c("set1", "set2"), each = 3L), nrow = 1L)
  names(map) <- unlist(tr$item_sets)
  expect_error(rasch_efrm(d, item_sets = map, groups = "group"),
               "item-to-set vector")
  X <- as.matrix(d[unlist(tr$item_sets)])
  expect_error(rasch_efrm(
    X, item_sets = tr$item_sets,
    groups = matrix(d$group, ncol = 1L), boot_reps = 0),
    "plain vector")
  expect_error(rasch_efrm(
    X, item_sets = tr$item_sets, groups = d$group,
    factors = matrix(rep(1:2, length.out = nrow(d)), ncol = 1L),
    boot_reps = 0), "plain vector")
})

test_that("MFRM role selectors cannot acquire dimensions silently", {
  d <- simulate_mfrm(12, 3, 3, seed = 734)
  expect_error(
    rasch_mfrm(d, person = "person", item = "item", score = "score",
               facets = matrix("rater", 1L)),
    "facet")
  expect_error(
    rasch_mfrm(d, person = "person", item = "item", score = "score",
               facets = "rater", interaction = matrix("rater", 1L)),
    "interaction")
  expect_error(
    rasch_mfrm(d, person = "person", item = "item", score = "score",
               facets = "rater", factors = matrix("group", 1L)),
    "plain character")
  expect_error(
    rasch_mfrm(d, person = "person", item = "item", score = "score",
               facets = "rater", factors = matrix(seq_len(nrow(d)), ncol = 1L)),
    "plain character")
})

test_that("by-value Rasch and CJ factors must be plain vectors", {
  X <- matrix(rbinom(80 * 4, 1, 0.5), 80, 4,
              dimnames = list(NULL, paste0("I", 1:4)))
  expect_error(rasch(X, factors = matrix(rep(1:2, each = 40), ncol = 1L)),
               "plain vector")

  d <- simulate_btl(n_objects = 5, n_judges = 12, reps_per_pair = 3,
                    seed = 735)
  fit <- btl(d, "object_a", "object_b", winner = "winner", judge = "judge")
  expect_error(btl_dif(fit, matrix(rep(1:2, length.out = nrow(fit$comparisons)),
                                   ncol = 1L)),
               "plain vector")
})
