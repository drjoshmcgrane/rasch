sim_lltm <- function(n = 600, seed = 4101) {
  set.seed(seed)
  predictors <- data.frame(
    item = paste0("I", 1:10),
    operation = rep(0:1, each = 5),
    format = rep(c("A", "B"), 5),
    stringsAsFactors = FALSE)
  delta <- 0.75 * predictors$operation +
    0.35 * (predictors$format == "B")
  theta <- rnorm(n)
  X <- sapply(delta, function(d) rbinom(n, 1, plogis(theta - d)))
  colnames(X) <- predictors$item
  list(X = X, predictors = predictors, theta = theta, delta = delta)
}

sim_lpcm_explanatory <- function(n = 700, seed = 4102) {
  set.seed(seed)
  predictors <- data.frame(
    item = paste0("P", 1:10),
    format = rep(c("A", "B"), 5),
    stringsAsFactors = FALSE)
  theta <- rnorm(n)
  X <- matrix(0L, n, nrow(predictors))
  for (i in seq_len(ncol(X))) {
    fb <- as.integer(predictors$format[i] == "B")
    tau <- c(-0.7 + 0.25 * fb, 0.65 + 0.55 * fb)
    X[, i] <- vapply(theta, function(b)
      sample.int(3L, 1L, prob = item_moments(b, tau)$P) - 1L, integer(1))
  }
  colnames(X) <- predictors$item
  list(X = X, predictors = predictors)
}

test_that("Kent calibration refuses singular or indefinite uncertainty", {
  C <- diag(2)
  good <- rasch:::.kent_calibration(4, C, 2 * diag(2), diag(2))
  expect_equal(good$lambda, c(2, 2))
  expect_equal(good$chisq, 2)
  expect_equal(good$p, pchisq(2, 2, lower.tail = FALSE))

  singular <- rasch:::.kent_calibration(4, C, diag(2), matrix(0, 2, 2))
  expect_true(is.na(singular$chisq))
  expect_true(is.na(singular$p))
  expect_length(singular$lambda, 0L)

  indefinite <- rasch:::.kent_calibration(4, C, diag(c(1, -1)), diag(2))
  expect_true(is.na(indefinite$chisq))
  expect_true(is.na(indefinite$p))

  asymmetric <- matrix(c(1, 0.5, 0, 1), 2L)
  expect_true(is.na(
    rasch:::.kent_calibration(4, C, asymmetric, diag(2))$chisq))
  expect_true(is.na(
    rasch:::.kent_calibration(4, C, diag(2), asymmetric)$chisq))
})

test_that("failed explanatory calibration withholds covariance inference", {
  d <- sim_lltm(n = 250, seed = 4110)
  M <- ncol(d$X)
  B <- rbind(diag(M - 1L), rep(-1, M - 1L))
  colnames(B) <- paste0("b", seq_len(ncol(B)))
  real_solve <- rasch:::.pcml_solve
  z <- testthat::with_mocked_bindings(
    rasch:::.pcml_design(d$X, B, colnames(B)),
    .pcml_solve = function(...) {
      ans <- real_solve(...)
      ans$converged <- FALSE
      ans
    },
    .package = "rasch")
  expect_false(z$converged)
  expect_true(all(is.na(z$thr$se)))
  expect_true(all(is.na(z$cov_tau)))
  expect_true(all(is.na(z$cov_beta)))
  expect_true(all(is.na(z$coefficients$se)))
  expect_true(all(is.na(z$coefficients$p)))
  expect_true(all(is.na(z$coefficients$p_adj)))
})

test_that("explanatory coefficients use supported person-cluster inference", {
  set.seed(991)
  n_cluster <- 24L
  repeats <- 3L
  predictors <- data.frame(
    item = paste0("I", 1:10),
    operation = rep(0:1, each = 5),
    format = rep(c("A", "B"), 5),
    stringsAsFactors = FALSE)
  theta <- rep(rnorm(n_cluster), each = repeats)
  delta <- -0.5 + 0.7 * predictors$operation +
    0.35 * (predictors$format == "B")
  X <- sapply(delta, function(d)
    rbinom(length(theta), 1, plogis(theta - d)))
  colnames(X) <- predictors$item
  id <- rep(sprintf("P%02d", seq_len(n_cluster)), each = repeats)

  fit <- rasch_explanatory(
    X, predictors, ~ operation + format, id = id)
  cf <- fit$est$coefficients
  expect_true(fit$est$cluster_support$repeated)
  expect_true(fit$est$cluster_inference)
  expect_identical(fit$est$cluster_support$n, n_cluster)
  expect_equal(fit$est$cluster_support$cr1,
               n_cluster / (n_cluster - 1L))
  expect_identical(fit$est$cluster_support$correction,
                   "linearised delete-one-person jackknife")
  expect_equal(cf$df, rep(n_cluster - 1L, nrow(cf)))
  expect_equal(cf$p, 2 * pt(-abs(cf$t), df = n_cluster - 1L),
               tolerance = 1e-12)

  # A fixed-departure refit repeats the calibration with the same independent
  # person clusters and therefore retains the same coefficient reference.
  relaxed <- relax_explanatory(fit, "I1", "location")
  expect_identical(relaxed$est$cluster_support$n, n_cluster)
  expect_true(relaxed$est$cluster_inference)
  expect_identical(relaxed$est$cluster_support$correction,
                   "linearised delete-one-person jackknife")
  expect_equal(relaxed$est$coefficients$df,
               rep(n_cluster - 1L, nrow(relaxed$est$coefficients)))

  unsupported_id <- rep(sprintf("C%d", 1:6), length.out = nrow(X))
  unsupported <- rasch_explanatory(
    X, predictors, ~ operation + format, id = unsupported_id)
  expect_false(unsupported$est$cluster_inference)
  expect_identical(unsupported$est$cluster_support$n, 6L)
  expect_true(all(is.na(unsupported$est$coefficients[c(
    "se", "t", "df", "p", "p_adj")])))
  expect_true(any(grepl("fewer than 10 person clusters",
                        unsupported$est$notes, fixed = TRUE)))

  low_X <- as.matrix(expand.grid(
    I1 = 0:1, I2 = 0:1, I3 = 0:1, I4 = 0:1))[2:9, , drop = FALSE]
  low_predictors <- data.frame(
    item = colnames(low_X), x = c(-1, -0.3, 0.4, 1))
  low <- rasch_explanatory(low_X, low_predictors, ~ x)
  expect_false(low$est$cluster_support$repeated)
  expect_false(low$est$cluster_inference)
  expect_true(all(is.na(low$est$coefficients[c(
    "se", "t", "df", "p", "p_adj")])))
})

test_that("one row per person keeps the finite coefficient reference", {
  # The sandwich rests on the same finite count of independent units whether
  # or not identifiers repeat. Reading single-row units as a limiting normal
  # rejects a true null far more often than the printed probability states.
  set.seed(4471)
  n_person <- 30L
  predictors <- data.frame(
    item = paste0("I", 1:10),
    operation = rep(0:1, each = 5),
    format = rep(c("A", "B"), 5),
    stringsAsFactors = FALSE)
  theta <- rnorm(n_person)
  delta <- -0.5 + 0.7 * predictors$operation +
    0.35 * (predictors$format == "B")
  X <- sapply(delta, function(d)
    rbinom(n_person, 1, plogis(theta - d)))
  colnames(X) <- predictors$item

  fit <- rasch_explanatory(X, predictors, ~ operation + format,
                           id = sprintf("P%02d", seq_len(n_person)))
  cf <- fit$est$coefficients
  n_units <- fit$est$cluster_support$n
  expect_false(fit$est$cluster_support$repeated)
  expect_true(fit$est$cluster_inference)
  expect_equal(cf$df, rep(n_units - 1L, nrow(cf)))
  expect_equal(cf$p, 2 * pt(-abs(cf$t), df = n_units - 1L),
               tolerance = 1e-12)
  # strictly larger than the limiting-normal probabilities it replaces
  expect_true(all(cf$p > 2 * pnorm(-abs(cf$t))))
})

test_that("explanatory CJ accepts predictors for set-aside boundary objects", {
  pr <- t(utils::combn(LETTERS[1:4], 2L))
  core <- do.call(rbind, lapply(seq_len(nrow(pr)), function(i)
    data.frame(a = rep(pr[i, 1L], 2L), b = rep(pr[i, 2L], 2L),
               winner = pr[i, ])))
  # E and F occur in the observed comparisons but form an isolated
  # undefeated/winless pair.  Neither can be extrapolated against the fitted
  # core, so neither appears in the reported object table.
  d <- rbind(core,
             data.frame(a = rep("E", 4L), b = rep("F", 4L),
                        winner = rep("E", 4L)))
  predictors <- data.frame(object = LETTERS[1:6], x = seq_len(6L))
  fit <- btl_explanatory(d, predictors, ~ x, "a", "b", winner = "winner")
  expect_s3_class(fit, "rasch_btl_explanatory")
  expect_setequal(fit$objects$object, LETTERS[1:4])
  expect_true(all(c("E", "F") %in%
    c(fit$observed_comparisons$object_a,
      fit$observed_comparisons$object_b)))
})

test_that("explanatory CJ rejects predictors seen only in unusable rows", {
  pairs <- rbind(c("A", "B"), c("A", "D"), c("B", "D"))
  d <- do.call(rbind, lapply(seq_len(nrow(pairs)), function(i)
    data.frame(a = rep(pairs[i, 1L], 20), b = rep(pairs[i, 2L], 20),
               winner = rep(pairs[i, ], 10))))
  d <- rbind(d, data.frame(a = "C", b = "A", winner = NA_character_))
  predictors <- data.frame(object = c("A", "B", "C", "D"), x = 1:4)

  expect_error(
    btl_explanatory(d, predictors, ~ x, "a", "b", winner = "winner"),
    "not present in the comparisons")
})

test_that("LLTM recovers item-feature effects and retains Rasch scoring", {
  d <- sim_lltm()
  f <- rasch_explanatory(d$X, d$predictors,
                         ~ operation + format, level = "item")

  expect_s3_class(f, "rasch_explanatory")
  expect_identical(f$explanatory_model, "LLTM")
  expect_true(f$est$converged)
  expect_equal(f$est$n_parameters, 2L)
  expect_lt(abs(f$est$coefficients["operation", "estimate"] - 0.75), .25)
  expect_lt(abs(f$est$coefficients["formatB", "estimate"] - 0.35), .25)
  # independent rows are still a finite number of sampling units
  n_units <- f$est$cluster_support$n
  expect_equal(f$est$coefficients$df, rep(n_units - 1L, 2L))
  expect_equal(f$est$coefficients$p,
               2 * pt(-abs(f$est$coefficients$t), df = n_units - 1L),
               tolerance = 1e-12)
  expect_equal(f$person$theta, f$score_table$theta[f$person$raw + 1L],
               tolerance = 1e-10)

  z <- explanatory_test(f)
  expect_equal(z$df, f$reference_fit$est$n_parameters - f$est$n_parameters)
  expect_true(is.finite(z$p_kent))
  expect_identical(z$p, z$p_kent)
  expect_true(is.finite(z$p_naive))
  free_tau <- f$reference_fit$est$thr$tau
  active_tau <- f$est$thr$tau
  departure <- free_tau - active_tau
  expected_r2 <- 1 - sum((departure - mean(departure))^2) /
    sum((free_tau - mean(free_tau))^2)
  expect_equal(z$r_squared, expected_r2)
  full_rank <- qr(cbind(1, f$est$B))$rank
  expect_equal(z$r_squared_adj,
               1 - (1 - z$r_squared) * (length(free_tau) - 1) /
                 (length(free_tau) - full_rank))
  expect_lte(z$r_squared_adj, z$r_squared)
  expect_identical(z$r2_basis, "threshold calibration")
  expect_lte(z$r_squared, 1)
  weak_fit <- f
  weak_fit$reference_fit$est$thr$weak[1] <- TRUE
  weak_fit$est$thr$weak[1] <- TRUE
  weak_z <- explanatory_test(weak_fit)
  keep <- seq_along(free_tau) != 1L
  weak_departure <- free_tau[keep] - active_tau[keep]
  weak_expected <- 1 -
    sum((weak_departure - mean(weak_departure))^2) /
    sum((free_tau[keep] - mean(free_tau[keep]))^2)
  expect_equal(weak_z$r_squared, weak_expected)
  weak_ok <- !weak_fit$reference_fit$est$thr$weak
  weak_n <- sum(weak_ok)
  weak_df <- weak_n - qr(cbind(1, weak_fit$est$B[weak_ok, , drop = FALSE]))$rank
  expect_equal(weak_z$r_squared_adj,
               1 - (1 - weak_z$r_squared) * (weak_n - 1) / weak_df)

  # an exclusion that removes a level's only support reduces the retained
  # design's rank, so the residual dimension comes from that rank rather
  # than from the full comparison's parameter count
  rare <- d$predictors
  rare$format <- factor(c("C", rep(c("A", "B"),
                                   length.out = nrow(rare) - 1L)))
  rf <- rasch_explanatory(d$X, rare, ~ operation + format, level = "item")
  wr <- rf
  wr$reference_fit$est$thr$weak[1] <- TRUE     # the only "C" item
  wz <- explanatory_test(wr)
  ok_r <- !wr$reference_fit$est$thr$weak
  rank_r <- qr(cbind(1, wr$est$B[ok_r, , drop = FALSE]))$rank
  expect_lt(rank_r, qr(cbind(1, wr$est$B))$rank)   # rank genuinely fell

  # with every calibration excluded there is nothing to adjust: the
  # coefficient is unavailable silently rather than after a warning
  all_weak <- rf
  all_weak$reference_fit$est$thr$weak[] <- TRUE
  expect_no_warning(all_z <- explanatory_test(all_weak))
  expect_true(is.na(all_z$r_squared_adj))
  expect_equal(wz$r_squared_adj,
               1 - (1 - wz$r_squared) * (sum(ok_r) - 1) /
                 (sum(ok_r) - rank_r))
})

test_that("LPCM accepts threshold effects and selected interactions", {
  d <- sim_lpcm_explanatory()
  f <- rasch_explanatory(d$X, d$predictors,
    ~ format + threshold + format:threshold, level = "item")

  expect_identical(f$explanatory_model, "LPCM")
  expect_true(f$est$converged)
  expect_true(all(c("formatB", "threshold2", "formatB:threshold2") %in%
                    f$est$coefficients$term))
  expect_true(all(is.finite(f$est$coefficients$se)))
  expect_equal(nrow(f$explanatory$metadata), sum(f$m))
})

test_that("a partly represented threshold departure adds only new directions", {
  d <- simulate_rasch(500, 6, model = "PCM", n_categories = 4, seed = 91)
  items <- sprintf("I%02d", 1:6)
  X <- as.data.frame(d[items])
  predictors <- expand.grid(
    threshold = 1:3, item = items, KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE)
  predictors <- predictors[
    order(match(predictors$item, items), predictors$threshold), ]
  predictors$x <- 0
  predictors$x[predictors$item == "I01"] <- c(1, 0, -1)
  predictors$trend <- rep(seq_along(items), each = 3)
  fit <- rasch_explanatory(
    X, predictors, ~ x + trend, level = "threshold")

  full <- rasch:::.explanatory_candidate(fit, "I01", "thresholds")
  addition <- rasch:::.explanatory_addition(fit$est$B, full)
  expect_equal(ncol(full), 2L)
  expect_equal(ncol(addition), 1L)
  expect_equal(qr(cbind(fit$est$B, addition), tol = 1e-10)$rank,
               ncol(fit$est$B) + 1L)

  diagnostics <- explanatory_diagnostics(fit)
  row <- diagnostics$item == "I01" &
    diagnostics$component == "Threshold structure"
  expect_true(any(row))
  expect_equal(diagnostics$parameters_added[row], 1L)
  relaxed <- relax_explanatory(fit, "I01", "thresholds")
  expect_s3_class(relaxed, "rasch_explanatory")
  expect_equal(tail(relaxed$explanatory$relaxations$parameters_added, 1L), 1L)
})

test_that("large-sample explanatory convergence is judged on parameter scale", {
  # This seed leaves an absolute projected score of about 2.1e-4 although
  # the remaining Newton move is only 2.1e-8 and the estimates are stable.
  set.seed(8400072L)
  x <- rep(c(-1, -0.5, 0, 0.5, 1), 2L)
  predictors <- data.frame(item = paste0("I", 1:10), x = x)
  theta <- rnorm(2000)
  X <- sapply(seq_along(x), function(i) {
    tau <- c(-0.8, 0.2, 1) + 0.45 * x[i]
    vapply(theta, function(th)
      sample.int(4L, 1L, prob = item_moments(th, tau)$P) - 1L, integer(1))
  })
  storage.mode(X) <- "integer"
  colnames(X) <- predictors$item

  fit <- rasch_explanatory(X, predictors, ~ x + threshold,
                           maxit = 80L, tol = 1e-8)
  expect_true(fit$est$converged)
  expect_true(fit$reference_fit$est$converged)
  expect_equal(fit$est$coefficients["x", "estimate"], 0.4452621,
               tolerance = 1e-6)
})

test_that("a stable small-sample LPCM fit is not rejected at a score boundary", {
  # The remaining Newton move is about 1.0e-7 although the extensive score
  # lies just above its absolute cutoff.
  set.seed(8200119L)
  x <- rep(c(-1, -0.3, 0.3, 1), 3L)
  predictors <- data.frame(item = paste0("I", 1:12), x = x)
  theta <- rnorm(300)
  X <- sapply(seq_along(x), function(i) {
    tau <- c(-0.8, 0.2, 1) + 0.45 * x[i]
    vapply(theta, function(th)
      sample.int(4L, 1L, prob = item_moments(th, tau)$P) - 1L, integer(1))
  })
  storage.mode(X) <- "integer"
  colnames(X) <- predictors$item

  fit <- rasch_explanatory(X, predictors, ~ x + threshold,
                           maxit = 80L, tol = 1e-8)
  expect_true(fit$est$converged)
  expect_true(fit$reference_fit$est$converged)
})

test_that("Rasch explanatory predictors may be continuous, categorical or ordinal", {
  set.seed(4104)
  p <- expand.grid(category = c("A", "B"),
                   ordinal = c("low", "moderate", "high"),
                   replicate = 1:3, stringsAsFactors = FALSE)
  p$item <- paste0("T", seq_len(nrow(p)))
  p$continuous <- scale(rnorm(nrow(p)))[, 1]
  p$ordinal <- ordered(p$ordinal,
                       levels = c("low", "moderate", "high"))
  delta <- .35 * p$continuous + .4 * (p$category == "B") +
    .3 * (as.integer(p$ordinal) - 2)
  theta <- rnorm(700)
  X <- sapply(delta, function(d) rbinom(length(theta), 1,
                                         plogis(theta - d)))
  colnames(X) <- p$item
  fit <- rasch_explanatory(
    X, p, ~ continuous + category + ordinal + category:ordinal)

  terms <- fit$est$coefficients$term
  expect_true("continuous" %in% terms)
  expect_true("categoryB" %in% terms)
  ordinal_terms <- terms[grepl("^ordinal", terms) & !grepl(":", terms)]
  expect_equal(length(ordinal_terms), 2L)
  expect_false(any(grepl("ordinal\\.[LQCTest]", terms)))
  expect_equal(sum(grepl("ordinal", terms) & grepl("categoryB", terms)), 2L)
  expect_true(all(abs(fit$est$coefficients[ordinal_terms, "estimate"] - .3) < .2))
  expect_true(is.ordered(fit$explanatory$metadata$ordinal))
})

test_that("ordinal contrast names remain unique when level labels collide", {
  set.seed(4110)
  lev <- c("a b", "a-b", "a.b", "c", "d")
  ord <- ordered(lev, levels = lev)
  X <- matrix(rbinom(1000, 1, .5), 200, 5,
              dimnames = list(NULL, paste0("I", 1:5)))
  p <- data.frame(item = colnames(X), ord = ord, check.names = FALSE)
  rf <- rasch_explanatory(X, p, ~ ord)
  expect_s3_class(rf, "rasch_explanatory")
  expect_identical(anyDuplicated(rf$est$coefficients$term), 0L)
  expect_true(all(grepl("adjacent_[1-4]", rf$est$coefficients$term)))

  d <- simulate_btl(5, 20, reps_per_pair = 8, seed = 4111)
  bp <- data.frame(object = paste0("O", 1:5), ord = ord,
                   check.names = FALSE)
  bf <- btl_explanatory(
    d, bp, ~ ord, "object_a", "object_b", winner = "winner",
    judge = "judge")
  expect_s3_class(bf, "rasch_btl_explanatory")
  expect_identical(anyDuplicated(bf$object_coefficients$term), 0L)
  expect_true(all(grepl("adjacent_[1-4]", bf$object_coefficients$term)))
})

test_that("threshold metadata and unidentified formulae are checked", {
  d <- sim_lpcm_explanatory(n = 350)
  p <- merge(d$predictors,
             data.frame(threshold = 1:2), by = NULL)
  p <- p[order(match(p$item, d$predictors$item), p$threshold), ]
  p$scoring <- p$threshold * ifelse(p$format == "B", 2, 1)
  f <- rasch_explanatory(d$X, p, ~ format + scoring,
                         level = "threshold")
  expect_s3_class(f, "rasch_explanatory")

  expect_error(rasch_explanatory(d$X, p[-1, ], ~ scoring,
                                 level = "threshold"), "missing")
  q <- d$predictors; q$constant <- 1
  expect_error(rasch_explanatory(d$X, q, ~ constant),
               "no estimable variation")
  q$bad <- seq_len(nrow(q)); q$bad[1] <- Inf
  expect_error(rasch_explanatory(d$X, q, ~ bad), "non-finite")
})

test_that("diagnostics use one Holm family and fixed departures refit all outputs", {
  d <- sim_lltm(n = 500)
  f <- rasch_explanatory(d$X, d$predictors, ~ operation + format)
  dg <- explanatory_diagnostics(f)
  expect_equal(dg$p_adj, p.adjust(dg$p, "holm"))
  expect_equal(sort(unique(dg$component)), "Item location")

  before <- f$person$theta
  g <- relax_explanatory(f, paste0(" ", dg$item[1], " "), "location")
  expect_equal(nrow(g$explanatory$relaxations), 1L)
  expect_gt(g$est$n_parameters, f$est$n_parameters)
  expect_false(isTRUE(all.equal(before, g$person$theta)))
  expect_equal(g$residuals,
               (g$X - g$moments$E) / sqrt(g$moments$V),
               tolerance = 1e-10)
  expect_error(relax_explanatory(f, data.frame(item = dg$item[1])),
               "exactly one item")
  failed <- f
  failed$est$converged <- FALSE
  expect_output(print(failed), "comparison: unavailable")
  expect_error(explanatory_test(failed), "did not converge")
  expect_error(explanatory_diagnostics(failed), "did not converge")
  expect_error(relax_explanatory(failed, dg$item[1]), "did not converge")
  bad_reference <- f
  bad_reference$reference_fit$est$converged <- FALSE
  expect_error(explanatory_test(bad_reference), "reference fit did not converge")
})

test_that("non-convergent explanatory candidates are withheld, not dropped", {
  d <- sim_lltm(n = 300)
  f <- rasch_explanatory(d$X, d$predictors, ~ operation + format)
  real_design <- rasch:::.pcml_design
  testthat::local_mocked_bindings(
    .pcml_design = function(...) {
      out <- real_design(...)
      out$converged <- FALSE
      out
    },
    .package = "rasch")
  dg <- explanatory_diagnostics(f)
  expect_gt(nrow(dg), 0L)
  expect_true(all(!dg$converged))
  expect_true(all(is.na(dg$departure)))
  expect_true(all(is.na(dg$p)))
  expect_true(all(is.na(dg$p_adj)))
  expect_match(attr(dg, "note"), "non-convergent")
  expect_error(relax_explanatory(f, dg$item[1L]),
               "relaxed explanatory calibration did not converge")
})

test_that("item changes preserve and refit explanatory calibrations", {
  d <- sim_lltm(n = 450)
  groups <- data.frame(group = rep(c("A", "B"), length.out = nrow(d$X)))
  f <- rasch_explanatory(d$X, d$predictors, ~ operation + format,
                         factors = groups)

  dropped <- drop_items(f, "I10")
  expect_s3_class(dropped, "rasch_explanatory")
  expect_equal(ncol(dropped$X), 9L)

  split <- split_items(f, "I1", by = "group")
  expect_s3_class(split, "rasch_explanatory")
  expect_equal(ncol(split$X), 11L)
  expect_true(any(split$explanatory$relaxations$component == "Item location"))
  expect_false(isTRUE(all.equal(f$person$theta, split$person$theta)))

  combined <- combine_items(f, list(c("I1", "I2")))
  expect_s3_class(combined, "rasch_explanatory")
  expect_equal(ncol(combined$X), 9L)
  expect_true(nrow(combined$explanatory$relaxations) >= 1L)
  expect_false(isTRUE(all.equal(f$person$theta, combined$person$theta)))

  dependence <- dependence_magnitude(f, dependent = "I2", independent = "I1")
  expect_s3_class(dependence$refit, "rasch_explanatory")
  expect_true(all(grepl("I2\\|I1=", tail(dependence$refit$items$item, 2L))))
  expect_false(isTRUE(all.equal(f$person$theta, dependence$refit$person$theta)))
})

test_that("multiple-choice data remain scored and inspectable after item changes", {
  set.seed(4106)
  n <- 450
  pred <- data.frame(item = paste0("M", 1:10),
                     feature = rep(0:1, each = 5))
  theta <- rnorm(n)
  raw <- sapply(.65 * pred$feature, function(d) {
    correct <- rbinom(n, 1, plogis(theta - d)) == 1L
    ifelse(correct, "A", sample(c("B", "C", "D"), n, replace = TRUE))
  })
  colnames(raw) <- pred$item
  group <- factor(rep(c("A", "B"), length.out = n))
  fit <- rasch_explanatory(
    raw, pred, ~ feature, factors = data.frame(group = group),
    key = stats::setNames(rep("A", ncol(raw)), colnames(raw)))

  split <- split_items(fit, "M1", by = "group")
  expect_type(split$X, "integer")
  expect_true(all(c("M1 (A)", "M1 (B)") %in% colnames(split$mc$raw)))
  expect_true(all(is.na(split$mc$raw[group == "B", "M1 (A)"])))
  expect_true(all(is.na(split$mc$raw[group == "A", "M1 (B)"])))
  expect_s3_class(distractor_analysis(split, "M1 (A)"), "data.frame")

  relaxed <- relax_explanatory(split, "M2", "location")
  expect_identical(relaxed$mc$raw, split$mc$raw)

  dropped <- drop_items(fit, "M10")
  expect_equal(ncol(dropped$mc$raw), 9L)
  combined <- combine_items(fit, list(c("M1", "M2")))
  expect_equal(ncol(combined$mc$raw), 8L)
  expect_true(all(colnames(combined$mc$raw) %in% paste0("M", 3:10)))
})

sim_explanatory_cj <- function(reps = 30, polytomous = FALSE, seed = 4103) {
  set.seed(seed)
  predictors <- data.frame(
    object = LETTERS[1:8],
    domain = rep(0:1, each = 4),
    format = rep(c("A", "B"), 4), stringsAsFactors = FALSE)
  beta <- 0.8 * predictors$domain + 0.35 * (predictors$format == "B")
  names(beta) <- predictors$object
  pairs <- t(combn(predictors$object, 2L))
  d <- data.frame(a = rep(pairs[, 1], each = reps),
                  b = rep(pairs[, 2], each = reps),
                  stringsAsFactors = FALSE)
  if (!polytomous) {
    p <- plogis(beta[d$a] - beta[d$b])
    d$winner <- ifelse(runif(nrow(d)) < p, d$a, d$b)
  } else {
    tau <- c(-0.8, 0.8)
    d$response <- vapply(seq_len(nrow(d)), function(r)
      sample.int(3L, 1L,
        prob = item_moments(beta[d$a[r]] - beta[d$b[r]], tau)$P) - 1L,
      integer(1))
  }
  list(data = d, predictors = predictors)
}

test_that("explanatory CJ recovers object effects and compares with a free fit", {
  d <- sim_explanatory_cj(reps = 35)
  f <- btl_explanatory(d$data, d$predictors, ~ domain + format,
                       object_a = "a", object_b = "b", winner = "winner")
  expect_s3_class(f, "rasch_btl_explanatory")
  expect_true(f$converged)
  expect_lt(abs(f$object_coefficients["domain", "estimate"] - .8), .3)
  expect_lt(abs(f$object_coefficients["formatB", "estimate"] - .35), .3)
  z <- explanatory_test(f)
  expect_equal(z$df, 5L)
  expect_true(is.finite(z$p_kent))
  expect_identical(z$p, z$p_kent)
  expect_true(is.finite(z$p_naive))
  free_location <- f$reference_fit$objects$location[
    match(f$objects$object, f$reference_fit$objects$object)]
  departure <- free_location - f$objects$location
  expected_r2 <- 1 - sum((departure - mean(departure))^2) /
    sum((free_location - mean(free_location))^2)
  expect_equal(z$r_squared, expected_r2)
  expect_identical(z$r2_basis, "object calibration")
  expect_lte(z$r_squared, 1)

  dg <- explanatory_diagnostics(f)
  expect_equal(dg$p_adj, p.adjust(dg$p, "holm"))
  g <- relax_btl_explanatory(f, paste0(" ", dg$object[1], " "))
  expect_s3_class(g, "rasch_btl_explanatory")
  expect_equal(nrow(g$explanatory$relaxations), 1L)
  expect_false(isTRUE(all.equal(f$objects$location, g$objects$location)))
  expect_error(relax_btl_explanatory(
    f, data.frame(object = dg$object[1])), "exactly one object")
  failed <- f
  failed$converged <- FALSE
  expect_error(explanatory_test(failed), "did not converge")
  expect_error(explanatory_diagnostics(failed), "did not converge")
  expect_error(relax_btl_explanatory(failed, dg$object[1]), "did not converge")
  bad_reference <- f
  bad_reference$reference_fit$converged <- FALSE
  expect_error(explanatory_test(bad_reference), "reference fit did not converge")

  real_refit <- rasch:::.btl_explanatory_refit
  testthat::local_mocked_bindings(
    .btl_explanatory_refit = function(...) {
      out <- real_refit(...)
      out$converged <- FALSE
      out
    },
    .package = "rasch")
  expect_error(relax_btl_explanatory(f, dg$object[1]),
               "relaxed explanatory comparative judgement calibration did not converge")
})

test_that("ordered explanatory CJ retains its response-threshold model", {
  d <- sim_explanatory_cj(reps = 45, polytomous = TRUE)
  f <- btl_explanatory(d$data, d$predictors,
                       ~ domain + format + domain:format,
                       object_a = "a", object_b = "b",
                       response = "response", thresholds = "pc")
  expect_s3_class(f, "rasch_btl_explanatory")
  expect_equal(f$m, 2L)
  expect_identical(f$thr_structure, "pc")
  expect_true("domain:formatB" %in% f$object_coefficients$term)
  expect_equal(f$thresholds$tau, f$reference_fit$thresholds$tau,
               tolerance = .5)
})

test_that("comparative judgement predictors retain all three variable types", {
  set.seed(4105)
  p <- expand.grid(category = c("A", "B"),
                   ordinal = c("low", "moderate", "high"),
                   replicate = 1:2, stringsAsFactors = FALSE)
  p$object <- paste0("O", seq_len(nrow(p)))
  p$continuous <- scale(rnorm(nrow(p)))[, 1]
  p$ordinal <- ordered(p$ordinal,
                       levels = c("low", "moderate", "high"))
  beta <- .3 * p$continuous + .45 * (p$category == "B") +
    .25 * (as.integer(p$ordinal) - 2)
  names(beta) <- p$object
  pairs <- t(combn(p$object, 2L))
  d <- data.frame(a = rep(pairs[, 1], each = 18),
                  b = rep(pairs[, 2], each = 18))
  prob <- plogis(beta[d$a] - beta[d$b])
  d$winner <- ifelse(runif(nrow(d)) < prob, d$a, d$b)
  fit <- btl_explanatory(
    d, p, ~ continuous + category + ordinal + category:ordinal,
    "a", "b", winner = "winner")

  terms <- fit$object_coefficients$term
  expect_true(all(c("continuous", "categoryB") %in% terms))
  ordinal_terms <- terms[grepl("^ordinal", terms) & !grepl(":", terms)]
  expect_equal(length(ordinal_terms), 2L)
  expect_false(any(grepl("ordinal\\.[LQCTest]", terms)))
  expect_equal(sum(grepl("ordinal", terms) & grepl("categoryB", terms)), 2L)
  expect_true(all(abs(fit$object_coefficients[ordinal_terms, "estimate"] - .25) < .2))
  expect_true(is.ordered(fit$explanatory$metadata$ordinal))
})

test_that("saturated explanatory designs reproduce free calibrations", {
  d <- sim_lltm(n = 450)
  fr <- rasch(d$X)
  fe <- rasch_explanatory(d$X, d$predictors, ~ item)
  expect_equal(fe$thresholds$tau, fr$thresholds$tau, tolerance = 1e-7)
  expect_equal(explanatory_test(fe)$df, 0L)

  cj <- sim_explanatory_cj(reps = 30)
  br <- btl(cj$data, "a", "b", winner = "winner")
  be <- btl_explanatory(cj$data, cj$predictors, ~ object,
                        "a", "b", winner = "winner")
  expect_equal(be$objects$location, br$objects$location, tolerance = 1e-7)
  expect_equal(explanatory_test(be)$df, 0L)
})
