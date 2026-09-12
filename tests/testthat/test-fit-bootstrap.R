# The parametric bootstrap null for the item fit statistics: that the
# replicates travel the same road as the observed data, that each statistic
# is read in the tail or tails that mean misfit, that the p-values have the
# resolution their replicate count allows, and that the correction runs in
# the direction the miscalibration requires.

# the resolution-floor warning is correct at these deliberately small B --
# assert it once below, and keep every other call quiet
fbq <- function(...) suppressWarnings(fit_bootstrap(...))

simb <- function(N, d, seed = 1) {
  set.seed(seed)
  th <- rnorm(N)
  X <- matrix(rbinom(N * length(d), 1, plogis(outer(th, d, "-"))),
              N, length(d))
  colnames(X) <- paste0("I", seq_along(d))
  X
}

STATS <- c("chisq", "fit_resid", "infit_ms", "outfit_ms", "infit_z", "outfit_z")

test_that("every fit statistic is carried with a probability in range", {
  fit <- rasch(simb(300, seq(-1.5, 1.5, length.out = 6), seed = 3))
  bs <- fbq(fit, B = 29, seed = 1)
  expect_s3_class(bs, "rasch_fit_bootstrap")
  expect_identical(bs$model_kind, "rasch")
  expect_equal(nrow(bs$items), 6L)
  expect_identical(bs$items$item, fit$items$item)
  for (st in STATS) {
    expect_equal(bs$items[[st]], fit$items[[st]], info = st)
    p <- bs$items[[paste0(st, "_p_boot")]]
    expect_true(all(p > 0 & p <= 1), info = st)
    expect_true(all(p >= 1 / (1 + bs$B_used)), info = st)
    expect_true(all(bs$items[[paste0(st, "_p_boot_adj")]] >= p - 1e-12), info = st)
  }
  expect_equal(bs$B, 29L)
  expect_true(bs$B_used <= bs$B)
  expect_identical(bs$theta, "conditional")
  expect_setequal(names(bs$replicates), STATS)
  expect_true(all(vapply(bs$replicates, nrow, 0L) == bs$B_used))
  expect_equal(nrow(bs$persons), nrow(fit$person))
  expect_identical(bs$persons$id, fit$person$id)
  for (st in c("fit_resid", "infit_ms", "outfit_ms", "infit_z", "outfit_z")) {
    p <- bs$persons[[paste0(st, "_p_boot")]]
    adj <- bs$persons[[paste0(st, "_p_boot_adj")]]
    expect_true(all(p > 0 & p <= 1, na.rm = TRUE), info = st)
    expect_true(all(adj >= p, na.rm = TRUE), info = st)
  }
  expect_identical(bs$person_adjustment$method,
                   "single-step maximum-statistic bootstrap")
  shown <- capture.output(returned <- print(bs))
  expect_identical(returned, bs)
  expect_lte(length(shown), 5L)
  expect_true(any(grepl("Usable replicates", shown)))
  expect_false(any(grepl("replicates\\$|fit_signature", shown)))

  broken <- list(zero_used = bs, negative_p = bs,
                 missing_accounting = bs, dropped_item = bs)
  broken$zero_used$B_used <- 0L
  broken$negative_p$items$chisq_p_boot_adj[1L] <- -4
  broken$missing_accounting$B_nonconverged <- NULL
  broken$dropped_item$items <- broken$dropped_item$items[-1L, , drop = FALSE]
  for (nm in names(broken))
    expect_error(.validate_fit_bootstrap(broken[[nm]], fit),
                 "incomplete or internally inconsistent", info = nm)
})

test_that("integrity signatures are portable and retain legacy validation", {
  x <- list(a = 1L, b = c(NA_real_, pi), label = "caf\u00e9")
  current <- .fit_boot_md5(x)
  expect_match(current, "^xdr3:[0-9a-f]{32}$")
  expect_true(.fit_boot_hash_matches(current, x))
  legacy <- .fit_boot_md5_legacy_candidates(x)
  expect_length(legacy, 2L)
  expect_true(all(vapply(legacy, .fit_boot_hash_matches, logical(1), x = x)))
  expect_false(.fit_boot_hash_matches(current, within(x, a <- 2L)))
  e <- new.env(parent = emptyenv())
  e$a <- 1L
  f <- ~ predictor
  environment(f) <- e
  formula_hash <- .fit_boot_md5(list(formula = f))
  e$unrelated_mutation <- 2L
  expect_identical(.fit_boot_md5(list(formula = f)), formula_hash)

  fit <- rasch(simb(200, seq(-1, 1, length.out = 5), seed = 301))
  bs <- fbq(fit, B = 9, seed = 302)
  old <- bs
  old$fit_signature$fingerprint <-
    .fit_boot_md5_legacy_candidates(fit)[2L]
  unsigned <- unclass(old)
  unsigned$result_signature <- NULL
  old$result_signature <- .fit_boot_md5_legacy_candidates(unsigned)[2L]
  expect_no_error(.validate_fit_bootstrap(old, fit))
})

test_that("the chi-square is read in one tail and the others in two", {
  # a statistic far below its null is misfit for the fit residual (a steeper
  # item than the model predicts) but not for a discrepancy that can only be
  # too large
  null <- c(rep(1, 50), rep(2, 50))
  expect_equal(.boot_p(0, null, "upper"), 1)
  expect_lt(.boot_p(0, null, "two"), 0.05)
  expect_lt(.boot_p(99, null, "upper"), 0.05)
  expect_lt(.boot_p(99, null, "two"), 0.05)
  # two-sided doubles the smaller tail and is capped at 1
  expect_lte(.boot_p(1.5, null, "two"), 1)
  expect_equal(.boot_p(NA_real_, null, "upper"), NA_real_)
})

test_that("the whole-test readings carry their own nulls", {
  fit <- rasch(simb(400, seq(-1.5, 1.5, length.out = 6), seed = 16))
  bs <- fbq(fit, B = 29, seed = 3)
  expect_equal(bs$total$chisq, fit$total_chisq)
  expect_equal(bs$total$fit_resid_mean, fit$item_fit_summary$mean)
  expect_equal(bs$total$fit_resid_sd, fit$item_fit_summary$sd)
  for (nm in c("chisq_p_boot", "fit_resid_mean_p_boot", "fit_resid_sd_p_boot"))
    expect_true(bs$total[[nm]] > 0 && bs$total[[nm]] <= 1, info = nm)
})

test_that("a seed reproduces the bootstrap and leaves the caller's stream alone", {
  fit <- rasch(simb(250, seq(-1.2, 1.2, length.out = 5), seed = 4))
  a <- fbq(fit, B = 19, seed = 7)
  b <- fbq(fit, B = 19, seed = 7)
  expect_equal(a$items$chisq_p_boot, b$items$chisq_p_boot)
  expect_equal(a$items$fit_resid_p_boot, b$items$fit_resid_p_boot)
  expect_equal(a$replicates$chisq, b$replicates$chisq)
  set.seed(99)
  before <- runif(1)
  set.seed(99)
  invisible(fbq(fit, B = 9, seed = 3))
  expect_equal(runif(1), before)
})

test_that("the replicates are drawn before dispatch, not inside a worker", {
  # what makes the result independent of the worker count is that every
  # random draw a replicate depends on is made in the coordinating process.
  # These tests run serially by policy (setup-efrm-workers.R), so the
  # property is checked where it lives: the same seed must fix the draws
  # regardless of how many replicates are asked for, so a longer run must
  # extend a shorter one rather than redraw it.
  fit <- rasch(simb(250, seq(-1.2, 1.2, length.out = 5), seed = 13))
  set.seed(4); th_a <- .fit_theta("resample", fit$person$theta, fit$psi)
  set.seed(4); th_b <- .fit_theta("resample", fit$person$theta, fit$psi)
  expect_equal(th_a, th_b)
  expect_equal(fbq(fit, B = 19, seed = 4)$replicates$chisq,
               fbq(fit, B = 19, seed = 4)$replicates$chisq)
})

test_that("the schemes for person locations all run and stay in range", {
  fit <- rasch(simb(250, seq(-1.2, 1.2, length.out = 5), seed = 5))
  for (s in c("conditional", "resample", "fixed", "normal")) {
    bs <- suppressWarnings(fit_bootstrap(fit, B = 19, theta = s, seed = 2))
    expect_identical(bs$theta, s)
    expect_true(all(bs$items$chisq_p_boot > 0 & bs$items$chisq_p_boot <= 1))
    if (s == "conditional") expect_s3_class(bs$persons, "data.frame")
    else expect_null(bs$persons)
  }
})

test_that("fixed locations remain paired with their response rows", {
  X <- rbind(c(NA, NA), c(0, 1), c(1, 0))
  th <- c(NA_real_, -1, 1)
  prepared <- .fit_theta_source(th, X, "fixed")
  expect_equal(prepared, c(0, -1, 1))
  expect_identical(.fit_theta("fixed", prepared, list(), n = nrow(X)),
                   prepared)
  X[1, 1] <- 0
  expect_error(.fit_theta_source(th, X, "fixed"),
               "finite person location", class = "rasch_refusal")
  expect_equal(.fit_theta_source(th, X, "resample"), c(-1, 1))
})

test_that("ability resampling treats a singleton as a value, not a sequence", {
  # base sample(3, n, replace = TRUE) samples 1:3; sampling positions is
  # robust if an internal or future caller supplies a singleton population.
  expect_equal(.fit_theta("resample", 3, list(), n = 8), rep(3, 8))
  expect_equal(.fit_theta("normal", 3,
                          list(var_theta = 0, mean_error_var = 0), n = 8),
               rep(3, 8))
})

test_that("source-tree bootstrap workers do not load another installation", {
  expect_warning(out <- with_mocked_bindings(
    .rasch_boot_apply(3L, identity, workers = 2L, label = "test bootstrap"),
    .rasch_namespace_is_installed = function() FALSE,
    .package = "rasch"), "source tree")
  expect_equal(out, as.list(1:3))
})

test_that("the maxT reference uses the joint null and never lowers raw p", {
  null <- cbind(a = 1:40, b = 40:1, c = c(1:20, 20:1))
  z <- .boot_maxt(c(40, 40, 20), null, "upper", min_success = 30L)
  expect_equal(z$family_n, 3L)
  expect_equal(z$family_boot, 40L)
  expect_true(all(z$p_adj >= z$p, na.rm = TRUE))
  expect_true(all(z$p_adj <= 1, na.rm = TRUE))
  expect_true(all(is.finite(z$p_adj)))

  # A constant member remains in the declared family. Degenerate null
  # columns are compared exactly rather than silently dropping that test.
  constant <- null; constant[, 3] <- 20
  zc <- .boot_maxt(c(40, 40, 20), constant, "upper", min_success = 30L)
  expect_true(all(is.finite(zc$p_adj)))
  expect_equal(zc$family_n, 3L)
  expect_equal(zc$adjusted_n, 3L)
  expect_true(all(is.finite(.boot_maxt(c(40, 40, 20), constant, "upper",
    min_success = 30L, mode = "centred")$p_adj)))
  expect_true(all(is.finite(.boot_maxt(c(40, 40, 20), constant, "upper",
    min_success = 30L, mode = "raw")$p_adj)))

  null[1:20, 3] <- NA
  z2 <- .boot_maxt(c(40, 40, 20), null, "upper", min_success = 30L)
  expect_true(is.na(z2$p[3]))
  expect_equal(z2$family_n, 3L)
  expect_true(all(is.na(z2$p_adj)))
})

test_that("studentised maxT externally standardises every null row", {
  set.seed(914)
  B <- 40L; K <- 5L
  null <- matrix(rnorm(B * K), B, K)
  observed <- rnorm(K)
  got <- .boot_maxt(observed, null, "two", min_success = 36L)

  Z <- matrix(NA_real_, B, K)
  for (i in seq_len(B)) {
    training <- null[-i, , drop = FALSE]
    Z[i, ] <- (null[i, ] - colMeans(training)) /
      apply(training, 2L, sd)
  }
  z_obs <- abs((observed - colMeans(null)) / apply(null, 2L, sd))
  max_null <- apply(abs(Z), 1L, max)
  expected <- vapply(z_obs, function(z)
    (1 + sum(max_null >= z)) / (B + 1), numeric(1))
  expect_equal(got$p_adj, pmax(got$p, expected), tolerance = 1e-14)
})

test_that("transformed maxT references exclude non-finite joint rows", {
  null <- cbind(a = seq(0.5, 2.5, length.out = 40),
                b = seq(2.5, 0.5, length.out = 40))
  null[1L, 1L] <- 0                         # log(0) is -Inf, not a usable row
  z <- .boot_maxt(c(2, 2), null, "two", min_success = 30L,
                  mode = "studentised", transform = log)
  expect_equal(z$family_boot, 39L)
  expect_true(all(is.finite(z$p_adj)))
  expect_true(all(z$p_adj >= z$p))
})

test_that("tiny maxT runs return marginal probabilities without invalid studentisation", {
  z1 <- .boot_maxt(1, matrix(0, 1, 1), "upper", min_success = 1L)
  z2 <- .boot_maxt(c(1, 2), matrix(c(0, 1, 1, 0), 2, 2), "upper",
                   min_success = 2L)
  expect_true(all(is.finite(z1$p)))
  expect_true(all(is.na(z1$p_adj)))
  expect_true(all(is.finite(z2$p)))
  expect_true(all(is.na(z2$p_adj)))
})

test_that("report transport attributes do not change fitted-model identity", {
  fit <- rasch(simb(100, seq(-1, 1, length.out = 5), seed = 915))
  signature <- .fit_boot_signature(fit)
  attr(fit, "report_dimensionality") <- list(marker = 1L)
  attr(fit, "report_invariance") <- list(marker = 2L)
  expect_true(.fit_boot_signature_matches(signature, fit))
})

test_that("the bootstrap statistic carries the bias the asymptotic df ignores", {
  # under a true model the replicated statistic sits well above its nominal
  # degrees of freedom: that gap is the miscalibration being corrected, and
  # the corrected p-values are accordingly larger than the asymptotic ones
  fit <- rasch(simb(1500, seq(-1.5, 1.5, length.out = 6), seed = 6))
  bs <- fbq(fit, B = 49, seed = 11)
  expect_true(mean(colMeans(bs$replicates$chisq)) >
                1.3 * stats::median(fit$items$df))
  expect_true(mean(bs$items$chisq_p_boot > fit$items$p) >= 5 / 6)
})

test_that("missing data are carried into every replicate", {
  X <- simb(300, seq(-1.5, 1.5, length.out = 6), seed = 8)
  X[cbind(sample(300, 60), sample(6, 60, replace = TRUE))] <- NA
  fit <- rasch(X)
  bs <- fbq(fit, B = 19, seed = 1)
  expect_true(all(bs$items$n_boot > 0))
  expect_true(all(is.finite(bs$items$chisq_p_boot)))
})

test_that("familywise adjustment covers each statistic's own family", {
  fit <- rasch(simb(300, seq(-1.5, 1.5, length.out = 6), seed = 9))
  bs <- fbq(fit, B = 19, seed = 5)
  for (st in STATS) {
    p <- bs$items[[paste0(st, "_p_boot")]]
    adj <- bs$items[[paste0(st, "_p_boot_adj")]]
    expect_equal(adj, stats::p.adjust(p, method = "holm"), info = st)
    expect_true(all(adj <= 1), info = st)
  }
})

test_that("the bootstrap refuses designs it cannot generate from", {
  fit <- rasch(simb(200, seq(-1, 1, length.out = 5), seed = 10))
  expect_error(fbq(list(items = 1)), "fitted model from rasch")
  for (cl in c("rasch_efrm", "rasch_mfrm", "rasch_explanatory")) {
    f2 <- fit; class(f2) <- c(cl, "rasch")
    expect_error(fbq(f2, B = 5), "generating structure",
                 class = "rasch_refusal")
  }
  f3 <- fit; f3$disc <- c(1, 1, 1, 1, 2)
  expect_error(fbq(f3, B = 5), "frame units", class = "rasch_refusal")
})

test_that("the bootstrap validates its own controls", {
  fit <- rasch(simb(200, seq(-1, 1, length.out = 5), seed = 12))
  expect_error(fbq(fit, B = 0), "whole number")
  expect_error(fbq(fit, B = c(10, 20)), "whole number")
  expect_error(fbq(fit, B = 2.5), "whole number")
  expect_error(fbq(fit, B = 5, workers = 0), "whole number")
  expect_error(fbq(fit, B = 5, seed = c(1, 2)), "whole number")
  expect_error(fbq(fit, B = 5, seed = 1.5), "whole number")
  expect_error(fbq(fit, B = 5, seed = -1), "whole number")
  expect_error(fbq(fit, B = 5, seed = .Machine$integer.max + 1),
               "integer range")
  expect_error(fbq(fit, B = 5, theta = "wishful"), "should be one of")

  failed <- fit; failed$est$converged <- FALSE
  expect_error(fbq(failed, B = 5), "observed Rasch fit did not converge",
               class = "rasch_refusal")
})

test_that("a non-converged replicate is not admitted to the null", {
  d <- simulate_rasch(250, 8, model = "PCM", n_categories = 4, seed = 1)
  f <- suppressWarnings(rasch(d, id = "id", model = "PCM"))
  z <- suppressWarnings(.fit_refit(
    f$X, f$model, f$n_groups, anchors = NULL, expected_m = f$m,
    maxit = 1L, tol = 1e-12))
  expect_identical(.fit_boot_status(z), "nonconverged")
})

test_that("a replicate with a shortened item scale is not admitted", {
  d <- simulate_rasch(30, 4, model = "PCM", n_categories = 3, seed = 2)
  fit <- rasch(d, id = "id", model = "PCM")
  set.seed(12001)
  xb <- .fit_gen_conditional(fit$X, fit$tau_list, NULL)
  expect_false(identical(
    apply(xb, 2L, max, na.rm = TRUE), as.integer(fit$m)
  ))
  z <- suppressWarnings(.fit_refit(
    xb, fit$model, fit$n_groups, anchors = NULL, expected_m = fit$m,
    maxit = 60L, tol = 1e-8
  ))
  expect_identical(.fit_boot_status(z), "error")
})


test_that("an anchored fit is bootstrapped under its own anchors", {
  X <- simb(300, seq(-1.5, 1.5, length.out = 6), seed = 15)
  free <- rasch(X)
  anc <- data.frame(item = c("I1", "I2"), k = 1L,
                    tau = free$thresholds$tau[free$thresholds$item %in% 1:2])
  fit <- rasch(X, anchors = anc)
  bs <- fbq(fit, B = 19, seed = 8)
  expect_true(all(bs$items$chisq_p_boot > 0 & bs$items$chisq_p_boot <= 1))
  expect_true(all(bs$items$n_boot > 0))
})

test_that("polytomous fits retain an adequate same-scale bootstrap null", {
  for (mod in c("PCM", "RSM")) {
    d <- simulate_rasch(400, 6, model = mod, n_categories = 4, seed = 21)
    # A short exploratory run may lose a replicate whose category scale is no
    # longer the fitted scale. It must account for those draws and retain the
    # documented minimum rather than admitting the shorter model.
    bs <- suppressWarnings(
      fbq(rasch(d, id = "id", model = mod), B = 29, seed = 2))
    expect_true(all(bs$items$chisq_p_boot > 0 & bs$items$chisq_p_boot <= 1))
    expect_true(all(is.finite(bs$items$fit_resid_p_boot)))
    expect_equal(bs$B_used + bs$B_failed, bs$B)
    expect_gte(bs$B_used, bs$minimum_usable)
  }
})

test_that("a thinned null is reported rather than left to be inferred", {
  d <- simulate_rasch(400, 6, model = "RSM", n_categories = 4, seed = 21)
  fit <- rasch(d, id = "id", model = "RSM")
  # fbq() would also swallow the warning under test; suppress only the
  # resolution-floor warning by matching this one specifically
  w <- character(0)
  withCallingHandlers(
    fit_bootstrap(fit, B = 29, theta = "normal", seed = 2),
    warning = function(cnd) {
      w <<- c(w, conditionMessage(cnd)); invokeRestart("muffleWarning")
    })
  expect_true(any(grepl("were unusable", w)))
})

test_that("normal bootstrap generation does not revert to noisy estimates", {
  theta <- c(-2, -1, 0, 1, 2)
  zero <- rasch:::.fit_theta(
    "normal", theta,
    list(var_theta = stats::var(theta),
         mean_error_var = stats::var(theta) + 1),
    n = 40L)
  expect_identical(zero, rep(mean(theta), 40L))
  expect_error(
    rasch:::.fit_theta(
      "normal", theta,
      list(var_theta = NA_real_, mean_error_var = NA_real_), n = 5L),
    "error-corrected person variance is unavailable")
})

test_that("a depleted or statistic-specific null is withheld", {
  fit <- rasch(simb(200, seq(-1, 1, length.out = 5), seed = 18))
  good <- lapply(STATS, function(st) fit$items[[st]])
  names(good) <- STATS
  good$persons <- lapply(c("fit_resid", "infit_ms", "outfit_ms",
                           "infit_z", "outfit_z"), function(st)
    fit$person[[st]])
  names(good$persons) <- c("fit_resid", "infit_ms", "outfit_ms",
                           "infit_z", "outfit_z")

  expect_error(
    suppressWarnings(with_mocked_bindings(
      fit_bootstrap(fit, B = 40, seed = 1),
      .rasch_boot_apply = function(n, fun, ...) c(
        rep(list(good), 29L), rep(list(NULL), n - 29L)),
      .package = "rasch")),
    "at least 36", class = "rasch_refusal")

  refused <- tryCatch(suppressWarnings(with_mocked_bindings(
    fit_bootstrap(fit, B = 40, seed = 1),
    .rasch_boot_apply = function(n, fun, ...) c(
      rep(list(good), 29L),
      rep(list(.fit_boot_failure("nonconverged")), 6L),
      rep(list(.fit_boot_failure("error")), 5L)),
    .package = "rasch")), error = identity)
  expect_s3_class(refused, "rasch_fit_bootstrap_refusal")
  expect_equal(refused$B_used, 29L)
  expect_equal(refused$B_nonconverged, 6L)
  expect_equal(refused$B_errors, 5L)

  expect_error(suppressWarnings(with_mocked_bindings(
    fit_bootstrap(fit, B = 40, seed = 1),
    .rasch_boot_apply = function(n, fun, ...) c(
      rep(list(good), 35L), rep(list(NULL), n - 35L)),
    .package = "rasch")), "at least 36", class = "rasch_refusal")

  bs <- suppressWarnings(with_mocked_bindings(
    fit_bootstrap(fit, B = 40, seed = 1),
    .rasch_boot_apply = function(n, fun, ...) c(
      rep(list(good), 36L), rep(list(NULL), n - 36L)),
    .package = "rasch"))
  expect_equal(bs$B_used, 36L)
  expect_equal(bs$B_failed, 4L)
  expect_equal(bs$B_nonconverged, 0L)
  expect_equal(bs$B_errors, 4L)
  expect_equal(bs$minimum_usable, 36L)
  expect_true(all(bs$items$n_boot_chisq == 36L))

  bad_stat <- good
  bad_stat$outfit_z[1L] <- NA_real_
  w <- character(0)
  bs2 <- withCallingHandlers(with_mocked_bindings(
    fit_bootstrap(fit, B = 40, seed = 1),
    .rasch_boot_apply = function(n, fun, ...) rep(list(bad_stat), n),
    .package = "rasch"), warning = function(cnd) {
      w <<- c(w, conditionMessage(cnd)); invokeRestart("muffleWarning")
    })
  expect_equal(bs2$items$n_boot_outfit_z[1L], 0L)
  expect_true(is.na(bs2$items$outfit_z_p_boot[1L]))
  ok <- is.finite(bs2$items$outfit_z_p_boot)
  expect_equal(bs2$items$outfit_z_p_boot_adj[ok],
               p.adjust(bs2$items$outfit_z_p_boot[ok], method = "holm",
                        n = nrow(bs2$items)))
  expect_true(any(grepl("too few usable replicated statistics for outfit_z", w)))

  bad_fr <- good
  bad_fr$fit_resid[1L] <- NA_real_
  bs_fr <- suppressWarnings(with_mocked_bindings(
    fit_bootstrap(fit, B = 40, seed = 1),
    .rasch_boot_apply = function(n, fun, ...) rep(list(bad_fr), n),
    .package = "rasch"))
  expect_equal(bs_fr$total$n_boot_fit_resid_mean, 0L)
  expect_equal(bs_fr$total$n_boot_fit_resid_sd, 0L)
  expect_true(is.na(bs_fr$total$fit_resid_mean_p_boot))
  expect_true(is.na(bs_fr$total$fit_resid_sd_p_boot))

  bs3 <- suppressWarnings(with_mocked_bindings(
    fit_bootstrap(fit, B = 40, seed = 1),
    .rasch_boot_apply = function(n, fun, ...) c(
      rep(list(good), 36L),
      rep(list(.fit_boot_failure("nonconverged")), 2L),
      rep(list(.fit_boot_failure("error")), 2L)),
    .package = "rasch"))
  expect_equal(bs3$B_used, 36L)
  expect_equal(bs3$B_nonconverged, 2L)
  expect_equal(bs3$B_errors, 2L)
})

test_that("the conditional generator preserves scores, masks, and theory", {
  d <- simulate_rasch(250, 8, model = "PCM", n_categories = 4,
                      missing = 0.08, seed = 11)
  f <- rasch(d, id = "id")
  nm <- is.na(f$X)
  set.seed(2)
  Xb <- .fit_gen_conditional(f$X, f$tau_list, nm)
  expect_identical(is.na(Xb), nm)
  expect_equal(rowSums(Xb, na.rm = TRUE), rowSums(f$X, na.rm = TRUE))
  expect_gt(mean(Xb != f$X, na.rm = TRUE), 0.1)
  # enumerable case: three dichotomous items, every score 1; the conditional
  # category probabilities are the scaled item weights
  tau3 <- list(-1, 0, 1)
  Xs <- matrix(0L, 4000, 3); Xs[, 1] <- 1L
  set.seed(3)
  Xg <- .fit_gen_conditional(Xs, tau3, NULL)
  w <- exp(-c(-1, 0, 1))
  expect_lt(max(abs(colMeans(Xg) - w / sum(w))), 0.03)

  # Finite extreme thresholds remain finite on the log scale and still
  # preserve the conditioning score exactly.
  Xe <- matrix(c(1L, 0L, 0L,
                 0L, 1L, 0L,
                 0L, 0L, 1L), 3, 3, byrow = TRUE)
  set.seed(4)
  Xext <- .fit_gen_conditional(Xe, list(-1000, 0, 1000), NULL)
  expect_true(all(is.finite(Xext)))
  expect_equal(rowSums(Xext), rowSums(Xe))
  Xm <- rbind(Xe, c(NA_integer_, NA_integer_, NA_integer_))
  Xmiss <- .fit_gen_conditional(Xm, list(-1000, 0, 1000), is.na(Xm))
  expect_true(all(is.na(Xmiss[4L, ])))
})

test_that("the conditional bootstrap preserves a booklet's missingness link", {
  set.seed(5)
  N <- 240
  grp <- rep(c("low", "high"), each = N / 2)
  th <- rnorm(N, ifelse(grp == "low", -0.75, 0.75), 1)
  X <- sapply(seq(-2, 2, length.out = 12),
              function(dd) rbinom(N, 1, plogis(th - dd)))
  colnames(X) <- sprintf("I%02d", 1:12)
  X[grp == "low", 9:12] <- NA
  X[grp == "high", 1:4] <- NA
  f <- rasch(X)
  obs <- mean(f$person$theta[grp == "high"], na.rm = TRUE) -
    mean(f$person$theta[grp == "low"], na.rm = TRUE)
  set.seed(9)
  Xb <- .fit_gen_conditional(f$X, f$tau_list, is.na(f$X))
  fb <- rasch(Xb)
  rep_diff <- mean(fb$person$theta[grp == "high"], na.rm = TRUE) -
    mean(fb$person$theta[grp == "low"], na.rm = TRUE)
  # abilities resampled independently of the mask would put this near zero
  expect_gt(obs, 0.8)
  expect_lt(abs(rep_diff - obs), 0.4)
})

test_that("bootstrap acceptance needs 90 percent for inferential runs", {
  expect_identical(.fit_min_boot_success(5), 3L)
  expect_identical(.fit_min_boot_success(29), 15L)
  expect_identical(.fit_min_boot_success(30), 30L)
  expect_identical(.fit_min_boot_success(40), 36L)
  expect_identical(.fit_min_boot_success(100), 90L)
  expect_identical(.rasch_min_boot_success(30), 30L)
  expect_identical(.rasch_min_boot_success(40), 30L)
  expect_identical(.rasch_min_boot_success(58), 30L)
  expect_identical(.rasch_min_boot_success(100), 51L)
  expect_identical(.rasch_min_boot_success(300), 151L)
  # the covariance guard also clears the number of independently estimable
  # directions supplied by its caller
  expect_identical(.rasch_min_boot_success(100, n_quantities = 80), 81L)
  expect_identical(.rasch_min_boot_success(40, n_quantities = 10), 30L)
})

test_that("a B below the adjusted-inference floor warns, naming the remedy", {
  fit <- rasch(simb(200, seq(-1, 1, length.out = 5), seed = 17))
  expect_warning(fit_bootstrap(fit, B = 99, seed = 1),
                 "cannot reach Holm-adjusted significance")
  expect_no_warning(fbq(fit, B = 199, seed = 1))
  expect_silent(capture.output(suppressMessages(
    invisible(fit_bootstrap(fit, B = 200, seed = 1)))))
})

test_that("an automatic interval count is held across replicates", {
  # The automatic rule reads the non-extreme count, which the ability-
  # sampling generators move across a 50-person boundary from replicate to
  # replicate. A null mixing interval counts, and so degrees of freedom, the
  # observed chi-square never used is not that statistic's null.
  d <- simulate_rasch(160, 6, seed = 8)
  f_auto <- rasch(d, id = "id")
  f_fix <- rasch(d, id = "id", n_groups = f_auto$n_groups)
  expect_equal(f_auto$n_groups, 3)
  expect_equal(f_auto$item_trait$chisq, f_fix$item_trait$chisq)

  fg <- .fit_boot_realised_groups(f_auto)
  expect_equal(fg$common, f_auto$n_groups)
  expect_null(fg$item)
  # a replicate that loses persons to the extremes would otherwise be scored
  # on two intervals instead of the observed three
  th <- f_auto$person$theta[1:120]
  ex <- rep(FALSE, 120)
  expect_equal(attr(.class_intervals(th, ex, NULL), "n_groups"), 2)
  expect_equal(attr(.class_intervals(th, ex, fg$common), "n_groups"), 3)

  # with the count held, the automatic fit and the same count requested
  # explicitly share one null
  b_auto <- fbq(f_auto, B = 99, theta = "resample", seed = 1, workers = 1)
  b_fix <- fbq(f_fix, B = 99, theta = "resample", seed = 1, workers = 1)
  expect_equal(b_auto$items$chisq_p_boot, b_fix$items$chisq_p_boot)
})

test_that("per-item interval counts are resolved on the observed responders", {
  X <- simb(400, seq(-1, 1, length.out = 8), seed = 744)
  X[121:400, 7:8] <- NA_integer_
  fit <- rasch(X)
  fg <- .fit_boot_realised_groups(fit)
  # the automatic rule applies per item with missing responses, so the held
  # counts are per item too
  expect_equal(fg$item,
               vapply(fit$ci_item,
                      function(g) as.integer(max(g, na.rm = TRUE)), 1L))
  expect_false(fg$item[7] == fg$item[1])
  # holding those counts reproduces the observed statistics exactly
  held <- .fit_refit(fit$X, fit$model, .refit_n_groups(fit),
                     fit$refit_spec$anchors, fit$m, 60, 1e-8,
                     fixed_groups = fg)
  expect_equal(unname(held$chisq), fit$items$chisq, tolerance = 1e-10)
})

test_that("fit_bootstrap hands the replicates the count it resolved", {
  # The interval count a replicate is scored on now travels in fixed_groups,
  # not in n_groups, so pin that argument rather than the request alone: an
  # automatic request carries the counts realised by the observed fit, an
  # explicit one carries none because the request is already fixed. Every
  # generator is pinned, including the default: its non-extreme set never
  # moves, so holding the counts leaves its null untouched, but the value
  # handed down is what decides the intervals and is pinned as such.
  original <- .fit_refit
  seen <- list()
  local_mocked_bindings(
    .fit_refit = function(X, model, n_groups, ..., fixed_groups = NULL) {
      seen[[length(seen) + 1L]] <<- list(n_groups = n_groups,
                                         fixed_groups = fixed_groups)
      original(X, model, n_groups, ..., fixed_groups = fixed_groups)
    })
  X <- simb(400, seq(-1, 1, length.out = 8), seed = 744)
  X[121:400, 7:8] <- NA_integer_
  for (gen in c("conditional", "resample")) for (ng in list(NULL, 4L)) {
    fit <- rasch(X, n_groups = ng)
    seen <- list()
    bs <- fbq(fit, B = 20, theta = gen, seed = 745, workers = 1)
    expect_length(seen, 20L)
    expect_true(all(vapply(seen, function(s) identical(s$n_groups, ng),
                           logical(1))))
    expected <- if (is.null(ng)) .fit_boot_realised_groups(fit) else NULL
    expect_true(all(vapply(seen, function(s)
      identical(s$fixed_groups, expected), logical(1))))
    if (is.null(ng)) {
      expect_equal(expected$common, fit$n_groups)
      expect_equal(expected$item,
                   vapply(fit$ci_item,
                          function(g) as.integer(max(g, na.rm = TRUE)), 1L))
    }
    expect_equal(bs$B_used, 20)
  }
})

test_that("an ability-sampling null saved under the old rule is refused", {
  # The held interval count changed what the resample/fixed/normal nulls
  # mean, so those stored results must be recomputed rather than reused. The
  # conditional generator retains every raw score, so its non-extreme set --
  # and its automatic interval count -- was already constant across
  # replicates; those results are unchanged and stay current.
  fit <- rasch(simb(154, seq(-1.2, 1.2, length.out = 6), seed = 11))
  cond <- fbq(fit, B = 20, seed = 5, workers = 1)
  samp <- fbq(fit, B = 20, theta = "resample", seed = 5, workers = 1)
  expect_identical(cond$algorithm, "loo-maxt-2")
  expect_identical(samp$algorithm, "loo-maxt-3")
  expect_no_error(.validate_fit_bootstrap(cond, fit))
  expect_no_error(.validate_fit_bootstrap(samp, fit))

  old <- unclass(samp)
  old$algorithm <- "loo-maxt-2"
  old$result_signature <- NULL
  old$result_signature <- .fit_boot_md5(old)
  class(old) <- class(samp)
  expect_error(.validate_fit_bootstrap(old, fit), "recompute")
})

test_that("paired-comparison fit bootstrap follows the fitted design", {
  d <- simulate_btl(5, 20, reps_per_pair = 8, seed = 41)
  fit <- btl(d, "object_a", "object_b", winner = "winner", judge = "judge")
  bs <- suppressWarnings(fit_bootstrap(fit, B = 19, workers = 1, seed = 4))
  expect_s3_class(bs, "rasch_fit_bootstrap")
  expect_identical(bs$model_kind, "btl")
  expect_identical(bs$model, "paired comparisons")
  expect_equal(bs$B_used, 19L)
  expect_identical(bs$objects$object, fit$objects$object)
  expect_identical(bs$judges$judge, fit$judges$judge)
  expect_equal(bs$pairs[, c("object_a", "object_b")],
               fit$pairs[, c("object_a", "object_b")])
  expect_true(bs$total$chisq_p_boot > 0 && bs$total$chisq_p_boot <= 1)
  for (tab in c("pairs", "objects", "judges")) {
    pcols <- grep("_p_boot$", names(bs[[tab]]), value = TRUE)
    for (nm in pcols) {
      adj <- bs[[tab]][[paste0(nm, "_adj")]]
      expect_true(all(adj >= bs[[tab]][[nm]], na.rm = TRUE),
                  info = paste(tab, nm))
    }
  }
  broken <- bs
  broken$objects$fit_resid_p_boot_adj[1L] <- -2
  expect_error(.validate_fit_bootstrap(broken, fit),
               "incomplete or internally inconsistent")
})

test_that("non-converged observed paired-comparison fits are refused", {
  d <- simulate_btl(5, 20, reps_per_pair = 4, seed = 47)
  fit <- btl(d, "object_a", "object_b", winner = "winner", judge = "judge")
  fit$converged <- FALSE
  expect_error(fit_bootstrap(fit, B = 5, seed = 1),
               "observed paired-comparison fit did not converge",
               class = "rasch_refusal")
})

test_that("paired-comparison fit bootstrap does not invent judges", {
  d <- simulate_btl(5, 20, reps_per_pair = 6, seed = 45)
  fit <- btl(d, "object_a", "object_b", winner = "winner")
  bs <- suppressWarnings(fit_bootstrap(fit, B = 5, workers = 1, seed = 7))
  expect_null(bs$judges)
  expect_null(bs$adjustment$judges)
  expect_null(bs$replicates$judges)
})

test_that("paired-comparison fit bootstrap retains explanatory restrictions", {
  set.seed(46)
  predictors <- data.frame(object = LETTERS[1:6], x = seq(-1, 1, length.out = 6))
  pr <- t(combn(predictors$object, 2))
  d <- data.frame(object_a = rep(pr[, 1], each = 20),
                  object_b = rep(pr[, 2], each = 20))
  bx <- setNames(predictors$x, predictors$object)
  d$winner <- ifelse(runif(nrow(d)) <
                       plogis(bx[d$object_a] - bx[d$object_b]),
                     d$object_a, d$object_b)
  fit <- btl_explanatory(d, predictors, ~ x,
                         object_a = "object_a", object_b = "object_b",
                         winner = "winner")
  bs <- suppressWarnings(fit_bootstrap(fit, B = 5, workers = 1, seed = 8))
  expect_identical(bs$objects$object, fit$objects$object)
  expect_equal(bs$B_used, 5L)
})

test_that("paired-comparison bootstrap covers ordered categories and history", {
  d <- simulate_btl(5, 24, reps_per_pair = 6, model = "polytomous",
                    n_categories = 4, seed = 42)
  fit <- btl(d, "object_a", "object_b", response = "response",
             judge = "judge", thresholds = "pc")
  generated <- .btl_boot_data(fit)
  expect_true(is.ordered(generated$response))
  expect_identical(levels(generated$response), as.character(0:fit$m))
  expect_no_error(suppressWarnings(
    fit_bootstrap(fit, B = 9, workers = 1, seed = 5)))
  stamped <- .new_fit_bootstrap(
    list(pairs = fit$pairs, objects = fit$objects, total = list(),
         B = 1L, B_used = 1L), fit, "btl")
  changed_threshold <- fit
  changed_threshold$thresholds$tau[1L] <-
    changed_threshold$thresholds$tau[1L] + 0.5
  expect_error(.validate_fit_bootstrap(stamped, changed_threshold),
               "different fitted model")

  h <- simulate_btl(5, 30, reps_per_pair = 6,
                    dependence = list(exposure = .2, carry_over = .1),
                    seed = 43)
  fh <- btl(h, "object_a", "object_b", winner = "winner",
            judge = "judge", order = "order")
  bh <- suppressWarnings(fit_bootstrap(fh, B = 9, workers = 1, seed = 6))
  expect_gte(bh$B_used, bh$minimum_usable)
  expect_true(all(vapply(bh$replicates$objects, nrow, 0L) == bh$B_used))
  changed_history <- fh
  changed_history$dependence$estimate[1L] <-
    changed_history$dependence$estimate[1L] + 0.5
  expect_error(.validate_fit_bootstrap(bh, changed_history),
               "different fitted model")

  dropped <- fh
  dropped$dependence <- dropped$dependence[-1L, , drop = FALSE]
  expect_error(.btl_boot_data(dropped), "dropped.*effect",
               class = "rasch_refusal")
})

test_that("paired-comparison bootstrap refuses designs it cannot regenerate", {
  d <- simulate_btl(5, 20, reps_per_pair = 6, seed = 44)
  ordinary <- btl(d, "object_a", "object_b", winner = "winner",
                  judge = "judge")
  oversized <- ordinary
  oversized$comparisons$weight[1L] <- .Machine$integer.max + 1
  expect_error(.btl_boot_data(oversized), "exceeds the integer range",
               class = "rasch_refusal")
  d$winner[1:4] <- "tie"
  half <- btl(d, "object_a", "object_b", winner = "winner", judge = "judge",
              ties = "half")
  expect_error(fit_bootstrap(half, B = 5, seed = 1), "half-weighted ties",
               class = "rasch_refusal")
  expect_error(fit_bootstrap(btl(d[-(1:4), ], "object_a", "object_b",
                                 winner = "winner", judge = "judge"),
                             B = 5, theta = "fixed", seed = 1),
               "theta.*not paired comparisons")
})
