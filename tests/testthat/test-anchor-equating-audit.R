test_that("a weak location anchor retains its exact item location", {
  set.seed(303)
  N <- 300L
  sim_one <- function(theta, tau) {
    score <- 0:length(tau)
    p <- exp(score * theta - c(0, cumsum(tau)))
    sample(score, 1L, prob = p)
  }
  theta <- rnorm(N)
  X <- sapply(seq_len(5L), function(i)
    vapply(theta, sim_one, integer(1), tau = c(-1, 0, 1)))
  colnames(X) <- paste0("I", seq_len(ncol(X)))
  category_one <- which(X[, 1L] == 1L)
  X[setdiff(category_one, category_one[seq_len(2L)]), 1L] <- 0L

  fit <- rasch(X, anchors = data.frame(item = "I1", k = NA_real_, tau = 0.5))
  expect_true(any(fit$thresholds$weak[fit$thresholds$item == 1L]))
  expect_true(all(is.na(fit$thresholds$se[
    fit$thresholds$item == 1L & fit$thresholds$weak
  ])))
  expect_equal(fit$items$location[1L], 0.5, tolerance = 1e-12)
  expect_equal(fit$items$se[1L], 0, tolerance = 1e-10)

  bank <- fit$items[1L, c("item", "location", "se", "max")]
  linked <- equate_tests(fit, bank, shift = "none")
  expect_identical(linked$n_common, 1L)
})

test_that("ordinary equating refuses incompatible exact anchors", {
  set.seed(502)
  X <- matrix(rbinom(1200L, 1L, 0.5), 200L, 6L,
              dimnames = list(NULL, paste0("I", 1:6)))
  fit <- rasch(X, anchors = data.frame(
    item = c("I1", "I2"), k = 1L, tau = c(0, 1)
  ))
  bank <- fit$items[, c("item", "location", "se", "max")]
  bank$se[] <- 0
  bank$location[bank$item == "I2"] <- 2

  expect_error(equate_tests(fit, bank, independent = TRUE),
               "zero pooled uncertainty.*incompatible origin shifts")
  expect_error(equate_tests(fit, bank[bank$item == "I2", ], shift = "none"),
               "zero pooled uncertainty.*declared fixed origin")
})

test_that("a numerical zero SE is not reclassified as an exact anchor", {
  set.seed(504)
  X <- matrix(rbinom(1200L, 1L, 0.5), 200L, 6L,
              dimnames = list(NULL, paste0("I", 1:6)))
  fit <- rasch(X)
  fit$items$se[] <- 0
  fit$est$cov_tau[,] <- 0
  bank <- fit$items[, c("item", "location", "se", "max")]
  bank$se[] <- 0
  bank$location[2L] <- bank$location[2L] + 0.5
  expect_no_error(eq <- equate_tests(fit, bank, independent = TRUE))
  expect_false(eq$inferential)
  expect_true(all(is.na(eq$table$p)))

  # A single spurious zero must not receive all of the link weight when the
  # remaining common items have usable estimated uncertainty.
  mixed <- rasch(X)
  mixed$items$se[1L] <- 0
  mixed$est$cov_tau[mixed$est$thr$item == 1L, ] <- 0
  mixed$est$cov_tau[, mixed$est$thr$item == 1L] <- 0
  mixed_bank <- mixed$items[, c("item", "location", "se", "max")]
  mixed_bank$location <- mixed$items$location - 0.25
  mixed_bank$location[1L] <- mixed$items$location[1L] - 5
  mixed_bank$se <- mixed$items$se
  mixed_link <- equate_tests(mixed, mixed_bank, independent = TRUE)
  expect_equal(mixed_link$shift, 0.25, tolerance = 1e-10)
  expect_identical(mixed_link$n, nrow(mixed_bank) - 1L)
  expect_match(mixed_link$note, "zero pooled estimated variance")
})

test_that("paired-comparison equating refuses incompatible exact anchors", {
  objects <- paste0("O", 1:4)
  fit <- structure(list(
    objects = data.frame(object = objects,
                         location = c(-1, -0.2, 0.3, 0.9),
                         se = c(0, 0, 0.2, 0.2)),
    cov_beta = diag(c(0, 0, 0.04, 0.04)),
    anchors = c(O1 = -1, O2 = -0.2),
    converged = TRUE, m = 1L, categories = 0:1,
    thr_structure = "none", clustered = FALSE
  ), class = "rasch_btl")
  bank <- data.frame(object = objects,
                     location = c(-1, 0.2, 0.3, 0.9),
                     se = c(0, 0, 0.2, 0.2))
  attr(bank, "cov_location") <- diag(bank$se^2)

  expect_error(btl_equate(fit, bank, independent = TRUE),
               "zero pooled uncertainty.*incompatible origin shifts")
  one <- bank[2L, ]
  attr(one, "cov_location") <- NULL
  expect_error(btl_equate(fit, one, shift = "none"),
               "zero pooled uncertainty.*declared fixed origin")
})

test_that("a zero BTL sandwich diagonal is not called a fixed anchor", {
  d <- simulate_btl(5, 20, reps_per_pair = 5, seed = 505)
  fit <- btl(d, "object_a", "object_b", "winner", judge = "judge")
  fit$objects$se[] <- 0
  fit$cov_beta[,] <- 0
  bank <- fit$objects[, c("object", "location", "se")]
  bank$se[] <- 0
  bank$location[2L] <- bank$location[2L] + 0.5
  attr(bank, "cov_location") <- matrix(0, nrow(bank), nrow(bank))
  expect_no_error(eq <- btl_equate(fit, bank, independent = TRUE))
  expect_false(eq$inferential)
  expect_true(all(is.na(eq$table$p)))

  mixed <- btl(d, "object_a", "object_b", "winner", judge = "judge")
  mixed$objects$se[1L] <- 0
  mixed$cov_beta[1L, ] <- 0
  mixed$cov_beta[, 1L] <- 0
  mixed_bank <- mixed$objects[, c("object", "location", "se")]
  mixed_bank$location <- mixed$objects$location - 0.25
  mixed_bank$location[1L] <- mixed$objects$location[1L] - 5
  mixed_bank$se <- mixed$objects$se
  mixed_link <- btl_equate(mixed, mixed_bank, independent = TRUE)
  expect_equal(mixed_link$shift, 0.25, tolerance = 1e-10)
  expect_match(paste(mixed_link$notes, collapse = " "),
               "zero pooled estimated variance")
})

test_that("an exact BTL anchor survives withheld judge covariance", {
  d <- simulate_btl(6, 8, reps_per_pair = 5, seed = 501)
  fit <- btl(d, "object_a", "object_b", "winner", judge = "judge",
             anchors = setNames(0.3, "O1"))
  expect_false(fit$cl$inference_available)
  expect_equal(fit$objects$se[fit$objects$object == "O1"], 0)
  expect_true(all(is.na(fit$objects$se[fit$objects$object != "O1"])))

  bank <- fit$objects[, c("object", "location", "se")]
  bank$se[] <- 0
  linked <- btl_equate(fit, bank, shift = "none")
  expect_true(all(is.na(linked$table$p_adj)))
  expect_identical(sum(is.finite(linked$table$se_1)), 1L)
})

test_that("unsupported person clustering withholds estimated uncertainty", {
  set.seed(503)
  X <- matrix(sample(0:2, 1000L, replace = TRUE), 200L, 5L,
              dimnames = list(NULL, paste0("I", 1:5)))

  one_cluster <- rasch(X, id = rep("P1", nrow(X)))
  expect_true(one_cluster$est$converged)
  expect_false(one_cluster$est$cluster_inference)
  expect_identical(one_cluster$est$cluster_support$n, 1L)
  expect_true(all(is.finite(one_cluster$items$location)))
  expect_true(all(is.na(one_cluster$items$se)))
  expect_true(all(is.na(one_cluster$est$cov_tau)))
  expect_match(paste(one_cluster$notes, collapse = " "),
               "uncertainty withheld.*fewer than 10 person clusters")

  anchored <- rasch(X, id = rep("P1", nrow(X)), anchors = data.frame(
    item = "I1", k = NA_real_, tau = 0.25
  ))
  expect_false(anchored$est$cluster_inference)
  expect_equal(anchored$items$location[anchored$items$item == "I1"],
               0.25, tolerance = 1e-12)
  expect_equal(anchored$items$se[anchored$items$item == "I1"], 0)
  expect_true(all(is.na(anchored$items$se[anchored$items$item != "I1"])))
  # The free thresholds within a location-anchored polytomous item remain
  # estimated quantities. The exact zero belongs only to their fixed mean.
  expect_true(all(is.na(anchored$thresholds$se[
    anchored$thresholds$item == 1L
  ])))

  bank <- anchored$items[, c("item", "location", "se", "max")]
  bank$se[] <- 0
  linked <- equate_tests(anchored, bank, shift = "none")
  expect_identical(linked$n_common, 5L)
  expect_identical(linked$n, 1L)
  expect_true(all(is.na(linked$table$p_adj)))
  expect_match(linked$note,
               "unavailable or unusable locations or uncertainty: I2, I3, I4, I5")

  # The unsupported fit contributes no estimated variance to its fixed item
  # mean. That one contrast can still use uncertainty supplied by an
  # independent reference; free items from the unsupported fit cannot.
  uncertain_bank <- bank[bank$item == "I1", , drop = FALSE]
  uncertain_bank$location <- uncertain_bank$location - 0.2
  uncertain_bank$se <- 0.1
  exact_test <- equate_tests(anchored, uncertain_bank, shift = "none")
  expect_true(exact_test$inferential)
  expect_equal(exact_test$table$se_diff, 0.1)
  expect_equal(exact_test$table$df, Inf)
  expect_equal(exact_test$table$p,
               2 * pnorm(-abs(exact_test$table$t)), tolerance = 1e-12)
})

test_that("ordinary equating uses finite person-cluster degrees of freedom", {
  make_fit <- function(seed) {
    set.seed(seed)
    G <- 12L; repeats <- 10L; L <- 5L
    theta <- rep(rnorm(G), each = repeats)
    X <- matrix(rbinom(length(theta) * L, 1L,
      plogis(outer(theta, seq(-1, 1, length.out = L), "-"))),
      length(theta), L,
      dimnames = list(NULL, paste0("I", seq_len(L))))
    rasch(X, id = rep(sprintf("P%02d", seq_len(G)), each = repeats))
  }
  fit1 <- make_fit(504)
  fit2 <- make_fit(505)
  expect_true(fit1$est$cluster_inference)
  expect_true(fit2$est$cluster_inference)

  fixed <- equate_tests(fit1, fit2, shift = "none", independent = TRUE)
  v1 <- fit1$items$se^2
  v2 <- fit2$items$se^2
  expected_df <- (v1 + v2)^2 / (v1^2 / 11 + v2^2 / 11)
  expected_t <- (fit1$items$location - fit2$items$location) /
    sqrt(v1 + v2)
  expect_equal(fixed$table$df, expected_df, tolerance = 1e-10)
  expect_equal(fixed$table$t, expected_t, tolerance = 1e-10)
  expect_equal(fixed$table$p,
               2 * pt(-abs(expected_t), df = expected_df),
               tolerance = 1e-12)
  expect_equal(fixed$table$p_adj,
               p.adjust(fixed$table$p, "holm", n = nrow(fixed$table)))
  expect_match(fixed$note, "Welch-Satterthwaite")

  linked <- equate_tests(fit1, fit2, independent = TRUE)
  w <- rasch:::.inverse_variance_weights(v1 + v2)
  u <- w / sum(w)
  Hc <- diag(length(u)) - matrix(u, length(u), length(u), byrow = TRUE)
  S1 <- rasch:::.equate_loc_cov(fit1, fit1$items$item)
  S2 <- rasch:::.equate_loc_cov(fit2, fit2$items$item)
  cv1 <- pmax(diag(Hc %*% S1 %*% t(Hc)), 0)
  cv2 <- pmax(diag(Hc %*% S2 %*% t(Hc)), 0)
  expected_link_df <- (cv1 + cv2)^2 / (cv1^2 / 11 + cv2^2 / 11)
  expect_equal(linked$table$df, expected_link_df, tolerance = 1e-9)
  expect_equal(linked$table$p,
               2 * pt(-abs(linked$table$t), df = expected_link_df),
               tolerance = 1e-12)

  bank <- fit2$items[, c("item", "location", "se", "max")]
  attr(bank, "cov_location") <- S2
  attr(bank, "df_location") <- 7
  bank_link <- equate_tests(fit1, bank, independent = TRUE)
  bank_df <- (cv1 + cv2)^2 / (cv1^2 / 11 + cv2^2 / 7)
  expect_equal(bank_link$table$df, bank_df, tolerance = 1e-9)
  bad_bank <- bank
  attr(bad_bank, "df_location") <- 0
  expect_error(equate_tests(fit1, bad_bank), "positive numeric degree")

  half <- rasch:::.equate_interval_halfwidth(data.frame(
    se_diff = c(1, 1, NA_real_), df = c(9, Inf, 9)
  ))
  expect_equal(half[1:2], c(qt(0.975, 9), qnorm(0.975)),
               tolerance = 1e-12)
  expect_true(is.na(half[3]))
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  expect_no_error(plot_equate(fit1, fit2, shift = "none",
                              independent = TRUE))
})

test_that("uninformative rows cannot manufacture person-cluster support", {
  make_data <- function(G, seed) {
    set.seed(seed)
    repeats <- 4L; L <- 5L
    theta <- rep(rnorm(G), each = repeats)
    X <- matrix(rbinom(length(theta) * L, 1L,
      plogis(outer(theta, seq(-1, 1, length.out = L), "-"))),
      length(theta), L,
      dimnames = list(NULL, paste0("I", seq_len(L))))
    list(X = X,
         id = rep(sprintf("P%02d", seq_len(G)), each = repeats))
  }
  augment <- function(z, n = 20L) {
    list(X = rbind(z$X, matrix(NA_integer_, n, ncol(z$X),
                              dimnames = list(NULL, colnames(z$X)))),
         id = c(z$id, sprintf("empty%02d", seq_len(n))))
  }

  supported <- make_data(12L, 506)
  extra <- augment(supported)
  f0 <- rasch(supported$X, id = supported$id)
  f1 <- rasch(extra$X, id = extra$id)
  expect_true(f0$est$cluster_inference)
  expect_true(f1$est$cluster_inference)
  expect_equal(f1$est$cluster_support, f0$est$cluster_support)
  expect_equal(f1$est$cov_tau, f0$est$cov_tau, tolerance = 1e-12)

  unsupported <- make_data(5L, 507)
  extra_small <- augment(unsupported, 40L)
  s0 <- rasch(unsupported$X, id = unsupported$id)
  s1 <- rasch(extra_small$X, id = extra_small$id)
  expect_false(s0$est$cluster_inference)
  expect_false(s1$est$cluster_inference)
  expect_identical(s1$est$cluster_support$n, 5L)
  expect_equal(s1$est$cluster_support$effective,
               s0$est$cluster_support$effective, tolerance = 1e-12)
  expect_true(all(is.na(s1$est$cov_tau)))

  # Repeated IDs only on non-contributing rows do not turn an ordinary
  # row-independent covariance into a finite-df clustered covariance.
  set.seed(508)
  X <- matrix(rbinom(500L, 1L, 0.5), 100L, 5L,
              dimnames = list(NULL, paste0("I", 1:5)))
  ordinary <- rasch(X, id = sprintf("U%03d", seq_len(nrow(X))))
  malformed <- ordinary
  malformed$est$cluster_support$repeated <- NULL
  expect_true(is.na(rasch:::.equate_cov_df(malformed)))
  X_empty <- rbind(X, matrix(NA_integer_, 2L, ncol(X),
                             dimnames = list(NULL, colnames(X))))
  empty_repeat <- rasch(X_empty,
    id = c(sprintf("U%03d", seq_len(nrow(X))), "empty", "empty"))
  expect_false(empty_repeat$est$cluster_support$repeated)
  expect_false(empty_repeat$repeated_ids)
  expect_equal(empty_repeat$est$cov_tau, ordinary$est$cov_tau,
               tolerance = 1e-12)
  expect_equal(empty_repeat$total_chisq_p, ordinary$total_chisq_p,
               tolerance = 1e-12)
  expect_equal(empty_repeat$items$p_adj, ordinary$items$p_adj,
               tolerance = 1e-12)
  bank <- ordinary$items[, c("item", "location", "se", "max")]
  bank$location <- bank$location + 0.1
  linked <- equate_tests(empty_repeat, bank, shift = "none")
  expect_true(all(linked$table$df == Inf))
})

test_that("independent-person sandwich needs support for every fitted direction", {
  set.seed(509)
  N <- 16L; L <- 30L
  X <- vapply(seq_len(L), function(j)
    sample(rep(0:1, each = N / 2L)), integer(N))
  colnames(X) <- paste0("I", seq_len(L))

  fit <- rasch(X)
  expect_true(fit$est$converged)
  expect_false(fit$est$cluster_support$repeated)
  expect_identical(fit$est$cluster_support$n, N)
  expect_false(fit$est$cluster_inference)
  expect_true(all(is.na(fit$est$cov_beta)))
  expect_true(all(is.na(fit$est$cov_tau)))
  expect_true(all(is.na(fit$items$se)))
  expect_match(paste(fit$notes, collapse = " "),
               "independent person units.*rank-deficient")
  bank <- fit$items[, c("item", "location", "se", "max")]
  bank$se[] <- 0
  linked <- equate_tests(fit, bank, shift = "none")
  expect_false(linked$inferential)
  expect_true(all(is.na(linked$table$p_adj)))

  low <- pcml(X)
  expect_false(low$cluster_inference)
  expect_true(all(is.na(low$cov_beta)))
  expect_true(all(is.na(low$cov_tau)))

  set.seed(510)
  adequate <- matrix(rbinom(1200L, 1L, 0.5), 120L, 10L,
                     dimnames = list(NULL, paste0("A", 1:10)))
  supported <- rasch(adequate)
  expect_false(supported$est$cluster_support$repeated)
  expect_true(supported$est$cluster_inference)
  expect_true(all(is.finite(supported$items$se)))
})

test_that("EFRM stage-one covariance needs independent-person support", {
  set.seed(511)
  N <- 16L; L <- 30L
  X <- vapply(seq_len(L), function(j)
    sample(rep(0:1, each = N / 2L)), integer(N))
  colnames(X) <- paste0("I", seq_len(L))
  dat <- data.frame(X, g = "all")

  descriptive <- rasch_efrm(dat, item_sets = list(core = colnames(X)),
                            groups = "g", boot_reps = 0)
  expect_true(descriptive$est$stage1_converged)
  expect_false(descriptive$est$cluster_inference)
  expect_identical(descriptive$est$cluster_support$n, N)
  expect_true(all(is.na(descriptive$est$cov_tau)))
  expect_true(all(is.na(descriptive$items$se)))
  expect_true(all(is.na(descriptive$unit_cov$cov_dtilde)))
  expect_match(paste(descriptive$notes, collapse = " "),
               "EFRM conditional calibration.*rank-deficient")

  set.seed(512)
  adequate_X <- matrix(rbinom(1200L, 1L, 0.5), 120L, 10L,
                       dimnames = list(NULL, paste0("A", 1:10)))
  adequate <- rasch_efrm(data.frame(adequate_X, g = "all"),
                         item_sets = list(core = colnames(adequate_X)),
                         groups = "g", boot_reps = 0)
  expect_true(adequate$est$cluster_inference)
  expect_true(all(is.finite(adequate$items$se)))
  expect_true(all(is.finite(adequate$unit_cov$cov_dtilde)))

  two_sets <- list(A = colnames(X)[seq_len(L / 2L)],
                   B = colnames(X)[L / 2L + seq_len(L / 2L)])
  for (method in c("hybrid", "bootstrap"))
    expect_error(rasch_efrm(dat, item_sets = two_sets, groups = "g",
                            se_method = method, boot_reps = 30, workers = 1),
                 "insufficient independent-person support.*stage-one covariance")
})

test_that("a single anchor changes only the origin and propagates downstream", {
  set.seed(601)
  N <- 500L
  theta <- rnorm(N)
  truth <- list(c(-1, 0.2), c(-0.6, 0.5), c(-0.2, 0.9),
                c(0.2, 1.2), c(0.5, 1.5))
  draw_item <- function(th, tau) {
    score <- 0:length(tau)
    prob <- exp(score * th - c(0, cumsum(tau)))
    sample(score, 1L, prob = prob)
  }
  X <- sapply(truth, function(tau)
    vapply(theta, draw_item, integer(1), tau = tau))
  colnames(X) <- paste0("I", 1:5)
  # Separate refits use different coordinates. Fit tightly, but allow for
  # numerical differences across BLAS/LAPACK builds when comparing them.
  # Exact anchor values are checked separately at stricter tolerance.
  fit_tol <- 1e-12
  comparison_tol <- 1e-8
  free <- rasch(X, tol = fit_tol)
  expect_true(free$est$converged)
  shift <- 0.8
  anchored <- rasch(X, anchors = data.frame(
    item = "I2", k = NA_real_, tau = free$items$location[2] + shift
  ), tol = fit_tol)
  expect_true(anchored$est$converged)
  expect_equal(anchored$items$location[2], free$items$location[2] + shift,
               tolerance = 1e-12)
  expect_equal(unlist(anchored$tau_list),
               unlist(free$tau_list) + shift, tolerance = comparison_tol)
  ok <- is.finite(free$person$theta) & is.finite(anchored$person$theta)
  expect_equal(anchored$person$theta[ok], free$person$theta[ok] + shift,
               tolerance = comparison_tol)
  expect_equal(anchored$person$se, free$person$se, tolerance = comparison_tol)
  expect_equal(anchored$residuals, free$residuals, tolerance = comparison_tol)
  q <- stats::setNames(seq(0.5, 1.5, length.out = ncol(X)), colnames(X))
  weighted_free <- weighted_person_estimates(free, q)
  weighted_anchored <- weighted_person_estimates(anchored, q)
  weighted_ok <- is.finite(weighted_free$theta) &
    is.finite(weighted_anchored$theta)
  expect_equal(weighted_anchored$theta[weighted_ok],
               weighted_free$theta[weighted_ok] + shift,
               tolerance = comparison_tol)
  expect_equal(weighted_anchored$se, weighted_free$se,
               tolerance = comparison_tol)
  expect_equal(weighted_anchored$weighted_score,
               weighted_free$weighted_score, tolerance = 1e-12)
  M <- nrow(free$thresholds)
  anchor_mean <- numeric(M)
  rows <- which(free$thresholds$item == 2L)
  anchor_mean[rows] <- 1 / length(rows)
  H <- diag(M) - matrix(1, M, 1L) %*% t(anchor_mean)
  expect_equal(anchored$est$cov_tau,
               H %*% free$est$cov_tau %*% t(H), tolerance = comparison_tol)

  # Fixing one numerical threshold also identifies only the origin. Its
  # covariance is the free covariance translated by that threshold, while
  # person measures move with the item scale and fitted residuals do not move.
  anchor_row <- which(free$thresholds$item == 3L &
                        free$thresholds$k == 1L)
  threshold_anchored <- rasch(X, anchors = data.frame(
    item = "I3", k = 1L,
    tau = free$thresholds$tau[anchor_row] + shift
  ), tol = fit_tol)
  expect_true(threshold_anchored$est$converged)
  expect_equal(threshold_anchored$thresholds$tau[anchor_row],
               free$thresholds$tau[anchor_row] + shift, tolerance = 1e-12)
  H_threshold <- diag(M)
  H_threshold <- H_threshold -
    matrix(1, M, 1L) %*% H_threshold[anchor_row, , drop = FALSE]
  expect_equal(threshold_anchored$thresholds$tau,
               free$thresholds$tau + shift, tolerance = comparison_tol)
  expect_equal(threshold_anchored$est$cov_tau,
               H_threshold %*% free$est$cov_tau %*% t(H_threshold),
               tolerance = comparison_tol)
  ok_threshold <- is.finite(free$person$theta) &
    is.finite(threshold_anchored$person$theta)
  expect_equal(threshold_anchored$person$theta[ok_threshold],
               free$person$theta[ok_threshold] + shift,
               tolerance = comparison_tol)
  expect_equal(threshold_anchored$person$se, free$person$se,
               tolerance = comparison_tol)
  expect_equal(threshold_anchored$residuals, free$residuals,
               tolerance = comparison_tol)

  d <- simulate_btl(6, 20, 5, seed = 602)
  b0 <- btl(d, "object_a", "object_b", "winner")
  beta0 <- setNames(b0$objects$location, b0$objects$object)
  target <- setNames(beta0[["O2"]] + 0.7, "O2")
  b1 <- btl(d, "object_a", "object_b", "winner", anchors = target)
  beta1 <- setNames(b1$objects$location, b1$objects$object)
  bshift <- target[[1L]] - beta0[["O2"]]
  expect_equal(beta1, beta0 + bshift, tolerance = 1e-10)
  expect_equal(b1$fitted_prob, b0$fitted_prob, tolerance = 1e-10)
  e <- match("O2", rownames(b0$cov_beta))
  Hb <- diag(nrow(b0$cov_beta)) -
    matrix(1, nrow(b0$cov_beta), 1L) %*%
    diag(nrow(b0$cov_beta))[e, , drop = FALSE]
  expect_equal(unname(b1$cov_beta),
               unname(Hb %*% b0$cov_beta %*% t(Hb)), tolerance = 1e-10)
  expect_equal(b1$objects$se[b1$objects$object == "O2"], 0)
})

test_that("DIF pairing follows residual rather than calibration support", {
  d <- simulate_rasch(120, 6, n_categories = 3, seed = 993)
  item_names <- grep("^I", names(d), value = TRUE)
  X <- as.matrix(d[item_names])
  n_second <- 30L
  second <- matrix(NA_integer_, n_second, ncol(X),
                   dimnames = list(NULL, colnames(X)))
  set.seed(994)
  second[, 1L] <- sample(0:2, n_second, replace = TRUE)
  dat <- data.frame(
    rbind(X, second),
    occasion = c(rep("T1", nrow(X)), rep("T2", n_second)))
  fit <- rasch(dat, factors = "occasion",
               id = c(d$id, d$id[seq_len(n_second)]))

  # These second occasions carry fitted residuals but no conditional item
  # pair. They define the DIF error stratum without changing the calibration
  # covariance's independent-unit interpretation.
  expect_false(fit$repeated_ids)
  expect_true(fit$repeated_residual_ids)
  out <- dif_anova(fit)
  expect_identical(out$within, "occasion")
  expect_true(any(out$summary$item == item_names[1L] &
                  out$summary$term == "occasion"))
})

test_that("paired DIF contrasts retain every planned between-person stratum", {
  person <- rep(paste0("P", 1:4), each = 2L)
  factors <- data.frame(
    sex = rep(c("A", "A", "B", "B"), each = 2L),
    occasion = rep(c("T1", "T2"), 4L))
  grp <- interaction(factors$sex, factors$occasion, sep = ":",
                     drop = TRUE)
  cellmap <- unique(data.frame(cell = as.character(grp), factors,
                               stringsAsFactors = FALSE))
  cellmap <- cellmap[match(levels(grp), cellmap$cell), , drop = FALSE]
  weights <- c("A:T1" = -0.5, "B:T1" = -0.5,
               "A:T2" = 0.5, "B:T2" = 0.5)
  z <- c(0.0, 0.4, 0.1, 0.6, -0.2, 0.1, 0.2, 0.7)

  complete <- rasch:::.dif_paired_cell_contrast(
    z, factors, grp, person, "occasion", cellmap, weights)
  expect_true(is.list(complete))

  # Removing one B person's second occasion leaves the B contrast mean but
  # not its sampling variance. Silently dropping B would change the planned
  # equal-sex contrast to the A-only contrast.
  z[8L] <- NA_real_
  expect_null(rasch:::.dif_paired_cell_contrast(
    z, factors, grp, person, "occasion", cellmap, weights))
})
