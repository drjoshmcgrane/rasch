test_that("full EFRM bootstrap excludes unidentified units from every covariance", {
  d <- simulate_efrm(150, 5, n_sets = 1, n_groups = 2, seed = 9207)
  sets <- attr(d, "truth")$item_sets
  original <- rasch:::.efrm_solve
  solves <- list()
  n_bad <- 5L
  bad_observed <- FALSE
  testthat::local_mocked_bindings(.efrm_solve = function(...) {
    z <- original(...)
    r <- length(solves)
    if (bad_observed || (r > 0L && r <= n_bad)) z$phi_unident[] <- TRUE
    solves[[r + 1L]] <<- z
    z
  }, .package = "rasch")
  run_fit <- function() rasch_efrm(
    d, sets, "group", id = "id", se_method = "bootstrap",
    boot_reps = 40, workers = 1, seed = 889)

  expect_no_warning(f <- run_fit())
  expect_length(solves, 41L)
  expect_true(all(vapply(solves, function(z) isTRUE(z$converged), TRUE)))
  expect_false(any(solves[[1L]]$phi_unident))
  draws <- solves[-1L]
  keep <- !vapply(draws, function(z) any(z$phi_unident), TRUE)
  expect_equal(sum(!keep), 5L)
  expect_identical(f$se_method, "bootstrap")
  expect_identical(f$full_boot_reps_used, 35L)
  expect_identical(f$full_boot_reps_failed, 5L)
  expect_equal(f$boot_reps_used, 35L)
  expect_equal(f$boot_reps_failed, 5L)
  accepted <- do.call(rbind, lapply(draws[keep], function(z)
    c(z$dtilde, log(z$phi))))
  C <- cov(accepted)
  nd <- length(solves[[1L]]$dtilde)
  g <- nd + seq_along(solves[[1L]]$phi)
  expect_equal(unname(f$unit_cov$cov_joint), unname(C), tolerance = 1e-12)
  expect_equal(unname(f$unit_cov$cov_dtilde),
               unname(C[seq_len(nd), seq_len(nd)]), tolerance = 1e-12)
  expect_equal(unname(f$unit_cov$cov_log_phi), unname(C[g, g]),
               tolerance = 1e-12)
  expect_equal(unname(f$phi_table$se_log_phi), sqrt(unname(diag(C)[g])),
               tolerance = 1e-12)

  # Too few identified draws must not produce a bootstrap covariance.
  n_bad <- 40L
  solves <- list()
  expect_warning(fallback <- run_fit(), "0 of 40 replicates were usable")
  expect_identical(fallback$se_method, "hybrid")
  # With one item set there is no linking stage; the hybrid fallback is
  # exactly the conditional stage-one covariance.
  expect_identical(fallback$unit_cov$method, "stage-one")
  expect_identical(fallback$full_boot_reps_requested, 40L)
  expect_identical(fallback$full_boot_reps_attempted, 40L)
  expect_identical(fallback$full_boot_reps_used, 0L)
  expect_identical(fallback$full_boot_reps_failed, 40L)
  expect_equal(fallback$unit_cov$cov_joint, solves[[1L]]$cov_joint)
  expect_true(all(is.finite(fallback$phi_table$se_log_phi)))

  # The same structural flag already refuses the observed estimator.
  bad_observed <- TRUE
  expect_error(rasch_efrm(d, sets, "group", id = "id", boot_reps = 0),
               "are unidentified")
})

test_that("bootstrap point estimates do not require per-draw analytic SEs", {
  d <- simulate_efrm(150, 5, n_sets = 1, n_groups = 2, seed = 9207)
  original <- rasch:::.efrm_solve
  calls <- 0L
  testthat::local_mocked_bindings(.efrm_solve = function(...) {
    z <- original(...)
    calls <<- calls + 1L
    if (calls > 1L) {
      z$se_log_phi[] <- NA_real_
      z$cov_joint[,] <- z$cov_log_phi[,] <- z$cov_dtilde[,] <- NA_real_
      z$cluster_inference <- FALSE
    }
    z
  }, .package = "rasch")
  expect_no_warning(f <- rasch_efrm(
    d, attr(d, "truth")$item_sets, "group", id = "id",
    se_method = "bootstrap", boot_reps = 40, workers = 1, seed = 889))
  expect_identical(f$se_method, "bootstrap")
  expect_identical(f$full_boot_reps_used, 40L)
  expect_identical(f$full_boot_reps_failed, 0L)
  expect_true(all(is.finite(f$unit_cov$cov_joint)))
  expect_true(all(is.finite(f$phi_table$se_log_phi)))
})

test_that("actual unidentified small-group resamples are counted as failures", {
  d <- simulate_efrm(120, 5, n_sets = 1, n_groups = 2,
                     group_unit_ratio = 1, seed = 2)
  sets <- attr(d, "truth")$item_sets
  d <- d[c(1:120, 121:123), ]
  original <- rasch:::.efrm_solve
  solves <- list()
  testthat::local_mocked_bindings(.efrm_solve = function(...) {
    z <- original(...)
    solves[[length(solves) + 1L]] <<- z
    z
  }, .package = "rasch")
  expect_no_warning(f <- rasch_efrm(
    d, sets, "group", id = "id", se_method = "bootstrap",
    boot_reps = 40, workers = 1, seed = 719))
  expect_length(solves, 41L)
  expect_false(any(solves[[1L]]$phi_unident))
  draws <- solves[-1L]
  unidentified <- vapply(draws, function(z)
    isTRUE(z$converged) && any(z$phi_unident), TRUE)
  expect_gt(sum(unidentified), 0L)
  keep <- vapply(draws, function(z)
    isTRUE(z$converged) && !any(z$phi_unident) &&
      all(is.finite(c(z$dtilde, log(z$phi)))), TRUE)
  expect_identical(f$se_method, "bootstrap")
  expect_equal(f$full_boot_reps_used, sum(keep))
  expect_equal(f$full_boot_reps_failed, sum(!keep))
  C <- cov(do.call(rbind, lapply(draws[keep], function(z)
    c(z$dtilde, log(z$phi)))))
  expect_equal(unname(f$unit_cov$cov_joint), unname(C), tolerance = 1e-12)
  expect_false(f$unit_support$phi_inference)
})
