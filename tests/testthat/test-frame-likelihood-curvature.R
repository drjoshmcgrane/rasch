test_that("likelihood curvature distinguishes maxima, saddles and flat directions", {
  good <- -matrix(c(2, 0.5, 0.5, 1), 2L)
  saddle <- matrix(c(-1, 2, 2, -1), 2L)
  for (sc in list(c(1, 1), c(1e-80, 1e80), c(7, 0.03))) {
    D <- diag(sc)
    expect_true(.likelihood_curvature_ok(D %*% good %*% D))
    expect_false(.likelihood_curvature_ok(D %*% saddle %*% D))
    expect_false(.likelihood_curvature_ok(D %*% diag(c(-1, 0)) %*% D))
  }
  expect_false(.likelihood_curvature_ok(matrix(NA_real_, 2L, 2L)))
  expect_false(.likelihood_curvature_ok(diag(2L)))
  expect_false(.likelihood_curvature_ok(matrix(c(-1, 0, 1, -1), 2L)))
})

test_that("EFRM refuses a stationary saddle and excludes it from bootstrap draws", {
  set.seed(1)
  ability <- rnorm(300)
  difficulties <- rep(c(-0.5, 0, 0.5), each = 2) + rep(c(-1.5, 1.5), 3)
  X <- sapply(difficulties, function(d) rbinom(300, 1, plogis(ability - d)))
  colnames(X) <- paste0("I", 1:6)
  d <- data.frame(rbind(X, X[, c(2, 1, 4, 3, 6, 5)]),
                  group = rep(c("a", "b"), each = 300))
  original <- .efrm_solve
  curvature <- .likelihood_curvature_ok
  H <- NULL
  sol <- NULL
  inject <- FALSE
  calls <- 0L
  draws <- list()
  testthat::local_mocked_bindings(
    .likelihood_curvature_ok = function(h, ...) {
      H <<- h
      curvature(h, ...)
    },
    .efrm_solve = function(...) {
      calls <<- calls + 1L
      z <- if (inject && calls > 1L && calls <= 6L) sol else original(...)
      if (!inject) sol <<- z else draws[[calls]] <<- z
      z
    }, .package = "rasch")
  expect_error(rasch_efrm(d, list(all = colnames(X)), "group", boot_reps = 0),
               "did not reach an identified local maximum")
  expect_false(sol$maximum_ok)
  expect_false(sol$converged)
  expect_true(sol$cluster_inference) # This is not a cluster-support failure.
  expect_true(all(is.na(sol$cov_joint)))
  expect_true(all(is.na(sol$se_log_phi)))

  # Independent binary pair-conditional likelihood, in the solver's free
  # threshold coordinates and reference-group unit parameterisation.
  A <- rbind(diag(5L), rep(-1, 5L))
  z <- c((sol$dtilde * sol$phi[1L])[1:5], sol$phi[2L] / sol$phi[1L])
  ll <- function(z) {
    delta <- drop(A %*% z[1:5])
    phi <- ifelse(d$group == "a", 1, z[6L])
    out <- 0
    for (i in 1:5) for (j in (i + 1L):6) {
      use <- d[[i]] + d[[j]] == 1
      eta <- phi[use] * (delta[j] - delta[i])
      yy <- d[[i]][use]
      out <- out + sum(yy * plogis(eta, log.p = TRUE) +
                        (1 - yy) * plogis(-eta, log.p = TRUE))
    }
    out
  }
  expect_equal(ll(z), sol$loglik, tolerance = 1e-10)
  expect_equal(unname(H), unname(stats::optimHess(z, ll)), tolerance = 1e-5)
  grad <- vapply(seq_along(z), function(j) {
    dz <- numeric(length(z)); dz[j] <- 1e-5
    (ll(z + dz) - ll(z - dz)) / 2e-5
  }, 0)
  expect_lt(max(abs(grad)), 1e-3)
  e <- eigen(H, symmetric = TRUE)
  expect_gt(e$values[1L], 1)
  gains <- vapply(c(-0.1, 0.1), function(h)
    ll(z + h * e$vectors[, 1L]) - ll(z), 0)
  expect_gt(max(gains), 0.1)

  # Feed the actual rejected solve to five otherwise ordinary bootstrap
  # replicates. It must be excluded before any transformation or covariance.
  inject <- TRUE
  calls <- 0L
  ordinary <- simulate_efrm(150, 5, n_sets = 1, n_groups = 2, seed = 9207)
  f <- rasch_efrm(ordinary, attr(ordinary, "truth")$item_sets, "group",
                  id = "id", se_method = "bootstrap", boot_reps = 40,
                  workers = 1, seed = 889)
  expect_true(f$est$converged)
  expect_equal(f$full_boot_reps_used, 35L)
  expect_equal(f$full_boot_reps_failed, 5L)
  keep <- draws[-seq_len(6L)]
  C <- cov(do.call(rbind, lapply(keep, function(x) c(x$dtilde, log(x$phi)))))
  expect_equal(unname(f$unit_cov$cov_joint), unname(C), tolerance = 1e-12)
})

test_that("BTL-EFRM rejects a full-Fisher-rank saddle using the exact Hessian", {
  set.seed(1)
  obj <- paste0("O", 1:6)
  beta <- rep(c(-0.5, 0, 0.5), each = 2) + rep(c(-1.5, 1.5), 3)
  pr <- t(combn(1:6, 2L))
  ix <- pr[rep(seq_len(nrow(pr)), each = 3L), , drop = FALSE]
  aa <- do.call(rbind, lapply(seq_len(40L), function(j) {
    y <- rbinom(nrow(ix), 1L, plogis(beta[ix[, 1L]] - beta[ix[, 2L]]))
    data.frame(a = ix[, 1L], b = ix[, 2L], y = y, judge = paste0("a", j),
               panel = "a")
  }))
  bb <- aa
  swap <- c(2, 1, 4, 3, 6, 5)
  bb$a <- swap[aa$a]; bb$b <- swap[aa$b]
  bb$judge <- sub("^a", "b", aa$judge); bb$panel <- "b"
  d <- rbind(aa, bb)
  H <- NULL
  original <- .likelihood_curvature_ok
  testthat::local_mocked_bindings(.likelihood_curvature_ok = function(h, ...) {
    H <<- h
    original(h, ...)
  }, .package = "rasch")
  s <- .btlef_stage1(d$a, d$b, d$y, d$panel, d$judge, 6L, 100L, 1e-7)
  expect_true(s$rank_ok)
  expect_true(s$cluster_ok)
  expect_false(s$maximum_ok)
  expect_false(s$converged)
  expect_true(all(is.na(s$cov_lrho)))
  expect_true(all(is.na(s$influence_lrho)))
  expect_true(all(is.na(s$se_beta)))
  A <- rbind(diag(5L), rep(-1, 5L))
  z <- c(s$beta[1:5], log(s$rho["b"]))
  ll <- function(z) {
    b <- drop(A %*% z[1:5])
    eta <- ifelse(d$panel == "a", 1, exp(z[6L])) * (b[d$a] - b[d$b])
    sum(d$y * plogis(eta, log.p = TRUE) +
          (1 - d$y) * plogis(-eta, log.p = TRUE))
  }
  expect_equal(ll(z), s$ll, tolerance = 1e-10)
  expect_equal(unname(H), unname(stats::optimHess(z, ll)), tolerance = 1e-5)
  grad <- vapply(seq_along(z), function(j) {
    dz <- numeric(length(z)); dz[j] <- 1e-5
    (ll(z + dz) - ll(z - dz)) / 2e-5
  }, 0)
  expect_lt(max(abs(grad)), 1e-3)
  e <- eigen(H, symmetric = TRUE)
  expect_gt(e$values[1L], 1)
  gains <- vapply(c(-0.1, 0.1), function(h)
    ll(z + h * e$vectors[, 1L]) - ll(z), 0)
  expect_gt(max(gains), 0.1)
  public <- data.frame(object_a = obj[d$a], object_b = obj[d$b],
    winner = obj[ifelse(d$y == 1, d$a, d$b)], judge = d$judge, panel = d$panel)
  for (method in c("conditional", "judge_bootstrap", "bootstrap"))
    expect_error(btl_efrm(public, "object_a", "object_b", "winner", "judge",
      panels = "panel", object_sets = list(all = obj), se_method = method,
      boot_reps = 30, workers = 1, seed = 71), "local-maximum curvature check")

  # Fixing the panel units removes the nonlinear coordinates. This is an
  # ordinary concave BTL likelihood and must remain estimable.
  fixed <- .btlef_stage1(d$a, d$b, d$y, d$panel, d$judge, 6L, 100L, 1e-7,
                         rho_fixed = c(a = 1, b = 1))
  expect_true(fixed$maximum_ok)
  expect_true(fixed$converged)
  expect_true(all(is.finite(fixed$se_beta)))
})

test_that("BTL frame curvature includes every nonlinear term away from convergence", {
  set.seed(8042)
  pairs <- t(combn(1:4, 2L))
  ix <- pairs[rep(seq_len(nrow(pairs)), 60L), , drop = FALSE]
  panel <- rep(c("a", "b", "c"), each = nrow(ix) / 3L)
  rho <- c(a = 0.65, b = 1, c = 1.5)[panel]
  beta <- c(-1.2, -0.4, 0.5, 1.1)
  y <- rbinom(nrow(ix), 1, plogis(rho * (beta[ix[, 1L]] - beta[ix[, 2L]])))
  H <- NULL
  original <- .likelihood_curvature_ok
  testthat::local_mocked_bindings(.likelihood_curvature_ok = function(h, ...) {
    H <<- h
    original(h, ...)
  }, .package = "rasch")
  s <- .btlef_stage1(ix[, 1L], ix[, 2L], y, panel,
                     paste0("j", rep(1:30, each = 12L)), 4L, 1L, 1e-7)
  expect_false(s$converged)
  A <- rbind(diag(3L), rep(-1, 3L))
  z <- c(s$beta[1:3], log(s$rho[s$free]))
  eta <- function(z) {
    b <- drop(A %*% z[1:3])
    r <- setNames(rep(1, length(s$panels)), s$panels)
    r[s$free] <- exp(z[3L + seq_along(s$free)])
    r[panel] * (b[ix[, 1L]] - b[ix[, 2L]])
  }
  ll <- function(z) sum(y * plogis(eta(z), log.p = TRUE) +
                         (1 - y) * plogis(-eta(z), log.p = TRUE))
  expect_equal(unname(H), unname(stats::optimHess(z, ll)), tolerance = 1e-5)
  # At this unfinished iterate the log-unit diagonal correction is nonzero;
  # a test only at an optimum would not pin that second derivative.
  for (h in seq_along(s$free)) {
    sel <- panel == s$free[h]
    correction <- sum((y[sel] - plogis(eta(z)[sel])) * eta(z)[sel])
    expect_gt(abs(correction), 0.1)
  }
})
