test_that("expected BTL residuals use common point totals and both corrections", {
  # Keep this reference deliberately independent of the package helper.  The
  # fitted value is an expected number of points for object_a, not a second
  # row count; both matrices therefore use the same pooled comparison total.
  expected_residual <- function(ia, ib, resp, w, K, expected_a, m = 1L) {
    S <- matrix(0, K, K)
    E <- matrix(0, K, K)
    for (r in seq_along(ia)) {
      S[ia[r], ib[r]] <- S[ia[r], ib[r]] + w[r] * resp[r]
      S[ib[r], ia[r]] <- S[ib[r], ia[r]] + w[r] * (m - resp[r])
      E[ia[r], ib[r]] <- E[ia[r], ib[r]] + w[r] * expected_a[r]
      E[ib[r], ia[r]] <- E[ib[r], ia[r]] + w[r] * (m - expected_a[r])
    }
    total <- S + t(S)
    observed <- qlogis((S + 0.5) / (total + 1))
    fitted <- qlogis((E + 0.5) / (total + 1))
    R <- observed - fitted
    R[total == 0] <- 0
    (R - t(R)) / 2
  }

  ia <- c(1L, 1L, 2L, 3L)
  ib <- c(2L, 3L, 3L, 4L)
  resp <- c(4, 0, 3, 1)
  w <- c(2, 3, 4, 5)
  expected_a <- c(3.2, 0.8, 2.4, 1.6)
  got <- .btl_resid_matrix_expected(ia, ib, resp, w, 4L,
                                    expected_a, m = 4L)
  expect_equal(got, expected_residual(ia, ib, resp, w, 4L,
                                      expected_a, m = 4L))

  # If the fitted expected points equal the observed points, continuity is
  # applied symmetrically and the residual is exactly zero, including at
  # ordered-response boundaries.
  expect_equal(
    .btl_resid_matrix_expected(ia, ib, resp, w, 4L, resp, m = 4L),
    matrix(0, 4L, 4L))
})

test_that("BTL expected residuals are invariant to orientation and row expansion", {
  ia <- c(1L, 1L, 2L)
  ib <- c(2L, 3L, 3L)
  resp <- c(1, 0, 1)
  expected_a <- c(0.65, 0.2, 0.8)
  w <- c(3, 2, 4)
  original <- .btl_resid_matrix_expected(ia, ib, resp, w, 3L,
                                         expected_a, m = 1L)

  reversed <- .btl_resid_matrix_expected(ib, ia, 1 - resp, w, 3L,
                                         1 - expected_a, m = 1L)
  expect_equal(reversed, original)

  ia_exp <- rep(ia, w)
  ib_exp <- rep(ib, w)
  resp_exp <- rep(resp, w)
  expected_exp <- rep(expected_a, w)
  expanded <- .btl_resid_matrix_expected(ia_exp, ib_exp, resp_exp,
                                         rep(1, length(ia_exp)), 3L,
                                         expected_exp, m = 1L)
  expect_equal(expanded, original)
})

test_that("dichotomous expected residual has the symmetric no-effect formula", {
  n <- 4
  ia <- rep(1L, n)
  ib <- rep(2L, n)
  resp <- rep(1, n)
  expected_a <- rep(0.5, n)
  got <- .btl_resid_matrix_expected(ia, ib, resp, rep(1, n), 2L,
                                    expected_a, m = 1L)
  want <- qlogis((n + 0.5) / (n + 1)) -
    qlogis((n * 0.5 + 0.5) / (n + 1))
  expect_equal(got[1, 2], want)
  expect_equal(got[2, 1], -want)
})

test_that("ordinary graded BTL dimensionality compares observed and fitted means", {
  objects <- c("A", "B", "C", "D", "E")
  beta <- setNames(c(-2.2, -0.9, -0.1, 0.8, 2.4), objects)
  tau <- c(-2.4, -0.35, 0.35, 2.4)
  m <- length(tau)
  pairs <- utils::combn(objects, 2L, simplify = FALSE)
  n <- 100000L
  d <- do.call(rbind, lapply(pairs, function(pair) {
    p <- item_moments(beta[pair[1L]] - beta[pair[2L]], tau)$P
    data.frame(a = pair[1L], b = pair[2L], response = 0:m,
               count = pmax(round(n * p), 1L))
  }))
  fit <- btl(d, "a", "b", response = "response", count = "count")
  result <- btl_dimensionality(fit, reps = 20L, seed = 81)
  expect_true(fit$converged)
  expect_lt(max(abs(result$residual_matrix)), 3e-4)

  cmp <- fit$comparisons
  expected <- .btl_fitted_score_moments(fit)$E
  want <- .btl_resid_matrix_expected(
    match(cmp$object_a, fit$objects$object),
    match(cmp$object_b, fit$objects$object), cmp$response, cmp$weight,
    nrow(fit$objects), expected, m = fit$m)
  expect_equal(result$residual_matrix, want)
})

test_that("binary position-adjusted BTL dimensionality uses fitted position means", {
  objects <- c("A", "B", "C", "D", "E")
  beta <- setNames(c(-2.2, -0.9, -0.1, 0.8, 2.4), objects)
  position <- 1.6
  n <- 100000L
  pairs <- utils::combn(objects, 2L, simplify = FALSE)
  d <- do.call(rbind, lapply(pairs, function(pair) do.call(rbind, lapply(
    list(pair, rev(pair)), function(ab) {
      p <- plogis(beta[ab[1L]] - beta[ab[2L]] + position)
      na <- round(n * p)
      data.frame(a = ab[1L], b = ab[2L], winner = c(ab[1L], ab[2L]),
                 count = c(na, n - na))
    }))))
  fit <- btl(d, "a", "b", winner = "winner", count = "count",
             position = TRUE)
  result <- btl_dimensionality(fit, reps = 20L, seed = 82)
  expect_true(fit$converged)
  expect_lt(max(abs(result$residual_matrix)), 3e-4)

  cmp <- fit$comparisons
  expected <- .btl_fitted_moments(fit, cmp)$E
  want <- .btl_resid_matrix_expected(
    match(cmp$object_a, fit$objects$object),
    match(cmp$object_b, fit$objects$object), cmp$response, cmp$weight,
    nrow(fit$objects), expected, m = fit$m)
  expect_equal(result$residual_matrix, want)
})

test_that("history references rebuild expectations for binary and ordered draws", {
  run_history <- function(generated, seed) {
    if ("winner" %in% names(generated)) {
      fit <- btl(generated, "object_a", "object_b", winner = "winner",
                 judge = "judge", order = "order")
    } else {
      fit <- btl(generated, "object_a", "object_b", response = "response",
                 judge = "judge", order = "order")
    }
    calls <- list()
    original <- get(".btl_resid_matrix_expected", asNamespace("rasch"))
    result <- testthat::with_mocked_bindings(
      btl_dimensionality(fit, reps = 20L, seed = seed,
                         independent_comparisons = TRUE),
      .btl_resid_matrix_expected = function(ia, ib, resp, w, K,
                                            expected_a, m = 1L) {
        calls[[length(calls) + 1L]] <<- list(
          ia = ia, ib = ib, resp = resp, w = w, K = K,
          expected_a = expected_a, m = m)
        original(ia, ib, resp, w, K, expected_a, m)
      },
      .package = "rasch")
    list(fit = fit, result = result, captured = calls)
  }

  check_history <- function(generated, seed) {
    run <- run_history(generated, seed)
    fit <- run$fit
    result <- run$result
    captured <- run$captured
    expect_length(captured, 21L) # observed matrix plus 20 replicates
    expect_length(result$reference$draws, 20L)
    expect_true(all(is.finite(result$reference$draws)))

    # The observed call uses fit$comparisons, whereas the history oracle is
    # most naturally rebuilt in fit$dependence_data's sorted row order.
    dd <- fit$dependence_data
    cmp <- fit$comparisons
    dd_key <- paste(dd$judge, dd$order, sep = "\r")
    cmp_key <- paste(cmp$judge, cmp$order, sep = "\r")
    observed_expected <- .btl_fitted_moments(fit, dd)$E
    observed_in_cmp_order <- observed_expected[match(cmp_key, dd_key)]
    expect_equal(captured[[1L]]$expected_a, observed_in_cmp_order)
    want <- .btl_resid_matrix_expected(
      match(cmp$object_a, fit$objects$object),
      match(cmp$object_b, fit$objects$object), cmp$response, cmp$weight,
      nrow(fit$objects), observed_in_cmp_order, m = fit$m)
    expect_equal(result$residual_matrix, want)

    # Rebuild every simulated row expectation from the generated response
    # history. This is an independent oracle for the nested sequential draw;
    # using the observed history would fail because the prior verdicts change.
    beta <- setNames(fit$objects$location, fit$objects$object)
    dep <- setNames(fit$dependence$estimate, fit$dependence$effect)
    simulated <- captured[-1L]
    differs_from_observed <- logical(length(simulated))
    for (i in seq_along(simulated)) {
      draw <- simulated[[i]]
      Z <- .btl_exposure(dd$object_a, dd$object_b, draw$resp, fit$m,
                         dd$judge, dd$order, dd$weight)
      lp <- unname(beta[dd$object_a] - beta[dd$object_b])
      if (length(dep))
        lp <- lp + drop(Z[, names(dep), drop = FALSE] %*% dep)
      oracle <- if (fit$m == 1L) stats::plogis(lp) else
        vapply(lp, function(z) item_moments(z, fit$thresholds$tau)$E, 0)
      expect_equal(draw$expected_a, oracle)
      differs_from_observed[i] <- any(abs(oracle - observed_expected) > 1e-12)
    }
    expect_true(any(differs_from_observed))
  }

  check_history(simulate_btl(
    5, 20, reps_per_pair = 10,
    dependence = list(exposure = 0.4, carry_over = 0.3), seed = 8301),
    8302)
  for (k in c(4L, 5L)) {
    check_history(simulate_btl(
      5, 20, reps_per_pair = 10, model = "polytomous",
      n_categories = k,
      dependence = list(exposure = 0.4, carry_over = 0.3),
      seed = 8310 + k), 8320 + k)
  }
})
