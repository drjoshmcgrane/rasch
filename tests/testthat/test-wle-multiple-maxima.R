test_that("WLE selects the greatest weighted likelihood in separated banks", {
  d <- c(rep(-6, 3), rep(6, 4))
  objective <- function(th, score = 3) {
    p <- plogis(th - d)
    score * th - sum(log1p(exp(th - d))) + log(sum(p * (1 - p))) / 2
  }
  best <- optimize(objective, c(2, 6), maximum = TRUE, tol = 1e-10)$maximum
  got <- person_wle(as.list(d))
  expect_equal(unname(got$theta["3"]), best, tolerance = 1e-7)
  expect_gt(objective(got$theta["3"]), objective(-4.053434) + .025)
  for (unit in c(1e-100, .5, 2, 1e100)) {
    other <- person_wle(as.list(d / unit + 20 / unit), disc = unit)
    expect_equal(other$theta * unit - 20, got$theta, tolerance = 1e-7)
    expect_equal(other$se * unit, got$se, tolerance = 1e-7)
  }
  X <- matrix(c(1, 1, 1, 0, 0, 0, 0), nrow = 1L)
  frame <- .efrm_person_estimates(X, as.list(d), rep(1, 7))
  expect_equal(frame$theta, unname(got$theta["3"]), tolerance = 1e-8)
  curve <- .person_wle_curve(as.list(d), rep(1, 7), rep(2, 7))
  expect_equal(.person_wle_maximum(curve, 3), frame$theta)
})

test_that("WLE handles equal maxima and disordered polytomous thresholds", {
  d <- c(rep(-6, 5), rep(6, 5))
  theta <- person_wle(as.list(d))$theta["5"]
  expect_lt(theta, -3)
  expect_equal(unname(theta), -3.601, tolerance = .001)
  tau <- list(c(12, -12), c(-8, -7), c(7, 8))
  objective <- function(th) {
    lp <- lapply(tau, function(tt) (0:length(tt)) * th - c(0, cumsum(tt)))
    psi <- sum(vapply(lp, function(z) max(z) + log(sum(exp(z - max(z)))), 0))
    info <- sum(vapply(tau, function(tt) item_moments(th, tt)$V, 0))
    3 * th - psi + log(info) / 2
  }
  z <- person_wle(tau)$theta["3"]
  grid <- seq(-20, 20, by = .02)
  expect_gte(objective(z) + 1e-9, max(vapply(grid, objective, 0)))
})

test_that("unequal external weights use the integral of their own score", {
  d <- c(rep(-6, 3), rep(6, 4))
  q <- c(rep(1, 6), 1.1); q <- q / max(q)
  response <- c(1, 1, 1, 0, 0, 0, 0)
  score <- function(th) {
    p <- plogis(th - d); V <- p * (1 - p)
    H <- sum(q * V); J <- sum(q^2 * V)
    sum(q * (response - p)) + J * sum(q * V * (1 - 2 * p)) / (2 * H^2)
  }
  grid <- seq(-15, 15, by = .02)
  g <- vapply(grid, score, 0)
  crossing <- which(g[-length(g)] > 0 & g[-1] < 0)
  roots <- vapply(crossing, function(i)
    uniroot(score, grid[c(i, i + 1L)], tol = 1e-10)$root, 0)
  values <- vapply(roots, function(r)
    integrate(Vectorize(score), roots[1L], r)$value, 0)
  expect_length(roots, 2L)
  curve <- .person_wle_curve(as.list(d), rep(1, 7), q)
  got <- .person_wle_maximum(curve, sum(q * response))
  expect_equal(got, roots[which.max(values)], tolerance = 1e-8)
})
