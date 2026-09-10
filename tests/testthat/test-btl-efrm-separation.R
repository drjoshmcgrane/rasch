.btlef_partial_separation_data <- function() {
  pair <- function(a, b, na, nb) data.frame(
    a = rep(a, na + nb), b = rep(b, na + nb),
    winner = c(rep(a, na), rep(b, nb)))
  within <- do.call(rbind, lapply(c("A", "B"), function(s)
    rbind(pair(paste0(s, 1), paste0(s, 2), 25, 75),
          pair(paste0(s, 2), paste0(s, 3), 25, 75),
          pair(paste0(s, 1), paste0(s, 3), 10, 90))))
  cross <- rbind(pair("B1", "A2", 0, 100), pair("B2", "A2", 50, 50),
                 pair("B3", "A2", 100, 0))
  d <- rbind(within, cross)
  set.seed(9102026)
  d$judge <- sample(rep(paste0("J", 1:20), length.out = nrow(d)))
  d$panel <- "P1"
  list(data = d, sets = list(A = paste0("A", 1:3), B = paste0("B", 1:3)))
}

test_that("cross-set separation allows zero-loading rows and positive units", {
  D <- cbind(alpha = c(-1, 0, 0, 1), kappa = 1)
  y <- c(0, 0, 1, 1)
  direction <- .btlef_separation_direction(D, y, 1L)
  expect_equal(direction, c(1, 0), tolerance = 1e-8)
  gain <- (2 * y - 1) * drop(D %*% direction)
  expect_equal(gain, c(1, 0, 0, 1), tolerance = 1e-8)
  # Reverse the hierarchy: increasing a positive alpha cannot fit it.
  # The finite boundary at alpha = 0 is handled by its own guard.
  expect_null(.btlef_separation_direction(D, 1 - y, 1L))
  expect_null(.btlef_separation_direction(rbind(D, D), c(y, 1 - y), 1L))
  expect_null(.btlef_separation_direction(matrix(0, 4, 2), y, 1L))

  # A separating origin can have either sign, with the unit held fixed.
  for (wins in c(0, 1)) {
    v <- .btlef_separation_direction(D, rep(wins, 4), 1L)
    expect_true(all((2 * wins - 1) * drop(D %*% v) >= -1e-8))
    expect_gt(mean((2 * wins - 1) * drop(D %*% v)), 0)
    expect_gte(v[1], 0)
  }
  # Scaling, row order, duplicating comparisons and reversing presentation
  # cannot change whether a valid recession direction exists.
  for (s in c(1e-6, 1, 1e6)) {
    Q <- sweep(D[c(4, 2, 1, 3), ], 2L, c(s, 2), "*")
    yy <- y[c(4, 2, 1, 3)]
    v <- .btlef_separation_direction(rbind(Q, -Q), c(yy, 1 - yy), 1L)
    expect_true(all((2 * yy - 1) * drop(Q %*% v) >= -1e-8))
    expect_gt(mean((2 * yy - 1) * drop(Q %*% v)), 0)
  }
})

test_that("partly separated public BTL-EFRM fits are refused at every iteration cap", {
  z <- .btlef_partial_separation_data()
  # Within-set locations are (-log(3), 0, log(3)); the cross likelihood has
  # derivative 200*log(3)*plogis(-alpha*log(3)) > 0 for every finite alpha.
  for (maxit in c(10L, 20L, 60L)) for (method in c("conditional", "judge_bootstrap"))
    expect_error(btl_efrm(z$data, "a", "b", "winner", "judge", "panel", z$sets,
      maxit = maxit, se_method = method, boot_reps = 30, workers = 1), "separated")
})

test_that("separation certificates respect mixed-outcome faces in higher dimensions", {
  set.seed(982715)
  for (p in 2:6) for (case in 1:10) {
    pos <- seq_len(ceiling(p / 2))
    v <- c(runif(length(pos), .5, 2), rnorm(floor(p / 2)))
    v <- v / sqrt(sum(v^2))
    orth <- matrix(rnorm((p - 1) * p), p - 1, p)
    orth <- orth - (orth %*% v) %*% t(v)
    positive <- matrix(rnorm(8 * p), 8, p)
    positive <- positive * drop(ifelse(positive %*% v >= 0, 1, -1))
    A <- rbind(orth, -orth, positive)
    # v is a known feasible direction. Opposite rows of orth force exact
    # zero loadings; the remaining rows have strictly positive loadings.
    out <- .btlef_separation_direction(A, rep(1, nrow(A)), pos)
    expect_false(is.null(out))
    if (is.null(out)) next
    expect_gte(min(A %*% out), -1e-7)
    expect_gt(mean(A %*% out), 1e-5)
    expect_true(all(out[pos] >= 0))
    expect_null(.btlef_separation_direction(rbind(A, A),
      c(rep(1, nrow(A)), rep(0, nrow(A))), pos))
  }
})

test_that("finite near-separated links and solver exhaustion are distinguished", {
  beta <- setNames(rep(c(-log(3), 0, log(3)), 2),
                   c(paste0("A", 1:3), paste0("B", 1:3)))
  set_of <- setNames(rep(c("A", "B"), each = 3), names(beta))
  a <- c(rep("B1", 20000), rep("B2", 100), rep("B3", 20000))
  y <- c(1, rep(0, 19999), rep(c(0, 1), 50), 0, rep(1, 19999))
  args <- list(a = a, b = rep("A2", length(a)), y = y,
    phg = rep(1, length(a)), sa = rep("B", length(a)), sb = rep("A", length(a)),
    bhat = beta, sets_u = c("A", "B"), tol = 1e-8, set_of = set_of)
  fit <- do.call(.btlef_stage2, c(args, list(maxit = 60)))
  expect_true(fit$converged)
  expect_true(fit$rank_ok)
  expect_equal(unname(fit$alpha["B"]), log(19999) / log(3), tolerance = 1e-7)
  expect_equal(unname(fit$kappa["B"]), 0, tolerance = 1e-7)
  expect_true(all(is.finite(fit$se_log_alpha)))
  expect_gt(mean(pmin(fit$p, 1 - fit$p) < 1e-4), .99)
  short <- do.call(.btlef_stage2, c(args, list(maxit = 1)))
  expect_identical(short$termination, "iteration_limit")
  expect_false(short$converged)
  expect_identical(fit$termination, "step")
})

test_that("both frame bootstraps exclude partly separated links", {
  z <- .btlef_partial_separation_data()
  d <- z$data
  d$winner[which(d$a == "B1" & d$b == "A2")[1:25]] <- "B1"
  d$winner[which(d$a == "B3" & d$b == "A2")[1:25]] <- "A2"
  sep <- z$data[substr(z$data$a, 1, 1) != substr(z$data$b, 1, 1), ]
  beta <- setNames(rep(c(-log(3), 0, log(3)), 2),
                   c(paste0("A", 1:3), paste0("B", 1:3)))
  original <- .btlef_stage2
  calls <- 0L
  testthat::local_mocked_bindings(.btlef_stage2 = function(...) {
    calls <<- calls + 1L
    if (calls == 1L) return(original(...))
    # Feed the real detector a fixed separated draw; never forge its status.
    original(sep$a, sep$b, as.integer(sep$winner == sep$a), rep(1, nrow(sep)),
      substr(sep$a, 1, 1), substr(sep$b, 1, 1), beta, c("A", "B"), 60, 1e-8)
  }, .package = "rasch")
  for (method in c("bootstrap", "judge_bootstrap")) {
    calls <- 0L
    expect_error(btl_efrm(d, "a", "b", "winner", "judge", "panel", z$sets,
      se_method = method, boot_reps = 30, workers = 1, seed = 812), "only 0 usable fits")
    expect_equal(calls, 31L)
  }
})
