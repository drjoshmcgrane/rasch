test_that("EFRM links a connected chain without direct overlap at its ends", {
  d <- simulate_efrm(n_per_group = 600, items_per_set = 5, n_sets = 3,
                     n_groups = 1, seed = 27109)
  sets <- attr(d, "truth")$item_sets
  complete <- d
  d[1:300, sets[[3]]] <- NA
  d[301:600, sets[[1]]] <- NA
  f <- rasch_efrm(d, item_sets = sets, groups = "group", id = "id",
                  boot_reps = 0, workers = 1)
  expect_true(f$est$converged)
  edges <- f$linking$alpha_edges
  expect_equal(nrow(edges), 2L)
  expect_setequal(paste(edges$set_a, edges$set_b),
                   c("set1 set2", "set2 set3"))
  expect_true(all(edges$converged %in% TRUE))
  expect_true(all(edges$n >= 30))
  expect_equal(prod(f$alpha_table$alpha), 1, tolerance = 1e-10)
  alpha <- setNames(f$alpha_table$alpha, f$alpha_table$set)
  expect_equal(unname(log(alpha[edges$set_b] / alpha[edges$set_a])),
                edges$log_slope, tolerance = 1e-10)
  expect_equal(unname(.efrm_source_matrix(f)),
                unname(as.matrix(d[, unlist(sets)])))

  # A small, redundant direct overlap does not identify an additional edge
  # and must not invalidate the two supported links.
  d[1:10, sets[[3]]] <- complete[1:10, sets[[3]]]
  sparse <- rasch_efrm(d, item_sets = sets, groups = "group", id = "id",
                       boot_reps = 0, workers = 1)
  expect_true(sparse$est$converged)
  expect_equal(nrow(sparse$linking$alpha_edges), 2L)

  # Removing the central-to-terminal overlap really disconnects the design.
  d[301:600, sets[[2]]] <- NA
  expect_error(rasch_efrm(d, item_sets = sets, groups = "group", id = "id",
                          boot_reps = 0, workers = 1),
               "not linked|not connected")
})

test_that("NPML distinguishes absent support from numerical link failure", {
  set.seed(124)
  n <- 120L
  Xm <- matrix(rbinom(n * 8, 1, .5), n, 8)
  vm <- data.frame(set = rep(c("a", "b"), each = 4), group = "g")
  pair <- function(X) .efrm_npml_pair(
    X, vm, tau_v = rep(list(0), 8), disc_v = rep(1, 8),
    sets_u = c("a", "b"), a = 1L, b = 2L, idx = seq_len(nrow(X)),
    init_log_ratio = 0, init_offset = 0, min_link_persons = 30L,
    report_support = TRUE)
  missing <- Xm
  missing[, 5:8] <- NA
  absent <- pair(missing)
  expect_s3_class(absent, "rasch_efrm_unlinked_pair")
  expect_identical(absent$reason, "too few common persons")
  short <- Xm
  short[, c(1, 5)] <- NA
  inadequate <- pair(short)
  expect_s3_class(inadequate, "rasch_efrm_unlinked_pair")
  expect_identical(inadequate$reason,
                    "too few informative common score patterns")

  # Force the optimiser to its finite search boundary on otherwise supported
  # response patterns. This must remain a failure, not an absent graph edge.
  local_mocked_bindings(
    efrm_negloglik_cpp = function(z, ...) (z[1] - log(20))^2 + z[2]^2,
    .package = "rasch")
  withr::local_options(rasch.efrm_cpp = TRUE)
  expect_null(pair(Xm))
})

test_that("supported direct set links retain their numerical estimates", {
  d <- simulate_efrm(n_per_group = 300, items_per_set = 5, n_sets = 2,
                     n_groups = 1, seed = 7321)
  sets <- attr(d, "truth")$item_sets
  current <- rasch_efrm(d, item_sets = sets, groups = "group", id = "id",
                        boot_reps = 0, workers = 1)
  original_pair <- .efrm_npml_pair
  local_mocked_bindings(
    .efrm_npml_pair = function(..., report_support = FALSE)
      original_pair(..., report_support = FALSE), .package = "rasch")
  previous <- rasch_efrm(d, item_sets = sets, groups = "group", id = "id",
                         boot_reps = 0, workers = 1)
  expect_true(current$est$converged)
  expect_true(previous$est$converged)
  expect_equal(current$alpha_table, previous$alpha_table, tolerance = 0)
  expect_equal(current$set_table, previous$set_table, tolerance = 0)
  expect_equal(current$thresholds, previous$thresholds, tolerance = 0)
  expect_equal(current$person, previous$person, tolerance = 0)
})

test_that("link graphs skip absent edges in observed and bootstrap fits", {
  n <- 80L
  u <- w <- g <- matrix(NA_real_, n, 3)
  no_edge <- structure(list(reason = "too few common persons"),
                       class = "rasch_efrm_unlinked_pair")
  calls <- 0L
  edge <- function(a, b, idx, ...) {
    calls <<- calls + 1L
    if (a == 1L && b == 3L) return(no_edge)
    # Resampling changes the link estimates, while retaining the chain.
    z <- mean(idx) / n
    list(log_ratio = z / 10, offset = -z / 5, n = n,
         converged = TRUE, edge_mass = 0, loglik = -100)
  }
  set.seed(9)
  linked <- suppressWarnings(.efrm_link_sets(
    u, w, g, c("A", "B", "C"), min_link_persons = 30L,
    pair_link = edge, boot_reps = 30, workers = 1))
  expect_equal(calls, 3L * 31L)
  expect_equal(linked$boot_reps_used, 30L)
  expect_equal(linked$boot_reps_failed, 0L)
  expect_equal(nrow(linked$edges), 2L)
  expect_true(all(is.finite(linked$cov_link)))

  # A numerical failure is not silently deleted, even if the remaining
  # graph would connect. The observed fit must keep that failure visible.
  bad_edge <- function(a, b, idx, ...) {
    if (a == 1L && b == 3L) return(NULL)
    edge(a, b, idx)
  }
  expect_error(.efrm_link_sets(
    u, w, g, c("A", "B", "C"), min_link_persons = 30L,
    pair_link = bad_edge, boot_reps = 0),
    "semiparametric likelihood link failed.*A.*C")
  disconnected <- function(a, b, idx, ...) {
    if (b == 3L) return(no_edge)
    edge(a, b, idx)
  }
  expect_error(.efrm_link_sets(
    u, w, g, c("A", "B", "C"), min_link_persons = 30L,
    pair_link = disconnected, boot_reps = 0), "not linked")
})

test_that("a redundant failed link invalidates the whole bootstrap replicate", {
  n <- 80L
  B <- 40L
  u <- w <- g <- matrix(NA_real_, n, 3)
  set.seed(9201)
  indices <- replicate(B, sample.int(n, n, replace = TRUE), simplify = FALSE)
  expected <- t(vapply(indices, function(idx) {
    z <- mean(idx) / n
    ab <- .2 + z; ac <- .8 + 2 * z; bc <- .2 - z
    c(-(ab + ac) / 3, (ab - bc) / 3, (ac + bc) / 3, 0, 0, 0)
  }, numeric(6)))
  for (failure in c("numerical", "nonconvergence")) {
    replicate_id <- -1L
    fail_ids <- seq_len(8L)
    edge <- function(a, b, idx, ...) {
      if (a == 1L && b == 2L) replicate_id <<- replicate_id + 1L
      z <- mean(idx) / n
      bad <- a == 1L && b == 3L && replicate_id %in% fail_ids
      if (bad && failure == "numerical") return(NULL)
      list(log_ratio = if (a == 1L && b == 2L) .2 + z else
             if (a == 1L && b == 3L) .8 + 2 * z else .2 - z,
           offset = 0, n = n, converged = !bad,
           edge_mass = 0, loglik = -100)
    }
    set.seed(9201)
    linked <- suppressWarnings(.efrm_link_sets(
      u, w, g, c("A", "B", "C"), min_link_persons = 30L,
      pair_link = edge, boot_reps = B, workers = 1))
    expect_equal(nrow(linked$edges), 3L)
    expect_true(all(linked$edges$converged %in% TRUE))
    expect_equal(linked$boot_reps_requested, B)
    expect_equal(linked$boot_reps_used, 32L)
    expect_equal(linked$boot_reps_failed, 8L)
    expect_equal(linked$link_reps, expected[-fail_ids, ], tolerance = 1e-12)
    expect_equal(linked$cov_link, cov(expected[-fail_ids, ]),
                 tolerance = 1e-12)

    # With 29 valid draws the absolute floor is not met, despite the two
    # surviving edges connecting all three sets in each failed draw.
    replicate_id <- -1L
    fail_ids <- seq_len(11L)
    set.seed(9201)
    expect_error(.efrm_link_sets(
      u, w, g, c("A", "B", "C"), min_link_persons = 30L,
      pair_link = edge, boot_reps = B, workers = 1),
      "29 of 40; at least 30")
  }
})

test_that("invalid supported moment links also fail the whole replicate", {
  x <- seq(-2, 2, length.out = 30L)
  u <- matrix(NA_real_, 90L, 3L)
  u[1:30, 1:2] <- x
  u[31:60, c(1, 3)] <- x
  u[61:90, 2:3] <- x
  w <- ifelse(is.finite(u), .01, NA_real_)
  g <- ifelse(is.finite(u), 1, NA_real_)
  # Only A--B loses its corrected variance. A--C--B remains connected, but
  # cannot replace the three-edge moment estimator in that replicate.
  regen <- function() {
    bad_w <- w
    bad_w[1:30, 1] <- 1e6
    list(u = u, w = bad_w, g = g)
  }
  point <- .efrm_link_sets(u, w, g, c("A", "B", "C"),
                           min_link_persons = 10L, boot_reps = 0)
  expect_equal(nrow(point$edges), 3L)
  set.seed(8194)
  expect_error(.efrm_link_sets(u, w, g, c("A", "B", "C"),
    min_link_persons = 10L, boot_reps = 30, regen = regen, workers = 1),
    "0 of 30; at least 30")
})

test_that("connected EFRM chains support hybrid and full-bootstrap uncertainty", {
  skip_on_cran()
  d <- simulate_efrm(n_per_group = 600, items_per_set = 5, n_sets = 3,
                     n_groups = 1, seed = 27109)
  sets <- attr(d, "truth")$item_sets
  d[1:300, sets[[3]]] <- NA
  d[301:600, sets[[1]]] <- NA
  for (method in c("hybrid", "bootstrap")) {
    f <- rasch_efrm(d, item_sets = sets, groups = "group", id = "id",
                    se_method = method, boot_reps = 40, workers = 1,
                    seed = 9814)
    expect_true(f$est$converged)
    expect_identical(f$se_method, method)
    expect_equal(nrow(f$linking$alpha_edges), 2L)
    expect_equal(f$boot_reps_requested, 40L)
    expect_gte(f$boot_reps_used, 30L)
    expect_equal(f$boot_reps_used + f$boot_reps_failed, 40L)
    expect_true(all(is.finite(f$alpha_table$se_log_alpha)))
    expect_true(all(f$alpha_table$se_log_alpha > 0))
  }
})
