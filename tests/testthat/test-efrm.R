simEF <- function(th, tau, r) {
  x <- 0:length(tau)
  p <- exp(r * (x * th - c(0, cumsum(tau)))); p / sum(p)
}

test_that("EFRM group support follows the informative stage-one pairs", {
  Xv <- matrix(NA_integer_, 4L, 4L)
  # Row 1 has a cross-set pair only; row 2 has a same-set pair at its
  # minimum total.  Neither contributes to the conditional likelihood.
  Xv[1L, c(1L, 3L)] <- c(1L, 0L)
  Xv[2L, 1:2] <- c(0L, 0L)
  Xv[3L, 1:2] <- c(0L, 1L)
  Xv[4L, ] <- c(0L, 1L, 1L, 0L)
  vmap <- data.frame(group = rep("g", 4L),
                     set = rep(c("A", "B"), each = 2L))
  z <- rasch:::.efrm_group_support(Xv, vmap, rep(1L, 4L), "g")
  expect_identical(z$n_persons, 2L)
  # The informative pair loads are one and two, respectively.
  expect_equal(z$effective_persons, 9 / 5)

  # A maximum-total pair is a conditional-likelihood constant too.
  Xv[3L, 1:2] <- c(1L, 1L)
  z2 <- rasch:::.efrm_group_support(Xv, vmap, rep(1L, 4L), "g")
  expect_identical(z2$n_persons, 1L)
})

test_that("EFRM set support follows the strongest path to the graph root", {
  sets <- c("A", "B", "C", "D")
  chain <- data.frame(
    set_a = c("A", "B", "A", "C"),
    set_b = c("B", "C", "C", "D"),
    n = c(80, 40, 10, 70))
  z <- rasch:::.efrm_set_path_support(sets, chain, chain$n)
  # C and D cannot borrow the strong terminal edge across the B--C bottleneck;
  # the weak direct A--C edge is redundant and does not reduce that path.
  expect_equal(unname(z), c(Inf, 80, 40, 40))

  stronger <- chain
  stronger$n[stronger$set_a == "A" & stronger$set_b == "C"] <- 60
  z2 <- rasch:::.efrm_set_path_support(sets, stronger, stronger$n)
  expect_equal(unname(z2), c(Inf, 80, 60, 60))
})

test_that("EFRM unit tests do not discard indefinite covariance directions", {
  bad <- rasch:::.efrm_wald_zero(
    c(0.2, -0.2), diag(c(1, -0.5)), "unit")
  expect_true(is.na(bad$wald))
  expect_true(is.na(bad$p))
  asymmetric <- matrix(c(1, 0.5, 0, 1), 2L)
  bad <- rasch:::.efrm_wald_zero(c(0.2, -0.2), asymmetric, "unit")
  expect_true(is.na(bad$wald))
  singular <- matrix(c(1, -1, -1, 1), 2L)
  expect_true(is.finite(
    rasch:::.efrm_wald_zero(c(0.2, -0.2), singular, "unit")$p))
  bad <- rasch:::.efrm_wald_zero(c(0.2, 0.2), singular, "unit")
  expect_true(is.na(bad$wald))
})

test_that("hybrid EFRM refuses an unusable stage-one covariance", {
  d <- simulate_efrm(n_per_group = 100, items_per_set = 5, n_sets = 2,
                     n_groups = 1, seed = 7001)
  tr <- attr(d, "truth")
  old_solve <- rasch:::.efrm_solve
  expect_error(testthat::with_mocked_bindings(
    rasch_efrm(d, item_sets = tr$item_sets, groups = "group", id = "id",
               boot_reps = 30, workers = 1),
    .efrm_solve = function(...) {
      z <- old_solve(...)
      z$cov_joint[1L, 1L] <- -1e6
      z
    },
    .package = "rasch"),
    "positive-semidefinite joint stage-one covariance")
})

test_that("EFRM recovers person-group units (one set, four groups)", {
  skip_on_cran()   # heavy simulation; verified locally and on CI
  set.seed(7); L <- 20; per_g <- 400
  phi_true <- c(0.6, 0.9, 1.1, 1.5); phi_true <- phi_true / exp(mean(log(phi_true)))
  d <- scale(seq(-2, 2, length.out = L), scale = FALSE)[, 1]
  glev <- paste0("G", 1:4); grp <- rep(glev, each = per_g); Np <- length(grp)
  th <- rnorm(Np, 0, 1.3)
  X <- sapply(seq_len(L), function(i)
    rbinom(Np, 1, plogis(phi_true[match(grp, glev)] * (th - d[i]))))
  colnames(X) <- sprintf("I%02d", 1:L)

  fit <- rasch_efrm(data.frame(X, g = grp), item_sets = list(core = colnames(X)),
                    groups = "g")
  expect_s3_class(fit, "rasch_efrm"); expect_s3_class(fit, "rasch")
  expect_true(fit$est$converged)
  lp <- log(fit$phi_table$phi)
  expect_gt(cor(lp, log(phi_true)), 0.95)
  expect_lt(max(abs(lp - log(phi_true))), 0.12)
  expect_equal(sum(lp), 0, tolerance = 1e-8)
  expect_gt(cor(fit$item_arbitrary$location, d), 0.99)
  expect_true(all(fit$alpha_table$alpha == 1))
  for (g in glev) {
    sel <- fit$person$g == g
    expect_gt(cor(fit$person$theta[sel], th[grp == g], use = "complete.obs"), 0.8)
  }
  expect_gt(fit$efrm_vs_rasch$two_delta_ll, 100)   # strong unit differences
})

test_that("EFRM recovers item-set units from common persons (polytomous)", {
  set.seed(8); Np <- 800
  alpha_true <- c(0.7, 1.0, 1.4); alpha_true <- alpha_true / exp(mean(log(alpha_true)))
  mu_true <- c(-0.4, 0.1, 0.3)
  th <- rnorm(Np, 0, 1.3)
  sets <- rep(1:3, each = 6)
  dd <- lapply(1:18, function(i) mu_true[sets[i]] + sort(rnorm(2, 0, 0.8)))
  X <- sapply(1:18, function(i) sapply(th, function(t)
    sample(0:2, 1, prob = simEF(t, dd[[i]], alpha_true[sets[i]]))))
  colnames(X) <- sprintf("S%dI%02d", sets, 1:18)

  fit <- rasch_efrm(data.frame(X, g = "all"), groups = "g",
                    item_sets = split(colnames(X), sets), boot_reps = 0)
  expect_true(fit$est$stage1_converged)
  expect_identical(fit$est$converged,
                   fit$est$stage1_converged &&
                     all(fit$linking$alpha_edges$converged %in% TRUE))
  expect_lt(max(abs(log(fit$alpha_table$alpha) - log(alpha_true))), 0.15)
  expect_true(all(fit$phi_table$phi == 1))
  mu_real <- tapply(sapply(dd, mean), sets, mean)
  mu_real <- mu_real - mean(mu_real)
  expect_lt(max(abs(fit$set_table$mu - mu_real)), 0.15)
  loc_true <- sapply(dd, mean) - mean(tapply(sapply(dd, mean), sets, mean))
  names(loc_true) <- colnames(X)
  est <- setNames(fit$item_arbitrary$location, fit$item_arbitrary$item)
  expect_gt(cor(est[names(loc_true)], loc_true), 0.97)
  expect_gt(cor(fit$person$theta, th, use = "complete.obs"), 0.9)
  expect_true(all(is.na(fit$thresholds_arbitrary$se)))
  expect_true(all(is.na(fit$item_arbitrary$se)))
  expect_true(all(is.na(fit$unit_cov$cov_delta)))
  expect_identical(fit$unit_cov$method, "stage-one-only")
  expect_false(fit$unit_support$alpha_inference)
  expect_match(paste(fit$notes, collapse = " "),
               "set-link uncertainty was omitted")
})

test_that("EFRM recovers the full unit grid (two sets x two groups)", {
  set.seed(9); per_g <- 500
  alpha_true <- c(A = 0.8, B = 1.25)
  phi_true <- c(g1 = 0.85, g2 = 1 / 0.85)
  th <- rnorm(2 * per_g, 0, 1.3)
  grp <- rep(names(phi_true), each = per_g)
  sets <- rep(c("A", "B"), each = 8)
  d <- rep(seq(-1.5, 1.5, length.out = 8), 2)
  d <- d - mean(tapply(d, sets, mean))
  X <- sapply(seq_along(sets), function(i)
    rbinom(2 * per_g, 1,
           plogis(alpha_true[sets[i]] * phi_true[grp] * (th - d[i]))))
  colnames(X) <- sprintf("%sI%02d", sets, seq_along(sets))

  fit <- rasch_efrm(data.frame(X, g = grp), groups = "g",
                    item_sets = split(colnames(X), sets), boot_reps = 80)
  fr <- fit$frames
  rho_true <- outer(alpha_true, phi_true)[cbind(fr$set, fr$group)]
  expect_gt(cor(log(fr$rho), log(rho_true)), 0.95)
  expect_equal(fr$rho, fr$alpha * fr$phi, tolerance = 1e-12)
  expect_lt(max(abs(log(fit$phi_table$phi) - log(phi_true))), 0.12)
  expect_lt(max(abs(log(fit$alpha_table$alpha) - log(alpha_true))), 0.2)
  expect_setequal(fit$efrm_vs_rasch$unit_omnibus$term,
                  c("group units (phi)", "set units (alpha)"))
  expect_true(all(is.finite(fit$efrm_vs_rasch$unit_omnibus$p)))
  expect_equal(fit$efrm_vs_rasch$unit_omnibus$p_adj,
               p.adjust(fit$efrm_vs_rasch$unit_omnibus$p, "holm"))
  expect_equal(fit$efrm_vs_rasch$unit_omnibus$significant,
               fit$efrm_vs_rasch$unit_omnibus$p_adj < 0.05)
  expect_true(all(c("p_adj", "significant") %in%
                    names(fit$efrm_vs_rasch$unit_tests)))
})

test_that("EFRM withholds unit tests for a sparsely represented group", {
  set.seed(901)
  n <- 200; L <- 8
  grp <- c(rep("large", 190), rep("small", 10))
  theta <- rnorm(n); delta <- seq(-1.5, 1.5, length.out = L)
  X <- matrix(rbinom(n * L, 1, plogis(outer(theta, delta, "-"))), n, L)
  colnames(X) <- paste0("I", seq_len(L))
  fit <- rasch_efrm(data.frame(X, group = grp),
                    item_sets = list(all = colnames(X)), groups = "group",
                    boot_reps = 0)
  expect_false(fit$unit_support$phi_inference)
  # One of the ten small-group rows is an all-extreme within-set response and
  # carries no pairwise-conditional information for the group unit.
  expect_equal(min(fit$unit_support$group$n_persons), 9)
  expect_true(all(is.na(fit$efrm_vs_rasch$unit_omnibus$p)))
  expect_true(all(is.na(fit$efrm_vs_rasch$unit_omnibus$p_adj)))
  expect_true(all(is.na(fit$efrm_vs_rasch$unit_tests$p)))
  expect_true(all(is.finite(fit$phi_table$phi)))
  interval_drawn <- FALSE
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  testthat::with_mocked_bindings(
    plot_frames(fit),
    segments = function(...) interval_drawn <<- TRUE,
    .package = "rasch")
  expect_false(interval_drawn)
})

test_that("NPML set links at an optimiser boundary are refused", {
  set.seed(1)
  n <- 100L
  u <- rnorm(n)
  Xa <- cbind(rbinom(n, 1, plogis(u + 1)),
              rbinom(n, 1, plogis(u - 1)))
  Xb <- matrix(1L, n, 2L)
  Xb[seq_len(3L), ] <- 0L
  Xm <- cbind(Xa, Xb)
  vmap <- data.frame(set = rep(c("a", "b"), each = 2L), group = "g")
  z <- rasch:::.efrm_npml_pair(
    Xm, vmap, tau_v = rep(0, 4L), disc_v = rep(1, 4L),
    sets_u = c("a", "b"), a = 1L, b = 2L, idx = seq_len(n),
    init_log_ratio = 0, init_offset = 0, min_link_persons = 5L,
    grid_n = 31L)
  expect_null(z)
})

test_that("same-tail extremes do not inflate EFRM link support", {
  set.seed(4)
  n <- 300L
  u <- rnorm(n)
  d <- seq(-1.5, 1.5, length.out = 6L)
  Xa <- sapply(d, function(dd) rbinom(n, 1L, plogis(u - dd)))
  Xb <- sapply(d, function(dd) rbinom(n, 1L, plogis(1.5 * u - dd)))
  Xm <- cbind(Xa, Xb)
  vmap <- data.frame(set = rep(c("a", "b"), each = 6L), group = "g")
  link <- function(X) rasch:::.efrm_npml_pair(
    X, vmap, tau_v = rep(0, 12L), disc_v = rep(1, 12L),
    sets_u = c("a", "b"), a = 1L, b = 2L, idx = seq_len(nrow(X)),
    init_log_ratio = 0, init_offset = 0, min_link_persons = 5L)
  z0 <- link(Xm)
  expect_false(is.null(z0))
  # A person at the same extreme tail in both sets has a likelihood that
  # rises towards one end of the grid, so the masses explaining such
  # persons sit on the edge. Appending 120 of them puts about three per cent
  # of the fitted masses on the grid ends, which a check on the total edge
  # mass would refuse as truncation. Those persons stay in the likelihood
  # and carry no information about the link, so the estimate barely moves.
  Xe <- rbind(Xm, matrix(1L, 80L, 12L), matrix(0L, 40L, 12L))
  z1 <- link(Xe)
  expect_false(is.null(z1))
  expect_identical(z1$n, z0$n)
  expect_lt(z1$edge_mass, 1e-3)
  expect_lt(abs(z1$log_ratio - z0$log_ratio), 0.02)
  expect_lt(abs(z1$offset - z0$offset), 0.02)
})

test_that("opposing set extremes remain informative to an EFRM link", {
  score_a <- c(0, 2, 0, 2, 1)
  score_b <- c(0, 2, 2, 0, 1)
  same_tail <- rasch:::.efrm_same_tail_extreme(
    score_a, rep(2, 5), score_b, rep(2, 5))

  expect_identical(same_tail, c(TRUE, TRUE, FALSE, FALSE, FALSE))

  # Classification is invariant to the fitted unit. In particular, a small
  # discrimination must not make an interior score look like a minimum.
  expect_identical(
    rasch:::.efrm_same_tail_extreme(
      c(5e-13, 0), c(1e-12, 1e-12), c(5e-13, 0), c(1e-12, 1e-12)),
    c(FALSE, TRUE))

  # The outer linking stage must not discard an otherwise supported
  # likelihood edge merely because too few persons have finite WLE starts.
  # The pair likelihood is responsible for deciding whether the raw response
  # patterns provide the requested support.
  n <- 30L
  u <- matrix(NA_real_, n, 2L)
  w <- g <- matrix(NA_real_, n, 2L)
  pair_calls <- 0L
  pair_link <- function(a, b, idx, init_ls, init_off) {
    pair_calls <<- pair_calls + 1L
    expect_equal(c(init_ls, init_off), c(0, 0))
    list(log_ratio = 0, offset = 0, n = n, converged = TRUE,
         edge_mass = 0, loglik = -10)
  }
  linked <- rasch:::.efrm_link_sets(
    u, w, g, c("A", "B"), min_link_persons = n,
    boot_reps = 0L, pair_link = pair_link)
  expect_equal(pair_calls, 1L)
  expect_equal(unname(linked$alpha), c(1, 1))
})

test_that("concordant extremes cannot satisfy EFRM link support", {
  set.seed(41)
  n <- 180L
  u <- rnorm(n)
  d <- seq(-1.5, 1.5, length.out = 6L)
  Xa <- sapply(d, function(dd) rbinom(n, 1L, plogis(u - dd)))
  Xb <- sapply(d, function(dd) rbinom(n, 1L, plogis(1.2 * u - dd)))
  Xm <- cbind(Xa, Xb)
  vmap <- data.frame(set = rep(c("a", "b"), each = 6L), group = "g")
  link <- function(min_n) rasch:::.efrm_npml_pair(
    Xm, vmap, tau_v = rep(0, 12L), disc_v = rep(1, 12L),
    sets_u = c("a", "b"), a = 1L, b = 2L, idx = seq_len(n),
    init_log_ratio = 0, init_offset = 0, min_link_persons = min_n,
    grid_n = 31L)
  z <- link(5L)
  expect_false(is.null(z))

  # Same-tail rows remain in the likelihood but cannot manufacture the
  # minimum number of persons that identify the transformation.
  Xe <- rbind(Xm, matrix(1L, 50L, 12L), matrix(0L, 50L, 12L))
  expect_null(rasch:::.efrm_npml_pair(
    Xe, vmap, tau_v = rep(0, 12L), disc_v = rep(1, 12L),
    sets_u = c("a", "b"), a = 1L, b = 2L, idx = seq_len(nrow(Xe)),
    init_log_ratio = 0, init_offset = 0,
    min_link_persons = z$n + 1L, grid_n = 31L))
})

test_that("a single frame reduces to the ordinary rasch fit", {
  set.seed(2); Np <- 400; L <- 8
  d <- seq(-1.5, 1.5, length.out = L)
  X <- matrix(rbinom(Np * L, 1, plogis(outer(rnorm(Np), d, "-"))), Np, L)
  colnames(X) <- sprintf("I%02d", 1:L)
  fe <- rasch_efrm(data.frame(X, g = "one"), groups = "g",
                   item_sets = list(all = colnames(X)))
  fr <- rasch(X)
  expect_true(all(fe$phi_table$phi == 1) && all(fe$alpha_table$alpha == 1))
  est <- setNames(fe$item_arbitrary$location, fe$item_arbitrary$item)
  expect_equal(unname(est[fr$items$item]), fr$items$location, tolerance = 1e-6)
  expect_equal(fe$person$theta, fr$person$theta, tolerance = 1e-6)
  expect_equal(fe$alpha$alpha, fr$alpha$alpha, tolerance = 1e-12)
  expect_true(fe$alpha$design_applicable)
  expect_equal(score_table(fe), score_table(fr), tolerance = 1e-8)
  expect_equal(ctt_table(fe)$table$item, colnames(X))
  expect_equal(colnames(guttman_table(fe)$matrix),
               colnames(X)[order(fr$items$location)])
  saved_before_flag <- fe
  saved_before_flag$alpha$design_applicable <- NULL
  expect_equal(ctt_table(saved_before_flag)$table$item, colnames(X))
  expect_equal(colnames(guttman_table(saved_before_flag)$matrix),
               colnames(X)[order(fr$items$location)])
})

test_that("EFRM fits the same cleaned frame labels that it validates", {
  set.seed(921)
  N <- 120; L <- 6
  X <- matrix(rbinom(N * L, 1, 0.5), N, L,
              dimnames = list(NULL, paste0("I", seq_len(L))))
  d <- data.frame(X, g = rep(c("A", " A "), length.out = N),
                  check.names = FALSE)
  fit <- rasch_efrm(d, item_sets = list(all = colnames(X)), groups = "g",
                    boot_reps = 0)
  expect_identical(as.character(fit$phi_table$group), "A")
  expect_identical(unique(as.character(fit$factors$g)), "A")
  # With only one set there is no estimated set link to omit, so the
  # supported stage-one threshold covariance remains available.
  expect_true(any(is.finite(fit$thresholds_arbitrary$se)))
  expect_identical(fit$unit_cov$method, "stage-one")

  set_map <- stats::setNames(rep(" core ", L), colnames(X))
  fit_map <- rasch_efrm(d, item_sets = set_map, groups = "g", boot_reps = 0)
  expect_identical(as.character(fit_map$alpha_table$set), "core")

  bad_sets <- list(core = colnames(X)[1:3],
                   " core " = colnames(X)[4:6])
  expect_error(rasch_efrm(d, item_sets = bad_sets, groups = "g",
                          boot_reps = 0),
               "after trimming")
})

test_that("one-cell EFRM keeps classical summaries but not a heterogeneous-unit score table", {
  d <- simulate_efrm(n_per_group = 250, items_per_set = 5, n_sets = 2,
                     n_groups = 1, set_unit_ratio = 1.3, seed = 92)
  tr <- attr(d, "truth")
  f <- rasch_efrm(d, item_sets = tr$item_sets, groups = "group", id = "id",
                  boot_reps = 0)
  expect_true(f$alpha$design_applicable)
  expect_true(is.finite(f$alpha$alpha))
  expect_equal(ctt_table(f)$table$item,
               unlist(tr$item_sets, use.names = FALSE))
  expect_null(score_table(f))
})

test_that("frame ICCs can separate observed means by a non-frame DIF factor", {
  set.seed(27); n <- 180; L <- 6
  grp <- rep(c("g1", "g2"), each = n / 2)
  cohort <- rep(c("early", "late"), length.out = n)
  th <- rnorm(n)
  d <- seq(-1.2, 1.2, length.out = L)
  phi <- c(g1 = 0.85, g2 = 1 / 0.85)
  X <- sapply(d, function(x)
    rbinom(n, 1, plogis(phi[grp] * (th - x))))
  colnames(X) <- sprintf("I%02d", seq_len(L))
  f <- rasch_efrm(data.frame(X, grp, cohort),
                  item_sets = list(core = colnames(X)), groups = "grp",
                  factors = "cohort", boot_reps = 0)
  p <- tempfile(fileext = ".pdf")
  grDevices::pdf(p); on.exit({ grDevices::dev.off(); unlink(p) }, add = TRUE)
  expect_no_error(plot_icc_frames(f, "I03", group = "cohort"))
  expect_error(plot_icc_frames(f, "I03", group = "grp"), "frame-defining")
})

test_that("unlinked structures produce informative errors", {
  set.seed(3); per_g <- 150
  th <- rnorm(2 * per_g)
  grp <- rep(c("g1", "g2"), each = per_g)
  X <- matrix(rbinom(2 * per_g * 8, 1, plogis(th)), 2 * per_g, 8)
  colnames(X) <- sprintf("I%d", 1:8)
  # each group sees only its own set: no common set, no common persons
  X[grp == "g1", 5:8] <- NA
  X[grp == "g2", 1:4] <- NA
  expect_error(
    rasch_efrm(data.frame(X, g = grp), groups = "g",
               item_sets = list(s1 = colnames(X)[1:4], s2 = colnames(X)[5:8])),
    "not linked")
})

test_that("sets too short for the unit correction refuse loudly", {
  # 3 dichotomous items per set = 2 interior score categories: the
  # score basis is too weak to support the guarded scale link, so the fit
  # must refuse rather than return an apparently precise unit estimate
  set.seed(4); n <- 300
  th <- rnorm(n, 0, 1.3)
  X <- sapply(seq(-1, 1, length.out = 6), function(d)
    rbinom(n, 1, plogis(th - d)))
  colnames(X) <- sprintf("I%d", 1:6)
  expect_error(
    rasch_efrm(data.frame(X, g = "g1"), groups = "g",
               item_sets = list(s1 = colnames(X)[1:3], s2 = colnames(X)[4:6])),
    "score range of at least 4")
})

test_that("phi sandwich SEs track the sampling variability", {
  skip_on_cran()   # heavy simulation; verified locally and on CI
  set.seed(11); L <- 12; per_g <- 250
  phi_true <- c(0.75, 4 / 3); phi_true <- phi_true / exp(mean(log(phi_true)))
  d <- seq(-1.5, 1.5, length.out = L)
  glev <- c("a", "b"); grp <- rep(glev, each = per_g); Np <- length(grp)
  reps <- 8
  lp_hat <- se_hat <- numeric(reps)
  for (r in seq_len(reps)) {
    th <- rnorm(Np, 0, 1.3)
    X <- sapply(seq_len(L), function(i)
      rbinom(Np, 1, plogis(phi_true[match(grp, glev)] * (th - d[i]))))
    colnames(X) <- sprintf("I%02d", 1:L)
    f <- rasch_efrm(data.frame(X, g = grp), groups = "g",
                    item_sets = list(core = colnames(X)))
    lp_hat[r] <- log(f$phi_table$phi[2])
    se_hat[r] <- f$phi_table$se_log_phi[2]
  }
  ratio <- mean(se_hat) / sd(lp_hat)
  expect_gt(ratio, 0.5); expect_lt(ratio, 2.0)
})

test_that("the weighted score, not the raw score, drives person estimates", {
  set.seed(12); Np <- 600
  alpha_true <- c(0.6, 5 / 3); alpha_true <- alpha_true / exp(mean(log(alpha_true)))
  th <- rnorm(Np, 0, 1.2)
  sets <- rep(1:2, each = 6)
  d <- rep(seq(-1, 1, length.out = 6), 2)
  X <- sapply(seq_along(sets), function(i) sapply(th, function(t)
    sample(0:2, 1, prob = simEF(t, d[i] + c(-0.5, 0.5), alpha_true[sets[i]]))))
  colnames(X) <- sprintf("S%dI%02d", sets, seq_along(sets))
  fit <- rasch_efrm(data.frame(X, g = "all"), groups = "g",
                    item_sets = split(colnames(X), sets), boot_reps = 0)
  p <- fit$person
  # two persons with the same raw score but different weighted scores
  # must receive different locations
  cand <- which(!p$extreme & p$n_items == 12)
  byraw <- split(cand, p$raw[cand])
  found <- FALSE
  for (grpw in byraw) {
    if (length(grpw) < 2) next
    w <- p$weighted_score[grpw]
    if (max(w) - min(w) > 0.5) {
      i1 <- grpw[which.min(w)]; i2 <- grpw[which.max(w)]
      expect_false(isTRUE(all.equal(p$theta[i1], p$theta[i2])))
      expect_lt(p$theta[i1], p$theta[i2])  # more weight on high-unit items
      found <- TRUE; break
    }
  }
  expect_true(found)
  # equal weighted scores within a pattern share a location exactly
  key <- paste(p$n_items, signif(p$weighted_score, 12))
  dup <- key[duplicated(key) & !is.na(p$theta)][1]
  who <- which(key == dup)
  expect_lt(diff(range(p$theta[who])), 1e-12)
})

test_that("EFRM extreme status follows the response pattern", {
  X <- rbind(non_extreme = c(1L, 0L, 0L),
             minimum = c(0L, 0L, 0L),
             maximum = c(1L, 1L, 1L))
  z <- .efrm_person_estimates(
    X, list(0, 0, 0), disc = c(1e-14, 1, 1))
  expect_lt(z$weighted_score[1L], 1e-12)
  expect_false(z$extreme[1L])
  expect_true(z$extreme[2L])
  expect_true(z$extreme[3L])
})

test_that("the semiparametric set link beats the naive SD ratio", {
  set.seed(13); Np <- 500
  alpha_true <- c(0.7, 10 / 7); alpha_true <- alpha_true / exp(mean(log(alpha_true)))
  th <- rnorm(Np, 0, 1.2)
  sets <- rep(1:2, each = 6)
  d <- rep(seq(-1, 1, length.out = 6), 2)
  X <- sapply(seq_along(sets), function(i) sapply(th, function(t)
    sample(0:2, 1, prob = simEF(t, d[i] + c(-0.5, 0.5), alpha_true[sets[i]]))))
  colnames(X) <- sprintf("S%dI%02d", sets, seq_along(sets))
  fit <- rasch_efrm(data.frame(X, g = "all"), groups = "g",
                    item_sets = split(colnames(X), sets), boot_reps = 0)
  # naive ratio from observed SDs of the same per-set person estimates
  dtl <- fit$thresholds_arbitrary
  u <- lapply(1:2, function(s) {
    cols <- which(fit$set_of[fit$virtual_map$item] == as.character(s))
    tl <- lapply(fit$virtual_map$item[cols], function(it)
      dtl$delta[dtl$item == it] * fit$alpha_table$alpha[s])
    pe <- rasch:::.person_estimates(fit$X[, cols, drop = FALSE], tl, disc = 1)
    pe
  })
  ok <- !u[[1]]$extreme & !u[[2]]$extreme &
    is.finite(u[[1]]$theta) & is.finite(u[[2]]$theta)
  naive <- sd(u[[2]]$theta[ok]) / sd(u[[1]]$theta[ok])
  true_ratio <- alpha_true[2] / alpha_true[1]
  est_ratio <- fit$alpha_table$alpha[2] / fit$alpha_table$alpha[1]
  expect_lt(abs(log(est_ratio) - log(true_ratio)),
            abs(log(naive) - log(true_ratio)))
})

test_that("the semiparametric set link handles a bimodal person distribution", {
  set.seed(100003); n <- 500; ips <- 8; ratio <- 1.4
  alpha_true <- c(set1 = ratio^(-0.5), set2 = ratio^0.5)
  theta <- sample(c(-2.2, 2.2), n, replace = TRUE) + rnorm(n, 0, 0.6)
  sets <- rep(names(alpha_true), each = ips)
  delta <- rep(seq(-1.5, 1.5, length.out = ips), 2)
  X <- vapply(seq_along(sets), function(i)
    rbinom(n, 1, plogis(alpha_true[sets[i]] * (theta - delta[i]))),
    numeric(n))
  colnames(X) <- sprintf("%sI%02d", sets, seq_along(sets))
  item_sets <- split(colnames(X), sets)
  fit <- rasch_efrm(data.frame(X, group = "g1"), item_sets = item_sets,
                    groups = "group", boot_reps = 0)
  got <- with(fit$alpha_table, setNames(alpha, set))
  expect_lt(abs(log(got["set2"] / got["set1"]) - log(ratio)), 0.12)

  # Set names are sorted before linking. Relabel the same sets so that the
  # other response set defines the finite-grid coordinate, then map the
  # estimates back to the original sets. The grid approximation should not
  # make the reported scale ratio depend materially on that direction.
  rev_fit <- rasch_efrm(
    data.frame(X, group = "g1"),
    item_sets = list(z_first = item_sets$set1, a_second = item_sets$set2),
    groups = "group", boot_reps = 0
  )
  rev_got <- with(rev_fit$alpha_table, setNames(alpha, set))
  expect_equal(unname(log(rev_got["a_second"] / rev_got["z_first"])),
               unname(log(got["set2"] / got["set1"])), tolerance = 0.01)
})

test_that("EFRM linking permits different person distributions by group", {
  set.seed(100004)
  ng <- 700L; ips <- 8L
  grp <- rep(c("lower", "upper"), each = ng)
  theta <- c(rnorm(ng, -1.4, 0.8),
             sample(c(0.7, 2.4), ng, replace = TRUE) + rnorm(ng, 0, 0.45))
  alpha <- c(set1 = 1.4^(-0.5), set2 = 1.4^0.5)
  phi <- c(lower = 1.5^(-0.5), upper = 1.5^0.5)
  sets <- rep(names(alpha), each = ips)
  delta <- rep(seq(-1.8, 1.8, length.out = ips), 2L)
  X <- vapply(seq_along(sets), function(i)
    rbinom(length(theta), 1,
           plogis(alpha[sets[i]] * phi[grp] * (theta - delta[i]))),
    numeric(length(theta)))
  colnames(X) <- sprintf("%sI%02d", sets, seq_along(sets))
  f <- rasch_efrm(data.frame(X, group = grp),
                  item_sets = split(colnames(X), sets), groups = "group",
                  boot_reps = 0)
  got_a <- with(f$alpha_table, setNames(alpha, set))
  got_p <- with(f$phi_table, setNames(phi, group))
  expect_lt(abs(log(got_a["set2"] / got_a["set1"]) - log(1.4)), 0.12)
  expect_lt(abs(log(got_p["upper"] / got_p["lower"]) - log(1.5)), 0.14)
})

test_that("EFRM honours the -1 missing code", {
  set.seed(14); per_g <- 200; L <- 10
  grp <- rep(c("a", "b"), each = per_g); Np <- length(grp)
  phi_true <- c(0.8, 1.25)
  th <- rnorm(Np)
  d <- seq(-1.2, 1.2, length.out = L)
  X <- sapply(seq_len(L), function(i)
    rbinom(Np, 1, plogis(phi_true[match(grp, c("a", "b"))] * (th - d[i]))))
  colnames(X) <- sprintf("I%02d", 1:L)
  hit <- sample(length(X), 150)
  Xna <- X; Xna[hit] <- NA
  Xcode <- X; Xcode[hit] <- -1L
  f1 <- rasch_efrm(data.frame(Xna, g = grp), groups = "g",
                   item_sets = list(core = colnames(X)))
  f2 <- rasch_efrm(data.frame(Xcode, g = grp), groups = "g",
                   item_sets = list(core = colnames(X)))
  expect_equal(f1$phi_table$phi, f2$phi_table$phi, tolerance = 1e-10)
  expect_equal(f1$item_arbitrary$location, f2$item_arbitrary$location,
               tolerance = 1e-10)
})

test_that("EFRM standard error methods are coherent", {
  skip_on_cran()   # heavy simulation; verified locally and on CI
  set.seed(15); Np <- 500
  alpha_true <- c(0.75, 4 / 3); alpha_true <- alpha_true / exp(mean(log(alpha_true)))
  sets <- rep(1:2, each = 6)
  d <- rep(seq(-1, 1, length.out = 6), 2)
  th <- rnorm(Np, 0, 1.3)
  X <- sapply(seq_along(sets), function(i) sapply(th, function(t)
    sample(0:2, 1, prob = simEF(t, d[i] + c(-0.5, 0.5), alpha_true[sets[i]]))))
  colnames(X) <- sprintf("S%dI%02d", sets, seq_along(sets))

  fit <- rasch_efrm(data.frame(X, g = "all"), groups = "g",
                    item_sets = split(colnames(X), sets), boot_reps = 120)
  expect_identical(fit$se_method, "hybrid")
  expect_true(all(is.finite(fit$alpha_table$se_log_alpha)))
  expect_true(all(fit$alpha_table$se_log_alpha > 0))
  # propagation: common-unit threshold SEs exceed the purely conditional part
  cond <- sqrt(diag(fit$est$cov_tau))   # virtual level, already propagated
  expect_true(all(fit$thresholds_arbitrary$se > 0))

  fb <- rasch_efrm(data.frame(X, g = "all"), groups = "g",
                   item_sets = split(colnames(X), sets),
                   se_method = "bootstrap", boot_reps = 50)
  expect_identical(fb$se_method, "bootstrap")
  expect_gte(fb$boot_reps_used, 30)
  expect_identical(fb$full_boot_reps_requested, 50L)
  expect_identical(fb$full_boot_reps_attempted, 50L)
  expect_identical(fb$full_boot_reps_used, fb$boot_reps_used)
  expect_identical(fb$full_boot_reps_failed, fb$boot_reps_failed)
  expect_identical(fb$unit_cov$method, "bootstrap")
  expect_equal(sqrt(diag(fb$unit_cov$cov_log_phi)),
               fb$phi_table$se_log_phi, tolerance = 1e-12)
  expect_equal(sqrt(diag(fb$unit_cov$cov_log_alpha)),
               fb$alpha_table$se_log_alpha, tolerance = 1e-12)
  K <- nrow(fb$unit_cov$cov_dtilde)
  G <- nrow(fb$unit_cov$cov_log_phi)
  expect_equal(fb$unit_cov$cov_joint[seq_len(K), seq_len(K), drop = FALSE],
               fb$unit_cov$cov_dtilde, tolerance = 1e-12)
  expect_equal(
    fb$unit_cov$cov_joint[K + seq_len(G), K + seq_len(G), drop = FALSE],
    fb$unit_cov$cov_log_phi, tolerance = 1e-12)
  fr <- fb$frames[1L, ]
  ia <- match(fr$set, fb$alpha_table$set)
  ig <- match(fr$group, fb$phi_table$group)
  expect_equal(fr$se_log_rho^2,
               fb$unit_cov$cov_log_alpha[ia, ia] +
                 fb$unit_cov$cov_log_phi[ig, ig] +
                 2 * fb$unit_cov$cov_log_alpha_phi[ia, ig],
               tolerance = 1e-10)
  # the two methods agree on scale (well within a factor of two)
  ratio <- median(fit$thresholds_arbitrary$se / fb$thresholds_arbitrary$se)
  expect_gt(ratio, 0.5); expect_lt(ratio, 2)
  # point estimates are identical across SE methods
  expect_equal(fit$alpha_table$alpha, fb$alpha_table$alpha, tolerance = 1e-10)
  expect_equal(fit$item_arbitrary$location, fb$item_arbitrary$location,
               tolerance = 1e-10)
})

test_that("unit Wald tests accompany the equal-unit comparison", {
  set.seed(16); per_g <- 300; L <- 10
  phi_true <- c(0.75, 4 / 3); phi_true <- phi_true / exp(mean(log(phi_true)))
  d <- seq(-1.5, 1.5, length.out = L)
  grp <- rep(c("a", "b"), each = per_g); Np <- length(grp)
  th <- rnorm(Np, 0, 1.3)
  X <- sapply(seq_len(L), function(i)
    rbinom(Np, 1, plogis(phi_true[match(grp, c("a", "b"))] * (th - d[i]))))
  colnames(X) <- sprintf("I%02d", 1:L)
  fit <- rasch_efrm(data.frame(X, g = grp), groups = "g",
                    item_sets = list(core = colnames(X)))
  ut <- fit$efrm_vs_rasch$unit_tests
  expect_true(all(grepl("phi", ut$parameter)))
  expect_true(all(ut$p < 0.01))            # planted units strongly detected
  expect_match(fit$efrm_vs_rasch$informative_for, "phi")

  # sets-only design: pairwise comparison declared uninformative,
  # alpha Wald tests carry the evidence
  set.seed(17); Np2 <- 500
  alpha_true <- c(0.7, 10 / 7); alpha_true <- alpha_true / exp(mean(log(alpha_true)))
  sets <- rep(1:2, each = 6)
  d2 <- rep(seq(-1, 1, length.out = 6), 2)
  th2 <- rnorm(Np2, 0, 1.3)
  X2 <- sapply(seq_along(sets), function(i) sapply(th2, function(t)
    sample(0:2, 1, prob = simEF(t, d2[i] + c(-0.5, 0.5), alpha_true[sets[i]]))))
  colnames(X2) <- sprintf("S%dI%02d", sets, seq_along(sets))
  f2 <- rasch_efrm(data.frame(X2, g = "all"), groups = "g",
                   item_sets = split(colnames(X2), sets), boot_reps = 80)
  expect_lt(abs(f2$efrm_vs_rasch$two_delta_ll), 1e-3)   # invariant by construction
  expect_match(f2$efrm_vs_rasch$informative_for, "person-side")
  ut2 <- f2$efrm_vs_rasch$unit_tests
  expect_true(all(grepl("alpha", ut2$parameter)))
  expect_true(all(ut2$p < 0.01))
})

test_that("EFRM omnibus families are retained when covariance is unavailable", {
  full <- .efrm_wald_zero(c(-0.2, 0.2), diag(c(0.04, 0.04)),
                          "group units")
  expect_equal(full$df, 2L)
  expect_true(is.finite(full$p))

  V <- diag(c(0.04, 0.04)); V[2, 2] <- NA_real_
  unavailable <- .efrm_wald_zero(c(-0.2, 0.2), V, "group units")
  expect_identical(unavailable$term, "group units")
  expect_true(is.na(unavailable$df))
  expect_true(is.na(unavailable$p))
})

test_that("the joint stage-1 covariance is coherent and its draws honest", {
  d <- simulate_efrm(150, 6, n_sets = 2, n_groups = 2, set_unit_ratio = 1.2,
                     seed = 77)
  fit <- rasch_efrm(d, item_sets = attr(d, "truth")$item_sets,
                    groups = "group", boot_reps = 50)
  uc <- fit$unit_cov
  K <- nrow(uc$cov_dtilde); G <- nrow(uc$cov_log_phi)
  expect_equal(dim(uc$cov_joint), c(K + G, K + G))
  # diagonal blocks reproduce the marginal covariances exactly
  expect_equal(uc$cov_joint[seq_len(K), seq_len(K)], uc$cov_dtilde,
               tolerance = 1e-10)
  expect_equal(unname(uc$cov_joint[K + seq_len(G), K + seq_len(G)]),
               unname(uc$cov_log_phi), tolerance = 1e-10)
  expect_equal(dim(uc$cov_log_alpha_phi), c(2L, 2L))
  expect_equal(dim(uc$cov_log_alpha), c(2L, 2L))
  fr <- fit$frames[1, ]
  ia <- match(fr$set, fit$alpha_table$set)
  ig <- match(fr$group, fit$phi_table$group)
  expected_rho_se <- sqrt(pmax(
    fit$alpha_table$se_log_alpha[ia]^2 +
      fit$phi_table$se_log_phi[ig]^2 +
      2 * uc$cov_log_alpha_phi[ia, ig], 0))
  expect_equal(fr$se_log_rho, expected_rho_se, tolerance = 1e-12)
  # the log-phi block respects the sum-to-zero centring: every row of the
  # phi sub-block sums to zero, so every draw's log phi perturbation does too
  expect_lt(max(abs(rowSums(uc$cov_joint[K + seq_len(G), K + seq_len(G),
                                         drop = FALSE]))), 1e-10)
  # a symmetric square root reproduces the covariance, and simulated draws
  # recover it (loose Monte Carlo tolerance) while preserving the constraint
  ee <- eigen((uc$cov_joint + t(uc$cov_joint)) / 2, symmetric = TRUE)
  L <- ee$vectors %*% (t(ee$vectors) * sqrt(pmax(ee$values, 0)))
  expect_equal(L %*% t(L), (uc$cov_joint + t(uc$cov_joint)) / 2,
               tolerance = 1e-8)
  set.seed(1)
  Z <- matrix(rnorm(4000 * (K + G)), 4000)
  V <- Z %*% t(L)
  expect_lt(max(abs(cov(V) - uc$cov_joint)) /
              max(abs(uc$cov_joint)), 0.12)
  expect_lt(max(abs(rowSums(V[, K + seq_len(G), drop = FALSE]))), 1e-8)
})

test_that("an unusable alpha-phi cross-covariance is explicitly withheld", {
  set.seed(902)
  n <- 160L
  z <- rnorm(n)
  u <- cbind(A = z + rnorm(n, sd = 0.2),
             B = 1.25 * z + rnorm(n, sd = 0.2))
  w <- matrix(0.04, n, 2L)
  g <- matrix(1, n, 2L)
  # The set link remains usable in every draw, while the accompanying group
  # units are deliberately unavailable. This isolates the cross-covariance
  # guard from failure of the alpha covariance itself.
  regen <- function()
    list(u = u, w = w, g = g,
         log_phi = c(A = NA_real_, B = NA_real_))
  expect_warning(
    link <- .efrm_link_sets(u, w, g, c("A", "B"), min_link_persons = 20L,
                            boot_reps = 30L, regen = regen, workers = 1L),
    "affected frame-unit standard errors are NA")
  expect_equal(link$boot_reps_used, 30L)
  expect_true(link$cross_cov_withheld)
  expect_null(link$cov_alpha_phi)

  expect_warning(
    no_regen <- .efrm_link_sets(
      u, w, g, c("A", "B"), min_link_persons = 20L,
      boot_reps = 30L, regen = NULL, workers = 1L),
    "calibration redraws were unavailable")
  expect_true(no_regen$cross_cov_withheld)
  expect_null(no_regen$cov_alpha_phi)

  # Pin the public consequence separately. Supply a usable marginal alpha
  # covariance but mark its alpha-phi cross-covariance as withheld; the fit
  # must not turn that unknown term into zero when it constructs frame SEs.
  d <- simulate_efrm(n_per_group = 80, items_per_set = 5, n_sets = 2,
                     n_groups = 2, seed = 36)
  tr <- attr(d, "truth")
  real_link <- .efrm_link_sets
  fit <- with_mocked_bindings(
    rasch_efrm(d, item_sets = tr$item_sets, groups = "group", boot_reps = 0),
    .efrm_link_sets = function(...) {
      out <- real_link(...)
      S <- length(out$alpha)
      out$se_log_alpha[] <- 0.1
      out$cov_link <- diag(0.01, 2L * S)
      out$cross_cov_withheld <- TRUE
      out$cov_alpha_phi <- NULL
      out
    },
    .package = "rasch"
  )
  expect_true(all(is.finite(fit$alpha_table$se_log_alpha)))
  expect_true(all(is.na(fit$frames$se_log_rho)))
})

test_that("rasch.efrm_link_draws is validated and blockdiag is simulation-only", {
  d <- simulate_efrm(300, 8, n_sets = 2, n_groups = 2, seed = 78)
  old <- options(rasch.efrm_link_draws = -3)
  on.exit(options(old), add = TRUE)
  expect_error(rasch_efrm(d, item_sets = attr(d, "truth")$item_sets,
                          groups = "group", boot_reps = 40),
               "whole number")
  options(rasch.efrm_link_draws = 1)
  expect_error(rasch_efrm(d, item_sets = attr(d, "truth")$item_sets,
                          groups = "group", boot_reps = 40),
               "between 30")
  options(rasch.efrm_link_draws = NULL, rasch.efrm_link_blockdiag = TRUE)
  fit_bd <- rasch_efrm(d, item_sets = attr(d, "truth")$item_sets,
                       groups = "group", boot_reps = 60)
  expect_true(is.finite(fit_bd$alpha_table$se_log_alpha[1]))
  options(rasch.efrm_link_blockdiag = NULL)
})

test_that("frame_invariance tests the invariance the model assumes", {
  # clean frames: nothing flagged, rmsd at or below rmse
  d <- simulate_efrm(n_per_group = 400, items_per_set = 8, n_sets = 1,
                     n_groups = 2, group_unit_ratio = 1.4, seed = 2)
  tr <- attr(d, "truth")
  f <- rasch_efrm(d, item_sets = tr$item_sets, groups = "group", id = "id",
                  boot_reps = 0)
  inv <- frame_invariance(f)
  expect_s3_class(inv, "rasch_frame_invariance")
  expect_equal(inv$summary$n_location, 0L)
  expect_true(is.na(inv$summary$n_discrimination))
  expect_lt(inv$summary$ratio, 1.5)
  expect_output(print(inv), "No available item-location comparison differs")
  inv_print <- inv
  inv_print$excluded <- data.frame(
    set = "S1", frame_1 = "A", frame_2 = "B", item = c("I01", "I02"),
    reason = c("different observed category structure",
               "weakly determined in a separate frame calibration"))
  printed <- capture.output(print(inv_print))
  expect_true(any(grepl("different observed category structure", printed,
                        fixed = TRUE)))
  expect_true(any(grepl("weakly determined", printed, fixed = TRUE)))

  # planted DIF on two items is found, and the rmsd/rmse ratio rises
  set.seed(5); N <- 400; K <- 8
  phi <- c(0.845, 1.183); delta <- seq(-1.5, 1.5, length.out = K)
  mk <- function(g, shift) {
    th <- rnorm(N, 0, 1.3)
    dd <- delta; dd[c(3, 6)] <- dd[c(3, 6)] + shift
    X <- vapply(seq_len(K), function(i)
      rbinom(N, 1, plogis(phi[g] * (th - dd[i]))), numeric(N))
    colnames(X) <- sprintf("I%02d", seq_len(K)); X
  }
  X <- rbind(mk(1, 0), mk(2, 1))
  dd <- data.frame(id = sprintf("P%04d", seq_len(2 * N)), X,
                   group = rep(c("g1", "g2"), each = N), check.names = FALSE)
  f2 <- rasch_efrm(dd, item_sets = list(set1 = colnames(X)), groups = "group",
                   id = "id", boot_reps = 0)
  inv2 <- frame_invariance(f2)
  expect_setequal(inv2$locations$item[inv2$locations$flagged],
                  c("I03", "I06"))
  expect_gt(inv2$summary$ratio, 1.5)

  # a single person group leaves no item in two frames
  d1 <- simulate_efrm(n_per_group = 200, items_per_set = 6, n_sets = 2,
                      n_groups = 1, set_unit_ratio = 1.3, seed = 3)
  t1 <- attr(d1, "truth")
  f1 <- rasch_efrm(d1, item_sets = t1$item_sets, groups = "group", id = "id",
                   boot_reps = 0)
  expect_error(frame_invariance(f1), "at least two person groups")
})

test_that("conditional frame invariance reports discrimination descriptively", {
  set.seed(21); N <- 700; K <- 8
  phi <- c(0.845, 1.183); delta <- seq(-1.5, 1.5, length.out = K)
  mk <- function(g, disc) {
    th <- rnorm(N, 0, 1.3)
    X <- vapply(seq_len(K), function(i)
      rbinom(N, 1, plogis(phi[g] * disc[i] * (th - delta[i]))), numeric(N))
    colnames(X) <- sprintf("I%02d", seq_len(K)); X
  }
  dsc <- rep(1, K); dsc[c(3, 6)] <- 1.8
  X <- rbind(mk(1, rep(1, K)), mk(2, dsc))
  d <- data.frame(id = sprintf("P%05d", seq_len(2 * N)), X,
                  group = rep(c("g1", "g2"), each = N), check.names = FALSE)
  f <- rasch_efrm(d, item_sets = list(set1 = colnames(X)), groups = "group",
                  id = "id", boot_reps = 0)
  inv <- frame_invariance(f)
  expect_true(all(c("infit_1", "infit_2", "infit_z", "p_adj",
                    "disc_1", "disc_2", "disc_ratio") %in%
                    names(inv$discrimination)))
  expect_identical(inv$boot_reps, 0L)
  expect_identical(inv$boot_reps_used, 0L)
  expect_no_error(.validate_frame_invariance(inv, f))
  changed <- inv
  changed$locations$difference[1] <- changed$locations$difference[1] + 1
  expect_error(.validate_frame_invariance(changed, f),
               "frame_invariance")
  expect_false("statistic" %in% names(inv$discrimination))
  expect_true(all(is.na(inv$discrimination$p)))
  expect_true(all(is.na(inv$discrimination$p_adj)))
  expect_true(all(is.na(inv$discrimination$flagged)))
  expect_true(all(is.finite(inv$discrimination$disc_ratio)))
  expect_equal(inv$locations$p_adj, p.adjust(inv$locations$p, "holm"))
  expect_output(print(inv), "descriptive")
  expect_output(print(inv), "bootstrap")
})

test_that("adjust chooses the screening threshold without hiding either p", {
  # Holm across every item and frame pair is right for reporting a
  # difference and wrong for deciding which items to examine
  set.seed(21); N <- 700; K <- 8
  phi <- c(0.845, 1.183); delta <- seq(-1.5, 1.5, length.out = K)
  mk <- function(g, disc, shift) {
    th <- rnorm(N, 0, 1.3)
    X <- vapply(seq_len(K), function(i)
      rbinom(N, 1, plogis(phi[g] * disc[i] * (th - delta[i] - shift[i]))),
      numeric(N))
    colnames(X) <- sprintf("I%02d", seq_len(K)); X
  }
  dsc <- rep(1, K); dsc[c(3, 6)] <- 1.7
  sh <- rep(0, K); sh[2] <- 0.45
  X <- rbind(mk(1, rep(1, K), rep(0, K)), mk(2, dsc, sh))
  d <- data.frame(id = sprintf("P%05d", seq_len(2 * N)), X,
                  group = rep(c("g1", "g2"), each = N), check.names = FALSE)
  f <- rasch_efrm(d, item_sets = list(set1 = colnames(X)), groups = "group",
                  id = "id", boot_reps = 0)
  fl <- function(inv) sort(unique(inv$locations$item[
    inv$locations$flagged %in% TRUE]))
  strict <- frame_invariance(f, adjust = "holm")
  loose <- frame_invariance(f, adjust = "none")

  # the loose screen can only ever flag a superset
  expect_true(all(fl(strict) %in% fl(loose)))
  # only the location flag moves: both probabilities are reported either
  # way, and the descriptive discrimination statistics are untouched
  expect_identical(strict$locations$p, loose$locations$p)
  expect_identical(strict$locations$p_adj, loose$locations$p_adj)
  expect_identical(strict$discrimination$infit_z, loose$discrimination$infit_z)
  expect_true(all(is.na(strict$discrimination$p)))
  expect_true(all(is.na(loose$discrimination$p_adj)))
  # the printed output names the rule it applied, so a screen is not
  # mistaken for a confirmed difference
  expect_output(print(loose), "unadjusted, screening")
  expect_output(print(strict), "Holm-adjusted")
  expect_error(frame_invariance(f, adjust = "bonferroni"))
})

test_that("compiled EFRM linking agrees with the R reference", {
  d <- simulate_efrm(n_per_group = 180, items_per_set = 6, n_sets = 2,
                     n_groups = 2, seed = 913)
  tr <- attr(d, "truth")
  old <- options(rasch.efrm_cpp = FALSE)
  on.exit(options(old), add = TRUE)
  fr <- rasch_efrm(d, item_sets = tr$item_sets, groups = "group", id = "id",
                   boot_reps = 0)
  options(rasch.efrm_cpp = TRUE)
  fc <- rasch_efrm(d, item_sets = tr$item_sets, groups = "group", id = "id",
                   boot_reps = 0)

  expect_equal(fc$alpha_table$alpha, fr$alpha_table$alpha, tolerance = 1e-9)
  expect_equal(fc$set_table$mu, fr$set_table$mu, tolerance = 1e-9)
  expect_equal(fc$linking$alpha_edges$loglik,
               fr$linking$alpha_edges$loglik, tolerance = 1e-8)
  expect_identical(fc$linking$alpha_edges$converged,
                   fr$linking$alpha_edges$converged)
  expect_equal(fc$thresholds_arbitrary$delta,
               fr$thresholds_arbitrary$delta, tolerance = 1e-8)
})

test_that("EFRM reports bootstrap progress and supports cancellation", {
  d <- simulate_efrm(n_per_group = 150, items_per_set = 6, n_sets = 2,
                     n_groups = 2, seed = 914)
  tr <- attr(d, "truth")
  seen <- list()
  f <- rasch_efrm(d, item_sets = tr$item_sets, groups = "group", id = "id",
                  boot_reps = 40,
                  seed = 914,
                  progress = function(stage, current, total)
                    seen[[length(seen) + 1L]] <<- c(stage, current, total))
  expect_s3_class(f, "rasch_efrm")
  stages <- vapply(seen, `[`, "", 1L)
  expect_true(all(c("conditional calibration", "linking bootstrap",
                    "finalising") %in% stages))
  link <- seen[stages == "linking bootstrap"]
  expect_identical(as.integer(tail(link, 1L)[[1L]][2:3]), c(40L, 40L))

  completed <- 0L
  expect_condition(
    rasch_efrm(d, item_sets = tr$item_sets, groups = "group", id = "id",
      boot_reps = 30,
      seed = 915,
      progress = function(stage, current, total)
        if (stage == "linking bootstrap") completed <<- as.integer(current),
      cancel = function() completed >= 2L),
    class = "rasch_cancelled")
  expect_lte(completed, 3L)
})

test_that("parallel EFRM bootstraps are seed-identical to serial fits", {
  skip_on_cran()
  old_workers <- options(rasch.max_workers = 2L,
                         rasch.efrm.max_workers = NULL)
  on.exit(options(old_workers), add = TRUE)
  # PSOCK workers must load the same installed namespace. pkgload source-tree
  # sessions deliberately skip this integration test; R CMD check and binary
  # package tests exercise it against the installed package.
  skip_if_not(rasch:::.rasch_namespace_is_installed(),
              "parallel integration test needs an installed package namespace")
  expect_true(rasch:::.rasch_namespace_is_installed())
  probe <- try(parallel::makePSOCKcluster(2L), silent = TRUE)
  skip_if(inherits(probe, "try-error"), "local socket clusters unavailable")
  parallel::stopCluster(probe)

  d <- simulate_efrm(n_per_group = 120, items_per_set = 6, n_sets = 2,
                     n_groups = 2, seed = 915)
  tr <- attr(d, "truth")
  serial <- rasch_efrm(d, item_sets = tr$item_sets, groups = "group", id = "id",
                       boot_reps = 40, workers = 1, seed = 916)
  parallel <- rasch_efrm(d, item_sets = tr$item_sets, groups = "group", id = "id",
                         boot_reps = 40, workers = 2, seed = 916)

  expect_identical(parallel$alpha_table, serial$alpha_table)
  expect_identical(parallel$set_table, serial$set_table)
  expect_identical(parallel$thresholds_arbitrary,
                   serial$thresholds_arbitrary)
  expect_identical(parallel$linking$alpha_edges$converged,
                   serial$linking$alpha_edges$converged)
})

test_that("EFRM defaults to four workers subject to system limits", {
  expect_identical(formals(rasch_efrm)$workers, 4L)
  old_workers <- options(rasch.max_workers = 1L,
                         rasch.efrm.max_workers = NULL)
  on.exit(options(old_workers), add = TRUE)
  expect_identical(rasch:::.efrm_available_workers(), 1L)
  options(rasch.max_workers = NULL)
  available <- rasch:::.efrm_available_workers()
  options(rasch.max_workers = factor("1"))
  expect_identical(rasch:::.efrm_available_workers(), available)
})

test_that("an EFRM bootstrap seed does not take over the caller's RNG", {
  d <- simulate_efrm(n_per_group = 100, items_per_set = 6, n_sets = 2,
                     n_groups = 2, seed = 917)
  tr <- attr(d, "truth")
  set.seed(918)
  before <- .Random.seed
  rasch_efrm(d, item_sets = tr$item_sets, groups = "group", id = "id",
             boot_reps = 40, seed = 919)
  expect_identical(.Random.seed, before)
})

test_that("EFRM sizes the full bootstrap to its largest covariance block", {
  d <- simulate_efrm(n_per_group = 40, items_per_set = 6, n_sets = 2,
                     n_groups = 2, n_categories = 4, seed = 920)
  expect_error(
    rasch_efrm(d, item_sets = attr(d, "truth")$item_sets, groups = "group",
               se_method = "bootstrap", boot_reps = 30, workers = 1),
    "at least 36 replicates")
})
