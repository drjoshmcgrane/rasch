# Numerical validation of the conditional estimation machinery: the
# pairwise conditional score and information against finite differences,
# and the extended-frames stage-1 solution as a constrained maximum of the
# within-frame pairwise conditional likelihood.

test_that("the pairwise conditional gradient and Hessian match finite differences", {
  set.seed(8)
  tau <- list(c(-1.2, 0.3), c(-0.4, 0.9), c(-0.8), c(0.5), c(-0.2, 0.6))
  X <- sapply(tau, function(tt) vapply(rnorm(180), function(b)
    sample(0:length(tt), 1, prob = item_moments(b, tt)$P), 0L))
  colnames(X) <- paste0("I", 1:5)
  m <- vapply(tau, length, 1L)
  thr <- threshold_index(m)
  pairs <- .pair_counts(X, m)
  t0 <- unlist(tau) + rnorm(nrow(thr), 0, 0.15)   # a non-stationary point

  f <- function(tv) .pcml_glh(tv, thr, pairs, m)$ll
  glh <- .pcml_glh(t0, thr, pairs, m)
  eps <- 1e-6
  g_num <- vapply(seq_along(t0), function(k) {
    e <- rep(0, length(t0)); e[k] <- eps
    (f(t0 + e) - f(t0 - e)) / (2 * eps)
  }, 0)
  expect_equal(unname(glh$g), g_num, tolerance = 1e-5)

  # Hessian columns against finite differences of the gradient
  for (k in c(1, 4, 8)) {
    e <- rep(0, length(t0)); e[k] <- eps
    h_num <- (.pcml_glh(t0 + e, thr, pairs, m)$g -
              .pcml_glh(t0 - e, thr, pairs, m)$g) / (2 * eps)
    expect_equal(unname(glh$H[, k]), unname(h_num), tolerance = 1e-4)
  }
})

test_that("the EFRM stage-1 estimate maximises the within-frame conditional likelihood", {
  set.seed(12)
  # 2 groups x 2 sets, both groups take both sets: units identified
  N <- 260
  grp <- rep(c("g1", "g2"), each = N / 2)
  phi <- c(g1 = 1, g2 = 1.5)
  tauA <- list(c(-1, 0.6), c(-0.3, 0.8), c(-0.9, 0.2))
  tauB <- list(c(-0.7, 0.5), c(0, 0.9), c(-1.1, 0.4))
  th <- rnorm(N)
  gen <- function(tt, rho) vapply(seq_len(N), function(n) {
    p <- item_moments(th[n], tt, disc = rho[n])$P
    sample(0:length(tt), 1, prob = p)
  }, 0L)
  rho <- phi[grp]
  X <- cbind(sapply(tauA, gen, rho = rho), sapply(tauB, gen, rho = rho))
  colnames(X) <- c(paste0("A", 1:3), paste0("B", 1:3))
  fit <- rasch_efrm(X, item_sets = list(A = paste0("A", 1:3),
                                        B = paste0("B", 1:3)),
                    groups = grp, se_method = "hybrid")

  # rebuild the stage-1 objective: the same-set pairwise conditional
  # log-likelihood over the virtual items, as a function of (dtilde, phi)
  vmap <- fit$virtual_map
  Xv <- fit$X
  m_v <- fit$m
  thr_v <- threshold_index(m_v)
  pairs <- .efrm_filter_pairs(.pair_counts(Xv, m_v), vmap)
  # cross-set exclusion is structural: every retained pair shares a set
  for (pc in pairs) expect_identical(vmap$set[pc$i], vmap$set[pc$j])

  # the fitted virtual thresholds as the likelihood sees them: the fit
  # reports common-unit thresholds with the frame units carried in disc,
  # so the effective thresholds are disc * tau
  tau_hat <- fit$disc[thr_v$item] * unlist(fit$tau_list)
  ll_hat <- .pcml_glh(tau_hat, thr_v, pairs, m_v)$ll
  expect_true(is.finite(ll_hat))

  # perturb along random FEASIBLE directions of the model structure:
  # tau_v = phi_g(v) x dtilde, where dtilde is COMMON to the group copies
  # of each set threshold, sum-zero within each set, and prod(phi) = 1.
  # No such perturbation may increase the conditional likelihood.
  phi_hat <- setNames(fit$phi_table$phi, fit$phi_table$group)
  gv <- vmap$group[thr_v$item]
  sv <- vmap$set[thr_v$item]
  key <- paste(vmap$item[thr_v$item], thr_v$k)  # the shared dtilde parameter
  # the reported thresholds live on the common arbitrary-unit scale, which
  # adds a per-frame origin constant; the pairwise likelihood is invariant
  # to exactly those constants, so centring within each frame recovers the
  # stage-1 dtilde
  frame <- paste(sv, gv)
  dtilde_hat <- tau_hat / phi_hat[gv]
  dtilde_hat <- dtilde_hat - ave(dtilde_hat, frame)
  # confirm the group copies agree exactly (the structural map is bilinear)
  expect_true(all(abs(tapply(dtilde_hat, key, function(v)
    diff(range(v)))) < 1e-6))
  # per-frame origin constants cancel in the pairwise likelihood, so the
  # centred reconstruction reproduces the fitted value
  ll_centred <- .pcml_glh(phi_hat[gv] * dtilde_hat, thr_v, pairs, m_v)$ll
  expect_equal(ll_centred, ll_hat, tolerance = 1e-8)
  ll_hat <- ll_centred
  ukey <- unique(key)
  set_of_key <- sv[match(ukey, key)]
  set.seed(99)
  worse <- 0L
  for (rep in 1:12) {
    du <- rnorm(length(ukey))
    for (s in unique(set_of_key))                 # sum-zero per set
      du[set_of_key == s] <- du[set_of_key == s] - mean(du[set_of_key == s])
    d_dir <- du[match(key, ukey)]                 # broadcast to copies
    lphi_dir <- rnorm(length(phi_hat)); lphi_dir <- lphi_dir - mean(lphi_dir)
    eps <- 5e-3
    phi_p <- phi_hat * exp(eps * lphi_dir[match(gv, names(phi_hat))])
    tau_p <- phi_p * (dtilde_hat + eps * d_dir)
    ll_p <- .pcml_glh(tau_p, thr_v, pairs, m_v)$ll
    if (ll_p <= ll_hat + 1e-8) worse <- worse + 1L
  }
  expect_equal(worse, 12L)

  # and the units were recovered
  expect_lt(abs(fit$phi_table$phi[fit$phi_table$group == "g2"] /
                fit$phi_table$phi[fit$phi_table$group == "g1"] - 1.5), 0.35)
})
