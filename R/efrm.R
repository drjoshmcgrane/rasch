# rasch :: extended frame of reference model
# ===========================================================================
# Humphry's extended frame of reference model (Humphry 2005; Humphry &
# Andrich 2008). A frame of reference is a class of persons responding to a
# class of items in a well-defined response context (Humphry & Andrich 2008);
# here a frame F_sg is one item-set by person-group cell, with scale
# parameter rho_sg = alpha_s * phi_g -- in the paper's terms (eq. 15) the
# RATIO of the reference unit to the frame's own unit, so a frame with
# rho > 1 has the smaller natural unit and the steeper curves:
#
#   P(X_ni = x) prop exp( rho_sg * ( x*theta_n - sum_{h<=x} delta_ih ) )
#
# Within a frame all curves are parallel, so the partial credit model holds
# in the frame's natural unit and the pairwise conditional logic
# (Zwinderman 1995) applies unchanged within frames.
#
# Humphry (2005) states the model for dichotomous responses; the polytomous
# form above is this package's extension of it. It is characterised by the
# properties the model's logic requires: the unit multiplies the whole
# exponent, so (i) within every frame the partial credit model holds in the
# frame's natural unit (natural thresholds rho*(delta - c), parallel
# curves), which is what makes the pairwise conditional cancellation valid;
# (ii) the weighted score sum_i rho_i x_i remains sufficient for the person
# parameter; and (iii) it reduces exactly to the dichotomous statement when
# every item has two categories and to the ordinary PCM when every unit is
# one. Per-threshold discriminations would destroy (i) and with it Rasch
# measurement; a unit on only part of the exponent would destroy the
# natural-unit reading. A consequence worth noting when interpreting
# results: category widths in natural units scale with the frame unit.
#
# Estimation is in two routes, each assigned to the data structure that
# identifies it:
#
# Route 1 (person-free, within-frame pairwise conditional ML): the centred,
# alpha-absorbed set thresholds dtilde (sum-zero per set) and the person
# group units phi_g, identified because a set taken by two groups shows the
# same threshold pattern at two scales. The structural map
# tau_v = phi_g(v) * dtilde is bilinear, so estimation alternates two exact
# linear Newton steps and finishes with a joint polish for exact sandwich
# standard errors. Frame origins are NOT pairwise-identifiable (an additive
# constant per frame cancels at every pair total) and item-set units
# alpha_s are exactly confounded with the spread of the set's thresholds,
# so neither is a parameter of this stage.
#
# Route 2 (person-side linking): alpha_s and set locations mu_s from persons
# common to two sets. In set s a person's natural coordinate is
# u_s = alpha_s * (theta - mu_s), hence u_b = r*u_a + c for a linked set pair.
# The response likelihood is integrated over masses estimated jointly with
# the link on a grid in u_a. This semiparametric mixing distribution can
# represent skewed, heavy-tailed and multimodal person populations without a
# prespecified normal shape.
# The pairwise log ratios and offsets are reconciled over the set-linking graph
# by weighted least squares. Truncated-score moments provide stable starting
# values and screen links whose score range is too weak to identify a scale.
#
# Identification: sum dtilde = 0 per set; sum_g log phi_g = 0;
# sum_s log alpha_s = 0; mean_s mu_s = 0. Together these give the
# product-of-units constraint over the frame grid with the arbitrary-unit
# origin at the mean set location. Cross-set item pairs carry no usable
# conditional information when units differ (the person parameter does not
# cancel) and are excluded from the pairwise stage.
# ===========================================================================

# Same-set pair filter: .pair_counts only returns pairs with overlapping
# persons, so cross-group pairs are already absent; cross-set pairs within a
# group must be removed because the person parameter does not cancel there.
.efrm_filter_pairs <- function(pairs, vmap) {
  Filter(function(pc) vmap$set[pc$i] == vmap$set[pc$j], pairs)
}

# Connected components by union-find; x is a list of integer edge pairs.
.efrm_components <- function(n, edges) {
  parent <- seq_len(n)
  find <- function(a) { while (parent[a] != a) a <- parent[a]; a }
  for (e in edges) {
    ra <- find(e[1]); rb <- find(e[2])
    if (ra != rb) parent[ra] <- rb
  }
  vapply(seq_len(n), find, 1L)
}

# Stage 1: block-coordinate Newton on the bilinear map tau = phi_g * dtilde,
# with a joint polish step and Godambe sandwich covariance.
.efrm_solve <- function(Xv, thr_v, m_v, vmap, pairs, drow, A_D,
                        maxit = 50, tol = 1e-7) {
  Mv <- nrow(thr_v)
  glevs <- sort(unique(vmap$group))
  G <- length(glevs)
  gidx <- match(vmap$group[thr_v$item], glevs)   # group of each virtual threshold
  Pd <- ncol(A_D)

  # one inner Newton block on a linear design tau = off + B beta
  newton_block <- function(B, off, beta, positive = FALSE, inner = 4) {
    glh <- .pcml_glh(drop(off + B %*% beta), thr_v, pairs, m_v)
    for (it in seq_len(inner)) {
      gb <- drop(crossprod(B, glh$g))
      Hb <- crossprod(B, glh$H %*% B)
      step <- tryCatch(solve(Hb, gb), error = function(e)
        solve(Hb - diag(1e-8, nrow(Hb)), gb))
      lam <- 1; ok <- FALSE; g2 <- glh
      for (half in 1:30) {
        cand <- beta - lam * step
        if (positive && any(cand <= 0)) { lam <- lam / 2; next }
        g2 <- .pcml_glh(drop(off + B %*% cand), thr_v, pairs, m_v)
        if (is.finite(g2$ll) && g2$ll >= glh$ll - 1e-12) { ok <- TRUE; break }
        lam <- lam / 2
      }
      if (!ok) break
      beta <- cand; glh <- g2
      if (max(abs(lam * step)) < tol) break
    }
    list(beta = beta, ll = glh$ll, glh = glh)
  }

  # start values: log-ratio least squares pooled over groups, phi = 1
  st <- .start_tau(Xv, thr_v)
  dt0 <- vapply(seq_len(max(drow)), function(d) mean(st[drow == d]), 0)
  # recentre per set
  for (s in unique(vmap$set)) {
    rows <- unique(drow[vmap$set[thr_v$item] == s])
    dt0[rows] <- dt0[rows] - mean(dt0[rows])
  }
  beta_d <- qr.coef(qr(A_D), dt0); beta_d[is.na(beta_d)] <- 0
  phi <- rep(1, G)

  ll <- -Inf
  for (outer in seq_len(maxit)) {
    # delta step: tau = B_d beta_d with B_d = phi_g(v) * A_D[drow(v), ]
    B_d <- A_D[drow, , drop = FALSE] * phi[gidx]
    res_d <- newton_block(B_d, 0, beta_d)
    beta_d <- res_d$beta
    dtil <- drop(A_D %*% beta_d)
    # phi step (skip when G = 1): tau = off + B_r phi[-1]
    if (G > 1L) {
      off <- dtil[drow] * (gidx == 1L)
      B_r <- matrix(0, Mv, G - 1L)
      for (g in 2:G) B_r[gidx == g, g - 1L] <- dtil[drow][gidx == g]
      res_r <- newton_block(B_r, off, phi[-1], positive = TRUE)
      phi <- c(1, res_r$beta)
      ll_new <- res_r$ll
    } else ll_new <- res_d$ll
    if (is.finite(ll) && abs(ll_new - ll) < tol * (abs(ll_new) + 1)) {
      ll <- ll_new; break
    }
    ll <- ll_new
  }

  # joint polish with the exact bilinear Hessian, then sandwich covariance
  tau_hat <- phi[gidx] * dtil[drow]
  glh <- .pcml_glh(tau_hat, thr_v, pairs, m_v)
  build_J <- function(dtil, phi) {
    J <- matrix(0, Mv, Pd + G - 1L)
    J[, seq_len(Pd)] <- A_D[drow, , drop = FALSE] * phi[gidx]
    if (G > 1L) for (g in 2:G)
      J[gidx == g, Pd + g - 1L] <- dtil[drow][gidx == g]
    J
  }
  struct_H <- function(J, glh) {
    H <- crossprod(J, glh$H %*% J)
    if (G > 1L) for (g in 2:G) {
      rows <- which(gidx == g)
      cb <- drop(crossprod(A_D[drow[rows], , drop = FALSE], glh$g[rows]))
      H[seq_len(Pd), Pd + g - 1L] <- H[seq_len(Pd), Pd + g - 1L] + cb
      H[Pd + g - 1L, seq_len(Pd)] <- H[Pd + g - 1L, seq_len(Pd)] + cb
    }
    H
  }
  for (polish in 1:10) {
    J <- build_J(dtil, phi)
    H <- struct_H(J, glh)
    g_full <- drop(crossprod(J, glh$g))
    step <- tryCatch(solve(H, g_full), error = function(e)
      solve(H - diag(1e-8, nrow(H)), g_full))
    cand_d <- beta_d - step[seq_len(Pd)]
    cand_p <- if (G > 1L) phi[-1] - step[Pd + seq_len(G - 1L)] else numeric(0)
    if (G > 1L && any(cand_p <= 0)) break
    dt_c <- drop(A_D %*% cand_d); phi_c <- c(1, cand_p)
    glh_c <- .pcml_glh(phi_c[gidx] * dt_c[drow], thr_v, pairs, m_v)
    if (!is.finite(glh_c$ll) || glh_c$ll < glh$ll - 1e-10) break
    beta_d <- cand_d; dtil <- dt_c; phi <- phi_c; glh <- glh_c
    if (max(abs(step)) < tol) break
  }
  tau_hat <- phi[gidx] * dtil[drow]

  J <- build_J(dtil, phi)
  H <- struct_H(J, glh)
  # identification of the group units: a (near-)null direction of the
  # UNRIDGED joint information that loads on a group's unit parameter
  # means the data cannot identify that unit -- the ridged solve would
  # land anywhere on the flat manifold and still show a small gradient.
  # (H is the negative expected Hessian here, so its eigenvalues are
  # nonnegative up to the bilinear cross term.)
  phi_unident <- setNames(rep(FALSE, G), glevs)
  if (G > 1L) {
    eh <- eigen(H, symmetric = TRUE)
    flat <- abs(eh$values) < max(abs(eh$values)) * 1e-10
    if (any(flat)) {
      V <- eh$vectors[, flat, drop = FALSE]
      lp_load <- sqrt(rowSums(V[Pd + seq_len(G - 1L), , drop = FALSE]^2))
      phi_unident[-1L][lp_load > 1e-2] <- TRUE
    }
  }
  Hinv <- tryCatch(solve(H), error = function(e) solve(H - diag(1e-8, nrow(H))))
  Jt <- .pcml_sandwich(Xv, thr_v, m_v, tau_hat, pairs)
  Jb <- crossprod(J, Jt %*% J)
  support <- attr(Jt, "cluster_support", exact = TRUE)
  support_guard <- .pcml_covariance_support(Jb, ncol(J), support)
  covb <- Hinv %*% Jb %*% Hinv
  if (!support_guard$inference) covb[,] <- NA_real_
  conv <- max(abs(drop(crossprod(J, glh$g)))) < 1e-3

  # recentre to sum_g log phi = 0; dtilde absorbs the constant
  cc <- mean(log(phi))
  phi_c <- phi * exp(-cc); dtil_c <- dtil * exp(cc)
  # covariance of log phi (ref group has zero variance), centred
  cov_lp <- matrix(0, G, G)
  if (G > 1L) {
    idx <- Pd + seq_len(G - 1L)
    Dl <- diag(1 / phi[-1], G - 1L)
    cov_lp[2:G, 2:G] <- Dl %*% covb[idx, idx, drop = FALSE] %*% Dl
  }
  Ac <- diag(G) - matrix(1 / G, G, G)
  cov_lp <- Ac %*% cov_lp %*% t(Ac)
  # dtilde_c = exp(cc) A_D beta_d, where
  # cc = mean(log(phi)).  Because cc is estimated jointly with beta_d, the
  # delta-method Jacobian must include its derivatives with respect to every
  # free phi as well as the beta--phi cross-covariance. Treating cc as fixed
  # can materially misstate the common-unit item standard errors.
  J_dt <- matrix(0, nrow(A_D), Pd + G - 1L)
  J_dt[, seq_len(Pd)] <- exp(cc) * A_D
  if (G > 1L)
    J_dt[, Pd + seq_len(G - 1L)] <-
      outer(dtil_c, 1 / (G * phi[-1L]))
  cov_dt <- J_dt %*% covb %*% t(J_dt)

  # joint covariance of (dtilde, log phi): both are functions of the same
  # jointly estimated (beta_d, phi), so downstream propagation (the
  # linking-stage calibration redraws) must draw them together -- the
  # cross-covariance is part of the estimate, not an approximation choice
  J_lp <- matrix(0, G, Pd + G - 1L)
  if (G > 1L) J_lp[2:G, Pd + seq_len(G - 1L)] <- diag(1 / phi[-1], G - 1L)
  J_lp <- Ac %*% J_lp
  J_all <- rbind(J_dt, J_lp)
  cov_joint <- J_all %*% covb %*% t(J_all)

  dimnames(cov_lp) <- list(glevs, glevs)
  list(dtilde = dtil_c, phi = setNames(phi_c, glevs),
       se_log_phi = setNames(sqrt(pmax(diag(cov_lp), 0)), glevs),
       cov_log_phi = cov_lp, phi_unident = phi_unident,
       cov_dtilde = cov_dt, cov_joint = cov_joint,
       loglik = glh$ll, iterations = outer,
       converged = conv, gidx = gidx,
       cluster_inference = support_guard$inference,
       cluster_support = support, cluster_note = support_guard$note)
}

# Refit the within-frame conditional likelihood under equal group units.  This
# is a separate optimiser from the EFRM calibration, so its likelihood is a
# usable comparison only when this solve itself reaches an identified optimum.
.efrm_equal_fit <- function(dtil, A_D, drow, thr_v, pairs, m_v,
                            maxit = 25L, tol = 1e-7) {
  B <- A_D[drow, , drop = FALSE]
  beta <- qr.coef(qr(A_D), dtil)
  beta[is.na(beta)] <- 0
  glh <- .pcml_glh(drop(B %*% beta), thr_v, pairs, m_v)
  iterations <- 0L
  for (it in seq_len(maxit)) {
    iterations <- it
    gb <- drop(crossprod(B, glh$g))
    Hb <- crossprod(B, glh$H %*% B)
    step <- tryCatch(solve(Hb, gb), error = function(e)
      tryCatch(solve(Hb - diag(1e-8, nrow(Hb)), gb),
               error = function(e2) NULL))
    if (is.null(step) || any(!is.finite(step))) break
    lam <- 1
    moved <- FALSE
    for (half in seq_len(30L)) {
      candidate <- beta - lam * step
      next_glh <- .pcml_glh(drop(B %*% candidate), thr_v, pairs, m_v)
      if (is.finite(next_glh$ll) && next_glh$ll >= glh$ll - 1e-12) {
        beta <- candidate
        glh <- next_glh
        moved <- TRUE
        break
      }
      lam <- lam / 2
    }
    if (!moved || max(abs(lam * step)) < tol) break
  }

  gb <- drop(crossprod(B, glh$g))
  Hb <- crossprod(B, glh$H %*% B)
  rc <- tryCatch(rcond(Hb), error = function(e) 0)
  rank_ok <- is.finite(rc) && rc > 1e-12
  move <- if (rank_ok)
    tryCatch(drop(solve(Hb, gb)), error = function(e) rep(Inf, length(gb)))
  else rep(Inf, length(gb))
  move_tol <- min(20 * tol, 1e-6)
  converged <- is.finite(glh$ll) && rank_ok && all(is.finite(gb)) &&
    (max(abs(gb)) < 1e-4 ||
       (all(is.finite(move)) && max(abs(move)) < move_tol))
  list(loglik = glh$ll, converged = converged, rank_ok = rank_ok,
       iterations = iterations, score = gb, newton_move = move)
}

# Per-person correction moments for the set-unit linking. For each
# missing-data pattern the WLE score map W and the model score
# distribution are exact, so the mean g(u) and variance w(u) of W(R) over
# the non-extreme scores are exact functions of u given the thresholds.
# The linking needs E[w(u)] and E[g'(u)] over the (unobserved) person
# distribution; because the raw score is sufficient, score-weight vectors
# phi with E[phi(R) | u] ~ target(u) on a grid make mean(phi(R_i)) an
# estimate of E[target(u)] over a broad class of person distributions. These
# moments now provide starting values and identification checks for the
# semiparametric likelihood link rather than the final set-unit estimate.
# Returns
# per-person phi_w(R_i) and phi_g(R_i), NA for extreme or empty rows.
.person_link_moments <- function(X, tau_list, disc = 1) {
  N <- nrow(X)
  obs <- !is.na(X)
  pat <- apply(obs, 1, function(z) paste(which(z), collapse = ","))
  w_i <- g_i <- rep(NA_real_, N)
  for (key in unique(pat)) {
    cols <- as.integer(strsplit(key, ",", fixed = TRUE)[[1]])
    if (length(cols) < 2L) next
    tl <- tau_list[cols]
    maxr <- sum(vapply(tl, length, 1L))
    # a pattern needs a score range of at least 4 (>= 3 interior score
    # categories) to provide a stable scale start and guard check. With only
    # 2 interior categories the score-map slope is too weak to support the
    # flexible distributional link reliably.
    if (maxr < 4L) next
    W <- person_wle(tl, disc = disc)$theta[as.character(0:maxr)]
    Wi <- W[-c(1L, maxr + 1L)]
    if (any(!is.finite(Wi))) next
    grid <- seq(min(Wi) - 3, max(Wi) + 3, length.out = 121L)
    # truncated score distribution at each grid point, by direct
    # convolution of the item category probabilities (same model as
    # person_wle; direct is faster than FFT at these lengths)
    B <- vapply(grid, function(u) {
      p <- 1
      for (tau in tl) {
        pi <- item_moments(u, tau, disc = disc)$P
        np <- numeric(length(p) + length(pi) - 1L)
        for (x in seq_along(pi)) {
          ix <- x:(x + length(p) - 1L)
          np[ix] <- np[ix] + p * pi[x]
        }
        p <- np
      }
      p <- pmax(p, 0)[-c(1L, maxr + 1L)]
      p / sum(p)
    }, numeric(maxr - 1L))
    gm <- colSums(B * Wi)
    wm <- colSums(B * Wi^2) - gm^2
    h <- grid[2L] - grid[1L]
    gp <- c(gm[2L] - gm[1L], (gm[-c(1L, 2L)] - gm[-c(120L, 121L)]) / 2,
            gm[121L] - gm[120L]) / h
    A2 <- t(B)
    M <- crossprod(A2) + diag(1e-6, maxr - 1L)
    phi_w <- drop(solve(M, crossprod(A2, wm)))
    phi_g <- drop(solve(M, crossprod(A2, gp)))
    sel <- which(pat == key)
    r <- rowSums(X[sel, cols, drop = FALSE])
    inner <- r >= 1 & r <= maxr - 1L
    w_i[sel[inner]] <- phi_w[r[inner]]
    g_i[sel[inner]] <- phi_g[r[inner]]
  }
  list(w = w_i, g = g_i)
}

# Reference kernels retained for numerical parity tests and platforms on which
# the compiled path is deliberately disabled with options(rasch.efrm_cpp=FALSE).
.efrm_npml_likelihood_r <- function(u, obs, score, taus, discs) {
  L <- outer(score, u)
  patterns <- apply(obs, 1L, paste0, collapse = "")
  for (pattern in unique(patterns)) {
    rows <- which(patterns == pattern)
    jj <- which(obs[rows[1L], ])
    den <- numeric(length(u))
    for (j in jj) {
      tt <- taus[[j]]; x <- 0:length(tt)
      lp <- discs[j] * (outer(u, x) -
        matrix(c(0, cumsum(tt)), nrow = length(u),
               ncol = length(x), byrow = TRUE))
      mx <- apply(lp, 1L, max)
      den <- den + mx + log(rowSums(exp(lp - mx)))
    }
    L[rows, ] <- L[rows, , drop = FALSE] -
      matrix(den, nrow = length(rows), ncol = length(u), byrow = TRUE)
  }
  L
}

.efrm_npml_fit_weights_r <- function(L, logw, mix_idx, count,
                                     maxit = 100L, tol = 1e-7) {
  conv <- FALSE; step <- Inf; ll_step <- Inf; ll <- -Inf; iterations <- 0L
  check_convergence <- is.finite(tol) && tol > 0
  observed_loglik <- function(current) {
    A <- L + current[mix_idx, , drop = FALSE]
    mx <- apply(A, 1L, max)
    sum(count * (mx + log(rowSums(exp(A - mx)))))
  }
  for (it in seq_len(maxit)) {
    A <- L + logw[mix_idx, , drop = FALSE]
    mx <- apply(A, 1L, max)
    post <- exp(A - mx); post <- post / rowSums(post)
    w <- matrix(NA_real_, nrow(logw), ncol(logw))
    for (h in seq_len(nrow(logw))) {
      sel <- mix_idx == h
      wh <- pmax(colSums(post[sel, , drop = FALSE] * count[sel]) /
                   sum(count[sel]), 1e-10)
      w[h, ] <- wh / sum(wh)
    }
    next_logw <- log(w)
    step <- max(abs(w - exp(logw)))
    logw <- next_logw
    iterations <- it
    if (check_convergence) {
      ll_new <- observed_loglik(logw)
      ll_step <- ll_new - ll
      ll <- ll_new
      # Adjacent grid masses can move slowly along an almost flat likelihood
      # ridge. Convergence of the observed likelihood is the relevant EM
      # criterion; requiring every mass itself to stop moving falsely labels a
      # stable set transformation as unconverged.
      if (is.finite(ll_step) && abs(ll_step) <= tol * (1 + abs(ll))) {
        conv <- TRUE
        break
      }
    }
  }
  # Coordinate-ascent calls deliberately use tol = 0 for a fixed number of
  # EM updates. Their intermediate likelihood is not used, so calculate it
  # once after the updates rather than once per iteration.
  if (!check_convergence) ll <- observed_loglik(logw)
  list(logw = logw, converged = conv, step = step,
       loglik = ll, loglik_step = ll_step, iterations = iterations)
}

# Text used only to identify repeated sufficient statistics before the NPML
# likelihood is evaluated. Seventeen significant digits round-trip an R
# double; fewer can pool distinct weighted scores and thereby change the
# likelihood when two fitted frame units are close.
.efrm_score_key <- function(x) sprintf("%.17g", x)

# A joint response is unbounded on the linking grid only when both sets point
# to the same tail. Opposing extremes (minimum in one set, maximum in the
# other) have a finite compromise location and therefore remain informative
# about the set transformation.
.efrm_same_tail_extreme <- function(score_a, max_a, score_b, max_b,
                                    tol = 1e-8) {
  # Work on the attainable-score proportion, rather than an absolute score
  # tolerance. Frame discriminations can be arbitrarily small or large; an
  # absolute tolerance can therefore call an interior response pattern an
  # extreme solely because of its fitted unit.
  valid_a <- is.finite(score_a) & is.finite(max_a) & max_a > 0
  valid_b <- is.finite(score_b) & is.finite(max_b) & max_b > 0
  prop_a <- prop_b <- rep(NA_real_, length(score_a))
  prop_a[valid_a] <- score_a[valid_a] / max_a[valid_a]
  prop_b[valid_b] <- score_b[valid_b] / max_b[valid_b]
  low_a <- valid_a & prop_a <= tol
  low_b <- valid_b & prop_b <= tol
  high_a <- valid_a & prop_a >= 1 - tol
  high_b <- valid_b & prop_b >= 1 - tol
  (low_a & low_b) | (high_a & high_b)
}

# Semiparametric likelihood link for two item sets. The common persons'
# distribution is represented by jointly estimated masses on a fixed grid in
# the first set's natural unit; the second set is evaluated at r*u + c. This
# retains the Rasch response model without prescribing a normal or unimodal
# person distribution. Response rows are compressed by
# observed-item pattern and weighted score before the EM iterations, so the
# link remains practical inside the person bootstrap.
.efrm_npml_pair <- function(Xm, vmap, tau_v, disc_v, sets_u, a, b, idx,
                            init_log_ratio, init_offset,
                            min_link_persons, grid_n = 61L,
                            report_support = FALSE) {
  # A missing edge does not make a connected set graph unidentified. Keep
  # insufficient response support distinct from a numerical link failure so
  # the graph solver can omit only the former. Direct numerical callers
  # retain the historical NULL result for an unavailable pair.
  unsupported <- function(reason) {
    if (!report_support) return(NULL)
    structure(list(reason = reason), class = "rasch_efrm_unlinked_pair")
  }
  ca <- which(vmap$set == sets_u[a]); cb <- which(vmap$set == sets_u[b])
  Xa <- Xm[idx, ca, drop = FALSE]; Xb <- Xm[idx, cb, drop = FALSE]
  oa <- !is.na(Xa); ob <- !is.na(Xb)
  keep <- rowSums(oa) > 0L & rowSums(ob) > 0L
  Xa <- Xa[keep, , drop = FALSE]; Xb <- Xb[keep, , drop = FALSE]
  oa <- oa[keep, , drop = FALSE]; ob <- ob[keep, , drop = FALSE]
  if (nrow(Xa) < min_link_persons)
    return(unsupported("too few common persons"))

  da <- disc_v[ca]; db <- disc_v[cb]
  score_a <- rowSums(sweep(Xa, 2L, da, "*"), na.rm = TRUE)
  score_b <- rowSums(sweep(Xb, 2L, db, "*"), na.rm = TRUE)
  pa <- apply(oa, 1L, paste0, collapse = "")
  pb <- apply(ob, 1L, paste0, collapse = "")
  key <- paste(pa, .efrm_score_key(score_a),
               pb, .efrm_score_key(score_b), sep = "\r")
  lev <- unique(key); grp <- match(key, lev)
  reps <- match(lev, key); count <- tabulate(grp, nbins = length(lev))
  oa <- oa[reps, , drop = FALSE]; ob <- ob[reps, , drop = FALSE]
  score_a <- score_a[reps]; score_b <- score_b[reps]
  # The flexible link needs an adequate attainable-score range in both sets.
  # This is a property of the observed-item pattern, not of the fitted frame
  # unit. Persons at opposite extremes of two adequate patterns remain
  # informative, whereas concordant extremes do not identify their relative
  # transformation.
  range_a <- rowSums(sweep(oa, 2L, lengths(tau_v[ca]), "*"))
  range_b <- rowSums(sweep(ob, 2L, lengths(tau_v[cb]), "*"))
  max_a <- rowSums(sweep(oa, 2L, da * lengths(tau_v[ca]), "*"))
  max_b <- rowSums(sweep(ob, 2L, db * lengths(tau_v[cb]), "*"))
  extreme_both <- .efrm_same_tail_extreme(
    score_a, max_a, score_b, max_b, 1e-8)
  informative <- count * (range_a >= 4L & range_b >= 4L & !extreme_both)
  n_informative <- sum(informative)
  if (n_informative < min_link_persons)
    return(unsupported("too few informative common score patterns"))
  # Group membership is observed and the person populations need not have the
  # same location, spread or shape. Keep one nonparametric margin per observed
  # group while estimating a common set transformation. A single pooled margin
  # would impose an avoidable population-distribution restriction whenever the
  # fitted group units differ.
  mix_group <- vapply(seq_len(nrow(oa)), function(i) {
    ja <- which(oa[i, ])
    jb <- which(ob[i, ])
    ga <- if (length(ja)) unique(vmap$group[ca[ja]]) else character(0)
    gb <- if (length(jb)) unique(vmap$group[cb[jb]]) else character(0)
    gg <- unique(c(ga, gb))
    if (length(gg) != 1L) stop("internal EFRM link mixes person groups")
    gg
  }, character(1L))
  mix_levels <- sort(unique(mix_group))
  mix_idx <- match(mix_group, mix_levels)

  # Log likelihood up to response-pattern constants, which cancel from the
  # posterior masses and do not depend on r or c. The compiled kernel is the
  # production path; the R implementation above remains an executable
  # numerical reference.
  use_cpp <- !identical(getOption("rasch.efrm_cpp", TRUE), FALSE)
  likelihood <- if (use_cpp) efrm_likelihood_cpp else .efrm_npml_likelihood_r

  grid <- seq(-8, 8, length.out = as.integer(grid_n))
  La <- likelihood(grid, oa, score_a, tau_v[ca], da)
  par <- c(init_log_ratio, init_offset)
  if (any(!is.finite(par))) par <- c(0, 0)
  # Joint NPML by coordinate ascent. For a fixed link, the EM step estimates
  # the grid masses from both sets; for fixed masses, a bounded likelihood
  # step estimates the log scale ratio and offset. Using both sets in the
  # mixing-distribution step matters: estimating the masses from the first
  # set alone would turn the second into a plug-in prediction problem rather
  # than maximise the likelihood of the linked responses.
  add_weights <- function(L, logw)
    L + logw[mix_idx, , drop = FALSE]
  fit_weights <- function(L, logw, maxit = 100L, tol = 0) {
    if (use_cpp)
      efrm_fit_weights_cpp(L, logw, mix_idx, count, maxit, tol)
    else
      .efrm_npml_fit_weights_r(L, logw, mix_idx, count, maxit, tol)
  }
  logw <- matrix(-log(length(grid)), length(mix_levels), length(grid),
                 dimnames = list(mix_levels, NULL))
  converged <- FALSE
  op <- NULL
  last_step <- Inf
  for (outer in seq_len(20L)) {
    old <- par
    Lb <- likelihood(exp(par[1L]) * grid + par[2L], ob, score_b,
                     tau_v[cb], db)
    ew <- fit_weights(La + Lb, logw)
    logw <- ew$logw
    objective <- if (use_cpp) function(z)
      efrm_negloglik_cpp(z, grid, ob, score_b, tau_v[cb], db,
                         La, logw, mix_idx, count)
    else function(z) {
      Lbz <- likelihood(exp(z[1L]) * grid + z[2L], ob, score_b,
                        tau_v[cb], db)
      A <- add_weights(La + Lbz, logw); mx <- apply(A, 1L, max)
      -sum(count * (mx + log(rowSums(exp(A - mx)))))
    }
    lower <- c(log(0.1), -20); upper <- c(log(10), 20)
    op <- stats::optim(par, objective, method = "L-BFGS-B",
                       lower = lower, upper = upper,
                       control = list(maxit = 30L, factr = 1e7))
    par <- op$par
    last_step <- max(abs(par - old))
    if (op$convergence == 0L && last_step < 1e-3) {
      converged <- TRUE; break
    }
  }
  # Refresh the nuisance distribution at the reported link and evaluate the
  # corresponding joint log likelihood.
  Lb <- likelihood(exp(par[1L]) * grid + par[2L], ob, score_b,
                   tau_v[cb], db)
  final_mass_maxit <- getOption("rasch.efrm_npml_mass_maxit", 500L)
  final_mass_maxit <- .check_whole(
    final_mass_maxit, "option rasch.efrm_npml_mass_maxit", 1)
  ew <- fit_weights(La + Lb, logw, final_mass_maxit, 1e-7)
  logw <- ew$logw; w <- exp(logw)
  A <- add_weights(La + Lb, logw); mx <- apply(A, 1L, max)
  ll <- sum(count * (mx + log(rowSums(exp(A - mx)))))
  converged <- (converged || (!is.null(op) && op$convergence == 0L &&
                              last_step < 1e-3)) && isTRUE(ew$converged)
  # Mass piled on the ends of the grid signals a truncated person
  # distribution only when persons with a finite location put it there. A
  # person at the same extreme tail in both sets has a likelihood that rises
  # monotonically towards one end of the grid, so the masses explaining such
  # persons sit on the edge of any finite box and carry no information about
  # the link. Opposing extremes have a finite compromise location and remain
  # informative. The truncation check therefore uses the posterior edge mass
  # of the remaining persons. The extreme persons stay in the
  # likelihood: dropping them would condition on the realised responses
  # without including that selection in the likelihood.
  post <- exp(A - mx)
  post_edge <- (post[, 1L] + post[, ncol(post)]) / rowSums(post)
  informative_by_group <- rowsum(informative, mix_idx)
  edge_mass_by_group <- rowsum(informative * post_edge, mix_idx) /
    informative_by_group
  # A group represented only by concordant extreme responses carries no
  # information about this set link. It must not manufacture a zero edge
  # mass, but another observed group can still identify the common link.
  edge_mass_by_group <- edge_mass_by_group[drop(informative_by_group) > 0, ,
                                           drop = FALSE]
  edge_mass <- max(edge_mass_by_group)
  if (isTRUE(getOption("rasch.efrm_link_debug", FALSE)))
    message("EFRM NPML: step=", signif(last_step, 4),
            ", mass step=", signif(ew$step, 4),
            ", mass loglik step=", signif(ew$loglik_step, 4),
            ", mass iterations=", ew$iterations,
            ", optim=", op$convergence,
            ", max group edge mass=", signif(edge_mass, 4),
            ", total edge mass=",
            signif(max(w[, 1L] + w[, ncol(w)]), 4),
            ", persons at the same extreme tail in both sets=",
            sum(count[extreme_both]))
  # L-BFGS-B regards an optimum on its artificial search boundary as a
  # successful fit. Such a value says only that the data want a ratio or
  # offset beyond the numerical box; it is not a finite interior estimate.
  boundary_tol <- 1e-6
  on_boundary <- any(par <= lower + boundary_tol |
                       par >= upper - boundary_tol)
  if (!all(is.finite(par)) || !is.finite(ll) || edge_mass > 0.02 ||
      on_boundary) return(NULL)
  list(log_ratio = par[1L], offset = par[2L], n = n_informative,
       converged = converged, edge_mass = edge_mass, loglik = ll)
}

# Stage 2: alpha_s and set locations mu_s from persons common to set pairs,
# with semiparametric person-distribution linking. Score moments provide
# stable starting values and guard checks; the pair likelihood supplies the
# reported scale ratio and offset.
# Standard errors and the joint covariance of (log alpha, mu) come from a
# person bootstrap of the whole linking stage: each replicate resamples
# persons, refits the nuisance masses and link for each set pair, and
# re-solves the linking least squares. This captures the nonlinearity of the
# semiparametric link and its reconciliation over the linking graph.
.rasch_cancelled <- function(label = "Estimation") {
  cond <- structure(
    list(message = paste(label, "cancelled"), call = NULL),
    class = c("rasch_cancelled", "error", "condition"))
  stop(cond)
}
.efrm_cancelled <- function() .rasch_cancelled("EFRM estimation")

# A covariance estimate needs three things at once: at least 30 usable
# draws (the documented contract -- at the permitted minimum of boot_reps =
# 30 that means every replicate, which is what "at least 30 are required"
# promises), a majority of the draws requested (a badly failing design must
# not look adequately sampled by asking for more), and more draws than the
# covariance has independent directions (fewer cannot span its effective
# rank). An earlier reading
# argued the 30-floor away as too strict at the minimum; 16 draws pricing a
# covariance is what that reasoning permitted.
.rasch_min_boot_success <- function(boot_reps, n_quantities = 0L) {
  # the documented contract (rasch_btl_efrm, rasch_efrm): inference needs at
  # least 30 successful replicates AND more than half of those requested. A
  # bare majority of 30 would let 16 draws price a covariance, and a
  # standard deviation from 16 draws is noise wearing a number. A covariance
  # additionally needs more draws than its effective rank, or the sample
  # covariance cannot span its free directions -- so the floor also clears
  # that count.
  as.integer(max(30L, floor(boot_reps / 2) + 1L, n_quantities + 1L))
}
.efrm_min_boot_success <- function(boot_reps, n_quantities = 0L) {
  .rasch_min_boot_success(boot_reps, n_quantities)
}

# Apply bootstrap jobs either serially or on a persistent socket cluster.
# The caller supplies random draws or per-replicate seeds before dispatch;
# workers use the coordinator's RNG settings when generating from those seeds.
# Small batches retain useful progress and cancellation checkpoints without
# repeatedly starting worker processes.
.rasch_available_workers <- function() {
  cores <- suppressWarnings(parallel::detectCores(logical = FALSE))
  if (length(cores) != 1L || !is.finite(cores) || cores < 1L)
    cores <- suppressWarnings(parallel::detectCores(logical = TRUE))
  if (length(cores) != 1L || !is.finite(cores) || cores < 1L) cores <- 4L

  limits <- as.integer(cores)
  env_limits <- Sys.getenv(c("SLURM_CPUS_PER_TASK", "PBS_NP", "NSLOTS",
                             "OMP_THREAD_LIMIT"), unset = NA_character_)
  env_limits <- suppressWarnings(as.integer(env_limits))
  limits <- c(limits, env_limits[is.finite(env_limits) & env_limits > 0L])

  worker_limit <- function(x) {
    if (length(x) != 1L || !is.numeric(x) || is.complex(x) ||
        !is.null(dim(x)) || !is.null(oldClass(x)) || !is.finite(x) ||
        x < 1L || x != floor(x) || x > .Machine$integer.max)
      return(NA_integer_)
    as.integer(x)
  }
  option_limits <- vapply(
    list(getOption("rasch.max_workers", NA_integer_),
         getOption("rasch.efrm.max_workers", NA_integer_)),
    worker_limit, integer(1))
  limits <- c(limits,
              option_limits[is.finite(option_limits) & option_limits >= 1L])

  check_limit <- tolower(Sys.getenv("_R_CHECK_LIMIT_CORES_", unset = ""))
  if (nzchar(check_limit) && !check_limit %in% c("false", "no", "0"))
    limits <- c(limits, 2L)

  max(1L, min(limits))
}
.efrm_available_workers <- function() .rasch_available_workers()

.rasch_namespace_is_installed <- function() {
  # A sourced copy can coexist with an older installed release. In that case
  # socket workers must not load the installed namespace and run different
  # code from their coordinator.
  env <- environment(.rasch_boot_apply)
  if (!isNamespace(env) ||
      !identical(unname(getNamespaceName(env)), "rasch"))
    return(FALSE)
  package_dir <- system.file(package = "rasch")
  nzchar(package_dir) && file.exists(file.path(package_dir, "DESCRIPTION"))
}

.rasch_boot_apply <- function(n, fun, workers = 1L, progress = NULL,
                              cancel = NULL, label = "Bootstrap estimation") {
  out <- vector("list", n)
  package_dir <- system.file(package = "rasch")
  if (workers > 1L && !.rasch_namespace_is_installed()) {
    # Under pkgload::load_all(), system.file() points at source/inst rather
    # than an installed package root. A socket worker cannot load that source
    # namespace and would otherwise find any older installed rasch version on
    # .libPaths(), making the coordinator and workers run different code.
    warning(label, " is running serially because the current rasch namespace ",
            "was loaded from a source tree; install this tree to validate ",
            "parallel execution", call. = FALSE)
    workers <- 1L
  }
  # Two jobs per worker reduces coordination overhead for cheap refits while
  # keeping the batches small enough that an unusually slow resample does not
  # leave the other workers idle for long. It also retains frequent progress
  # and cancellation checks.
  batch_size <- 2L * workers
  batches <- split(seq_len(n), ceiling(seq_len(n) / batch_size))
  if (workers == 1L) {
    for (ids in batches) {
      if (is.function(cancel) && isTRUE(cancel())) .rasch_cancelled(label)
      out[ids] <- lapply(ids, fun)
      if (is.function(progress)) progress(max(ids), n)
    }
    return(out)
  }

  cl <- tryCatch(parallel::makePSOCKcluster(workers),
                 error = function(e) e)
  if (inherits(cl, "error")) {
    # A machine that cannot open sockets can still run the same jobs. The
    # supplied draws or replicate seeds give the same serial result.
    warning("could not start ", label, " workers (",
            conditionMessage(cl), "); running serially", call. = FALSE)
    for (ids in batches) {
      if (is.function(cancel) && isTRUE(cancel())) .rasch_cancelled(label)
      out[ids] <- lapply(ids, fun)
      if (is.function(progress)) progress(max(ids), n)
    }
    return(out)
  }
  on.exit(try(parallel::stopCluster(cl), silent = TRUE), add = TRUE)
  paths <- unique(c(dirname(package_dir), .libPaths()))
  # A replicate's integer seed identifies the same draws only under the same
  # uniform, normal and sampling generators. Socket workers start with R's
  # defaults rather than inheriting RNGkind() from the coordinator.
  rng_kind <- RNGkind()
  setup_worker <- function(paths, rng_kind) {
    .libPaths(paths)
    ns <- loadNamespace("rasch")
    if (!identical(normalizePath(getNamespaceInfo(ns, "path")),
                   normalizePath(file.path(paths[1L], "rasch"))))
      stop("bootstrap worker loaded a different rasch installation; ",
           "restart R with the analysis library first on .libPaths()")
    do.call(RNGkind, as.list(rng_kind))
    invisible(NULL)
  }
  # Kept source references can carry the package namespace into the socket
  # message. Strip them before dispatch, or unserialisation may load another
  # installed rasch before this callback can set the coordinator's library.
  setup_worker <- utils::removeSource(setup_worker)
  environment(setup_worker) <- baseenv()
  parallel::clusterCall(cl, setup_worker, paths, rng_kind)
  holder <- list2env(list(.rasch_boot_fun = fun), parent = emptyenv())
  parallel::clusterExport(cl, ".rasch_boot_fun", envir = holder)
  worker_call <- function(i)
    get(".rasch_boot_fun", envir = .GlobalEnv, inherits = FALSE)(i)
  environment(worker_call) <- baseenv()

  for (ids in batches) {
    if (is.function(cancel) && isTRUE(cancel())) .rasch_cancelled(label)
    out[ids] <- parallel::parLapply(cl, ids, worker_call)
    if (is.function(progress)) progress(max(ids), n)
  }
  out
}

.efrm_boot_apply <- function(n, fun, workers = 1L, progress = NULL,
                             cancel = NULL)
  .rasch_boot_apply(n, fun, workers, progress, cancel,
                    label = "EFRM bootstrap")

.efrm_link_sets <- function(u_mat, w_mat, g_mat, sets_u, min_link_persons,
                            boot_reps = 300, regen = NULL, pair_link = NULL,
                            progress = NULL, cancel = NULL, workers = 1L) {
  S <- ncol(u_mat)
  A <- rbind(diag(S - 1L), rep(-1, S - 1L))

  # one pass over a person index vector: per-edge corrected slopes and
  # offsets, then the two weighted least-squares solves. um/wm/gm default
  # to the point-estimate person matrices; the bootstrap can pass
  # regenerated ones (see below).
  link_once <- function(idx, um = u_mat, wm = w_mat, gm = g_mat,
                        hard = FALSE, pair_fun = pair_link) {
    edges <- list(); ls_est <- off_est <- off_n <- numeric(0)
    link_converged <- link_edge_mass <- link_loglik <- numeric(0)
    for (a in seq_len(S - 1)) for (b in (a + 1):S) {
      # Non-extreme scores are needed only for the corrected-moment start and
      # the weak-link screen. The likelihood itself must retain every common
      # person, including an extreme score in either set: dropping them would
      # condition on the realised responses without including that selection
      # in the likelihood.
      ok_start <- idx[is.finite(um[idx, a]) & is.finite(um[idx, b]) &
                        is.finite(wm[idx, a]) & is.finite(wm[idx, b]) &
                        is.finite(gm[idx, a]) & is.finite(gm[idx, b])]
      if (is.null(pair_fun) && length(ok_start) < min_link_persons) next
      # Corrected score moments are a stable start for the likelihood when
      # available. They are not its support criterion: opposing extreme
      # scores have no finite WLE but do identify a finite set
      # transformation. In that case the likelihood starts at the neutral
      # transformation and applies its own support and boundary checks.
      ls <- off <- 0
      moment_ok <- length(ok_start) >= 2L
      moment_failure <- "variance"
      if (moment_ok) {
        u1 <- um[ok_start, a]; u2 <- um[ok_start, b]
        d1 <- mean(gm[ok_start, a]); d2 <- mean(gm[ok_start, b])
        moment_ok <- is.finite(d1) && is.finite(d2) && d1 >= 0.1 && d2 >= 0.1
        if (!moment_ok) moment_failure <- "slope"
        if (moment_ok) {
          v1 <- (var(u1) - mean(wm[ok_start, a])) / d1^2
          v2 <- (var(u2) - mean(wm[ok_start, b])) / d2^2
          moment_ok <- is.finite(v1) && is.finite(v2) && v1 > 0 && v2 > 0
          if (moment_ok) {
            ls <- 0.5 * (log(v2) - log(v1))          # log(alpha_b / alpha_a)
            off <- mean(u2) - exp(ls) * mean(u1)
          }
        }
      }
      if (is.null(pair_fun) && !moment_ok) {
        if (hard) {
          if (identical(moment_failure, "slope"))
            stop("the score maps of sets '", sets_u[a], "' and '", sets_u[b],
                 "' are too flat to link (degenerate map slope)")
          stop("too little true person variance to link sets '", sets_u[a],
               "' and '", sets_u[b], "'")
        }
        return(NULL)
      }
      np <- NULL
      if (!is.null(pair_fun)) {
        np <- pair_fun(a, b, idx, ls, off)
        if (inherits(np, "rasch_efrm_unlinked_pair")) next
        if (is.null(np)) {
          if (hard)
            stop("the semiparametric likelihood link failed for sets '",
                 sets_u[a], "' and '", sets_u[b],
                 "': there were too few informative common patterns (a ",
                 "score range of at least 4 is required within each set), ",
                 "or the finite transformation did not pass its numerical ",
                 "and grid-boundary checks")
          return(NULL)
        }
        # Retain a usable point estimate with an explicit warning, but do not
        # admit a replicate with a failed supported link to the covariance.
        # A redundant edge still belongs to the fitted statistic; deleting it
        # would silently substitute a different linking graph for that draw.
        if (!hard && !isTRUE(np$converged)) return(NULL)
        ls <- np$log_ratio; off <- np$offset
      }
      edges[[length(edges) + 1L]] <- c(a, b)
      ls_est <- c(ls_est, ls)
      off_est <- c(off_est, off)  # alpha_b (mu_a - mu_b)
      off_n <- c(off_n, if (is.null(pair_fun)) length(ok_start) else np$n)
      link_converged <- c(link_converged,
                          if (is.null(np)) NA else np$converged)
      link_edge_mass <- c(link_edge_mass,
                          if (is.null(np)) NA else np$edge_mass)
      link_loglik <- c(link_loglik, if (is.null(np)) NA else np$loglik)
    }
    if (!length(edges)) {
      if (hard) stop("no set pairs share enough persons with informative ",
                     "score patterns to link the units: each person's ",
                     "pattern needs a score range of at least 4 within each ",
                     "set (at least four dichotomous items, or fewer ",
                     "polytomous ones), and concordant extreme scores do not ",
                     "count as support for semiparametric set linking")
      return(NULL)
    }
    comp <- utils::getFromNamespace(".efrm_components", "rasch")(S, edges)
    if (length(unique(comp)) > 1L) {
      if (hard)
        stop("item sets are not linked by common persons; relative units (alpha) ",
             "are unidentified between: ",
             paste(tapply(sets_u, comp, paste, collapse = "+"), collapse = " | "))
      return(NULL)
    }
    C <- matrix(0, length(edges), S)
    for (e in seq_along(edges)) { C[e, edges[[e]][2]] <- 1; C[e, edges[[e]][1]] <- -1 }
    sw <- sqrt(off_n); M <- C %*% A
    la <- drop(A %*% qr.coef(qr(M * sw), ls_est * sw)); la[is.na(la)] <- 0
    alpha <- exp(la)
    dmu <- off_est / alpha[vapply(edges, `[`, 1L, 2)]
    Cm <- matrix(0, length(edges), S)
    for (e in seq_along(edges)) { Cm[e, edges[[e]][1]] <- 1; Cm[e, edges[[e]][2]] <- -1 }
    mu <- drop(A %*% qr.coef(qr((Cm %*% A) * sw), dmu * sw)); mu[is.na(mu)] <- 0
    list(la = la, mu = mu, edges = edges, ls_est = ls_est, off_n = off_n,
         link_converged = link_converged, link_edge_mass = link_edge_mass,
         link_loglik = link_loglik)
  }

  N <- nrow(u_mat)
  point <- link_once(seq_len(N), hard = TRUE)
  # The corrected-moment link has no iterative edge optimiser and therefore
  # records NA, not FALSE, in link_converged. Only the optional NPML link has
  # a numerical convergence condition. Treating those structural NAs as
  # failures silently disabled every moment-link bootstrap.
  point_converged <- is.null(pair_link) ||
    all(point$link_converged %in% TRUE)

  # person bootstrap of the linking stage (skipped inside an outer
  # bootstrap). When a regen closure is supplied, each replicate also
  # redraws the within-frame thresholds and group units from their
  # estimated covariances and rebuilds the person estimates, so the
  # set-unit uncertainty carries the calibration noise that person
  # resampling alone cannot see: the set unit is a scale, and error in the
  # estimated threshold spread moves every person estimate's variance
  # coherently. Without this the hybrid log-alpha standard error
  # understates by ~20% in simulation and the unit tests reject at ~9-10%
  # instead of 5%.
  cov_link <- NULL
  cov_alpha_phi <- NULL
  cross_withheld <- FALSE
  link_reps <- dtilde_reps <- NULL
  if (boot_reps > 0 && point_converged) {
    # a pool of parameter draws cycled across the person resamples: each
    # replicate pairs one draw with one resample, so the replicate spread
    # still carries both variance components, at a fraction of the regen
    # cost of one draw per replicate
    n_draws <- if (is.null(regen)) 0L else {
      nd <- getOption("rasch.efrm_link_draws", max(50L, boot_reps %/% 5L))
      .check_whole(nd, "option rasch.efrm_link_draws", 30)
    }
    draws <- if (n_draws > 0L) vector("list", n_draws)
    draw_id <- integer(boot_reps)
    boot_idx <- vector("list", boot_reps)
    if (is.function(progress)) progress(0L, boot_reps)
    # Generate randomness in the same order for every worker count. Neither
    # link_once() nor the numerical optimiser draws random numbers.
    for (r in seq_len(boot_reps)) {
      if (is.function(cancel) && isTRUE(cancel())) .efrm_cancelled()
      if (n_draws > 0L) {
        di <- ((r - 1L) %% n_draws) + 1L
        if (is.null(draws[[di]])) draws[[di]] <- regen()
        draw_id[r] <- di
      }
      boot_idx[[r]] <- sample.int(N, N, replace = TRUE)
    }
    reps <- matrix(NA_real_, boot_reps, 2L * S)
    phi_reps <- NULL
    one_boot <- function(r) {
      mats <- if (n_draws == 0L) NULL else draws[[draw_id[r]]]
      b <- if (is.null(mats)) link_once(boot_idx[[r]])
           else link_once(boot_idx[[r]],
                          um = mats$u, wm = mats$w, gm = mats$g,
                          pair_fun = if (is.null(mats$pair_link))
                            pair_link else mats$pair_link)
      list(link = b,
           log_phi = if (is.null(mats)) NULL else mats$log_phi,
           dtilde = if (is.null(mats)) NULL else mats$dtilde)
    }
    ans <- .efrm_boot_apply(boot_reps, one_boot, workers,
                            progress = progress, cancel = cancel)
    for (r in seq_len(boot_reps)) {
      mats <- ans[[r]]
      b <- mats$link
      if (!is.null(b)) reps[r, ] <- c(b$la, b$mu)
      if (!is.null(mats$log_phi)) {
        if (is.null(phi_reps))
          phi_reps <- matrix(NA_real_, boot_reps, length(mats$log_phi),
                             dimnames = list(NULL, names(mats$log_phi)))
        phi_reps[r, ] <- mats$log_phi
      }
      if (!is.null(mats$dtilde)) {
        if (is.null(dtilde_reps))
          dtilde_reps <- matrix(NA_real_, boot_reps, length(mats$dtilde))
        dtilde_reps[r, ] <- mats$dtilde
      }
    }
    # Infinite log-units are failed link estimates, not usable bootstrap
    # draws. complete.cases() admits them and would let one row turn the
    # covariance into NaN.
    complete <- rowSums(is.finite(reps)) == ncol(reps)
    if (!is.null(dtilde_reps))
      complete <- complete &
        rowSums(is.finite(dtilde_reps)) == ncol(dtilde_reps)
    reps_ok <- reps[complete, , drop = FALSE]
    # log(alpha) and mu each obey a sum-zero identification constraint, so
    # their joint covariance has 2(S - 1) independent directions rather than
    # the 2S columns used to store it.
    min_success <- .efrm_min_boot_success(boot_reps, 2L * (S - 1L))
    if (nrow(reps_ok) < min_success)
      stop("too few unit-linking bootstrap replicates were usable for a ",
           "stable alpha covariance (", nrow(reps_ok), " of ", boot_reps,
           "; at least ", min_success, " are required). Raise `boot_reps`, ",
           "or strengthen the linking design")
    cov_link <- stats::cov(reps_ok)
    link_reps <- reps_ok
    if (!is.null(dtilde_reps))
      dtilde_reps <- dtilde_reps[complete, , drop = FALSE]
    if (!is.null(phi_reps)) {
      joint_ok <- complete &
        rowSums(is.finite(phi_reps)) == ncol(phi_reps)
      # The cross-covariance feeds frame-unit standard errors, so it meets the
      # same absolute/majority floor as the other bootstrap quantities.
      # Only the alpha-by-phi cross block is retained here; no joint
      # (S + G)-dimensional covariance is inverted.  Its sampling guard is
      # therefore the absolute/majority rule, while the alpha covariance above
      # is guarded against the 2(S - 1) free directions in its own block.
      joint_min <- .rasch_min_boot_success(boot_reps)
      if (sum(joint_ok) >= joint_min) {
        C <- stats::cov(cbind(reps[joint_ok, seq_len(S), drop = FALSE],
                              phi_reps[joint_ok, , drop = FALSE]))
        cov_alpha_phi <- C[seq_len(S), S + seq_len(ncol(phi_reps)),
                           drop = FALSE]
        dimnames(cov_alpha_phi) <- list(sets_u, colnames(phi_reps))
      } else {
        cross_withheld <- TRUE
        warning("only ", sum(joint_ok), " joint alpha-phi bootstrap draws ",
                "were usable (", joint_min, " required); the alpha-phi ",
                "cross-covariance is withheld and the affected frame-unit ",
                "standard errors are NA", call. = FALSE)
      }
    } else {
      # A person-only link can estimate alpha, but it contains no joint draw
      # with log(phi).  That missing cross-stage covariance is unknown, not
      # evidence that it is zero.
      cross_withheld <- TRUE
      warning("calibration redraws were unavailable, so the alpha-phi ",
              "cross-covariance is withheld and the affected frame-unit ",
              "standard errors are NA", call. = FALSE)
    }
  }

  list(alpha = setNames(exp(point$la), sets_u),
       se_log_alpha = if (is.null(cov_link)) setNames(rep(NA_real_, S), sets_u)
         else setNames(sqrt(pmax(diag(cov_link)[seq_len(S)], 0)), sets_u),
       mu = setNames(point$mu, sets_u),
       cov_link = cov_link,
       cov_alpha_phi = cov_alpha_phi,
       cross_cov_withheld = cross_withheld,
       link_reps = link_reps,
       dtilde_reps = dtilde_reps,
       boot_reps_requested = as.integer(boot_reps),
       boot_reps_used = if (is.null(link_reps)) 0L else nrow(link_reps),
       boot_reps_failed = as.integer(boot_reps) -
         (if (is.null(link_reps)) 0L else nrow(link_reps)),
       edges = data.frame(set_a = sets_u[vapply(point$edges, `[`, 1L, 1)],
                          set_b = sets_u[vapply(point$edges, `[`, 1L, 2)],
                          n = point$off_n, log_slope = point$ls_est,
                          converged = point$link_converged,
                          edge_mass = point$link_edge_mass,
                          loglik = point$link_loglik))
}

# Person estimation under unequal frame units: the weighted score
# W = sum_i rho_i x_i is sufficient, so roots are solved once per
# missing-data pattern and unique weighted score.
.efrm_person_estimates <- function(X, tau_list, disc) {
  N <- nrow(X)
  obs <- !is.na(X)
  m <- vapply(tau_list, length, 1L)
  pat <- apply(obs, 1, function(z) paste(which(z), collapse = ","))
  theta <- se <- rep(NA_real_, N)
  raw <- rowSums(X, na.rm = TRUE); raw[rowSums(obs) == 0L] <- NA
  max_raw <- as.numeric(obs %*% m)
  Xw <- sweep(X, 2, disc, "*")
  W <- rowSums(Xw, na.rm = TRUE); W[rowSums(obs) == 0L] <- NA
  # Extreme status belongs to the response pattern, not to the numerical
  # size of its sufficient statistic. With a small frame unit, a positive
  # response can have W below any fixed tolerance without being a zero score.
  away_from_zero <- obs & X != 0
  away_from_max <- obs & sweep(X, 2L, m, `!=`)
  extreme <- rowSums(obs) > 0L &
    (rowSums(away_from_zero, na.rm = TRUE) == 0L |
       rowSums(away_from_max, na.rm = TRUE) == 0L)

  for (key in unique(pat)) {
    cols <- as.integer(strsplit(key, ",", fixed = TRUE)[[1]])
    if (!length(cols)) next
    sel <- which(pat == key)
    r <- disc[cols]; tl <- tau_list[cols]
    curve <- .person_wle_curve(tl, r)
    # Work with relative units in the score equation and its moments;
    # probabilities retain the fitted units. This removes a common scale
    # factor without changing the root, and avoids powers of extreme units.
    unit_scale <- max(r)
    rp <- r / unit_scale
    pattern_score <- as.numeric(X[sel, cols, drop = FALSE] %*% rp)
    # Cache genuinely equal sufficient statistics only. Rounding can merge
    # different response patterns when frame units are highly unequal.
    for (Wu in unique(pattern_score)) {
      who <- sel[pattern_score == Wu]
      root <- .person_wle_maximum(curve, Wu)
      theta[who] <- root
      if (!is.na(root)) {
        V <- vapply(seq_along(cols), function(j)
          item_moments(root, tl[[j]], disc = r[j])$V, 0)
        se[who] <- (1 / sqrt(sum(rp^2 * V))) / unit_scale
      }
    }
  }
  data.frame(n_items = rowSums(obs), raw = raw, max_raw = max_raw,
             weighted_score = W, theta = theta, se = se, extreme = extreme)
}

# Person support for group-unit inference.  The stage-one likelihood uses
# same-set pairs only, and a pair at its minimum or maximum total has just one
# feasible allocation.  Count precisely those within-set comparisons that
# can inform the group unit; counting all co-observed virtual cells would
# include both cross-set pairs and conditional-likelihood constants.
.efrm_group_support <- function(Xv, vmap, m_v, glevs) {
  do.call(rbind, lapply(glevs, function(g) {
    load <- numeric(nrow(Xv))
    for (s in unique(vmap$set[vmap$group == g])) {
      cc <- which(vmap$group == g & vmap$set == s)
      if (length(cc) < 2L) next
      ij <- utils::combn(cc, 2L)
      for (k in seq_len(ncol(ij))) {
        i <- ij[1L, k]; j <- ij[2L, k]
        xi <- Xv[, i]; xj <- Xv[, j]
        informative <- !is.na(xi) & !is.na(xj) &
          !(xi == 0 & xj == 0) &
          !(xi == m_v[i] & xj == m_v[j])
        informative[is.na(informative)] <- FALSE
        load <- load + informative
      }
    }
    ww <- load[load > 0]
    data.frame(group = g, n_persons = length(ww),
               effective_persons = if (length(ww))
                 sum(ww)^2 / sum(ww^2) else 0,
               stringsAsFactors = FALSE)
  }))
}

# Support for a linked set parameter is limited by the weakest edge on its
# route to the reference set.  With more than one route, retain the route
# whose weakest edge is strongest.  Looking only at incident edges both lends
# strong terminal-edge support across a weak upstream link and lets a weak,
# redundant edge suppress a set that has a well-supported alternative path.
.efrm_set_path_support <- function(sets, edges, weight) {
  S <- length(sets)
  if (S == 1L) return(setNames(Inf, sets))
  cap <- setNames(rep(0, S), sets)
  cap[sets[1L]] <- Inf
  done <- setNames(rep(FALSE, S), sets)
  repeat {
    candidates <- which(!done)
    if (!length(candidates)) break
    u <- candidates[which.max(cap[candidates])]
    if (cap[u] <= 0) break
    done[u] <- TRUE
    rows <- which(edges$set_a == sets[u] | edges$set_b == sets[u])
    for (e in rows) {
      v_name <- if (edges$set_a[e] == sets[u]) edges$set_b[e] else
        edges$set_a[e]
      v <- match(v_name, sets)
      if (is.na(v)) next
      candidate <- min(cap[u], weight[e])
      if (!done[v] && candidate > cap[v]) cap[v] <- candidate
    }
  }
  cap
}

.efrm_wald_zero <- function(est, Sigma, term) {
  if (length(est) < 2L) return(NULL)
  unavailable <- function() data.frame(
    term = term, df = NA_integer_, wald = NA_real_, p = NA_real_)
  if (is.null(Sigma) || !is.matrix(Sigma) ||
      nrow(Sigma) != length(est) || ncol(Sigma) != length(est) ||
      any(!is.finite(est)) || any(!is.finite(Sigma)) ||
      !.covariance_is_symmetric(Sigma))
    return(unavailable())
  ee <- eigen((Sigma + t(Sigma)) / 2, symmetric = TRUE)
  cut <- max(abs(ee$values)) * 1e-8
  if (!is.finite(cut) || cut == 0 || min(ee$values) < -cut)
    return(unavailable())
  use <- ee$values > cut
  if (!any(use)) return(unavailable())
  estimable <- ee$vectors[, use, drop = FALSE]
  omitted <- est - drop(estimable %*% crossprod(estimable, est))
  if (sqrt(sum(omitted^2)) >
      1e-7 * max(1, sqrt(sum(est^2)))) return(unavailable())
  Sinv <- estimable %*% (t(estimable) / ee$values[use])
  W <- drop(t(est) %*% Sinv %*% est)
  data.frame(term = term, df = sum(use), wald = W,
             p = stats::pchisq(W, sum(use), lower.tail = FALSE))
}

#' Fit the extended frame of reference model
#'
#' Fits Humphry's extended frame of reference model, in which the unit can
#' differ across item-set by person-group frames. For item \eqn{i} in set
#' \eqn{s} and person \eqn{n} in group \eqn{g},
#' \deqn{P(X_{ni}=x)=\frac{\exp\{\rho_{sg}[x\theta_n-
#'   \sum_{k=1}^{x}\delta_{ik}]\}}
#'   {\sum_{y=0}^{m_i}\exp\{\rho_{sg}[y\theta_n-
#'   \sum_{k=1}^{y}\delta_{ik}]\}},\qquad
#'   \rho_{sg}=\alpha_s\phi_g.}
#'
#' @details
#' The partial credit model holds within each frame in its natural unit.
#' \eqn{\phi_g} and \eqn{\alpha_s} are unit \emph{ratios} in the sense of
#' Humphry and Andrich (2008, eq. 15): each is the common reference unit over
#' the frame's own unit. The identification constraints set the geometric mean
#' of the group units and of the set units to one; no observed group or set is
#' the reference level. A value above one therefore denotes a finer natural
#' unit than the corresponding geometric-mean unit and steeper curves on the
#' common scale. Ratios between two observed levels are obtained directly, for
#' example as \eqn{\alpha_s/\alpha_t}.
#' Person-group ratios \eqn{\phi_g} are identified from common item
#' thresholds across groups. Item sets partition the items, so set ratios
#' \eqn{\alpha_s} are identified instead from persons observed in more than
#' one set. The
#' set-linking graph and the group-by-set frame graph must each connect to a
#' common scale. Direct overlap between every pair of item sets is not
#' required: sets can be linked through intermediate sets. Pairs without
#' enough informative common persons contribute no edge; the remaining
#' graph must still connect all sets. A bootstrap replicate is unusable if
#' any supported link fails numerically or does not converge, even when
#' other links still connect the sets.
#'
#' Set units use a semiparametric likelihood for persons observed in each
#' linked pair of sets. For sets \eqn{a} and \eqn{b}, it maximises
#' \deqn{\prod_n\int P(X_{na}\mid u)P(X_{nb}\mid ru+c)\,dF_{g(n)}(u),}
#' where the masses of each observed group's \eqn{F_g}, the scale ratio
#' \eqn{r} and the offset \eqn{c} are estimated jointly on a fixed grid. This
#' avoids prescribing a normal or common person distribution across groups.
#' A link whose scale or offset reaches the numerical search boundary is
#' refused rather than reported as a finite estimate. So is a link whose
#' grid truncates the person distribution: persons with a finite location
#' in at least one set may place no more than two per cent of their
#' posterior mass on the ends of the grid. Persons at the same extreme tail
#' in both sets remain in the likelihood but do not count towards this check
#' or the link's inferential support; their likelihood rises towards one end
#' of any finite grid and carries no information about the link. Opposing
#' extremes retain a finite compromise location and do contribute.
#' The conditional thresholds and group units are held fixed in this step;
#' only \eqn{r}, \eqn{c}, and the nuisance masses are estimated. The linked
#' parameters are then
#' \deqn{\delta_{ik}=\widetilde\delta_{ik}/\alpha_s+\mu_s,
#' \qquad \rho_{sg}=\alpha_s\phi_g.}
#' Score moments supply starting values and screen weak links. Response
#' patterns must span a score range of at least four within a set. Overlapping
#' item sets are not permitted. The public
#' convergence flag covers the conditional calibration, the set-link
#' transformation and its nonparametric nuisance masses;
#' \code{stage1_converged} records the conditional stage separately.
#'
#' The empirical Godambe covariance from the conditional stage requires at
#' least ten informative persons, at least eight effective persons, more
#' effective persons than fitted stage-one directions, and full rank in the
#' projected score covariance. If these conditions fail, point estimates can
#' be returned with \code{boot_reps = 0}, but unit uncertainty is withheld.
#' A hybrid or full-bootstrap fit is refused because its linking and unit
#' covariance would otherwise inherit an unsupported stage-one covariance.
#'
#' The hybrid covariance combines the pairwise Godambe covariance with a
#' person bootstrap for set linking. Each replicate jointly redraws the
#' within-frame thresholds and group units, then rebuilds the link. The joint
#' draws retain covariance among common-scale thresholds, set units and group
#' units. With \code{se_method = "bootstrap"}, the complete model is refitted
#' to each person resample. Refits that do not converge or have unidentified
#' group units are discarded and counted as failed replicates.
#'
#' The \code{efrm_vs_rasch} component records the within-frame composite
#' log-likelihood comparison between group-dependent and equal group units.
#' This difference is descriptive and contains no information about set units,
#' which are identified at the linking stage. The accompanying Wald omnibus
#' tests provide inference for the group- and set-unit families. Their
#' probabilities are Holm-adjusted as one omnibus family; the individual
#' unit contrasts form a second Holm-adjusted follow-up family. An unavailable
#' probability remains in its declared family. Unit estimates
#' are retained for sparse designs, but probabilities require at least 50
#' persons or effective persons in every group. Set-unit inference requires
#' at least 50 informative common persons on the strongest bottleneck path
#' from every set to the first set, which is used only as the support graph's
#' bookkeeping root. Thus a weak upstream link limits a terminal set, while a
#' weak redundant edge does not suppress a stronger route.
#' Group-unit and dependent set-unit probabilities are withheld when any
#' group unit has a reported standard error above 5 log units. The estimates
#' and covariance remain descriptive. This check uses the returned uncertainty
#' method, including the full bootstrap when available.
#'
#' The model assumes that an item retains its location and discrimination
#' across the frames in which it appears, apart from the frame unit.
#' \code{\link{frame_invariance}} examines this assumption by separate frame
#' calibrations. Misfit concentrated within one item set can also distort its
#' estimated unit; inspect item fit and targeting before interpreting unit
#' differences. \code{\link{drop_items}} and \code{\link{resolve_frames}}
#' provide refitted sensitivity analyses.
#'
#' The dichotomous model follows Humphry (2005) and Humphry and Andrich
#' (2008). The polytomous, multigroup and crossed-frame forms are extensions
#' implemented in this package. The discrete nonparametric margin follows the
#' Rasch estimation approach of Follmann (1988); its use for linked item-set
#' units is an extension implemented here.
#' @param data Persons-by-items data (matrix or data frame, like
#'   \code{\link{rasch}}), plus a person-group column.
#' @param item_sets A named list mapping set names to item-column names, or
#'   a named character vector mapping every analysed item exactly once to a
#'   set. The vector cannot name items outside the analysis. Items not
#'   mentioned form their own set \code{"(rest)"} when a list is given.
#' @param groups Name of the person-group column in \code{data}, or a vector
#'   with one entry per person.
#'   Several columns define crossed group cells. Their units are returned in
#'   \code{phi_table}; \code{phi_factorial} and
#'   \code{phi_factorial_tests} contain the GLS factorial decomposition and
#'   omnibus Wald tests. Raw probabilities are retained in \code{p}; decisions
#'   use \code{p_adj}, Holm-adjusted across the factorial terms. Structurally
#'   unidentified units are refused. Very
#'   imprecise but identified units are retained with a warning.
#' @param id Person identifier, either a column name or one value per row.
#'   EFRM data require one response row per person, so identifiers must be
#'   unique when supplied. Missing or blank identifiers are treated as
#'   different unknown persons.
#' @param factors,items,n_groups,na_codes As in \code{\link{rasch}}.
#' @param maxit,tol Outer iteration cap and convergence tolerance of the
#'   bilinear pairwise stage.
#' @param min_link_persons Minimum number of common persons required for a
#'   set pair to contribute to the unit linking.
#' @param se_method \code{"hybrid"} (sandwich + linking bootstrap + delta
#'   propagation; default) or \code{"bootstrap"} (full person
#'   bootstrap of all stages).
#' @param boot_reps Bootstrap replicates; defaults to 300 for the linking
#'   bootstrap and 200 for the full bootstrap. Use zero to omit set-link
#'   uncertainty; in a multi-set fit, common-unit item and threshold standard
#'   errors are then unavailable. Otherwise at least 30 replicates are
#'   required. A bootstrap covariance is reported only when more than half of
#'   the requested replicates are usable.
#'   Inference is returned only when at least 30 replicates succeed, a
#'   majority of those requested, and the requested count exceeds the number
#'   of independent directions in the largest covariance block used by the
#'   fit. The fit stops if the linking covariance cannot meet that rule; an
#'   unsuccessful full bootstrap falls back to hybrid standard errors with a
#'   warning and retains its replicate accounting.
#' @param progress Optional function called as \code{progress(stage, current,
#'   total)} during long uncertainty calculations. It is intended for
#'   interfaces and batch logging and does not alter estimation.
#' @param cancel Optional zero-argument function checked between bootstrap
#'   batches. Returning \code{TRUE} stops with a \code{rasch_cancelled}
#'   condition. A serial fit uses one replicate per batch.
#' @param workers Number of parallel bootstrap workers. The default is four,
#'   reduced when fewer physical cores are available or the R process has a
#'   lower system limit. Random samples are generated before distribution, so
#'   a fixed seed gives the same result for any worker count. Every worker
#'   holds its own copy of the bootstrap state.
#' @param seed Optional bootstrap seed. The caller's random-number state is
#'   restored when estimation finishes; see \code{\link{rasch_rng}} for
#'   generator support.
#' @return An object of classes \code{"rasch_efrm"} and \code{"rasch"}.
#'   Model-specific components include \code{frames}, \code{phi_table},
#'   \code{alpha_table}, \code{set_table}, common-unit item and threshold
#'   tables, group-specific \code{score_curves} (expected weighted sufficient
#'   score and conditional standard error by person location and exact
#'   observed-item pattern), \code{efrm_vs_rasch}, and
#'   \code{linking}, the person support used for unit inference in
#'   \code{unit_support}, and the active covariance blocks in
#'   \code{unit_cov}. For a full-bootstrap fit, all blocks in
#'   \code{unit_cov} are calculated from the same usable person resamples;
#'   otherwise they are the analytic within-frame and, when requested, hybrid
#'   linking covariances used by the reported tests. With several item sets
#'   and \code{boot_reps = 0}, \code{cov_delta} and the corresponding
#'   common-unit standard errors are unavailable because set-link uncertainty
#'   has not been estimated. The requested, usable and failed
#'   uncertainty
#'   replicates used by the returned uncertainty method are reported as
#'   \code{boot_reps_requested}, \code{boot_reps_used} and
#'   \code{boot_reps_failed}; the hybrid set-link counts are repeated inside
#'   \code{linking}. When a full bootstrap was requested, its requested,
#'   attempted, usable and failed counts are retained separately in the
#'   corresponding \code{full_boot_reps_*} components, including when the fit
#'   falls back to hybrid standard errors. See the extended frame of reference
#'   vignette for their interpretation. If the within-frame calibration does
#'   not converge, its covariance blocks, standard errors and all later
#'   inferential probabilities are withheld. Failure of only a set link does
#'   not invalidate the already
#'   converged within-frame calibration or group-unit estimates, but
#'   common-unit item, frame and person uncertainty is withheld because it
#'   depends on that link, including the standard errors in \code{score_curves}.
#' @references
#' Andrich, D. (1982). An extension of the Rasch model for ratings providing
#' both location and dispersion parameters. Psychometrika, 47(1), 105--113.
#'
#' Andrich, D. and Luo, G. (2003). Conditional pairwise estimation in the
#' Rasch model for ordered response categories using principal components.
#' Journal of Applied Measurement, 4(3), 205--221.
#'
#' Andrich, D. and Marais, I. (2019). A Course in Rasch Measurement Theory:
#' Measuring in the Educational, Social and Health Sciences. Springer.
#'
#' Follmann, D. (1988). Consistent estimation in the Rasch model based on
#' nonparametric margins. Psychometrika, 53, 553--562.
#' \doi{10.1007/BF02294407}
#'
#' Humphry, S. M. (2005). Maintaining a Common Arbitrary Unit in Social
#' Measurement. PhD thesis, Murdoch University.
#'
#' Humphry, S. M. (2010). Modeling the effects of person group factors on
#' discrimination. Educational and Psychological Measurement, 70(2),
#' 215--231.
#'
#' Humphry, S. M. (2012). Item set discrimination and the unit in the Rasch
#' model. Journal of Applied Measurement, 13(2), 165--180.
#'
#' Montuoro, P. and Humphry, S. M. (2024). Modeling the effect of reading
#' item clarity on item discrimination. Journal of Applied Measurement,
#' 24(3/4), 121--132.
#'
#' Humphry, S. M. and Andrich, D. (2008). Understanding the unit in the Rasch
#' model. Journal of Applied Measurement, 9(3), 249--264.
#' @seealso \code{\link{frame_invariance}}, which tests the item invariance
#'   this model assumes rather than imposing it, and \code{\link{drop_items}},
#'   which removes an item the test flags and refits. Also
#'   \code{\link{rasch}}, \code{\link{rasch_mfrm}},
#'   \code{\link{test_information}}, and \code{\link{simulate_efrm}}.
#' @examples
#' \donttest{
#' set.seed(1); Np <- 400
#' simP <- function(th, tau, r) { x <- 0:length(tau)
#'   p <- exp(r * (x * th - c(0, cumsum(tau)))); p / sum(p) }
#' grp <- rep(c("A", "B"), each = Np / 2)
#' phi <- c(A = 0.8, B = 1.25)
#' d <- seq(-1.5, 1.5, length.out = 10)
#' theta <- rnorm(Np)
#' X <- sapply(seq_along(d), function(i) sapply(seq_len(Np), function(n)
#'   sample(0:1, 1, prob = simP(theta[n], d[i], phi[grp[n]]))))
#' colnames(X) <- sprintf("I%02d", seq_along(d))
#' fit <- rasch_efrm(data.frame(X, grp = grp), item_sets = list(core = colnames(X)),
#'                   groups = "grp")
#' fit$phi_table
#' }
#' @export
rasch_efrm <- function(data, item_sets, groups, id = NULL, factors = NULL,
                       items = NULL, n_groups = NULL,
                       na_codes = -1, maxit = 50, tol = 1e-7,
                       min_link_persons = 30,
                       se_method = c("hybrid", "bootstrap"),
                       boot_reps = NULL, progress = NULL, cancel = NULL,
                       workers = 4L, seed = NULL) {
  .factors_sym <- substitute(factors)
  .factors_label <- if (is.name(.factors_sym))
    as.character(.factors_sym) else "factor"
  .check_column_names(data)
  .check_controls(maxit, tol)
  if (!is.null(id) && (!is.atomic(id) || !is.null(dim(id))))
    stop("`id` must name one data column or be a plain vector with one value per row",
         call. = FALSE)
  if (!is.null(items) &&
      (!(is.character(items) || is.numeric(items)) || is.complex(items) ||
       !is.null(dim(items)) || !is.null(oldClass(items)) || !length(items) ||
       anyNA(items)))
    stop("`items` must be a non-empty plain vector of item names or indices",
         call. = FALSE)
  if (!is.atomic(groups) || !length(groups) || !is.null(dim(groups)))
    stop("`groups` must name data columns or be a plain vector with one value per person",
         call. = FALSE)
  if (!is.null(factors) && !is.data.frame(factors) &&
      (!is.atomic(factors) || !is.null(dim(factors))))
    stop("`factors` must be a data frame, column names, or a plain vector with one value per person",
         call. = FALSE)
  min_link_persons <- .check_whole(min_link_persons, "min_link_persons", 1)
  if (!is.null(n_groups))
    n_groups <- .check_whole(n_groups, "n_groups", 2)
  n_groups_requested <- n_groups
  se_method <- match.arg(se_method)
  if (is.null(boot_reps)) boot_reps <- if (se_method == "hybrid") 300L else 200L
  boot_reps <- .check_whole(boot_reps, "boot_reps", 0)
  if (boot_reps > 0L && boot_reps < 30L)
    stop("EFRM uncertainty needs either zero or at least 30 bootstrap replicates")
  if (!is.null(progress) && !is.function(progress))
    stop("progress must be NULL or a function")
  if (!is.null(cancel) && !is.function(cancel))
    stop("cancel must be NULL or a function")
  workers <- .check_whole(workers, "workers", 1)
  workers <- min(workers, .efrm_available_workers(), max(1L, boot_reps))
  if (workers > 1L &&
      !file.exists(system.file("DESCRIPTION", package = "rasch")))
    stop("parallel EFRM workers require an installed package; install rasch ",
         "before using workers above one")
  if (!is.null(seed)) {
    seed <- .check_whole(seed, "seed", 0)
    old_seed <- .sim_seed_capture()
    on.exit(.sim_seed_restore(old_seed), add = TRUE)
    set.seed(seed)
  }
  report <- function(stage, current, total) {
    if (is.function(cancel) && isTRUE(cancel())) .efrm_cancelled()
    if (is.function(progress)) progress(stage, current, total)
    invisible(NULL)
  }
  report("conditional calibration", 0L, 1L)
  # Simulation-only grid override used to verify that reported links are not
  # artefacts of the 61-point default. It is deliberately not a public model
  # option: changing the grid is a validation exercise, not an analyst choice.
  link_grid_n <- getOption("rasch.efrm_link_grid_n", 61L)
  link_grid_n <- .check_whole(link_grid_n,
                              "option rasch.efrm_link_grid_n", 21)
  # --- roles ----------------------------------------------------------------
  id_vec <- NULL; fac_df <- NULL; grp <- NULL; grp_name <- "group"
  grp_components <- NULL
  if (is.data.frame(data)) {
    nm <- names(data)
    if (.role_columns(groups, nm, nrow(data))) {
      # column name(s); a character vector of length nrow(data) is the
      # group values themselves and is handled below
      miss <- setdiff(groups, nm)
      if (length(miss))
        stop("group column(s) not found in the data: ",
             paste(miss, collapse = ", "))
      # a repeated grouping adds a redundant crossed factor and destroys the
      # factorial decomposition of the frame units
      if (anyDuplicated(groups))
        stop("group column(s) named more than once: ",
             paste(unique(groups[duplicated(groups)]), collapse = ", "))
      if (length(groups) == 1L) {
        grp <- data[[groups]]; grp_name <- groups
      } else {
        # SEVERAL frame-defining factors: the frames are their crossed
        # cells, and a factorial decomposition of the cell units is
        # reported in phi_factorial
        grp_components <- data[groups]
        grp_components[] <- lapply(grp_components, .role_text_values)
        # each source grouping must be checked before the cells are crossed:
        # a blank component survives inside a crossed label ("g1: " is not
        # blank) and would be estimated as a frame of its own
        blank_c <- vapply(grp_components, function(v)
          any(!is.na(v) & !nzchar(v)), TRUE)
        if (any(blank_c))
          stop("blank value(s) in frame group column(s): ",
               paste(groups[blank_c], collapse = ", "),
               "; a whitespace-only label is not a group")
        grp <- .factor_cells(grp_components, sep = ":")
        grp_name <- paste(groups, collapse = ":")
      }
    }
    if (is.character(id) && length(id) == 1L) {
      if (!id %in% nm) stop("id column '", id, "' not found in the data")
      id_vec <- data[[id]]
    } else if (!is.null(id)) {
      if (length(id) != nrow(data))
        stop("`id` has ", length(id), " entries but the data has ",
             nrow(data), " rows")
      id_vec <- id
    }
    factors_are_cols <- .role_columns(factors, nm, nrow(data))
    if (factors_are_cols) {
      missf <- setdiff(factors, nm)
      if (length(missf))
        stop("factor column(s) not found in the data: ",
             paste(missf, collapse = ", "))
      if (anyDuplicated(factors))
        stop("factor column(s) named more than once: ",
             paste(unique(factors[duplicated(factors)]), collapse = ", "))
      fac_df <- data[, factors, drop = FALSE]
    } else if (is.data.frame(factors)) {
      if (nrow(factors) != nrow(data))
        stop("`factors` has ", nrow(factors), " rows but the data has ",
             nrow(data))
      if (anyDuplicated(names(factors)))
        stop("duplicate factor column name(s): ",
             paste(unique(names(factors)[duplicated(names(factors))]),
                   collapse = ", "))
      clash <- intersect(names(factors), nm)
      different <- clash[!vapply(clash, function(cn)
        .same_role_values(factors[[cn]], data[[cn]]), logical(1))]
      if (length(different) && is.null(items))
        stop("external factor column(s) share item-data names but contain ",
             "different values: ", paste(different, collapse = ", "),
             ". Rename the external factor column(s), or name the item ",
             "columns explicitly with items=", call. = FALSE)
      fac_df <- factors
    } else if (!is.null(factors) && is.atomic(factors)) {
      if (length(factors) != nrow(data))
        stop("`factors` must be column name(s) in the data, a data frame ",
             "with one row per person, or a vector with one entry per row")
      fac_df <- stats::setNames(data.frame(factors, stringsAsFactors = FALSE),
                                .factors_label)
    } else if (!is.null(factors))
      stop("`factors` must be column name(s) in the data, a data frame ",
           "with one row per person, or a vector with one entry per row")
    # Frame-defining columns are already stored below as part of the frame
    # structure. Repeating them as ordinary factors creates duplicate names
    # and can make a later DIF or refit select the wrong column.
    if (!is.null(fac_df) && .role_columns(groups, nm, nrow(data)))
      fac_df <- fac_df[, !names(fac_df) %in% groups, drop = FALSE]
    # a data column whose values are identical to a by-value role vector is
    # almost certainly that same variable: exclude it so it is not also
    # scored as a numeric item (the rule rasch() applies)
    val_of <- function(v) if (is.null(v)) NULL else
      nm[vapply(data, function(col)
        length(col) == length(v) &&
          .same_role_values(col, v), logical(1))]
    groups_are_cols <- .role_columns(groups, nm, nrow(data))
    # a data column whose values are identical to a by-value role vector
    # may be that same variable, or a genuine item that happens to agree.
    # Deciding silently risks fitting the wrong analysis either way, so an
    # ambiguous match is refused unless `items` states the item columns
    id_by_value <- !is.null(id) && !(is.character(id) && length(id) == 1L)
    val_matched <- c(if (id_by_value) val_of(id),
                     if (!is.null(factors) && !factors_are_cols &&
                         !is.data.frame(factors)) val_of(factors),
                     if (!groups_are_cols) val_of(groups))
    if (is.null(items) && length(val_matched))
      stop("data column(s) identical to a supplied role vector: ",
           paste(unique(val_matched), collapse = ", "),
           ". If they are the same variable, drop them from the data or ",
           "name the item columns with items=; a genuine item identical ",
           "to a role must be listed in items=")
    drop_cols <- c(if (is.character(id) && length(id) == 1L) id,
                   if (factors_are_cols) factors
                   # With explicit items=, an external factor data frame is
                   # separate from `data`, so a shared name is not a role
                   # collision.  Actual in-data role columns remain guarded.
                   else if (is.data.frame(factors) && is.null(items))
                     intersect(names(factors), nm),
                   if (groups_are_cols) groups)
    item_cols <- if (is.null(items)) setdiff(nm, drop_cols)
    else if (is.character(items)) {
      miss <- setdiff(items, nm)
      if (length(miss))
        stop("item column(s) not found in the data: ",
             paste(miss, collapse = ", "))
      items
    } else {
      if (!is.numeric(items) || any(!is.finite(items)) ||
          any(items != floor(items)) || any(items < 1) ||
          any(items > length(nm)))
        stop("numeric `items` indices must be whole numbers between 1 and ",
             length(nm))
      nm[as.integer(items)]
    }
    if (anyDuplicated(item_cols))
      stop("item column(s) named more than once: ",
           paste(unique(item_cols[duplicated(item_cols)]), collapse = ", "))
    clash <- intersect(item_cols, drop_cols)
    if (length(clash))
      stop("`items` includes id, group, or factor column(s): ",
           paste(clash, collapse = ", "),
           " -- name only item columns")
    X <- as.matrix(data[, item_cols, drop = FALSE])
  } else {
    X <- as.matrix(data)
    if (is.character(id) && length(id) == 1L)
      stop("id column '", id, "' cannot be looked up in matrix input; ",
           "supply the id values as a vector")
    if (!is.null(id)) {
      if (length(id) != nrow(X))
        stop("`id` has ", length(id), " entries but the data has ",
             nrow(X), " rows")
      id_vec <- id
    }
    if (is.data.frame(factors)) {
      if (nrow(factors) != nrow(X))
        stop("`factors` has ", nrow(factors), " rows but the data has ",
             nrow(X))
      if (anyDuplicated(names(factors)))
        stop("duplicate factor column name(s): ",
             paste(unique(names(factors)[duplicated(names(factors))]),
                   collapse = ", "))
      fac_df <- factors
    } else if (!is.null(factors) && is.atomic(factors)) {
      if (length(factors) != nrow(X))
        stop("`factors` must be a data frame with one row per person or a ",
             "vector with one entry per row (matrix input has no columns ",
             "to look up)")
      fac_df <- stats::setNames(data.frame(factors, stringsAsFactors = FALSE),
                                .factors_label)
    } else if (!is.null(factors))
      stop("`factors` must be a data frame with one row per person or a ",
           "vector with one entry per row")
    if (!is.null(items)) {
      if (anyDuplicated(items))
        stop("item column(s) named more than once: ",
             paste(unique(items[duplicated(items)]), collapse = ", "))
      if (is.character(items)) {
        miss <- setdiff(items, colnames(X))
        if (length(miss))
          stop("item column(s) not found in the data: ",
               paste(miss, collapse = ", "))
        X <- X[, items, drop = FALSE]
      } else {
        if (!is.numeric(items) || any(!is.finite(items)) ||
            any(items != floor(items)) || any(items < 1) ||
            any(items > ncol(X)))
          stop("numeric `items` indices must be whole numbers between 1 and ",
               ncol(X))
        X <- X[, as.integer(items), drop = FALSE]
      }
    }
  }
  if (is.null(grp)) {
    if (length(groups) != nrow(X))
      stop("'groups' must name a column of data or give one value per person")
    grp <- groups
  }
  # a person with no frame group cannot be placed in any set: frame
  # expansion would leave their row entirely missing and the fit would
  # complete without recording the loss. A whitespace-only label is not a
  # group either -- it would become a frame of its own
  gchr <- .role_text_values(grp)
  bad_grp <- is.na(grp) | !nzchar(gchr)
  if (any(bad_grp))
    stop(sum(bad_grp), " person(s) have a missing or blank frame group; ",
         "their responses would be dropped from every set -- assign a group ",
         "or remove those rows")
  # Fit exactly the labels that were validated. Otherwise leading or trailing
  # whitespace can split one substantive group into several frame units.
  grp <- factor(gchr)
  if (nlevels(grp) < 1L) stop("no person groups found")
  if (is.null(id_vec)) id_vec <- seq_len(nrow(X))
  id_text <- .role_text_values(id_vec)
  present_id <- !is.na(id_text) & nzchar(id_text)
  if (anyDuplicated(id_text[present_id]))
    stop("rasch_efrm needs one response row per person; duplicate identifiers ",
         "cannot be treated as independent people by the set-link likelihood ",
         "or its person bootstrap")
  # The crossed-cell column is internal metadata. Keep its readable name
  # unless it would collide with a component or an ordinary person factor.
  .check_factor_frame(fac_df)
  taken_factor_names <- unique(c(names(grp_components), names(fac_df)))
  if (grp_name %in% taken_factor_names) {
    grp_name <- ".frame_group"
    while (grp_name %in% taken_factor_names)
      grp_name <- paste0(grp_name, ".")
  }

  prep <- .prepare_X(X, na_codes = na_codes); X <- prep$X
  notes <- prep$notes
  m_item <- apply(X, 2, max, na.rm = TRUE); L <- ncol(X)

  # --- item sets --------------------------------------------------------------
  if (is.list(item_sets)) {
    if (is.data.frame(item_sets) || !is.null(dim(item_sets)) ||
        !all(vapply(item_sets, function(s)
          (is.character(s) || is.factor(s)) && is.null(dim(s)), logical(1))))
      stop("item_sets must be a named list of plain item-name vectors, or a named item-to-set vector",
           call. = FALSE)
    if (is.null(names(item_sets)) || anyNA(names(item_sets)) ||
        any(!nzchar(trimws(names(item_sets)))))
      stop("item_sets must be a NAMED list (each element a set of item ",
           "names); a blank name is not a set")
    names(item_sets) <- trimws(names(item_sets))
    if (anyDuplicated(names(item_sets)))
      stop("duplicate set name(s) in item_sets after trimming: ",
           paste(unique(names(item_sets)[duplicated(names(item_sets))]),
                 collapse = ", "))
    # an empty set is a frame the design cannot carry: fitting without it
    # answers a different question from the one that was asked
    empty <- vapply(item_sets, function(s)
      !length(s) || all(is.na(s)) || all(!nzchar(trimws(as.character(s)))),
      TRUE)
    if (any(empty))
      stop("item set(s) with no items: ",
           paste(names(item_sets)[empty], collapse = ", "),
           "; every set needs at least one item name")
    ov <- unlist(item_sets)
    if (anyDuplicated(ov))
      stop("item(s) assigned to more than one set: ",
           paste(unique(ov[duplicated(ov)]), collapse = ", "))
    unknown <- setdiff(ov, colnames(X))
    if (length(unknown))
      stop("item_sets name item(s) not in the data: ",
           paste(unknown, collapse = ", "))
    set_of <- rep(NA_character_, L); names(set_of) <- colnames(X)
    for (s in names(item_sets)) {
      hit <- intersect(item_sets[[s]], colnames(X))
      set_of[hit] <- s
    }
    if (anyNA(set_of)) {
      # the generated name for unlisted items must not be one the caller
      # already used, or two different sets would become one
      if ("(rest)" %in% names(item_sets))
        stop("item_sets already contains a set named '(rest)', which is the ",
             "generated name for items no set lists; rename it, or list ",
             "every item explicitly")
      set_of[is.na(set_of)] <- "(rest)"
    }
  } else {
    if (!(is.character(item_sets) || is.factor(item_sets)) ||
        !is.null(dim(item_sets)))
      stop("item_sets must be a named list of plain item-name vectors, or a named item-to-set vector",
           call. = FALSE)
    if (is.null(names(item_sets)) || anyNA(names(item_sets)) ||
        any(!nzchar(trimws(names(item_sets)))))
      stop("item_sets must be a named list or named vector; a blank item ",
           "name is not an item")
    sv <- .role_text_values(item_sets)
    if (anyNA(item_sets) || any(!nzchar(sv)))
      stop("item_sets maps item(s) to a blank set name: ",
           paste(names(item_sets)[is.na(item_sets) | !nzchar(sv)],
                 collapse = ", "))
    if (anyDuplicated(names(item_sets)))
      stop("duplicate item(s) in the item_sets map: ",
           paste(unique(names(item_sets)[duplicated(names(item_sets))]),
                 collapse = ", "),
           "; each item may be assigned to one set")
    extra <- setdiff(names(item_sets), colnames(X))
    if (length(extra))
      stop("item_sets map item(s) not in the data: ",
           paste(extra, collapse = ", "))
    set_of <- sv[match(colnames(X), names(item_sets))]
    if (anyNA(set_of))
      stop("item(s) missing from the item_sets map: ",
           paste(colnames(X)[is.na(set_of)], collapse = ", "))
    names(set_of) <- colnames(X)
  }
  sets_u <- sort(unique(set_of)); S <- length(sets_u)
  # STRING-SORTED group levels, matching .efrm_solve's internal ordering
  # exactly: interaction() orders levels first-factor-fastest, the solver
  # sorts labels as strings, and attaching the solver's phi positionally
  # to a differently ordered glevs swapped crossed cells (g2:N took
  # g1:S's unit in a 2-by-2, manufacturing a spurious factor effect)
  glevs <- sort(levels(grp)); G <- length(glevs)

  # --- virtual columns: item x group ------------------------------------------
  vmap <- expand.grid(item = colnames(X), group = glevs,
                      KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  vmap$set <- set_of[vmap$item]
  vmap$vkey <- as.character(.factor_cells(vmap[c("item", "group")],
                                           sep = ":"))
  Xv <- matrix(NA_integer_, nrow(X), nrow(vmap),
               dimnames = list(NULL, vmap$vkey))
  for (v in seq_len(nrow(vmap))) {
    sel <- which(as.character(grp) == vmap$group[v])
    Xv[sel, v] <- X[sel, vmap$item[v]]
  }
  keep <- colSums(!is.na(Xv)) > 0L
  Xv <- Xv[, keep, drop = FALSE]; vmap <- vmap[keep, , drop = FALSE]
  rownames(vmap) <- NULL
  m_v <- m_item[vmap$item]
  thr_v <- threshold_index(m_v)

  # delta enumeration over underlying items (ordered by set, then item)
  iord <- order(match(set_of, sets_u), colnames(X))
  items_o <- colnames(X)[iord]
  thr_items <- threshold_index(m_item[items_o])
  Md <- nrow(thr_items)
  # Bootstrap rank is determined by the free directions after the set, group
  # and origin constraints, not by the length of the concatenated output.  The
  # largest full-bootstrap block is either the common-unit threshold block or
  # the joint (dtilde, log phi) block exposed in unit_cov.
  link_cov_rank <- 2L * max(S - 1L, 0L)
  stage1_joint_rank <- max(Md - S, 0L) + max(G - 1L, 0L)
  full_cov_rank <- max(Md - 1L, stage1_joint_rank, S - 1L, 0L)
  bootstrap_cov_dim <- if (se_method == "bootstrap")
    max(link_cov_rank, full_cov_rank) else link_cov_rank
  if (boot_reps > 0L && boot_reps <= bootstrap_cov_dim)
    stop("EFRM uncertainty for this design needs at least ",
         bootstrap_cov_dim + 1L, " replicates to span the free directions in ",
         "its largest covariance block; increase `boot_reps`, or use zero to ",
         "omit unit uncertainty")
  # map each virtual threshold row to its delta row
  drow <- vapply(seq_len(nrow(thr_v)), function(r) {
    it <- vmap$item[thr_v$item[r]]
    which(thr_items$item == match(it, items_o) & thr_items$k == thr_v$k[r])
  }, 1L)
  # per-set sum-zero blocks for dtilde
  set_of_drow <- set_of[items_o][thr_items$item]
  A_D <- matrix(0, Md, Md - S)
  cursor <- 0L
  for (s in sets_u) {
    rows <- which(set_of_drow == s); Ms <- length(rows)
    if (Ms < 2L) stop("set '", s, "' needs at least two thresholds")
    A_D[rows, cursor + seq_len(Ms - 1L)] <- rbind(diag(Ms - 1L), rep(-1, Ms - 1L))
    cursor <- cursor + Ms - 1L
  }

  # --- pairwise stage ----------------------------------------------------------
  pairs <- .efrm_filter_pairs(.pair_counts(Xv, m_v), vmap)
  if (!length(pairs)) stop("no informative within-frame item pairs")

  # phi-link check: two groups are joined only when they SHARE at least two
  # observed items of a common set. Sharing a set label alone is not enough:
  # the unit ratio phi_g/phi_h is identified through the common structural
  # thresholds of items both groups answered, and a spread comparison needs
  # at least two of them (disjoint item subsets of the same set leave the
  # ratio unidentified even though the old set-level check passed).
  has_data <- colSums(!is.na(Xv)) > 0L
  edges_g <- list()
  for (s in unique(vmap$set)) {
    by_grp <- lapply(glevs, function(g) {
      sel <- vmap$set == s & vmap$group == g & has_data
      unique(vmap$item[sel])
    })
    for (g1 in seq_along(glevs)) for (g2 in seq_len(g1 - 1L)) {
      if (length(intersect(by_grp[[g1]], by_grp[[g2]])) >= 2L)
        edges_g[[length(edges_g) + 1L]] <- c(g1, g2)
    }
  }
  comp_g <- .efrm_components(G, edges_g)
  if (length(unique(comp_g)) > 1L)
    stop("person groups are not linked by any common item set; relative units ",
         "(phi) are unidentified between: ",
         paste(tapply(glevs, comp_g, paste, collapse = "+"), collapse = " | "))

  sol <- .efrm_solve(Xv, thr_v, m_v, vmap, pairs, drow, A_D,
                     maxit = maxit, tol = tol)
  if (!isTRUE(sol$converged))
    warning("EFRM estimation did NOT converge in ", sol$iterations,
            " iterations; increase maxit or inspect the frame design",
            call. = FALSE)
  if (!isTRUE(sol$cluster_inference)) {
    notes <- c(notes, paste("EFRM conditional calibration:",
                            sol$cluster_note))
    if (boot_reps > 0L)
      stop(paste(
        "EFRM unit uncertainty is unavailable because the conditional",
        "calibration has insufficient independent-person support for its",
        "stage-one covariance; use boot_reps = 0 for a descriptive fit, or",
        "collect more informative persons"), call. = FALSE)
  }
  # STRUCTURAL non-identification (a flat direction of the information
  # along a unit) is an error: every common-unit quantity would silently
  # depend on the arbitrary point the optimiser stopped at. PRACTICAL
  # weakness -- a reported SE above 5 log-units, i.e. the unit uncertain
  # beyond a factor of exp(5) ~ 150 -- keeps its estimate for sensitivity
  # work but warns loudly and is noted: the SE already says the data
  # carry essentially no unit information. That check follows the choice
  # between analytic and full-bootstrap uncertainty below.
  bad_g <- sol$phi_unident |
    (isTRUE(sol$cluster_inference) & !is.finite(sol$se_log_phi))
  if (length(glevs) > 1L && any(bad_g))
    stop("the unit(s) of group(s) ", paste(glevs[bad_g], collapse = ", "),
         " are unidentified: the thresholds in their frames carry no ",
         "usable spread for phi to scale (information rank/conditioning ",
         "failure) -- refit without these groups, or include items whose ",
         "difficulties differ within the sets they answer")
  phi <- sol$phi; dtil <- sol$dtilde

  # Equal-unit comparison on the same conditional information.  It has its
  # own optimiser and therefore its own convergence/rank decision; a finite
  # objective left at an unfinished iterate is not a maximised likelihood.
  equal_fit <- .efrm_equal_fit(dtil, A_D, drow, thr_v, pairs, m_v,
                               maxit = 25L, tol = tol)
  if (!isTRUE(equal_fit$converged))
    notes <- c(notes, paste(
      "the equal-group-unit conditional refit did not converge or was",
      "rank-deficient, so its descriptive likelihood difference is unavailable"))

  # --- person-side linking (alpha, mu) ----------------------------------------
  # one builder for every linking path (point estimate, calibration-redraw
  # pool, outer bootstrap): person estimates u plus the correction moments
  # w = phi_w(R) and g = phi_g(R) per set (.person_link_moments), kept
  # coherent so every bootstrap replicate re-runs the corrected computation
  person_mats <- function(Xm, dt, phiv) {
    u <- w <- g <- matrix(NA_real_, nrow(Xm), S,
                          dimnames = list(NULL, sets_u))
    for (si in seq_len(S)) for (gl in glevs) {
      cols <- which(vmap$set == sets_u[si] & vmap$group == gl)
      if (!length(cols)) next
      # thresholds in dtilde units for these virtual columns
      tl <- lapply(cols, function(v) dt[drow[which(thr_v$item == v)]])
      Xs <- Xm[, cols, drop = FALSE]
      pe <- utils::getFromNamespace(".person_estimates", "rasch")(
        Xs, tl, disc = phiv[gl])
      sel <- which(pe$n_items > 0 & !pe$extreme)
      u[sel, si] <- pe$theta[sel]
      lm <- utils::getFromNamespace(".person_link_moments", "rasch")(
        Xs, tl, disc = phiv[gl])
      w[sel, si] <- lm$w[sel]; g[sel, si] <- lm$g[sel]
    }
    list(u = u, w = w, g = g)
  }
  make_pair_link <- function(Xm, dt, phiv) {
    tau_v_link <- lapply(seq_len(ncol(Xm)), function(v)
      dt[drow[which(thr_v$item == v)]])
    disc_v_link <- unname(phiv[vmap$group])
    function(a, b, idx, init_ls, init_off)
      utils::getFromNamespace(".efrm_npml_pair", "rasch")(
        Xm, vmap, tau_v_link, disc_v_link, sets_u,
        a, b, idx, init_ls, init_off, min_link_persons,
        grid_n = link_grid_n, report_support = TRUE)
  }
  report("conditional calibration", 1L, 1L)
  if (S > 1L) {
    pm <- person_mats(Xv, dtil, phi)
    pair_link <- make_pair_link(Xv, dtil, phi)
    # calibration-noise propagation into the linking bootstrap: a
    # symmetric square root of each stage-1 covariance, used to redraw
    # (dtilde, log phi) per bootstrap replicate under the documented
    # item-side/person-side independence treatment
    mat_sqrt <- function(V) {
      if (!.covariance_supports_wald(V)) return(NULL)
      ee <- eigen((V + t(V)) / 2, symmetric = TRUE)
      ee$vectors %*% (t(ee$vectors) * sqrt(pmax(ee$values, 0)))
    }
    K_dt <- length(dtil)
    # simulation-only escape hatches for the validation studies:
    # rasch.efrm_link_blockdiag = TRUE reproduces the marginal
    # (cross-covariance-free) draw for paired comparison with the joint
    # draw; rasch.efrm_link_draws overrides the parameter-draw pool size
    cj <- sol$cov_joint
    if (isTRUE(getOption("rasch.efrm_link_blockdiag", FALSE)) &&
        !is.null(cj)) {
      cj[seq_len(K_dt), K_dt + seq_along(phi)] <- 0
      cj[K_dt + seq_along(phi), seq_len(K_dt)] <- 0
    }
    L_j <- mat_sqrt(cj)
    if (boot_reps > 0L && is.null(L_j))
      stop(paste(
        "hybrid set-unit inference requires a finite, symmetric",
        "positive-semidefinite joint stage-one covariance; use boot_reps = 0",
        "for a descriptive fit"), call. = FALSE)
    regen <- if (!is.null(L_j)) function() {
      v <- drop(L_j %*% stats::rnorm(ncol(L_j)))
      phi_draw <- phi * exp(v[K_dt + seq_along(phi)])
      dt_draw <- dtil + v[seq_len(K_dt)]
      c(person_mats(Xv, dt_draw, phi_draw),
        list(log_phi = log(phi_draw),
             dtilde = dt_draw,
             pair_link = make_pair_link(Xv, dt_draw, phi_draw)))
    } else NULL
    link <- .efrm_link_sets(pm$u, pm$w, pm$g, sets_u, min_link_persons,
                            boot_reps = boot_reps, regen = regen,
                            pair_link = pair_link,
                            progress = function(current, total)
                              report("linking bootstrap", current, total),
                            cancel = cancel, workers = workers)
    alpha <- link$alpha; mu <- link$mu
    if (!all(link$edges$converged %in% TRUE)) {
      warning("one or more semiparametric set links stopped before the scale ",
              "transformation and nuisance masses met their convergence ",
              "tolerances; inspect fit$linking$alpha_edges",
              call. = FALSE)
      notes <- c(notes, paste(
        "one or more semiparametric set links stopped before the scale",
        "transformation and nuisance masses met their convergence tolerances"))
    }
  } else {
    alpha <- setNames(1, sets_u); mu <- setNames(0, sets_u)
    link <- list(alpha = alpha, se_log_alpha = setNames(0, sets_u), mu = mu,
                 cov_link = NULL, cov_alpha_phi = NULL,
                 cross_cov_withheld = FALSE,
                 boot_reps_requested = 0L, boot_reps_used = 0L,
                 boot_reps_failed = 0L, edges = data.frame())
  }

  # --- optional full person bootstrap of all stages ----------------------------
  boot <- NULL
  full_boot_reps_requested <- if (se_method == "bootstrap") boot_reps else 0L
  full_boot_reps_attempted <- 0L
  full_boot_reps_used <- 0L
  full_boot_reps_failed <- 0L
  if (se_method == "bootstrap" && boot_reps > 0L &&
      all(link$edges$converged %in% TRUE)) {
    full_boot_reps_attempted <- boot_reps
    Npers <- nrow(Xv)
    boot_replicate <- function(idx) {
      Xb <- Xv[idx, , drop = FALSE]
      pb <- utils::getFromNamespace(".efrm_filter_pairs", "rasch")(
        utils::getFromNamespace(".pair_counts", "rasch")(Xb, m_v), vmap)
      if (!length(pb)) return(NULL)
      sb <- utils::getFromNamespace(".efrm_solve", "rasch")(
        Xb, thr_v, m_v, vmap, pb, drow, A_D,
        maxit = maxit, tol = tol)
      # A small gradient does not establish identification: the ridged
      # solver can converge in a flat group-unit direction. Apply the same
      # structural refusal as the observed fit before linking or collecting
      # this draw. Its analytic SEs are not needed for bootstrap covariance.
      if (!isTRUE(sb$converged) || any(sb$phi_unident)) return(NULL)
      if (S > 1L) {
        pm_b <- person_mats(Xb, sb$dtilde, sb$phi)
        lb <- utils::getFromNamespace(".efrm_link_sets", "rasch")(
          pm_b$u, pm_b$w, pm_b$g, sets_u,
          min_link_persons, boot_reps = 0,
          pair_link = make_pair_link(Xb, sb$dtilde, sb$phi))
        if (!all(lb$edges$converged %in% TRUE)) return(NULL)
        ab <- lb$alpha; mb <- lb$mu
      } else { ab <- setNames(1, sets_u); mb <- setNames(0, sets_u) }
      db <- sb$dtilde / ab[set_of_drow] + mb[set_of_drow]
      # Retain both the stage-one thresholds and their common-unit
      # transformation.  The former are needed to return a joint
      # (dtilde, log phi) covariance from the same successful full-bootstrap
      # replicates as every other reported unit covariance.
      c(log(sb$phi), log(ab), mb, sb$dtilde, db)
    }
    collect <- matrix(NA_real_, boot_reps, G + 2L * S + 2L * Md)
    report("full person bootstrap", 0L, boot_reps)
    boot_idx <- lapply(seq_len(boot_reps), function(r)
      sample.int(Npers, Npers, replace = TRUE))
    ans <- .efrm_boot_apply(boot_reps, function(r)
      tryCatch(boot_replicate(boot_idx[[r]]), error = function(e) NULL),
      workers = workers,
      progress = function(current, total)
        report("full person bootstrap", current, total),
      cancel = cancel)
    for (r in seq_len(boot_reps)) {
      res <- ans[[r]]
      if (!is.null(res)) collect[r, ] <- res
    }
    collect <- collect[
      rowSums(is.finite(collect)) == ncol(collect), , drop = FALSE]
    full_boot_reps_used <- nrow(collect)
    full_boot_reps_failed <- boot_reps - full_boot_reps_used
    min_success <- .rasch_min_boot_success(boot_reps, full_cov_rank)
    if (full_boot_reps_used < min_success) {
      fallback_note <- sprintf(
        paste0("full person bootstrap: %d of %d replicates were usable; ",
               "at least %d are required, so hybrid standard errors were returned"),
        full_boot_reps_used, boot_reps, min_success)
      warning(fallback_note, call. = FALSE)
      notes <- c(notes, fallback_note)
    } else boot <- collect
  } else if (se_method == "bootstrap" && boot_reps > 0L) {
    notes <- c(notes, paste(
      "full person bootstrap was not attempted because the fitted set link",
      "did not meet its convergence criterion; hybrid standard errors were returned"))
  }

  # --- assembly in arbitrary units ----------------------------------------------
  delta <- dtil / alpha[set_of_drow] + mu[set_of_drow]
  rho_v <- alpha[vmap$set] * phi[vmap$group]
  thr_v$tau <- delta[drow]

  # Covariance of delta = dtilde/alpha + mu. The full bootstrap supplies
  # delta directly. In the hybrid bootstrap, each calibration draw is used to
  # rebuild the set link; retaining that draw beside its linked alpha and mu
  # preserves their covariance. Adding the two variance components as if
  # independent would discard precisely this shared-calibration term.
  cov_delta <- sol$cov_dtilde / tcrossprod(alpha[set_of_drow])
  unit_cov_method <- "stage-one"
  if (!is.null(boot)) {
    cov_delta <- stats::cov(
      boot[, G + 2L * S + Md + seq_len(Md), drop = FALSE])
    unit_cov_method <- "bootstrap"
  } else if (S > 1L && !is.null(link$link_reps) &&
             !is.null(link$dtilde_reps)) {
    lr <- link$link_reps
    dr <- link$dtilde_reps
    delta_reps <- matrix(NA_real_, nrow(lr), Md)
    for (r in seq_len(nrow(lr))) {
      ar <- exp(lr[r, seq_len(S)])
      mr <- lr[r, S + seq_len(S)]
      delta_reps[r, ] <- dr[r, ] / ar[match(set_of_drow, sets_u)] +
        mr[match(set_of_drow, sets_u)]
    }
    cov_delta <- stats::cov(delta_reps)
    unit_cov_method <- "hybrid"
  } else if (S > 1L && !is.null(link$cov_link)) {
    Sidx <- match(set_of_drow, sets_u)
    Caa <- link$cov_link[seq_len(S), seq_len(S), drop = FALSE]
    Cam <- link$cov_link[seq_len(S), S + seq_len(S), drop = FALSE]
    Cmm <- link$cov_link[S + seq_len(S), S + seq_len(S), drop = FALSE]
    g1 <- -(dtil / alpha[set_of_drow])     # d delta / d log alpha
    T2 <- matrix(g1, Md, Md) * Cam[Sidx, Sidx, drop = FALSE]
    cov_delta <- cov_delta + tcrossprod(g1) * Caa[Sidx, Sidx, drop = FALSE] +
      T2 + t(T2) + Cmm[Sidx, Sidx, drop = FALSE]
    unit_cov_method <- "hybrid"
  } else if (S > 1L) {
    # The transformation to the common unit contains estimated set scales
    # and origins.  Treating them as fixed would report only the conditional
    # calibration component and materially understate uncertainty.
    cov_delta[,] <- NA_real_
    unit_cov_method <- "stage-one-only"
    notes <- c(notes, paste(
      "common-unit item and threshold standard errors are withheld because",
      "set-link uncertainty was omitted (boot_reps = 0)"))
  }
  cov_tau <- cov_delta[drow, drow, drop = FALSE]
  thr_v$se <- sqrt(pmax(diag(cov_tau), 0))
  thr_v$anchored <- FALSE
  # a virtual item threshold on a near-empty category is a boundary
  # artefact here too: flag it and report NA rather than a covariance-based
  # number, the same honesty rasch()/pcml()/rasch_mfrm() apply
  weak <- .pcml_weak_thresholds(Xv, m_v, thr_v, colnames(Xv))
  thr_v$weak <- weak$flag
  thr_v$se[weak$flag] <- NA_real_
  if (length(weak$notes)) notes <- c(notes, weak$notes)
  link_converged <- S == 1L || all(link$edges$converged %in% TRUE)
  est <- list(model = "EFRM", thr = thr_v, cov_tau = cov_tau,
              loglik = sol$loglik, iterations = sol$iterations,
              converged = isTRUE(sol$converged) && link_converged,
              stage1_converged = sol$converged,
              m = m_v, anchors = NULL,
              n_parameters = (Md - S) + (G - 1L) + 2L * (S - 1L),
              cluster_inference = sol$cluster_inference,
              cluster_support = sol$cluster_support,
              cluster_note = sol$cluster_note)

  fac_all <- data.frame(g = as.character(grp), stringsAsFactors = FALSE)
  names(fac_all) <- grp_name
  if (!is.null(grp_components)) {
    gc_chr <- as.data.frame(lapply(grp_components, as.character),
                            stringsAsFactors = FALSE, check.names = FALSE)
    names(gc_chr) <- names(grp_components)
    fac_all <- cbind(fac_all, gc_chr)
  }
  if (!is.null(fac_df)) fac_all <- cbind(fac_all, fac_df)
  fit <- .assemble_fit("EFRM", Xv, est, id_vec, fac_all, n_groups,
                       notes, disc = rho_v)
  # Several frame-specific copies of one item are calibration cells, not
  # additional administered items. Alpha and one raw-score conversion over
  # the expanded columns are therefore undefined. Retain both for the
  # one-frame reduction, where the fitted columns are the items and
  # the model is the ordinary Rasch model.
  expanded_cells <- any(table(vmap$item) > 1L)
  if (expanded_cells) {
    fit$alpha <- list(
      alpha = NA_real_, n = NA_integer_, applicable = FALSE,
      design_applicable = FALSE,
      reason = "not applicable when an item has several frame response cells")
    fit$score_table <- NULL
    fit$notes <- unique(c(fit$notes, paste(
      "a universal raw-score conversion is not defined across the expanded",
      "frame response cells; use score_curves and design-specific information")))
  } else fit$alpha$design_applicable <- TRUE
  # every frame-defining factor (the crossed cell and its components) is
  # frame structure, not a testable DIF factor
  fit$frame_group <- c(grp_name,
                       if (!is.null(grp_components)) names(grp_components))

  # --- structural tables -----------------------------------------------------------
  se_lp <- if (!is.null(boot)) apply(boot[, seq_len(G), drop = FALSE], 2, sd)
           else unname(sol$se_log_phi)
  se_la <- if (!is.null(boot)) apply(boot[, G + seq_len(S), drop = FALSE], 2, sd)
           else unname(link$se_log_alpha)
  fit$phi_table <- data.frame(group = glevs, phi = unname(phi),
                              se_log_phi = se_lp)
  # Use the uncertainty method actually returned. A successful full
  # bootstrap supersedes the analytic stage-one SE, so an intermediate
  # analytic weak-unit flag must not override its final SE.
  weak_g <- isTRUE(sol$converged) & is.finite(se_lp) & se_lp > 5
  if (G > 1L && any(weak_g)) {
    warning("the unit(s) of group(s) ",
            paste(glevs[weak_g], collapse = ", "), " are only weakly ",
            "identified (SE of log phi above 5): the estimates are kept ",
            "for sensitivity work, but the data carry essentially no ",
            "information about these units", call. = FALSE)
  }
  # Unit tests need support for the unit being tested. For phi, each person
  # contributes through within-frame item pairs; for alpha, support is the
  # number of informative common persons on the weakest edge of the strongest
  # path from that set to the reference set. Concordant extreme pairs remain
  # in the likelihood but contain no information about the transformation and
  # are not support for its Wald tests.
  # The estimates remain available below the boundary, but normal/Wald
  # probabilities are not reported. Sparse-null simulations showed material
  # size inflation at 10--30 persons and nominal behaviour at 100; 50 is the
  # prespecified minimum for inferential use.
  min_unit_persons <- 50L
  group_support <- .efrm_group_support(Xv, vmap, m_v, glevs)
  phi_support_ok <- all(group_support$n_persons >= min_unit_persons &
    group_support$effective_persons >=
      min_unit_persons - sqrt(.Machine$double.eps))
  set_support <- data.frame(set = sets_u, n_common_persons = Inf,
                            stringsAsFactors = FALSE)
  if (S > 1L)
    set_support$n_common_persons <- unname(.efrm_set_path_support(
      sets_u, link$edges, link$edges$n))
  alpha_support_ok <- link_converged &&
    (S == 1L || all(set_support$n_common_persons >= min_unit_persons))
  # A sufficient head count is not itself inferential availability. The
  # stage-one score covariance must support phi, and a multi-set alpha needs
  # the set-link covariance that boot_reps = 0 deliberately omits.
  phi_ok <- isTRUE(sol$converged) && !any(weak_g) &&
    phi_support_ok && isTRUE(sol$cluster_inference)
  alpha_cov_available <- S == 1L || !is.null(boot) ||
    !is.null(link$cov_link)
  # With one set alpha is fixed at one; otherwise its link uses phi.
  alpha_ok <- isTRUE(sol$converged) && (S == 1L || !any(weak_g)) &&
    alpha_support_ok && alpha_cov_available
  fit$unit_support <- list(group = group_support, set = set_support,
                           minimum_persons = min_unit_persons,
                           phi_inference = phi_ok,
                           alpha_inference = alpha_ok)
  if (!phi_ok) fit$notes <- unique(c(fit$notes,
    if (!isTRUE(sol$converged)) paste(
      "group-unit probabilities are withheld because the conditional",
      "calibration did not converge")
    else if (any(weak_g)) paste(
      "group-unit probabilities are withheld because at least one group",
      "unit has a standard error above 5 log units")
    else if (!isTRUE(sol$cluster_inference)) paste(
      "group-unit probabilities are withheld because the conditional",
      "calibration covariance lacks sufficient independent-person support")
    else paste0(
      "group-unit probabilities are withheld because at least one group has ",
      "fewer than 50 persons or effective persons contributing within-frame pairs")))
  if (!alpha_ok) fit$notes <- unique(c(fit$notes,
    if (!isTRUE(sol$converged)) paste(
      "set-unit probabilities are withheld because the conditional",
      "calibration did not converge")
    else if (any(weak_g)) paste(
      "set-unit probabilities are withheld because their links depend on",
      "a group unit with a standard error above 5 log units")
    else if (!link_converged)
      paste("set-unit probabilities and common-unit uncertainty are withheld",
            "because at least one fitted set link did not converge")
    else if (!alpha_cov_available)
      paste("set-unit probabilities are withheld because set-link uncertainty",
            "was omitted (boot_reps = 0)")
    else paste0("set-unit probabilities are withheld because at least one ",
                "set has fewer than 50 informative common persons on the ",
                "strongest path to the reference set")))
  # factorial decomposition of the cell units: generalised least squares
  # of log phi_cell on sum-coded main effects (and the interaction when
  # every cell is observed), using the JOINT covariance of the cell
  # log-units -- from the bootstrap replicates when available, else the
  # solver's centred analytic covariance. The centring constraint
  # (sum log phi = 0) makes the cells dependent and the covariance
  # singular along the constant, so a spectral pseudo-inverse is used.
  # Coefficients are reported descriptively (estimate, SE); inference is
  # carried by MULTI-DEGREE-OF-FREEDOM Wald tests per term in
  # phi_factorial_tests (coefficient-wise z would not be an omnibus test
  # for factors with more than two levels).
  if (!is.null(grp_components)) {
    # Recover component values from the recorded design rather than parsing
    # the composite label. Component levels may themselves contain colons.
    first <- match(glevs, as.character(grp))
    dd <- as.data.frame(lapply(grp_components, function(x)
      as.character(x[first])), stringsAsFactors = FALSE, check.names = FALSE)
    names(dd) <- names(grp_components)
    for (cn in names(dd)) dd[[cn]] <- factor(dd[[cn]])
    estimable <- all(vapply(dd, nlevels, 1L) >= 2L) &&
      nrow(dd) > sum(vapply(dd, nlevels, 1L) - 1L) &&
      !anyNA(fit$phi_table$phi)
    if (estimable) {
      original_names <- names(dd)
      safe_names <- sprintf("efrm_factor_%04d", seq_along(original_names))
      names(dd) <- safe_names
      ctr <- stats::setNames(rep(list("contr.sum"), ncol(dd)), safe_names)
      full <- nrow(dd) >= prod(vapply(dd, nlevels, 1L))
      fml <- stats::as.formula(paste("~", paste(safe_names, collapse =
        if (full && ncol(dd) > 1L) " * " else " + ")))
      Xf <- stats::model.matrix(fml, dd, contrasts.arg = ctr)
      display_names <- vapply(original_names, function(x) {
        if (grepl("[:`]", x))
          paste0("`", gsub("`", "``", x, fixed = TRUE), "`") else x
      }, "")
      relabel_factorial <- function(x) {
        vapply(strsplit(x, ":", fixed = TRUE), function(parts) {
          parts <- vapply(parts, function(part) {
            hit <- which(startsWith(part, safe_names))
            if (length(hit) != 1L) return(part)
            paste0(display_names[hit],
                   substring(part, nchar(safe_names[hit]) + 1L))
          }, "")
          paste(parts, collapse = ":")
        }, "")
      }
      # the centring constraint (sum log phi = 0) makes the intercept
      # inestimable -- and it lies in the null space of the centred
      # covariance, so keeping it would make the GLS cross-product
      # singular. The sum-coded effects describe the centred cells fully.
      asg <- attr(Xf, "assign")
      Xf <- Xf[, asg != 0L, drop = FALSE]
      asg <- asg[asg != 0L]
      Sig <- if (!is.null(boot))
        stats::cov(boot[, seq_len(G), drop = FALSE])
      else sol$cov_log_phi
      if (!.covariance_supports_wald(Sig, G)) {
        fit$notes <- unique(c(fit$notes, paste(
          "the crossed group-unit decomposition is unavailable because its",
          "covariance is asymmetric or not positive semidefinite")))
      } else {
        eg <- eigen((Sig + t(Sig)) / 2, symmetric = TRUE)
        pos <- eg$values > max(eg$values) * 1e-8
        Sinv <- eg$vectors[, pos, drop = FALSE] %*%
          (t(eg$vectors[, pos, drop = FALSE]) / eg$values[pos])
        XtS <- t(Xf) %*% Sinv
        V <- tryCatch(solve(XtS %*% Xf), error = function(e) NULL)
        if (!is.null(V) && .covariance_supports_wald(V, ncol(Xf))) {
          cf <- drop(V %*% (XtS %*% log(fit$phi_table$phi)))
          fit$phi_factorial <- data.frame(
            term = relabel_factorial(colnames(Xf)), log_unit = cf,
            se = sqrt(pmax(diag(V), 0)), stringsAsFactors = FALSE)
          tls <- attr(stats::terms(fml), "term.labels")
          tests <- list()
          for (tno in seq_along(tls)) {
            ii <- which(asg == tno)
            if (!length(ii)) next
            Vt <- V[ii, ii, drop = FALSE]
            W <- if (.covariance_is_psd(Vt))
              tryCatch(drop(t(cf[ii]) %*% solve(Vt) %*% cf[ii]),
                       error = function(e) NA_real_) else NA_real_
            if (is.finite(W) && W < 0) W <- NA_real_
            tests[[length(tests) + 1L]] <- data.frame(
              term = relabel_factorial(tls[tno]), df = length(ii), wald = W,
              p = stats::pchisq(W, length(ii), lower.tail = FALSE),
              stringsAsFactors = FALSE)
          }
          fit$phi_factorial_tests <- do.call(rbind, tests)
          fit$phi_factorial_tests$p_adj <- NA_real_
          usable <- is.finite(fit$phi_factorial_tests$p)
          fit$phi_factorial_tests$p_adj[usable] <- stats::p.adjust(
            fit$phi_factorial_tests$p[usable], method = "holm",
            n = nrow(fit$phi_factorial_tests))
          fit$phi_factorial_tests$significant <- ifelse(
            is.finite(fit$phi_factorial_tests$p_adj),
            fit$phi_factorial_tests$p_adj < 0.05, NA)
          if (!phi_ok) {
            fit$phi_factorial_tests$p <- NA_real_
            fit$phi_factorial_tests$p_adj <- NA_real_
            fit$phi_factorial_tests$significant <- NA
          }
        } else {
          fit$notes <- unique(c(fit$notes, paste(
            "the crossed group-unit decomposition is unavailable because its",
            "coefficient covariance could not be estimated reliably")))
        }
      }
    }
  }
  fit$alpha_table <- data.frame(set = sets_u, alpha = unname(alpha),
                                se_log_alpha = se_la)
  fit$set_table <- data.frame(set = sets_u, mu = unname(mu),
                              alpha = unname(alpha),
                              n_items = as.integer(table(set_of)[sets_u]))
  # a distinct threshold is weak if any virtual threshold folded into it
  # rests on a near-empty category; carry that to the common-unit tables
  weak_d <- logical(length(delta))
  agg_w <- tapply(weak$flag, drow, any)
  weak_d[as.integer(names(agg_w))] <- as.logical(agg_w)
  se_arb <- sqrt(pmax(diag(cov_delta), 0)); se_arb[weak_d] <- NA_real_
  fit$thresholds_arbitrary <- data.frame(item = items_o[thr_items$item],
                                         set = set_of_drow,
                                         k = thr_items$k, delta = delta,
                                         se = se_arb, weak = weak_d)
  fit$item_arbitrary <- do.call(rbind, lapply(seq_along(items_o), function(i) {
    rows <- which(thr_items$item == i)
    weak_i <- any(weak_d[rows])
    data.frame(item = items_o[i], set = set_of[items_o[i]],
               location = mean(delta[rows]),
               se = if (weak_i) NA_real_ else
                 sqrt(max(mean(cov_delta[rows, rows, drop = FALSE]), 0)),
               weak = weak_i)
  }))
  rownames(fit$item_arbitrary) <- NULL

  # frame table with pooled fit
  fr <- unique(vmap[, c("set", "group")])
  fit$frames <- do.call(rbind, lapply(seq_len(nrow(fr)), function(j) {
    cols <- which(vmap$set == fr$set[j] & vmap$group == fr$group[j])
    npers <- sum(rowSums(!is.na(Xv[, cols, drop = FALSE])) > 0)
    gf <- .group_col_fit(fit$residuals, fit$moments, cols, disc = rho_v,
                         extreme = fit$person$extreme,
                         f_cell = fit$summary_stats$df_factor)
    data.frame(set = fr$set[j], group = fr$group[j],
               n_persons = npers, n_items = length(cols),
               alpha = unname(alpha[fr$set[j]]), phi = unname(phi[fr$group[j]]),
               rho = unname(alpha[fr$set[j]] * phi[fr$group[j]]),
               # with bootstrap replicates the SE of log rho comes from
               # the JOINT draws (log phi + log alpha per replicate), so
               # cross-stage dependence is captured. The hybrid calculation
               # uses the matching covariance from its joint calibration
               # redraws rather than adding the two marginal variances.
               se_log_rho = if (!is.null(boot))
                 stats::sd(boot[, match(fr$group[j], glevs)] +
                           boot[, G + match(fr$set[j], sets_u)])
               else {
                 is <- match(fr$set[j], sets_u)
                 ig <- match(fr$group[j], glevs)
                 # a cross-covariance withheld for too few joint draws is
                 # unknown, not zero: adding the marginal variances as
                 # though it were zero is not necessarily conservative,
                 # so the standard error is NA rather than a guess
                 if (isTRUE(link$cross_cov_withheld)) NA_real_ else {
                   cross <- if (is.null(link$cov_alpha_phi)) 0 else
                     link$cov_alpha_phi[is, ig]
                   sqrt(pmax(fit$alpha_table$se_log_alpha[is]^2 +
                               fit$phi_table$se_log_phi[ig]^2 + 2 * cross, 0))
                 }
               },
               origin = unname(mu[fr$set[j]]),
               infit_ms = gf$infit_ms, outfit_ms = gf$outfit_ms,
               fit_resid = gf$fit_resid, n_responses = gf$n)
  }))
  rownames(fit$frames) <- NULL

  # equal-unit comparison and score curves. The pairwise comparison is only
  # informative about the group units (phi): the within-frame likelihood is
  # invariant to the set units (alpha), which are identified person-side, so
  # their evidence is the Wald test on log alpha.
  Sig_phi <- if (!is.null(boot))
    stats::cov(boot[, seq_len(G), drop = FALSE]) else sol$cov_log_phi
  Sig_alpha <- if (S > 1L) {
    if (!is.null(boot)) stats::cov(boot[, G + seq_len(S), drop = FALSE])
    else if (!is.null(link$cov_link))
      link$cov_link[seq_len(S), seq_len(S), drop = FALSE]
    else NULL
  } else NULL
  unit_omnibus <- do.call(rbind, Filter(Negate(is.null), list(
    if (G > 1L) .efrm_wald_zero(log(phi), Sig_phi, "group units (phi)"),
    if (S > 1L) .efrm_wald_zero(log(alpha), Sig_alpha, "set units (alpha)"))))
  if (!is.null(unit_omnibus)) {
    unit_omnibus$p[unit_omnibus$term == "group units (phi)" & !phi_ok] <- NA_real_
    unit_omnibus$p[unit_omnibus$term == "set units (alpha)" & !alpha_ok] <- NA_real_
    unit_omnibus$p_adj <- NA_real_
    usable <- is.finite(unit_omnibus$p)
    unit_omnibus$p_adj[usable] <- stats::p.adjust(
      unit_omnibus$p[usable], method = "holm", n = nrow(unit_omnibus))
    unit_omnibus$significant <- ifelse(
      is.finite(unit_omnibus$p_adj), unit_omnibus$p_adj < 0.05, NA)
  }

  ut <- rbind(
    if (G > 1L) data.frame(parameter = paste0("log phi[", glevs, "]"),
                           estimate = log(fit$phi_table$phi),
                           se = fit$phi_table$se_log_phi,
                           family = "phi"),
    if (S > 1L) data.frame(parameter = paste0("log alpha[", sets_u, "]"),
                           estimate = log(alpha),
                           se = fit$alpha_table$se_log_alpha,
                           family = "alpha"))
  if (!is.null(ut)) {
    ut$z <- .wald_ratio(ut$estimate, ut$se)
    ut$p <- 2 * pnorm(-abs(ut$z))
    ut$p[ut$family == "phi" & !phi_ok] <- NA_real_
    ut$p[ut$family == "alpha" & !alpha_ok] <- NA_real_
    ut$p_adj <- NA_real_
    usable <- is.finite(ut$p)
    ut$p_adj[usable] <- stats::p.adjust(
      ut$p[usable], method = "holm", n = nrow(ut))
    ut$significant <- ifelse(is.finite(ut$p_adj), ut$p_adj < 0.05, NA)
    ut$family <- NULL
    rownames(ut) <- NULL
  }
  # Expose one coherent covariance method.  In a full-bootstrap fit all
  # blocks below use the same retained person resamples; otherwise they are
  # the corresponding analytic/hybrid blocks used by the reported tests.
  if (!is.null(boot)) {
    i_phi <- seq_len(G)
    i_alpha <- G + seq_len(S)
    i_dtilde <- G + 2L * S + seq_len(Md)
    cov_dtilde_unit <- stats::cov(boot[, i_dtilde, drop = FALSE])
    cov_joint_unit <- stats::cov(cbind(
      boot[, i_dtilde, drop = FALSE], boot[, i_phi, drop = FALSE]))
    cov_alpha_phi_unit <- if (S > 1L)
      stats::cov(boot[, i_alpha, drop = FALSE],
                 boot[, i_phi, drop = FALSE]) else NULL
  } else {
    cov_dtilde_unit <- sol$cov_dtilde
    cov_joint_unit <- sol$cov_joint
    cov_alpha_phi_unit <- link$cov_alpha_phi
  }
  fit$unit_cov <- list(cov_dtilde = cov_dtilde_unit,
                       cov_log_phi = Sig_phi,
                       cov_log_alpha = Sig_alpha,
                       cov_joint = cov_joint_unit,
                       cov_log_alpha_phi = cov_alpha_phi_unit,
                       cov_delta = cov_delta,
                       method = unit_cov_method)
  ll_equal <- if (isTRUE(equal_fit$converged)) equal_fit$loglik else NA_real_
  fit$efrm_vs_rasch <- list(ll_efrm = sol$loglik, ll_equal = ll_equal,
                            two_delta_ll = if (is.finite(ll_equal))
                              2 * (sol$loglik - ll_equal) else NA_real_,
                            equal_converged = equal_fit$converged,
                            extra_parameters = G - 1L,
                            informative_for = if (G > 1L) "group units (phi)"
                              else "nothing: single group, set units are identified person-side",
                            unit_omnibus = unit_omnibus,
                            unit_tests = ut)
  grid <- .default_model_grid(fit, half_width = 6, by = 0.1)
  # A score curve belongs to an administration, not only to a group: two
  # people in one group who were given different items have different maximum
  # weighted scores and different expected totals.  The response matrix has
  # no separate availability indicator, so its exact non-missing item pattern
  # is the only observable administration pattern.  Do not enlarge that
  # pattern to the whole set: doing so adds unobserved items to the weighted
  # sufficient score and its information.
  obs <- !is.na(X)
  gch <- as.character(grp)
  # The fixed-width bit key is collision-free for the ordered item columns.
  # Human-readable labels are never parsed back into item membership, so item
  # and set names may themselves contain the display separators.
  key <- paste0(gch, "\r", apply(obs, 1L, function(v)
    paste0(as.integer(v), collapse = "")))
  reps <- which(!duplicated(key) & rowSums(obs) > 0L)
  design_label <- function(v) {
    parts <- vapply(sets_u, function(s) {
      ii <- which(set_of == s)
      jj <- ii[v[ii]]
      if (!length(jj)) return(NA_character_)
      if (length(jj) == length(ii)) return(s)
      paste0(s, " [", paste(colnames(X)[jj], collapse = ", "), "]")
    }, "")
    paste(parts[!is.na(parts)], collapse = " + ")
  }
  fit$score_curves <- do.call(rbind, lapply(reps, function(p) {
    g <- gch[p]
    cols <- which(obs[p, ])
    r_i <- alpha[set_of] * phi[g]
    ew <- vapply(grid, function(th) sum(vapply(cols, function(i)
      r_i[i] * item_moments(th, delta[thr_items$item == match(colnames(X)[i], items_o)],
                            disc = r_i[i])$E, 0)), 0)
    info <- vapply(grid, function(th) sum(vapply(cols, function(i)
      r_i[i]^2 * item_moments(th, delta[thr_items$item == match(colnames(X)[i], items_o)],
                              disc = r_i[i])$V, 0)), 0)
    data.frame(group = g, design = design_label(obs[p, ]),
               n_persons = sum(key == key[p]),
               theta = grid, expected_score = ew, sem = 1 / sqrt(info))
  }))
  fit$linking <- list(
    phi_edges = edges_g,
    alpha_edges = link$edges,
    boot_reps_requested = link$boot_reps_requested,
    boot_reps_used = link$boot_reps_used,
    boot_reps_failed = link$boot_reps_failed,
    alpha_method = if (S > 1L)
      "finite-grid semiparametric maximum likelihood" else "not applicable",
    alpha_grid = if (S > 1L)
      c(lower = -8, upper = 8, points = link_grid_n) else NULL)
  report("finalising", 1L, 1L)
  fit$se_method <- if (!is.null(boot)) "bootstrap" else "hybrid"
  fit$boot_reps_requested <- if (!is.null(boot)) boot_reps else
    link$boot_reps_requested
  fit$boot_reps_used <- if (!is.null(boot)) nrow(boot) else
    link$boot_reps_used
  fit$boot_reps_failed <- fit$boot_reps_requested - fit$boot_reps_used
  fit$full_boot_reps_requested <- full_boot_reps_requested
  fit$full_boot_reps_attempted <- full_boot_reps_attempted
  fit$full_boot_reps_used <- full_boot_reps_used
  fit$full_boot_reps_failed <- full_boot_reps_failed
  fit$workers <- workers
  fit$seed <- seed
  fit$virtual_map <- vmap
  fit$set_of <- set_of
  fit$refit_spec <- list(
    groups = if (!is.null(grp_components)) names(grp_components) else grp_name,
    factors = setdiff(names(fit$factors), fit$frame_group),
    n_groups = n_groups_requested, na_codes = na_codes,
    maxit = maxit, tol = tol, min_link_persons = min_link_persons,
    se_method = se_method, boot_reps = boot_reps, workers = workers,
    seed = seed)
  if (!link_converged) {
    # The within-frame dtilde/phi calibration is still valid, but alpha and
    # mu define the transformation to the common unit.  A covariance using
    # only dtilde while that transformation failed would look precise by
    # silently treating the uncertain link as fixed.
    fit$est$thr$se[] <- NA_real_
    fit$est$cov_tau[,] <- NA_real_
    fit$thresholds_arbitrary$se[] <- NA_real_
    fit$item_arbitrary$se[] <- NA_real_
    fit$frames$se_log_rho[] <- NA_real_
    fit$score_curves$sem[] <- NA_real_
    for (nm in c("cov_log_alpha", "cov_log_alpha_phi", "cov_delta"))
      if (!is.null(fit$unit_cov[[nm]])) fit$unit_cov[[nm]][] <- NA_real_
  }
  if (!isTRUE(sol$converged)) {
    # A failed stage-one calibration invalidates every later uncertainty
    # calculation. Set-link non-convergence is different: the stage-one
    # group-unit and response-fit results remain interpretable on their own
    # scale and are therefore handled separately by the link guards.
    fit$est$thr$se[] <- NA_real_
    fit$est$cov_tau[,] <- NA_real_
    for (nm in setdiff(names(fit$unit_cov), "method"))
      if (!is.null(fit$unit_cov[[nm]])) fit$unit_cov[[nm]][] <- NA_real_
    fit$score_curves$sem[] <- NA_real_
    fit$phi_table$se_log_phi[] <- NA_real_
    fit$alpha_table$se_log_alpha[] <- NA_real_
    fit$thresholds_arbitrary$se[] <- NA_real_
    fit$item_arbitrary$se[] <- NA_real_
    fit$frames$se_log_rho[] <- NA_real_
    if (!is.null(fit$phi_factorial)) fit$phi_factorial$se[] <- NA_real_
    if (!is.null(fit$phi_factorial_tests)) {
      for (nm in intersect(c("wald", "p", "p_adj"),
                           names(fit$phi_factorial_tests)))
        fit$phi_factorial_tests[[nm]][] <- NA_real_
      fit$phi_factorial_tests$significant[] <- NA
    }
    for (part in c("unit_omnibus", "unit_tests")) {
      tab <- fit$efrm_vs_rasch[[part]]
      if (is.null(tab)) next
      for (nm in intersect(c("se", "z", "wald", "p", "p_adj"), names(tab)))
        tab[[nm]][] <- NA_real_
      if ("significant" %in% names(tab)) tab$significant[] <- NA
      fit$efrm_vs_rasch[[part]] <- tab
    }
    fit$efrm_vs_rasch$ll_efrm <- NA_real_
    fit$efrm_vs_rasch$two_delta_ll <- NA_real_
  }
  fit <- .tag_tables(fit)
  class(fit) <- c("rasch_efrm", "rasch")
  fit
}

#' @export
print.rasch_efrm <- function(x, ...) {
  separation_quality <- x$separation_quality %||% x$power_of_fit %||%
    .separation_quality(x$psi$PSI)
  cat(sprintf("rasch extended frame of reference analysis: %d items in %d set(s) x %d group(s) = %d frames, %d persons\n",
              length(x$set_of), nrow(x$alpha_table), nrow(x$phi_table),
              nrow(x$frames), nrow(x$X)))
  cat(sprintf("Within-frame pairwise conditional ML: %s in %d iterations\n",
              if (x$est$converged) "converged" else "NOT converged",
              x$est$iterations))
  cat(sprintf("PSI %.3f, separation quality: %s\n", x$psi$PSI,
              separation_quality))
  cat("\nPerson group units (phi):\n")
  print(x$phi_table, digits = 3, row.names = FALSE)
  cat("\nItem set units (alpha) and locations:\n")
  print(.fmt_df(merge(x$alpha_table, x$set_table[, c("set", "mu", "n_items")],
                      by = "set", sort = FALSE)), row.names = FALSE)
  cat(sprintf("\nEqual-unit comparison: 2(ll_EFRM - ll_equal) = %.3f with %d extra unit parameter(s)\n",
              x$efrm_vs_rasch$two_delta_ll, x$efrm_vs_rasch$extra_parameters))
  cat("(composite likelihood: descriptive; informative for ",
      x$efrm_vs_rasch$informative_for, ")\n", sep = "")
  if (!is.null(x$efrm_vs_rasch$unit_omnibus)) {
    cat("Omnibus Wald tests of equal units (Holm-adjusted family):\n")
    print(.fmt_df(x$efrm_vs_rasch$unit_omnibus), row.names = FALSE)
  }
  if (!is.null(x$efrm_vs_rasch$unit_tests)) {
    cat("Holm-adjusted exploratory unit contrasts (H0: unit = 1):\n")
    print(.fmt_df(x$efrm_vs_rasch$unit_tests), row.names = FALSE)
  }
  if (length(x$notes)) cat(sprintf("\nNotes: %s\n", paste(x$notes, collapse = "; ")))
  invisible(x)
}

#' Plot frame units
#'
#' Caterpillar plot of the frame units \code{rho_sg = alpha_s phi_g} on the
#' log scale, grouped by item set and coloured by person group, with 95 per
#' cent error bars; frames with pooled fit residuals beyond the band are
#' highlighted. Error bars are omitted when either fitted unit family does
#' not meet its inferential support conditions.
#'
#' @param fit A fitted object from \code{\link{rasch_efrm}}.
#' @param band Pooled fit residual band beyond which a frame is highlighted.
#' @return Called for its plotting side effect; invisibly \code{NULL}.
#' @examples
#' \donttest{
#' # see ?rasch_efrm for a complete simulated example
#' }
#' @export
plot_frames <- function(fit, band = 2.5) {
  if (!inherits(fit, "rasch_efrm")) stop("plot_frames needs a rasch_efrm fit")
  .check_response_display_fit(fit, "frame-unit plots")
  .check_band(band)
  fr <- fit$frames[order(fit$frames$set, fit$frames$rho), ]
  n <- nrow(fr)
  lr <- log(fr$rho)
  # A finite sandwich standard error can remain useful descriptively below
  # the calibrated support boundary, but drawing it as a 95 per cent interval
  # would reinstate inference that the fit explicitly withholds. A frame unit
  # depends on every non-structural unit family in the fit, so show intervals
  # only when those families and the row itself support them.
  inferential <-
    (nrow(fit$phi_table) <= 1L || isTRUE(fit$unit_support$phi_inference)) &&
    (nrow(fit$alpha_table) <= 1L || isTRUE(fit$unit_support$alpha_inference))
  have_se <- inferential && any(is.finite(fr$se_log_rho))
  lo <- if (have_se) lr - 1.96 * fr$se_log_rho else lr
  hi <- if (have_se) lr + 1.96 * fr$se_log_rho else lr
  labs <- paste0(fr$set, " \u00d7 ", fr$group)
  glev <- sort(unique(fr$group))
  colr <- .rr$pal[(match(fr$group, glev) - 1L) %% length(.rr$pal) + 1L]
  op <- par(mar = c(4.2, 9, 3.2, 1.5), mgp = c(2.5, 0.7, 0), tcl = -0.25,
            las = 1, col.axis = .rr$ink, col.lab = .rr$ink, col.main = .rr$ink,
            font.main = 2, cex.main = 1.15)
  on.exit(par(op))
  plot(NA, xlim = range(c(lo, hi, 0), na.rm = TRUE) + c(-0.1, 0.1),
       ylim = c(0.5, n + 0.5),
       xlab = "log unit (log rho)", ylab = "", axes = FALSE, main = "")
  if (!have_se)
    mtext("unit intervals unavailable for this fit", side = 3,
          line = 0.2, adj = 1, cex = 0.72, col = .rr$soft, font = 3)
  abline(h = seq_len(n), col = .rr$grid, lwd = 0.8)
  abline(v = 0, lty = 2, col = .rr$soft)
  .rr_axis(1)
  axis(2, at = seq_len(n), labels = labs, cex.axis = 0.75,
       col = .rr$grid, col.ticks = NA)
  misfit <- !is.na(fr$fit_resid) & abs(fr$fit_resid) > band
  if (have_se)
    segments(lo, seq_len(n), hi, seq_len(n), lwd = 2.2,
             col = ifelse(misfit, .rr$red, .rr$soft))
  points(lr, seq_len(n), pch = 21, cex = 1.5,
         bg = ifelse(misfit, .rr$red, colr), col = "white", lwd = 1.2)
  .rr_legend("bottomright", glev, pch = 21,
             pt.bg = .rr$pal[seq_along(glev)], col = "white", pt.cex = 1.2)
  if (any(misfit))
    mtext(sprintf("%d frame(s) with |pooled fit residual| > %.1f", sum(misfit), band),
          side = 3, line = 0.2, adj = 0, cex = 0.8, col = .rr$red)
  invisible(NULL)
}

#' Plot an item's characteristic curves across frames
#'
#' Plots the model expected-score curve for one item in each frame,
#' with observed class-interval means overlaid. Differences between the model
#' curves reflect the fitted frame units. A nominated non-frame person factor
#' separates the observed means within each frame for a DIF display.
#'
#' @param fit A fitted object from \code{\link{rasch_efrm}}.
#' @param item Underlying item name.
#' @param n_groups Number of class intervals for the observed means. By
#'   default, use the fit's count, with at least two requested intervals.
#'   Tied locations may produce fewer intervals.
#' @param grid Logit grid.
#' @param group Optional person grouping vector, or one or more names of
#'   non-frame factors nominated in the fit. Several names define their
#'   factor-combination cells.
#' @return Called for its plotting side effect; invisibly \code{NULL}.
#' @examples
#' \donttest{
#' # see ?rasch_efrm for a complete simulated example
#' }
#' @export
plot_icc_frames <- function(fit, item, n_groups = NULL,
                            grid = NULL, group = NULL) {
  if (!inherits(fit, "rasch_efrm")) stop("plot_icc_frames needs a rasch_efrm fit")
  if (!is.atomic(item) || !is.null(dim(item)) || length(item) != 1L ||
      is.na(item))
    stop("`item` must name exactly one item")
  item <- .role_text_values(item)
  if (!nzchar(item)) stop("`item` must name exactly one item")
  n_groups <- if (is.null(n_groups)) max(2L, fit$n_groups)
              else .check_whole(n_groups, "n_groups", 2)
  grid <- .model_grid(fit, grid, half_width = 5)
  vm <- fit$virtual_map
  rows <- which(vm$item == item)
  if (!length(rows)) stop("no such item: ", item)
  group_label <- NULL
  if (.role_columns(group,
                    if (is.null(fit$factors)) character(0) else names(fit$factors),
                    nrow(fit$X))) {
    bad <- intersect(group, fit$frame_group %||% character(0))
    if (length(bad))
      stop("frame-defining factors belong to the frame curves, not the DIF overlay: ",
           paste(bad, collapse = ", "))
    if (is.null(fit$factors) || !all(group %in% names(fit$factors)))
      stop("every named group must be a person factor in the fit")
    group_label <- paste(group, collapse = " x ")
    group <- if (length(group) == 1L) fit$factors[[group]] else
      .factor_cells(fit$factors[group], sep = ":")
  }
  if (!is.null(group) && length(group) != nrow(fit$X))
    stop("group must have one value per person")
  gf <- if (is.null(group)) NULL else droplevels(factor(group))
  if (!is.null(gf) && nlevels(gf) < 2L)
    stop("group must contain at least two observed levels")
  thr <- fit$thresholds_arbitrary
  tau_i <- thr$delta[thr$item == item]
  mmax <- length(tau_i)
  op <- .rr_canvas(range(grid), c(0, mmax),
                   "Person location (logits, common unit)", "Expected score",
                   if (is.null(group_label)) item else
                     sprintf("%s by %s", item, group_label))
  on.exit(par(op))
  th_all <- fit$person$theta
  frame_cols <- .rr$pal[(seq_along(rows) - 1L) %% length(.rr$pal) + 1L]
  group_pch <- if (is.null(gf)) 21 else
    rep(c(21:25, 0:20), length.out = nlevels(gf))
  for (j in seq_along(rows)) {
    v <- rows[j]
    rho <- fit$disc[v]
    colr <- frame_cols[j]
    Ecurve <- vapply(grid, function(t) item_moments(t, tau_i, disc = rho)$E, 0)
    lines(grid, Ecurve, lwd = 2.6, col = colr)
    x <- fit$X[, v]
    candidate <- !is.na(th_all) & !is.na(x) & !fit$person$extreme
    if (sum(candidate) >= 2 * n_groups) {
      ng <- min(n_groups, max(2, floor(sum(candidate) / 15)))
      ci <- .class_intervals(ifelse(is.na(x), NA_real_, th_all),
                             fit$person$extreme, ng)
      ok <- !is.na(ci)
      if (is.null(gf)) {
        points(tapply(th_all[ok], ci[ok], mean), tapply(x[ok], ci[ok], mean),
               pch = 21, bg = colr, col = "white", cex = 1.3, lwd = 1)
      } else {
        g <- gf[ok]
        for (k in seq_len(nlevels(gf))) {
          take <- g == levels(gf)[k]
          if (sum(take) < 2L) next
          points(tapply(th_all[ok][take], ci[ok][take], mean),
                 tapply(x[ok][take], ci[ok][take], mean),
                 pch = group_pch[k], bg = colr, col = colr,
                 cex = 1.15, lwd = 1)
        }
      }
    }
  }
  .rr_legend("topleft",
             sprintf("%s (rho %.3f)", vm$group[rows], fit$disc[rows]),
             lwd = 2.6, col = frame_cols)
  if (is.null(gf)) {
    .rr_legend("bottomright", c("Model", "Observed"),
               lwd = c(2.6, NA), pch = c(NA, 21),
               pt.bg = c(NA, .rr$ink), col = c(.rr$ink, "white"),
               pt.cex = 1.15)
  } else {
    .rr_legend("bottomright", levels(gf), pch = group_pch,
               pt.bg = "white", col = .rr$ink, pt.cex = 1.05)
  }
  invisible(NULL)
}
