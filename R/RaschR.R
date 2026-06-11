# RaschR :: a pairwise Rasch analysis engine in R
# ===========================================================================
# Clean-room from published theory: Choppin (1968, 1985); Andrich & Luo (1993)
# pairwise estimation; Warm (1989) weighted likelihood; Andrich item-trait
# interaction chi-square; Saaty principal eigenvector. No commercial source.
# Australian English; no em dashes by house style.
#
# Conditional pairwise identity underpinning estimation: for items i, j and
# thresholds k (on i), l (on j), patterns A = (X_i=k-1, X_j=l) and
# B = (X_i=k, X_j=l-1) share pair-total k+l-1, so theta cancels and
#   P(A)/P(B) = exp(tau_ik - tau_jl).
# Dichotomous data is the special case m_i = 1.
# ===========================================================================

# ---------------------------------------------------------------------------
# ESTIMATION: pairwise comparisons and two solvers
# ---------------------------------------------------------------------------

#' Enumerate item-category thresholds
#'
#' Builds the index mapping each item-category threshold to a global id, given
#' the maximum score of each item.
#'
#' @param m Integer vector of maximum scores per item (1 for dichotomous items).
#' @return A data frame with columns \code{id}, \code{item}, and \code{k} (the
#'   within-item threshold number).
#' @examples
#' threshold_index(c(1, 3, 2))
#' @export
threshold_index <- function(m) {
  thr <- do.call(rbind, lapply(seq_along(m), function(i)
    if (m[i] >= 1) data.frame(item = i, k = seq_len(m[i])) else NULL))
  thr$id <- seq_len(nrow(thr)); thr[, c("id", "item", "k")]
}

# Pairwise threshold-difference matrix D[p,q] ~ tau_p - tau_q (cross-item only)
# with inverse-variance weights W[p,q] = (cA*cB)/(cA+cB)  (var of log ratio ~ 1/cA+1/cB).
#' Pairwise threshold-difference matrix
#'
#' Computes the conditional pairwise log-difference between every pair of
#' item-category thresholds, with inverse-variance weights. The person
#' parameter cancels in each comparison, giving sample-free estimates.
#'
#' @param X Persons-by-items integer score matrix (categories from 0).
#' @param thr Threshold index from \code{\link{threshold_index}}.
#' @param cont Continuity constant added only to comparisons containing an empty
#'   cell.
#' @return A list with the difference matrix \code{D} and weight matrix \code{W}.
#' @examples
#' set.seed(1)
#' X <- matrix(rbinom(300 * 4, 1, 0.5), 300, 4)
#' poly_diff_matrix(X, threshold_index(rep(1, 4)))$D
#' @export
poly_diff_matrix <- function(X, thr, cont = 0.5) {
  X <- as.matrix(X); M <- nrow(thr)
  D <- matrix(NA_real_, M, M); W <- matrix(0, M, M)
  for (p in seq_len(M)) {
    i <- thr$item[p]; k <- thr$k[p]
    for (q in seq_len(M)) {
      j <- thr$item[q]; l <- thr$k[q]
      if (i == j) next
      both <- !is.na(X[, i]) & !is.na(X[, j])
      cA <- sum(X[both, i] == (k - 1L) & X[both, j] == l)
      cB <- sum(X[both, i] == k        & X[both, j] == (l - 1L))
      if (cA + cB > 0) {
        a <- cA; b <- cB
        if (cA == 0 || cB == 0) { a <- cA + cont; b <- cB + cont }  # only if empty
        D[p, q] <- log(a / b)
        W[p, q] <- (a * b) / (a + b)                                # inverse-variance
      }
    }
  }
  list(D = D, W = W)
}

# (a) Weighted least squares on the difference graph (inverse-variance weights).
#' Solve the pairwise threshold system by weighted least squares
#'
#' @param D Threshold-difference matrix from \code{\link{poly_diff_matrix}}.
#' @param W Optional inverse-variance weight matrix.
#' @return Numeric vector of threshold locations, centred at zero.
#' @examples
#' set.seed(1)
#' X <- matrix(rbinom(300 * 4, 1, 0.5), 300, 4)
#' dm <- poly_diff_matrix(X, threshold_index(rep(1, 4)))
#' solve_LS(dm$D, dm$W)
#' @export
solve_LS <- function(D, W = NULL) {
  M <- nrow(D); rows <- which(!is.na(D) & upper.tri(D), arr.ind = TRUE)
  C <- matrix(0, nrow(rows), M); d <- numeric(nrow(rows)); wt <- numeric(nrow(rows))
  for (r in seq_len(nrow(rows))) {
    p <- rows[r, 1]; q <- rows[r, 2]
    C[r, p] <- 1; C[r, q] <- -1; d[r] <- D[p, q]
    wt[r] <- if (is.null(W)) 1 else W[p, q]
  }
  sw <- sqrt(wt)
  tau <- c(qr.coef(qr((C[, -M, drop = FALSE]) * sw), d * sw), 0); tau[is.na(tau)] <- 0
  tau - mean(tau)
}
# (b) Weighted reciprocal averaging: iterative form of the principal-eigenvector
#     solution; coincides with (a) for incomplete (polytomous) designs.
#' Solve the pairwise threshold system by weighted reciprocal averaging
#'
#' The iterative form of the principal-eigenvector solution; it coincides with
#' \code{\link{solve_LS}} once the comparison design is structurally incomplete.
#'
#' @param D Threshold-difference matrix from \code{\link{poly_diff_matrix}}.
#' @param W Optional inverse-variance weight matrix.
#' @param tol Convergence tolerance.
#' @param maxit Maximum number of iterations.
#' @return Numeric vector of threshold locations, centred at zero.
#' @examples
#' set.seed(1)
#' X <- matrix(rbinom(300 * 4, 1, 0.5), 300, 4)
#' dm <- poly_diff_matrix(X, threshold_index(rep(1, 4)))
#' solve_RA(dm$D, dm$W)
#' @export
solve_RA <- function(D, W = NULL, tol = 1e-11, maxit = 8000) {
  M <- nrow(D); tau <- rep(0, M); def <- !is.na(D)
  if (is.null(W)) W <- matrix(1, M, M)
  for (it in seq_len(maxit)) {
    new <- vapply(seq_len(M), function(p) {
      q <- which(def[p, ]); if (!length(q)) return(tau[p])
      sum(W[p, q] * (D[p, q] + tau[q])) / sum(W[p, q])
    }, numeric(1))
    new <- new - mean(new)
    if (max(abs(new - tau)) < tol) return(new)
    tau <- new
  }
  tau
}

# Dichotomous principal-eigenvector solver (clean closed form, m_i = 1).
#' Dichotomous item difficulties by the principal eigenvector
#'
#' Estimates dichotomous Rasch item difficulties as the Perron eigenvector of
#' the reciprocal pairwise comparison matrix (the Saaty solution).
#'
#' @param X Persons-by-items 0/1 matrix.
#' @param cont Continuity constant for the comparison counts.
#' @return Named numeric vector of item difficulties, centred at zero.
#' @examples
#' set.seed(1)
#' d <- seq(-2, 2, length.out = 6)
#' X <- matrix(rbinom(500 * 6, 1, plogis(outer(rnorm(500), d, "-"))), 500, 6)
#' colnames(X) <- paste0("I", 1:6)
#' est_eigen_dich(X)
#' @export
est_eigen_dich <- function(X, cont = 0.5) {
  X <- as.matrix(X); L <- ncol(X)
  N <- matrix(0, L, L)
  for (i in seq_len(L)) for (j in seq_len(L)) if (i != j) {
    both <- !is.na(X[, i]) & !is.na(X[, j])
    N[i, j] <- sum(X[both, i] == 1L & X[both, j] == 0L)
  }
  Nc <- N + cont; diag(Nc) <- 0
  R <- matrix(1, L, L)
  for (i in seq_len(L)) for (j in seq_len(L)) if (i != j) R[i, j] <- Nc[i, j] / Nc[j, i]
  ev <- eigen(R, symmetric = FALSE)
  v  <- abs(Re(ev$vectors[, which.max(Re(ev$values))]))
  d  <- -log(v / sum(v)); setNames(d - mean(d), colnames(X))
}

# ---------------------------------------------------------------------------
# CATEGORY MOMENTS (E, variance, 3rd and 4th central moments)
# ---------------------------------------------------------------------------
#' Category-score moments for a polytomous item
#'
#' @param theta Person location, in logits.
#' @param tau_i Numeric vector of the item's threshold parameters.
#' @return A list with category probabilities \code{P}, expected score \code{E},
#'   variance \code{V}, and third and fourth central moments \code{mu3},
#'   \code{mu4}.
#' @examples
#' item_moments(0.5, c(-1, 0, 1))
#' @export
item_moments <- function(theta, tau_i) {
  m <- length(tau_i); x <- 0:m
  num <- exp(x * theta - c(0, cumsum(tau_i))); P <- num / sum(num)
  E <- sum(x * P); d <- x - E
  list(P = P, E = E, V = sum(d^2 * P), mu3 = sum(d^3 * P), mu4 = sum(d^4 * P))
}

# ---------------------------------------------------------------------------
# PERSON ESTIMATION: Warm's WLE (finite at extreme scores) + standard errors
# WLE solves:  (R - sum E_i) + (sum mu3_i) / (2 sum V_i) = 0
# ---------------------------------------------------------------------------
#' Warm's weighted likelihood person estimates
#'
#' Computes the weighted likelihood estimate (WLE) of person location for every
#' possible raw score, with standard errors. WLE estimates are finite at the
#' extreme (zero and maximum) scores, unlike the maximum likelihood estimate.
#'
#' @param tau_list List of per-item threshold vectors.
#' @return A list with \code{theta} and \code{se}, each named by raw score.
#' @examples
#' person_wle(list(c(-1, 0), c(-0.5, 0.5), c(0, 1)))
#' @export
person_wle <- function(tau_list) {
  Smax <- sum(vapply(tau_list, length, 1L))
  theta <- se <- setNames(rep(NA_real_, Smax + 1L), as.character(0:Smax))
  for (R in 0:Smax) {
    g <- function(th) {
      mo <- lapply(tau_list, item_moments, theta = th)
      E  <- sum(vapply(mo, `[[`, 0, "E"));  V <- sum(vapply(mo, `[[`, 0, "V"))
      m3 <- sum(vapply(mo, `[[`, 0, "mu3"))
      (R - E) + m3 / (2 * V)
    }
    root <- tryCatch(uniroot(g, c(-20, 20), tol = 1e-9)$root, error = function(e) NA_real_)
    theta[as.character(R)] <- root
    if (!is.na(root)) {
      V <- sum(vapply(lapply(tau_list, item_moments, theta = root), `[[`, 0, "V"))
      se[as.character(R)] <- 1 / sqrt(V)               # WLE SE ~ 1/sqrt(information)
    }
  }
  list(theta = theta, se = se)
}

# ---------------------------------------------------------------------------
# FIT: standardised residuals, infit/outfit mean squares, standardised
# fit residual (Wilson-Hilferty normalising transform, ~ N(0,1) under fit).
# The exact scaling constant is a placeholder until calibrated against
# output from established Rasch software.
# ---------------------------------------------------------------------------
.wh <- function(ms, q) (ms^(1/3) - 1) * (3 / q) + (q / 3)   # mean square -> z

fit_from_moments <- function(X, score, Emat, Vmat, M4mat) {
  L <- ncol(X)
  out <- data.frame(item = colnames(X), infit_ms = NA_real_, outfit_ms = NA_real_,
                    infit_z = NA_real_, outfit_z = NA_real_, fit_resid = NA_real_)
  for (i in seq_len(L)) {
    obs <- !is.na(X[, i])
    E <- Emat[cbind(score[obs] + 1L, i)]; V <- Vmat[cbind(score[obs] + 1L, i)]
    C4 <- M4mat[cbind(score[obs] + 1L, i)]
    res <- X[obs, i] - E; n <- length(res)
    outfit <- mean(res^2 / V); infit <- sum(res^2) / sum(V)
    qo <- sqrt(max(sum(C4 / V^2) / n^2 - 1 / n, 1e-8))
    qi <- sqrt(max(sum(C4 - V^2) / sum(V)^2, 1e-8))
    out$outfit_ms[i] <- outfit; out$infit_ms[i] <- infit
    out$outfit_z[i]  <- .wh(outfit, qo); out$infit_z[i] <- .wh(infit, qi)
    out$fit_resid[i] <- sum(res / sqrt(V)) / sqrt(n)   # standardised aggregate
  }
  out
}

# ---------------------------------------------------------------------------
# ITEM-TRAIT INTERACTION CHI-SQUARE over class intervals (+ Bonferroni)
# chi2_i = sum_g n_gi (Obar_gi - Ebar_gi)^2 / Vbar_gi ;  df_i = G - 1
# ---------------------------------------------------------------------------
item_trait <- function(X, score, theta_person, Emat, Vmat, n_groups = 10,
                       adjust_N = NA) {
  L <- ncol(X)
  rk <- rank(theta_person, ties.method = "first")
  brk <- unique(quantile(rk, probs = seq(0, 1, length.out = n_groups + 1)))
  g <- cut(rk, breaks = brk, include.lowest = TRUE, labels = FALSE)
  G <- max(g, na.rm = TRUE)
  chi <- setNames(numeric(L), colnames(X))
  for (i in seq_len(L)) for (gg in seq_len(G)) {
    sel <- which(g == gg & !is.na(X[, i]))
    if (length(sel) < 2) next
    Obar <- mean(X[sel, i])
    Ebar <- mean(Emat[cbind(score[sel] + 1L, i)])
    Vbar <- mean(Vmat[cbind(score[sel] + 1L, i)])
    chi[i] <- chi[i] + length(sel) * (Obar - Ebar)^2 / Vbar
  }
  df_i <- G - 1
  # Sample-size adjustment: rescale chi-square to a reference N so
  # that probabilities are comparable across studies of different size.
  if (!is.na(adjust_N)) chi <- chi * (adjust_N / nrow(X))
  p <- pchisq(chi, df_i, lower.tail = FALSE)
  data.frame(item = colnames(X), chisq = chi, df = df_i, p = p,
             p_bonf = pmin(p * L, 1), misfit_bonf = p < 0.05 / L)
}

# ---------------------------------------------------------------------------
# PERSON SEPARATION INDEX (separation reliability)
# ---------------------------------------------------------------------------
psi_index <- function(theta_person, se_person) {
  ok <- !is.na(theta_person) & !is.na(se_person)
  vt <- var(theta_person[ok]); mse <- mean(se_person[ok]^2)
  psi <- max((vt - mse) / vt, 0)
  sep <- if (psi < 1) sqrt(psi / (1 - psi)) else Inf
  list(PSI = psi, separation = sep, var_theta = vt, mean_error_var = mse)
}

# ---------------------------------------------------------------------------
# THRESHOLD / CATEGORY DIAGNOSTICS: ordering, reversals, modal categories
# ---------------------------------------------------------------------------
threshold_diag <- function(X, thr, grid = seq(-8, 8, by = 0.05)) {
  items <- sort(unique(thr$item)); res <- list()
  for (i in items) {
    tau_i <- thr$tau[thr$item == i]
    modal <- unique(vapply(grid, function(th) which.max(item_moments(th, tau_i)$P) - 1L, 1L))
    res[[as.character(i)]] <- list(
      item = colnames(X)[i], thresholds = tau_i,
      ordered = all(diff(tau_i) > 0),
      reversed_at = which(diff(tau_i) <= 0) + 1L,
      never_modal_categories = setdiff(0:length(tau_i), modal),
      category_counts = as.integer(table(factor(X[, i], levels = 0:length(tau_i)))))
  }
  res
}

# ---------------------------------------------------------------------------
# CONSTRAINED RSM FIT: tau_ik = location_i + step_k, estimated jointly from
# the pairwise comparisons (removes post-hoc decomposition shrinkage).
# Each comparison: D[p,q] = (lam_i - lam_j) + (kap_k - kap_l).
# ---------------------------------------------------------------------------
#' Fit the rating scale model under constraint
#'
#' Estimates a rating scale model directly from the pairwise comparisons,
#' decomposing each threshold into an item location plus a common step.
#'
#' @param X Persons-by-items integer score matrix with equal maximum score
#'   across items.
#' @param cont Continuity constant.
#' @return A list with item \code{location}s, common \code{step}s, and the
#'   implied item-by-threshold matrix \code{tau}.
#' @examples
#' set.seed(1)
#' simP <- function(th, t) { x <- 0:length(t); p <- exp(x * th - c(0, cumsum(t))); p / sum(p) }
#' loc <- seq(-1, 1, length.out = 6); step <- c(-0.8, -0.2, 0.3, 0.7); th <- rnorm(400)
#' X <- sapply(loc, function(b) sapply(th, function(t) sample(0:4, 1, prob = simP(t, b + step))))
#' colnames(X) <- paste0("R", 1:6)
#' rasch_rsm(X)$location
#' @export
rasch_rsm <- function(X, cont = 0.5) {
  X <- as.matrix(X); storage.mode(X) <- "integer"
  m <- apply(X, 2, max, na.rm = TRUE)
  if (length(unique(m)) != 1L) stop("RSM requires equal max score across items")
  mm <- m[1]; L <- ncol(X); thr <- threshold_index(m)
  dm <- poly_diff_matrix(X, thr, cont); D <- dm$D; W <- dm$W
  rows <- which(!is.na(D) & upper.tri(D), arr.ind = TRUE)
  nC <- L + mm; C <- matrix(0, nrow(rows), nC); d <- numeric(nrow(rows)); wt <- numeric(nrow(rows))
  for (r in seq_len(nrow(rows))) {
    p <- rows[r, 1]; q <- rows[r, 2]
    i <- thr$item[p]; k <- thr$k[p]; j <- thr$item[q]; l <- thr$k[q]
    C[r, i] <- C[r, i] + 1; C[r, j] <- C[r, j] - 1
    C[r, L + k] <- C[r, L + k] + 1; C[r, L + l] <- C[r, L + l] - 1
    d[r] <- D[p, q]; wt[r] <- W[p, q]
  }
  keep <- setdiff(seq_len(nC), c(L, L + mm)); sw <- sqrt(wt)
  beta <- numeric(nC); beta[keep] <- qr.coef(qr((C[, keep, drop = FALSE]) * sw), d * sw)
  beta[is.na(beta)] <- 0
  lam <- beta[1:L]; kap <- beta[(L + 1):(L + mm)]
  ck <- mean(kap); kap <- kap - ck; lam <- lam + ck; lam <- lam - mean(lam)
  list(location = setNames(lam, colnames(X)), step = kap,
       tau = outer(lam, kap, "+"))
}

# ---------------------------------------------------------------------------
# TOP-LEVEL ANALYSIS
# ---------------------------------------------------------------------------
#' Fit and diagnose a Rasch model
#'
#' Runs a full Rasch analysis: pairwise conditional item estimation, Warm
#' weighted likelihood person estimates, fit statistics, the item-trait
#' interaction chi-square, person separation, and threshold diagnostics.
#'
#' @param X Persons-by-items integer score matrix, categories starting at 0.
#'   Missing values are allowed (pairwise deletion in estimation).
#' @param model Either \code{"PCM"} (partial credit) or \code{"RSM"} (rating
#'   scale, fitted under constraint).
#' @param solver Pairwise solver for PCM: \code{"LS"} (weighted least squares)
#'   or \code{"RA"} (reciprocal averaging).
#' @param n_groups Number of class intervals for the item-trait chi-square.
#' @param cont Continuity constant for the pairwise comparisons.
#' @param adjust_N Optional reference sample size; if supplied, item-trait
#'   chi-squares are rescaled to this size.
#' @return A list containing the item summary (\code{items}), \code{thresholds},
#'   person estimates and standard errors, item-trait results, person separation
#'   index (\code{psi}), threshold diagnostics, and the standardised residual
#'   matrix.
#' @examples
#' set.seed(1)
#' d <- seq(-2, 2, length.out = 8)
#' X <- matrix(rbinom(500 * 8, 1, plogis(outer(rnorm(500), d, "-"))), 500, 8)
#' colnames(X) <- paste0("I", 1:8)
#' fit <- rasch(X, model = "PCM")
#' fit$items
#' fit$psi$PSI
#' @export
rasch <- function(X, model = c("PCM", "RSM"), solver = c("LS", "RA"),
                         n_groups = 10, cont = 0.5, adjust_N = NA) {
  model <- match.arg(model); solver <- match.arg(solver)
  X <- as.matrix(X); storage.mode(X) <- "integer"
  if (is.null(colnames(X))) colnames(X) <- sprintf("I%02d", seq_len(ncol(X)))
  m <- apply(X, 2, max, na.rm = TRUE); L <- ncol(X)

  thr <- threshold_index(m)
  if (model == "RSM") {
    rs <- rasch_rsm(X, cont)
    tau_list <- lapply(seq_len(L), function(i) rs$tau[i, ])
    thr$tau <- unlist(tau_list)
  } else {
    dm <- poly_diff_matrix(X, thr, cont)
    thr$tau <- if (solver == "LS") solve_LS(dm$D, dm$W) else solve_RA(dm$D, dm$W)
    tau_list <- lapply(seq_len(L), function(i) thr$tau[thr$item == i])
  }

  pe <- person_wle(tau_list); Smax <- sum(m)
  Emat <- Vmat <- M4mat <- matrix(NA_real_, Smax + 1L, L)
  for (s in 0:Smax) if (!is.na(pe$theta[as.character(s)])) {
    mo <- lapply(tau_list, item_moments, theta = pe$theta[as.character(s)])
    Emat[s + 1L, ]  <- vapply(mo, `[[`, 0, "E")
    Vmat[s + 1L, ]  <- vapply(mo, `[[`, 0, "V")
    M4mat[s + 1L, ] <- vapply(mo, `[[`, 0, "mu4")
  }
  raw <- rowSums(X, na.rm = TRUE)
  theta_person <- pe$theta[as.character(raw)]; se_person <- pe$se[as.character(raw)]

  # standardised residual matrix Z[n,i] = (x - E) / sqrt(V) at each person's fit
  Eobs <- Emat[raw + 1L, , drop = FALSE]; Vobs <- Vmat[raw + 1L, , drop = FALSE]
  Z <- (X - Eobs) / sqrt(Vobs)

  fit <- fit_from_moments(X, raw, Emat, Vmat, M4mat)
  it  <- item_trait(X, raw, theta_person, Emat, Vmat, n_groups = n_groups, adjust_N = adjust_N)
  ps  <- psi_index(theta_person, se_person)
  td  <- threshold_diag(X, thr)

  items <- data.frame(item = colnames(X), max = m,
                      location = vapply(tau_list, mean, 0),
                      infit_ms = fit$infit_ms, outfit_ms = fit$outfit_ms,
                      fit_resid = fit$fit_resid,
                      chisq = it$chisq, chisq_p = it$p, misfit = it$misfit_bonf)
  list(model = model, X = X, items = items, thresholds = thr, tau_list = tau_list,
       theta_by_score = pe$theta, se_by_score = pe$se,
       theta_person = theta_person, se_person = se_person, raw = raw,
       Emat = Emat, Vmat = Vmat, residuals = Z, n_groups = n_groups,
       item_trait = it, psi = ps, thresholds_diag = td,
       total_chisq = sum(it$chisq), total_df = L * max(it$df))
}

# ===========================================================================
# DIMENSIONALITY AND LOCAL DEPENDENCE
# ===========================================================================

# Residual correlation matrix. Under unidimensionality + local independence the
# off-diagonal correlations sit near -1/(L-1); large positive values flag local
# dependence between item pairs.
#' Residual correlations for local dependence
#'
#' @param fit A fitted object from \code{\link{rasch}}.
#' @param flag Excess-correlation threshold above the average for flagging a
#'   dependent item pair.
#' @return A list with the residual correlation \code{matrix}, the \code{average}
#'   off-diagonal correlation, and a table of \code{flagged} pairs.
#' @examples
#' set.seed(1)
#' d <- seq(-2, 2, length.out = 8)
#' X <- matrix(rbinom(500 * 8, 1, plogis(outer(rnorm(500), d, "-"))), 500, 8)
#' colnames(X) <- paste0("I", 1:8)
#' residual_correlations(rasch(X))$average
#' @export
residual_correlations <- function(fit, flag = 0.2) {
  Z <- fit$residuals
  R <- cor(Z, use = "pairwise.complete.obs")
  off <- R[upper.tri(R)]; avg <- mean(off, na.rm = TRUE)
  idx <- which(upper.tri(R) & (R - avg) > flag, arr.ind = TRUE)
  pairs <- if (nrow(idx))
    data.frame(item_a = colnames(Z)[idx[, 1]], item_b = colnames(Z)[idx[, 2]],
               resid_cor = R[idx], excess = R[idx] - avg) else
    data.frame(item_a = character(), item_b = character(),
               resid_cor = numeric(), excess = numeric())
  list(matrix = R, average = avg, flagged = pairs[order(-pairs$resid_cor), ])
}

# PCA of the residual correlation matrix. The first residual contrast (PC1)
# carries any second dimension; items with opposing loadings define the split.
#' Principal components of the residual correlations
#'
#' @param fit A fitted object from \code{\link{rasch}}.
#' @return A list with the residual \code{eigenvalues}, their \code{prop}ortions,
#'   the first-contrast \code{loadings}, and the \code{first_eigen}value.
#' @examples
#' set.seed(1)
#' d <- seq(-2, 2, length.out = 8)
#' X <- matrix(rbinom(500 * 8, 1, plogis(outer(rnorm(500), d, "-"))), 500, 8)
#' colnames(X) <- paste0("I", 1:8)
#' residual_pca(rasch(X))$first_eigen
#' @export
residual_pca <- function(fit) {
  R <- cor(fit$residuals, use = "pairwise.complete.obs")
  ev <- eigen(R, symmetric = TRUE)
  loadings <- ev$vectors[, 1] * sqrt(pmax(ev$values[1], 0))
  data.frame(item = colnames(fit$residuals), pc1_loading = loadings) -> ld
  list(eigenvalues = ev$values, prop = ev$values / sum(ev$values),
       loadings = ld[order(-ld$pc1_loading), ], first_eigen = ev$values[1])
}

# Dimensionality test (Smith 2002): split items by PC1 loading sign, estimate
# each person on each subset, then a paired t per person. The proportion of
# |t| > 1.96 beyond 5% (binomial CI excluding 0.05) signals multidimensionality.
#' Residual-component dimensionality test
#'
#' Splits items by the sign of their first residual-contrast loading, estimates
#' each person on each subset, and compares the two sets of estimates with a
#' per-person t-test (Smith 2002). A significant-test proportion whose lower
#' confidence bound exceeds \code{alpha} signals multidimensionality.
#'
#' @param fit A fitted object from \code{\link{rasch}}.
#' @param alpha Nominal significance level for the per-person t-tests.
#' @return A list with the proportion of significant tests, its confidence
#'   interval, the item split, and a \code{multidimensional} verdict.
#' @examples
#' set.seed(1)
#' d <- seq(-2, 2, length.out = 8)
#' X <- matrix(rbinom(500 * 8, 1, plogis(outer(rnorm(500), d, "-"))), 500, 8)
#' colnames(X) <- paste0("I", 1:8)
#' dimensionality_test(rasch(X))$multidimensional
#' @export
dimensionality_test <- function(fit, alpha = 0.05) {
  pca <- residual_pca(fit); X <- fit$X
  pos <- which(pca$loadings$pc1_loading[match(colnames(X), pca$loadings$item)] > 0)
  neg <- setdiff(seq_len(ncol(X)), pos)
  if (length(pos) < 2 || length(neg) < 2)
    return(list(note = "need >= 2 items each side of the contrast"))
  subset_theta <- function(cols) {
    tl <- fit$tau_list[cols]; pe <- person_wle(tl)
    r <- rowSums(X[, cols, drop = FALSE], na.rm = TRUE)
    list(theta = pe$theta[as.character(r)], se = pe$se[as.character(r)])
  }
  a <- subset_theta(pos); b <- subset_theta(neg)
  ok <- !is.na(a$theta) & !is.na(b$theta) & !is.na(a$se) & !is.na(b$se)
  t <- (a$theta[ok] - b$theta[ok]) / sqrt(a$se[ok]^2 + b$se[ok]^2)
  prop <- mean(abs(t) > qnorm(1 - alpha / 2)); n <- sum(ok)
  ci <- prop + c(-1, 1) * 1.96 * sqrt(prop * (1 - prop) / n)
  list(prop_significant = prop, ci = ci, n = n,
       multidimensional = ci[1] > alpha,
       items_positive = colnames(X)[pos], items_negative = colnames(X)[neg],
       first_eigenvalue = pca$first_eigen)
}

# DIF by two-way residual ANOVA: residual ~ group * class-interval, per item.
# group main effect = uniform DIF; group:class interaction = non-uniform DIF.
#' Differential item functioning by two-way residual ANOVA
#'
#' For each item, analyses the standardised residuals by person group and trait
#' class interval. The group main effect indicates uniform DIF and the
#' group-by-interval interaction indicates non-uniform DIF, with Bonferroni
#' flagging.
#'
#' @param fit A fitted object from \code{\link{rasch}}.
#' @param group Grouping variable (factor or vector) of length equal to the
#'   number of persons.
#' @param n_groups Number of trait class intervals; defaults to the value used
#'   in the fit.
#' @return A data frame of uniform and non-uniform DIF statistics per item.
#' @examples
#' set.seed(1); n <- 600
#' d <- seq(-2, 2, length.out = 8); g <- rep(c("a", "b"), each = n / 2)
#' sh <- matrix(0, n, 8); sh[g == "b", 3] <- 1
#' X <- matrix(rbinom(n * 8, 1, plogis(outer(rnorm(n), d, "-") - sh)), n, 8)
#' colnames(X) <- paste0("I", 1:8)
#' dif_anova(rasch(X), group = g)
#' @export
dif_anova <- function(fit, group, n_groups = NULL) {
  Z <- fit$residuals; L <- ncol(Z)
  if (is.null(n_groups)) n_groups <- fit$n_groups
  rk <- rank(fit$theta_person, ties.method = "first")
  brk <- unique(quantile(rk, probs = seq(0, 1, length.out = n_groups + 1)))
  ci <- factor(cut(rk, breaks = brk, include.lowest = TRUE, labels = FALSE))
  group <- factor(group)
  out <- data.frame(item = colnames(Z), F_uniform = NA_real_, p_uniform = NA_real_,
                    F_nonuniform = NA_real_, p_nonuniform = NA_real_)
  for (i in seq_len(L)) {
    d <- data.frame(z = Z[, i], g = group, ci = ci)
    d <- d[stats::complete.cases(d), ]
    a <- tryCatch(stats::anova(stats::lm(z ~ g * ci, data = d)), error = function(e) NULL)
    if (is.null(a)) next
    rn <- rownames(a)
    if ("g" %in% rn)    { out$F_uniform[i] <- a["g", "F value"];    out$p_uniform[i] <- a["g", "Pr(>F)"] }
    if ("g:ci" %in% rn) { out$F_nonuniform[i] <- a["g:ci", "F value"]; out$p_nonuniform[i] <- a["g:ci", "Pr(>F)"] }
  }
  out$uniform_DIF    <- out$p_uniform    < 0.05 / L   # Bonferroni
  out$nonuniform_DIF <- out$p_nonuniform < 0.05 / L
  out
}
