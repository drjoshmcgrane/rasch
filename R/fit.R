# RaschR :: fit statistics
# ===========================================================================
# Test-of-fit statistics for the pairwise analysis: standardised residuals;
# item and person
# fit residuals (Wilson-Hilferty standardised mean squares, approximately
# N(0,1) under fit, negative for overfit and positive for underfit); infit
# and outfit mean squares; the item-trait interaction chi-square and the
# class-interval ANOVA F; the person separation index with and without
# extremes; Cronbach's alpha; targeting; and the test information function.
# ===========================================================================

.wh <- function(ms, q) (ms^(1/3) - 1) * (3 / q) + (q / 3)   # mean square -> z

# Because each person's location is estimated from their own responses, the
# expected squared standardised residual in cell (n, i) is close to
# 1 - V_ni / sum_j V_nj rather than 1; mean squares are rescaled by this
# information share so they centre on 1 under fit.
.z2_expectation <- function(mo, Z, disc = NULL) {
  Vobs <- mo$V; Vobs[is.na(Z)] <- NA
  if (!is.null(disc)) Vobs <- sweep(Vobs, 2, disc^2, "*")
  share <- Vobs / rowSums(Vobs, na.rm = TRUE)
  pmax(1 - share, 1e-4)
}

# Item fit from per-person model moments (observed cells only).
.item_fit <- function(X, Z, mo, disc = NULL) {
  L <- ncol(X)
  E2 <- .z2_expectation(mo, Z, disc)
  out <- data.frame(item = colnames(X), infit_ms = NA_real_, outfit_ms = NA_real_,
                    infit_z = NA_real_, outfit_z = NA_real_, fit_resid = NA_real_,
                    n = NA_integer_)
  for (i in seq_len(L)) {
    ok <- which(!is.na(Z[, i]))
    if (length(ok) < 3) next
    z2 <- Z[ok, i]^2
    V <- mo$V[ok, i]; C4 <- mo$M4[ok, i]; n <- length(ok)
    e2 <- E2[ok, i]
    outfit <- sum(z2) / sum(e2)
    infit  <- sum(z2 * V) / sum(e2 * V)
    qo <- sqrt(max(sum(C4 / V^2) / n^2 - 1 / n, 1e-8))
    qi <- sqrt(max(sum(C4 - V^2) / sum(V)^2, 1e-8))
    out$outfit_ms[i] <- outfit; out$infit_ms[i] <- infit
    out$outfit_z[i] <- .wh(outfit, qo); out$infit_z[i] <- .wh(infit, qi)
    # fit residual: the information-weighted form is well calibrated
    # (mean near 0, SD near 1 on model-true data)
    out$fit_resid[i] <- out$infit_z[i]
    out$n[i] <- n
  }
  out
}

# Person fit residuals: each person's standardised residuals across their
# observed items, summarised exactly as for items.
.person_fit <- function(X, Z, mo, disc = NULL) {
  N <- nrow(X)
  E2 <- .z2_expectation(mo, Z, disc)
  out <- data.frame(infit_ms = rep(NA_real_, N), outfit_ms = NA_real_,
                    fit_resid = NA_real_)
  for (n in seq_len(N)) {
    ok <- which(!is.na(Z[n, ]))
    if (length(ok) < 3) next
    z2 <- Z[n, ok]^2
    V <- mo$V[n, ok]; C4 <- mo$M4[n, ok]; k <- length(ok)
    e2 <- E2[n, ok]
    outfit <- sum(z2) / sum(e2)
    out$outfit_ms[n] <- outfit
    out$infit_ms[n] <- sum(z2 * V) / sum(e2 * V)
    qo <- sqrt(max(sum(C4 / V^2) / k^2 - 1 / k, 1e-8))
    out$fit_resid[n] <- .wh(outfit, qo)
  }
  out
}

# Class intervals over non-extreme person locations.
.class_intervals <- function(theta, extreme, n_groups) {
  g <- rep(NA_integer_, length(theta))
  use <- which(!is.na(theta) & !extreme)
  rk <- rank(theta[use], ties.method = "first")
  brk <- unique(quantile(rk, probs = seq(0, 1, length.out = n_groups + 1)))
  g[use] <- cut(rk, breaks = brk, include.lowest = TRUE, labels = FALSE)
  g
}

# Item-trait interaction chi-square and class-interval ANOVA F per item.
.item_trait <- function(X, Z, mo, ci, adjust_N = NA) {
  L <- ncol(X); G <- max(ci, na.rm = TRUE)
  chi <- setNames(numeric(L), colnames(X))
  Fv <- pF <- rep(NA_real_, L)
  for (i in seq_len(L)) {
    for (gg in seq_len(G)) {
      sel <- which(ci == gg & !is.na(X[, i]))
      if (length(sel) < 2) next
      Obar <- mean(X[sel, i])
      Ebar <- mean(mo$E[sel, i]); Vbar <- mean(mo$V[sel, i])
      chi[i] <- chi[i] + length(sel) * (Obar - Ebar)^2 / Vbar
    }
    d <- data.frame(z = Z[, i], g = factor(ci))
    d <- d[stats::complete.cases(d), ]
    a <- tryCatch(stats::anova(stats::lm(z ~ g, data = d)), error = function(e) NULL)
    if (!is.null(a) && "g" %in% rownames(a)) {
      Fv[i] <- a["g", "F value"]; pF[i] <- a["g", "Pr(>F)"]
    }
  }
  df_i <- G - 1
  n_used <- sum(!is.na(ci))
  if (!is.na(adjust_N)) chi <- chi * (adjust_N / n_used)
  p <- pchisq(chi, df_i, lower.tail = FALSE)
  data.frame(item = colnames(X), chisq = chi, df = df_i, p = p,
             p_bonf = pmin(p * L, 1), misfit_bonf = p < 0.05 / L,
             F_value = Fv, p_F = pF)
}

# Person separation index (separation reliability; Andrich 1982).
.psi <- function(theta, se, keep = TRUE) {
  ok <- !is.na(theta) & !is.na(se) & keep
  if (sum(ok) < 3) return(list(PSI = NA_real_, separation = NA_real_,
                               var_theta = NA_real_, mean_error_var = NA_real_,
                               n = sum(ok)))
  vt <- var(theta[ok]); mse <- mean(se[ok]^2)
  psi <- max((vt - mse) / vt, 0)
  sep <- if (psi < 1) sqrt(psi / (1 - psi)) else Inf
  list(PSI = psi, separation = sep, var_theta = vt, mean_error_var = mse,
       n = sum(ok))
}

# Cronbach's alpha on complete cases, reported alongside the PSI.
.alpha <- function(X) {
  Xc <- X[stats::complete.cases(X), , drop = FALSE]
  if (nrow(Xc) < 3 || ncol(Xc) < 2) return(list(alpha = NA_real_, n = nrow(Xc)))
  L <- ncol(Xc); vi <- apply(Xc, 2, var); vt <- var(rowSums(Xc))
  list(alpha = L / (L - 1) * (1 - sum(vi) / vt), n = nrow(Xc))
}

# Qualitative power-of-test-of-fit assessment, driven by the PSI.
.fit_power <- function(psi) {
  if (is.na(psi)) "unknown"
  else if (psi >= 0.9) "excellent"
  else if (psi >= 0.8) "good"
  else if (psi >= 0.7) "reasonable"
  else if (psi >= 0.5) "low"
  else "too low"
}

# Targeting summary: how well item thresholds cover the person distribution.
.targeting <- function(person, thresholds) {
  ok <- !is.na(person$theta)
  th <- person$theta[ok]
  list(person_mean = mean(th), person_sd = sd(th),
       person_mean_noext = mean(person$theta[ok & !person$extreme]),
       item_mean = 0,
       threshold_range = range(thresholds$tau),
       prop_below = mean(th < min(thresholds$tau)),
       prop_above = mean(th > max(thresholds$tau)))
}

#' Test information function
#'
#' Fisher information of the whole test over a grid of person locations, with
#' the corresponding standard error of measurement.
#'
#' @param fit A fitted object from \code{\link{rasch}}.
#' @param grid Logit grid over which to evaluate the information.
#' @return A data frame with \code{theta}, \code{info}, and \code{sem}.
#' @examples
#' set.seed(1)
#' d <- seq(-1.5, 1.5, length.out = 6)
#' X <- matrix(rbinom(300 * 6, 1, plogis(outer(rnorm(300), d, "-"))), 300, 6)
#' colnames(X) <- paste0("I", 1:6)
#' head(test_information(rasch(X)))
#' @export
test_information <- function(fit, grid = seq(-6, 6, by = 0.1)) {
  L <- length(fit$tau_list)
  disc <- if (is.null(fit$disc)) rep(1, L) else fit$disc
  info <- vapply(grid, function(th)
    sum(vapply(seq_len(L), function(i)
      disc[i]^2 * item_moments(th, fit$tau_list[[i]], disc = disc[i])$V, 0)), 0)
  data.frame(theta = grid, info = info, sem = 1 / sqrt(info))
}
