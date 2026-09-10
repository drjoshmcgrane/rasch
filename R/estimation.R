# rasch :: estimation
# ===========================================================================
# Pairwise conditional maximum likelihood after Zwinderman (1995); the
# principal-components parameterisation in pcml_pc follows Andrich & Luo
# (2003). For items i, j with maximum scores m_i,
# m_j, the distribution of X_i given the pair total X_i + X_j = r is free of
# the person parameter:
#
#   P(X_i = k | X_i + X_j = r) = exp(-L_i(k) - L_j(r-k)) / sum_k' exp(...)
#
# where L_i(k) = sum_{h<=k} tau_ih is the cumulative threshold sum. The
# pairwise conditional log-likelihood, summed over all item pairs and pair
# totals, is maximised by Newton-Raphson. Standard errors use the Godambe
# sandwich covariance of the pairwise composite likelihood. The rating scale model is
# the same likelihood under the constraint tau_ik = delta_i + kappa_k,
# imposed through the design matrix. Dichotomous data is the special case
# m_i = 1. Australian English; no em dashes by house style.
# ===========================================================================


# Adjust the finite probabilities while retaining every member of the
# declared family. stats::p.adjust() otherwise reduces n when p contains NA,
# which makes an incomplete family less stringent than the question asked.
.p_adjust_family <- function(p, method = "holm", n = length(p)) {
  if (length(n) != 1L || !is.numeric(n) || is.complex(n) ||
      !is.null(dim(n)) || !is.null(oldClass(n)) || is.na(n) ||
      !is.finite(n) || n < length(p) || n > .Machine$integer.max ||
      n != floor(n))
    stop("`n` must be one whole number at least as large as `p`",
         call. = FALSE)
  out <- rep(NA_real_, length(p))
  usable <- is.finite(p)
  if (any(usable))
    out[usable] <- stats::p.adjust(p[usable], method = method,
                                   n = as.integer(n))
  out
}

# A zero estimated standard error is not evidence of an infinitely precise
# effect. It usually marks a singular sandwich, constant bootstrap column or
# fixed parameter. Withhold the Wald statistic instead of manufacturing an
# infinite value and p = 0. Any positive finite standard error remains
# testable: comparing it with the coefficient would make the answer depend on
# the units in which a predictor or contrast happened to be expressed.
.wald_ratio <- function(estimate, se) {
  n <- max(length(estimate), length(se))
  if (!n) return(numeric(0))
  if (!length(estimate) || !length(se) ||
      !length(estimate) %in% c(1L, n) || !length(se) %in% c(1L, n))
    stop("Wald estimates and standard errors have incompatible lengths",
         call. = FALSE)
  out <- rep(NA_real_, n)
  estimate <- rep_len(estimate, length(out))
  se <- rep_len(se, length(out))
  usable <- is.finite(estimate) & is.finite(se) & se > 0
  out[usable] <- estimate[usable] / se[usable]
  out
}

# Relative inverse-variance weights. Dividing by the smallest positive
# variance before inversion gives the same weighted mean as 1 / v without an
# arbitrary absolute floor or overflow. An exactly zero variance is a fixed
# contrast and therefore carries all the weight; several such contrasts are
# weighted equally.
.inverse_variance_weights <- function(v) {
  if (!is.numeric(v) || is.complex(v) || !length(v) ||
      any(!is.finite(v)) || any(v < 0))
    stop("variances must be finite non-negative numeric values", call. = FALSE)
  if (any(v == 0)) return(as.numeric(v == 0))
  rel <- v / min(v)
  1 / rel
}

# A bank's stated marginal standard errors and attached covariance must refer
# to the same calibration. Use relative agreement at their own scale: an
# absolute floor would accept materially different uncertainties merely
# because both happen to be expressed in small units.
.se_covariance_agree <- function(se, cov_se, tolerance = 1e-6) {
  if (length(se) != length(cov_se)) return(FALSE)
  scale <- pmax(abs(se), abs(cov_se))
  abs(se - cov_se) <= tolerance * scale
}

# Positive-semidefiniteness is also scale free. Covariance matrices may be
# singular because of identification constraints, so retain eigenvalues that
# are negative only at the stated relative numerical tolerance.
.covariance_is_psd <- function(C, tolerance = 1e-8) {
  if (!is.matrix(C) || !is.numeric(C) || is.complex(C) ||
      nrow(C) != ncol(C) ||
      any(!is.finite(C))) return(FALSE)
  if (!nrow(C)) return(TRUE)
  cscale <- max(abs(C))
  if (cscale > 0 && max(abs(C - t(C))) > tolerance * cscale) return(FALSE)
  ev <- eigen((C + t(C)) / 2, symmetric = TRUE, only.values = TRUE)$values
  scale <- max(abs(ev))
  !length(ev) || scale == 0 || min(ev) >= -tolerance * scale
}

.covariance_is_symmetric <- function(C, tolerance = 1e-8) {
  if (!is.matrix(C) || !is.numeric(C) || is.complex(C) ||
      nrow(C) != ncol(C) ||
      any(!is.finite(C))) return(FALSE)
  if (!nrow(C)) return(TRUE)
  scale <- max(abs(C))
  scale == 0 || max(abs(C - t(C))) <= tolerance * scale
}

# Local-maximum check in free likelihood coordinates. Expected information
# alone misses the residual-weighted curvature of a nonlinear parameter map.
# Standardise the diagonal before testing definiteness so changing a
# parameter's units does not change the decision. Flat directions do not
# establish an identified maximum, even when the score is zero.
.likelihood_curvature_ok <- function(H, tolerance = 1e-10) {
  if (!.covariance_is_symmetric(H)) return(FALSE)
  if (!nrow(H)) return(TRUE)
  if (any(diag(H) >= 0)) return(FALSE)
  sc <- sqrt(-diag(H))
  information <- sweep(sweep(-H, 1L, sc, "/"), 2L, sc, "/")
  if (any(!is.finite(information))) return(FALSE)
  ev <- eigen((information + t(information)) / 2,
              symmetric = TRUE, only.values = TRUE)$values
  min(ev) > tolerance * max(abs(ev))
}

# A covariance used for Wald inference must be more than present. This common
# gate prevents a malformed or materially indefinite matrix from becoming a
# zero variance through downstream pmax(..., 0) guards. Singular matrices are
# allowed: identifying constraints commonly make a full parameter covariance
# positive semidefinite rather than positive definite.
.covariance_supports_wald <- function(C, n = NULL) {
  if (!is.matrix(C) || !is.numeric(C) || is.complex(C) ||
      nrow(C) != ncol(C) ||
      any(!is.finite(C))) return(FALSE)
  if (!is.null(n) && !identical(dim(C), c(as.integer(n), as.integer(n))))
    return(FALSE)
  .covariance_is_symmetric(C) && .covariance_is_psd(C)
}

# First-order Kent calibration for a composite likelihood-ratio statistic.
# The calibration is inferential, so an indefinite or singular estimated
# sensitivity in the tested directions makes it unavailable. In particular,
# do not let a generic solve() error escape or turn numerical eigenvalues into
# a plausible-looking probability.
.kent_calibration <- function(statistic, C, covariance,
                              sensitivity_inverse,
                              tolerance = sqrt(.Machine$double.eps)) {
  unavailable <- list(chisq = NA_real_, p = NA_real_, lambda = numeric(0))
  if (length(statistic) != 1L || !is.numeric(statistic) ||
      is.complex(statistic) || !is.null(dim(statistic)) ||
      !is.finite(statistic) || statistic < 0 || !is.matrix(C) ||
      !is.numeric(C) || is.complex(C) || !ncol(C) || any(!is.finite(C)) ||
      !is.matrix(covariance) || !is.numeric(covariance) ||
      is.complex(covariance) || !is.matrix(sensitivity_inverse) ||
      !is.numeric(sensitivity_inverse) || is.complex(sensitivity_inverse))
    return(unavailable)
  p <- nrow(C)
  if (!identical(dim(covariance), c(p, p)) ||
      !identical(dim(sensitivity_inverse), c(p, p)) ||
      any(!is.finite(covariance)) || any(!is.finite(sensitivity_inverse)) ||
      !.covariance_supports_wald(covariance, p) ||
      !.covariance_supports_wald(sensitivity_inverse, p))
    return(unavailable)
  num <- crossprod(C, covariance %*% C)
  den <- crossprod(C, sensitivity_inverse %*% C)
  num <- (num + t(num)) / 2
  den <- (den + t(den)) / 2
  ed <- tryCatch(eigen(den, symmetric = TRUE), error = function(e) NULL)
  if (is.null(ed) || any(!is.finite(ed$values))) return(unavailable)
  dscale <- max(abs(ed$values))
  if (!is.finite(dscale) || dscale <= 0 ||
      min(ed$values) <= tolerance * dscale)
    return(unavailable)
  inv_root <- ed$vectors %*%
    (diag(1 / sqrt(ed$values), nrow = length(ed$values))) %*%
    t(ed$vectors)
  ratio <- inv_root %*% num %*% inv_root
  lambda <- tryCatch(eigen((ratio + t(ratio)) / 2, symmetric = TRUE,
                           only.values = TRUE)$values,
                     error = function(e) NULL)
  if (is.null(lambda) || any(!is.finite(lambda))) return(unavailable)
  lscale <- max(abs(lambda))
  if (!is.finite(lscale) || lscale <= 0 ||
      min(lambda) <= tolerance * lscale)
    return(unavailable)
  adjusted <- statistic * length(lambda) / sum(lambda)
  if (!is.finite(adjusted) || adjusted < 0) return(unavailable)
  list(chisq = adjusted,
       p = stats::pchisq(adjusted, length(lambda), lower.tail = FALSE),
       lambda = lambda)
}


# The columns the person table generates for itself. A factor sharing one of
# these names would be bound into the same table and silently replace the
# calculated column -- a factor called class_interval would stand in for the
# fitted intervals every fit statistic is computed over.
.person_reserved <- c("id", "n_items", "raw", "max_raw", "weighted_score",
                      "max_weighted_score",
                      "theta", "se", "extreme", "infit_ms", "outfit_ms",
                      "infit_z", "outfit_z", "fit_resid", "natural_resid",
                      "df_fit", "class_interval")

# A person-factor frame must carry unique, non-empty column names that no
# generated column already claims; an ambiguous factor structure poisons
# every downstream DIF family.
.check_factor_frame <- function(fac_df) {
  if (is.null(fac_df)) return(invisible(NULL))
  nms <- names(fac_df)
  if (is.null(nms) || any(is.na(nms)) || any(!nzchar(trimws(nms))))
    stop("every person factor needs a non-empty name (not whitespace-only)",
         call. = FALSE)
  if (anyDuplicated(nms))
    stop("duplicate factor column name(s): ",
         paste(unique(nms[duplicated(nms)]), collapse = ", "), call. = FALSE)
  clash <- intersect(nms, .person_reserved)
  if (length(clash))
    stop("person factor name(s) reserved for the fitted person table: ",
         paste(clash, collapse = ", "),
         "; rename them, or the calculated column would be replaced",
         call. = FALSE)
  invisible(NULL)
}

# One whole finite number within a range; NULL allowed only when said so.
.check_whole <- function(x, name, min = 1, max = Inf, null_ok = FALSE) {
  if (is.null(x)) {
    if (null_ok) return(invisible(NULL))
    stop("`", name, "` must be supplied", call. = FALSE)
  }
  upper <- min(max, .Machine$integer.max)
  if (length(x) != 1L || !is.numeric(x) || is.complex(x) ||
      !is.null(dim(x)) || !is.null(oldClass(x)) || !is.finite(x) ||
      x != floor(x) || x < min || x > upper)
    stop("`", name, "` must be one whole number between ", min, " and ",
         if (is.finite(max)) max else "the integer range", call. = FALSE)
  invisible(as.integer(x))
}

# One probability strictly inside (0, 1).
.check_prob <- function(x, name) {
  if (length(x) != 1L || !is.numeric(x) || is.complex(x) ||
      !is.null(dim(x)) || !is.null(oldClass(x)) || !is.finite(x) ||
      x <= 0 || x >= 1)
    stop("`", name, "` must be one probability strictly between 0 and 1",
         call. = FALSE)
  invisible(x)
}

# Shared validation for the Newton controls every estimator accepts.
.check_controls <- function(maxit, tol) {
  if (length(maxit) != 1L || !is.numeric(maxit) || is.complex(maxit) ||
      !is.null(dim(maxit)) || !is.null(oldClass(maxit)) ||
      !is.finite(maxit) ||
      maxit != floor(maxit) || maxit < 1 || maxit > .Machine$integer.max)
    stop("`maxit` must be one whole positive iteration cap", call. = FALSE)
  if (length(tol) != 1L || !is.numeric(tol) || is.complex(tol) ||
      !is.null(dim(tol)) || !is.null(oldClass(tol)) || !is.finite(tol) ||
      tol <= 0)
    stop("`tol` must be one positive finite tolerance", call. = FALSE)
  invisible(NULL)
}

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
  if (!is.numeric(m) || is.complex(m) || !is.null(dim(m)) ||
      !is.null(oldClass(m)) || !length(m) || any(!is.finite(m)) ||
      any(m != floor(m)) || any(m < 0) || any(m > .Machine$integer.max))
    stop("`m` must hold at least one whole non-negative maximum score")
  thr <- do.call(rbind, lapply(seq_along(m), function(i)
    if (m[i] >= 1) data.frame(item = i, k = seq_len(m[i])) else NULL))
  if (is.null(thr))
    return(data.frame(id = integer(0), item = integer(0), k = integer(0)))
  thr$id <- seq_len(nrow(thr)); thr[, c("id", "item", "k")]
}

# Cross-tabulated pair counts: for every item pair i < j, the matrix of
# joint category counts over persons observed on both items.
.pair_counts <- function(X, m) {
  L <- ncol(X); out <- vector("list", 0L)
  for (i in seq_len(L - 1)) for (j in (i + 1):L) {
    both <- !is.na(X[, i]) & !is.na(X[, j])
    if (!any(both)) next
    idx <- X[both, i] * (m[j] + 1L) + X[both, j] + 1L
    n <- matrix(tabulate(idx, nbins = (m[i] + 1L) * (m[j] + 1L)),
                nrow = m[i] + 1L, byrow = TRUE)
    out[[length(out) + 1L]] <- list(i = i, j = j, n = n)
  }
  out
}

# ---------------------------------------------------------------------------
# Structural identification checks, run before solving. The pairwise
# conditional likelihood factorises over item pairs, so relative locations
# between two blocks of items are identified only if some person answered
# items in both (the item-pair graph is connected); on a disconnected design
# the likelihood is flat in the between-block shift and Newton lands wherever
# the ridge sends it -- garbage that must be an error, not a result. Anchors
# rescue a block: a block containing an anchored item has its origin fixed.
# ---------------------------------------------------------------------------
# Refuse fractional scores BEFORE integer coercion truncates them: every
# public estimator must give the same answer as rasch() here (which already
# rejects), not silently fit altered data.
.check_integer_scores <- function(X, what = "X") {
  # factors must be read through their LABELS: as.numeric(factor) returns
  # level codes, which let factor("1.9") slip past as the integer 3
  xc <- as.character(if (is.factor(X)) as.character(X) else X)
  obs <- !is.na(xc)
  Xn <- suppressWarnings(as.numeric(xc))
  nonnum <- obs & is.na(Xn)
  if (any(nonnum))
    stop("non-numeric score(s) in ", what, " (e.g. '", xc[nonnum][1],
         "'); scores must be integer counts", call. = FALSE)
  noninf <- obs & !is.finite(Xn)
  if (any(noninf))
    stop("non-finite score(s) in ", what, " (e.g. ", xc[noninf][1],
         "); scores must be integer counts", call. = FALSE)
  outside <- obs & is.finite(Xn) &
    (Xn > .Machine$integer.max | Xn < -.Machine$integer.max)
  if (any(outside))
    stop("score(s) outside the supported integer range in ", what,
         " (e.g. ", format(Xn[outside][1]),
         "); rescore the response categories before analysis", call. = FALSE)
  bad <- obs & Xn != round(Xn)
  if (any(bad))
    stop("non-integer score(s) in ", what, " (e.g. ",
         format(Xn[bad][1]), "); Rasch categories are integer counts -- ",
         "round or rescore explicitly before analysis", call. = FALSE)
  invisible(TRUE)
}

# The low-level estimators do not perform rasch()'s category preparation.
# Their pair tables index scores directly as 0, ..., m, so negative values are
# silently omitted by tabulate(), a gap creates a threshold for a category that
# was never observed, and a constant/all-missing column has no estimable item
# parameter. Refuse those inputs before integer storage or pair construction.
.pcml_score_matrix <- function(X) {
  X <- as.matrix(X)
  if (length(dim(X)) != 2L || nrow(X) < 1L || ncol(X) < 2L)
    stop("`X` must contain at least one person and two item columns",
         call. = FALSE)
  nm <- colnames(X)
  if (!is.null(nm) && (anyNA(nm) || any(!nzchar(trimws(nm)))))
    stop("item column names must be non-missing and non-empty (not whitespace-only)",
         call. = FALSE)
  .check_integer_scores(X, "the score matrix")
  storage.mode(X) <- "integer"
  for (j in seq_len(ncol(X))) {
    z <- sort(unique(X[!is.na(X[, j]), j]))
    lab <- if (is.null(colnames(X))) paste0("column ", j) else
      paste0("item ", colnames(X)[j])
    if (!length(z))
      stop(lab, " has no observed scores", call. = FALSE)
    if (length(z) < 2L)
      stop(lab, " is constant; item parameters require at least two observed scores",
           call. = FALSE)
    if (!identical(z, seq.int(0L, max(z))))
      stop(lab, " must use consecutive integer categories from 0; observed: ",
           paste(z, collapse = ", "), call. = FALSE)
  }
  X
}

.pcml_check_connected <- function(pairs, L, item_names, anchored = integer(0)) {
  # a co-observed pair whose every observed total is 0 or the maximum has a
  # single feasible conditional allocation and carries NO information: it
  # must not count as a link, or one respondent scoring (0, 0) across two
  # blocks would "connect" them while the likelihood stays flat between
  # them (every intermediate total has at least two allocations, so any
  # response off the two extreme-total corners is a real link)
  informative <- vapply(pairs, function(p) {
    sum(p$n) - p$n[1L, 1L] - p$n[nrow(p$n), ncol(p$n)] > 0
  }, TRUE)
  pairs <- pairs[informative]
  edges <- if (length(pairs))
    do.call(rbind, lapply(pairs, function(p) c(p$i, p$j)))
  else matrix(integer(0), 0L, 2L)
  comp <- .btlef_components(L, edges)
  if (length(unique(comp)) == 1L) return(invisible(comp))
  if (length(anchored)) {
    bad <- setdiff(unique(comp), unique(comp[anchored]))
    if (!length(bad)) return(invisible(comp))
    blocks <- vapply(bad, function(cc)
      paste(item_names[comp == cc], collapse = ", "), "")
    stop("the item-pair graph is not connected, and block(s) without an ",
         "anchored item have no identified origin: ",
         paste0("{", blocks, "}", collapse = " "),
         "; link the blocks through common persons or anchor an item in ",
         "every block", call. = FALSE)
  }
  blocks <- tapply(item_names, comp, paste, collapse = ", ")
  stop("the item-pair graph is not connected: no person answered items in ",
       "more than one of the blocks ",
       paste0("{", blocks, "}", collapse = " | "),
       "; relative locations between the blocks are unidentified -- link ",
       "them through common items or persons, or anchor item(s) in every ",
       "block", call. = FALSE)
}

# A threshold k is informed by the persons observed in its two adjacent
# categories; when either count is tiny the conditional estimate can run
# away (a category with one response sends its threshold toward the
# boundary) while the ridged covariance reports a spuriously small standard
# error. Flag such thresholds so the caller can report the estimate with an
# NA standard error and a note naming the cause -- an honest answer, not a
# manufactured one. Categories with zero responses are the caller's problem
# (rasch() rescores them away); the danger zone handled here is 1-2.
.pcml_weak_thresholds <- function(X, m, thr, item_names, min_count = 3L,
                                  min_item_count = 8L) {
  flag <- logical(nrow(thr)); notes <- character(0)
  for (i in seq_len(ncol(X))) {
    cnt <- tabulate(X[, i] + 1L, nbins = m[i] + 1L)
    weak_k <- which(pmin(cnt[-length(cnt)], cnt[-1]) < min_count)
    if (length(weak_k)) {
      flag[thr$item == i & thr$k %in% weak_k] <- TRUE
      kc <- which(cnt < min_count) - 1L
      notes <- c(notes, sprintf(
        "item %s: only %s response(s) in category %s; threshold(s) %s and the item location are weakly determined (SE reported as NA) -- consider pc_components or collapsing categories",
        item_names[i], paste(cnt[kc + 1L], collapse = "/"),
        paste(kc, collapse = "/"), paste(weak_k, collapse = "/")))
    }
    # an item's thresholds are estimated JOINTLY, so a critically sparse
    # category destabilises its siblings too, not only the adjacent
    # threshold: in simulation, a ~4-response category left a sibling
    # threshold's reported SE understated four-fold while its own local
    # counts looked healthy (~7 responses gave ~1.7x). Flag the whole item
    # once any category falls below min_item_count.
    kc_it <- which(cnt < min_item_count) - 1L
    if (length(kc_it) && any(!flag[thr$item == i])) {
      already <- all(flag[thr$item == i])
      flag[thr$item == i] <- TRUE
      if (!already) notes <- c(notes, sprintf(
        "item %s: category %s has only %s response(s); all of the item's jointly estimated thresholds are unreliable at this sparsity (SEs reported as NA) -- consider pc_components or collapsing categories",
        item_names[i], paste(kc_it, collapse = "/"),
        paste(cnt[kc_it + 1L], collapse = "/")))
    }
  }
  list(flag = flag, notes = notes)
}

# ---------------------------------------------------------------------------
# Starting values: weighted least squares on the pairwise log-ratios. Used
# only to seed Newton-Raphson; the returned estimates always come from the
# conditional likelihood itself.
# ---------------------------------------------------------------------------
.start_tau <- function(X, thr, cont = 0.5) {
  M <- nrow(thr)
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
        if (cA == 0 || cB == 0) { a <- cA + cont; b <- cB + cont }
        D[p, q] <- log(a / b)
        W[p, q] <- (a * b) / (a + b)
      }
    }
  }
  rows <- which(!is.na(D) & upper.tri(D), arr.ind = TRUE)
  if (!nrow(rows)) return(rep(0, M))
  C <- matrix(0, nrow(rows), M); d <- numeric(nrow(rows)); wt <- numeric(nrow(rows))
  for (r in seq_len(nrow(rows))) {
    p <- rows[r, 1]; q <- rows[r, 2]
    C[r, p] <- 1; C[r, q] <- -1; d[r] <- D[p, q]; wt[r] <- W[p, q]
  }
  sw <- sqrt(wt)
  tau <- c(qr.coef(qr((C[, -M, drop = FALSE]) * sw), d * sw), 0)
  tau[is.na(tau)] <- 0
  tau - mean(tau)
}

# Pseudo log-likelihood, gradient, and Hessian over the full threshold vector.
.pcml_glh <- function(tau, thr, pairs, m) {
  M <- nrow(thr)
  g <- numeric(M); H <- matrix(0, M, M); ll <- 0
  cum <- lapply(seq_along(m), function(i) cumsum(tau[thr$item == i]))
  ids <- lapply(seq_along(m), function(i) thr$id[thr$item == i])
  for (pc in pairs) {
    i <- pc$i; j <- pc$j; n <- pc$n
    Li <- c(0, cum[[i]]); Lj <- c(0, cum[[j]])
    idx <- c(ids[[i]], ids[[j]]); mi <- m[i]; mj <- m[j]
    for (r in seq_len(mi + mj - 1L)) {
      ks <- max(0L, r - mj):min(mi, r)
      if (length(ks) < 2L) next
      nk <- n[cbind(ks + 1L, r - ks + 1L)]
      N <- sum(nk)
      if (N == 0) next
      lp <- -(Li[ks + 1L] + Lj[r - ks + 1L])
      lp <- lp - max(lp); elp <- exp(lp); p <- elp / sum(elp)
      ll <- ll + sum(nk * (lp - log(sum(elp))))
      # local coefficients d lp_k / d tau, item i columns then item j columns
      U <- cbind(
        -outer(ks, seq_len(mi), ">="),
        -outer(r - ks, seq_len(mj), ">="))
      storage.mode(U) <- "double"
      g[idx] <- g[idx] + drop(crossprod(U, nk - N * p))
      Ep <- drop(crossprod(U, p))
      S  <- crossprod(U, p * U)
      H[idx, idx] <- H[idx, idx] - N * (S - tcrossprod(Ep))
    }
  }
  list(ll = ll, g = g, H = H)
}

# Convert response-row identifiers to cluster indices. Missing identifiers are
# unknown people, not one shared cluster, so each receives its own index.
.pcml_cluster_index <- function(cluster, N) {
  if (!is.atomic(cluster) || !is.null(dim(cluster)) || length(cluster) != N)
    stop("internal calibration clusters must give one plain identifier per response row",
         call. = FALSE)
  z <- .role_text_values(cluster)
  missing <- is.na(z) | !nzchar(z)
  known <- unique(z[!missing])
  group <- match(z, known)
  if (any(missing))
    group[missing] <- length(known) + seq_len(sum(missing))
  as.integer(group)
}

# Godambe sandwich covariance for the pairwise pseudo-likelihood. The naive
# inverse information overstates precision because every response enters
# L - 1 overlapping pairs; the sandwich H^-1 J H^-1 with J the empirical
# covariance of the per-person scores corrects this.
.pcml_sandwich <- function(X, thr, m, tau, pairs, cluster = NULL) {
  M <- nrow(thr); N <- nrow(X)
  cum <- lapply(seq_along(m), function(i) cumsum(tau[thr$item == i]))
  ids <- lapply(seq_along(m), function(i) thr$id[thr$item == i])
  S <- matrix(0, N, M)
  # Count only conditional item-pair observations whose total admits at
  # least two response allocations. This is design information, rather than
  # the realised score magnitude: an informative cluster may legitimately
  # have a score vector that cancels to zero at the fitted parameters.
  pair_load <- integer(N)
  for (pc in pairs) {
    i <- pc$i; j <- pc$j; mi <- m[i]; mj <- m[j]
    Li <- c(0, cum[[i]]); Lj <- c(0, cum[[j]])
    idx <- c(ids[[i]], ids[[j]])
    # per-cell score vector u_k - ubar_r, indexed by cell (k, l)
    V <- matrix(0, (mi + 1L) * (mj + 1L), mi + mj)
    for (r in seq_len(mi + mj - 1L)) {
      ks <- max(0L, r - mj):min(mi, r)
      if (length(ks) < 2L) next
      lp <- -(Li[ks + 1L] + Lj[r - ks + 1L])
      lp <- lp - max(lp); p <- exp(lp) / sum(exp(lp))
      U <- cbind(-outer(ks, seq_len(mi), ">="),
                 -outer(r - ks, seq_len(mj), ">="))
      storage.mode(U) <- "double"
      ub <- drop(crossprod(U, p))
      V[ks * (mj + 1L) + (r - ks) + 1L, ] <- sweep(U, 2, ub)
    }
    both <- which(!is.na(X[, i]) & !is.na(X[, j]))
    if (!length(both)) next
    total <- X[both, i] + X[both, j]
    informative <- pmin(mi, total) > pmax(0, total - mj)
    if (any(informative))
      pair_load[both[informative]] <- pair_load[both[informative]] + 1L
    cell <- X[both, i] * (mj + 1L) + X[both, j] + 1L
    S[both, idx] <- S[both, idx] + V[cell, , drop = FALSE]
  }
  # Rows are the independent units unless an ID joins repeated rows from the
  # same person. Build support for both cases: even without repeated IDs, a
  # sandwich over fewer independent people than fitted directions is singular
  # and cannot provide inferential covariance.
  group_all <- if (is.null(cluster)) seq_len(N) else
    .pcml_cluster_index(cluster, N)
  contributes <- pair_load > 0L
  group_raw <- group_all[contributes]
  if (length(group_raw)) {
    group <- match(group_raw, unique(group_raw))
    cluster_load <- as.numeric(rowsum(pair_load[contributes], group,
                                      reorder = FALSE))
    S <- rowsum(S[contributes, , drop = FALSE], group = group,
                reorder = FALSE)
  } else {
    cluster_load <- numeric(0)
    S <- matrix(0, 0L, ncol(S))
  }
  repeated <- anyDuplicated(group_raw) > 0L
  # With repeated rows, the empirical meat is formed from a finite number of
  # person clusters. Apply the same leading CR1 correction used for clustered
  # comparative judgements. The t reference used by explanatory coefficient
  # tests addresses the reference distribution; it does not replace this
  # degrees-of-freedom correction to the covariance itself.
  cr1 <- if (repeated && length(cluster_load) > 1L)
    length(cluster_load) / (length(cluster_load) - 1L) else 1
  cluster_support <- list(
    repeated = repeated,
    n = length(cluster_load),
    effective = if (length(cluster_load))
      sum(cluster_load)^2 / sum(cluster_load^2) else 0,
    informative_rows = sum(contributes),
    pair_contributions = sum(pair_load),
    cr1 = cr1
  )
  out <- cr1 * crossprod(S)
  attr(out, "cluster_support") <- cluster_support
  out
}

# Linearised delete-one-cluster covariance for a restricted PCML fit.  The
# ordinary cluster sandwich is first-order unbiased only as the number of
# clusters grows.  In a small repeated-person explanatory analysis, clusters
# with different numbers of response occasions can have appreciably different
# leverage.  For cluster g, solve (A - A_g)d_g = s_g, where A_g and s_g are its
# sensitivity and score contributions at the full-data estimate, then apply
# the usual delete-one-cluster jackknife covariance to the d_g.  This is the
# CR3/bias-reduced linearisation of an actual leave-one-person refit, without
# fitting the model G additional times.
.pcml_linearised_cluster_cov <- function(X, thr, m, tau, pairs, B, H_beta,
                                         cluster) {
  parameter_scale <- .design_column_scale(B)
  scale_outer <- outer(parameter_scale, parameter_scale)
  B <- sweep(B, 2L, parameter_scale, `/`)
  H_beta <- H_beta / scale_outer
  N <- nrow(X); P <- ncol(B)
  group_all <- .pcml_cluster_index(cluster, N)
  G_all <- max(group_all)
  S <- matrix(0, N, P)
  A_g <- array(0, c(G_all, P, P))
  pair_load <- integer(N)
  cum <- lapply(seq_along(m), function(i) cumsum(tau[thr$item == i]))
  ids <- lapply(seq_along(m), function(i) thr$id[thr$item == i])

  for (pc in pairs) {
    i <- pc$i; j <- pc$j; mi <- m[i]; mj <- m[j]
    Li <- c(0, cum[[i]]); Lj <- c(0, cum[[j]])
    idx <- c(ids[[i]], ids[[j]])
    B_pair <- B[idx, , drop = FALSE]
    score_cell <- matrix(0, (mi + 1L) * (mj + 1L), P)
    info_total <- vector("list", mi + mj + 1L)
    for (r in seq_len(mi + mj - 1L)) {
      ks <- max(0L, r - mj):min(mi, r)
      if (length(ks) < 2L) next
      lp <- -(Li[ks + 1L] + Lj[r - ks + 1L])
      lp <- lp - max(lp); prob <- exp(lp) / sum(exp(lp))
      U <- cbind(-outer(ks, seq_len(mi), ">="),
                 -outer(r - ks, seq_len(mj), ">="))
      storage.mode(U) <- "double"
      Ub <- U %*% B_pair
      mean_u <- drop(crossprod(Ub, prob))
      centred <- sweep(Ub, 2L, mean_u)
      score_cell[ks * (mj + 1L) + (r - ks) + 1L, ] <- centred
      info_total[[r + 1L]] <- crossprod(centred, prob * centred)
    }
    both <- which(!is.na(X[, i]) & !is.na(X[, j]))
    if (!length(both)) next
    total <- X[both, i] + X[both, j]
    informative <- pmin(mi, total) > pmax(0, total - mj)
    if (any(informative)) {
      rows <- both[informative]
      pair_load[rows] <- pair_load[rows] + 1L
      for (h in which(informative)) {
        g <- group_all[both[h]]
        A_g[g, , ] <- A_g[g, , ] + info_total[[total[h] + 1L]]
      }
    }
    cell <- X[both, i] * (mj + 1L) + X[both, j] + 1L
    S[both, ] <- S[both, , drop = FALSE] +
      score_cell[cell, , drop = FALSE]
  }

  contributes <- pair_load > 0L
  used <- unique(group_all[contributes])
  G <- length(used)
  if (G < 2L) return(NULL)
  group <- match(group_all[contributes], used)
  S_g <- rowsum(S[contributes, , drop = FALSE], group, reorder = FALSE)
  A_g <- A_g[used, , , drop = FALSE]
  A <- apply(A_g, c(2L, 3L), sum)
  if (!is.matrix(A)) A <- matrix(A, P, P)
  target <- -H_beta
  scale <- max(1, max(abs(target)))
  if (any(!is.finite(A)) || max(abs(A - target)) > 1e-7 * scale)
    stop("internal cluster sensitivity decomposition did not reproduce the fitted information",
         call. = FALSE)

  delta <- matrix(NA_real_, G, P)
  for (g in seq_len(G)) {
    Ag <- matrix(A_g[g, , , drop = FALSE], P, P)
    keep_A <- A - Ag
    rc <- tryCatch(rcond(keep_A), error = function(e) 0)
    if (!is.finite(rc) || rc <= 1e-12) return(NULL)
    delta[g, ] <- tryCatch(solve(keep_A, S_g[g, ]),
                            error = function(e) rep(NA_real_, P))
  }
  if (any(!is.finite(delta))) return(NULL)
  delta <- sweep(delta, 2L, colMeans(delta))
  (G - 1) / G * crossprod(delta) / scale_outer
}

# Decide whether an empirical composite-score meat can support every fitted
# parameter direction. The same rules apply wherever the PCML score sandwich
# is used: ordinary calibrations, explanatory restrictions and EFRM stage one.
# The point estimate comes from the conditional likelihood and remains useful
# when this guard withholds its empirical covariance.
.pcml_covariance_support <- function(Jb, p, support) {
  if (is.null(support))
    return(list(inference = TRUE, note = NULL, rank_deficient = FALSE))
  rank_j <- tryCatch(qr(Jb, tol = 1e-10)$rank,
                     error = function(e) NA_integer_)
  rank_deficient <- !is.finite(rank_j) || rank_j < p
  inference <- is.numeric(support$n) && length(support$n) == 1L &&
    is.finite(support$n) && is.numeric(support$effective) &&
    length(support$effective) == 1L && is.finite(support$effective) &&
    support$n >= 10L && support$effective >= 8 &&
    support$effective > p && !rank_deficient
  if (inference)
    return(list(inference = TRUE, note = NULL,
                rank_deficient = FALSE))

  unit <- if (isTRUE(support$repeated))
    "person clusters" else "independent person units"
  reason <- if (!is.numeric(support$n) || length(support$n) != 1L ||
      !is.finite(support$n) || support$n < 10L)
    paste("fewer than 10", unit)
  else if (rank_deficient)
    paste("a rank-deficient", if (isTRUE(support$repeated))
      "person-cluster" else "person-unit", "score covariance")
  else if (!is.numeric(support$effective) ||
      length(support$effective) != 1L || !is.finite(support$effective) ||
      support$effective < 8)
    paste("fewer than 8 effective", unit)
  else
    paste0("the effective count of ", unit,
           " does not exceed the number of fitted parameters")
  note <- sprintf(paste0(
    "item-parameter uncertainty withheld: %d %s (%.1f effective) ",
    "support %d fitted parameters, with %s; point estimates remain available"),
    support$n, unit, support$effective, p, reason)
  list(inference = FALSE, note = note, rank_deficient = rank_deficient)
}

# A PCM threshold next to a category that no informative response pattern
# observes has no finite pairwise conditional estimate: the solver runs it
# towards infinity and the projected information goes singular, which the
# rank check then reports as a connectivity problem. Name the category
# instead. A threshold fixed by an anchor needs no data, so an item whose
# thresholds are all fixed passes.
.pcml_check_uninformative <- function(X, m, inames, fixed = integer(0),
                                      thr) {
  cc <- .conditional_categories(X)
  bad <- character(0)
  for (i in seq_along(m)) {
    if (all(thr$id[thr$item == i] %in% fixed)) next
    lost <- setdiff(seq.int(0L, m[i]), cc[[i]])
    if (length(lost))
      bad <- c(bad, sprintf("%s (categor%s %s)", inames[i],
                            if (length(lost) > 1L) "ies" else "y",
                            paste(lost, collapse = ", ")))
  }
  if (length(bad))
    stop("categories observed only in extreme response patterns (every ",
         "answered item at its minimum or at its maximum) carry no pairwise ",
         "conditional information, so their thresholds cannot be estimated: ",
         paste(bad, collapse = "; "),
         "; collapse those categories before calling pcml() -- rasch() ",
         "does so automatically", call. = FALSE)
  invisible(NULL)
}

# Bound design columns without squaring their entries to obtain the scale.
.design_column_scale <- function(B) {
  z <- apply(abs(B), 2L, max)
  if (any(!is.finite(z))) stop("the design contains non-finite values")
  z[z == 0] <- 1
  if (any(!is.finite(z^2) | z^2 == 0))
    stop("the predictor units are outside the representable covariance range; ",
         "rescale the predictors", call. = FALSE)
  z
}

# Newton-Raphson on tau = offset + B beta, where B removes the location
# indeterminacy, imposes the rating scale or facet structure, or restricts
# estimation to the unanchored thresholds (offset carrying the anchors).
.pcml_solve <- function(X, thr, m, B, beta0, offset = 0, maxit = 60, tol = 1e-8,
                        pairs = NULL, cluster = NULL) {
  if (is.null(pairs)) pairs <- .pair_counts(X, m)
  if (!length(pairs)) stop("no informative item pairs: check the data")
  # Work in bounded design columns. Predictor units must not determine the
  # Newton stopping rule, information rank or covariance support. Restore the
  # caller's coefficient coordinates, including both covariance axes, below.
  parameter_scale <- .design_column_scale(B)
  B <- sweep(B, 2L, parameter_scale, `/`)
  beta <- beta0 * parameter_scale
  glh <- .pcml_glh(drop(offset + B %*% beta), thr, pairs, m)
  it <- 0L
  for (it in seq_len(maxit)) {
    gb <- drop(crossprod(B, glh$g))
    Hb <- crossprod(B, glh$H %*% B)
    step <- tryCatch(solve(Hb, gb), error = function(e)
      solve(Hb - diag(1e-8, nrow(Hb)), gb))
    # step halving on the pseudo-likelihood
    lam <- 1; ok <- FALSE; g2 <- glh
    for (half in 1:30) {
      cand <- beta - lam * step
      g2 <- .pcml_glh(drop(offset + B %*% cand), thr, pairs, m)
      if (is.finite(g2$ll) && g2$ll >= glh$ll - 1e-12) { ok <- TRUE; break }
      lam <- lam / 2
    }
    if (!ok) break
    done <- max(abs(B %*% (lam * step))) < tol
    beta <- cand; glh <- g2
    if (done) break
  }
  Hb <- crossprod(B, glh$H %*% B)
  # the projected information must have full rank at the solution: a
  # singular Hb means some parameter direction is unidentified (e.g.
  # blocks linked only through non-informative extreme-total pairs slip
  # past graph checks), and the ridged inverse would report plausible or
  # even zero standard errors for a direction the data never determined
  rc <- tryCatch(rcond(Hb), error = function(e) 0)
  if (!(is.finite(rc) && rc > 1e-12))
    stop("the projected information matrix is singular (reciprocal ",
         "condition number ", format(rc, digits = 3), "): some parameter ",
         "direction is not identified by the data -- typically blocks of ",
         "items linked only through responses with no conditional ",
         "information (all at the minimum or maximum)", call. = FALSE)
  Hinv <- tryCatch(solve(Hb), error = function(e)
    solve(Hb - diag(1e-8, nrow(Hb))))
  gb_final <- drop(crossprod(B, glh$g))
  # The projected score is an extensive quantity and therefore grows with
  # sample size. At large N it can remain just above a fixed absolute cutoff
  # after the parameter estimates and log likelihood have stopped changing.
  # Accept either a small score in the bounded design or a small full Newton
  # move on the threshold scale, independent of predictor units. The information
  # rank check above prevents a small move in an unidentified direction from
  # being mistaken for convergence. Cap that allowance so an extremely loose
  # user tolerance cannot certify a visibly unfinished fit.
  newton_move <- drop(Hinv %*% gb_final)
  move_tol <- min(20 * tol, 1e-6)
  converged <- max(abs(gb_final)) < 1e-4 ||
    max(abs(B %*% newton_move)) < move_tol
  J  <- .pcml_sandwich(X, thr, m, drop(offset + B %*% beta), pairs,
                       cluster = cluster)
  Jb <- crossprod(B, J %*% B)
  covb <- Hinv %*% Jb %*% Hinv
  covt <- B %*% covb %*% t(B)
  cluster_support <- attr(J, "cluster_support", exact = TRUE)
  support_guard <- .pcml_covariance_support(Jb, ncol(B), cluster_support)
  cluster_inference <- support_guard$inference
  cluster_note <- support_guard$note
  if (!cluster_inference) {
      # A cluster sandwich with insufficient independent support can collapse
      # to zero (one cluster) or be rank-deficient while still looking
      # numerical. Keep exact fixed thresholds at zero and withhold every
      # estimated covariance entry instead of presenting false precision.
      fixed <- rowSums(abs(B)) == 0
      covb[,] <- NA_real_
      covt[,] <- NA_real_
      if (any(fixed)) {
        covt[fixed, ] <- 0
        covt[, fixed] <- 0
      }
  }
  scale_outer <- outer(parameter_scale, parameter_scale)
  list(tau = drop(offset + B %*% beta), beta = beta / parameter_scale,
       cov_beta = covb / scale_outer,
       cov_tau = covt, se_tau = sqrt(pmax(diag(covt), 0)),
       H_beta = Hb * scale_outer,
       loglik = glh$ll, iterations = it,
       converged = converged, cluster_inference = cluster_inference,
       cluster_support = cluster_support, cluster_note = cluster_note)
}

.pcml_anchor_k <- function(k) {
  if (!is.null(dim(k)) || !is.null(oldClass(k)))
    stop("anchor `k` must be a plain vector", call. = FALSE)
  if (is.logical(k) && all(is.na(k))) return(as.numeric(k))
  if (!is.numeric(k) || any(is.nan(k)) || any(!is.na(k) &
      (!is.finite(k) | k != floor(k) | k < 1)))
    stop("anchor `k` values must be positive whole threshold numbers or NA for item locations",
         call. = FALSE)
  as.numeric(k)
}

.pcml_anchor_columns <- function(anchors) {
  n <- nrow(anchors)
  item <- anchors$item
  item_kind <- is.character(item) || is.factor(item) ||
    (is.numeric(item) && !is.complex(item) && is.null(oldClass(item)))
  if (!item_kind || !is.null(dim(item)) || length(item) != n)
    stop("anchor `item` must be one plain item name or index per row",
         call. = FALSE)
  if (is.character(item) || is.factor(item)) {
    text <- as.character(item)
    if (anyNA(text) || any(!nzchar(trimws(text))))
      stop("anchor item names must be non-missing and non-blank",
           call. = FALSE)
  } else if (any(!is.finite(item)) || any(item != floor(item))) {
    stop("numeric anchor item indices must be finite whole numbers",
         call. = FALSE)
  }

  anchors$k <- .pcml_anchor_k(anchors$k)
  if (length(anchors$k) != n)
    stop("anchor `k` must contain one value per row", call. = FALSE)

  tau <- anchors$tau
  if (!is.numeric(tau) || is.complex(tau) || !is.null(dim(tau)) ||
      !is.null(oldClass(tau)) || length(tau) != n || any(!is.finite(tau)))
    stop("anchor `tau` must contain one plain finite numeric value per row",
         call. = FALSE)

  if ("average" %in% names(anchors)) {
    average <- anchors$average
    if (!is.logical(average) || !is.null(dim(average)) ||
        !is.null(oldClass(average)) || length(average) != n || anyNA(average))
      stop("the anchors `average` column must contain one plain TRUE or FALSE per row",
           call. = FALSE)
  }
  anchors
}

#' Estimate Rasch thresholds by pairwise conditional maximum likelihood
#'
#' Estimates PCM or RSM thresholds by Newton--Raphson maximisation of the
#' pairwise conditional likelihood (Zwinderman 1995).
#'
#' @details
#' For the PCM, the adjacent-category log odds are
#' \deqn{\log\{P(X_{ni}=k)/P(X_{ni}=k-1)\}=\theta_n-\delta_{ik}.}
#' Conditioning on the score for an item pair removes \eqn{\theta_n}. The PCM
#' estimates each \eqn{\delta_{ik}}; the RSM imposes
#' \eqn{\delta_{ik}=\beta_i+\tau_k} through a design matrix.
#'
#' @param X Persons-by-items integer score matrix. Each item must have at least
#'   two observed categories, numbered consecutively from 0. Missing values are
#'   handled by pairwise deletion, so linked booklet designs and
#'   random missingness estimate without imputation; the item-pair graph must
#'   be connected (some person answering items in both of any two blocks),
#'   otherwise relative locations between blocks are unidentified and the fit
#'   stops with an error naming the blocks -- unless \code{anchors} fix an
#'   item in every block, the disjoint-form equating case.
#' @param model \code{"PCM"} or \code{"RSM"}.
#' @param anchors Optional anchor table for equating: a data frame with
#'   columns \code{item} (name or column index), \code{k}, and \code{tau}
#'   (the anchor value). A numeric \code{k} fixes that single threshold
#'   (individual anchoring); \code{k = NA} fixes the item's mean location at
#'   \code{tau} while its thresholds remain free (location anchoring). The
#'   remaining parameters are estimated on the anchored scale and no
#'   recentring is applied. An optional logical column \code{average},
#'   \code{TRUE} on every row, selects RUMM2030's average item anchoring
#'   instead: every parameter is estimated free, and the calibration is
#'   shifted so that the mean location of the anchor items equals the mean
#'   of their \code{tau} values (one row per item, \code{k = NA}). Only the
#'   origin changes, so every item, the anchors included, keeps its
#'   estimated position relative to the others. Numeric-threshold and
#'   item-location anchor values are treated as fixed constants: their
#'   uncertainty is not included in the returned covariance or standard
#'   errors. With several numeric-threshold anchors, their stated relative
#'   spacing is a model constraint as well as a choice of origin. Column names
#'   must be unique. PCM only.
#' @param maxit,tol Newton-Raphson iteration cap and convergence tolerance.
#' @return A list containing the threshold table \code{thr}, covariance matrix
#'   \code{cov_tau}, pairwise conditional log-likelihood, iteration count,
#'   convergence flag, notes, and maximum scores \code{m}. In \code{thr},
#'   \code{weak} marks all thresholds of an item with fewer than eight
#'   responses in any category, or a threshold adjacent to a category with
#'   fewer than three responses. Standard errors for weak thresholds are
#'   reported as \code{NA}. If estimation does not converge, the function
#'   warns and all standard errors and covariance entries are \code{NA}.
#'   Sandwich uncertainty is also withheld when fewer than 10 independent
#'   persons, fewer than 8 effective persons, no more effective persons than
#'   fitted parameters, or a rank-deficient person-score covariance cannot
#'   support inference. Effective support is based on informative conditional
#'   item-pair contributions; point estimates and exact anchors remain.
#' @references
#' Zwinderman, A. H. (1995). Pairwise parameter estimation in Rasch models.
#' Applied Psychological Measurement, 19(4), 369--375.
#' @examples
#' set.seed(1)
#' d <- seq(-1.5, 1.5, length.out = 6)
#' X <- matrix(rbinom(400 * 6, 1, plogis(outer(rnorm(400), d, "-"))), 400, 6)
#' colnames(X) <- paste0("I", 1:6)
#' pcml(X)$thr
#' # anchor two items at fixed values (equating)
#' anchors <- data.frame(item = c("I1", "I6"), k = 1,
#'                       tau = c(-1.5, 1.5))
#' pcml(X, anchors = anchors)$thr
#' @export
pcml <- function(X, model = c("PCM", "RSM"), anchors = NULL,
                 maxit = 60, tol = 1e-8) {
  out <- .pcml_fit(X, model = model, anchors = anchors,
                   maxit = maxit, tol = tol)
  if (!isTRUE(out$converged))
    warning("estimation did NOT converge in ", out$iterations,
            " iterations; standard errors and covariance are unavailable",
            call. = FALSE)
  out
}

# Internal entry used by rasch() when several response rows belong to the
# same person. The public low-level pcml() interface continues to regard its
# matrix rows as independent units.
.pcml_fit <- function(X, model = c("PCM", "RSM"), anchors = NULL,
                      maxit = 60, tol = 1e-8, cluster = NULL) {
  model <- match.arg(model)
  .check_controls(maxit, tol)
  if (!is.null(colnames(X)) && anyDuplicated(colnames(X)))
    stop("item column names must be unique: ",
         paste(unique(colnames(X)[duplicated(colnames(X))]), collapse = ", "))
  X <- .pcml_score_matrix(X)
  m <- apply(X, 2, max, na.rm = TRUE); L <- ncol(X)
  thr <- threshold_index(m); M <- nrow(thr)
  inames <- if (is.null(colnames(X))) paste0("V", seq_len(L)) else colnames(X)
  pairs <- .pair_counts(X, m)
  weak <- .pcml_weak_thresholds(X, m, thr, inames)

  average <- FALSE
  if (!is.null(anchors)) {
    if (!is.data.frame(anchors))
      stop("anchors must be a data frame with columns item, k, tau")
    .check_column_names(anchors)
    if (model != "PCM") stop("anchoring is supported for the PCM only")
    if (!all(c("item", "k", "tau") %in% names(anchors)))
      stop("anchors needs columns item, k, tau")
    if (!nrow(anchors))
      stop("anchors must contain at least one fixed threshold or item mean")
    anchors <- .pcml_anchor_columns(anchors)
    a_item <- if (is.character(anchors$item) || is.factor(anchors$item))
      match(as.character(anchors$item), colnames(X))
    else {
      ai <- anchors$item
      if (!is.numeric(ai) || any(!is.finite(ai)) || any(ai != floor(ai)) ||
          any(ai < 1) || any(ai > ncol(X)))
        stop("numeric anchor item indices must be whole numbers between 1 and ",
             ncol(X))
      as.integer(ai)
    }
    if (anyNA(a_item)) stop("anchor item(s) not found among the item columns")
    # average = TRUE selects RUMM's average item anchoring: the calibration
    # is estimated free and then shifted so the mean location of the anchor
    # items equals the mean of their anchor values. No item is fixed, so it
    # is handled by the free branches below with a different origin
    if ("average" %in% names(anchors)) {
      avg <- anchors$average
      if (any(avg) && !all(avg))
        stop("average anchoring applies to the whole anchor set: `average` ",
             "must be TRUE on every row, or absent")
      average <- all(avg)
    }
    if (average) {
      if (!all(is.na(anchors$k)))
        stop("average anchoring shifts the calibration by item locations: ",
             "give one row per anchor item with k = NA")
      if (anyDuplicated(a_item))
        stop("duplicate average anchor(s) for an item")
    }
  }
  if (!is.null(anchors) && !average) {
    # k = NA anchors the item's mean location (location anchoring) with its
    # thresholds free; a numeric k fixes that single threshold. For a
    # dichotomous item the two coincide.
    is_mean <- is.na(anchors$k)
    conv <- is_mean & m[a_item] == 1L
    anchors$k[conv] <- 1L; is_mean[conv] <- FALSE

    ft_item <- a_item[!is_mean]
    a_id <- thr$id[match(paste(ft_item, anchors$k[!is_mean]),
                         paste(thr$item, thr$k))]
    if (anyNA(a_id)) stop("anchor threshold number(s) out of range for the item")
    if (anyDuplicated(a_id)) stop("duplicate anchor threshold(s)")
    mean_items <- a_item[is_mean]; mean_tau <- anchors$tau[is_mean]
    if (anyDuplicated(mean_items)) stop("duplicate location anchor(s) for an item")
    if (length(intersect(mean_items, ft_item)))
      stop("an item cannot carry both a location anchor and threshold anchors")

    offset <- numeric(M)
    offset[a_id] <- anchors$tau[!is_mean]
    for (j in seq_along(mean_items))
      offset[thr$item == mean_items[j]] <- mean_tau[j]

    # start values shifted onto the anchored scale
    st <- .start_tau(X, thr)
    shifts <- c(anchors$tau[!is_mean] - st[a_id],
                vapply(seq_along(mean_items), function(j)
                  mean_tau[j] - mean(st[thr$item == mean_items[j]]), 0))
    st <- st + mean(shifts)

    # design: identity columns for plain free thresholds; a sum-zero spread
    # block for each location-anchored item; nothing for fixed thresholds
    plain <- which(!(seq_len(M) %in% a_id) & !(thr$item %in% mean_items))
    blocks <- list()
    if (length(plain)) blocks$plain <- list(B = diag(M)[, plain, drop = FALSE],
                                            beta0 = st[plain])
    for (j in seq_along(mean_items)) {
      rows <- which(thr$item == mean_items[j]); mi <- length(rows)
      A <- matrix(0, M, mi - 1L)
      A[rows, ] <- rbind(diag(mi - 1L), rep(-1, mi - 1L))
      s <- st[rows] - mean(st[rows])
      blocks[[paste0("mean", j)]] <- list(B = A, beta0 = s[-mi])
    }
    if (!length(blocks)) stop("at least one parameter must remain free")
    B <- do.call(cbind, lapply(blocks, `[[`, "B"))
    beta0 <- unlist(lapply(blocks, `[[`, "beta0"), use.names = FALSE)

    .pcml_check_uninformative(X, m, inames, fixed = a_id, thr = thr)
    .pcml_check_connected(pairs, L, inames,
                          anchored = unique(c(ft_item, mean_items)))
    sol <- .pcml_solve(X, thr, m, B, beta0, offset = offset,
                       maxit = maxit, tol = tol, pairs = pairs,
                       cluster = cluster)
    thr$tau <- sol$tau; thr$se <- sol$se_tau; thr$se[a_id] <- 0
    thr$anchored <- seq_len(M) %in% a_id | thr$item %in% mean_items
    # location anchoring (k = NA) fixes only an item's MEAN location; its
    # individual thresholds stay free and estimated. Only genuinely fixed
    # thresholds (a_id) may suppress the weak-category flag -- a
    # location-anchored item's free threshold sitting on a near-empty category
    # is still a boundary artefact and must keep weak = TRUE / se = NA
    thr$weak <- weak$flag & !(seq_len(M) %in% a_id)
    thr$se[thr$weak] <- NA_real_
    if (!isTRUE(sol$converged)) {
      thr$se[] <- NA_real_
      thr$se[a_id] <- 0
      sol$cov_tau[,] <- NA_real_
      if (length(a_id)) {
        sol$cov_tau[a_id, ] <- 0
        sol$cov_tau[, a_id] <- 0
      }
      sol$cov_beta[,] <- NA_real_
    }
    return(list(model = model, thr = thr, cov_tau = sol$cov_tau,
                loglik = sol$loglik, iterations = sol$iterations,
                converged = sol$converged, m = m, anchors = anchors,
                n_parameters = ncol(B), B = B, cov_beta = sol$cov_beta,
                H_beta = sol$H_beta,
                notes = c(weak$notes, sol$cluster_note),
                cluster_inference = sol$cluster_inference,
                cluster_support = sol$cluster_support))
  }

  if (model == "RSM") {
    if (length(unique(m)) != 1L) stop("RSM requires equal max score across items")
    mm <- m[1]
    # parameters: delta_1..delta_{L-1}, kappa_1..kappa_{mm-1}; sum-zero each
    P <- (L - 1L) + (mm - 1L)
    B <- matrix(0, M, P)
    for (row in seq_len(M)) {
      i <- thr$item[row]; k <- thr$k[row]
      if (i < L) B[row, i] <- 1 else B[row, seq_len(L - 1L)] <- -1
      if (mm > 1L) {
        if (k < mm) B[row, L - 1L + k] <- B[row, L - 1L + k] + 1
        else B[row, L:(L + mm - 2L)] <- B[row, L:(L + mm - 2L)] - 1
      }
    }
    st <- .start_tau(X, thr)
    del <- vapply(seq_len(L), function(i) mean(st[thr$item == i]), 0)
    kap <- vapply(seq_len(mm), function(k) mean(st[thr$k == k] -
                    del[thr$item[thr$k == k]]), 0)
    del <- del - mean(del); kap <- kap - mean(kap)
    beta0 <- c(del[-L], if (mm > 1L) kap[-mm] else numeric(0))
  } else {
    # sum-zero over all thresholds during estimation; recentred afterwards
    B <- rbind(diag(M - 1L), rep(-1, M - 1L))
    st <- .start_tau(X, thr)
    beta0 <- st[-M] - mean(st)
  }

  if (model == "PCM") .pcml_check_uninformative(X, m, inames, thr = thr)
  .pcml_check_connected(pairs, L, inames)
  sol <- .pcml_solve(X, thr, m, B, beta0, maxit = maxit, tol = tol,
                     pairs = pairs, cluster = cluster)
  # fix the origin -- mean item location zero, or under average anchoring
  # the mean location of the anchor items at the mean of their anchor
  # values -- and move the covariance to that parameterisation with it. The
  # origin is the linear functional c = a'tau with a = 1/(L m_i) on item i's
  # thresholds (1/(A m_i) over the A anchor items); tau_new = (I - 1 a') tau
  # + target, so cov_new = (I - 1 a') cov (I - a 1'). Under equal max scores
  # the recentring constant is identically zero on the estimation constraint
  # (sum of all thresholds zero) and the transform is a no-op, which is why
  # equal-m calibration checks could not see its absence; with MIXED max
  # scores the untransformed covariance mis-states every threshold and
  # item-location SE (verified by simulation: per-item empSD/SE 0.90-1.21
  # before, 0.96-1.09 after).
  if (average) {
    a_c <- numeric(M); in_anchor <- thr$item %in% a_item
    a_c[in_anchor] <- 1 / (length(a_item) * m[thr$item[in_anchor]])
    target <- mean(anchors$tau)
  } else {
    a_c <- 1 / (L * m[thr$item]); target <- 0
  }
  sol$tau <- sol$tau - sum(a_c * sol$tau) + target
  A_c <- diag(M) - matrix(1, M, 1) %*% t(a_c)
  sol$cov_tau <- A_c %*% sol$cov_tau %*% t(A_c)
  sol$se_tau <- sqrt(pmax(diag(sol$cov_tau), 0))
  if (!isTRUE(sol$converged)) {
    sol$se_tau[] <- NA_real_
    sol$cov_tau[,] <- NA_real_
    sol$cov_beta[,] <- NA_real_
  }
  thr$tau <- sol$tau; thr$se <- sol$se_tau; thr$anchored <- FALSE
  thr$weak <- weak$flag
  thr$se[thr$weak] <- NA_real_

  list(model = model, thr = thr, cov_tau = sol$cov_tau,
       loglik = sol$loglik, iterations = sol$iterations,
       converged = sol$converged, m = m,
       anchors = if (average) anchors else NULL,
       n_parameters = ncol(B), B = B, cov_beta = sol$cov_beta,
       H_beta = sol$H_beta, notes = c(weak$notes, sol$cluster_note),
       cluster_inference = sol$cluster_inference,
       cluster_support = sol$cluster_support)
}

# ---------------------------------------------------------------------------
# Andrich principal-components reparameterisation (Andrich 1978, 1985; Pedler
# 1987): an optional alternative to the free-threshold pcml() above, useful
# when some categories are sparsely populated. Each item's mi thresholds are
# re-expressed as up to four orthogonal-polynomial components in the
# category score x = 0, ..., mi: location (linear), spread (quadratic),
# skewness (cubic), and kurtosis (quartic),
#
#   L_i(x) = x.omega_1i - x(mi-x).omega_2i - x(mi-x)(2x-mi).omega_3i - ...
#
# with tau_ik = L_i(k) - L_i(k-1) (Andrich 1985, eqs 1.6-1.7). Every
# component's coefficient differences telescope to zero across an item's
# thresholds except location's (which differences to a constant 1), so
# location is exactly the item's mean threshold and carries the same
# across-item additive-shift redundancy as in .start_tau; it is sum-zero
# constrained across items the same way pcml()'s RSM delta is. The
# higher-order components are free per item, capped at an item's own
# threshold count (a dichotomous item has location only). The family stops
# at the quartic (kurtosis) term. With four components it spans the free PCM
# when every item has at most four thresholds. Longer scales are restricted
# to the four-component polynomial trend.
# ---------------------------------------------------------------------------
.pc_gcoefs <- function(mi) {
  x <- 0:mi
  G <- cbind(x,
             -x * (mi - x),
             -x * (mi - x) * (2 * x - mi),
             -x * (mi - x) * (5 * x^2 - 5 * x * mi + mi^2 + 1))
  G[-1, , drop = FALSE] - G[-nrow(G), , drop = FALSE]
}

# Retain only components supported by the item's threshold count. Incremental
# rank checking also protects against numerical degeneracy in the design.
.pc_select <- function(G, ncomp, tol = 1e-8) {
  if (ncomp < 2L) return(integer(0))
  base <- G[, 1, drop = FALSE]; r0 <- qr(base, tol = tol)$rank
  keep <- integer(0)
  for (l in 2:ncomp) {
    test <- cbind(base, G[, l])
    r1 <- qr(test, tol = tol)$rank
    if (r1 > r0) { keep <- c(keep, l); base <- test; r0 <- r1 }
  }
  keep
}

#' Estimate Rasch thresholds using a principal-component parameterisation
#'
#' Re-expresses each item's thresholds as orthogonal polynomial components:
#' location, spread, skewness, and kurtosis (Andrich and Luo 2003).
#' Estimation uses the same pairwise conditional likelihood as
#' \code{\link{pcml}}. With at most four thresholds per item the full
#' parameterisation is exact. Items with five or more thresholds are fitted by
#' a reduced-rank polynomial trend, which can stabilise sparse categories at
#' the cost of restricting the threshold pattern.
#'
#' @details For scores \eqn{x=0,\ldots,m}, the cumulative threshold function is
#' \deqn{C(x)=x\omega_1-x(m-x)\omega_2-x(m-x)(2x-m)\omega_3
#'       -x(m-x)(5x^2-5mx+m^2+1)\omega_4.}
#' Threshold \eqn{k} is \eqn{C(k)-C(k-1)}. The four coefficients are
#' location, spread, skewness and kurtosis, respectively.
#'
#' @param X Persons-by-items integer score matrix. Each item must have at least
#'   two observed categories, numbered consecutively from 0. Missing values are
#'   handled by pairwise deletion.
#' @param n_components Maximum number of components per item: 1 (location
#'   only) up to 4 (location, spread, skewness, kurtosis).
#'   Capped per item at its own number of thresholds. Kurtosis requires
#'   at least four thresholds (five response categories).
#' @param maxit,tol Newton-Raphson iteration cap and convergence tolerance.
#' @return A list with the threshold table \code{thr} (columns \code{id},
#'   \code{item}, \code{k}, \code{tau}, \code{se}), the component table
#'   \code{components} (one row per item, with \code{location},
#'   \code{spread}, \code{skewness}, \code{kurtosis} and their standard
#'   errors, \code{NA} where an item's rank does not support that
#'   component), the threshold covariance matrix \code{cov_tau}, the
#'   pairwise conditional log-likelihood, the iteration count, a convergence
#'   flag, and the max-score vector \code{m}. If estimation does not converge,
#'   the function warns and all standard errors and covariance entries are
#'   \code{NA}. The independent-person and effective-support conditions
#'   described for \code{\link{pcml}} also apply.
#' @references
#' Andrich, D. and Luo, G. (2003). Conditional pairwise estimation in the
#' Rasch model for ordered response categories using principal components.
#' Journal of Applied Measurement, 4(3), 205--221.
#'
#' Zwinderman, A. H. (1995). Pairwise parameter estimation in Rasch models.
#' Applied Psychological Measurement, 19(4), 369--375.
#'
#' Andrich, D. (1978). A rating formulation for ordered response categories.
#' Psychometrika, 43(4), 561--573.
#'
#' Andrich, D. (1985). An elaboration of Guttman scaling with Rasch models
#' for measurement. In N. B. Tuma (Ed.), Sociological Methodology 1985
#' (pp. 33--80). Jossey-Bass.
#'
#' Pedler, P. J. (1987). Accounting for psychometric dependence with a class
#' of latent trait models. PhD thesis, University of Western Australia.
#' @examples
#' set.seed(1)
#' d <- seq(-1.5, 1.5, length.out = 6)
#' X <- matrix(rbinom(400 * 6, 1, plogis(outer(rnorm(400), d, "-"))), 400, 6)
#' colnames(X) <- paste0("I", 1:6)
#' pcml_pc(X)$components
#' @export
pcml_pc <- function(X, n_components = 4, maxit = 60, tol = 1e-8) {
  out <- .pcml_pc_fit(X, n_components = n_components,
                      maxit = maxit, tol = tol)
  if (!isTRUE(out$converged))
    warning("estimation did NOT converge in ", out$iterations,
            " iterations; standard errors and covariance are unavailable",
            call. = FALSE)
  out
}

.pcml_pc_fit <- function(X, n_components = 4, maxit = 60, tol = 1e-8,
                         cluster = NULL) {
  .check_controls(maxit, tol)
  n_components <- .check_whole(n_components, "n_components", 1, 4)
  if (!is.null(colnames(X)) && anyDuplicated(colnames(X)))
    stop("item column names must be unique: ",
         paste(unique(colnames(X)[duplicated(colnames(X))]), collapse = ", "))
  X <- .pcml_score_matrix(X)
  m <- apply(X, 2, max, na.rm = TRUE); L <- ncol(X)
  thr <- threshold_index(m); M <- nrow(thr)
  ncomp <- pmin(m, n_components, 4L)

  Gs <- lapply(m, .pc_gcoefs)
  keep <- Map(.pc_select, Gs, ncomp)
  extra <- lengths(keep)

  # column layout: (L - 1) sum-zero location columns, then one free column
  # per item for each higher-order component actually identified for it
  P <- (L - 1L) + sum(extra)
  B <- matrix(0, M, P)
  off <- (L - 1L) + cumsum(c(0L, extra))[seq_len(L)]

  for (i in seq_len(L)) {
    rows <- which(thr$item == i)
    G <- Gs[[i]]
    if (i < L) B[rows, i] <- G[, 1] else B[rows, seq_len(L - 1L)] <- -G[, 1]
    if (extra[i] > 0L)
      B[rows, off[i] + seq_len(extra[i])] <- G[, keep[[i]], drop = FALSE]
  }

  st <- .start_tau(X, thr)
  loc <- vapply(seq_len(L), function(i) mean(st[thr$item == i]), 0)
  loc <- loc - mean(loc)
  beta0 <- c(loc[-L], numeric(sum(extra)))

  inames <- if (is.null(colnames(X))) paste0("V", seq_len(L)) else colnames(X)
  pairs <- .pair_counts(X, m)
  .pcml_check_connected(pairs, L, inames)
  sol <- .pcml_solve(X, thr, m, B, beta0, maxit = maxit, tol = tol,
                     pairs = pairs, cluster = cluster)
  thr$tau <- sol$tau; thr$se <- sol$se_tau

  labs <- c("spread", "skewness", "kurtosis")
  comp <- data.frame(item = inames,
                     location = c(sol$beta[seq_len(L - 1L)],
                                  -sum(sol$beta[seq_len(L - 1L)])),
                     location_se = NA_real_,
                     spread = NA_real_, spread_se = NA_real_,
                     skewness = NA_real_, skewness_se = NA_real_,
                     kurtosis = NA_real_, kurtosis_se = NA_real_)
  Bloc <- rbind(diag(L - 1L), rep(-1, L - 1L))
  cov_loc <- Bloc %*% sol$cov_beta[seq_len(L - 1L), seq_len(L - 1L), drop = FALSE] %*% t(Bloc)
  comp$location_se <- sqrt(pmax(diag(cov_loc), 0))
  for (i in seq_len(L)) if (extra[i] > 0L) {
    cols <- off[i] + seq_len(extra[i])
    for (j in seq_along(keep[[i]])) {
      lab <- labs[keep[[i]][j] - 1L]
      comp[[lab]][i] <- sol$beta[cols[j]]
      comp[[paste0(lab, "_se")]][i] <- sqrt(pmax(sol$cov_beta[cols[j], cols[j]], 0))
    }
  }
  if (!isTRUE(sol$converged)) {
    thr$se[] <- NA_real_
    se_cols <- grep("_se$", names(comp), value = TRUE)
    for (nm in se_cols) comp[[nm]][] <- NA_real_
    sol$cov_tau[,] <- NA_real_
    sol$cov_beta[,] <- NA_real_
  }

  list(model = "PCM", n_components = n_components,
       pc_algorithm = "guttman-four-1", thr = thr,
       components = comp, cov_tau = sol$cov_tau, loglik = sol$loglik,
       iterations = sol$iterations, converged = sol$converged, m = m,
       anchors = NULL, n_parameters = ncol(B), B = B,
       cov_beta = sol$cov_beta, H_beta = sol$H_beta,
       notes = sol$cluster_note,
       cluster_inference = sol$cluster_inference,
       cluster_support = sol$cluster_support)
}
