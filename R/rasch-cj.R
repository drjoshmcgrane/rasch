# rasch :: item responses and comparative judgements in one calibration
# ===========================================================================
# Items can be calibrated from two kinds of evidence: persons answering them,
# and judges comparing them. Each kind is a frame of reference in the sense of
# Humphry and Andrich (2008), and each has its own unit. This estimator writes
# one set of item parameters into both, scales the judgement frames by a unit
# relative to the response frame, and maximises the sum of three likelihoods
# that are each free of person and judge parameters:
#
#   responses    conditional likelihood of each person's pattern given the
#                raw score (CML with elementary symmetric functions), unit 1;
#   comparisons  Bradley-Terry, P(a beats b) = logistic(alpha (lambda_a -
#                lambda_b)), the judge's own location cancelling in the
#                difference;
#   rankings     Plackett-Luce, a ranking of n objects as n - 1 sequential
#                choices with P(next is a) = exp(kappa lambda_a) / sum over
#                the objects still unplaced, which reduces to Bradley-Terry
#                when n = 2.
#
# The item parameters are the thresholds psi_ik of the partial credit model,
# k = 1..m_i, a dichotomous item having the single threshold psi_i1. What a
# judge compares is an object: an item, whose location lambda is the mean of
# its thresholds, or a single threshold of an item, whose location is that
# threshold. The judgement blocks therefore act on lambda = O psi for an
# object matrix O, and the response block on psi directly.
#
# Because every block is a likelihood the sum is one too: there is no weight
# to choose, the inverse observed information is the covariance, and the
# likelihood ratio test of the model against separate locations per frame is
# a test of the invariance the model assumes. The pairwise alternative,
# turning each response pattern into item comparisons and pooling them with
# the judgements, is a composite likelihood that over-counts the derived
# pairs; it gives the same locations in expectation but needs a Godambe
# weight and a sandwich to be efficient and honest. This function is the
# full-likelihood route and leaves the pairwise estimators untouched.
#
# The response frame fixes the unit; alpha and kappa are the units of the
# judgement frames relative to it. An alpha of 0.5 means judges separate the
# objects at half the discrimination of the test, so a logit on the judges'
# scale is half a logit on the test's. Units are estimated on the log scale
# so positivity is automatic, and reported on the natural scale.
#
# Identification: the mean item location is zero. Every item must be reached
# by at least one frame, and the item graph over all frames must be
# connected, or the relative location of a block is not identified.
# ===========================================================================

# Category terms of the partial credit model for thresholds psi: eta[[i]][x
# + 1] = exp(-(psi_i1 + ... + psi_ix)), with eta[[i]][1] = 1.
.cj_eta <- function(psi, o)
  lapply(o$idx, function(ix) c(1, exp(-cumsum(psi[ix]))))

# Person-pattern elementary symmetric functions. Row n of A marks the items
# person n answered; G[n, r + 1] is gamma_r over those items, the sum over
# their response patterns with raw score r of the product of category terms,
# so persons with different answered sets get their own gamma. `skip` omits
# items, which is how the derivative terms gamma^(i) and gamma^(ij) are
# formed. Each item multiplies the running polynomial by its own.
.cj_esf_rows <- function(eta, A, skip = integer(0)) {
  N <- nrow(A); R <- sum(lengths(eta) - 1L)
  G <- matrix(0, N, R + 1L); G[, 1L] <- 1
  for (i in seq_along(eta)) {
    if (i %in% skip) next
    rows <- A[, i]
    if (!any(rows)) next
    e <- eta[[i]]; mi <- length(e) - 1L
    Gi <- G[rows, , drop = FALSE]; Gn <- Gi
    for (x in seq_len(mi))
      Gn[, (x + 1L):(R + 1L)] <- Gn[, (x + 1L):(R + 1L), drop = FALSE] +
        e[x + 1L] * Gi[, seq_len(R + 1L - x), drop = FALSE]
    G[rows, ] <- Gn
  }
  G
}

# Collapse the response matrix to the distinct (answered set, raw score)
# patterns that the conditional likelihood depends on, with a count each.
# Persons with every answered item at its minimum or maximum contribute a
# constant and are dropped; persons answering fewer than two items likewise.
# `m` gives the number of thresholds of each item; T counts, for each
# threshold, the kept persons who scored at or above it.
.cj_cml_prep <- function(X, m = NULL) {
  X <- as.matrix(X)
  if (is.null(m)) m <- pmax(apply(X, 2L, function(v) max(c(0L, v), na.rm = TRUE)), 1L)
  m <- as.integer(m)
  A <- !is.na(X)
  Xz <- X; Xz[!A] <- 0L
  r <- as.integer(rowSums(Xz))
  n_obs <- rowSums(A)
  max_r <- as.vector(A %*% m)
  keep <- n_obs >= 2L & r > 0L & r < max_r
  idx <- split(seq_len(sum(m)), rep(seq_along(m), m))
  T_ <- unlist(lapply(seq_along(m), function(i)
    vapply(seq_len(m[i]), function(k) sum(Xz[keep, i] >= k), 0)), use.names = FALSE)
  key <- paste(apply(A[keep, , drop = FALSE], 1L, paste, collapse = ""), r[keep])
  first <- !duplicated(key)
  list(A = A[keep, , drop = FALSE][first, , drop = FALSE],
       r = r[keep][first],
       w = as.numeric(table(factor(key, levels = key[first]))),
       T = T_, m = m, idx = idx, P = sum(m), I = ncol(X), n_persons = sum(keep))
}

.cj_cml_ll <- function(psi, o) {
  if (!nrow(o$A)) return(0)
  eta <- .cj_eta(psi, o)
  G <- .cj_esf_rows(eta, o$A)
  g <- G[cbind(seq_len(nrow(o$A)), o$r + 1L)]
  if (any(!is.finite(g)) || any(g <= 0)) return(-Inf)
  -sum(psi * o$T) - sum(o$w * log(g))
}

# Conditional category probabilities for each pattern row. Q[[i]][n, x + 1]
# is P(x_ni = x | r_n), zero where unanswered; C[[i]][n, k] is P(x_ni >= k |
# r_n), the conditional expectation of the threshold indicator.
.cj_cml_p <- function(psi, o) {
  eta <- .cj_eta(psi, o); N <- nrow(o$A)
  G <- .cj_esf_rows(eta, o$A)
  gr <- G[cbind(seq_len(N), o$r + 1L)]
  Q <- vector("list", o$I); C <- vector("list", o$I)
  for (i in seq_len(o$I)) {
    mi <- o$m[i]
    Gi <- .cj_esf_rows(eta, o$A, skip = i)
    Qi <- matrix(0, N, mi + 1L)
    for (x in 0:mi) {
      ok <- o$A[, i] & o$r >= x
      if (any(ok))
        Qi[ok, x + 1L] <- eta[[i]][x + 1L] * Gi[cbind(which(ok), o$r[ok] - x + 1L)] / gr[ok]
    }
    Ci <- matrix(0, N, mi)
    for (k in seq_len(mi)) Ci[, k] <- rowSums(Qi[, (k + 1L):(mi + 1L), drop = FALSE])
    Q[[i]] <- Qi; C[[i]] <- Ci
  }
  list(Q = Q, C = C, gr = gr, eta = eta)
}

.cj_cml_grad <- function(psi, o, pp = NULL) {
  if (!nrow(o$A)) return(numeric(o$P))
  if (is.null(pp)) pp <- .cj_cml_p(psi, o)
  -o$T + unlist(lapply(pp$C, function(Ci) colSums(o$w * Ci)), use.names = FALSE)
}

# Minus the conditional covariance of the threshold indicators, summed over
# patterns. Within an item the indicators are nested, so the joint
# probability of two is the tail of the higher; across items it needs the
# pattern functions with both items removed.
.cj_cml_hess <- function(psi, o, pp = NULL) {
  if (!nrow(o$A)) return(matrix(0, o$P, o$P))
  if (is.null(pp)) pp <- .cj_cml_p(psi, o)
  eta <- pp$eta; gr <- pp$gr; C <- pp$C; N <- nrow(o$A); I <- o$I
  H <- matrix(0, o$P, o$P)
  for (i in seq_len(I)) {
    ix <- o$idx[[i]]; mi <- o$m[i]; Ci <- C[[i]]
    for (k in seq_len(mi)) for (l in seq_len(mi))
      H[ix[k], ix[l]] <- -sum(o$w * (Ci[, max(k, l)] - Ci[, k] * Ci[, l]))
  }
  for (i in seq_len(I - 1L)) for (j in (i + 1L):I) {
    mi <- o$m[i]; mj <- o$m[j]
    Gij <- .cj_esf_rows(eta, o$A, skip = c(i, j))
    both <- o$A[, i] & o$A[, j]
    # joint conditional distribution of the two responses, then its upper
    # tails, which are the joint probabilities of the threshold indicators
    J <- array(0, c(N, mi + 1L, mj + 1L))
    for (x in 0:mi) for (y in 0:mj) {
      ok <- both & o$r >= x + y
      if (any(ok))
        J[ok, x + 1L, y + 1L] <- eta[[i]][x + 1L] * eta[[j]][y + 1L] *
          Gij[cbind(which(ok), o$r[ok] - x - y + 1L)] / gr[ok]
    }
    for (x in rev(seq_len(mi))) J[, x, ] <- J[, x, ] + J[, x + 1L, ]
    for (y in rev(seq_len(mj))) J[, , y] <- J[, , y] + J[, , y + 1L]
    ix <- o$idx[[i]]; jx <- o$idx[[j]]
    for (k in seq_len(mi)) for (l in seq_len(mj)) {
      s <- -sum(o$w * (J[, k + 1L, l + 1L] - C[[i]][, k] * C[[j]][, l]))
      H[ix[k], jx[l]] <- s; H[jx[l], ix[k]] <- s
    }
  }
  H
}

# Bradley-Terry block from the ordered count matrix W, W[i, j] the number of
# comparisons i won against j, with unit u on the difference.
.cj_bt_ll <- function(delta, u, W) {
  D <- u * outer(delta, delta, "-")
  sum(W * stats::plogis(D, log.p = TRUE))
}

# Gradient and Hessian in (delta, log u).
.cj_bt_parts <- function(delta, u, W) {
  Z <- outer(delta, delta, "-")
  P <- stats::plogis(u * Z)
  M <- W + t(W)
  V <- M * P * (1 - P); diag(V) <- 0
  R <- W - M * P; diag(R) <- 0
  g_delta <- u * rowSums(R)
  H_dd <- u^2 * V; diag(H_dd) <- -rowSums(H_dd)
  g_u <- sum(W * (1 - P) * Z)
  H_uu <- -sum(W * P * (1 - P) * Z^2)
  A <- W * ((1 - P) - u * P * (1 - P) * Z); diag(A) <- 0
  H_du <- rowSums(A) - colSums(A)
  # chain rule to the log unit: d/dlog u = u d/du
  list(g = c(g_delta, u * g_u),
       H = rbind(cbind(H_dd, u * H_du), c(u * H_du, u^2 * H_uu + u * g_u)))
}

# Plackett-Luce block. `rk` is a list of integer vectors, each a ranking of
# object indices from first to last. Each ranking contributes one choice per
# position but the last: the object placed, against everything still
# unplaced.
.cj_pl_ll <- function(delta, u, rk) {
  ll <- 0
  for (v in rk) {
    s <- u * delta[v]; n <- length(s)
    for (k in seq_len(n - 1L)) {
      rest <- s[k:n]; m <- max(rest)
      ll <- ll + s[k] - m - log(sum(exp(rest - m)))
    }
  }
  ll
}

.cj_pl_parts <- function(delta, u, rk) {
  I <- length(delta)
  g_delta <- numeric(I); H_dd <- matrix(0, I, I)
  g_u <- 0; H_uu <- 0; H_du <- numeric(I)
  for (v in rk) {
    n <- length(v)
    for (k in seq_len(n - 1L)) {
      idx <- v[k:n]; d <- delta[idx]
      s <- u * d; p <- exp(s - max(s)); p <- p / sum(p)
      e <- numeric(length(idx)); e[1L] <- 1
      Vp <- diag(p, nrow = length(p)) - tcrossprod(p)      # cov of the choice
      g_delta[idx] <- g_delta[idx] + u * (e - p)
      H_dd[idx, idx] <- H_dd[idx, idx] - u^2 * Vp
      g_u <- g_u + d[1L] - sum(p * d)
      H_uu <- H_uu - drop(crossprod(d, Vp %*% d))
      H_du[idx] <- H_du[idx] + (e - p) - u * drop(Vp %*% d)
    }
  }
  list(g = c(g_delta, u * g_u),
       H = rbind(cbind(H_dd, u * H_du), c(u * H_du, u^2 * H_uu + u * g_u)))
}

# Newton ascent with step halving. The joint likelihood in the free units is
# not concave everywhere, so a Hessian that is not negative definite has its
# eigenvalues reflected before the step is taken. Convergence is judged on
# the parameter scale, by the size of the remaining Newton step.
.cj_newton <- function(par0, fn, gr, he, maxit = 200L, tol = 1e-8) {
  ascent_step <- function(H, g) {
    st <- tryCatch(solve(H, -g), error = function(e) NULL)
    if (!is.null(st) && all(is.finite(st)) && sum(st * g) > 0) return(st)
    ev <- eigen((H + t(H)) / 2, symmetric = TRUE)
    lam <- -pmax(abs(ev$values), 1e-8)
    as.vector(ev$vectors %*% (crossprod(ev$vectors, -g) / lam))
  }
  par <- par0; f <- fn(par)
  if (!is.finite(f)) stop("the log-likelihood is not finite at the starting values",
                          call. = FALSE)
  if (!length(par)) return(list(par = par, ll = f, H = matrix(0, 0L, 0L),
                                converged = TRUE, iterations = 0L))
  converged <- FALSE; iterations <- 0L
  for (it in seq_len(maxit)) {
    iterations <- it
    g <- gr(par); H <- he(par)
    step <- ascent_step(H, g)
    if (any(!is.finite(step))) step <- g / max(1, max(abs(g)))
    if (max(abs(step)) < tol) { converged <- TRUE; break }
    s <- 1; moved <- FALSE
    for (k in seq_len(60L)) {
      cand <- par + s * step; f_new <- fn(cand)
      if (is.finite(f_new) && f_new > f) { moved <- TRUE; break }
      s <- s / 2
    }
    if (!moved) { converged <- max(abs(step)) < sqrt(tol); break }
    par <- cand; f <- f_new
  }
  H <- he(par)
  list(par = par, ll = f, H = H, converged = converged, iterations = iterations)
}

# Contrast basis for a weighted sum-zero constraint: par = B beta with beta
# free and sum(c * par) = 0. Equal weights give the usual sum-zero basis.
.cj_basis <- function(n, c = rep(1, n)) {
  if (n < 2L) return(matrix(0, n, 0L))
  rbind(diag(n - 1L), -c[-n] / c[n])
}

# Basis with one sum-zero constraint per group, for a frame whose objects
# fall into disconnected blocks: block-diagonal in the group bases.
.cj_basis_groups <- function(group) {
  n <- length(group); B <- matrix(0, n, 0L)
  for (g in unique(group)) {
    rows <- which(group == g); Bg <- .cj_basis(length(rows))
    if (!ncol(Bg)) next
    add <- matrix(0, n, ncol(Bg)); add[rows, ] <- Bg
    B <- cbind(B, add)
  }
  B
}

# Connected components of an adjacency matrix, as an integer label per node.
.cj_components <- function(adj) {
  n <- nrow(adj); comp <- integer(n); id <- 0L
  for (s in seq_len(n)) {
    if (comp[s]) next
    id <- id + 1L
    reach <- seq_len(n) == s
    repeat {
      new <- reach | colSums(adj[reach, , drop = FALSE]) > 0
      if (all(new == reach)) break
      reach <- new
    }
    comp[reach] <- id
  }
  comp
}

# Fit one frame on its own, unit fixed at one, over the basis B. Used for
# starting values and for the separate calibrations the invariance test
# compares.
.cj_fit_alone <- function(ll, parts, B, start = NULL, maxit = 200L, tol = 1e-8) {
  nb <- ncol(B)
  b0 <- if (is.null(start) || !nb) rep(0, nb) else
    as.vector(solve(crossprod(B), crossprod(B, start)))
  fn <- function(b) ll(as.vector(B %*% b))
  gr <- function(b) as.vector(crossprod(B, parts(as.vector(B %*% b))$g))
  he <- function(b) crossprod(B, parts(as.vector(B %*% b))$H %*% B)
  fit <- .cj_newton(b0, fn, gr, he, maxit, tol)
  covb <- if (nb) tryCatch(solve(-fit$H), error = function(e) matrix(NA_real_, nb, nb))
    else matrix(0, 0L, 0L)
  cov <- B %*% covb %*% t(B)
  list(par = as.vector(B %*% fit$par), cov = cov,
       se = sqrt(pmax(diag(cov), 0)), ll = fit$ll, converged = fit$converged,
       iterations = fit$iterations, n_free = nb)
}

# Resolve item and threshold columns from a judgement source to object keys:
# the item name for an item location, "item:k" for its k-th threshold.
.cj_object_keys <- function(it, thr, items, m, what, noun = "items") {
  unknown <- setdiff(unique(it), items)
  if (length(unknown))
    stop("`", what, "` names objects that are not ", noun, ": ",
         paste(sQuote(unknown, FALSE), collapse = ", "), call. = FALSE)
  if (is.null(thr)) return(it)
  thr <- suppressWarnings(as.numeric(thr))
  has <- !is.na(thr)
  if (!any(has)) return(it)
  if (is.null(m))
    stop("threshold-level judgements need response data to define the thresholds",
         call. = FALSE)
  mi <- m[match(it, items)]
  bad <- has & (thr != floor(thr) | thr < 1 | thr > mi)
  if (any(bad))
    stop("`", what, "` names thresholds an item does not have: ",
         paste(unique(paste0(it[bad], ":", thr[bad])), collapse = ", "),
         call. = FALSE)
  ifelse(has, paste0(it, ":", thr), it)
}

# Resolve comparison rows to (winner key, loser key).
.cj_comparison_keys <- function(comparisons, object_a, object_b, winner,
                                threshold_a, threshold_b, items, m, notes,
                                noun = "items") {
  need <- c(object_a, object_b, winner)
  miss <- setdiff(need, names(comparisons))
  if (length(miss))
    stop("`comparisons` lacks column(s): ", paste(miss, collapse = ", "),
         call. = FALSE)
  a <- .role_text_values(comparisons[[object_a]])
  b <- .role_text_values(comparisons[[object_b]])
  w <- .role_text_values(comparisons[[winner]])
  if (anyNA(a) || anyNA(b) || anyNA(w))
    stop("`comparisons` has missing objects or winners", call. = FALSE)
  ka <- .cj_object_keys(a, comparisons[[threshold_a]], items, m, "comparisons", noun)
  kb <- .cj_object_keys(b, comparisons[[threshold_b]], items, m, "comparisons", noun)
  if (any(ka == kb))
    stop("`comparisons` compares an object with itself", call. = FALSE)
  tie <- tolower(w) %in% c("tie", "draw")
  if (any(tie)) {
    notes <- c(notes, sprintf("%d tied comparison(s) dropped", sum(tie)))
    a <- a[!tie]; b <- b[!tie]; w <- w[!tie]; ka <- ka[!tie]; kb <- kb[!tie]
  }
  # the winner names an object key, or the bare item when only one side is
  # that item
  win_a <- w == ka | (w == a & a != b)
  win_b <- w == kb | (w == b & a != b)
  if (any(!(win_a | win_b)))
    stop("`comparisons` has winners that are neither compared object",
         call. = FALSE)
  list(win = ifelse(win_a, ka, kb), lose = ifelse(win_a, kb, ka),
       n = length(w), notes = notes)
}

# Resolve ranking rows to a list of object-key vectors, one per ranking.
.cj_ranking_keys <- function(rankings, ranking, item, rank, threshold,
                             items, m, notes, noun = "items") {
  need <- c(ranking, item, rank)
  miss <- setdiff(need, names(rankings))
  if (length(miss))
    stop("`rankings` lacks column(s): ", paste(miss, collapse = ", "),
         call. = FALSE)
  id <- .role_text_values(rankings[[ranking]])
  it <- .role_text_values(rankings[[item]])
  rk <- suppressWarnings(as.numeric(rankings[[rank]]))
  if (anyNA(id) || anyNA(it) || anyNA(rk))
    stop("`rankings` has missing ranking identifiers, items or ranks",
         call. = FALSE)
  key <- .cj_object_keys(it, rankings[[threshold]], items, m, "rankings", noun)
  out <- vector("list", 0L); dropped <- 0L
  for (g in split(seq_along(id), id)) {
    o <- order(rk[g])
    r_sorted <- rk[g][o]; v <- key[g][o]
    if (anyDuplicated(v) || anyDuplicated(r_sorted) || length(v) < 2L) {
      dropped <- dropped + 1L
      next
    }
    out[[length(out) + 1L]] <- v
  }
  if (dropped)
    notes <- c(notes, sprintf(paste0("%d ranking(s) dropped: fewer than two ",
                                     "objects, a repeated object, or tied ranks"),
                              dropped))
  list(rk = out, n = length(out), notes = notes)
}

#' Calibrate items from responses and comparative judgements together
#'
#' Fits one set of item thresholds to persons' responses and to judges'
#' comparisons or rankings of the same items, giving each judgement frame
#' its own unit relative to the response frame, and tests whether the frames
#' agree about the items.
#'
#' @details
#' The response frame contributes the conditional likelihood of each
#' person's pattern given the raw score under the partial credit model, so
#' persons are eliminated exactly as in \code{\link{rasch}}; a dichotomous
#' item has one threshold, its location. Comparisons contribute the
#' Bradley-Terry likelihood \eqn{P(a \succ b) = \mathrm{logistic}\{\alpha
#' (\lambda_a - \lambda_b)\}}, in which a judge's own location cancels.
#' Rankings contribute the Plackett-Luce likelihood, a ranking of \eqn{n}
#' objects being \eqn{n - 1} sequential choices with \eqn{P(\text{next} = a)
#' \propto \exp(\kappa \lambda_a)} over the objects still unplaced; with
#' \eqn{n = 2} it is Bradley-Terry. The three blocks share the thresholds
#' and are summed, so the estimator is a full likelihood: standard errors
#' are from the inverse observed information and no weighting of the sources
#' is required.
#'
#' What a judge compares is an object. By default an object is an item and
#' its location \eqn{\lambda} is the mean of the item's thresholds, so
#' judgements inform the location of a polytomous item and leave the spacing
#' of its thresholds to the responses. A threshold column in the judgement
#' data makes the object a single threshold of the item, \eqn{\lambda} being
#' that threshold: a judge then says that reaching category \eqn{k} of one
#' item is harder than reaching category \eqn{l} of another. Both kinds of
#' object may appear in one source, and two thresholds of one item may be
#' compared with each other.
#'
#' \eqn{\alpha} and \eqn{\kappa} are the units of the judgement frames
#' relative to the response frame (Humphry and Andrich 2008). A unit of 0.5
#' means the judges separate the objects at half the discrimination the test
#' does. They are estimated on the log scale; \code{units} fixes either at
#' one instead, which asserts that the frame shares the test's unit.
#'
#' The model assumes each object has one location across the frames it
#' appears in. \code{invariance} tests that assumption in two ways. The
#' likelihood ratio compares the fitted model with separate locations per
#' frame (under which the units are absorbed into the locations and drop
#' out); each judgement frame then has one free location per object it
#' judges, less one per connected block of its own design, and the degrees
#' of freedom are the total of those less the number of free units. The
#' object table places each frame's separate calibration on the reference
#' scale by dividing by the fitted unit, centres it within each block, and
#' tests each object's difference from the reference calibration by a Wald
#' statistic that conditions on the fitted unit; probabilities are
#' Holm-adjusted across objects. An object that fails is one the judges and
#' the test-takers disagree about, and combining the sources moves its
#' estimate toward whichever source has the more information.
#'
#' Items are calibrated when at least one frame reaches them, and a
#' disconnected item graph is refused. An item every person answered the
#' same way carries no conditional information; its responses are left out
#' and it is located, as a single threshold, by the judgements, which must
#' reach it.
#'
#' The fitted thresholds are returned as an anchor table in the
#' \code{anchors} component. Passing it to \code{\link{rasch}} runs item and
#' person fit, targeting and the other response diagnostics on the combined
#' calibration; \code{rasch} keeps at least one item free, so pass the rows
#' of every item but one and that item is re-estimated on the anchored
#' scale.
#'
#' With \code{data = NULL} and both judgement sources supplied, the
#' comparisons become the reference frame and the rankings' unit is relative
#' to theirs; the objects are the items the two sources name, and threshold
#' columns are refused because no responses define the thresholds.
#'
#' Ties are dropped with a note; judge clustering and judge fit are not
#' provided.
#'
#' \strong{Measuring persons.} With \code{objects = "persons"} the judges
#' compare the persons' work rather than the items, and the function
#' locates each person from their responses and from those judgements
#' together. The items must be calibrated already: \code{anchors} gives
#' their thresholds, from \code{\link{rasch}}, from an item-mode
#' \code{rasch_cj} fit, or as a data frame. The response block is then the
#' partial credit likelihood of each person's responses given the anchored
#' thresholds, one location per person, and the judgement blocks are as
#' above with the persons as objects and units relative to the test. The
#' items are held, so a person's information grows with their own items and
#' judgements, and the standard errors are from the observed information.
#' The judgement tables name persons by the row names of \code{data}, by
#' the \code{id} column or vector, or by \code{"P1"}, \code{"P2"}, ...
#' when there are none; a person judged but absent from \code{data} has no
#' responses and is placed by the judgements. Threshold columns have no
#' meaning here and are ignored.
#'
#' A person has a finite location when the evidence points both ways: a
#' score inside its range, or at least one win and one loss (a rank above
#' someone and a rank below someone) among the persons being estimated.
#' A person lacking either is set aside as extreme and reported at the Warm
#' estimate from their responses, as \code{rasch} reports extreme persons,
#' or with no location if they have no responses; a judged group none of
#' whose members has responses is not on the test scale and is refused. The
#' invariance table compares each frame's separate placement of every
#' person, divided by the fitted unit and centred within each connected
#' block of the design over the persons both frames locate, with the
#' response placement, by a Wald contrast with Holm adjustment. There is no
#' likelihood ratio test in this mode: the separate model has a location
#' per person per frame, so its parameters grow with the persons and the
#' ratio is not chi-square. A person who fails is one whose judged work
#' does not match their responses.
#'
#' The unit of a judgement frame is estimated alongside one location per
#' person, so it carries the incidental-parameter bias of joint estimation
#' rather than the consistency of the item mode, where persons are
#' conditioned out. In simulations with 300 persons, 8 to 40 dichotomous
#' items and 4 to 20 comparisons per person the unit was 5 to 10 percent
#' high and its standard error understated, while the person locations and
#' their standard errors were calibrated. Fix the unit with \code{units}
#' when it is known. Too few judgements per person leave the unit without
#' a maximum, because some ordering of the persons the responses allow
#' agrees with every judgement; the fit then reports a runaway unit.
#'
#' @param data Persons-by-items response data, dichotomous or polytomous, as
#'   for \code{\link{rasch}}, or \code{NULL} to combine comparisons and
#'   rankings without responses. In the person mode, responses to the
#'   anchored items, coded from 0 to the item's top category.
#' @param comparisons Optional data frame of paired comparisons, one row
#'   each.
#' @param object_a,object_b,winner Column names in \code{comparisons}, as in
#'   \code{\link{btl}}. The winner is the object judged to have the higher
#'   location, for difficulty the harder one; it is given as the item name,
#'   or as \code{"item:k"} when both sides are thresholds of the same item.
#'   Winners equal to \code{"tie"} or \code{"draw"} are dropped.
#' @param threshold_a,threshold_b Optional column names in
#'   \code{comparisons} giving the threshold number of each side, \code{NA}
#'   for the item's location. Used only when the columns are present.
#' @param rankings Optional data frame of rankings in long form, one row per
#'   object placed: a ranking identifier, the item, and its rank within that
#'   ranking. Rankings may cover any subset of two or more objects.
#' @param ranking,item,rank Column names in \code{rankings}. Rank 1 is the
#'   object with the highest location.
#' @param threshold Optional column name in \code{rankings} giving the
#'   threshold number of each row, \code{NA} for the item's location. Used
#'   only when the column is present.
#' @param units Named vector giving the unit of the \code{comparisons} and
#'   \code{rankings} frames: \code{NA} (the default) to estimate, or \code{1}
#'   to fix at the response unit.
#' @param items Optional item columns to analyse, as in \code{\link{rasch}}.
#' @param na_codes As in \code{\link{rasch}}.
#' @param objects What the judges compare: \code{"items"} (the default),
#'   or \code{"persons"} to measure the persons from their responses and
#'   judgements of their work.
#' @param anchors Person mode only: the calibrated item thresholds, as a
#'   \code{\link{rasch}} fit, an item-mode \code{rasch_cj} fit, or a data
#'   frame with columns \code{item}, \code{k} and \code{tau} giving every
#'   threshold of every item in \code{data}.
#' @param id Person mode only: the name of a column of \code{data} holding
#'   the person identifiers the judgement tables use, or a vector of them,
#'   one per row. Defaults to the row names.
#' @param maxit,tol Newton iteration cap and convergence tolerance on the
#'   parameter scale.
#' @return An object of class \code{"rasch_cj"} with components
#'   \code{items} (item, location, se, and the separate calibration of each
#'   item-level object from each frame on the reference scale),
#'   \code{thresholds} (item, k, threshold, se), \code{objects} (the judged
#'   objects and the frames reaching them), \code{units} (frame, unit, se,
#'   and whether it was estimated), \code{invariance} (a list with the
#'   likelihood ratio test and the per-object table comparing each judgement
#'   frame with the reference frame), \code{anchors} (a data frame ready for
#'   \code{\link{rasch}}), \code{cov} (covariance of the thresholds),
#'   \code{cov_items} (covariance of the item locations), \code{loglik},
#'   \code{converged}, \code{iterations}, frame sizes in \code{n}, the
#'   \code{reference} frame, and \code{notes}.
#'
#'   In the person mode, \code{mode} is \code{"persons"} and the object
#'   holds \code{persons} (person, n_items, raw, max_raw, the combined
#'   location and se, the location and se from responses alone,
#'   \code{Inf} or \code{-Inf} at an extreme score, the centred location
#'   from each judgement frame alone on the test scale, whether each frame
#'   reaches the person, and whether the person was set aside as extreme),
#'   \code{units}, \code{invariance} (a list with the per-person table
#'   \code{persons}), \code{anchors} (the thresholds used), \code{cov}
#'   (covariance of the estimated locations), \code{loglik},
#'   \code{converged}, \code{iterations}, \code{n} and \code{notes}.
#' @references Bradley, R. A. and Terry, M. E. (1952). Rank analysis of
#'   incomplete block designs: I. The method of paired comparisons.
#'   Biometrika, 39, 324--345.
#'
#'   Humphry, S. M. and Andrich, D. (2008). Understanding the unit in the
#'   Rasch model. Journal of Applied Measurement, 9(3), 249--264.
#'
#'   Plackett, R. L. (1975). The analysis of permutations. Applied
#'   Statistics, 24(2), 193--202.
#' @seealso \code{\link{rasch}}, \code{\link{btl}},
#'   \code{\link{frame_invariance}}.
#' @examples
#' set.seed(1)
#' delta <- seq(-1.5, 1.5, length.out = 8)
#' names(delta) <- sprintf("I%02d", 1:8)
#' theta <- rnorm(300)
#' X <- sapply(delta, function(d) as.integer(runif(300) < plogis(theta - d)))
#' pairs <- t(combn(names(delta), 2))[sample(28, 200, replace = TRUE), ]
#' p_a <- plogis(0.6 * (delta[pairs[, 1]] - delta[pairs[, 2]]))
#' cj <- data.frame(a = pairs[, 1], b = pairs[, 2],
#'                  winner = ifelse(runif(200) < p_a, pairs[, 1], pairs[, 2]))
#' fit <- rasch_cj(X, comparisons = cj, object_a = "a", object_b = "b",
#'                 winner = "winner")
#' fit$units
#' fit$invariance$lr
#'
#' # persons from their responses and judgements of their work
#' pairs <- cbind(sample(300, 1200, replace = TRUE), sample(300, 1200, replace = TRUE))
#' pairs <- pairs[pairs[, 1] != pairs[, 2], ]
#' p_a <- plogis(0.8 * (theta[pairs[, 1]] - theta[pairs[, 2]]))
#' work <- data.frame(a = paste0("P", pairs[, 1]), b = paste0("P", pairs[, 2]))
#' work$winner <- ifelse(runif(nrow(work)) < p_a, work$a, work$b)
#' pfit <- rasch_cj(X, comparisons = work, object_a = "a", object_b = "b",
#'                  winner = "winner", objects = "persons", anchors = fit)
#' head(pfit$persons)
#' @export
rasch_cj <- function(data, comparisons = NULL, object_a = "object_a",
                     object_b = "object_b", winner = "winner",
                     threshold_a = "threshold_a", threshold_b = "threshold_b",
                     rankings = NULL, ranking = "ranking", item = "item",
                     rank = "rank", threshold = "threshold",
                     units = c(comparisons = NA, rankings = NA),
                     items = NULL, na_codes = -1,
                     objects = c("items", "persons"), anchors = NULL,
                     id = NULL, maxit = 200, tol = 1e-8) {
  objects <- match.arg(objects)
  has_resp <- !is.null(data)
  if (has_resp) .check_column_names(data)
  .check_controls(maxit, tol)
  if (is.null(comparisons) && is.null(rankings))
    stop("supply `comparisons`, `rankings`, or both; with responses alone use rasch()",
         call. = FALSE)
  if (objects == "items" && !is.null(anchors))
    stop("`anchors` is for the person mode: set `objects = \"persons\"`",
         call. = FALSE)
  if (objects == "items" && !is.null(id))
    stop("`id` is for the person mode: set `objects = \"persons\"`", call. = FALSE)
  if (!has_resp && (is.null(comparisons) || is.null(rankings)))
    stop("without responses, supply both `comparisons` and `rankings`; ",
         "for one judgement source alone use btl()", call. = FALSE)
  if (!is.null(comparisons) && !is.data.frame(comparisons))
    stop("`comparisons` must be a data frame", call. = FALSE)
  if (!is.null(rankings) && !is.data.frame(rankings))
    stop("`rankings` must be a data frame", call. = FALSE)
  if (!is.numeric(units) && !is.logical(units))
    stop("`units` must be a named vector of NA or 1", call. = FALSE)
  u_spec <- c(comparisons = NA_real_, rankings = NA_real_)
  if (!is.null(names(units))) {
    bad <- setdiff(names(units), names(u_spec))
    if (length(bad))
      stop("`units` names must be \"comparisons\" or \"rankings\"", call. = FALSE)
    u_spec[names(units)] <- as.numeric(units)
  } else if (length(units) == 2L) {
    u_spec[] <- as.numeric(units)
  } else stop("`units` must be a named vector of NA or 1", call. = FALSE)
  if (any(!is.na(u_spec) & u_spec != 1))
    stop("a fixed unit must be 1; other values rescale the reference frame",
         call. = FALSE)
  if (objects == "persons")
    return(.cj_persons(data, anchors, id, comparisons, object_a, object_b,
                       winner, threshold_a, threshold_b, rankings, ranking,
                       item, rank, threshold, u_spec, items, na_codes, maxit,
                       tol, match.call()))
  notes <- character(0)
  ref <- if (has_resp) "responses" else "comparisons"
  if (!has_resp) u_spec[["comparisons"]] <- 1

  if (!has_resp) {
    # judgement-only fit: the comparisons are the reference frame and the
    # items are whatever the two sources name, each a single location
    item_names <- unique(c(.role_text_values(comparisons[[object_a]]),
                           .role_text_values(comparisons[[object_b]]),
                           .role_text_values(rankings[[item]])))
    item_names <- item_names[!is.na(item_names)]
    X <- matrix(NA_integer_, 0L, length(item_names),
                dimnames = list(NULL, item_names))
    constant <- character(0)
    m <- rep(1L, length(item_names))
  } else {
    X <- if (is.data.frame(data)) data else as.data.frame(data)
    if (!is.null(items)) X <- X[, items, drop = FALSE]
    prep <- .prepare_X(X, na_codes = na_codes)
    # rasch() drops an item every person answered the same way, because the
    # conditional likelihood cannot place it. Here the judgements can, so
    # such items are kept with their responses set aside, and required to be
    # reached by a judgement frame below.
    raw_names <- colnames(as.matrix(X))
    if (is.null(raw_names)) raw_names <- sprintf("I%02d", seq_len(ncol(X)))
    constant <- setdiff(raw_names, colnames(prep$X))
    Xp <- prep$X
    notes <- c(notes, grep("^dropped constant item", prep$notes, value = TRUE,
                           invert = TRUE))
    if (length(constant)) {
      for (nm in constant) {
        Xp <- cbind(Xp, NA_integer_)
        colnames(Xp)[ncol(Xp)] <- nm
      }
      Xp <- Xp[, raw_names[raw_names %in% colnames(Xp)], drop = FALSE]
      notes <- c(notes, paste0("item(s) with no response variation, located ",
                               "by the judgements alone: ",
                               paste(constant, collapse = ", ")))
    }
    X <- Xp
    item_names <- colnames(X)
    # .prepare_X has already rescored gaps to consecutive categories from 0
    m <- vapply(seq_len(ncol(X)), function(i) {
      z <- X[!is.na(X[, i]), i]
      if (!length(z)) 1L else as.integer(max(z))
    }, 1L)
  }
  I <- length(item_names)
  if (I < 2L) stop("at least two items are required", call. = FALSE)
  o <- .cj_cml_prep(X, m)
  P <- o$P; idx <- o$idx
  n_resp <- o$n_persons

  # judgement sources resolved to object keys, then to one object table
  cmp <- if (!is.null(comparisons))
    .cj_comparison_keys(comparisons, object_a, object_b, winner, threshold_a,
                        threshold_b, item_names, if (has_resp) m else NULL, notes)
  else list(win = character(0), lose = character(0), n = 0L, notes = notes)
  notes <- cmp$notes
  rkl <- if (!is.null(rankings))
    .cj_ranking_keys(rankings, ranking, item, rank, threshold, item_names,
                     if (has_resp) m else NULL, notes)
  else list(rk = list(), n = 0L, notes = notes)
  notes <- rkl$notes
  has_cmp <- cmp$n > 0L; has_rk <- rkl$n > 0L
  if (!is.null(comparisons) && !has_cmp)
    stop("`comparisons` has no usable rows", call. = FALSE)
  if (!is.null(rankings) && !has_rk)
    stop("`rankings` has no usable rows", call. = FALSE)

  keys <- unique(c(cmp$win, cmp$lose, unlist(rkl$rk, use.names = FALSE)))
  is_thr <- grepl(":[0-9]+$", keys)
  obj_item <- ifelse(is_thr, sub(":[0-9]+$", "", keys), keys)
  obj_thr <- rep(NA_integer_, length(keys))
  obj_thr[is_thr] <- as.integer(sub("^.*:", "", keys[is_thr]))
  ord <- order(match(obj_item, item_names), obj_thr, na.last = FALSE)
  keys <- keys[ord]; obj_item <- obj_item[ord]; obj_thr <- obj_thr[ord]
  n_obj <- length(keys)
  # object matrix: an item's location is the mean of its thresholds
  O <- matrix(0, n_obj, P)
  for (j in seq_len(n_obj)) {
    i <- match(obj_item[j], item_names)
    if (is.na(obj_thr[j])) O[j, idx[[i]]] <- 1 / m[i]
    else O[j, idx[[i]][obj_thr[j]]] <- 1
  }
  W <- matrix(0, n_obj, n_obj)
  if (has_cmp)
    W <- matrix(as.numeric(table(factor(match(cmp$win, keys), levels = seq_len(n_obj)),
                                 factor(match(cmp$lose, keys), levels = seq_len(n_obj)))),
                n_obj, n_obj)
  rk <- lapply(rkl$rk, function(v) match(v, keys))

  # connectivity over all frames: two items are linked when a pattern
  # answers both, a comparison pairs their objects, or a ranking places both
  obj_of_item <- match(obj_item, item_names)
  adj <- matrix(FALSE, I, I)
  if (nrow(o$A)) adj <- adj | (crossprod(o$A * 1) > 0)
  adj_obj <- W + t(W) > 0
  for (v in rk) adj_obj[v, v] <- TRUE
  for (j in which(rowSums(adj_obj) > 0))
    adj[obj_of_item[j], obj_of_item[adj_obj[j, ]]] <- TRUE
  adj <- adj | t(adj); diag(adj) <- TRUE
  if (length(constant)) {
    unplaced <- constant[!(constant %in% obj_item)]
    if (length(unplaced))
      stop("item(s) with no response variation and no judgements cannot be ",
           "located: ", paste(unplaced, collapse = ", "), call. = FALSE)
  }
  comp_items <- .cj_components(adj)
  if (max(comp_items) > 1L)
    stop("the items are not connected across the frames; no source links ",
         paste(item_names[comp_items != 1L], collapse = ", "), " to ",
         paste(item_names[comp_items == 1L], collapse = ", "), call. = FALSE)

  # frames and their unit parameters: theta = (beta, log units of free frames)
  frames <- c(if (has_resp) "responses", if (has_cmp) "comparisons",
              if (has_rk) "rankings")
  free <- c(if (has_cmp) is.na(u_spec[["comparisons"]]),
            if (has_rk) is.na(u_spec[["rankings"]]))
  names(free) <- setdiff(frames, "responses")
  cw <- rep(1 / m, m)                      # mean item location is zero
  B <- .cj_basis(P, cw); nb <- ncol(B); nu <- sum(free)
  u_pos <- integer(0)
  if (nu) u_pos <- nb + seq_len(nu)
  names(u_pos) <- names(free)[free]
  unit_of <- function(th, f) if (isTRUE(free[[f]])) exp(th[u_pos[[f]]]) else 1

  bt_ll <- function(lam, u) .cj_bt_ll(lam, u, W)
  bt_parts <- function(lam, u) .cj_bt_parts(lam, u, W)
  pl_ll <- function(lam, u) .cj_pl_ll(lam, u, rk)
  pl_parts <- function(lam, u) .cj_pl_parts(lam, u, rk)

  fn <- function(th) {
    psi <- as.vector(B %*% th[seq_len(nb)]); lam <- as.vector(O %*% psi)
    ll <- if (has_resp) .cj_cml_ll(psi, o) else 0
    if (has_cmp) ll <- ll + bt_ll(lam, unit_of(th, "comparisons"))
    if (has_rk) ll <- ll + pl_ll(lam, unit_of(th, "rankings"))
    ll
  }
  # assemble gradient and Hessian in theta from the per-frame blocks, the
  # judgement blocks given in (lambda, log unit) and mapped through O,
  # dropping the unit row when the unit is fixed
  assemble <- function(th) {
    psi <- as.vector(B %*% th[seq_len(nb)]); lam <- as.vector(O %*% psi)
    np <- nb + nu
    g <- numeric(np); H <- matrix(0, np, np)
    g_p <- numeric(P); H_pp <- matrix(0, P, P)
    if (has_resp) {
      pp <- .cj_cml_p(psi, o)
      g_p <- .cj_cml_grad(psi, o, pp); H_pp <- .cj_cml_hess(psi, o, pp)
    }
    add <- function(parts, f) {
      jo <- seq_len(n_obj)
      g_p <<- g_p + as.vector(crossprod(O, parts$g[jo]))
      H_pp <<- H_pp + crossprod(O, parts$H[jo, jo] %*% O)
      if (isTRUE(free[[f]])) {
        k <- u_pos[[f]]
        g[k] <<- parts$g[n_obj + 1L]
        cross <- as.vector(crossprod(B, crossprod(O, parts$H[jo, n_obj + 1L])))
        H[seq_len(nb), k] <<- H[seq_len(nb), k] + cross
        H[k, seq_len(nb)] <<- H[k, seq_len(nb)] + cross
        H[k, k] <<- parts$H[n_obj + 1L, n_obj + 1L]
      }
    }
    if (has_cmp) add(bt_parts(lam, unit_of(th, "comparisons")), "comparisons")
    if (has_rk) add(pl_parts(lam, unit_of(th, "rankings")), "rankings")
    g[seq_len(nb)] <- as.vector(crossprod(B, g_p))
    H[seq_len(nb), seq_len(nb)] <- H[seq_len(nb), seq_len(nb)] +
      crossprod(B, H_pp %*% B)
    list(g = g, H = H)
  }
  gr <- function(th) assemble(th)$g
  he <- function(th) assemble(th)$H

  # separate calibrations: starting values, and the invariance alternative.
  # A judgement frame is fitted over the objects it reaches, with one
  # constraint per connected block of its own design.
  sep <- list(); judged <- list(); blocks <- list()
  # the response frame alone cannot place a constant item, so its
  # thresholds are held at zero there and its objects sit out the contrasts
  inf_thr <- which(rep(!(item_names %in% constant), m))
  B_r <- matrix(0, P, max(length(inf_thr) - 1L, 0L))
  B_r[inf_thr, ] <- .cj_basis(length(inf_thr), cw[inf_thr])
  sep$responses <- if (has_resp && nrow(o$A))
    .cj_fit_alone(function(p) .cj_cml_ll(p, o),
                  function(p) { pp <- .cj_cml_p(p, o)
                    list(g = .cj_cml_grad(p, o, pp), H = .cj_cml_hess(p, o, pp)) },
                  B_r, maxit = maxit, tol = tol)
  else NULL
  frame_alone <- function(ll, parts, adj_f) {
    jo <- which(rowSums(adj_f) > 0)
    comp <- .cj_components(adj_f[jo, jo, drop = FALSE])
    Bf <- .cj_basis_groups(comp)
    fit <- .cj_fit_alone(function(l) { z <- numeric(n_obj); z[jo] <- l; ll(z, 1) },
                         function(l) { z <- numeric(n_obj); z[jo] <- l
                           p <- parts(z, 1)
                           list(g = p$g[jo], H = p$H[jo, jo, drop = FALSE]) },
                         Bf, maxit = maxit, tol = tol)
    list(fit = fit, jo = jo, comp = comp)
  }
  if (has_cmp) {
    a <- frame_alone(bt_ll, bt_parts, W + t(W) > 0)
    sep$comparisons <- a$fit; judged$comparisons <- a$jo; blocks$comparisons <- a$comp
  }
  if (has_rk) {
    adj_r <- matrix(FALSE, n_obj, n_obj)
    for (v in rk) adj_r[v, v] <- TRUE
    a <- frame_alone(pl_ll, pl_parts, adj_r)
    sep$rankings <- a$fit; judged$rankings <- a$jo; blocks$rankings <- a$comp
  }

  # start from the response calibration where it exists and is finite; a
  # separate calibration can be infinite for a threshold no person
  # discriminates, and such thresholds start at zero for the judgements to
  # place. Without responses, start from the reference judgement frame.
  start_p <- rep(0, P)
  if (!is.null(sep$responses) && sep$responses$converged) {
    start_p <- sep$responses$par
    start_p[!is.finite(start_p)] <- 0
  } else if (!has_resp) {
    start_p[obj_of_item[judged[[ref]]]] <- sep[[ref]]$par
  }
  start_p <- start_p - sum(cw * start_p) / sum(cw)
  th0 <- c(as.vector(solve(crossprod(B), crossprod(B, start_p))), rep(0, nu))
  fit <- .cj_newton(th0, fn, gr, he, maxit, tol)
  th <- fit$par
  psi <- as.vector(B %*% th[seq_len(nb)])
  covth <- tryCatch(solve(-fit$H), error = function(e) matrix(NA_real_, nb + nu, nb + nu))
  cov_p <- B %*% covth[seq_len(nb), seq_len(nb), drop = FALSE] %*% t(B)
  L <- matrix(0, I, P)                     # item locations from thresholds
  for (i in seq_len(I)) L[i, idx[[i]]] <- 1 / m[i]
  delta <- as.vector(L %*% psi); names(delta) <- item_names
  cov_d <- L %*% cov_p %*% t(L)
  dimnames(cov_d) <- list(item_names, item_names)
  thr_names <- unlist(lapply(seq_len(I), function(i)
    paste0(item_names[i], ":", seq_len(m[i]))), use.names = FALSE)
  dimnames(cov_p) <- list(thr_names, thr_names)
  se_d <- sqrt(pmax(diag(cov_d), 0)); se_p <- sqrt(pmax(diag(cov_p), 0))
  if (!fit$converged) {
    notes <- c(notes, "estimation did not converge; standard errors and tests withheld")
    se_d[] <- NA_real_; cov_d[] <- NA_real_; se_p[] <- NA_real_; cov_p[] <- NA_real_
  }
  for (f in names(free)[free]) {
    u <- unit_of(th, f)
    if (u < 1e-3)
      notes <- c(notes, sprintf(paste0("the %s unit collapsed towards zero: ",
        "the judgements carry no information about the objects, or their ",
        "orientation is reversed (the winner or first rank should be the ",
        "object with the higher location)"), f))
  }

  units_tab <- data.frame(frame = frames, unit = 1, se = NA_real_,
                          estimated = FALSE, stringsAsFactors = FALSE)
  for (f in names(free)) {
    r <- match(f, units_tab$frame)
    units_tab$unit[r] <- unit_of(th, f)
    units_tab$estimated[r] <- isTRUE(free[[f]])
    if (isTRUE(free[[f]]) && fit$converged) {
      k <- u_pos[[f]]
      # delta method from the log scale
      units_tab$se[r] <- units_tab$unit[r] * sqrt(max(covth[k, k], 0))
    }
  }

  # invariance: likelihood ratio against separate locations per frame, and
  # per-object Wald contrasts of each judgement frame with the reference
  ll_sep <- sum(vapply(sep, function(s) s$ll, 0))
  df <- sum(vapply(setdiff(names(sep), ref), function(f) sep[[f]]$n_free, 0)) - nu
  lr <- data.frame(statistic = 2 * (ll_sep - fit$ll), df = as.integer(df))
  lr$p <- if (fit$converged && lr$df > 0 && is.finite(lr$statistic))
    stats::pchisq(max(lr$statistic, 0), lr$df, lower.tail = FALSE) else NA_real_

  item_tab <- data.frame(item = item_names, location = delta, se = se_d,
                         stringsAsFactors = FALSE, row.names = NULL)
  thr_tab <- data.frame(item = rep(item_names, m),
                        k = unlist(lapply(m, seq_len), use.names = FALSE),
                        threshold = psi, se = se_p,
                        stringsAsFactors = FALSE, row.names = NULL)
  # the reference calibration of every object, on the object scale
  ref_lam <- NULL
  if (ref == "responses" && !is.null(sep$responses)) {
    ref_lam <- as.vector(O %*% sep$responses$par)
    ref_cov <- O %*% sep$responses$cov %*% t(O)
    ref_jo <- which(!(obj_item %in% constant))
  } else if (ref != "responses") {
    ref_lam <- numeric(n_obj); ref_lam[judged[[ref]]] <- sep[[ref]]$par
    ref_cov <- matrix(0, n_obj, n_obj)
    ref_cov[judged[[ref]], judged[[ref]]] <- sep[[ref]]$cov
    ref_jo <- judged[[ref]]
  }
  inv_tab <- NULL
  for (f in setdiff(names(free), ref)) {
    u <- unit_of(th, f); s <- sep[[f]]; jo <- judged[[f]]; comp <- blocks[[f]]
    # centre within each block of the frame's design
    C <- matrix(0, length(jo), length(jo))
    for (g in unique(comp)) { rows <- comp == g; C[rows, rows] <- -1 / sum(rows) }
    diag(C) <- diag(C) + 1
    l_f <- as.vector(C %*% (s$par / u))
    at_item <- is.na(obj_thr[jo])
    col <- rep(NA_real_, I)
    col[match(obj_item[jo][at_item], item_names)] <- l_f[at_item]
    item_tab[[paste0("location_", f)]] <- col
    if (!is.null(ref_lam) && fit$converged) {
      both <- jo %in% ref_jo
      l_r <- as.vector(C %*% ref_lam[jo])
      diff <- l_f - l_r
      V <- C %*% (ref_cov[jo, jo, drop = FALSE] + s$cov / u^2) %*% t(C)
      z <- diff / sqrt(pmax(diag(V), 0))
      z[!is.finite(z) | !both] <- NA_real_
      p <- 2 * stats::pnorm(-abs(z))
      inv_tab <- rbind(inv_tab, data.frame(
        frame = f, item = obj_item[jo], threshold = obj_thr[jo],
        reference = l_r, judgements = l_f,
        difference = diff, se = sqrt(pmax(diag(V), 0)), z = z, p = p,
        p_adj = stats::p.adjust(p, "holm"), stringsAsFactors = FALSE,
        row.names = NULL))
    }
  }
  if (ref == "responses" && !is.null(sep$responses)) {
    item_tab$location_responses <- as.vector(L %*% sep$responses$par)
  } else if (ref != "responses") {
    col <- rep(NA_real_, I)
    col[match(obj_item[judged[[ref]]], item_names)] <- sep[[ref]]$par
    item_tab[[paste0("location_", ref)]] <- col
  }

  obj_tab <- data.frame(object = keys, item = obj_item, threshold = obj_thr,
                        location = as.vector(O %*% psi),
                        comparisons = seq_len(n_obj) %in% judged$comparisons,
                        rankings = seq_len(n_obj) %in% judged$rankings,
                        stringsAsFactors = FALSE, row.names = NULL)
  anchors <- data.frame(item = thr_tab$item, k = thr_tab$k, tau = psi,
                        stringsAsFactors = FALSE, row.names = NULL)

  structure(list(items = item_tab, thresholds = thr_tab, objects = obj_tab,
                 units = units_tab,
                 invariance = list(lr = lr, items = inv_tab),
                 anchors = anchors, cov = cov_p, cov_items = cov_d,
                 loglik = fit$ll, loglik_separate = ll_sep,
                 converged = fit$converged, iterations = fit$iterations,
                 n = c(persons = n_resp, comparisons = cmp$n, rankings = rkl$n),
                 reference = ref, notes = notes, call = match.call()),
            class = "rasch_cj")
}

#' @export
print.rasch_cj <- function(x, ...) {
  if (identical(x$mode, "persons")) {
    ps <- x$persons
    cat(sprintf(paste0("Combined person measurement: %d persons (%d with ",
                       "responses to %d anchored items, %d extreme) from %d ",
                       "comparisons, %d rankings\n"),
                nrow(ps), x$n[["persons"]], length(unique(x$anchors$item)),
                sum(ps$extreme), x$n[["comparisons"]], x$n[["rankings"]]))
  } else {
    n_thr <- if (is.null(x$thresholds)) nrow(x$items) else nrow(x$thresholds)
    cat(sprintf(paste0("Combined calibration: %d items%s from %d informative ",
                       "persons, %d comparisons, %d rankings\n"),
                nrow(x$items),
                if (n_thr > nrow(x$items)) sprintf(" (%d thresholds)", n_thr) else "",
                x$n[["persons"]], x$n[["comparisons"]], x$n[["rankings"]]))
  }
  cat(sprintf("Full likelihood: %s in %d iterations; log-likelihood %.2f\n",
              if (x$converged) "converged" else "NOT converged", x$iterations,
              x$loglik))
  ref <- if (is.null(x$reference)) "responses" else x$reference
  for (r in which(x$units$frame != ref)) {
    cat(sprintf("Unit of %s relative to %s: %.3f%s\n",
                x$units$frame[r], ref, x$units$unit[r],
                if (x$units$estimated[r] && is.finite(x$units$se[r]))
                  sprintf(" (se %.3f)", x$units$se[r]) else if (!x$units$estimated[r])
                  " (fixed)" else ""))
  }
  lr <- x$invariance$lr
  if (!is.null(lr) && is.finite(lr$p))
    cat(sprintf("Invariance across frames: LR %.2f on %d df, p = %s\n",
                lr$statistic, lr$df, .fmt_p(lr$p)))
  if (identical(x$mode, "persons")) {
    tab <- x$invariance$persons
    if (!is.null(tab)) {
      for (f in unique(tab$frame)) {
        r <- tab$frame == f & !is.na(tab$p)
        hit <- r & tab$p_adj < 0.05
        cat(sprintf(paste0("Persons placed differently by the %s (%d with a ",
                           "contrast): %d at p < 0.05, %d at Holm p < 0.05%s\n"),
                    f, sum(r), sum(tab$p[r] < 0.05), sum(hit),
                    if (any(hit)) paste0(": ", paste(tab$person[hit], collapse = ", "))
                    else ""))
      }
    }
  } else {
    tab <- x$invariance$items
    if (!is.null(tab)) {
      hit <- !is.na(tab$p_adj) & tab$p_adj < 0.05
      lab <- if ("threshold" %in% names(tab))
        ifelse(is.na(tab$threshold), tab$item, paste0(tab$item, ":", tab$threshold))
      else tab$item
      cat(sprintf("Objects differing between frames (Holm p < 0.05): %s\n",
                  if (any(hit)) paste(unique(lab[hit]), collapse = ", ") else "none"))
    }
  }
  for (n in x$notes) cat("Note: ", n, "\n", sep = "")
  invisible(x)
}
