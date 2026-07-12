# rasch :: extended frame of reference for paired comparisons
# ===========================================================================
# The extended frame of reference model of Humphry (2005) and Humphry and
# Andrich (2008) states that the unit of the latent scale is a property of
# the frame of measurement, not a universal constant. rasch_efrm() fits that
# model for the persons-by-items design; btl_efrm() is this package's
# extension of the same idea to the Bradley-Terry-Luce family of paired
# comparisons (Bradley and Terry 1952; Luce 1959), where the frame is a
# judge-panel by object-set cell.
#
# Objects k are partitioned into sets s(k) (S sets), and judges j into panels
# g(j) (G panels). Each object has a true common-scale value
#
#     v_k = alpha_{s(k)} * beta_k + kappa_{s(k)},
#
# with beta_k the within-set (frame-unit) calibration location, alpha_s > 0
# the set unit, and kappa_s the set origin. A comparison judged in panel g
# carries the panel unit phi_g, and:
#
#   * WITHIN a set s (both objects in s):
#         logit P(a beats b) = phi_g * (beta_a - beta_b).
#     The origin kappa_s cancels, and the set unit alpha_s is confounded with
#     the spread of beta -- exactly as the item-set unit is absorbed in the
#     within-frame stage of the Rasch EFRM -- so it is not a parameter here.
#
#   * ACROSS sets (a in A, b in B):
#         logit P(a beats b) = phi_g * (v_a - v_b)
#                            = phi_g * (alpha_A beta_a - alpha_B beta_b
#                                        + kappa_A - kappa_B).
#     The cross-set comparisons place the two sets on one common scale and so
#     identify alpha and kappa.
#
# The lineage of a frame-dependent unit for comparative judgement is
# Thurstone's (1927) varying discriminal dispersion and the varying-precision
# paired-comparison models catalogued by David (1988); the measurement-unit
# reading of it is Humphry's. The paired-comparison form fitted here is this
# package's extension, stated for dichotomous winner data.
#
# Estimation is in two conditional stages, mirroring rasch_efrm():
#
# Stage 1 (within frames): for each set the within-set comparisons, pooled
# over panels, fit the bilinear model logit = rho_{gs} (b_a - b_b) with
# b sum-zero and one reference panel fixed at rho = 1. This is the same
# constrained bilinear maximisation the Rasch EFRM performs in its stage 1;
# the ratios rho_{gs} = phi_g / phi_{ref(s)} estimate the panel units up to
# the set's reference panel, and are reconciled across sets by a
# precision-weighted least squares over the panel-by-set linking graph, then
# normalised to geometric mean one over panels. The reconciled reference-panel
# units put every set's b on the common panel scale, giving beta.
#
# Stage 2 (linking sets): with beta-hat and phi-hat fixed, the cross-set
# comparisons are a low-dimensional maximum likelihood in (log alpha, kappa)
# for the non-reference sets, solved by Newton with the analytic gradient and
# Hessian. The key theoretical point, and the reason this design is worth
# stating, is that the linking uses only comparison OUTCOMES: no distributional
# assumption about the objects is made, so the set units are identified WITHIN
# the conditional (person-free) framework. This is unlike the persons-by-items
# EFRM, whose item-set units are identified only from the person side (their
# distribution), a genuinely distributional step. The paired-comparison design
# supplies its own conditional link and needs no such assumption.
# ===========================================================================

# connected components of an undirected graph on 1..n given a two-column
# integer edge matrix; union-find, as used across the package's linking code
.btlef_components <- function(n, edges) {
  parent <- seq_len(n)
  find <- function(a) { while (parent[a] != a) a <- parent[a]; a }
  if (length(edges)) for (r in seq_len(nrow(edges))) {
    ra <- find(edges[r, 1]); rb <- find(edges[r, 2])
    if (ra != rb) parent[ra] <- rb
  }
  vapply(seq_len(n), find, 1L)
}

# Stage 1 bilinear solve for ONE frame set (or for the pooled single-unit
# model, called with one panel). Objects are indexed 1..K; panel is a
# character vector per comparison; the most-used panel is the reference
# (rho = 1). Parameters are the sum-zero location contrasts and the free
# panels' log discrimination. Fisher scoring with a step-halving line search
# gives the point estimate; a judge-clustered Godambe sandwich gives the
# covariance, exactly as in btl().
.btlef_stage1 <- function(ia, ib, y, panel, jd, K, maxit, tol) {
  R <- length(ia)
  pcount <- table(panel)
  present <- names(pcount)
  ref <- present[which.max(as.integer(pcount))]     # reference panel: rho = 1
  free <- setdiff(present, ref); Gf <- length(free)
  pf <- match(panel, free)                           # free-panel index, NA on ref
  B <- rbind(diag(K - 1L), rep(-1, K - 1L))          # K x (K-1) sum-zero map
  Bd <- B[ia, , drop = FALSE] - B[ib, , drop = FALSE]
  np <- (K - 1L) + Gf

  eval_th <- function(th) {
    bfree <- th[seq_len(K - 1L)]
    beta <- as.numeric(B %*% bfree)
    rho <- rep(1, R)
    if (Gf) {
      rr <- exp(th[(K - 1L) + seq_len(Gf)])
      ok <- !is.na(pf); rho[ok] <- rr[pf[ok]]
    }
    d <- beta[ia] - beta[ib]
    p <- plogis(rho * d)
    ll <- sum(ifelse(y == 1, log(pmax(p, 1e-300)), log(pmax(1 - p, 1e-300))))
    list(beta = beta, rho = rho, d = d, p = p, ll = ll)
  }
  design <- function(cur) {
    J <- Bd * cur$rho
    if (Gf) {
      Jr <- matrix(0, R, Gf)
      for (h in seq_len(Gf)) {
        sel <- which(pf == h); Jr[sel, h] <- cur$rho[sel] * cur$d[sel]
      }
      J <- cbind(J, Jr)
    }
    J
  }

  theta <- numeric(np); cur <- eval_th(theta)
  for (it in seq_len(maxit)) {
    J <- design(cur); u <- y - cur$p; av <- cur$p * (1 - cur$p)
    g <- crossprod(J, u); Fi <- crossprod(J, J * av)
    step <- tryCatch(solve(Fi, g),
                     error = function(e) solve(Fi + diag(1e-8, np), g))
    lam <- 1; moved <- FALSE
    for (half in 1:30) {
      cand <- theta + lam * as.numeric(step); c2 <- eval_th(cand)
      if (is.finite(c2$ll) && c2$ll >= cur$ll - 1e-12) {
        theta <- cand; cur <- c2; moved <- TRUE; break
      }
      lam <- lam / 2
    }
    if (!moved) break
    if (max(abs(lam * step)) < tol) break
  }

  # judge-clustered Godambe sandwich (unclustered when every judge appears once)
  J <- design(cur); u <- y - cur$p; av <- cur$p * (1 - cur$p)
  Fi <- crossprod(J, J * av)
  bread <- tryCatch(solve(Fi), error = function(e) solve(Fi + diag(1e-8, np)))
  Sr <- J * u
  Sc <- rowsum(Sr, jd)
  cov_theta <- bread %*% crossprod(Sc) %*% bread
  conv <- max(abs(crossprod(J, u))) < 1e-4

  cov_bb <- B %*% cov_theta[seq_len(K - 1L), seq_len(K - 1L), drop = FALSE] %*% t(B)
  se_beta <- sqrt(pmax(diag(cov_bb), 0))
  rho_p <- setNames(rep(1, length(present)), present)
  cov_lrho <- matrix(0, Gf, Gf, dimnames = list(free, free))
  if (Gf) {
    li <- (K - 1L) + seq_len(Gf)
    rho_p[free] <- exp(theta[li])
    cov_lrho <- cov_theta[li, li, drop = FALSE]
    dimnames(cov_lrho) <- list(free, free)
  }
  list(beta = cur$beta, se_beta = se_beta, p = cur$p, ll = cur$ll,
       ref = ref, panels = present, rho = rho_p, free = free,
       cov_lrho = cov_lrho, converged = conv)
}

# Reconcile the per-set panel ratios into one set of panel units phi with
# geometric mean one, by generalised least squares on the panel graph. Each
# set contributes the observations log rho_{gs} = log phi_g - log phi_{ref(s)}
# for its free panels g, with the within-set covariance of those log-ratios
# carried across from stage 1; the observations are independent between sets,
# so the full observation covariance is block-diagonal by set. GLS on this
# gives both the precision-weighted point estimate (the reconciliation over
# the sets where a panel appears) and correctly correlated standard errors.
# Errors informatively when the panels are not connected through shared sets.
.btlef_reconcile_phi <- function(panels_u, blocks) {
  G <- length(panels_u)
  if (G == 1L)
    return(list(phi = setNames(1, panels_u),
                se_log_phi = setNames(NA_real_, panels_u),
                lphi = setNames(0, panels_u)))
  # flatten the per-set blocks into one observation vector, design and
  # block-diagonal covariance
  y <- numeric(0); pan <- ref <- character(0)
  Cov <- matrix(0, 0, 0)
  for (bk in blocks) {
    if (!length(bk$free)) next
    idx <- length(y) + seq_along(bk$free)
    y <- c(y, bk$lrho[bk$free]); pan <- c(pan, bk$free)
    ref <- c(ref, rep(bk$ref, length(bk$free)))
    cb <- bk$cov[bk$free, bk$free, drop = FALSE]
    cb <- cb + diag(1e-10, nrow(cb))                       # numerical floor
    Z <- matrix(0, nrow(Cov) + nrow(cb), ncol(Cov) + ncol(cb))
    if (nrow(Cov)) Z[seq_len(nrow(Cov)), seq_len(ncol(Cov))] <- Cov
    Z[nrow(Cov) + seq_len(nrow(cb)), ncol(Cov) + seq_len(ncol(cb))] <- cb
    Cov <- Z
  }
  if (!length(y))
    stop("the panels cannot be linked: no set contains comparisons from more ",
         "than one panel, so the panel units phi are unidentified")
  ei <- cbind(match(pan, panels_u), match(ref, panels_u))
  comp <- .btlef_components(G, ei)
  if (length(unique(comp)) > 1L)
    stop("the panel-by-set graph is not connected; panel units (phi) are ",
         "unidentified between: ",
         paste(tapply(panels_u, comp, paste, collapse = "+"), collapse = " | "))

  g0 <- panels_u[1]                                        # arbitrary anchor
  cols <- setdiff(panels_u, g0)
  X <- matrix(0, length(y), length(cols))
  for (r in seq_along(y)) {
    cp <- match(pan[r], cols); if (!is.na(cp)) X[r, cp] <- X[r, cp] + 1
    cr <- match(ref[r], cols); if (!is.na(cr)) X[r, cr] <- X[r, cr] - 1
  }
  W <- solve(Cov)                                          # GLS weight
  XtW <- t(X) %*% W
  covred <- solve(XtW %*% X)
  bred <- covred %*% (XtW %*% y)
  lphi <- setNames(numeric(G), panels_u); lphi[cols] <- bred
  cov_full <- matrix(0, G, G, dimnames = list(panels_u, panels_u))
  cov_full[cols, cols] <- covred
  A <- diag(G) - matrix(1 / G, G, G)                       # centre to geo-mean 1
  lphi_c <- as.numeric(A %*% lphi)
  cov_c <- A %*% cov_full %*% t(A)
  list(phi = setNames(exp(lphi_c), panels_u),
       se_log_phi = setNames(sqrt(pmax(diag(cov_c), 0)), panels_u),
       lphi = setNames(lphi_c, panels_u))
}

# Stage 2: cross-set linking. With the frame locations beta and panel units
# phi held fixed, estimate (log alpha, kappa) for the non-reference sets by
# Newton on the cross-set comparison likelihood. Standard errors are the
# inverse observed information, conditional on stage 1 (the stage-1
# uncertainty is not propagated -- see the roxygen note).
.btlef_stage2 <- function(a, b, y, phg, sa, sb, bhat, sets_u, maxit, tol) {
  S <- length(sets_u); free <- sets_u[-1L]; nf <- S - 1L; np <- 2L * nf
  ba <- bhat[a]; bb <- bhat[b]
  fa <- match(sa, free); fb <- match(sb, free)             # NA on the reference set
  R <- length(y)

  eval_th <- function(th) {
    la <- th[seq_len(nf)]; kap <- th[nf + seq_len(nf)]
    alpha <- setNames(rep(1, S), sets_u); kappa <- setNames(rep(0, S), sets_u)
    alpha[free] <- exp(la); kappa[free] <- kap
    va <- alpha[sa] * ba + kappa[sa]; vb <- alpha[sb] * bb + kappa[sb]
    p <- plogis(phg * (va - vb))
    ll <- sum(ifelse(y == 1, log(pmax(p, 1e-300)), log(pmax(1 - p, 1e-300))))
    list(alpha = alpha, kappa = kappa, p = p, ll = ll)
  }
  # derivative design D (R x np): columns log alpha (free sets), then kappa
  design <- function(cur) {
    D <- matrix(0, R, np)
    for (j in seq_len(nf)) {
      s <- free[j]; aj <- cur$alpha[[s]]
      onA <- !is.na(fa) & fa == j; onB <- !is.na(fb) & fb == j
      D[, j] <- phg * (onA * aj * ba - onB * aj * bb)       # d eta / d log alpha_s
      D[, nf + j] <- phg * (onA - onB)                      # d eta / d kappa_s
    }
    D
  }

  theta <- numeric(np); cur <- eval_th(theta)
  for (it in seq_len(maxit)) {
    D <- design(cur); u <- y - cur$p; av <- cur$p * (1 - cur$p)
    g <- crossprod(D, u); Fi <- crossprod(D, D * av)
    step <- tryCatch(solve(Fi, g),
                     error = function(e) solve(Fi + diag(1e-8, np), g))
    lam <- 1; moved <- FALSE
    for (half in 1:30) {
      cand <- theta + lam * as.numeric(step); c2 <- eval_th(cand)
      if (is.finite(c2$ll) && c2$ll >= cur$ll - 1e-12) {
        theta <- cand; cur <- c2; moved <- TRUE; break
      }
      lam <- lam / 2
    }
    if (!moved) break
    if (max(abs(lam * step)) < tol) break
  }

  # observed information: the la-diagonal carries the curvature d2 eta / d la^2
  D <- design(cur); u <- y - cur$p; av <- cur$p * (1 - cur$p)
  H <- crossprod(D, D * av)
  for (j in seq_len(nf)) {
    s <- free[j]; aj <- cur$alpha[[s]]
    onA <- !is.na(fa) & fa == j; onB <- !is.na(fb) & fb == j
    curv <- phg * (onA * aj * ba - onB * aj * bb)
    H[j, j] <- H[j, j] - sum(u * curv)
  }
  cov <- tryCatch(solve(H), error = function(e)
    solve(crossprod(D, D * av) + diag(1e-8, np)))
  conv <- max(abs(crossprod(D, u))) < 1e-4
  se <- sqrt(pmax(diag(cov), 0))
  list(alpha = cur$alpha, kappa = cur$kappa, p = cur$p, ll = cur$ll,
       cov = cov, se_log_alpha = setNames(se[seq_len(nf)], free),
       se_kappa = setNames(se[nf + seq_len(nf)], free),
       free = free, converged = conv)
}

# pooled log-of-mean-square fit residual over a set of comparisons, using the
# frame model's fitted probabilities (the paired-comparison residual logic of
# btl(): z = (y - p) / sqrt(p(1-p)), Andrich and Marais 2019 ch. 23)
.btlef_frame_fit <- function(y, p) {
  n <- length(y)
  if (n < 3L) return(NA_real_)
  V <- pmax(p * (1 - p), 1e-12)
  z2 <- (y - p)^2 / V
  mu4 <- p * (1 - p)^4 + (1 - p) * p^4
  c4v <- mu4 / V^2 - 1
  y2 <- sum(z2); f <- n; v <- sum(c4v)
  if (v > 1e-8 && y2 > 0) f * (log(y2) - log(f)) / sqrt(v) else NA_real_
}

#' Fit the extended frame of reference model for paired comparisons
#'
#' Estimates object locations from paired comparisons when the unit of the
#' latent scale differs across frames -- judge-panel by object-set cells -- an
#' extension of Humphry's (2005) extended frame of reference model to the
#' Bradley-Terry-Luce family. Objects are partitioned into sets and judges into
#' panels. Each object has a common-scale value \code{v_k = alpha_s beta_k +
#' kappa_s}, with \code{beta_k} the within-set calibration location,
#' \code{alpha_s > 0} the set unit and \code{kappa_s} the set origin; a
#' comparison judged in panel \code{g} carries the panel unit \code{phi_g}.
#' Within a set the comparison logit is \code{phi_g (beta_a - beta_b)} (the set
#' origin cancels and the set unit is confounded with the spread of
#' \code{beta}, so neither is identified within a set); across sets it is
#' \code{phi_g (v_a - v_b)}, which places the sets on one common scale.
#'
#' Estimation follows \code{\link{rasch_efrm}} in two conditional stages. In
#' stage one the within-set comparisons of each set, pooled over panels, fit
#' the bilinear model \code{logit = rho_{gs} (b_a - b_b)} with \code{b}
#' sum-zero and the most-used panel fixed at \code{rho = 1}; the ratios
#' \code{rho_{gs} = phi_g / phi_{ref(s)}} estimate the panel units up to each
#' set's reference panel and are reconciled across sets by a precision-weighted
#' least squares over the panel-by-set linking graph, then normalised to
#' geometric mean one. In stage two the cross-set comparisons are a
#' low-dimensional maximum likelihood in \code{(log alpha, kappa)} for the
#' non-reference sets, with \code{alpha_1 = 1} and \code{kappa_1 = 0} fixing
#' the reference set (the first, alphabetically).
#'
#' The set units are identified here WITHIN the conditional framework: the
#' cross-set linking uses only comparison outcomes and makes no distributional
#' assumption about the objects. This is the substantive difference from the
#' persons-by-items EFRM, whose item-set units can only be identified from the
#' person side -- that is, from the distribution of the persons over the linked
#' sets. The paired-comparison design supplies its own conditional link, so no
#' such distributional step is needed. See Humphry (2005) and Humphry and
#' Andrich (2008) for the theory of the unit, and Thurstone (1927) and David
#' (1988) for the varying-discriminal-dispersion lineage from which the
#' frame-dependent unit descends. The paired-comparison form is this package's
#' extension of Humphry's model.
#'
#' Standard errors are staged and their honesty is explicit. The frame
#' locations and log panel ratios carry judge-clustered Godambe sandwich
#' standard errors from stage one; the panel units \code{phi} carry standard
#' errors from the reconciliation. The set units \code{log alpha} and origins
#' \code{kappa} carry inverse-observed-information standard errors from stage
#' two that are CONDITIONAL on stage one -- the stage-one uncertainty in
#' \code{beta} and \code{phi} is not propagated into them, exactly as the
#' hybrid standard errors of \code{\link{rasch_efrm}} do not propagate the
#' pairwise-stage uncertainty into the linking stage. The common-scale values
#' \code{v} and their standard errors combine the conditional pieces by the
#' delta method.
#'
#' A single set (\code{S = 1}) reduces the model to panel units alone; stage
#' two is skipped and the print states the panel-units model. When
#' additionally \code{G = 1} the fit reduces exactly to \code{\link{btl}} on
#' the same data. The equal-unit (single-unit) comparison refits plain
#' \code{\link{btl}} on all comparisons pooled and reports the descriptive
#' composite log-likelihood difference against the frames model; because that
#' comparison is a composite likelihood, the inference on the units is carried
#' by the Wald tests on \code{log phi_g} and \code{log alpha_s} in
#' \code{phi_table} and \code{alpha_table}.
#'
#' @param data A data frame with one comparison per row.
#' @param object_a,object_b Names of the columns holding the two compared
#'   objects.
#' @param winner Name of the column holding the winner; its value must equal
#'   the row's \code{object_a} or \code{object_b} entry, \code{"tie"} or
#'   \code{"draw"} marks a tie, and anything else is treated as missing.
#' @param judge Name of the judge column (clusters the stage-one standard
#'   errors and defines the panels when \code{panels} is a judge attribute).
#' @param panels Either the name of a judge-attribute column in \code{data} or
#'   a named vector mapping judge to panel.
#' @param object_sets A named list mapping set names to character vectors of
#'   object names; every compared object must belong to exactly one set.
#' @param response Not supported: this first implementation fits dichotomous
#'   winner data only. Supplying it raises an informative error.
#' @param ties \code{"drop"} (default, removed with a note) or \code{"error"}.
#' @param min_link Minimum number of cross-set comparisons a set pair must
#'   supply to be used for linking; sets not reachable from the reference set
#'   through sufficient cross-set pairs raise an error.
#' @param maxit,tol Newton iteration cap and convergence tolerance.
#' @return An object of class \code{"rasch_btl_efrm"}: \code{objects} (object,
#'   set, \code{beta_set} and its standard error, common-scale \code{v} and its
#'   standard error), \code{phi_table} (panel units with Wald tests against
#'   \code{log phi = 0}), \code{alpha_table} and \code{kappa_table} (set units
#'   and origins with Wald tests; the reference set carries \code{alpha = 1},
#'   \code{kappa = 0} with no standard error), \code{frames} (panel by set:
#'   unit \code{rho = phi alpha}, comparison count, pooled fit residual),
#'   \code{equal_unit} (the descriptive single-unit comparison), \code{n_cross}
#'   (cross-set comparison counts per set pair), \code{notes} and
#'   \code{converged}.
#' @references Bradley, R. A. and Terry, M. E. (1952). Rank analysis of
#'   incomplete block designs: I. The method of paired comparisons.
#'   Biometrika, 39, 324-345.
#'
#'   David, H. A. (1988). The Method of Paired Comparisons (2nd ed.). Griffin.
#'
#'   Humphry, S. M. (2005). Maintaining a common arbitrary unit in social
#'   measurement. PhD thesis, Murdoch University.
#'
#'   Humphry, S. M. and Andrich, D. (2008). Understanding the unit in the
#'   Rasch model. Journal of Applied Measurement, 9(3), 249-264.
#'
#'   Luce, R. D. (1959). Individual Choice Behavior. Wiley.
#'
#'   Thurstone, L. L. (1927). A law of comparative judgment. Psychological
#'   Review, 34, 273-286.
#' @examples
#' \donttest{
#' d <- simulate_btl_efrm(n_objects_per_set = 6, n_sets = 2, n_panels = 2,
#'                        set_units = c(1, 1.4), set_origins = c(0, 0.8),
#'                        seed = 1)
#' fit <- btl_efrm(d, "object_a", "object_b", winner = "winner",
#'                 judge = "judge", panels = "panel",
#'                 object_sets = attr(d, "truth")$object_sets)
#' fit$alpha_table
#' }
#' @export
btl_efrm <- function(data, object_a, object_b, winner, judge, panels,
                     object_sets, response = NULL,
                     ties = c("drop", "error"), min_link = 20,
                     maxit = 60, tol = 1e-8) {
  ties <- match.arg(ties)
  if (!is.null(response))
    stop("btl_efrm fits dichotomous winner data only in this first ",
         "implementation; a graded `response` is not supported. Reduce the ",
         "graded margins to a winner, or use btl() for a single-frame graded ",
         "analysis.")
  data <- as.data.frame(data)
  for (col in c(object_a, object_b, winner, judge))
    if (!col %in% names(data)) stop("column not found: ", col)
  a <- trimws(as.character(data[[object_a]]))
  b <- trimws(as.character(data[[object_b]]))
  wn <- trimws(as.character(data[[winner]]))
  jd <- as.character(data[[judge]])

  # panels: a judge-attribute column, or a named judge -> panel vector
  if (length(panels) == 1L && is.character(panels) && panels %in% names(data)) {
    pan <- as.character(data[[panels]])
  } else if (!is.null(names(panels)) && all(nzchar(names(panels)))) {
    pan <- unname(as.character(panels)[match(jd, names(panels))])
  } else {
    stop("`panels` must name a column of `data` or be a named vector ",
         "mapping judge to panel")
  }
  notes <- character(0)

  keep <- !is.na(a) & !is.na(b) & !is.na(wn) & !is.na(jd) & !is.na(pan) & a != b
  if (any(!keep)) {
    notes <- c(notes, sprintf("%d row(s) dropped (missing or self-comparison)",
                              sum(!keep)))
    a <- a[keep]; b <- b[keep]; wn <- wn[keep]; jd <- jd[keep]; pan <- pan[keep]
  }
  if (!length(a)) stop("no usable comparisons")

  y <- ifelse(wn == a, 1L, ifelse(wn == b, 0L, NA_integer_))
  is_tie <- is.na(y) & tolower(wn) %in% c("tie", "draw")
  miss <- is.na(y) & !is_tie
  if (any(miss)) {
    notes <- c(notes, sprintf(
      "%d row(s) with winner matching neither object treated as missing", sum(miss)))
    sel <- !miss
    a <- a[sel]; b <- b[sel]; y <- y[sel]; jd <- jd[sel]; pan <- pan[sel]
  }
  if (anyNA(y)) {
    nt <- sum(is.na(y))
    if (ties == "error") stop(nt, " tie(s) present; set ties = 'drop'")
    notes <- c(notes, sprintf("%d tie(s) dropped", nt))
    sel <- !is.na(y)
    a <- a[sel]; b <- b[sel]; y <- y[sel]; jd <- jd[sel]; pan <- pan[sel]
  }
  if (!length(a)) stop("no usable comparisons after cleaning")

  # --- object sets ----------------------------------------------------------
  if (!is.list(object_sets) || is.null(names(object_sets)) ||
      any(!nzchar(names(object_sets))))
    stop("`object_sets` must be a named list: set name -> object names")
  objs_all <- sort(unique(c(a, b)))
  set_of <- setNames(rep(NA_character_, length(objs_all)), objs_all)
  multi <- character(0)
  for (s in names(object_sets)) {
    hit <- intersect(as.character(object_sets[[s]]), objs_all)
    for (o in hit) {
      if (!is.na(set_of[o])) multi <- c(multi, o)
      set_of[o] <- s
    }
  }
  if (length(multi))
    stop("object(s) assigned to more than one set: ",
         paste(unique(multi), collapse = ", "))
  if (anyNA(set_of))
    stop("object(s) in the data not found in `object_sets` (every compared ",
         "object must belong to exactly one set): ",
         paste(objs_all[is.na(set_of)], collapse = ", "))
  sets_u <- sort(unique(set_of)); S <- length(sets_u)
  panels_u <- sort(unique(pan)); G <- length(panels_u)
  sa <- set_of[a]; sb <- set_of[b]
  within <- sa == sb

  if (any(table(set_of) < 2L))
    stop("every set needs at least two objects; offending set(s): ",
         paste(names(which(table(set_of) < 2L)), collapse = ", "))

  # --- stage 1: within-set bilinear solves ----------------------------------
  bhat <- setNames(rep(NA_real_, length(objs_all)), objs_all)
  se_bhat <- bhat
  phi_ref_of_set <- setNames(rep(NA_real_, S), sets_u)
  ref_of_set <- setNames(rep(NA_character_, S), sets_u)
  within_p <- rep(NA_real_, length(a))                # frame-model fitted p
  blocks <- list()                                    # per-set panel-ratio blocks
  s1 <- vector("list", S); names(s1) <- sets_u
  ll_within <- 0
  for (s in sets_u) {
    rows <- which(within & sa == s)
    os <- sort(names(set_of)[set_of == s]); Ks <- length(os)
    ia <- match(a[rows], os); ib <- match(b[rows], os)
    if (anyNA(ia) || anyNA(ib) || any(is.na(match(os, unique(c(a[rows], b[rows]))))))
      stop("set '", s, "' has object(s) with no within-set comparison; ",
           "each object needs at least one comparison inside its own set")
    fit1 <- .btlef_stage1(ia, ib, y[rows], pan[rows], jd[rows], Ks, maxit, tol)
    s1[[s]] <- fit1
    bhat[os] <- fit1$beta; se_bhat[os] <- fit1$se_beta
    within_p[rows] <- fit1$p
    ref_of_set[s] <- fit1$ref
    ll_within <- ll_within + fit1$ll
    blocks[[s]] <- list(ref = fit1$ref, free = fit1$free,
                        lrho = setNames(log(fit1$rho[fit1$free]), fit1$free),
                        cov = fit1$cov_lrho)
  }

  # --- reconcile panel units and put beta on the common panel scale ---------
  rec <- .btlef_reconcile_phi(panels_u, blocks)
  phi <- rec$phi
  for (s in sets_u) {
    pr <- phi[[ref_of_set[s]]]; phi_ref_of_set[s] <- pr
    os <- names(set_of)[set_of == s]
    bhat[os] <- bhat[os] / pr; se_bhat[os] <- se_bhat[os] / pr
  }
  # recompute within-set fitted p on the common scale: logit = phi_g (bhat_a - bhat_b)
  within_p[within] <- plogis(phi[pan[within]] * (bhat[a[within]] - bhat[b[within]]))

  stage1_conv <- all(vapply(s1, function(z) z$converged, TRUE))

  # --- stage 2: cross-set linking -------------------------------------------
  cross <- which(!within)
  n_cross <- data.frame(set_a = character(0), set_b = character(0),
                        n = integer(0), stringsAsFactors = FALSE)
  alpha <- setNames(rep(1, S), sets_u); kappa <- setNames(rep(0, S), sets_u)
  se_log_alpha <- setNames(rep(NA_real_, S), sets_u)
  se_kappa <- setNames(rep(NA_real_, S), sets_u)
  cov2 <- NULL; stage2_conv <- TRUE; ll_cross <- 0; cross_p <- NULL
  if (S > 1L) {
    if (!length(cross))
      stop("no cross-set comparisons: the sets cannot be linked to a common ",
           "scale (set units alpha and origins kappa are unidentified)")
    # unordered set-pair counts
    key <- ifelse(sa[cross] < sb[cross], paste(sa[cross], sb[cross]),
                  paste(sb[cross], sa[cross]))
    tab <- table(key); parts <- do.call(rbind, strsplit(names(tab), " ", fixed = TRUE))
    n_cross <- data.frame(set_a = parts[, 1], set_b = parts[, 2],
                          n = as.integer(tab), stringsAsFactors = FALSE)
    rownames(n_cross) <- NULL
    used <- n_cross$n >= min_link
    edges <- cbind(match(n_cross$set_a[used], sets_u),
                   match(n_cross$set_b[used], sets_u))
    comp <- .btlef_components(S, edges)
    ref_comp <- comp[1]                                 # reference set = sets_u[1]
    if (any(comp != ref_comp))
      stop("set(s) not reachable from the reference set '", sets_u[1],
           "' through cross-set pairs with at least min_link = ", min_link,
           " comparisons: ", paste(sets_u[comp != ref_comp], collapse = ", "),
           " (increase the cross-set data or lower min_link)")
    st2 <- .btlef_stage2(a[cross], b[cross], y[cross], phi[pan[cross]],
                         sa[cross], sb[cross], bhat, sets_u, maxit, tol)
    alpha <- st2$alpha; kappa <- st2$kappa; cov2 <- st2$cov
    se_log_alpha[st2$free] <- st2$se_log_alpha
    se_kappa[st2$free] <- st2$se_kappa
    stage2_conv <- st2$converged; ll_cross <- st2$ll; cross_p <- st2$p
  }

  # --- common-scale values with delta-method standard errors ----------------
  v <- alpha[set_of[objs_all]] * bhat[objs_all] + kappa[set_of[objs_all]]
  se_v <- se_bhat[objs_all]                             # reference set: v = beta
  free <- sets_u[-1L]
  if (S > 1L) for (o in objs_all) {
    s <- set_of[[o]]; if (s == sets_u[1]) next
    j <- match(s, free); idx <- c(j, (S - 1L) + j)
    C2 <- cov2[idx, idx, drop = FALSE]
    gvec <- c(alpha[[s]] * bhat[[o]], 1)                # d v / d(log alpha, kappa)
    var_link <- drop(t(gvec) %*% C2 %*% gvec)
    se_v[[o]] <- sqrt(pmax(alpha[[s]]^2 * se_bhat[[o]]^2 + var_link, 0))
  }

  # --- equal-unit (single-unit) comparison ----------------------------------
  ll_frames <- ll_within + ll_cross
  single <- tryCatch(
    .btlef_stage1(match(a, objs_all), match(b, objs_all), y,
                  rep("all", length(a)), jd, length(objs_all), maxit, tol),
    error = function(e) NULL)
  ll_single <- if (is.null(single)) NA_real_ else single$ll
  equal_unit <- list(
    loglik_frames = ll_frames, loglik_single = ll_single,
    difference = if (is.na(ll_single)) NA_real_ else ll_frames - ll_single,
    note = paste("descriptive composite-likelihood difference;",
                 "the Wald tests on log phi and log alpha carry the inference"))

  # --- structural tables ----------------------------------------------------
  z_phi <- log(phi) / rec$se_log_phi
  phi_table <- data.frame(panel = panels_u, phi = unname(phi),
                          se_log_phi = unname(rec$se_log_phi),
                          z = unname(z_phi), p = unname(2 * pnorm(-abs(z_phi))),
                          stringsAsFactors = FALSE)
  z_al <- log(alpha) / se_log_alpha
  alpha_table <- data.frame(set = sets_u, alpha = unname(alpha),
                            se_log_alpha = unname(se_log_alpha),
                            z = unname(z_al), p = unname(2 * pnorm(-abs(z_al))),
                            stringsAsFactors = FALSE)
  z_ka <- kappa / se_kappa
  kappa_table <- data.frame(set = sets_u, kappa = unname(kappa),
                            se_kappa = unname(se_kappa),
                            z = unname(z_ka), p = unname(2 * pnorm(-abs(z_ka))),
                            stringsAsFactors = FALSE)

  objects <- data.frame(object = objs_all, set = unname(set_of[objs_all]),
                        beta_set = unname(bhat[objs_all]),
                        se_beta = unname(se_bhat[objs_all]),
                        v = unname(v), se_v = unname(se_v),
                        stringsAsFactors = FALSE)
  rownames(objects) <- NULL

  # frames: one row per panel-by-set cell holding within-set comparisons
  fr <- list()
  for (s in sets_u) for (g in panels_u) {
    rows <- which(within & sa == s & pan == g)
    if (!length(rows)) next
    fr[[length(fr) + 1L]] <- data.frame(
      panel = g, set = s, rho = unname(phi[[g]] * alpha[[s]]),
      n_comparisons = length(rows),
      fit_resid = .btlef_frame_fit(y[rows], within_p[rows]),
      stringsAsFactors = FALSE)
  }
  frames <- if (length(fr)) do.call(rbind, fr) else NULL
  if (!is.null(frames)) rownames(frames) <- NULL

  if (S == 1L)
    notes <- c(notes, "single set: panel-units model (set units alpha not estimated)")
  if (G == 1L && S == 1L)
    notes <- c(notes, "single panel and single set: reduces to btl()")

  out <- list(objects = objects, phi_table = phi_table,
              alpha_table = alpha_table, kappa_table = kappa_table,
              frames = frames, equal_unit = equal_unit, n_cross = n_cross,
              sets = sets_u, panels = panels_u, reference_set = sets_u[1],
              n_comparisons = length(a),
              converged = stage1_conv && stage2_conv,
              se_note = paste("stage-2 alpha and kappa standard errors are",
                              "conditional on stage 1; stage-1 uncertainty in",
                              "beta and phi is not propagated into them"),
              notes = notes)
  class(out) <- "rasch_btl_efrm"
  out
}

#' @export
print.rasch_btl_efrm <- function(x, ...) {
  cat(sprintf(paste0("Bradley-Terry-Luce extended frame of reference: ",
                     "%d objects in %d set(s) x %d panel(s), %d comparisons\n"),
              nrow(x$objects), nrow(x$alpha_table), nrow(x$phi_table),
              x$n_comparisons))
  cat(sprintf("Two-stage conditional ML: %s\n",
              if (x$converged) "converged" else "NOT converged"))
  if (nrow(x$alpha_table) == 1L)
    cat("Model: panel units only (single set; set units not estimated)\n")
  cat("\nPanel units (phi; Wald H0: log phi = 0):\n")
  print(.fmt_df(x$phi_table), row.names = FALSE)
  if (nrow(x$alpha_table) > 1L) {
    cat("\nSet units (alpha) and origins (kappa; reference set = ",
        x$reference_set, "):\n", sep = "")
    at <- merge(x$alpha_table[, c("set", "alpha", "se_log_alpha", "p")],
                x$kappa_table[, c("set", "kappa", "se_kappa")],
                by = "set", sort = FALSE)
    print(.fmt_df(at), row.names = FALSE)
  }
  eu <- x$equal_unit
  if (!is.na(eu$difference))
    cat(sprintf(paste0("\nEqual-unit comparison: ll_frames - ll_single = ",
                       "%.3f (%s)\n"), eu$difference, eu$note))
  if (length(x$notes)) cat(sprintf("Notes: %s\n", paste(x$notes, collapse = "; ")))
  invisible(x)
}

#' Plot the frame units of a paired-comparison EFRM fit
#'
#' Caterpillar plot of the estimated units on the log scale: one row per panel
#' unit \code{phi_g} and one per set unit \code{alpha_s}, with 95 per cent
#' intervals, the reference (unit one) marked, mirroring
#' \code{\link{plot_frames}} in the package's house style.
#'
#' @param fit A fitted object from \code{\link{btl_efrm}}.
#' @return Called for its plotting side effect; invisibly \code{NULL}.
#' @examples
#' \donttest{
#' # see ?btl_efrm for a complete simulated example
#' }
#' @export
plot_btl_units <- function(fit) {
  if (!inherits(fit, "rasch_btl_efrm"))
    stop("plot_btl_units needs a rasch_btl_efrm fit")
  ph <- fit$phi_table; al <- fit$alpha_table
  rows <- rbind(
    data.frame(label = paste0("panel: ", ph$panel), kind = "panel",
               est = log(ph$phi), se = ph$se_log_phi, stringsAsFactors = FALSE),
    if (nrow(al) > 1L)
      data.frame(label = paste0("set: ", al$set), kind = "set",
                 est = log(al$alpha), se = al$se_log_alpha,
                 stringsAsFactors = FALSE))
  rows$se[!is.finite(rows$se)] <- 0
  rows <- rows[order(rows$kind, rows$est), ]
  n <- nrow(rows)
  lo <- rows$est - 1.96 * rows$se; hi <- rows$est + 1.96 * rows$se
  colr <- ifelse(rows$kind == "panel", .rr$blue, .rr$purple)
  op <- par(mar = c(4.2, 9, 3.2, 1.5), mgp = c(2.5, 0.7, 0), tcl = -0.25,
            las = 1, col.axis = .rr$ink, col.lab = .rr$ink, col.main = .rr$ink,
            font.main = 2, cex.main = 1.15)
  on.exit(par(op))
  plot(NA, xlim = range(c(lo, hi, 0)) + c(-0.1, 0.1), ylim = c(0.5, n + 0.5),
       xlab = "log unit", ylab = "", axes = FALSE, main = "")
  abline(h = seq_len(n), col = .rr$grid, lwd = 0.8)
  abline(v = 0, lty = 2, col = .rr$soft)
  axis(1, col = .rr$grid, col.ticks = .rr$soft)
  axis(2, at = seq_len(n), labels = rows$label, cex.axis = 0.75,
       col = .rr$grid, col.ticks = NA)
  segments(lo, seq_len(n), hi, seq_len(n), lwd = 2.2, col = .rr$soft)
  points(rows$est, seq_len(n), pch = 21, cex = 1.5, bg = colr,
         col = "white", lwd = 1.2)
  .rr_legend("bottomright", c("panel unit (phi)", "set unit (alpha)"),
             pch = 21, pt.bg = c(.rr$blue, .rr$purple), col = "white",
             pt.cex = 1.2)
  invisible(NULL)
}
