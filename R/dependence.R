# rasch :: quantifying response dependence
# ===========================================================================
# Two magnitude estimators for local (response) dependence beyond the
# residual-correlation screen. (1) The magnitude estimate d of Andrich and
# Kreiner (2010), generalised to polytomous items by Andrich, Humphry and
# Marais (2012): the dependent item is resolved into one item per category
# of the independent item, both originals are deleted, and dependence shows
# as opposite shifts +/- d of the resolved thresholds. (2) The spread-
# parameter screen of Andrich (1985): a subtest formed from independent
# items cannot have a spread (half threshold distance) below the value the
# binomial distribution gives for its maximum score, so a spread estimate
# under that least upper bound indicates dependence among the members.
# ===========================================================================

#' Estimate the magnitude of response dependence between two items
#'
#' Quantifies how strongly a dependent item's response follows an
#' independent item's response, in logits, by the resolution method of
#' Andrich and Kreiner (2010; polytomous generalisation Andrich, Humphry
#' and Marais 2012), by resolution of the dependent item. The
#' dependent item is resolved into one item per category of the independent
#' item (each carrying the responses of the persons who gave that category),
#' both original items are removed, and the model refitted. Under
#' dependence of magnitude \eqn{d}, threshold \eqn{k} of the resolved item
#' for category \eqn{x_i} is shifted by \eqn{-d} when \eqn{k \le x_i} and
#' \eqn{+d} otherwise, so each threshold yields
#' \eqn{\hat d_k = (\hat\delta_{ji(k)}(x_i = k-1) - \hat\delta_{ji(k)}(x_i = k))/2}
#' and \eqn{\hat d} is their mean (eq. 24.7 of Andrich and Marais 2019).
#' The resolved threshold estimates share the calibration of the remaining
#' items and are therefore correlated. The standard error of \eqn{\hat d} is
#' calculated from their full sandwich covariance. If that covariance is
#' unavailable or not positive semidefinite, the point estimate is retained
#' as a descriptive magnitude and Wald inference is withheld.
#' When repeated person identifiers contribute to the refit, probabilities use
#' a t reference with the number of independent person clusters minus one
#' degree of freedom. Independent response rows retain the asymptotic normal
#' reference, represented by infinite degrees of freedom.
#'
#' Polytomous resolution requires an unconstrained partial credit model so
#' that each resolved threshold can move independently. A rating scale or
#' principal-component threshold constraint is therefore refused. The refit
#' otherwise retains the original fit grouping, keyed scoring, anchors on
#' items that remain, and optimisation controls. Both calibrations must
#' converge before the magnitude and its standard error are reported. MFRM
#' virtual items are resolved through this unconstrained PCM. EFRM virtual
#' frames are mutually exclusive and must first be reduced to an observable
#' frame or linked design block. For an explanatory fit, the remaining items
#' retain their explanatory restrictions and the resolved copies receive free
#' fixed departures.
#' Resolution is refused if any retained item or resolved copy loses its
#' original score categories during calibration.
#'
#' @param fit A fitted object from \code{\link{rasch}}.
#' @param dependent,independent Item names or indices: the item hypothesised
#'   to depend, and the item it depends on. Both must share the same maximum
#'   score (the formalisation requires it).
#' @return A list of class \code{"rasch_dependence"}: the estimate \code{d},
#'   its \code{se}, \code{t}, reference \code{df} and \code{p} for the
#'   hypothesis \eqn{d = 0},
#'   the per-threshold table \code{thresholds} (columns \code{k},
#'   \code{delta_lo}, \code{delta_hi}, \code{d_k}, \code{se_k}), and the
#'   resolved \code{refit}.
#' @references Andrich, D. and Kreiner, S. (2010). Quantifying response
#'   dependence between two dichotomous items using the Rasch model.
#'   Applied Psychological Measurement, 34, 181-192. Andrich, D., Humphry,
#'   S. M. and Marais, I. (2012). Quantifying local, response dependence
#'   between two polytomous items using the Rasch model. Applied
#'   Psychological Measurement, 36, 309-324.
#' @examples
#' set.seed(1); N <- 700
#' d0 <- seq(-1.5, 1.5, length.out = 8)
#' X <- matrix(rbinom(N * 8, 1, plogis(outer(rnorm(N), d0, "-"))), N, 8)
#' X[, 5] <- ifelse(runif(N) < 0.75, X[, 4], X[, 5])   # I5 follows I4
#' colnames(X) <- paste0("I", 1:8)
#' dependence_magnitude(rasch(X), dependent = "I5", independent = "I4")
#' @export
dependence_magnitude <- function(fit, dependent, independent) {
  if (length(dependent) != 1L || length(independent) != 1L)
    stop("`dependent` and `independent` must each name exactly one item")
  if (!inherits(fit, "rasch")) stop("dependence_magnitude needs a rasch fit")
  if (inherits(fit, "rasch_efrm"))
    .refuse("dependence magnitude is not defined across mutually exclusive ",
            "EFRM virtual frames; extract and calibrate an observable frame ",
            "or linked design block first")
  if (!isTRUE(fit$est$converged))
    stop("the fitted calibration did not converge; dependence magnitude is unavailable")
  j <- .item_idx(fit, dependent); i <- .item_idx(fit, independent)
  if (i == j) stop("dependent and independent must be different items")
  mi <- fit$m[i]; mj <- fit$m[j]
  if (mi != mj)
    stop("the two items must share the same maximum score (here ",
         mj, " and ", mi, "); see Andrich, Humphry and Marais (2012)")
  if (mi > 1L && (identical(fit$model, "RSM") ||
                  !is.null(fit$refit_spec$pc_components)))
    stop("polytomous dependence magnitude needs an unconstrained PCM ",
         "resolution; refit with model = \"PCM\" and no principal-component ",
         "threshold constraint")
  X <- fit$X
  nm_i <- colnames(X)[i]; nm_j <- colnames(X)[j]

  # resolve item j by the categories of item i; drop both originals
  keep <- setdiff(colnames(X), c(nm_i, nm_j))
  Xn <- X[, keep, drop = FALSE]
  res_names <- character(mi + 1)
  for (x in 0:mi) {
    col <- X[, j]
    col[is.na(X[, i]) | X[, i] != x] <- NA
    obs <- sort(unique(col[!is.na(col)]))
    if (!identical(obs, 0:mj))
      stop(sprintf(paste0("resolved item for %s = %d does not observe every ",
                          "category of %s (found: %s); too little data to ",
                          "resolve"), nm_i, x, nm_j,
                   paste(obs, collapse = ",")))
    Xn <- cbind(Xn, col)
    res_names[x + 1] <- sprintf("%s|%s=%d", nm_j, nm_i, x)
    if (res_names[x + 1] %in% colnames(Xn)[-ncol(Xn)])
      stop("generated resolved-item name already exists: ", res_names[x + 1])
    colnames(Xn)[ncol(Xn)] <- res_names[x + 1]
  }
  # The resolution method requires every resolved threshold to be free. This
  # also routes an MFRM virtual-item analysis through the ordinary PCM rather
  # than handing the structural model label to rasch().
  score_max <- c(stats::setNames(fit$m[match(keep, fit$items$item)], keep),
                 stats::setNames(rep(mj, length(res_names)), res_names))
  refit <- if (inherits(fit, "rasch_explanatory")) {
    inherit <- c(stats::setNames(keep, keep),
                 stats::setNames(rep(nm_j, length(res_names)), res_names))
    .explanatory_refit_modified(fit, Xn, inherit = inherit,
                                fully_relaxed = res_names)
  } else .rasch_refit(fit, Xn, model = "PCM", require_anchor = FALSE,
                      score_max = score_max)
  .require_fitted_score_structure(refit, score_max, "the dependence refit")
  if (!isTRUE(refit$est$converged))
    stop("the resolved calibration did not converge; dependence magnitude is unavailable")
  if (!all(res_names %in% refit$items$item))
    stop("a resolved item was dropped during re-analysis; too little data")

  thr <- refit$thresholds
  item_of <- refit$items$item[thr$item]
  cv <- refit$est$cov_tau
  cov_ok <- length(dim(cv)) == 2L && nrow(cv) == ncol(cv) &&
    nrow(cv) == nrow(thr) && all(is.finite(cv)) &&
    .covariance_is_symmetric(cv) && .covariance_is_psd(cv)
  w <- numeric(nrow(thr))
  tab <- data.frame(k = seq_len(mj), delta_lo = NA_real_, delta_hi = NA_real_,
                    d_k = NA_real_, se_k = NA_real_)
  for (k in seq_len(mj)) {
    lo <- thr[item_of == res_names[k]     & thr$k == k, ]  # x_i = k - 1: +d
    hi <- thr[item_of == res_names[k + 1] & thr$k == k, ]  # x_i = k    : -d
    tab$delta_lo[k] <- lo$tau; tab$delta_hi[k] <- hi$tau
    tab$d_k[k] <- (lo$tau - hi$tau) / 2
    if (cov_ok)
      tab$se_k[k] <- sqrt(pmax((cv[lo$id, lo$id] + cv[hi$id, hi$id] -
                                2 * cv[lo$id, hi$id]) / 4, 0))
    w[lo$id] <- w[lo$id] + 1; w[hi$id] <- w[hi$id] - 1
  }
  w <- w / (2 * mj)
  d <- mean(tab$d_k)
  se <- if (cov_ok) sqrt(pmax(drop(crossprod(w, cv %*% w)), 0)) else NA_real_
  statistic <- .wald_ratio(d, se)
  df <- .dif_refit_df(refit)
  # a magnitude built on weakly identified resolved thresholds keeps its
  # descriptive value, but the covariance behind its standard error is not
  # trustworthy; inference is withheld, as it is for the thresholds
  # themselves
  res_thr <- thr[item_of %in% res_names, ]
  weak_res <- isTRUE(any(res_thr$weak))
  note <- character(0)
  if (weak_res) {
    se <- statistic <- df <- NA_real_
    note <- c(note, paste("resolved thresholds are weakly identified (sparse",
                         "resolved categories); the magnitude is descriptive and",
                         "inference is withheld"))
    tab$se_k <- NA_real_
  }
  if (!cov_ok) {
    se <- statistic <- df <- NA_real_
    tab$se_k <- NA_real_
    note <- c(note, paste("the resolved-threshold covariance is unavailable or",
                         "not positive semidefinite; the magnitude is descriptive",
                         "and inference is withheld"))
  } else if (!weak_res && !is.finite(statistic)) {
    se <- statistic <- df <- NA_real_
    note <- c(note, paste("the resolved contrast has zero estimated uncertainty;",
                         "Wald inference is withheld"))
  } else if (!weak_res && is.na(df)) {
    se <- statistic <- NA_real_
    tab$se_k <- NA_real_
    note <- c(note, paste(
      "the resolved calibration does not carry enough independent",
      "person-cluster support for a reference distribution; inference is withheld"))
  }
  out <- list(d = d, se = se, t = statistic, df = df,
              p = if (!is.finite(statistic) || is.na(df)) NA_real_ else
                2 * stats::pt(-abs(statistic), df = df),
              thresholds = tab, dependent = nm_j, independent = nm_i,
              note = if (length(note)) paste(note, collapse = "; ") else NULL,
              refit = refit)
  out <- .tag_tables(out)
  class(out) <- "rasch_dependence"
  out
}

#' @export
print.rasch_dependence <- function(x, ...) {
  cat(sprintf("Response dependence of %s on %s (Andrich & Kreiner resolution)\n",
              x$dependent, x$independent))
  # Exact access matters for pre-`t` saved fits: they have `thresholds` but
  # no `t`, and `$t` would partially match that data frame.
  statistic <- x[["t"]] %||% x[["z"]]
  if (is.na(x$se) || is.na(statistic) || is.na(x$p))
    cat(sprintf("  d = %.3f logits (descriptive; inference withheld)\n", x$d))
  else {
    df <- x[["df"]] %||% Inf
    reference <- if (is.finite(df)) sprintf("t(%g)", df) else "z"
    cat(sprintf("  d = %.3f logits (se %.3f), %s = %.2f, p = %s\n",
                x$d, x$se, reference, statistic, .fmt_p(x$p)))
  }
  if (!is.null(x$note)) cat("  Note:", x$note, "\n")
  if (nrow(x$thresholds) > 1) {
    cat("  per threshold:\n")
    print(.fmt_df(x$thresholds), row.names = FALSE)
  }
  invisible(x)
}

# Andrich (1985) least upper bounds for the spread parameter of an
# independent subtest, from the binomial threshold structure (Andrich &
# Marais 2019, Table 24.1).
.spread_lub <- c(`2` = 0.69, `3` = 0.55, `4` = 0.41, `5` = 0.35,
                 `6` = 0.29, `7` = 0.25, `8` = 0.22)

#' Spread-parameter test for dependence within subtests
#'
#' Andrich's (1985) least-upper-bound screen: the spread component
#' \eqn{\lambda} of a polytomous item (half the distance between successive
#' thresholds in the principal-components parameterisation, estimated here
#' by \code{\link{pcml_pc}}) cannot fall below the value implied by the
#' binomial distribution when the item is a subtest of equally difficult,
#' independent dichotomous items; different difficulties only raise it.
#' Sampling uncertainty matters when the fitted spread lies near the bound.
#' The function therefore reports whether the point estimate is below the
#' bound separately from a one-sided test of \eqn{H_0: \lambda \geq \lambda_0}
#' against \eqn{H_1: \lambda < \lambda_0}. Evidence of dependence requires the
#' adjusted one-sided probability to be below \code{alpha}; a point estimate
#' below the bound alone is not treated as a verdict. The tabulated bounds
#' are exact for subtests of two and three items. For larger subtests the
#' binomial thresholds are no longer equally spaced, and the spread that
#' \code{\link{pcml_pc}} recovers from independent, equally difficult
#' components sits above the tabulated value (about 0.45, 0.41, 0.34, and
#' 0.27 for four, five, six, and eight items against 0.41, 0.35, 0.29, and
#' 0.22), so the screen is conservative there: a subtest of four or more
#' items needs a spread well below the tabulated bound before the test
#' reports dependence.
#' Applied to the superitems recorded by \code{\link{combine_items}}. The
#' binomial bound applies only when every component was dichotomous; a
#' composite containing a polytomous item is shown but its bound and verdict
#' are withheld. When person identifiers repeat, the spread refit's sandwich
#' covariance clusters the score contributions by person and probabilities use
#' a t reference with the number of independent person clusters minus one
#' degree of freedom. Independent rows retain the asymptotic normal reference.
#' The input calibration and the principal-components refit must both converge.
#'
#' @param fit A fitted object from \code{\link{rasch}}.
#' @param maxit,tol Passed to the \code{\link{pcml_pc}} refit.
#' @param alpha Significance level for the one-sided dependence screen.
#' @param p_adjust Multiplicity adjustment across the eligible superitems;
#'   one of \code{stats::p.adjust.methods}. An eligible superitem remains in
#'   the family if its probability is unavailable.
#' @return A data frame with one row per recorded superitem: \code{item},
#'   \code{m}, whether the binomial bound is \code{eligible}, the
#'   \code{spread} estimate and its \code{se}, the bound \code{lub}
#'   (available for dichotomous-component subtests with maximum scores 2 to
#'   8), \code{t} = (spread - lub)/se, reference \code{df}, the one-sided
#'   \code{p} and adjusted \code{p_adj}, \code{below_bound} for the
#'   point-estimate comparison, and
#'   \code{dependent} for adjusted evidence at \code{alpha}.
#'   Items not formed by \code{combine_items()} are omitted.
#'   The result retains \code{alpha} and \code{p_adjust} as attributes.
#' @references
#' Andrich, D. (1985). An elaboration of Guttman scaling with Rasch models
#' for measurement. In N. B. Tuma (Ed.), Sociological Methodology 1985
#' (pp. 33--80). Jossey-Bass.
#'
#' Andrich, D. and Marais, I. (2019). A Course in Rasch Measurement Theory:
#' Measuring in the Educational, Social and Health Sciences. Springer.
#' @examples
#' set.seed(1); N <- 600
#' d0 <- seq(-1.5, 1.5, length.out = 8)
#' X <- matrix(rbinom(N * 8, 1, plogis(outer(rnorm(N), d0, "-"))), N, 8)
#' X[, 5] <- ifelse(runif(N) < .85, X[, 4], 1 - X[, 4])
#' X[, 6] <- ifelse(runif(N) < .75, X[, 4], 1 - X[, 4]) # a dependent triple
#' colnames(X) <- paste0("I", 1:8)
#' fit2 <- combine_items(rasch(X), list(c("I4", "I5", "I6"), c("I1", "I2", "I3")))
#' spread_test(fit2)
#' @export
spread_test <- function(fit, maxit = 60, tol = 1e-8,
                        alpha = 0.05, p_adjust = "holm") {
  if (!inherits(fit, "rasch")) stop("spread_test needs a rasch fit")
  .check_prob(alpha, "alpha")
  if (!is.character(p_adjust) || length(p_adjust) != 1L ||
      !is.null(dim(p_adjust)) || !is.null(oldClass(p_adjust)) ||
      !p_adjust %in% stats::p.adjust.methods)
    stop("p_adjust must name a method in stats::p.adjust.methods")
  if (!isTRUE(fit$est$converged))
    stop("the fitted calibration did not converge; the spread test is unavailable")
  sub_names <- intersect(names(fit$subtest_map), fit$items$item)
  if (!length(sub_names))
    stop("no recorded superitems: form dichotomous items into subtests with combine_items() first")
  idx_fit <- match(sub_names, fit$items$item)
  eligible <- vapply(sub_names, function(it)
    isTRUE(fit$subtest_binary[[it]]) && fit$m[match(it, fit$items$item)] %in%
      as.integer(names(.spread_lub)), logical(1))
  if (!any(eligible)) {
    out <- data.frame(item = sub_names, m = fit$m[idx_fit], eligible = FALSE,
                      spread = NA_real_, se = NA_real_, lub = NA_real_,
                      t = NA_real_, df = NA_real_, p = NA_real_, p_adj = NA_real_,
                      below_bound = NA, dependent = NA)
    attr(out, "alpha") <- alpha
    attr(out, "p_adjust") <- p_adjust
    class(out) <- c("rasch_spread", "data.frame")
    return(out)
  }
  # pcml_pc() is also a public stand-alone estimator and therefore has no
  # person-ID argument. Here the source fit supplies those IDs: pass them to
  # the same internal fit so repeated observations contribute one clustered
  # sandwich score per person rather than masquerading as independent rows.
  pc <- .pcml_pc_fit(fit$X, maxit = maxit, tol = tol,
                     cluster = .dif_ids(fit$person$id))
  if (!isTRUE(pc$converged))
    stop("the principal-components refit did not converge; the spread test is unavailable")
  cmp <- pc$components
  idx_pc <- match(sub_names, cmp$item)
  if (anyNA(idx_pc))
    stop("a recorded superitem was dropped during the spread refit")
  out <- data.frame(item = cmp$item[idx_pc], m = fit$m[idx_fit],
                    eligible = eligible,
                    spread = cmp$spread[idx_pc], se = cmp$spread_se[idx_pc],
                    lub = ifelse(eligible,
                      unname(.spread_lub[as.character(fit$m[idx_fit])]), NA_real_))
  pc_refit <- list(est = pc, repeated_ids = fit$repeated_ids,
                   person = fit$person)
  reference_df <- .dif_refit_df(pc_refit)
  out$t <- .wald_ratio(out$spread - out$lub, out$se)
  out$df <- ifelse(out$eligible, reference_df, NA_real_)
  if (is.na(reference_df)) {
    out$se[] <- NA_real_
    out$t[] <- NA_real_
  }
  out$p <- ifelse(out$eligible & is.finite(out$t) & !is.na(out$df),
                  stats::pt(out$t, df = out$df), NA_real_)
  out$p_adj <- NA_real_
  use <- out$eligible & is.finite(out$p)
  out$p_adj[use] <- stats::p.adjust(
    out$p[use], method = p_adjust, n = sum(out$eligible))
  out$below_bound <- ifelse(out$eligible & is.finite(out$spread),
                            out$spread < out$lub, NA)
  out$dependent <- ifelse(use, out$p_adj < alpha, NA)
  rownames(out) <- NULL
  attr(out, "alpha") <- alpha
  attr(out, "p_adjust") <- p_adjust
  if (is.na(reference_df)) attr(out, "note") <- paste(
    "spread inference is withheld because the principal-components refit",
    "does not carry enough independent-person support for its sandwich",
    "covariance")
  class(out) <- c("rasch_spread", "data.frame")
  out
}

#' @export
print.rasch_spread <- function(x, ...) {
  cat(sprintf(
    "Spread-parameter screen (Andrich 1985): one-sided evidence below the binomial bound (%s adjustment; alpha %.3f)\n",
    attr(x, "p_adjust") %||% "holm", attr(x, "alpha") %||% 0.05))
  d <- as.data.frame(x)
  names(d)[names(d) == "lub"] <- "bound"
  print(.fmt_df(d), row.names = FALSE)
  if (!is.null(attr(x, "note"))) cat("Note:", attr(x, "note"), "\n")
  invisible(x)
}
