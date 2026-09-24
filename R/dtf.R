# ---------------------------------------------------------------------------
# Differential test and bundle functioning from a resolved calibration.
#
# After resolve_dif() (or split_items()) each group has its own copy of every
# split item and shares the unsplit items, which act as the anchors that put
# the groups' copies on one scale. Between a group and the reference the
# test therefore differs only through its split copies, and every quantity
# here is a function of the fitted thresholds: the item shifts are linear
# contrasts, the score-to-measure and expected-score differences are smooth
# functions of them, and the delta method with the calibration covariance
# gives every standard error.
# ---------------------------------------------------------------------------

# Expected test score at one location over an item set, with its slope in
# theta and its gradient in the thresholds. For the PCM the category
# probability is proportional to exp(x theta - sum_{k <= x} tau_k), so
# d E / d tau_k = -Cov(X, 1[X >= k]) and d E / d theta = Var(X). ids_list
# holds each item's threshold row ids in the calibration covariance.
.dtf_expected <- function(theta, tau_list, ids_list, n_par) {
  total <- 0; slope <- 0; grad <- numeric(n_par)
  for (i in seq_along(tau_list)) {
    mo <- item_moments(theta, tau_list[[i]])
    x <- seq_along(mo$P) - 1L
    total <- total + mo$E; slope <- slope + mo$V
    # tail sums over x >= k for k = 1..m
    tail_xp <- rev(cumsum(rev(x * mo$P)))[-1L]
    tail_p <- rev(cumsum(rev(mo$P)))[-1L]
    grad[ids_list[[i]]] <- -(tail_xp - mo$E * tail_p)
  }
  list(E = total, V = slope, grad = grad)
}

# Root of an expected-score equation over one item set.
.dtf_root <- function(f, interval) {
  stats::uniroot(f, interval, tol = 1e-10, maxiter = 1000L)$root
}

# The location a person needs on item set `b` to have the expected score
# that location `theta` gives on item set `a`, with the gradient of the
# harder-positive logit difference theta - theta_b in the thresholds.
.dtf_logit_shift <- function(theta, a, b, interval) {
  ea <- .dtf_expected(theta, a$tau, a$ids, a$n_par)
  th_b <- .dtf_root(function(t)
    .dtf_expected(t, b$tau, b$ids, b$n_par)$E - ea$E,
    range(c(interval, theta - 5, theta + 5)))
  eb <- .dtf_expected(th_b, b$tau, b$ids, b$n_par)
  # theta_b solves E_b(theta_b) = E_a(theta): d theta_b / d tau =
  # (grad E_a - grad E_b) / V_b, and the reported shift is theta - theta_b
  list(shift = theta - th_b, grad = (eb$grad - ea$grad) / eb$V)
}

.dtf_se <- function(grad, cv) {
  if (is.null(cv)) return(NA_real_)
  v <- drop(t(grad) %*% cv %*% grad)
  if (!is.finite(v)) NA_real_ else sqrt(max(v, 0))
}

.dtf_groups <- function(fit, by) {
  if (is.null(by)) stop("`by` must name the splitting factor", call. = FALSE)
  n <- nrow(fit$X)
  if (!is.null(dim(by)) || !is.atomic(by) || is.complex(by))
    stop("`by` must name person factor(s) nominated in the fit or be an ",
         "ordinary grouping vector", call. = FALSE)
  if (is.character(by) && length(by) < n && length(by) >= 1L) {
    if (length(by) == 1L && !by %in% names(fit$factors) && grepl(":", by))
      by <- strsplit(by, ":", fixed = TRUE)[[1]]
    unknown <- setdiff(by, names(fit$factors))
    if (length(unknown))
      stop("'", paste(unknown, collapse = "', '"),
           "' is not a person factor nominated in the fit", call. = FALSE)
    grp <- if (length(by) == 1L) fit$factors[[by]] else
      .factor_cells(fit$factors[by], sep = ":")
    label <- paste(by, collapse = ":")
  } else {
    if (length(by) != n)
      stop("`by` must have one entry per person (", n, ")", call. = FALSE)
    grp <- by; label <- "group"
  }
  if (is.character(grp) || is.factor(grp)) {
    grp <- .role_text_values(grp)
    grp[!is.na(grp) & !nzchar(grp)] <- NA_character_
  }
  grp <- factor(grp)
  if (nlevels(grp) < 2L)
    stop("the grouping needs at least two levels", call. = FALSE)
  list(grp = grp, label = label)
}

#' Differential test and bundle functioning
#'
#' Measures how much a test, or a bundle of its items, functions differently
#' for one group of persons than for a reference group once the differential
#' item functioning has been resolved. The item shifts a split leaves behind
#' are combined into test-level and bundle-level differences in logits and in
#' raw-score units, each with a standard error from the calibration
#' covariance.
#'
#' @details
#' A resolved calibration (\code{\link{resolve_dif}} or
#' \code{\link{split_items}}) carries one copy of each split item per group
#' and a single copy of every unsplit item. The unsplit items are the
#' anchors: they are the items the resolution judged invariant, and they
#' place the groups' copies on one scale. Every difference reported here is
#' relative to that anchoring. A shift is positive when the test is harder
#' for the group than for the reference.
#'
#' For group \eqn{g} and reference \eqn{h} with item sets \eqn{I_g} and
#' \eqn{I_h} (identical apart from the split copies), the function reports:
#' \describe{
#'   \item{item shifts}{\eqn{\delta_{ig}-\delta_{ih}} for each source item,
#'   zero by construction for an anchor, with a Wald test from the
#'   threshold covariance.}
#'   \item{mean shift}{the average item shift over the test, the amount by
#'   which the group's copy of the test is harder on average.}
#'   \item{score-to-measure differences}{\eqn{\theta_g(r)-\theta_h(r)},
#'   the difference between the measures the two item sets assign to the
#'   same raw score \eqn{r}, from the expected-score equations
#'   \eqn{\sum_{i\in I_g}E_i(\theta)=r}. This is the bias a person in the
#'   group would carry if scored on the reference calibration.}
#'   \item{curves}{at each location \eqn{\theta}, the expected-score
#'   difference \eqn{T_h(\theta)-T_g(\theta)} and the logit difference
#'   \eqn{\theta-\tilde\theta} where \eqn{T_h(\tilde\theta)=T_g(\theta)}.}
#'   \item{test summaries}{the curves averaged over the group's own
#'   persons at their estimated locations, signed (differences may cancel)
#'   and unsigned (they may not), in logits, in score units and as a
#'   percentage of the score range, with the largest score-unit difference
#'   among those persons.}
#' }
#' Standard errors use the delta method: an item shift and the mean shift
#' are linear in the thresholds, and for the score-to-measure and curve
#' quantities \eqn{\partial E_i/\partial\tau_{ik}=-\mathrm{Cov}(X_i,
#' 1[X_i\ge k])} and \eqn{\partial E_i/\partial\theta=\mathrm{Var}(X_i)}.
#' The unsigned summaries carry a delta-method standard error but no test:
#' a folded difference has no null distribution at zero.
#'
#' A bundle (Douglas, Roussos and Stout 1996) is a named set of source
#' items. Its table gives the mean of its members' shifts with a Wald test,
#' a homogeneity chi-square on one fewer degree of freedom than the bundle
#' has members (do the members shift by the same amount?), and the signed
#' and unsigned expected-score differences over the bundle alone, as a
#' percentage of the bundle's score range. \code{\link{dif_anova}} tests
#' whether a bundle functions differently; this function sizes the
#' difference.
#'
#' @param fit A \code{\link{resolve_dif}} result, a \code{\link{rasch}} fit
#'   whose items were split with \code{\link{split_items}}, or an ordinary
#'   fit together with \code{items} to split first.
#' @param by The splitting factor: a person factor name (or several, for
#'   their joint cells) or a grouping vector with one entry per person.
#'   Taken from the resolution when \code{fit} is a \code{resolve_dif}
#'   result that split by one factor.
#' @param items Items to split by \code{by} before measuring, for an
#'   ordinary fit.
#' @param bundles Optional named list of source-item vectors, each with at
#'   least two items.
#' @param reference The level of \code{by} the other groups are compared
#'   with; the first level by default.
#' @param p_adjust Multiplicity adjustment over the item shifts and,
#'   separately, the bundle shifts; a method of
#'   \code{\link[stats]{p.adjust}}.
#' @param alpha Significance level for the adjusted probabilities.
#' @param grid Number of locations on which the curves are evaluated,
#'   spanning the fitted persons' locations.
#' @return A list of class \code{"rasch_dtf"} with tables \code{test} (one
#'   row per group), \code{items} (one row per source item and group),
#'   \code{bundles} (\code{NULL} unless requested), \code{scores} (the
#'   score-to-measure differences) and \code{curves}; \code{anchors}, the
#'   unsplit items; \code{groups}, \code{reference}, \code{by} and
#'   \code{notes}.
#' @references
#' Chalmers, R. P. (2018). Model-based measures for detecting and
#' quantifying response bias. Psychometrika, 83(3), 696--732.
#'
#' Douglas, J. A., Roussos, L. A. and Stout, W. (1996). Item-bundle DIF
#' hypothesis testing: Identifying suspect bundles and assessing their
#' differential functioning. Journal of Educational Measurement, 33(4),
#' 465--484.
#'
#' Andrich, D. and Hagquist, C. (2015). Real and artificial differential
#' item functioning in polytomous items. Educational and Psychological
#' Measurement, 75(2), 185--207.
#' @examples
#' set.seed(1); n <- 600
#' d <- seq(-2, 2, length.out = 8); g <- rep(c("a", "b"), each = n / 2)
#' sh <- matrix(0, n, 8); sh[g == "b", 2:3] <- 0.8
#' X <- matrix(rbinom(n * 8, 1, plogis(outer(rnorm(n), d, "-") - sh)), n, 8)
#' colnames(X) <- paste0("I", 1:8)
#' fit <- rasch(data.frame(X, grp = g), factors = "grp")
#' d2 <- dtf(fit, by = "grp", items = c("I2", "I3"),
#'           bundles = list(pair = c("I2", "I3")))
#' d2
#' d2$scores
#' plot_dtf(d2)
#' @seealso \code{\link{dif_anova}} for detecting item and bundle DIF,
#'   \code{\link{resolve_dif}} for the resolution this function measures
#'   and \code{\link{plot_dtf}} for the curves.
#' @export
dtf <- function(fit, by = NULL, items = NULL, bundles = NULL,
                reference = NULL, p_adjust = "holm", alpha = 0.05,
                grid = 61L) {
  .check_dif_args(alpha, p_adjust)
  .check_whole(grid, "grid", 5)
  notes <- character(0)
  from_resolution <- inherits(fit, "rasch_resolve_dif")
  if (from_resolution) {
    res <- fit; fit <- res$fit
    if (is.null(by)) {
      fac <- unique(res$splits$factor)
      if (length(fac) > 1L)
        stop("the resolution split by several factors (",
             paste(fac, collapse = ", "), "); name the one to compare in `by`",
             call. = FALSE)
      if (length(fac) == 1L) by <- fac
      else if (!is.null(fit$factors) && ncol(fit$factors) == 1L)
        by <- names(fit$factors)
    }
    if (!is.null(items))
      stop("`items` applies to an unsplit fit; the resolution has already ",
           "split its items", call. = FALSE)
  }
  if (!inherits(fit, "rasch") ||
      inherits(fit, c("rasch_mfrm", "rasch_efrm", "rasch_btl")))
    stop("dtf needs an ordinary rasch fit or a resolve_dif result",
         call. = FALSE)
  if (!isTRUE(fit$est$converged))
    stop("the fitted calibration did not converge; DTF is unavailable",
         call. = FALSE)
  if (!is.null(fit$disc) && length(unique(fit$disc)) > 1L)
    stop("DTF needs items sharing one discrimination", call. = FALSE)
  if (is.null(by) && !is.null(fit$factors) && ncol(fit$factors) == 1L)
    by <- names(fit$factors)
  g <- .dtf_groups(fit, by)
  grp <- g$grp
  if (!is.null(items)) {
    if (!is.null(fit$split_map))
      stop("the fit already carries split items; pass `items` with an ",
           "unsplit fit", call. = FALSE)
    fit <- split_items(fit, items, by = grp)
  } else if (is.null(fit$split_map) && !from_resolution) {
    stop("no item is split, so the groups share every item and there is ",
         "nothing to compare; resolve DIF first or name `items` to split",
         call. = FALSE)
  }
  if (is.null(reference)) reference <- levels(grp)[1L]
  if (length(reference) != 1L || !reference %in% levels(grp))
    stop("`reference` must be one level of the grouping: ",
         paste(levels(grp), collapse = ", "), call. = FALSE)
  reference <- as.character(reference)
  focal <- setdiff(levels(grp), reference)

  # each source item must have exactly one calibrated copy answered by each
  # group; a source with no copy for some group (its level was omitted from
  # the split) leaves the comparison
  map <- .split_source_map(fit)
  sources <- unique(unname(map))
  X <- fit$X
  copy_of <- function(level) vapply(sources, function(s) {
    cand <- names(map)[map == s]
    rows <- !is.na(grp) & grp == level
    used <- cand[colSums(!is.na(X[rows, cand, drop = FALSE])) > 0L]
    if (length(used) > 1L)
      stop("item ", s, " has several copies answered by group '", level,
           "' (", paste(used, collapse = ", "), "); the fit is not split by ",
           "this grouping", call. = FALSE)
    if (!length(used)) NA_character_ else used
  }, "")
  copies <- lapply(c(reference, focal), copy_of)
  names(copies) <- c(reference, focal)
  present <- Reduce(`&`, lapply(copies, function(z) !is.na(z)))
  dropped <- sources[!present]
  if (length(dropped)) {
    notes <- c(notes, sprintf(paste(
      "item(s) without a calibrated copy for every group were left out of",
      "the comparison: %s"), paste(dropped, collapse = ", ")))
    sources <- sources[present]
    copies <- lapply(copies, function(z) z[present])
  }
  if (length(sources) < 2L)
    stop("fewer than two items are calibrated for every group", call. = FALSE)
  anchors <- sources[vapply(sources, function(s)
    length(unique(vapply(copies, `[[`, "", s))) == 1L, TRUE)]
  if (!length(anchors))
    stop("every item is split by the grouping, so no anchor places the ",
         "groups on one scale; keep some invariant items unsplit",
         call. = FALSE)

  thr <- fit$thresholds
  n_par <- nrow(thr)
  cv <- fit$est$cov_tau
  cov_ok <- is.matrix(cv) && identical(dim(cv), c(n_par, n_par)) &&
    all(is.finite(cv))
  if (!cov_ok) {
    cv <- NULL
    notes <- c(notes, paste(
      "the threshold covariance is unavailable, so standard errors and",
      "tests are withheld"))
  }
  item_index <- function(name) match(name, fit$items$item)
  ids_of <- function(name) thr$id[thr$item == item_index(name)]
  loc_of <- function(name) fit$items$location[item_index(name)]
  m_of <- function(name) fit$m[item_index(name)]
  # threshold-mean contrast for a location, as a vector over the covariance
  loc_vec <- function(name) {
    v <- numeric(n_par); ids <- ids_of(name); v[ids] <- 1 / length(ids); v
  }
  set_of <- function(level) {
    nm <- copies[[level]]
    list(tau = fit$tau_list[item_index(nm)], ids = lapply(nm, ids_of),
         n_par = n_par, names = nm)
  }
  ref_set <- set_of(reference)
  M <- sum(vapply(sources, function(s) m_of(copies[[reference]][s]), 1L))
  interval <- .person_root_interval(fit$tau_list, 1)
  theta_hat <- fit$person$theta
  grid_range <- range(theta_hat[is.finite(theta_hat)])
  if (!all(is.finite(grid_range))) grid_range <- range(unlist(fit$tau_list))
  grid_theta <- seq(grid_range[1] - 0.5, grid_range[2] + 0.5,
                    length.out = grid)

  bundles <- .dtf_bundles(bundles, sources)

  item_rows <- list(); test_rows <- list(); bundle_rows <- list()
  score_rows <- list(); curve_rows <- list()
  for (lv in focal) {
    cp <- copies[[lv]]; cr <- copies[[reference]]
    split <- cp != cr
    if (sum(vapply(sources, function(s) m_of(cp[s]), 1L)) != M)
      stop("internal error: the score range differs between groups")
    # item shifts as threshold-mean contrasts
    vecs <- lapply(sources, function(s)
      if (cp[s] == cr[s]) numeric(n_par) else loc_vec(cp[s]) - loc_vec(cr[s]))
    shifts <- vapply(sources, function(s) loc_of(cp[s]) - loc_of(cr[s]), 0)
    ses <- vapply(vecs, .dtf_se, 0, cv = cv)
    ses[!split] <- if (cov_ok) 0 else NA_real_
    z <- ifelse(split, .wald_ratio(shifts, ses), NA_real_)
    p <- 2 * stats::pnorm(-abs(z))
    item_rows[[lv]] <- data.frame(
      group = lv, item = sources, item_group = unname(cp),
      item_reference = unname(cr),
      location_group = vapply(cp, loc_of, 0),
      location_reference = vapply(cr, loc_of, 0),
      shift = shifts, se = ses, z = z, p = p, split = split,
      stringsAsFactors = FALSE, row.names = NULL)
    mean_vec <- Reduce(`+`, vecs) / length(sources)
    shift_mean <- mean(shifts)
    se_mean <- .dtf_se(mean_vec, cv)
    # score-to-measure: the measure each item set assigns to a raw score
    fs <- set_of(lv)
    sc <- lapply(seq_len(M - 1L), function(r) {
      th_r <- .dtf_root(function(t)
        .dtf_expected(t, ref_set$tau, ref_set$ids, n_par)$E - r, interval)
      th_g <- .dtf_root(function(t)
        .dtf_expected(t, fs$tau, fs$ids, n_par)$E - r, interval)
      er <- .dtf_expected(th_r, ref_set$tau, ref_set$ids, n_par)
      eg <- .dtf_expected(th_g, fs$tau, fs$ids, n_par)
      grad <- -eg$grad / eg$V + er$grad / er$V
      c(theta_reference = th_r, theta_group = th_g,
        shift = th_g - th_r, se = .dtf_se(grad, cv))
    })
    sc <- as.data.frame(do.call(rbind, sc))
    score_rows[[lv]] <- data.frame(group = lv, score = seq_len(M - 1L), sc,
                                   stringsAsFactors = FALSE, row.names = NULL)
    # curves over the location grid
    cu <- lapply(grid_theta, function(th) {
      eg <- .dtf_expected(th, fs$tau, fs$ids, n_par)
      er <- .dtf_expected(th, ref_set$tau, ref_set$ids, n_par)
      ls <- .dtf_logit_shift(th, fs, ref_set, interval)
      c(theta = th, expected_reference = er$E, expected_group = eg$E,
        shift_score = er$E - eg$E, se_score = .dtf_se(er$grad - eg$grad, cv),
        shift_logit = ls$shift, se_logit = .dtf_se(ls$grad, cv))
    })
    cu <- as.data.frame(do.call(rbind, cu))
    curve_rows[[lv]] <- data.frame(group = lv, cu, stringsAsFactors = FALSE,
                                   row.names = NULL)
    # summaries over the group's own persons at their estimated locations
    th_p <- theta_hat[!is.na(grp) & grp == lv & is.finite(theta_hat)]
    n_p <- length(th_p)
    summ <- .dtf_summaries(th_p, fs, ref_set, interval, cv, M,
                           function(th, a, b) .dtf_logit_shift(th, a, b, interval))
    if (!any(split)) {
      # the group shares every item with the reference: every difference is
      # zero by construction, and root-finding noise must not become a test
      summ[] <- lapply(summ, function(v) if (is.na(v)) v else 0)
      score_rows[[lv]]$shift <- 0
      score_rows[[lv]]$se <- if (cov_ok) 0 else NA_real_
      curve_rows[[lv]][c("shift_score", "shift_logit")] <- 0
      curve_rows[[lv]][c("se_score", "se_logit")] <- if (cov_ok) 0 else
        NA_real_
      notes <- c(notes, sprintf(paste(
        "group '%s' shares every item with the reference, so its",
        "differences are zero by construction"), lv))
    }
    if (n_p == 0L)
      notes <- c(notes, sprintf(paste(
        "group '%s' has no person with a finite location, so its",
        "person-weighted summaries are unavailable"), lv))
    test_rows[[lv]] <- data.frame(
      group = lv, n = n_p, shift_mean = shift_mean, se = se_mean,
      z = .wald_ratio(shift_mean, se_mean), summ,
      stringsAsFactors = FALSE, row.names = NULL)
    # bundles: mean shift, homogeneity, and the bundle's own expected score
    for (b in names(bundles)) {
      mem <- bundles[[b]]; k <- length(mem)
      idx <- match(mem, sources)
      bvec <- Reduce(`+`, vecs[idx]) / k
      b_shift <- mean(shifts[idx]); b_se <- .dtf_se(bvec, cv)
      hom <- c(chisq = NA_real_, df = k - 1, p = NA_real_)
      if (cov_ok && any(split[idx])) {
        Vb <- do.call(cbind, lapply(vecs[idx], function(v) v))
        S <- t(Vb) %*% cv %*% Vb
        C <- cbind(-1, diag(k - 1L))
        d <- C %*% shifts[idx]; W <- C %*% S %*% t(C)
        q <- tryCatch(drop(t(d) %*% solve(W, d)), error = function(e) NA_real_)
        if (is.finite(q))
          hom <- c(chisq = q, df = k - 1, p = stats::pchisq(q, k - 1,
                                                          lower.tail = FALSE))
      }
      bset <- list(tau = fs$tau[idx], ids = fs$ids[idx], n_par = n_par)
      rset <- list(tau = ref_set$tau[idx], ids = ref_set$ids[idx],
                   n_par = n_par)
      Mb <- sum(vapply(mem, function(s) m_of(cp[s]), 1L))
      bs <- .dtf_summaries(th_p, bset, rset, interval, cv, Mb, NULL)
      bundle_rows[[length(bundle_rows) + 1L]] <- data.frame(
        group = lv, bundle = b, n_items = k, shift_mean = b_shift,
        se = b_se, z = .wald_ratio(b_shift, b_se),
        chisq_hom = hom[["chisq"]], df_hom = hom[["df"]], p_hom = hom[["p"]],
        bs[c("sDBF_score", "se_sDBF_score", "uDBF_score", "se_uDBF_score",
             "sDBF_pct", "uDBF_pct")],
        stringsAsFactors = FALSE, row.names = NULL)
    }
  }
  items_tab <- do.call(rbind, item_rows); rownames(items_tab) <- NULL
  items_tab$p_adj <- NA_real_
  tested <- is.finite(items_tab$p)
  if (any(tested))
    items_tab$p_adj[tested] <- .p_adjust_family(items_tab$p[tested],
                                                method = p_adjust)
  items_tab$significant <- is.finite(items_tab$p_adj) & items_tab$p_adj < alpha
  test_tab <- do.call(rbind, test_rows); rownames(test_tab) <- NULL
  test_tab$p <- 2 * stats::pnorm(-abs(test_tab$z))
  test_tab$p_sDTF_logit <- 2 * stats::pnorm(-abs(
    .wald_ratio(test_tab$sDTF_logit, test_tab$se_sDTF_logit)))
  test_tab$p_sDTF_score <- 2 * stats::pnorm(-abs(
    .wald_ratio(test_tab$sDTF_score, test_tab$se_sDTF_score)))
  test_tab <- test_tab[, c("group", "n", "shift_mean", "se", "z", "p",
                           "sDTF_logit", "se_sDTF_logit", "p_sDTF_logit",
                           "uDTF_logit", "se_uDTF_logit",
                           "sDTF_score", "se_sDTF_score", "p_sDTF_score",
                           "uDTF_score", "se_uDTF_score",
                           "sDTF_pct", "uDTF_pct", "max_score")]
  bundles_tab <- NULL
  if (length(bundle_rows)) {
    bundles_tab <- do.call(rbind, bundle_rows); rownames(bundles_tab) <- NULL
    bundles_tab$p <- 2 * stats::pnorm(-abs(bundles_tab$z))
    bundles_tab$p_adj <- NA_real_
    tested <- is.finite(bundles_tab$p)
    if (any(tested))
      bundles_tab$p_adj[tested] <- .p_adjust_family(bundles_tab$p[tested],
                                                    method = p_adjust)
    bundles_tab$significant <- is.finite(bundles_tab$p_adj) &
      bundles_tab$p_adj < alpha
    bundles_tab$p_sDBF_score <- 2 * stats::pnorm(-abs(
      .wald_ratio(bundles_tab$sDBF_score, bundles_tab$se_sDBF_score)))
    bundles_tab <- bundles_tab[, c(
      "group", "bundle", "n_items", "shift_mean", "se", "z", "p", "p_adj",
      "significant", "chisq_hom", "df_hom", "p_hom",
      "sDBF_score", "se_sDBF_score", "p_sDBF_score", "uDBF_score",
      "se_uDBF_score", "sDBF_pct", "uDBF_pct")]
  }
  scores_tab <- do.call(rbind, score_rows); rownames(scores_tab) <- NULL
  curves_tab <- do.call(rbind, curve_rows); rownames(curves_tab) <- NULL
  out <- list(test = test_tab, items = items_tab, bundles = bundles_tab,
              scores = scores_tab, curves = curves_tab,
              anchors = anchors, dropped = dropped,
              by = g$label, reference = reference, groups = focal,
              max_score = M, n_items = length(sources),
              alpha = alpha, p_adjust = p_adjust, notes = notes)
  out <- .tag_tables(out)
  class(out) <- "rasch_dtf"
  out
}

# Person-weighted signed and unsigned differences between two item sets at
# the locations th_p, in score units (with the delta-method gradient of the
# signed and unsigned means) and, when a logit-shift function is supplied,
# in logits. Names follow the test-level labels; the bundle table renames.
.dtf_summaries <- function(th_p, a, b, interval, cv, M, logit_fun) {
  na <- c(sDTF_logit = NA_real_, se_sDTF_logit = NA_real_,
          uDTF_logit = NA_real_, se_uDTF_logit = NA_real_,
          sDTF_score = NA_real_, se_sDTF_score = NA_real_,
          uDTF_score = NA_real_, se_uDTF_score = NA_real_,
          sDTF_pct = NA_real_, uDTF_pct = NA_real_, max_score = NA_real_)
  rename <- function(x) {
    if (is.null(logit_fun)) {
      x <- x[c("sDTF_score", "se_sDTF_score", "uDTF_score", "se_uDTF_score",
               "sDTF_pct", "uDTF_pct")]
      names(x) <- sub("DTF", "DBF", names(x), fixed = TRUE)
    }
    as.data.frame(as.list(x))
  }
  if (!length(th_p)) return(rename(na))
  # persons sharing a location share every quantity below
  u <- sort(unique(th_p)); w <- as.numeric(table(factor(th_p, levels = u)))
  w <- w / sum(w)
  ds <- matrix(0, length(u), a$n_par); dl <- ds
  s_score <- s_logit <- numeric(length(u))
  for (j in seq_along(u)) {
    ea <- .dtf_expected(u[j], a$tau, a$ids, a$n_par)
    eb <- .dtf_expected(u[j], b$tau, b$ids, b$n_par)
    s_score[j] <- eb$E - ea$E; ds[j, ] <- eb$grad - ea$grad
    if (!is.null(logit_fun)) {
      ls <- logit_fun(u[j], a, b)
      s_logit[j] <- ls$shift; dl[j, ] <- ls$grad
    }
  }
  signed <- function(s, d) c(sum(w * s), .dtf_se(drop(w %*% d), cv))
  unsigned <- function(s, d) c(sum(w * abs(s)),
                               .dtf_se(drop((w * sign(s)) %*% d), cv))
  ss <- signed(s_score, ds); us <- unsigned(s_score, ds)
  out <- na
  out[c("sDTF_score", "se_sDTF_score")] <- ss
  out[c("uDTF_score", "se_uDTF_score")] <- us
  out[c("sDTF_pct", "uDTF_pct")] <- 100 * c(ss[1], us[1]) / M
  out["max_score"] <- max(abs(s_score))
  if (!is.null(logit_fun)) {
    sl <- signed(s_logit, dl); ul <- unsigned(s_logit, dl)
    out[c("sDTF_logit", "se_sDTF_logit")] <- sl
    out[c("uDTF_logit", "se_uDTF_logit")] <- ul
  }
  rename(out)
}

.dtf_bundles <- function(bundles, sources) {
  if (is.null(bundles)) return(list())
  if (!is.list(bundles) || !length(bundles) || is.data.frame(bundles))
    stop("`bundles` must be a named list of item-name vectors", call. = FALSE)
  nms <- names(bundles)
  if (is.null(nms) || anyNA(nms) || any(!nzchar(nms)))
    stop("every bundle needs a name", call. = FALSE)
  if (anyDuplicated(nms))
    stop("bundle name(s) used more than once: ",
         paste(unique(nms[duplicated(nms)]), collapse = ", "), call. = FALSE)
  out <- lapply(nms, function(b) {
    members <- bundles[[b]]
    if (!is.character(members) || anyNA(members) || !length(members))
      stop("bundle '", b, "' must name its items", call. = FALSE)
    members <- unique(members)
    missing <- setdiff(members, sources)
    if (length(missing))
      stop("bundle '", b, "' names item(s) not calibrated for every group: ",
           paste(missing, collapse = ", "), call. = FALSE)
    if (length(members) < 2L)
      stop("bundle '", b, "' needs at least two items", call. = FALSE)
    members
  })
  names(out) <- nms
  out
}

#' @export
print.rasch_dtf <- function(x, ...) {
  cat(sprintf(paste0("Differential test functioning by %s (reference: %s; ",
                     "%d of %d items anchor the groups)\n"),
              x$by, x$reference, length(x$anchors), x$n_items))
  cat("Positive values: harder for the group than for the reference.\n")
  show <- x$test[, c("group", "n", "shift_mean", "se", "p",
                     "sDTF_logit", "uDTF_logit", "sDTF_score", "uDTF_score",
                     "sDTF_pct", "uDTF_pct")]
  print(.fmt_df(show), row.names = FALSE)
  if (!is.null(x$bundles)) {
    cat("Bundles:\n")
    showb <- x$bundles[, c("group", "bundle", "n_items", "shift_mean", "se",
                           "p_adj", "significant", "chisq_hom", "p_hom",
                           "sDBF_score", "uDBF_score", "sDBF_pct")]
    print(.fmt_df(showb), row.names = FALSE)
  }
  sp <- x$items[x$items$split, , drop = FALSE]
  if (nrow(sp)) {
    cat("Split items:\n")
    print(.fmt_df(sp[, c("group", "item", "shift", "se", "p_adj",
                         "significant")]), row.names = FALSE)
  }
  if (length(x$notes)) cat("Notes:", paste(x$notes, collapse = "; "), "\n")
  invisible(x)
}

#' Plot differential test functioning
#'
#' Draws the expected-score curves of the reference and one other group
#' over the location scale, with the difference between them below.
#'
#' @param x A \code{\link{dtf}} result.
#' @param group The group to draw; the first compared group by default.
#' @param ... Further arguments passed to \code{\link{plot}}.
#' @return \code{x}, invisibly.
#' @examples
#' set.seed(2)
#' n <- 400
#' d <- seq(-1.5, 1.5, length.out = 8)
#' g <- rep(c("a", "b"), each = n / 2)
#' sh <- matrix(0, n, 8); sh[g == "b", 2:3] <- 0.8
#' X <- matrix(rbinom(n * 8, 1, plogis(outer(rnorm(n), d, "-") - sh)), n, 8)
#' colnames(X) <- paste0("I", 1:8)
#' fit <- rasch(data.frame(X, grp = g), factors = "grp")
#' plot_dtf(dtf(fit, by = "grp", items = c("I2", "I3")))
#' @export
plot_dtf <- function(x, group = NULL, ...) {
  if (is.null(group)) group <- x$groups[1L]
  if (length(group) != 1L || !group %in% x$groups)
    stop("`group` must be one of: ", paste(x$groups, collapse = ", "),
         call. = FALSE)
  cu <- x$curves[x$curves$group == group, , drop = FALSE]
  op <- graphics::par(mfrow = c(2, 1), mar = c(4, 4, 2, 1))
  on.exit(graphics::par(op))
  graphics::plot(cu$theta, cu$expected_reference, type = "l",
                 xlab = "Location (logits)", ylab = "Expected score",
                 ylim = c(0, x$max_score),
                 main = sprintf("Expected test score: %s and %s",
                                x$reference, group), ...)
  graphics::lines(cu$theta, cu$expected_group, lty = 2)
  graphics::legend("topleft", legend = c(x$reference, group), lty = 1:2,
                   bty = "n")
  band <- 1.96 * cu$se_score
  ylim <- range(c(cu$shift_score - band, cu$shift_score + band, 0),
                na.rm = TRUE)
  graphics::plot(cu$theta, cu$shift_score, type = "n",
                 xlab = "Location (logits)",
                 ylab = "Score difference (reference - group)", ylim = ylim)
  if (all(is.finite(band)))
    graphics::polygon(c(cu$theta, rev(cu$theta)),
                      c(cu$shift_score - band, rev(cu$shift_score + band)),
                      col = "grey85", border = NA)
  graphics::abline(h = 0, col = "grey60")
  graphics::lines(cu$theta, cu$shift_score)
  invisible(x)
}
