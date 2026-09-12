# rasch :: fit statistics
# ===========================================================================
# Test-of-fit statistics for the pairwise analysis: standardised residuals;
# the log-of-mean-square fit residual of Andrich & Marais (2019, ch. 23)
# for items and persons (with its
# untransformed "natural" form); infit and outfit mean squares with
# Wilson-Hilferty standardisations (the mean squares' null variance is
# about 2/N and the mean square reads as relative ICC slope -- Wu & Adams
# 2013, JAM 14(4), whose claimed N-independent null for the standardised
# forms does NOT survive parameter estimation; see fit_bootstrap());
# the class-interval ANOVA item-fit F;
# the item-trait interaction chi-square; the person separation index with
# and without extremes; Cronbach's alpha; targeting; and the test
# information function.
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

# ---------------------------------------------------------------------------
# The fit residual of Andrich & Marais (2019, ch. 23; see also Andrich 1988).
# The observed cells of non-extreme persons
# carry C - (N + P) model-testing degrees of freedom, where P item parameters
# and N person locations were estimated (with complete dichotomous data this
# is (N-1)(I-1) - (m-1)); apportioned equally, each cell carries
# f_cell = df_total / C. Summing z^2 over an item's (person's) observed
# cells gives Y^2 with E[Y^2] = f (the summed cell df) and model variance
# V[Y^2] = sum(C4/V^2 - 1) (the dichotomous form 2 tanh((b-d)/2) sinh(b-d)
# generalised to ordered categories). The reported fit residual is the
# symmetrising log-of-mean-square transform
#   T2 = f (ln Y^2 - ln f) / sqrt(V[Y^2])    (A&M 2019, eq. 23.14)
# with expectation 0 and variance 1 under fit: negative = too deterministic /
# over-discriminating (Guttman-like), positive = erratic /
# under-discriminating. The conventional flagging value is |T2| > 2.5
# (Andrich & Marais 2019, ch. 15). The untransformed statistic
#   T1 = (Y^2 - f) / sqrt(V[Y^2])
# is the "natural" fit residual and is kept alongside as natural_resid.
# Extreme persons are excluded throughout, exactly as they are set aside
# from the calibrating sample, so their person fit residual is NA and they
# contribute nothing to item fit or to the degrees of freedom.
# ---------------------------------------------------------------------------
.fitres_df <- function(Z, extreme, n_parameters) {
  obs <- !is.na(Z) & !extreme
  C <- sum(obs)
  N_ne <- sum(rowSums(obs) > 0)
  df_total <- C - N_ne - n_parameters
  f_cell <- if (C > 0 && df_total > 0) df_total / C else NA_real_
  list(obs = obs, f_cell = f_cell,
       f_item = f_cell * colSums(obs), f_person = f_cell * rowSums(obs))
}

# One margin's fit residuals from cell sums: y2 = sum z^2, v = sum of
# per-cell variances C4/V^2 - 1, f = summed cell df, n = cells used.
.fitres_transform <- function(y2, v, f, n, min_cells = 2L) {
  ok <- n >= min_cells & !is.na(f) & f > 0 & v > 1e-8 & y2 > 0
  natural <- fit_resid <- rep(NA_real_, length(y2))
  natural[ok] <- (y2[ok] - f[ok]) / sqrt(v[ok])
  fit_resid[ok] <- f[ok] * (log(y2[ok]) - log(f[ok])) / sqrt(v[ok])
  list(fit_resid = fit_resid, natural = natural, df = ifelse(ok, f, NA_real_))
}

.fitres <- function(Z, mo, extreme, n_parameters) {
  dfs <- .fitres_df(Z, extreme, n_parameters)
  z2 <- Z^2; z2[!dfs$obs] <- NA
  vcell <- mo$M4 / mo$V^2 - 1; vcell[!dfs$obs] <- NA
  it <- .fitres_transform(colSums(z2, na.rm = TRUE),
                        colSums(vcell, na.rm = TRUE),
                        dfs$f_item, colSums(dfs$obs))
  pe <- .fitres_transform(rowSums(z2, na.rm = TRUE),
                        rowSums(vcell, na.rm = TRUE),
                        dfs$f_person, rowSums(dfs$obs), min_cells = 3L)
  list(items = it, persons = pe, f_cell = dfs$f_cell)
}

# Class-interval ANOVA item fit (Andrich & Marais 2019, ch. 15): per item,
# a one-way analysis of variance of the standardised residuals
# over the class intervals. Under fit the interval means share a common
# zero mean, so the between-interval F on (G - 1, n - G) degrees of freedom
# tests the same item-trait interaction as the chi-square but through the
# ANOVA calibration. Reported with Holm and Bonferroni familywise
# adjustments across items.
.item_anova <- function(Z, ci, extreme, ci_list = NULL) {
  L <- ncol(Z)
  out <- data.frame(item = colnames(Z), F_anova = NA_real_, df1 = NA_integer_,
                    df2 = NA_integer_, p = NA_real_)
  for (i in seq_len(L)) {
    ci_i <- if (is.null(ci_list)) ci else ci_list[[i]]
    sel <- which(!is.na(Z[, i]) & !is.na(ci_i) & !extreme)
    if (!length(sel)) next
    g <- ci_i[sel]; z <- Z[sel, i]
    keep <- g %in% as.integer(names(which(table(g) >= 2)))
    g <- factor(g[keep]); z <- z[keep]
    G <- nlevels(g); n <- length(z)
    if (G < 2 || n - G < 1) next
    mg <- tapply(z, g, mean); ng <- tabulate(g)
    ssb <- sum(ng * (mg - mean(z))^2)
    ssw <- sum((z - mg[g])^2)
    if (ssw <= 0) next
    out$F_anova[i] <- (ssb / (G - 1)) / (ssw / (n - G))
    out$df1[i] <- G - 1L; out$df2[i] <- n - G
    out$p[i] <- pf(out$F_anova[i], G - 1, n - G, lower.tail = FALSE)
  }
  usable <- is.finite(out$p)
  out$p_adj <- out$p_bonf <- rep(NA_real_, nrow(out))
  out$p_adj[usable] <- p.adjust(out$p[usable], method = "holm", n = L)
  out$p_bonf[usable] <- p.adjust(out$p[usable], method = "bonferroni", n = L)
  out
}

# Item fit from per-person model moments (observed cells only). Extreme-score
# persons are excluded: their measures are at a boundary, so their residuals
# are structurally near zero and would deflate the mean-squares toward
# apparent fit -- the same convention the log-of-mean-square fit residual and
# the item-fit ANOVA already use, and standard in Rasch fit reporting.
.item_fit <- function(X, Z, mo, disc = NULL, extreme = NULL) {
  L <- ncol(X)
  if (is.null(extreme)) extreme <- rep(FALSE, nrow(X))
  E2 <- .z2_expectation(mo, Z, disc)
  out <- data.frame(item = colnames(X), infit_ms = NA_real_, outfit_ms = NA_real_,
                    infit_z = NA_real_, outfit_z = NA_real_, n = NA_integer_)
  for (i in seq_len(L)) {
    ok <- which(!is.na(Z[, i]) & !extreme)
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
    out$n[i] <- n
  }
  out
}

# Person fit residuals: each person's standardised residuals across their
# observed items, summarised exactly as for items. Items with an extreme
# (minimum- or maximum-possible) total are excluded, mirroring the exclusion
# of extreme persons from item fit: their locations are at a boundary and
# their residuals are structurally near zero.
.person_fit <- function(X, Z, mo, disc = NULL, item_extreme = NULL) {
  N <- nrow(X)
  if (is.null(item_extreme)) item_extreme <- rep(FALSE, ncol(X))
  E2 <- .z2_expectation(mo, Z, disc)
  out <- data.frame(infit_ms = rep(NA_real_, N), outfit_ms = NA_real_,
                    infit_z = NA_real_, outfit_z = NA_real_)
  for (n in seq_len(N)) {
    ok <- which(!is.na(Z[n, ]) & !item_extreme)
    if (length(ok) < 3) next
    z2 <- Z[n, ok]^2
    V <- mo$V[n, ok]; C4 <- mo$M4[n, ok]; k <- length(ok)
    e2 <- E2[n, ok]
    outfit <- sum(z2) / sum(e2)
    out$outfit_ms[n] <- outfit
    out$infit_ms[n] <- sum(z2 * V) / sum(e2 * V)
    qo <- sqrt(max(sum(C4 / V^2) / k^2 - 1 / k, 1e-8))
    qi <- sqrt(max(sum(C4 - V^2) / sum(V)^2, 1e-8))
    out$outfit_z[n] <- .wh(outfit, qo)
    out$infit_z[n] <- .wh(out$infit_ms[n], qi)
  }
  out
}

# The default number of class intervals (Andrich & Marais 2019, ch. 15):
# as many intervals of at least 50 persons as the non-extreme sample
# allows, at most 10, at least 2.
.default_n_groups <- function(n_ne) max(2L, min(10L, n_ne %/% 50L))

# Allocate locations to n_groups contiguous intervals, as equal-sized as
# possible WITHOUT splitting ties: persons sharing a location are
# indistinguishable (equal raw scores give equal measures), so they belong
# to the same interval and interval sizes are generally unequal, as in the
# worked class-interval tables of Andrich & Marais (2019, ch. 13; sizes
# such as 13/20/16). A boundary falls where the cumulative count comes
# closest to each equal-share target.
.ci_allocate <- function(th, n_groups) {
  ut <- sort(unique(th))
  if (length(ut) <= n_groups) return(match(th, ut))
  # tabulate(match()) not table(factor()): two distinct doubles can share a
  # printed representation, and factor(levels = ut) then dies on
  # "duplicated" levels (person-mean locations from tapply hit this)
  cnt <- tabulate(match(th, ut), nbins = length(ut))
  cum <- cumsum(cnt)
  n <- length(th)
  b <- integer(n_groups - 1L)
  lo <- 1L
  for (gg in seq_len(n_groups - 1L)) {
    target <- n * gg / n_groups
    cand <- seq(lo, length(ut) - (n_groups - gg))
    b[gg] <- cand[which.min(abs(cum[cand] - target))]
    lo <- b[gg] + 1L
  }
  grp_of_ut <- findInterval(seq_along(ut), b + 1L) + 1L
  grp_of_ut[match(th, ut)]
}

# Class intervals over non-extreme person locations. n_groups = NULL
# applies the default rule above.
.class_intervals <- function(theta, extreme, n_groups = NULL) {
  g <- rep(NA_integer_, length(theta))
  use <- which(!is.na(theta) & !extreme)
  if (!length(use)) {
    attr(g, "n_groups") <- 0L
    return(g)
  }
  if (is.null(n_groups)) n_groups <- .default_n_groups(length(use))
  g[use] <- .ci_allocate(theta[use], n_groups)
  attr(g, "n_groups") <- max(g[use], na.rm = TRUE)
  g
}

# Class intervals compiled per item, the automatic adjustment with missing
# data (Andrich & Marais 2019, ch. 15): each item allocates the persons who
# answered it into intervals of its own, with the group count from the same
# rule applied to that item's responders. Returns a list of allocation vectors, one per item.
.class_intervals_by_item <- function(X, theta, extreme, n_groups = NULL) {
  lapply(seq_len(ncol(X)), function(i) {
    obs <- !is.na(X[, i]) & !is.na(theta) & !extreme
    g <- rep(NA_integer_, length(theta))
    if (!any(obs)) return(g)
    ng <- if (is.null(n_groups)) .default_n_groups(sum(obs)) else n_groups
    g[obs] <- .ci_allocate(theta[obs], ng)
    g
  })
}

# Item-trait interaction chi-square over class intervals, per item. ci is
# the common allocation vector; ci_list, when supplied, gives each item its
# own allocation (per-item basis). The degrees of freedom are per item:
# (number of class intervals contributing at least 2 responders to that
# item) - 1, so items with missing data are tested on the intervals they
# actually reach.
.item_trait <- function(X, mo, ci, ci_list = NULL) {
  L <- ncol(X)
  chi <- setNames(numeric(L), colnames(X))
  used <- integer(L)
  invalid <- logical(L)
  for (i in seq_len(L)) {
    ci_i <- if (is.null(ci_list)) ci else ci_list[[i]]
    G <- suppressWarnings(max(ci_i, na.rm = TRUE))
    if (!is.finite(G)) next
    for (gg in seq_len(G)) {
      sel <- which(ci_i == gg & !is.na(X[, i]))
      if (length(sel) < 2) next
      Obar <- mean(X[sel, i])
      Ebar <- mean(mo$E[sel, i]); Vbar <- mean(mo$V[sel, i])
      if (!is.finite(Obar) || !is.finite(Ebar) ||
          !is.finite(Vbar) || Vbar <= 0) {
        invalid[i] <- TRUE
        next
      }
      chi[i] <- chi[i] + length(sel) * (Obar - Ebar)^2 / Vbar
      used[i] <- used[i] + 1L
    }
  }
  # an item whose responders fall in fewer than two class intervals has no
  # estimable item-by-trait interaction: its chi-square and df are NA, not
  # a manufactured df = 1 with a valid-looking p
  df_i <- ifelse(used >= 2L & !invalid, used - 1L, NA_integer_)
  chi[is.na(df_i)] <- NA_real_
  p <- pchisq(chi, df_i, lower.tail = FALSE)
  usable <- is.finite(p)
  p_adj <- p_bonf <- rep(NA_real_, length(p))
  p_adj[usable] <- p.adjust(p[usable], method = "holm", n = L)
  p_bonf[usable] <- p.adjust(p[usable], method = "bonferroni", n = L)
  data.frame(item = colnames(X), chisq = chi, df = df_i, p = p,
             p_adj = p_adj, p_bonf = p_bonf)
}

# Correlation that degrades to NA (rather than erroring) when fewer than 3
# complete pairs are available.
.safe_cor <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 3 || sd(x[ok]) == 0 || sd(y[ok]) == 0) return(NA_real_)
  cor(x[ok], y[ok])
}

# Distribution summary of a fit-statistic column: mean, SD, skewness, and
# (excess) kurtosis (the summary block of Andrich & Marais 2019, app. C).
.dist_stats <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 3) return(list(mean = NA_real_, sd = NA_real_,
                                 skewness = NA_real_, kurtosis = NA_real_))
  m <- mean(x); s <- sd(x); d <- x - m
  if (!is.finite(s) || s == 0)
    return(list(mean = m, sd = s, skewness = NA_real_, kurtosis = NA_real_))
  list(mean = m, sd = s,
       skewness = mean(d^3) / (mean(d^2))^1.5,
       kurtosis = mean(d^4) / (mean(d^2))^2 - 3)
}

#' Class-interval detail for one item's chi-square test of fit
#'
#' The per-class-interval breakdown behind an item's item-trait chi-square,
#' as dissected in Andrich and Marais (2019, ch. 13): for every class
#' interval the size, the maximum and mean
#' person location, the standardised residual between observed and expected
#' interval means, its squared chi-square component, the observed and
#' expected means (OM, EV), the sample-size-free effect size
#' ES = (OM - EV)/sqrt(mean V), and per response category the observed
#' proportion (OBS.P), the mean model probability (EST.P), and the observed
#' conditional threshold proportion (OBS.T), the proportion scoring k among
#' those scoring k - 1 or k.
#'
#' @param fit A fitted object from \code{\link{rasch}}.
#' @param item Item name or index.
#' @return A list with \code{item}, \code{location}, the \code{intervals}
#'   data frame, the \code{categories} data frame, the whole-sample observed
#'   mean \code{ave}, and the item's total \code{chisq}, \code{df}, and
#'   \code{p}. Intervals with fewer than 2 responders are shown but carry no
#'   chi-square contribution (\code{used = FALSE}), matching the item-trait
#'   computation. The same applies when the model variance for an interval
#'   is unavailable or zero. The probability is \code{NA} when person IDs
#'   repeat because the asymptotic reference counts response rows rather than
#'   independent persons; the interval summaries and chi-square remain
#'   descriptive.
#' @seealso \code{\link{fit_bootstrap}}, which refers the item's total, and
#'   every other item fit statistic, to a bootstrap null rather than to its
#'   asymptotic distribution.
#' @examples
#' set.seed(1)
#' d <- seq(-1.5, 1.5, length.out = 6)
#' X <- matrix(rbinom(400 * 6, 1, plogis(outer(rnorm(400), d, "-"))), 400, 6)
#' colnames(X) <- paste0("I", 1:6)
#' chisq_detail(rasch(X), "I3")$intervals
#' @export
chisq_detail <- function(fit, item) {
  if (!inherits(fit, "rasch") || inherits(fit, "rasch_btl"))
    stop("chisq_detail needs a response-data Rasch fit", call. = FALSE)
  if (!isTRUE(fit$est$converged))
    stop("the fitted calibration did not converge; item-trait detail is unavailable",
         call. = FALSE)
  if (!.efrm_link_converged(fit))
    stop("the fitted set-unit link did not converge; item-trait detail is unavailable",
         call. = FALSE)
  if (length(item) != 1L) stop("`item` must name exactly one item")
  i <- .item_idx(fit, item)
  # per-item interval allocation when the fit carries one (missing data)
  ci <- if (!is.null(fit$ci_item)) fit$ci_item[[i]] else fit$person$class_interval
  if (all(is.na(ci)))
    .refuse("item ", fit$items$item[i], " has no persons in any class ",
            "interval (only extreme or missing responders); the ",
            "class-interval detail is unavailable")
  th <- fit$person$theta
  x <- fit$X[, i]; E <- fit$moments$E[, i]; V <- fit$moments$V[, i]
  mi <- length(fit$tau_list[[i]])
  G <- max(ci, na.rm = TRUE)
  iv <- data.frame(interval = seq_len(G), n = 0L, theta_max = NA_real_,
                   theta_mean = NA_real_, obs_mean = NA_real_,
                   exp_value = NA_real_, residual = NA_real_,
                   chisq = NA_real_, es = NA_real_, used = FALSE)
  cats <- expand.grid(interval = seq_len(G), category = 0:mi)
  cats <- cats[order(cats$interval, cats$category), ]
  cats$obs_p <- cats$est_p <- cats$obs_t <- NA_real_
  disc_i <- if (is.null(fit$disc)) 1 else fit$disc[i]
  for (g in seq_len(G)) {
    sel <- which(!is.na(ci) & ci == g & !is.na(x))
    iv$n[g] <- length(sel)
    if (!length(sel)) next
    iv$theta_max[g] <- max(th[sel]); iv$theta_mean[g] <- mean(th[sel])
    OM <- mean(x[sel]); EV <- mean(E[sel]); Vbar <- mean(V[sel])
    iv$obs_mean[g] <- OM; iv$exp_value[g] <- EV
    valid_moment <- is.finite(OM) && is.finite(EV) &&
      is.finite(Vbar) && Vbar > 0
    if (valid_moment)
      iv$es[g] <- (OM - EV) / sqrt(Vbar)
    if (length(sel) >= 2 && valid_moment) {
      iv$residual[g] <- sqrt(length(sel)) * (OM - EV) / sqrt(Vbar)
      iv$chisq[g] <- iv$residual[g]^2
      iv$used[g] <- TRUE
    }
    P <- vapply(th[sel], function(b)
      item_moments(b, fit$tau_list[[i]], disc = disc_i)$P, numeric(mi + 1))
    est_p <- rowMeans(P)
    obs_n <- as.integer(table(factor(x[sel], levels = 0:mi)))
    rows <- cats$interval == g
    cats$obs_p[rows] <- obs_n / length(sel)
    cats$est_p[rows] <- est_p
    for (k in seq_len(mi)) {
      pair <- obs_n[k] + obs_n[k + 1]
      cats$obs_t[rows][k + 1] <- if (pair > 0) obs_n[k + 1] / pair else NA_real_
    }
  }
  it_row <- fit$item_trait[i, ]
  item_p <- if (.has_repeated_residual_units(fit)) NA_real_ else it_row$p
  list(item = fit$items$item[i], location = fit$items$location[i],
       intervals = iv, categories = cats, ave = mean(x, na.rm = TRUE),
       chisq = it_row$chisq, df = it_row$df, p = item_p)
}

# Person separation index (separation reliability; Andrich 1982), with the
# separation ratio and the number of distinct strata (4G + 1)/3 (Wright &
# Masters 1982): the count of statistically distinguishable performance
# levels the instrument supports.
.psi <- function(theta, se, keep = TRUE) {
  ok <- is.finite(theta) & is.finite(se) & se >= 0 & keep %in% TRUE
  if (sum(ok) < 3) return(list(PSI = NA_real_, separation = NA_real_,
                               strata = NA_real_, var_theta = NA_real_,
                               mean_error_var = NA_real_, n = sum(ok)))
  vt <- var(theta[ok]); mse <- mean(se[ok]^2)
  if (!is.finite(vt) || vt <= 0 || !is.finite(mse))
    return(list(PSI = NA_real_, separation = NA_real_, strata = NA_real_,
                var_theta = vt, mean_error_var = mse, n = sum(ok)))
  psi <- max((vt - mse) / vt, 0)
  sep <- if (psi < 1) sqrt(psi / (1 - psi)) else Inf
  strata <- if (is.finite(sep)) (4 * sep + 1) / 3 else Inf
  list(PSI = psi, separation = sep, strata = strata, var_theta = vt,
       mean_error_var = mse, n = sum(ok))
}

# Cronbach's alpha (Cronbach 1951) on complete cases, reported alongside
# the PSI. Alpha has no missing-data form, so the applicable flag carries
# that caveat (the complete-case value is still reported, with its n).
.alpha <- function(X) {
  Xc <- X[stats::complete.cases(X), , drop = FALSE]
  applicable <- nrow(Xc) == nrow(X)
  if (nrow(Xc) < 3 || ncol(Xc) < 2) return(list(
    alpha = NA_real_, n = nrow(Xc), applicable = applicable,
    design_applicable = TRUE))
  L <- ncol(Xc); vi <- apply(Xc, 2, var); vt <- var(rowSums(Xc))
  if (!is.finite(vt) || vt <= 0)          # constant total score: undefined
    return(list(alpha = NA_real_, n = nrow(Xc), applicable = applicable,
                design_applicable = TRUE))
  list(alpha = L / (L - 1) * (1 - sum(vi) / vt), n = nrow(Xc),
       applicable = applicable, design_applicable = TRUE)
}

# Whether a structural fit reduces to one observable response cell per item.
# Current fits record this in alpha$design_applicable. The map fallback keeps
# projects saved before that field was introduced safe: an absent flag is not
# permission unless the fitted response cells can be matched one-to-one to
# items.
.classical_design_applicable <- function(fit) {
  structural <- inherits(fit, c("rasch_mfrm", "rasch_efrm"))
  if (!structural) return(TRUE)
  vm <- fit$virtual_map
  if (is.null(vm) || !all(c("vkey", "item") %in% names(vm)) ||
      is.null(colnames(fit$X))) return(FALSE)
  item_names <- vm$item[match(colnames(fit$X), vm$vkey)]
  map_ok <- !anyNA(item_names) && !anyDuplicated(item_names)
  recorded <- fit$alpha$design_applicable
  if (length(recorded)) isTRUE(recorded) && map_ok else map_ok
}

# Qualitative description of person separation, driven by the PSI. This is not
# the statistical power of a test of fit, which also depends on sample size,
# the departure being tested, targeting, category use and the test statistic.
.separation_quality <- function(psi) {
  if (is.na(psi)) "unknown"
  else if (psi >= 0.9) "excellent"
  else if (psi >= 0.8) "good"
  else if (psi >= 0.7) "reasonable"
  else if (psi >= 0.5) "low"
  else "too low"
}

# Compatibility for fitted objects and downstream code created before the
# public label was corrected from power of fit to separation quality.
.fit_power <- .separation_quality

# Targeting summary: how well item thresholds cover the person distribution.
.targeting <- function(person, thresholds) {
  ok <- is.finite(person$theta)
  th <- person$theta[ok]
  tau_ok <- is.finite(thresholds$tau)
  tau <- thresholds$tau[tau_ok]
  item_location <- if (length(tau) && "item" %in% names(thresholds))
    unname(tapply(tau, thresholds$item[tau_ok], mean)) else numeric(0)
  tr <- if (length(tau)) range(tau) else c(NA_real_, NA_real_)
  noext <- ok & person$extreme %in% FALSE
  list(person_mean = if (length(th)) mean(th) else NA_real_,
       person_sd = if (length(th) > 1L) sd(th) else NA_real_,
       person_mean_noext = if (any(noext))
         mean(person$theta[noext]) else NA_real_,
       item_mean = if (length(item_location)) mean(item_location) else NA_real_,
       threshold_range = tr,
       prop_below = if (length(th) && length(tau)) mean(th < tr[1L]) else NA_real_,
       prop_above = if (length(th) && length(tau)) mean(th > tr[2L]) else NA_real_)
}

# One list of display names, joined so the joined string names exactly that
# list. A name may itself contain the separator, so a name that would
# otherwise read as several is quoted and a quote within a name is escaped:
# "a+b"+c is two sets, a+b+c is three. Labels are never parsed back into
# membership; the quoting is what stops two different designs reading alike.
.label_names <- function(x, sep = "+") {
  x <- as.character(x)
  q <- grepl(sep, x, fixed = TRUE) | grepl("\"", x, fixed = TRUE)
  x[q] <- encodeString(x[q], quote = "\"")
  paste(x, collapse = sep)
}

# Display label of one design block, built from the block's own columns, so a
# block restricted to an item selection never names items outside it.
.design_label <- function(fit, cols) {
  vm <- fit$virtual_map
  if (inherits(fit, "rasch_efrm")) {
    g <- vm$group[cols[1L]]
    gcols <- which(vm$group == g)
    sets_of_col <- vm$set[gcols]
    psets <- sort(unique(sets_of_col[gcols %in% cols]))
    # a set is named as a whole only when the block holds every column the
    # group has in it; otherwise the items themselves are named. The set
    # clause is unconditional: an administration is which items of which set
    # a group took, so a one-set frame still names its set, and the group
    # stays because the label alone identifies the design in
    # test_information() and in every curve legend.
    partial <- !all(gcols[sets_of_col %in% psets] %in% cols)
    paste0("group=", g, ", sets=", .label_names(psets),
      if (partial) paste0(", items=", .label_names(vm$item[cols])) else "")
  } else if (inherits(fit, "rasch_mfrm")) {
    fs <- fit$facet_spec
    cell <- .factor_keys(vm[, fs, drop = FALSE])
    cells <- unique(cell)
    active <- cells[cells %in% cell[cols]]
    parts <- vapply(active, function(k) {
      ii <- which(cell == k)
      jj <- ii[ii %in% cols]
      paste0(paste(paste0(fs, "=", unlist(vm[ii[1L], fs, drop = FALSE])),
                   collapse = ", "),
             if (length(jj) < length(ii))
               paste0(" [items=", .label_names(vm$item[jj]), "]")
             else "")
    }, "")
    paste(parts, collapse = " + ")
  } else "test"
}

# Administrable virtual-item blocks of a fit: one per design a person
# could actually take. Ordinary fits: the whole test. EFRM: one block per
# person group AND exact observed item pattern in that group. MFRM: one
# block per observed item-by-facet pattern for a person. Every pattern is
# taken at face value, including where item nonresponse leaves nearly
# every person a pattern of their own: an item a person left unanswered
# carries no information about where that person is, so the curve that
# describes them is the one over the items they answered. Merging their
# pattern into a fuller one -- let alone into a union of patterns nobody
# took -- would claim information no response supports and understate
# their SEM. More designs than a legend can carry is therefore a matter
# for the curve plots, not a reason to report information as though the
# unanswered items had been answered. Shared by test_information() and the
# test-level curve plots so they cannot disagree.
.design_blocks <- function(fit) {
  L <- length(fit$tau_list)
  blocks <- list(test = seq_len(L))
  if (inherits(fit, "rasch_efrm")) {
    # a group's virtual block can span item sets that were only PARTIALLY
    # administered within the group (a linking design: most persons take
    # one set, a linking subsample takes several). Summing all the group's
    # sets would describe a form nobody in the majority sub-population
    # ever took, understating their SEM -- so split each group by the
    # distinct item-administration patterns actually observed
    vm <- fit$virtual_map
    blocks <- list(); labels <- character(0)
    for (g in unique(vm$group)) {
      gcols <- which(vm$group == g)
      grows <- rowSums(!is.na(fit$X[, gcols, drop = FALSE])) > 0
      if (!any(grows)) next
      # Keep the exact observed items, including partial sets. Set membership
      # is used only for labels, not to enlarge a person's administration.
      answered <- !is.na(fit$X[grows, gcols, drop = FALSE])
      pat <- .factor_keys(as.data.frame(answered, check.names = FALSE))
      for (p in unique(pat)) {
        cols <- gcols[answered[match(p, pat), ]]
        blocks[[length(blocks) + 1L]] <- cols
        labels <- c(labels, .design_label(fit, cols))
      }
    }
    # Readable labels are for display only. If literal group or set names make
    # two labels look the same, retain both designs and mark them distinctly.
    lab <- labels
    if (anyDuplicated(lab)) {
      dup <- duplicated(lab) | duplicated(lab, fromLast = TRUE)
      lab[dup] <- paste0(lab[dup], " [design ", seq_along(lab)[dup], "]")
    }
    names(blocks) <- make.unique(lab)
  } else if (inherits(fit, "rasch_mfrm")) {
    observed <- !is.na(fit$X)
    observed <- observed[rowSums(observed) > 0L, , drop = FALSE]
    pat <- .factor_keys(as.data.frame(observed, check.names = FALSE))
    blocks <- lapply(unique(pat), function(p)
      which(observed[match(p, pat), ]))
    labs <- vapply(blocks, function(cols) .design_label(fit, cols), "")
    if (anyDuplicated(labs)) {
      dup <- duplicated(labs) | duplicated(labs, fromLast = TRUE)
      labs[dup] <- paste0(labs[dup], " [design ", seq_along(labs)[dup], "]")
    }
    names(blocks) <- make.unique(labs)
  }
  blocks
}

#' Test information function
#'
#' Fisher information over a grid of person locations, with the corresponding
#' standard error of measurement. Ordinary Rasch fits return one whole-test
#' curve. EFRM fits return one curve per person group and per item
#' administration pattern actually observed within that group (in a linking
#' design, persons who took only the core set get a core-only curve, and the
#' linking subsample gets the pooled one). MFRM fits return one curve per
#' observed item-by-facet pattern for a person, so ratings that jointly
#' inform the same person measure are added and mutually exclusive designs
#' remain separate. Partly answered sets or facet conditions contribute only
#' their observed items; a missing response is not treated as an administered
#' item when defining these patterns. Where item nonresponse leaves nearly
#' every person a pattern of their own, that is what these fits return:
#' an unanswered item carries no information about the person who left it,
#' so no pattern is merged into a fuller one and no curve of theirs is
#' drawn over a design nobody was administered.
#'
#' @details
#' For an administrable block \eqn{\mathcal A}, the information and standard
#' error of measurement are
#' \deqn{I(\theta)=\sum_{i\in\mathcal A}d_i^2
#' \operatorname{Var}(X_i\mid\theta),\qquad
#' \operatorname{SEM}(\theta)=I(\theta)^{-1/2},}
#' where \eqn{d_i} is the frame unit or discrimination multiplier. For an
#' ordinary Rasch fit, \eqn{d_i=1}.
#' Information is returned only for a converged calibration. Comparative
#' Judgement designs use \code{\link{btl_information}} instead.
#'
#' @param fit A fitted object from \code{\link{rasch}}.
#' @param grid Logit grid over which to evaluate the information.
#' @param items Optional item selection: item names or indices. Every design
#'   block is restricted to the named items, so a restricted person-item map
#'   can carry the information of its own selection. The \code{design} labels
#'   are those of the restricted blocks, and blocks that differ only outside
#'   the selection are returned once.
#' @return A data frame with \code{theta}, \code{info}, and \code{sem}. For
#'   EFRM and MFRM fits it also contains a \code{design} column identifying
#'   the administrable frame or facet design.
#' @references
#' Andrich, D. and Marais, I. (2019). A Course in Rasch Measurement Theory:
#' Measuring in the Educational, Social and Health Sciences. Springer.
#' @seealso \code{\link{targeting_table}} and \code{\link{plot_tif}}.
#' @examples
#' set.seed(1)
#' d <- seq(-1.5, 1.5, length.out = 6)
#' X <- matrix(rbinom(300 * 6, 1, plogis(outer(rnorm(300), d, "-"))), 300, 6)
#' colnames(X) <- paste0("I", 1:6)
#' head(test_information(rasch(X)))
#' @export
test_information <- function(fit, grid = NULL, items = NULL) {
  if (!inherits(fit, "rasch")) stop("test_information needs a rasch fit")
  if (inherits(fit, "rasch_btl"))
    stop("test_information is for response-data models; use btl_information() for a Comparative Judgement fit")
  converged <- if (!is.null(fit$est$converged)) fit$est$converged else
    fit$converged
  if (!isTRUE(converged))
    stop("the fitted calibration did not converge; test information is unavailable")
  if (!.efrm_link_converged(fit))
    stop("the fitted set-unit link did not converge; test information is unavailable")
  if (is.null(grid)) {
    grid <- .default_model_grid(fit, by = 0.1)
  } else if (!is.numeric(grid) || is.complex(grid) || !length(grid) ||
             !is.null(dim(grid)) || !is.null(oldClass(grid)) ||
             any(!is.finite(grid))) {
    stop("`grid` must contain one or more plain finite numeric locations",
         call. = FALSE)
  }
  L <- length(fit$tau_list)
  disc <- if (is.null(fit$disc)) rep(1, L) else fit$disc
  # an item subset restricts every design block to the named items, so a
  # restricted person-item map can carry the information of its own
  # selection rather than the whole instrument's
  if (!is.null(items) && (!is.atomic(items) || !is.null(dim(items))))
    stop("`items` must be a vector of item names or indices")
  keep <- if (is.null(items)) seq_len(L)
          else if (is.numeric(items) || is.complex(items)) {
            if (is.complex(items) || any(!is.finite(items)) ||
                any(items != floor(items)))
              stop("`items` indices must be whole numbers")
            if (any(items < 1L) || any(items > L))
              stop("`items` indices must lie in 1..", L)
            ki <- as.integer(items)
            if (anyDuplicated(ki))
              stop("item indices must not name the same item more than once")
            ki
          } else {
            nm <- as.character(items)
            if (anyNA(nm) || any(!nzchar(trimws(nm))))
              stop("`items` must contain non-missing, non-empty item names")
            # Resolve literal names before treating surrounding whitespace as
            # selector syntax, so deliberately spaced columns remain usable.
            direct <- match(nm, fit$items$item)
            canonical <- nm
            canonical[is.na(direct)] <-
              .role_text_values(nm[is.na(direct)])
            if (anyDuplicated(canonical))
              stop("item(s) named more than once: ",
                   paste(unique(canonical[duplicated(canonical)]),
                         collapse = ", "))
            nm <- canonical
            ki <- match(nm, fit$items$item)
            # an EFRM calibrates virtual item-by-group cells; an underlying
            # item name selects every cell it appears in
            if (anyNA(ki) && !is.null(fit$virtual_map)) {
              vm <- fit$virtual_map
              ki <- lapply(nm, function(x) {
                j <- match(x, fit$items$item)
                if (!is.na(j)) return(j)
                which(as.character(vm$item) == x)
              })
              bad <- vapply(ki, length, 0L) == 0L
              if (any(bad))
                stop("item(s) not in the fit: ",
                     paste(nm[bad], collapse = ", "))
              ki <- unique(unlist(ki))
            } else if (anyNA(ki)) {
              stop("item(s) not in the fit: ",
                   paste(nm[is.na(ki)], collapse = ", "))
            }
            ki
          }
  blocks <- .design_blocks(fit)
  blocks <- lapply(blocks, function(ii) intersect(ii, keep))
  blocks <- blocks[vapply(blocks, length, 0L) > 0L]
  if (!length(blocks))
    stop("the item selection leaves no items in any design block")
  if (!is.null(items)) {
    # the curve is the selection's, so its label must be too: blocks that
    # differ only outside the selection are now one design, and a retained
    # label must name only the items that produced the numbers
    blocks <- blocks[!duplicated(vapply(blocks, paste, "", collapse = "+"))]
    names(blocks) <- make.unique(vapply(blocks, function(ii)
      .design_label(fit, ii), ""))
  }
  ans <- lapply(seq_along(blocks), function(j) {
    ii <- blocks[[j]]
    info <- vapply(grid, function(th)
      sum(vapply(ii, function(i)
        disc[i]^2 * item_moments(th, fit$tau_list[[i]],
                                 disc = disc[i])$V, 0)), 0)
    out <- data.frame(theta = grid, info = info, sem = 1 / sqrt(info))
    if (length(blocks) > 1L || names(blocks)[j] != "test")
      out$design <- names(blocks)[j]
    out
  })
  out <- do.call(rbind, ans)
  rownames(out) <- NULL
  out
}

# ---------------------------------------------------------------------------
# Estimated item discrimination
# ===========================================================================
# The slope that maximises an item's own likelihood with the person measures
# and the item's thresholds held at the values the Rasch model gave them.
# One free parameter per item, fitted one item at a time, so it is a
# description of how steeply an item sorts the people the model has already
# located -- not a two-parameter estimate, which would relocate everything at
# once. Reported for polytomous items as well: item_moments() carries the
# discrimination through the partial credit structure, so the same slope
# multiplies every threshold of the item.
#
# Efficiency: the likelihood only needs the DISTINCT person measures, of which
# there are at most one per raw score under complete data, so each evaluation
# costs far less than a pass over the sample.
# ---------------------------------------------------------------------------
.item_discrim <- function(theta, X, tau_list, extreme, bounds = c(0.05, 5)) {
  ok <- !extreme & is.finite(theta)
  if (sum(ok) < 3L) return(rep(NA_real_, length(tau_list)))
  th <- theta[ok]
  Xo <- X[ok, , drop = FALSE]
  ut <- sort(unique(th))
  idx <- match(th, ut)
  vapply(seq_along(tau_list), function(i) {
    y <- Xo[, i]
    g <- !is.na(y)
    # an item with no variation among the non-extreme persons has no slope
    if (sum(g) < 3L || length(unique(y[g])) < 2L) return(NA_real_)
    tau <- tau_list[[i]]
    m <- length(tau)
    cnt <- table(factor(idx[g], levels = seq_along(ut)),
                 factor(y[g], levels = 0:m))
    keep <- rowSums(cnt) > 0
    cnt <- matrix(cnt[keep, ], nrow = sum(keep))
    uth <- ut[keep]
    nll <- function(a) {
      lp <- vapply(uth, function(t)
        log(pmax(item_moments(t, tau, disc = a)$P, 1e-300)), numeric(m + 1L))
      -sum(cnt * t(lp))
    }
    o <- tryCatch(stats::optimize(nll, bounds), error = function(e) NULL)
    if (is.null(o)) NA_real_ else o$minimum
  }, 0)
}
