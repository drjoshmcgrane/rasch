# rmt :: dimensionality and local dependence
# ===========================================================================
# Residual correlations for local dependence, principal components of the
# residuals, and the Smith (2002) paired t-test of unidimensionality on the
# first residual contrast.
# ===========================================================================

#' Residual correlations for local dependence (Yen's Q3)
#'
#' The pairwise correlations of the standardised response residuals are
#' Yen's (1984) Q3 statistics. Under unidimensionality and local
#' independence the off-diagonal values sit near \code{-1/(L-1)}; large
#' positive values flag local dependence between item pairs. Following
#' Christensen, Makransky and Horton (2017), each Q3 is also reported
#' relative to the average off-diagonal value (\code{q3_star}), and a pair
#' is flagged when that excess passes \code{flag}.
#'
#' @param fit A fitted object from \code{\link{rasch}}.
#' @param flag Excess above the average off-diagonal Q3 at which a pair is
#'   flagged as dependent.
#' @return A list with the Q3 \code{matrix}, the \code{average} off-diagonal
#'   value, \code{pairs} (every item pair with \code{q3}, \code{q3_star} and
#'   a \code{flagged} indicator, sorted by \code{q3}), and the subset of
#'   \code{flagged} pairs.
#' @references Yen, W. M. (1984). Effects of local item dependence on the
#'   fit and equating performance of the three-parameter logistic model.
#'   \emph{Applied Psychological Measurement}, 8(2), 125-145.
#'
#'   Christensen, K. B., Makransky, G., & Horton, M. (2017). Critical values
#'   for Yen's Q3: identification of local dependence in the Rasch model
#'   using residual correlations. \emph{Applied Psychological Measurement},
#'   41(3), 178-194.
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
  idx <- which(upper.tri(R), arr.ind = TRUE)
  pairs <- data.frame(item_a = colnames(Z)[idx[, 1]],
                      item_b = colnames(Z)[idx[, 2]],
                      q3 = R[idx], q3_star = R[idx] - avg,
                      flagged = (R[idx] - avg) > flag)
  pairs <- pairs[!is.na(pairs$q3), ]
  pairs <- pairs[order(-pairs$q3), ]
  rownames(pairs) <- NULL
  list(matrix = R, average = avg, pairs = pairs,
       flagged = pairs[pairs$flagged, ])
}

#' Principal components of the residual correlations
#'
#' The first residual contrast (PC1) carries any second dimension; items with
#' opposing loadings define the split used by the unidimensionality t-test.
#' Loadings for the leading components and the eigenvalue table support
#' inspection beyond the first contrast.
#'
#' @param fit A fitted object from \code{\link{rasch}}.
#' @param n_components Number of leading components to return loadings and
#'   eigenvalue rows for (capped at the number of items).
#' @return A list with the residual \code{eigenvalues}, their
#'   \code{prop}ortions, the first-contrast \code{loadings} (sorted), the
#'   \code{loadings_matrix} for the leading components, the
#'   \code{eigen_table} (component, eigenvalue, proportion, cumulative), and
#'   the \code{first_eigen}value.
#' @examples
#' set.seed(1)
#' d <- seq(-2, 2, length.out = 8)
#' X <- matrix(rbinom(500 * 8, 1, plogis(outer(rnorm(500), d, "-"))), 500, 8)
#' colnames(X) <- paste0("I", 1:8)
#' residual_pca(rasch(X))$first_eigen
#' @export
residual_pca <- function(fit, n_components = 10) {
  R <- cor(fit$residuals, use = "pairwise.complete.obs")
  # column pairs with no overlapping persons (for example item-by-group
  # columns from different groups) carry no dependence information
  R[is.na(R)] <- 0; diag(R) <- 1
  ev <- eigen(R, symmetric = TRUE)
  k <- min(n_components, ncol(R))
  loadings <- ev$vectors[, 1] * sqrt(pmax(ev$values[1], 0))
  ld <- data.frame(item = colnames(fit$residuals), pc1_loading = loadings)
  lm <- ev$vectors[, seq_len(k), drop = FALSE] %*%
    diag(sqrt(pmax(ev$values[seq_len(k)], 0)), k)
  colnames(lm) <- paste0("PC", seq_len(k))
  # a pairwise-complete correlation matrix need not be positive
  # semi-definite; proportions are taken over the positive eigenvalue mass
  # so the cumulative share cannot exceed one
  tot <- sum(pmax(ev$values, 0))
  list(eigenvalues = ev$values, prop = pmax(ev$values, 0) / tot,
       loadings = ld[order(-ld$pc1_loading), ],
       loadings_matrix = data.frame(item = colnames(fit$residuals), lm),
       eigen_table = data.frame(component = seq_len(k),
                                eigenvalue = ev$values[seq_len(k)],
                                proportion = pmax(ev$values[seq_len(k)], 0) / tot,
                                cumulative = cumsum(pmax(ev$values[seq_len(k)], 0)) / tot),
       first_eigen = ev$values[1])
}

#' Scree plot of the residual components with parallel analysis
#'
#' Eigenvalues of the residual correlation matrix for the leading components,
#' with a parallel-analysis reference: the mean eigenvalues of residual-sized
#' random normal matrices sharing the data's missingness pattern. Observed
#' eigenvalues above the reference suggest structure beyond chance.
#'
#' @param fit A fitted object from \code{\link{rasch}}.
#' @param n_components Number of leading components to display.
#' @param parallel Draw the parallel-analysis reference line.
#' @param reps Random replicates for the reference.
#' @return Called for its plotting side effect; invisibly the eigen table.
#' @examples
#' set.seed(1)
#' d <- seq(-2, 2, length.out = 8)
#' X <- matrix(rbinom(500 * 8, 1, plogis(outer(rnorm(500), d, "-"))), 500, 8)
#' colnames(X) <- paste0("I", 1:8)
#' plot_scree(rasch(X))
#' @export
plot_scree <- function(fit, n_components = 10, parallel = TRUE, reps = 20) {
  pc <- residual_pca(fit, n_components)
  k <- nrow(pc$eigen_table)
  obs <- pc$eigen_table$eigenvalue
  pa <- NULL
  if (parallel) {
    Z <- fit$residuals
    sim <- replicate(reps, {
      Zr <- matrix(rnorm(length(Z)), nrow(Z))
      Zr[is.na(Z)] <- NA
      Rr <- cor(Zr, use = "pairwise.complete.obs")
      Rr[is.na(Rr)] <- 0; diag(Rr) <- 1
      eigen(Rr, symmetric = TRUE, only.values = TRUE)$values[seq_len(k)]
    })
    pa <- rowMeans(sim)
  }
  ylim <- c(0, max(c(obs, pa, 1)) * 1.12)
  op <- .rr_canvas(c(0.5, k + 0.5), ylim, "Component", "Eigenvalue",
                   grid_x = FALSE)
  on.exit(par(op))
  abline(h = 1, lty = 3, col = .rr$soft)
  if (!is.null(pa)) lines(seq_len(k), pa, lwd = 2, lty = 5, col = .rr$red)
  lines(seq_len(k), obs, lwd = 2.6, col = .rr$blue)
  points(seq_len(k), obs, pch = 21, bg = .rr$blue, col = "white", cex = 1.5)
  axis(1, at = seq_len(k), col = NA, col.ticks = NA)
  .rr_legend("topright",
             if (is.null(pa)) "Observed" else c("Observed", "Parallel analysis"),
             lwd = c(2.6, if (!is.null(pa)) 2), lty = c(1, if (!is.null(pa)) 5),
             col = c(.rr$blue, if (!is.null(pa)) .rr$red))
  invisible(pc$eigen_table)
}

#' Residual-component test of unidimensionality
#'
#' Estimates each person separately on two item subsets and compares the two
#' estimates with a per-person t-test (Smith 2002). By default the subsets
#' are defined by the sign of the first residual-contrast loading; they can
#' also be nominated manually (for example, by content). Under
#' unidimensionality and local independence the two subset estimates are
#' independent given the person location, so
#' \code{t = (theta_A - theta_B) / sqrt(se_A^2 + se_B^2)} is approximately
#' standard normal and about \code{alpha} of the tests should reach
#' significance. Persons with an extreme score on either subset are excluded
#' (their weighted-likelihood estimates are most biased there). The
#' proportion of significant tests is reported with an exact
#' (Clopper-Pearson) binomial confidence interval; a lower bound above
#' \code{alpha} signals multidimensionality.
#'
#' @param fit A fitted object from \code{\link{rasch}}.
#' @param alpha Nominal significance level for the per-person t-tests.
#' @param items_positive,items_negative Optional character vectors naming the
#'   two item subsets; both must be given (disjoint, at least two items
#'   each), otherwise the first residual contrast defines the split.
#' @return A list with the proportion of significant tests, its exact
#'   confidence interval, the sample sizes (\code{n} used,
#'   \code{n_excluded_extreme}), the item split and its source, a
#'   \code{multidimensional} verdict, and \code{paired_t}, the paired t-test
#'   of the two subset means (the group-level comparison, which requires
#'   pairing because both estimates come from the same persons).
#' @examples
#' set.seed(1)
#' d <- seq(-2, 2, length.out = 8)
#' X <- matrix(rbinom(500 * 8, 1, plogis(outer(rnorm(500), d, "-"))), 500, 8)
#' colnames(X) <- paste0("I", 1:8)
#' dimensionality_test(rasch(X))$multidimensional
#' @export
dimensionality_test <- function(fit, alpha = 0.05, items_positive = NULL,
                                items_negative = NULL) {
  X <- fit$X
  manual <- !is.null(items_positive) || !is.null(items_negative)
  pca <- residual_pca(fit)
  if (manual) {
    if (is.null(items_positive) || is.null(items_negative))
      stop("supply both item subsets, or neither")
    pos <- match(items_positive, colnames(X))
    neg <- match(items_negative, colnames(X))
    if (anyNA(pos) || anyNA(neg))
      stop("subset item(s) not in the fit: ",
           paste(c(items_positive[is.na(pos)], items_negative[is.na(neg)]),
                 collapse = ", "))
    if (length(intersect(pos, neg)))
      stop("the two subsets must be disjoint")
    split_source <- "manual"
  } else {
    pos <- which(pca$loadings$pc1_loading[match(colnames(X), pca$loadings$item)] > 0)
    neg <- setdiff(seq_len(ncol(X)), pos)
    split_source <- "first residual contrast"
  }
  if (length(pos) < 2 || length(neg) < 2)
    return(list(note = "need >= 2 items in each subset"))
  est_sub <- function(cols) {
    d <- if (is.null(fit$disc)) rep(1, ncol(X)) else fit$disc
    if (length(unique(d[cols])) == 1L)
      .person_estimates(X[, cols, drop = FALSE], fit$tau_list[cols],
                        disc = d[cols][1])
    else .efrm_person_estimates(X[, cols, drop = FALSE], fit$tau_list[cols],
                                d[cols])
  }
  a <- est_sub(pos); b <- est_sub(neg)
  usable <- !is.na(a$theta) & !is.na(b$theta) & !is.na(a$se) & !is.na(b$se)
  ok <- usable & !a$extreme & !b$extreme
  t <- (a$theta[ok] - b$theta[ok]) / sqrt(a$se[ok]^2 + b$se[ok]^2)
  n <- sum(ok)
  if (n < 10) return(list(note = "fewer than 10 usable persons for the t-test"))
  n_sig <- sum(abs(t) > qnorm(1 - alpha / 2))
  bt <- stats::binom.test(n_sig, n, p = alpha)
  # paired t-test of the two subset means (the group-level comparison: the
  # two estimates come from the same persons, so the means need pairing;
  # Andrich & Marais 2019, ch. 24)
  dd <- a$theta[ok] - b$theta[ok]
  pt <- stats::t.test(dd)
  list(prop_significant = n_sig / n, ci = as.numeric(bt$conf.int), n = n,
       n_excluded_extreme = sum(usable) - n,
       multidimensional = bt$conf.int[1] > alpha,
       split = split_source,
       items_positive = colnames(X)[pos], items_negative = colnames(X)[neg],
       first_eigenvalue = pca$first_eigen,
       paired_t = list(mean_difference = mean(dd),
                       t = unname(pt$statistic), df = unname(pt$parameter),
                       p = pt$p.value))
}

#' Magnitude of multidimensionality from a subtest analysis
#'
#' Estimates how strongly two or more hypothesised subscales measure
#' distinct traits, by Andrich's (2016) comparison of two reliability
#' calculations: one treating all items as independent (which inflates
#' reliability under multidimensionality) and one on the subtest analysis
#' in which each subscale is combined into a single polytomous super-item
#' (which absorbs the unique subscale variance). Under the bifactor
#' formalisation \eqn{\beta_{ns} = \beta_n + c\,\beta'_{ns}} (Marais and
#' Andrich 2008), with \eqn{S} subscales of \eqn{K} items,
#' \deqn{c^2 = S\,(r_1/r_2 - 1) \frac{SK - 1}{S(K - 1)},}
#' the latent correlation between subscales is \eqn{\rho = 1/(1 + c^2)},
#' and \eqn{A = S/(S + c^2)} is the proportion of common (non-unique,
#' non-error) variance. Both the person separation index and coefficient
#' alpha versions are reported (Andrich and Marais 2019, ch. 24).
#'
#' @param fit A fitted object from \code{\link{rasch}}.
#' @param subtests A list of character vectors assigning \emph{every} item
#'   of the fit to one subscale (at least two subscales of two or more
#'   items). Unequal subscale sizes use their mean as \eqn{K}.
#' @return A list of class \code{"rmt_dim_magnitude"}: the comparison
#'   \code{table} (rows PSI and alpha; columns \code{run1}, \code{subtest},
#'   \code{c2}, \code{c}, \code{rho}, \code{A}), the subtest \code{refit},
#'   and the design constants \code{S} and \code{K}.
#' @references Andrich, D. (2016). Components of variance of scales with a
#'   bifactor structure from two calculations of coefficient alpha.
#'   Educational Measurement: Issues and Practice, 35(4), 25-30.
#' @examples
#' set.seed(1); N <- 500
#' common <- rnorm(N); u1 <- rnorm(N); u2 <- rnorm(N)
#' d <- rep(seq(-1, 1, length.out = 5), 2)
#' X <- sapply(1:10, function(i) rbinom(N, 1,
#'   plogis(common + 0.8 * (if (i <= 5) u1 else u2) - d[i])))
#' colnames(X) <- paste0("I", 1:10)
#' fit <- rasch(X)
#' dimensionality_magnitude(fit,
#'   list(paste0("I", 1:5), paste0("I", 6:10)))$table
#' @export
dimensionality_magnitude <- function(fit, subtests) {
  if (!inherits(fit, "rasch")) stop("dimensionality_magnitude needs a rasch fit")
  if (!is.list(subtests) || length(subtests) < 2)
    stop("supply a list of at least two subscales")
  allit <- unlist(subtests)
  if (!setequal(allit, fit$items$item) || anyDuplicated(allit))
    stop("subtests must assign every item of the fit to exactly one subscale")
  S <- length(subtests)
  K <- mean(lengths(subtests))
  g <- S * (K - 1) / (S * K - 1)
  refit <- combine_items(fit, subtests)
  ratio <- function(r1, r2) {
    if (any(!is.finite(c(r1, r2))) || r2 <= 0) return(rep(NA_real_, 4))
    c2 <- max(S * (r1 / r2 - 1) / g, 0)
    c(c2, sqrt(c2), 1 / (1 + c2), S / (S + c2))
  }
  psi_row <- ratio(fit$psi$PSI, refit$psi$PSI)
  alp_row <- ratio(fit$alpha$alpha, refit$alpha$alpha)
  tab <- data.frame(index = c("PSI", "alpha"),
                    run1 = c(fit$psi$PSI, fit$alpha$alpha),
                    subtest = c(refit$psi$PSI, refit$alpha$alpha),
                    c2 = c(psi_row[1], alp_row[1]),
                    c = c(psi_row[2], alp_row[2]),
                    rho = c(psi_row[3], alp_row[3]),
                    A = c(psi_row[4], alp_row[4]))
  out <- list(table = tab, refit = refit, S = S, K = K,
              alpha_applicable = fit$alpha$applicable)
  class(out) <- "rmt_dim_magnitude"
  out
}

#' @export
print.rmt_dim_magnitude <- function(x, ...) {
  cat(sprintf("Magnitude of multidimensionality (Andrich 2016): %d subscales, mean %.1f items\n",
              x$S, x$K))
  tab <- x$table
  num <- vapply(tab, is.numeric, TRUE)
  tab[num] <- lapply(tab[num], round, 3)
  print(tab, row.names = FALSE)
  cat("rho = latent correlation between subscales; A = proportion of common variance\n")
  if (isFALSE(x$alpha_applicable))
    cat("note: alpha computed on complete cases only (missing data present)\n")
  invisible(x)
}
