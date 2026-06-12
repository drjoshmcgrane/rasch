# RaschR :: dimensionality and local dependence
# ===========================================================================
# Residual correlations for local dependence, principal components of the
# residuals, and the Smith (2002) paired t-test of unidimensionality on the
# first residual contrast.
# ===========================================================================

#' Residual correlations for local dependence
#'
#' Under unidimensionality and local independence the off-diagonal residual
#' correlations sit near \code{-1/(L-1)}; large positive values flag local
#' dependence between item pairs (a common convention inspects values more
#' than 0.2 above the average off-diagonal correlation; Christensen,
#' Makransky and Horton 2017).
#'
#' @param fit A fitted object from \code{\link{rasch}}.
#' @param flag Excess-correlation threshold above the average for flagging a
#'   dependent item pair.
#' @return A list with the residual correlation \code{matrix}, the
#'   \code{average} off-diagonal correlation, and a table of \code{flagged}
#'   pairs.
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
  list(eigenvalues = ev$values, prop = ev$values / sum(ev$values),
       loadings = ld[order(-ld$pc1_loading), ],
       loadings_matrix = data.frame(item = colnames(fit$residuals), lm),
       eigen_table = data.frame(component = seq_len(k),
                                eigenvalue = ev$values[seq_len(k)],
                                proportion = ev$values[seq_len(k)] / sum(ev$values),
                                cumulative = cumsum(ev$values[seq_len(k)]) / sum(ev$values)),
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
                   "Scree of the residual components", grid_x = FALSE)
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
#'   \code{n_excluded_extreme}), the item split and its source, and a
#'   \code{multidimensional} verdict.
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
  list(prop_significant = n_sig / n, ci = as.numeric(bt$conf.int), n = n,
       n_excluded_extreme = sum(usable) - n,
       multidimensional = bt$conf.int[1] > alpha,
       split = split_source,
       items_positive = colnames(X)[pos], items_negative = colnames(X)[neg],
       first_eigenvalue = pca$first_eigen)
}
