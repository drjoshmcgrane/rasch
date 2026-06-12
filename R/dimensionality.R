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
#'
#' @param fit A fitted object from \code{\link{rasch}}.
#' @return A list with the residual \code{eigenvalues}, their
#'   \code{prop}ortions, the first-contrast \code{loadings}, and the
#'   \code{first_eigen}value.
#' @examples
#' set.seed(1)
#' d <- seq(-2, 2, length.out = 8)
#' X <- matrix(rbinom(500 * 8, 1, plogis(outer(rnorm(500), d, "-"))), 500, 8)
#' colnames(X) <- paste0("I", 1:8)
#' residual_pca(rasch(X))$first_eigen
#' @export
residual_pca <- function(fit) {
  R <- cor(fit$residuals, use = "pairwise.complete.obs")
  # column pairs with no overlapping persons (for example item-by-group
  # columns from different groups) carry no dependence information
  R[is.na(R)] <- 0; diag(R) <- 1
  ev <- eigen(R, symmetric = TRUE)
  loadings <- ev$vectors[, 1] * sqrt(pmax(ev$values[1], 0))
  ld <- data.frame(item = colnames(fit$residuals), pc1_loading = loadings)
  list(eigenvalues = ev$values, prop = ev$values / sum(ev$values),
       loadings = ld[order(-ld$pc1_loading), ], first_eigen = ev$values[1])
}

#' Residual-component test of unidimensionality
#'
#' Splits the items by the sign of their first residual-contrast loading,
#' estimates each person on each subset, and compares the two estimates with a
#' per-person t-test (Smith 2002). A proportion of significant tests whose
#' lower binomial confidence bound exceeds \code{alpha} signals
#' multidimensionality.
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
  est_sub <- function(cols) {
    d <- if (is.null(fit$disc)) rep(1, ncol(X)) else fit$disc
    if (length(unique(d[cols])) == 1L)
      .person_estimates(X[, cols, drop = FALSE], fit$tau_list[cols],
                        disc = d[cols][1])
    else .efrm_person_estimates(X[, cols, drop = FALSE], fit$tau_list[cols],
                                d[cols])
  }
  a <- est_sub(pos); b <- est_sub(neg)
  ok <- !is.na(a$theta) & !is.na(b$theta) & !is.na(a$se) & !is.na(b$se)
  t <- (a$theta[ok] - b$theta[ok]) / sqrt(a$se[ok]^2 + b$se[ok]^2)
  prop <- mean(abs(t) > qnorm(1 - alpha / 2)); n <- sum(ok)
  ci <- prop + c(-1, 1) * 1.96 * sqrt(prop * (1 - prop) / n)
  list(prop_significant = prop, ci = ci, n = n,
       multidimensional = ci[1] > alpha,
       items_positive = colnames(X)[pos], items_negative = colnames(X)[neg],
       first_eigenvalue = pca$first_eigen)
}
