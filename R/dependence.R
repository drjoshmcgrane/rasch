# rmt :: quantifying response dependence
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
#' Because the resolved items are answered by disjoint persons, the
#' estimates are independent, and the standard error pools the threshold
#' variances: \eqn{\hat\sigma^2_k = (\hat\sigma^2_{(k)(k-1)} +
#' \hat\sigma^2_{(k)(k)})/4}, with \eqn{V[\hat d] = \bar{\hat\sigma^2_k}/m}
#' (eqs. 24.9-24.11).
#'
#' @param fit A fitted object from \code{\link{rasch}}.
#' @param dependent,independent Item names or indices: the item hypothesised
#'   to depend, and the item it depends on. Both must share the same maximum
#'   score (the formalisation requires it).
#' @return A list of class \code{"rmt_dependence"}: the estimate \code{d},
#'   its \code{se}, \code{z} and \code{p} for the hypothesis \eqn{d = 0},
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
  if (!inherits(fit, "rasch")) stop("dependence_magnitude needs a rasch fit")
  j <- .item_idx(fit, dependent); i <- .item_idx(fit, independent)
  if (i == j) stop("dependent and independent must be different items")
  mi <- fit$m[i]; mj <- fit$m[j]
  if (mi != mj)
    stop("the two items must share the same maximum score (here ",
         mj, " and ", mi, "); see Andrich, Humphry and Marais (2012)")
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
    colnames(Xn)[ncol(Xn)] <- res_names[x + 1]
  }
  refit <- rasch(Xn, model = fit$model, id = fit$person$id,
                 factors = fit$factors, n_groups = fit$n_groups)
  if (!all(res_names %in% refit$items$item))
    stop("a resolved item was dropped during re-analysis; too little data")

  thr <- refit$thresholds
  item_of <- refit$items$item[thr$item]
  tab <- data.frame(k = seq_len(mj), delta_lo = NA_real_, delta_hi = NA_real_,
                    d_k = NA_real_, se_k = NA_real_)
  for (k in seq_len(mj)) {
    lo <- thr[item_of == res_names[k]     & thr$k == k, ]  # x_i = k - 1: +d
    hi <- thr[item_of == res_names[k + 1] & thr$k == k, ]  # x_i = k    : -d
    tab$delta_lo[k] <- lo$tau; tab$delta_hi[k] <- hi$tau
    tab$d_k[k] <- (lo$tau - hi$tau) / 2
    tab$se_k[k] <- sqrt((lo$se^2 + hi$se^2) / 4)
  }
  d <- mean(tab$d_k)
  se <- sqrt(mean(tab$se_k^2) / mj)
  z <- d / se
  out <- list(d = d, se = se, z = z, p = 2 * pnorm(-abs(z)),
              thresholds = tab, dependent = nm_j, independent = nm_i,
              refit = refit)
  class(out) <- "rmt_dependence"
  out
}

#' @export
print.rmt_dependence <- function(x, ...) {
  cat(sprintf("Response dependence of %s on %s (Andrich & Kreiner resolution)\n",
              x$dependent, x$independent))
  cat(sprintf("  d = %.3f logits (se %.3f), z = %.2f, p = %s\n",
              x$d, x$se, x$z, .fmt_p(x$p)))
  if (nrow(x$thresholds) > 1) {
    cat("  per threshold:\n")
    print(round(x$thresholds, 3), row.names = FALSE)
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
#' A spread estimate below the bound therefore indicates response
#' dependence among the members (Andrich and Marais 2019, Table 24.1).
#' Typically applied after \code{\link{combine_items}}, whose super-items
#' are exactly such subtests.
#'
#' @param fit A fitted object from \code{\link{rasch}}.
#' @param maxit,tol Passed to the \code{\link{pcml_pc}} refit.
#' @return A data frame with one row per polytomous item: \code{item},
#'   \code{m}, the \code{spread} estimate and its \code{se}, the bound
#'   \code{lub} (available for maximum scores 2 to 8), \code{z} =
#'   (spread - lub)/se, and \code{dependent} = spread below the bound.
#'   Dichotomous items carry no spread and are omitted.
#' @examples
#' set.seed(1); N <- 600
#' d0 <- seq(-1.5, 1.5, length.out = 8)
#' X <- matrix(rbinom(N * 8, 1, plogis(outer(rnorm(N), d0, "-"))), N, 8)
#' X[, 5] <- X[, 4]; X[, 6] <- X[, 4]                 # a dependent triple
#' colnames(X) <- paste0("I", 1:8)
#' fit2 <- combine_items(rasch(X), list(c("I4", "I5", "I6"), c("I1", "I2", "I3")))
#' spread_test(fit2)
#' @export
spread_test <- function(fit, maxit = 60, tol = 1e-8) {
  if (!inherits(fit, "rasch")) stop("spread_test needs a rasch fit")
  poly <- which(fit$m >= 2)
  if (!length(poly)) stop("no polytomous items: no spread parameters to test")
  pc <- pcml_pc(fit$X, maxit = maxit, tol = tol)
  cmp <- pc$components
  out <- data.frame(item = cmp$item[poly], m = fit$m[poly],
                    spread = cmp$spread[poly], se = cmp$spread_se[poly],
                    lub = unname(.spread_lub[as.character(fit$m[poly])]))
  out$z <- (out$spread - out$lub) / out$se
  out$dependent <- !is.na(out$lub) & !is.na(out$spread) & out$spread < out$lub
  rownames(out) <- NULL
  class(out) <- c("rmt_spread", "data.frame")
  out
}

#' @export
print.rmt_spread <- function(x, ...) {
  cat("Spread-parameter screen (Andrich 1985): spread below the binomial bound indicates dependence\n")
  d <- as.data.frame(x)
  names(d)[names(d) == "lub"] <- "bound"
  print(.fmt_df(d), row.names = FALSE)
  invisible(x)
}
