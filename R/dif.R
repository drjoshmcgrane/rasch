# RaschR :: differential item functioning
# ===========================================================================
# DIF by two-way residual analysis of variance (Hagquist & Andrich 2017):
# for each item
# and each person factor, standardised residuals are analysed by factor and
# trait class interval. The factor main effect indicates uniform DIF; the
# factor-by-interval interaction indicates non-uniform DIF. Any number of
# person factors can be analysed in one call.
# ===========================================================================

#' Differential item functioning by two-way residual ANOVA
#'
#' For each item and each person factor, analyses the standardised residuals
#' by factor group and trait class interval. The group main effect indicates
#' uniform DIF and the group-by-interval interaction indicates non-uniform
#' DIF, with Bonferroni flagging across items within each factor.
#'
#' @param fit A fitted object from \code{\link{rasch}}.
#' @param factors A vector (one factor), a data frame of person factors, or a
#'   character vector naming factor columns already nominated in the fit (via
#'   \code{rasch(..., factors = )}). Defaults to every factor stored in the
#'   fit.
#' @param n_groups Number of trait class intervals; defaults to the value used
#'   in the fit.
#' @return A data frame of uniform and non-uniform DIF statistics per item and
#'   factor, with mean squares, F statistics, probabilities, and Bonferroni
#'   flags.
#' @examples
#' set.seed(1); n <- 600
#' d <- seq(-2, 2, length.out = 8); g <- rep(c("a", "b"), each = n / 2)
#' sh <- matrix(0, n, 8); sh[g == "b", 3] <- 1
#' X <- matrix(rbinom(n * 8, 1, plogis(outer(rnorm(n), d, "-") - sh)), n, 8)
#' colnames(X) <- paste0("I", 1:8)
#' dif_anova(rasch(X), factors = data.frame(group = g))
#' @export
dif_anova <- function(fit, factors = NULL, n_groups = NULL) {
  Z <- fit$residuals; L <- ncol(Z)
  if (is.null(factors)) factors <- fit$factors
  if (is.null(factors)) stop("no person factors supplied or stored in the fit")
  if (is.character(factors) && !is.null(fit$factors) &&
      all(factors %in% names(fit$factors)))
    factors <- fit$factors[, factors, drop = FALSE]
  if (!is.data.frame(factors)) factors <- data.frame(group = factors)
  if (is.null(n_groups)) n_groups <- fit$n_groups

  ci <- fit$person$class_interval
  if (is.null(ci) || !identical(n_groups, fit$n_groups))
    ci <- .class_intervals(fit$person$theta, fit$person$extreme, n_groups)
  ci <- factor(ci)

  res <- list()
  for (fname in names(factors)) {
    grp <- factor(factors[[fname]])
    out <- data.frame(factor = fname, item = colnames(Z),
                      F_uniform = NA_real_, p_uniform = NA_real_,
                      F_nonuniform = NA_real_, p_nonuniform = NA_real_)
    for (i in seq_len(L)) {
      d <- data.frame(z = Z[, i], g = grp, ci = ci)
      d <- d[stats::complete.cases(d), ]
      if (nrow(d) < 10 || length(unique(d$g)) < 2) next
      a <- tryCatch(stats::anova(stats::lm(z ~ g * ci, data = d)),
                    error = function(e) NULL)
      if (is.null(a)) next
      rn <- rownames(a)
      if ("g" %in% rn) {
        out$F_uniform[i] <- a["g", "F value"]; out$p_uniform[i] <- a["g", "Pr(>F)"]
      }
      if ("g:ci" %in% rn) {
        out$F_nonuniform[i] <- a["g:ci", "F value"]
        out$p_nonuniform[i] <- a["g:ci", "Pr(>F)"]
      }
    }
    out$uniform_DIF    <- !is.na(out$p_uniform)    & out$p_uniform    < 0.05 / L
    out$nonuniform_DIF <- !is.na(out$p_nonuniform) & out$p_nonuniform < 0.05 / L
    res[[fname]] <- out
  }
  out <- do.call(rbind, res)
  rownames(out) <- NULL
  out
}
