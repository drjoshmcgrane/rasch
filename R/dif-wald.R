# rasch :: conditional Wald test of DIF on the resolved calibration
# ===========================================================================
# Each item is split by a person factor and the split copies are calibrated
# together with the unsplit items, which anchor the groups on one scale. A
# Wald test of the copies' locations against one another is the conditional
# analogue of Andersen's likelihood-ratio test for that item: no class
# intervals are formed, the raw score conditions the comparison exactly, and
# the test keeps its power after earlier splits, where a residual analysis
# by class interval loses it.
# ===========================================================================

#' Conditional Wald test of DIF on the resolved calibration
#'
#' Tests each item for uniform DIF by splitting it by a person factor,
#' recalibrating with the unsplit items as the anchor, and comparing the
#' locations of the split copies with a Wald test. The comparison conditions
#' on the raw score through the conditional likelihood, so no class
#' intervals are formed and no person location is estimated; it is the
#' item-level analogue of Andersen's (1973) likelihood-ratio test.
#'
#' A split copy of an item that persons in one level of the factor only
#' answered has no group contrast and is left out with a note. Levels with
#' fewer than \code{min_n} distinct responders to an item are dropped from
#' that item's comparison. Items whose resolution is unavailable (an
#' externally anchored item, a group that does not observe every score
#' category, a refit that fails, a location that rests on a near-empty
#' category) keep their row with \code{NA} statistics and a note.
#'
#' Every other item is treated as invariant while one item is tested. Under
#' pervasive DIF that assumption fails and the anchor carries artificial DIF
#' (Andrich and Hagquist 2012), so the test is best used iteratively, as
#' \code{\link{resolve_dif}} does with \code{criterion = "wald"}, or with a
#' set of anchors chosen on other grounds.
#'
#' @param fit A fitted object from \code{\link{rasch}}, including one that
#'   already carries splits.
#' @param factors Person factors to test, as in \code{\link{dif_anova}};
#'   defaults to every nominated factor. Each factor is tested on its own.
#' @param items Items to test; defaults to every item in the fit.
#' @param p_adjust Multiplicity adjustment over every item-by-factor test,
#'   as in \code{\link[stats]{p.adjust}}.
#' @param alpha Significance level for the adjusted probabilities.
#' @param min_n Levels with fewer distinct responders to an item are dropped
#'   from that item's comparison.
#' @return A list of class \code{"rasch_dif_wald"}. \code{summary} has one
#'   row per item and factor: \code{n_levels} compared, \code{shift} (for
#'   two levels the location of the second level minus the first, positive
#'   when the item is harder for the second level; for more levels the range
#'   of the locations), its \code{se} (two levels only), the \code{wald}
#'   statistic, its hypothesis \code{df} (levels minus one), the reference
#'   \code{ref_df} (infinite for independent response rows, the cluster
#'   count minus one for a supported repeated-person calibration),
#'   \code{p}, \code{p_adj} and \code{significant}. \code{levels} gives each
#'   compared level's resolved \code{location}, \code{se} and \code{n}.
#'   \code{notes} records what was left out or withheld and why.
#' @references
#' Andersen, E. B. (1973). A goodness of fit test for the Rasch model.
#' Psychometrika, 38(1), 123--140.
#'
#' Glas, C. A. W. and Verhelst, N. D. (1995). Testing the Rasch model. In
#' G. H. Fischer and I. W. Molenaar (eds), Rasch Models: Foundations, Recent
#' Developments, and Applications (pp. 69--95). Springer.
#'
#' Kopf, J., Zeileis, A. and Strobl, C. (2015). Anchor selection strategies
#' for DIF analysis: review, assessment, and new approaches. Educational and
#' Psychological Measurement, 75(1), 22--56.
#'
#' Andrich, D. and Hagquist, C. (2012). Real and artificial differential
#' item functioning. Journal of Educational and Behavioral Statistics,
#' 37(3), 387--416.
#' @seealso \code{\link{dif_anova}} for the residual analysis by class
#'   interval, \code{\link{dif_size}} for the pairwise magnitudes of one
#'   item and \code{\link{resolve_dif}} for iterative resolution.
#' @examples
#' set.seed(1); n <- 600
#' d <- seq(-2, 2, length.out = 8); g <- rep(c("a", "b"), each = n / 2)
#' sh <- matrix(0, n, 8); sh[g == "b", 3] <- 0.8
#' X <- matrix(rbinom(n * 8, 1, plogis(outer(rnorm(n), d, "-") - sh)), n, 8)
#' colnames(X) <- paste0("I", 1:8)
#' fit <- rasch(data.frame(X, grp = g), factors = "grp")
#' dif_wald(fit)
#' @export
dif_wald <- function(fit, factors = NULL, items = NULL, p_adjust = "holm",
                     alpha = 0.05, min_n = 20L) {
  .check_dif_args(alpha, p_adjust, min_n = min_n)
  if (!inherits(fit, "rasch") || inherits(fit, c("rasch_mfrm", "rasch_efrm")))
    stop("dif_wald needs an ordinary rasch fit with person factors")
  if (!isTRUE(fit$est$converged))
    stop("the fitted calibration did not converge; the Wald test is unavailable")
  fac <- .dif_factors(fit, factors)
  .check_dif_factor_levels(fac)
  all_items <- fit$items$item
  if (is.null(items)) items <- all_items
  if (!is.character(items) || !length(items) || anyNA(items))
    stop("`items` must name at least one item")
  bad <- setdiff(items, all_items)
  if (length(bad)) stop("item(s) not in the fit: ", paste(bad, collapse = ", "))
  items <- unique(items)
  notes <- character(0)
  rows <- list(); lev_rows <- list()
  for (f in names(fac)) {
    grp <- factor(fac[[f]])
    absent <- character(0)
    for (it in items) {
      i <- match(it, all_items)
      n_lev <- .dif_cell_n(grp, !is.na(fit$X[, i]), fit$person$id)
      # a copy answered in one level only has no contrast to test
      if (sum(n_lev > 0L) < 2L) { absent <- c(absent, it); next }
      rs <- .dif_resolve(fit, it, grp, min_n)
      notes <- c(notes, rs$notes)
      k <- length(rs$levs)
      usable <- k >= 2L && all(is.finite(rs$loc)) && !any(rs$weak) &&
        !identical(rs$score_compatible, FALSE) &&
        .covariance_supports_wald(rs$vloc, k)
      ref_df <- rs$df
      if (usable && (is.na(ref_df) || ref_df <= 0)) {
        usable <- FALSE
        notes <- c(notes, paste0(
          it, " [", f, "]: the resolved calibration does not carry enough ",
          "independent-person support for covariance-based inference; the ",
          "Wald test is withheld"))
      }
      shift <- se <- W <- p <- NA_real_
      q <- k - 1L
      if (usable) {
        C <- cbind(-1, diag(q))
        d <- drop(C %*% rs$loc)
        V <- C %*% rs$vloc %*% t(C)
        V <- (V + t(V)) / 2
        if (qr(V)$rank == q) {
          W <- drop(t(d) %*% solve(V, d))
          p <- stats::pf(W / q, q, ref_df, lower.tail = FALSE)
        } else notes <- c(notes, paste0(
          it, " [", f, "]: the contrast covariance is singular; the Wald ",
          "test is withheld"))
        shift <- if (k == 2L) d else diff(range(rs$loc))
        if (k == 2L) se <- sqrt(max(V[1L, 1L], 0))
      } else if (k >= 2L && all(is.finite(rs$loc)))
        shift <- if (k == 2L) diff(rs$loc) else diff(range(rs$loc))
      rows[[length(rows) + 1L]] <- data.frame(
        item = it, factor = f, n_levels = k, shift = shift, se = se,
        wald = W, df = as.numeric(q), ref_df = ref_df, p = p,
        stringsAsFactors = FALSE)
      lev_se <- if (usable) sqrt(pmax(diag(rs$vloc), 0)) else
        rep(NA_real_, k)
      lev_rows[[length(lev_rows) + 1L]] <- data.frame(
        item = it, factor = f, level = rs$levs, location = rs$loc,
        se = lev_se, n = as.integer(n_lev[rs$levs]),
        stringsAsFactors = FALSE)
    }
    if (length(absent)) notes <- c(notes, sprintf(
      "%s: not tested because persons in one level only answered: %s",
      f, paste(absent, collapse = ", ")))
  }
  summary <- if (length(rows)) do.call(rbind, rows) else data.frame(
    item = character(), factor = character(), n_levels = integer(),
    shift = numeric(), se = numeric(), wald = numeric(), df = numeric(),
    ref_df = numeric(), p = numeric(), stringsAsFactors = FALSE)
  summary$p_adj <- .p_adjust_family(summary$p, method = p_adjust)
  summary$significant <- ifelse(is.na(summary$p_adj), NA,
                                summary$p_adj < alpha)
  rownames(summary) <- NULL
  levels <- if (length(lev_rows)) do.call(rbind, lev_rows) else data.frame(
    item = character(), factor = character(), level = character(),
    location = numeric(), se = numeric(), n = integer(),
    stringsAsFactors = FALSE)
  rownames(levels) <- NULL
  out <- list(summary = summary, levels = levels, factors = names(fac),
              alpha = alpha, p_adjust = p_adjust, min_n = min_n,
              notes = unique(notes))
  out <- .tag_tables(out)
  class(out) <- "rasch_dif_wald"
  out
}

#' @export
print.rasch_dif_wald <- function(x, ...) {
  cat(sprintf("Conditional Wald test of DIF by %s (resolved locations, logits)\n",
              paste(x$factors, collapse = ", ")))
  s <- x$summary
  if (!nrow(s)) {
    cat("no item could be tested\n")
  } else {
    num <- vapply(s, is.numeric, TRUE)
    s[num] <- lapply(s[num], round, 3)
    s$significant <- ifelse(is.na(s$significant), "",
                            ifelse(s$significant, "*", ""))
    print(s, row.names = FALSE)
    ref <- unique(x$summary$ref_df[!is.na(x$summary$ref_df)])
    reference <- if (!length(ref)) "every test withheld" else
      if (all(is.infinite(ref))) "chi-square reference" else
      "F reference on the cluster degrees of freedom"
    cat(sprintf("p adjusted by %s over %d item-by-factor test(s); %s\n",
                x$p_adjust, nrow(s), reference))
    cat("shift: second level minus first for two levels, range of the",
        "locations otherwise; positive means harder for the second level\n")
  }
  if (length(x$notes)) cat("notes:", paste(x$notes, collapse = "; "), "\n")
  invisible(x)
}
