# rmt :: structure amendments (subtests and item splitting)
# ===========================================================================
# Two standard remedies applied by re-analysis. Local dependence: combine the
# dependent items into a single polytomous super-item (the subtest) whose
# score is the sum of its members. Differential item functioning: split the
# offending item into one item per group level, each observed only by that
# group, so every group receives its own location and the invariance
# violation is resolved within the model.
# ===========================================================================

#' Combine items into subtests and re-analyse
#'
#' Forms one polytomous super-item from each nominated group of items (its
#' score is the member sum; missing if any member is missing), keeps all other
#' items as they are, and refits the model with the same settings. The usual
#' treatment for item pairs flagged by \code{\link{residual_correlations}}.
#'
#' @param fit A fitted object from \code{\link{rasch}}.
#' @param groups A list of character vectors, each naming two or more items to
#'   combine; a single vector is also accepted.
#' @param model Model for the re-analysis; defaults to \code{"PCM"}, which is
#'   almost always required because subtests change the maximum scores.
#' @return A new \code{\link{rasch}} fit on the combined structure, with the
#'   combinations recorded in its notes.
#' @examples
#' set.seed(1); Np <- 500; L <- 8
#' d <- seq(-2, 2, length.out = L)
#' X <- matrix(rbinom(Np * L, 1, plogis(outer(rnorm(Np), d, "-"))), Np, L)
#' colnames(X) <- paste0("I", 1:L)
#' X[, 5] <- ifelse(runif(Np) < 0.9, X[, 4], X[, 5])   # dependent pair
#' fit <- rasch(X)
#' fit2 <- combine_items(fit, list(c("I4", "I5")))
#' fit2$items$item
#' @export
combine_items <- function(fit, groups, model = "PCM") {
  if (!inherits(fit, "rasch")) stop("combine_items needs a rasch fit")
  if (inherits(fit, "rasch_mfrm"))
    stop("combine items in the long-format data and refit rasch_mfrm instead")
  if (inherits(fit, "rasch_efrm"))
    stop("combine items in the source data and refit rasch_efrm instead")
  if (is.character(groups)) groups <- list(groups)
  X <- fit$X
  used <- unlist(groups)
  bad <- setdiff(used, colnames(X))
  if (length(bad)) stop("item(s) not in the fit: ", paste(bad, collapse = ", "))
  if (anyDuplicated(used)) stop("an item appears in more than one subtest")
  short <- vapply(groups, length, 1L) < 2
  if (any(short)) stop("each subtest needs at least two items")

  keep <- setdiff(colnames(X), used)
  # a super-item requires all its members: summing whatever happens to be
  # present would let the maximum vary by person, so partial responses stay
  # missing -- under linked designs this can empty a subtest entirely
  newcols <- lapply(groups, function(g) rowSums(X[, g, drop = FALSE]))
  empty <- vapply(newcols, function(v) sum(!is.na(v)) < 10, TRUE)
  if (any(empty))
    stop("subtest(s) with fewer than 10 persons answering every member ",
         "item (missing data): ",
         paste(vapply(groups[empty], paste, "", collapse = "+"),
               collapse = ", "),
         "; choose items answered by the same persons")
  Xn <- cbind(X[, keep, drop = FALSE], do.call(cbind, newcols))
  colnames(Xn) <- c(keep, vapply(groups, paste, "", collapse = "+"))

  refit <- rasch(Xn, model = model, id = fit$person$id, factors = fit$factors,
                 n_groups = fit$n_groups)
  refit$notes <- c(refit$notes,
                   vapply(groups, function(g)
                     sprintf("subtest formed from: %s", paste(g, collapse = ", ")), ""))
  refit
}

#' Split items by a person factor to resolve DIF
#'
#' Replaces each nominated item with one item per level of a person factor,
#' each carrying that level's responses only (other levels missing). Every
#' group then receives its own item location, which resolves the invariance
#' violation flagged by \code{\link{dif_anova}}; the distance between the
#' split locations estimates the DIF size. The model is refitted with the
#' same settings.
#'
#' @param fit A fitted object from \code{\link{rasch}}.
#' @param items Character vector naming the item(s) to split.
#' @param by The name of a person factor nominated in the fit, or a grouping
#'   vector with one entry per person.
#' @return A new \code{\link{rasch}} fit in which each split item appears as
#'   \code{"item (level)"}, with the splits recorded in its notes.
#' @examples
#' set.seed(1); n <- 600
#' d <- seq(-2, 2, length.out = 8); g <- rep(c("a", "b"), each = n / 2)
#' sh <- matrix(0, n, 8); sh[g == "b", 3] <- 1
#' X <- matrix(rbinom(n * 8, 1, plogis(outer(rnorm(n), d, "-") - sh)), n, 8)
#' colnames(X) <- paste0("I", 1:8)
#' fit <- rasch(data.frame(X, grp = g), factors = "grp")
#' fit2 <- split_items(fit, "I3", by = "grp")
#' fit2$items$item
#' @export
split_items <- function(fit, items, by) {
  if (!inherits(fit, "rasch")) stop("split_items needs a rasch fit")
  if (inherits(fit, "rasch_mfrm"))
    stop("split items in the long-format data and refit rasch_mfrm instead")
  if (inherits(fit, "rasch_efrm"))
    stop("amend the source data and refit rasch_efrm instead")
  X <- fit$X
  bad <- setdiff(items, colnames(X))
  if (length(bad)) stop("item(s) not in the fit: ", paste(bad, collapse = ", "))
  if (is.character(by) && length(by) == 1L) {
    if (is.null(fit$factors) || !by %in% names(fit$factors))
      stop("'", by, "' is not a person factor nominated in the fit")
    grp <- fit$factors[[by]]
  } else {
    if (length(by) != nrow(X)) stop("'by' must have one entry per person")
    grp <- by
  }
  grp <- factor(grp)
  if (nlevels(grp) < 2) stop("the splitting factor needs at least two levels")

  keep <- setdiff(colnames(X), items)
  Xn <- X[, keep, drop = FALSE]
  for (it in items) for (lv in levels(grp)) {
    col <- X[, it]
    col[is.na(grp) | grp != lv] <- NA
    Xn <- cbind(Xn, col)
    colnames(Xn)[ncol(Xn)] <- paste0(it, " (", lv, ")")
  }
  refit <- rasch(Xn, model = fit$model, id = fit$person$id,
                 factors = fit$factors, n_groups = fit$n_groups)
  refit$notes <- c(refit$notes,
                   sprintf("item %s split by group into: %s", items,
                           paste0(items, " (", paste(levels(grp), collapse = "/"), ")")))
  refit
}

#' Resolve differential item functioning by iterative item splitting
#'
#' Splits DIF items one at a time, largest effect first, refitting after
#' each split, until no item shows significant DIF (or the anchor set would
#' fall too low). Splitting the item with the largest real DIF first
#' removes the artificial DIF it induces on other items (Andrich & Hagquist
#' 2012, 2015), so the procedure resolves genuine DIF without chasing the
#' artificial DIF that a single simultaneous pass would flag. Each split
#' resolves the item into one copy per group (or per factor-combination
#' cell for an interaction), with independent locations and thresholds, so
#' both uniform and non-uniform DIF are resolved together.
#'
#' @param fit A fitted object from \code{\link{rasch}} carrying person
#'   factors.
#' @param factors Person factors to test, as in \code{\link{dif_anova}};
#'   defaults to every nominated factor.
#' @param alpha Significance level for the adjusted probabilities.
#' @param p_adjust Multiplicity adjustment across items each round.
#' @param min_anchors Minimum number of original items to leave unsplit; the
#'   procedure stops before the anchor set falls below this (pervasive DIF is
#'   not artificial DIF). Default \code{max(3, items / 4)}.
#' @param max_splits Hard cap on the number of splits. Default: the number
#'   of items.
#' @return A list of class \code{"rmt_resolve_dif"}: the final resolved
#'   \code{fit}, the \code{splits} performed (order, item, factor, partial
#'   eta-squared, DIF magnitude in logits), the \code{stopped} reason, and
#'   the residual \code{dif} table for the final fit.
#' @references Andrich, D., & Hagquist, C. (2012). Real and artificial
#'   differential item functioning. \emph{Journal of Educational and
#'   Behavioral Statistics}, 37(3), 387-416.
#' @examples
#' set.seed(1); n <- 600
#' d <- seq(-2, 2, length.out = 8); g <- rep(c("a", "b"), each = n / 2)
#' sh <- matrix(0, n, 8); sh[g == "b", 3] <- 1.2      # one strong DIF item
#' X <- matrix(rbinom(n * 8, 1, plogis(outer(rnorm(n), d, "-") - sh)), n, 8)
#' colnames(X) <- paste0("I", 1:8)
#' fit <- rasch(data.frame(X, grp = g), factors = "grp")
#' resolve_dif(fit)$splits
#' @export
resolve_dif <- function(fit, factors = NULL, alpha = 0.05, p_adjust = "BH",
                        min_anchors = NULL, max_splits = NULL) {
  if (!inherits(fit, "rasch") || inherits(fit, c("rasch_mfrm", "rasch_efrm")))
    stop("resolve_dif needs an ordinary rasch fit with person factors")
  fac0 <- .dif_factors(fit, factors)
  fnames <- names(fac0)
  L0 <- ncol(fit$X)
  if (is.null(min_anchors)) min_anchors <- max(3L, L0 %/% 4L)
  if (is.null(max_splits)) max_splits <- L0

  # significant, non-superseded group terms of the current fit, with the
  # factors to split by and the partial eta-squared to rank on
  flagged <- function(cur) {
    keep <- intersect(fnames, names(cur$factors))
    if (!length(keep)) return(NULL)
    s <- dif_anova(cur, factors = keep, p_adjust = p_adjust, alpha = alpha)$summary
    s <- s[(s$uniform_DIF | s$nonuniform_DIF) & !s$superseded, , drop = FALSE]
    if (!nrow(s)) return(NULL)
    data.frame(item = s$item,
               vars = vapply(s$term, function(t)
                 paste(.term_vars(t), collapse = "+"), ""),
               eta2 = pmax(ifelse(s$uniform_DIF, s$eta2_uniform, 0),
                           ifelse(s$nonuniform_DIF, s$eta2_nonuniform, 0)),
               stringsAsFactors = FALSE)
  }
  # the original item a (possibly resolved) column belongs to
  base_of <- function(nm) sub(" \\(.*$", "", nm)

  cur <- fit
  splits <- list(); done <- character(0); stopped <- "no significant DIF remains"
  repeat {
    fl <- flagged(cur)
    if (is.null(fl) || !nrow(fl)) break
    fl <- fl[order(-fl$eta2), , drop = FALSE]
    # first flagged item-factor not already split, ranked by effect size
    pick <- NULL
    for (r in seq_len(nrow(fl))) {
      key <- paste(fl$item[r], fl$vars[r])
      if (!key %in% done) { pick <- fl[r, ]; break }
    }
    if (is.null(pick)) { stopped <- "remaining DIF cannot be resolved further"; break }
    if (length(splits) >= max_splits) { stopped <- "reached the split cap"; break }
    n_anchor <- L0 - length(unique(base_of(unlist(lapply(splits, `[[`, "item")))))
    if (n_anchor <= min_anchors) {
      stopped <- sprintf("stopped to keep %d anchor items (pervasive DIF is not artificial DIF)",
                         min_anchors)
      break
    }
    by_vars <- strsplit(pick$vars, "+", fixed = TRUE)[[1]]
    grp <- if (length(by_vars) == 1L) cur$factors[[by_vars]] else
      interaction(cur$factors[by_vars], sep = ":", drop = TRUE)
    # DIF magnitude in logits for the item about to be resolved
    mag <- tryCatch(max(abs(dif_size(cur, pick$item, by = grp)$pairs$difference)),
                    error = function(e) NA_real_)
    refit <- tryCatch(split_items(cur, pick$item, by = grp),
                      error = function(e) NULL)
    if (is.null(refit)) { done <- c(done, paste(pick$item, pick$vars)); next }
    splits[[length(splits) + 1L]] <- list(
      order = length(splits) + 1L, item = pick$item, factor = pick$vars,
      eta2 = pick$eta2, magnitude = mag)
    done <- c(done, paste(pick$item, pick$vars))
    cur <- refit
  }
  split_df <- if (length(splits))
    do.call(rbind, lapply(splits, function(s) as.data.frame(s,
                                                            stringsAsFactors = FALSE))) else
    data.frame(order = integer(), item = character(), factor = character(),
               eta2 = numeric(), magnitude = numeric())
  rownames(split_df) <- NULL
  final_dif <- tryCatch(flagged(cur), error = function(e) NULL)
  out <- list(fit = cur, splits = split_df, n_splits = nrow(split_df),
              stopped = stopped,
              n_remaining_dif = if (is.null(final_dif)) 0L else nrow(final_dif))
  class(out) <- "rmt_resolve_dif"
  out
}

#' @export
print.rmt_resolve_dif <- function(x, ...) {
  cat(sprintf("Iterative DIF resolution: %d split(s); %s\n",
              x$n_splits, x$stopped))
  if (x$n_splits) {
    d <- x$splits; d$eta2 <- round(d$eta2, 3); d$magnitude <- round(d$magnitude, 3)
    print(d, row.names = FALSE)
  }
  cat(sprintf("Remaining items with significant DIF: %d\n", x$n_remaining_dif))
  invisible(x)
}
