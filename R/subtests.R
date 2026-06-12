# RaschR :: structure amendments (subtests and item splitting)
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
  newcols <- lapply(groups, function(g) rowSums(X[, g, drop = FALSE]))
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
