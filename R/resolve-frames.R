# rasch :: resolving an item that does not hold across frames
# ===========================================================================
# The frame model gives each item one location, shared across the frames it
# appears in and scaled by the frame unit. That tie is what makes the frames
# comparable, and it is an assumption. When frame_invariance() finds an item
# the tie does not hold for, there are two remedies. Removing the item takes
# it out of every frame, so it stops measuring anyone, including in the
# frames it behaved perfectly well in. Resolving it gives it a separate
# location in each frame: it stops linking the frames, which is what the
# diagnosis says is wrong with it, and goes on measuring the person within
# their own frame.
#
# Resolving costs a parameter per extra frame and removes the item from the
# link, so the group units then rest on the items that remain common. It is
# the better remedy when the item measures well within frames and the worse
# one when it does not measure at all.
# ===========================================================================

#' Resolve items that do not hold across frames
#'
#' Gives each named item a separate location in every frame in which it was
#' administered, then refits the EFRM.
#'
#' @details
#' A resolved item continues to contribute to person measurement within each
#' frame but no longer constrains the link between those frames. Its versions
#' are named \code{"item (frame)"}. The remaining common items and the linked
#' set design must still identify the frame units; otherwise the refit is
#' refused by the model's connectivity and rank checks.
#' Each set must also retain links between its groups of item versions to
#' identify their relative origins. A unit link supplied by another set
#' cannot replace these origin links.
#'
#' Resolve an item when its within-frame measurement remains defensible but
#' its cross-frame location does not. This refit does not estimate a separate
#' discrimination and therefore does not resolve a discrimination-only flag
#' from \code{\link{frame_invariance}}. Review or remove such an item instead.
#' Use \code{\link{drop_items}} when the item should no longer contribute to
#' measurement. Each resolved version must observe every score category of the
#' source item. If it does not, the refit is refused rather than renumbering
#' that frame's scores.
#'
#' @param fit A fitted object from \code{\link{rasch_efrm}}.
#' @param items Item names to resolve.
#' @param boot_reps Bootstrap replicates for the refit. The default retains
#'   the fitted specification; a number overrides it.
#' @return A refitted object of class \code{"rasch_efrm"}, carrying a note
#'   for each item resolved. The resolved versions appear in the item table
#'   as \code{"item (frame)"}.
#' @seealso \code{\link{frame_invariance}}, which identifies the items to
#'   resolve; \code{\link{drop_items}}, which removes an item instead; and
#'   \code{\link{split_items}}, the equivalent for an ordinary fit.
#' @examples
#' d <- simulate_efrm(n_per_group = 200, items_per_set = 6, n_sets = 2,
#'                    n_groups = 2, set_unit_ratio = 1.3, seed = 5)
#' tr <- attr(d, "truth")
#' fit <- rasch_efrm(d, item_sets = tr$item_sets, groups = "group",
#'                   id = "id", boot_reps = 0)
#' fit2 <- resolve_frames(fit, "S1I02", boot_reps = 0)
#' grep("S1I02", fit2$items$item, value = TRUE)
#' @export
resolve_frames <- function(fit, items, boot_reps = NULL) {
  if (!inherits(fit, "rasch"))
    stop("resolve_frames needs a fit from rasch_efrm()")
  if (!inherits(fit, "rasch_efrm"))
    stop("resolve_frames needs a frame model; for an ordinary fit with ",
         "person factors use split_items() instead")
  items <- .structural_item_names(items, "resolve")

  all_items <- names(fit$set_of)
  bad <- setdiff(items, all_items)
  if (length(bad))
    stop("item(s) not in the fit: ", paste(bad, collapse = ", "),
         "; the fit holds: ", paste(utils::head(all_items, 8), collapse = ", "),
         if (length(all_items) > 8) ", ..." else "")

  grp <- .frame_group_values(fit)
  glev <- levels(factor(grp))
  gvec <- as.character(grp)

  src <- .efrm_source_matrix(fit, all_items)
  source_max <- .efrm_source_maxima(fit)

  cols <- list()
  new_set <- character()
  source_base <- character()
  made <- list()
  for (it in all_items) {
    if (!it %in% items) {
      cols[[it]] <- src[, it]
      new_set[it] <- fit$set_of[[it]]
      source_base[it] <- it
      next
    }
    for (lv in glev) {
      v <- src[, it]
      v[is.na(gvec) | gvec != lv] <- NA
      if (all(is.na(v))) next          # the item was not taken in this frame
      nm <- paste0(it, " (", lv, ")")
      if (nm %in% all_items || nm %in% names(cols))
        stop("generated resolved-item name already exists: ", nm)
      cols[[nm]] <- v
      new_set[nm] <- fit$set_of[[it]]
      source_base[nm] <- it
      made[[it]] <- c(made[[it]], lv)
    }
  }

  source <- as.data.frame(cols, check.names = FALSE,
                         stringsAsFactors = FALSE)
  resolved_max <- stats::setNames(
    as.numeric(source_max[source_base[names(source)]]), names(source))
  .require_score_structure(source, resolved_max,
                           "the resolved-frame refit")
  refit <- .efrm_refit(fit, source, new_set, boot_reps = boot_reps,
                       score_max = resolved_max)
  if (!isTRUE(refit$est$converged))
    stop("the resolved frame calibration did not converge; the sensitivity analysis is unavailable")
  if (!.efrm_link_converged(refit))
    stop("the resolved set-unit link did not converge; the sensitivity analysis is unavailable")
  refit$notes <- c(refit$notes,
                   vapply(names(made), function(it)
                     sprintf("item %s resolved by frame into: %s", it,
                             paste0(it, " (",
                                    paste(made[[it]], collapse = "/"), ")")),
                     character(1)))
  refit
}
