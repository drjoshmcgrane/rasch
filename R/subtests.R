# rasch :: structure amendments (subtests and item splitting)
# ===========================================================================
# Two standard remedies applied by re-analysis. Local dependence: combine the
# dependent items into a single polytomous super-item (the subtest) whose
# score is the sum of its members. Differential item functioning: split the
# offending item into one item per group level, each observed only by that
# group, so every group receives its own location and the invariance
# violation is resolved within the model.
# ===========================================================================

# One source label per current calibration item. An unchanged item maps to
# itself; every group-specific copy of a split item maps to the same source.
# Structural refits use this record to protect the unsplit reference set.
.split_source_map <- function(fit) {
  items <- colnames(fit$X)
  map <- fit$split_map
  if (is.null(map)) return(stats::setNames(items, items))
  if (is.null(names(map)) || anyDuplicated(names(map)))
    stop("the split-item provenance is malformed; refit from the source data")
  missing <- setdiff(items, names(map))
  if (length(missing))
    stop("the split-item provenance is incomplete for: ",
         paste(missing, collapse = ", "), "; refit from the source data")
  map <- as.character(map[items])
  names(map) <- items
  if (anyNA(map) || any(!nzchar(map)))
    stop("the split-item provenance contains a missing source item; refit ",
         "from the source data")
  map
}

.n_unsplit_sources <- function(map) {
  by_source <- split(names(map), unname(map))
  sum(vapply(names(by_source), function(source) {
    current <- by_source[[source]]
    length(current) == 1L && identical(current, source)
  }, logical(1)))
}

.split_source_items <- function(items, map) {
  items <- unique(as.character(items))
  source <- unname(map[items])
  source[is.na(source)] <- items[is.na(source)]
  unique(source)
}

#' Combine items into subtests and re-analyse
#'
#' Replaces each nominated item group by a polytomous super-item whose score is
#' the sum of its members, then refits the model. The function is commonly used
#' to examine item groups identified by \code{\link{residual_correlations}}.
#' Every total from zero to the sum of the component maxima must be observed;
#' otherwise the refit is refused rather than renumbering the superitem score.
#' The refit is also refused if calibration merges an observed category that
#' lacks conditional information, or changes a retained item's scoring.
#'
#' @param fit A fitted object from \code{\link{rasch}}.
#' @param groups A list of character vectors, each naming two or more items to
#'   combine; a single vector is also accepted.
#' @param model Model for the re-analysis; defaults to \code{"PCM"}, which is
#'   almost always required because subtests change the maximum scores.
#'   Explanatory fits retain their predictor restrictions on unchanged items
#'   and freely estimate the new superitems. They require \code{"PCM"}; an
#'   RSM restriction cannot be added through this argument.
#' @return A new \code{\link{rasch}} fit on the combined structure, with the
#'   combinations recorded in its notes. Person and item estimates are
#'   recalculated. Fit grouping, external anchors on unchanged items, keyed
#'   scoring, PCM constraints and optimisation controls are retained.
#'   Anchored items cannot be combined because the resulting superitem has no
#'   corresponding external anchor. A group-specific copy produced by
#'   \code{split_items()} cannot be combined: form the subtest before applying
#'   a DIF split. Existing split-item provenance is retained for items not
#'   included in a new subtest.
#' @examples
#' set.seed(1); Np <- 500; L <- 8
#' d <- seq(-2, 2, length.out = L)
#' X <- matrix(rbinom(Np * L, 1, plogis(outer(rnorm(Np), d, "-"))), Np, L)
#' colnames(X) <- paste0("I", 1:L)
#' X[, 5] <- ifelse(runif(Np) < 0.9, X[, 4], X[, 5])   # dependent pair
#' fit <- rasch(X)
#' fit2 <- combine_items(fit, list(c("I4", "I5")))
#' fit2$items$item
#' @seealso \code{\link{drop_items}} to remove an item rather than combine
#'   it, and \code{\link{residual_correlations}} for the dependence that
#'   motivates combining.
#' @export
combine_items <- function(fit, groups, model = "PCM") {
  if (!inherits(fit, "rasch")) stop("combine_items needs a rasch fit")
  .require_refittable_calibration(fit)
  if (!is.character(model) || length(model) != 1L ||
      !is.null(dim(model)) || !is.null(oldClass(model)) || is.na(model))
    stop("`model` must name one model: \"PCM\" or \"RSM\"", call. = FALSE)
  model <- match.arg(model, c("PCM", "RSM"))
  if (inherits(fit, "rasch_explanatory") && model != "PCM")
    stop("combine_items() retains the explanatory design; model = \"RSM\" ",
         "is not supported for explanatory fits. Use model = \"PCM\", ",
         "or refit the combined scores explicitly with rasch(..., model = \"RSM\")",
         call. = FALSE)
  if (is.character(groups) && length(groups)) groups <- list(groups)
  if (!is.list(groups) || !length(groups) ||
      any(!vapply(groups, length, 0L)))
    stop("`groups` must be item names, or a non-empty list of item groups, ",
         "each naming at least one item")
  groups <- lapply(groups, function(g) {
    if (is.factor(g)) g <- as.character(g)
    if (!is.character(g) || !is.null(dim(g)) || anyNA(g) ||
        any(!nzchar(trimws(g))))
      stop("each subtest must be an ordinary vector of non-missing item names")
    g
  })
  if (inherits(fit, "rasch_mfrm"))
    stop("combine items in the long-format data and refit rasch_mfrm instead")
  if (inherits(fit, "rasch_efrm"))
    stop("combine items in the source data and refit rasch_efrm instead")
  X <- fit$X
  used <- unlist(groups)
  bad <- setdiff(used, colnames(X))
  if (length(bad)) stop("item(s) not in the fit: ", paste(bad, collapse = ", "))
  if (anyDuplicated(used)) stop("an item appears in more than one subtest")
  short <- vapply(groups, length, 1L) < 2
  if (any(short)) stop("each subtest needs at least two items")
  split_provenance <- NULL
  if (!is.null(fit$split_map)) {
    split_provenance <- .split_source_map(fit)
    split_members <- used[unname(split_provenance[used]) != used]
    if (length(split_members))
      stop("group-specific split item(s) cannot be combined into a subtest: ",
           paste(unique(split_members), collapse = ", "),
           "; form the subtest before applying DIF splits")
  }

  keep <- setdiff(colnames(X), used)
  anchors <- fit$refit_spec$anchors
  if (!is.null(anchors) && any(as.character(anchors$item) %in% used))
    stop("an anchored item cannot be combined into a subtest; refit with ",
         "anchors on items that remain unchanged")
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
  # a subtest total spans the SUM of its members' maxima, so a code that
  # could never collide with a single item's score can fall inside the
  # superitem's range: applying it there would delete complete responses
  # and renumber the categories that remain
  group_max <- vapply(groups, function(g)
    sum(fit$m[match(g, fit$items$item)]), 0)
  smax <- max(group_max)
  all_codes <- (fit$refit_spec %||% list())$na_codes %||% -1
  numeric_codes <- suppressWarnings(as.numeric(as.character(all_codes)))
  dropped_codes <- all_codes[is.finite(smax) & is.finite(numeric_codes) &
                             numeric_codes >= 0 & numeric_codes <= smax]
  super_names <- vapply(groups, paste, "", collapse = "+")
  if (anyDuplicated(super_names) || any(super_names %in% keep))
    stop("the generated subtest names duplicate an existing item name")
  if (!is.null(split_provenance)) {
    kept_sources <- unique(unname(split_provenance[keep]))
    if (any(super_names %in% kept_sources))
      stop("the generated subtest names duplicate an existing split-item ",
           "source name")
  }
  colnames(Xn) <- c(keep, super_names)
  score_max <- c(stats::setNames(fit$m[match(keep, fit$items$item)], keep),
                 stats::setNames(group_max, super_names))
  .require_score_structure(Xn, score_max, "the subtest refit")

  refit <- if (inherits(fit, "rasch_explanatory")) {
    inherit <- c(stats::setNames(keep, keep),
                 stats::setNames(vapply(groups, `[`, "", 1L), super_names))
    .explanatory_refit_modified(fit, Xn, inherit = inherit,
                                fully_relaxed = super_names)
  } else .rasch_refit(fit, Xn, model = model,
                      score_max = score_max)
  if (!isTRUE(refit$est$converged))
    stop("the subtest calibration did not converge; the combined analysis is unavailable")
  .require_fitted_score_structure(refit, score_max, "the subtest refit")
  if (length(dropped_codes))
    refit$notes <- c(refit$notes, paste0(
      "missing-data code(s) ", paste(dropped_codes, collapse = ", "),
      " not applied to the subtest scores: they fall inside a subtest's ",
      "score range, where they would delete complete responses"))
  old_map <- fit$subtest_map %||% list()
  old_binary <- fit$subtest_binary %||% logical(0)
  members <- lapply(groups, function(g) unlist(lapply(g, function(it)
    old_map[[it]] %||% it), use.names = FALSE))
  binary <- vapply(groups, function(g) all(vapply(g, function(it) {
    if (it %in% names(old_binary)) isTRUE(old_binary[[it]])
    else fit$m[match(it, fit$items$item)] == 1L
  }, logical(1))), logical(1))
  keep_map <- old_map[intersect(names(old_map), keep)]
  keep_binary <- old_binary[intersect(names(old_binary), keep)]
  refit$subtest_map <- c(keep_map, stats::setNames(members, super_names))
  refit$subtest_binary <- c(keep_binary,
                            stats::setNames(binary, super_names))
  if (!is.null(split_provenance)) {
    refit$split_map <- c(split_provenance[keep],
                         stats::setNames(super_names, super_names))
  }
  refit$notes <- c(refit$notes,
                   vapply(groups, function(g)
                     sprintf("subtest formed from: %s", paste(g, collapse = ", ")), ""))
  refit
}

#' Split items by a person factor to resolve DIF
#'
#' Replaces each nominated item with one item per estimable level of a person
#' factor, each carrying that level's responses only (other levels missing).
#' A level must contain every score from zero to the fitted item maximum;
#' levels missing a score are omitted and recorded in the notes. The split
#' is refused if calibration subsequently merges a category lacking conditional
#' information within a retained group. Every retained
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
#'   \code{"item (level)"}, with the splits recorded in its notes. Person and
#'   item estimates are recalculated with the fitted grouping, keyed scoring,
#'   anchors on unchanged items, PCM constraints and optimisation controls.
#'   An anchored item cannot be split.
#' @examples
#' set.seed(1); n <- 600
#' d <- seq(-2, 2, length.out = 8); g <- rep(c("a", "b"), each = n / 2)
#' sh <- matrix(0, n, 8); sh[g == "b", 3] <- 1
#' X <- matrix(rbinom(n * 8, 1, plogis(outer(rnorm(n), d, "-") - sh)), n, 8)
#' colnames(X) <- paste0("I", 1:8)
#' fit <- rasch(data.frame(X, grp = g), factors = "grp")
#' fit2 <- split_items(fit, "I3", by = "grp")
#' fit2$items$item
#' @seealso \code{\link{resolve_dif}}, which applies this iteratively;
#'   \code{\link{drop_items}}, which removes an item rather than resolving
#'   it; and \code{\link{dif_anova}}, which identifies the items to split.
#' @export
split_items <- function(fit, items, by) {
  if (!inherits(fit, "rasch")) stop("split_items needs a rasch fit")
  .require_refittable_calibration(fit)
  if (inherits(fit, "rasch_mfrm"))
    stop("split items in the long-format data and refit rasch_mfrm instead")
  if (inherits(fit, "rasch_efrm"))
    stop("amend the source data and refit rasch_efrm instead")
  items <- .structural_item_names(items, "split")
  X <- fit$X
  bad <- setdiff(items, colnames(X))
  if (length(bad)) stop("item(s) not in the fit: ", paste(bad, collapse = ", "))
  if (!is.null(dim(by)) || !is.atomic(by) || is.complex(by))
    stop("`by` must be one person-factor name or an ordinary grouping vector",
         call. = FALSE)
  if (is.character(by) && length(by) == 1L) {
    if (is.null(fit$factors) || !by %in% names(fit$factors))
      stop("'", by, "' is not a person factor nominated in the fit")
    grp <- fit$factors[[by]]
  } else {
    if (length(by) != nrow(X)) stop("'by' must have one entry per person")
    grp <- by
  }
  if (is.character(grp) || is.factor(grp)) {
    grp <- .role_text_values(grp)
    grp[!is.na(grp) & !nzchar(grp)] <- NA_character_
  }
  grp <- factor(grp)
  if (nlevels(grp) < 2) stop("the splitting factor needs at least two levels")

  keep <- setdiff(colnames(X), items)
  anchors <- fit$refit_spec$anchors
  if (!is.null(anchors) && any(as.character(anchors$item) %in% items))
    stop("an anchored item cannot be split; nominate invariant anchors and ",
         "refit before resolving DIF")
  Xn <- as.data.frame(X[, keep, drop = FALSE], check.names = FALSE)
  key_extra <- NULL
  made <- list()
  split_levels <- lapply(items, function(it) {
    mi <- fit$m[match(it, fit$items$item)]
    ok <- vapply(levels(grp), function(lv) {
      z <- X[!is.na(grp) & grp == lv, it]
      identical(as.integer(sort(unique(z[!is.na(z)]))), seq.int(0L, mi))
    }, logical(1))
    lv <- levels(grp)[ok]
    if (length(lv) < 2L)
      stop("item ", it, " cannot be split without changing its score ",
           "structure: fewer than two factor levels observe every score ",
           "from 0 to ", mi)
    lv
  })
  names(split_levels) <- items
  omitted_levels <- lapply(items, function(it)
    setdiff(levels(grp), split_levels[[it]]))
  names(omitted_levels) <- items
  for (it in items) for (lv in split_levels[[it]]) {
    col <- X[, it]
    col[is.na(grp) | grp != lv] <- NA
    nm <- paste0(it, " (", lv, ")")
    if (nm %in% names(Xn)) stop("generated split-item name already exists: ", nm)
    if (!inherits(fit, "rasch_explanatory") &&
        !is.null(fit$mc) && it %in% colnames(fit$mc$raw)) {
      col <- fit$mc$raw[, it]
      col[is.na(grp) | grp != lv] <- NA
      kr <- fit$refit_spec$key
      if (is.data.frame(kr)) {
        kr <- kr[as.character(kr$item) == it, , drop = FALSE]
        kr$item <- nm
        key_extra <- rbind(key_extra, kr)
      }
    }
    Xn[[nm]] <- col
    made[[it]] <- c(made[[it]], nm)
  }
  refit <- if (inherits(fit, "rasch_explanatory")) {
    inherit <- c(stats::setNames(keep, keep),
                 unlist(lapply(items, function(it)
                   stats::setNames(rep(it, length(made[[it]])), made[[it]]))))
    .explanatory_refit_modified(fit, Xn, inherit = inherit,
                                location_relaxed = unlist(made,
                                                          use.names = FALSE))
  } else .rasch_refit(fit, Xn, key_extra = key_extra)
  if (!isTRUE(refit$est$converged))
    stop("the split-item calibration did not converge; the resolved analysis is unavailable")
  # Observing every category before fitting is not sufficient: a category
  # confined to extreme patterns within one group can be merged in preparation.
  # The split must retain each copy's original score scale after calibration.
  expected_max <- c(stats::setNames(fit$m[match(keep, fit$items$item)], keep),
    unlist(lapply(items, function(it) stats::setNames(
      rep(fit$m[match(it, fit$items$item)], length(made[[it]])), made[[it]]))))
  retained_max <- refit$m[match(names(expected_max), refit$items$item)]
  changed <- is.na(retained_max) | retained_max != expected_max
  if (any(changed))
    stop("the split refit cannot preserve the fitted score structure for ",
         paste(names(expected_max)[changed], collapse = ", "),
         "; a category lacks conditional information within a group. ",
         "Pool the affected group or leave this split unresolved", call. = FALSE)
  if (length(fit$subtest_map)) {
    sm <- fit$subtest_map[intersect(names(fit$subtest_map), keep)]
    sb <- fit$subtest_binary[intersect(names(fit$subtest_binary), keep)]
    for (it in intersect(items, names(fit$subtest_map))) {
      sm[made[[it]]] <- rep(list(fit$subtest_map[[it]]), length(made[[it]]))
      sb[made[[it]]] <- rep(isTRUE(fit$subtest_binary[[it]]),
                            length(made[[it]]))
    }
    refit$subtest_map <- sm
    refit$subtest_binary <- sb
  }
  old_map <- .split_source_map(fit)
  refit$split_map <- c(old_map[keep], unlist(lapply(items, function(it)
    stats::setNames(rep(unname(old_map[it]), length(made[[it]])), made[[it]]))))
  refit$notes <- c(refit$notes,
                   vapply(items, function(it)
                     sprintf("item %s split by group into: %s", it,
                             paste(made[[it]], collapse = ", ")), ""),
                   unlist(lapply(items, function(it) {
                     lv <- omitted_levels[[it]]
                     if (!length(lv)) return(character(0))
                     sprintf(paste0("item %s was not split for unavailable ",
                                    "level(s): %s (not every fitted score ",
                                    "category was observed)"),
                             it, paste(lv, collapse = ", "))
                   }), use.names = FALSE))
  refit
}

#' Resolve differential item functioning by iterative item splitting
#'
#' Splits items with uniform DIF one at a time, beginning with the largest
#' estimated effect,
#' and refits after each split. This order addresses the artificial DIF that a
#' large departure can induce in otherwise invariant items (Andrich and
#' Hagquist 2012, 2015). Each split gives the item a separate location in every
#' factor cell. A PCM also estimates the split copies' thresholds separately;
#' an RSM retains its common rating-scale threshold structure. A location split does not model a
#' group-specific discrimination, so items with non-uniform DIF are left for
#' review rather than being made untestable by a split. The procedure stops
#' when no resolvable uniform DIF remains or the remaining unsplit reference
#' set reaches \code{min_anchors}. Items fixed by external anchors are not
#' split.
#'
#' @param fit A fitted object from \code{\link{rasch}}.
#' @param factors Person factors to test, as in \code{\link{dif_anova}};
#'   defaults to every nominated factor.
#' @param alpha Significance level for the adjusted probabilities.
#' @param p_adjust Multiplicity adjustment for the DIF tests in each round.
#' @param min_n Minimum distinct responders required in every item-by-factor
#'   cell before an automatic split is allowed. Repeated response rows from
#'   one person count once within a cell. The omnibus DIF test determines
#'   whether a split is needed; pairwise follow-ups describe where the
#'   difference lies but are not a second significance gate.
#' @param min_anchors Minimum number of original items to leave unsplit as the
#'   internal reference set. The procedure stops before this set becomes
#'   smaller; pervasive DIF is not artificial DIF. Default
#'   \code{max(3, items / 4)}.
#' @param max_splits Hard cap on the number of splits. Default: the number
#'   of items.
#' @param effects \code{"main"} fits the factors additively;
#'   \code{"factorial"} also tests their interactions. The same model is
#'   used at every round and in the final DIF assessment.
#' @return A list of class \code{"rasch_resolve_dif"}: the final resolved
#'   \code{fit}, the \code{splits} performed (order, item, factor, partial
#'   eta-squared, source item, DIF magnitude in logits), the \code{stopped}
#'   reason, the residual \code{dif} table, and the number of distinct source
#'   items that still show DIF in the final fit. \code{n_untested} counts the
#'   item-term tests the final assessment could not estimate although the
#'   design could answer them; those terms are reported as neither DIF nor no
#'   DIF, so the remaining-DIF count is a lower bound whenever
#'   \code{n_untested} is positive. A split copy answered in one level of its
#'   splitting factor only is not counted: its term is structurally absent,
#'   not lost. Both counts are \code{NA} when no item-term test was estimable
#'   at all. \code{effects} records the factor model used.
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
#' @seealso \code{\link{split_items}} for a single split,
#'   \code{\link{drop_items}} to remove an item instead, and
#'   \code{\link{dif_anova}} for the test it resolves.
#' @export
resolve_dif <- function(fit, factors = NULL, alpha = 0.05, p_adjust = "holm",
                        min_n = 20L, min_anchors = NULL, max_splits = NULL,
                        effects = c("main", "factorial")) {
  effects <- match.arg(effects)
  .check_dif_args(alpha, p_adjust, min_n = min_n)
  if (!inherits(fit, "rasch") || inherits(fit, c("rasch_mfrm", "rasch_efrm")))
    stop("resolve_dif needs an ordinary rasch fit with person factors")
  fac0 <- .dif_factors(fit, factors)
  fnames <- names(fac0)
  if (!length(fnames)) stop("`factors` must include at least one person factor")
  source_map <- .split_source_map(fit)
  L0 <- length(unique(unname(source_map)))
  if (L0 < 3L) stop("DIF resolution needs at least three fitted items")
  if (is.null(min_anchors))
    min_anchors <- min(L0 - 1L, max(3L, L0 %/% 4L))
  if (is.null(max_splits)) max_splits <- L0
  if (length(min_anchors) != 1L || !is.numeric(min_anchors) ||
      is.complex(min_anchors) || !is.null(dim(min_anchors)) ||
      !is.null(oldClass(min_anchors)) || !is.finite(min_anchors) ||
      min_anchors < 2L || min_anchors != floor(min_anchors) ||
      min_anchors >= L0)
    stop("min_anchors must be a whole number from 2 to one fewer than the ",
         "number of fitted items")
  if (length(max_splits) != 1L || !is.numeric(max_splits) ||
      is.complex(max_splits) || !is.null(dim(max_splits)) ||
      !is.null(oldClass(max_splits)) || !is.finite(max_splits) ||
      max_splits < 0L || max_splits != floor(max_splits) ||
      max_splits > .Machine$integer.max)
    stop("max_splits must be a non-negative whole number")

  # significant, non-superseded group terms of the current fit, with the
  # factors to split by and the partial eta-squared to rank on
  last_da <- NULL
  flagged <- function(cur, resolvable_only = TRUE) {
    # Splitting preserves response rows. Retain the supplied factor values,
    # including external metadata and replacements for a stored factor,
    # instead of silently reselecting columns from the original fit.
    da <- dif_anova(cur, factors = fac0, effects = effects,
                    p_adjust = p_adjust, alpha = alpha)
    # An item-term test that is not estimable is reported as NA and reads
    # as "not significant" here. Keep the assessment itself so the caller
    # can tell no DIF apart from no test.
    last_da <<- da
    s <- da$summary
    # Splitting supplies cell-specific locations and thresholds, not
    # cell-specific slopes. Restrict automatic resolution to uniform-only DIF;
    # otherwise the split removes the between-group overlap and merely makes a
    # non-uniform violation disappear from the DIF table.
    take <- if (resolvable_only)
      s$uniform_DIF & !s$nonuniform_DIF & !s$superseded
    else
      (s$uniform_DIF | s$nonuniform_DIF) & !s$superseded
    vars <- da$summary_factors[take]
    s <- s[take, , drop = FALSE]
    if (!nrow(s)) return(NULL)
    out <- data.frame(
      item = s$item, factor = vapply(vars, .dif_term_label, ""),
      eta2 = pmax(ifelse(s$uniform_DIF, s$eta2_uniform, 0),
                  ifelse(s$nonuniform_DIF, s$eta2_nonuniform, 0)),
      uniform = s$uniform_DIF, nonuniform = s$nonuniform_DIF,
      stringsAsFactors = FALSE)
    attr(out, "term_factors") <- vars
    out
  }
  cur <- fit
  splits <- list(); done <- character(0); skipped <- character(0)
  skipped_anchor <- character(0)
  stopped <- "no significant DIF remains"
  # TRUE while `stopped` is a verdict on the DIF that remains, which the
  # final assessment can contradict; FALSE once the loop stops for a reason
  # of its own that stays true whatever that assessment shows.
  stopped_is_verdict <- TRUE
  repeat {
    fl <- flagged(cur)
    if (is.null(fl) || !nrow(fl)) {
      remaining <- flagged(cur, resolvable_only = FALSE)
      if (!is.null(remaining) && any(remaining$nonuniform))
        stopped <- paste("no resolvable uniform DIF remains;",
                         "non-uniform DIF requires item review")
      break
    }
    vars <- attr(fl, "term_factors")
    ord <- order(-fl$eta2)
    fl <- fl[ord, , drop = FALSE]
    vars <- vars[ord]
    # first flagged item-factor not already split, ranked by effect size
    pick <- NULL
    for (r in seq_len(nrow(fl))) {
      vk <- paste0(nchar(vars[[r]], type = "bytes"), ":", vars[[r]], ";",
                   collapse = "")
      key <- paste0(nchar(fl$item[r], type = "bytes"), ":", fl$item[r], ";", vk)
      if (!key %in% done) { pick <- fl[r, ]; pick_vars <- vars[[r]]; break }
    }
    if (is.null(pick)) { stopped <- "remaining DIF cannot be resolved further"; break }
    if (length(splits) >= max_splits) {
      stopped <- "reached the split cap"; stopped_is_verdict <- FALSE; break
    }
    n_anchor <- .n_unsplit_sources(.split_source_map(cur))
    if (n_anchor <= min_anchors) {
      stopped <- sprintf("stopped to keep %d anchor items (pervasive DIF is not artificial DIF)",
                         min_anchors)
      stopped_is_verdict <- FALSE
      break
    }
    by_vars <- pick_vars
    grp <- if (length(by_vars) == 1L) fac0[[by_vars]] else
      .factor_cells(fac0[by_vars], sep = ":")
    grp <- factor(grp)
    external <- cur$refit_spec$anchors
    if (!is.null(external) &&
        pick$item %in% as.character(external$item)) {
      skipped_anchor <- c(skipped_anchor, paste(pick$item, pick$factor))
      done <- c(done, key)
      next
    }
    # The adjusted omnibus term is the decision to split. Pairwise follow-ups
    # answer a different question and cannot be imposed as a second rejection
    # gate: a multi-df omnibus may be significant even when no one of many
    # Holm-adjusted pairs is. Check only whether every level to be split has
    # adequate, non-boundary support and a common score-category structure.
    intended_levels <- levels(droplevels(grp[!is.na(grp)]))
    support <- tryCatch(.dif_resolve(cur, pick$item, grp, min_n),
                        error = function(e) NULL)
    support_ok <- !is.null(support) &&
      setequal(support$levs, intended_levels) &&
      !any(support$weak %in% TRUE) &&
      !identical(support$score_compatible, FALSE) &&
      all(is.finite(support$loc))
    if (!support_ok) {
      skipped <- c(skipped, paste(pick$item, pick$factor))
      done <- c(done, key)
      next
    }
    # Use the complete-design marginal follow-up for the reported logit
    # magnitude. Its significance does not decide whether the split proceeds.
    dp <- tryCatch(dif_posthoc(
      cur, pick$item, term = by_vars,
      factors = fac0,
      p_adjust = p_adjust, alpha = alpha, min_n = min_n),
      error = function(e) NULL)
    # DIF magnitude in logits, over the trustworthy (non-weak) pairs only
    mag <- if (!is.null(dp)) {
      d <- abs(dp$table$estimate[is.finite(dp$table$estimate) &
                                   (is.finite(dp$table$se) |
                                      is.finite(dp$table$statistic))])
      if (length(d)) max(d) else NA_real_
    } else NA_real_
    refit <- tryCatch(split_items(cur, pick$item, by = grp),
                      error = function(e) NULL)
    if (is.null(refit)) { done <- c(done, key); next }
    base_map <- .split_source_map(cur)
    splits[[length(splits) + 1L]] <- list(
      order = length(splits) + 1L, item = pick$item, factor = pick$factor,
      base_item = unname(base_map[pick$item]), eta2 = pick$eta2,
      magnitude = mag)
    done <- c(done, key)
    cur <- refit
  }
  split_df <- if (length(splits))
    do.call(rbind, lapply(splits, function(s) as.data.frame(s,
                                                            stringsAsFactors = FALSE))) else
    data.frame(order = integer(), item = character(), factor = character(),
               base_item = character(), eta2 = numeric(), magnitude = numeric())
  rownames(split_df) <- NULL
  final_dif <- tryCatch(flagged(cur, resolvable_only = FALSE),
    error = function(e) stop("the final DIF assessment failed: ",
                             conditionMessage(e), call. = FALSE))
  remaining_items <- if (is.null(final_dif) || !nrow(final_dif)) {
    character(0)
  } else {
    .split_source_items(final_dif$item, .split_source_map(cur))
  }
  # dif_anova reports a test it cannot estimate as an NA row rather than an
  # error, and the summary reads NA as not significant. An item term the
  # design could have answered is therefore silence, not a clean bill of
  # health. A split copy is the opposite case: it lives in one level of its
  # splitting factor only, so its item-term test is structurally absent
  # rather than lost, and says nothing about the DIF that remains. Count
  # only the tests that were both requested and answerable, and withhold
  # the counts outright when none of them was estimable.
  s_fin <- last_da$summary
  icol <- match(s_fin$item, colnames(cur$X))
  answerable <- vapply(seq_len(nrow(s_fin)), function(r) {
    if (is.na(icol[r])) return(TRUE)
    seen <- !is.na(cur$X[, icol[r]])
    all(vapply(last_da$summary_factors[[r]], function(v) {
      lv <- as.character(fac0[seen, v])
      length(unique(lv[!is.na(lv)])) >= 2L
    }, TRUE))
  }, TRUE)
  tested <- is.finite(s_fin$p_uniform_adj) | is.finite(s_fin$p_nonuniform_adj)
  n_untested <- sum(answerable & !tested)
  no_test <- !any(tested)
  if (no_test) {
    unknown <- paste("no item-term test in the final DIF assessment was",
                     "estimable; remaining DIF is unknown")
    stopped <- if (stopped_is_verdict) unknown else paste0(stopped, "; ", unknown)
  } else if (n_untested)
    stopped <- paste0(stopped, sprintf(
      paste("; %d of %d item-term test(s) in the final assessment were not",
            "estimable, so remaining DIF may be understated"),
      n_untested, sum(answerable)))
  notes <- character(0)
  # The assessment's own notes are kept verbatim, its count of terms it could
  # not estimate included: those terms stay in its adjusted-probability
  # family, so they make the verdict above conservative. That is a fact about
  # the family, not the claim about remaining DIF that n_untested withholds.
  if (length(last_da$notes))
    notes <- c(notes, paste("final DIF assessment:", last_da$notes))
  if (length(skipped)) notes <- c(notes,
    sprintf("%d flagged item-factor(s) not split because one or more cells had fewer than min_n distinct responders, weak boundary estimates, or incompatible response categories (%s)",
            length(skipped), paste(skipped, collapse = "; ")))
  if (length(skipped_anchor)) notes <- c(notes,
    sprintf("%d externally anchored item-factor(s) not split (%s)",
            length(skipped_anchor), paste(skipped_anchor, collapse = "; ")))
  out <- list(algorithm = "factor-design-resolution-2",
              fit = cur, splits = split_df, n_splits = nrow(split_df),
              stopped = stopped, dif = final_dif, notes = notes,
              effects = effects,
              n_remaining_dif = if (no_test) NA_integer_ else
                length(remaining_items),
              n_untested = as.integer(n_untested),
              n_nonuniform = if (no_test) NA_integer_ else
                if (is.null(final_dif)) 0L else
                sum(final_dif$nonuniform %in% TRUE))
  out <- .tag_tables(out)
  class(out) <- "rasch_resolve_dif"
  out
}

#' @export
print.rasch_resolve_dif <- function(x, ...) {
  cat(sprintf("Iterative DIF resolution: %d split(s); %s\n",
              x$n_splits, x$stopped))
  if (x$n_splits) {
    d <- x$splits; d$eta2 <- round(d$eta2, 3); d$magnitude <- round(d$magnitude, 3)
    print(d, row.names = FALSE)
  }
  cat(sprintf("Remaining items with significant DIF: %d\n", x$n_remaining_dif))
  if (isTRUE(x$n_untested > 0))
    cat(sprintf("Item-term tests not estimable: %d\n", x$n_untested))
  if (isTRUE(x$n_nonuniform > 0))
    cat(sprintf("Non-uniform item-factor findings requiring review: %d\n",
                x$n_nonuniform))
  if (length(x$notes)) cat("\n", paste(x$notes, collapse = "\n"), "\n", sep = "")
  invisible(x)
}
