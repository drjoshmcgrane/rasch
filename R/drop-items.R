# rasch :: removing items from a fitted analysis
# ===========================================================================

.frame_group_vars <- function(fit) {
  fg <- fit$frame_group
  if (is.null(fg) || !length(fg)) stop("the frame variables are not recorded")
  if (length(fg) > 1L) fg[-1L] else fg[1L]
}

.frame_group_values <- function(fit) {
  fg <- fit$frame_group[1L]
  if (is.null(fit$factors) || !fg %in% names(fit$factors))
    stop("the fitted frame membership is not available")
  fit$factors[[fg]]
}

.structural_item_names <- function(items, action = "change") {
  if (is.factor(items)) items <- as.character(items)
  if (!is.character(items) || !is.null(dim(items)) || !length(items) ||
      anyNA(items) ||
      any(!nzchar(trimws(items))))
    stop("`items` must contain at least one non-missing item name to ",
         action, call. = FALSE)
  if (anyDuplicated(items))
    stop("item(s) named more than once: ",
         paste(unique(items[duplicated(items)]), collapse = ", "),
         call. = FALSE)
  items
}

# Structural refits must not let ordinary data preparation renumber an item's
# scores. `maxima` is a named vector giving the intended maximum score for
# each column being checked. A missing intermediate category makes that
# parameterisation unestimable in the refit; refusing it is preferable to
# silently fitting a different scoring structure.
.require_score_structure <- function(source, maxima, context = "the refit") {
  if (is.null(maxima) || !length(maxima)) return(invisible(TRUE))
  source <- as.data.frame(source, check.names = FALSE,
                          stringsAsFactors = FALSE)
  maxima <- maxima[!is.na(names(maxima)) & names(maxima) %in% names(source)]
  bad <- character(0); found <- character(0)
  for (nm in names(maxima)) {
    v <- source[[nm]]
    z <- suppressWarnings(as.numeric(as.character(v)))
    obs <- sort(unique(z[is.finite(z) & z >= 0]))
    expected <- seq.int(0L, as.integer(maxima[[nm]]))
    if (!identical(as.integer(obs), expected)) {
      bad <- c(bad, nm)
      found <- c(found, paste(obs, collapse = ","))
    }
  }
  if (length(bad))
    stop(context, " cannot preserve the fitted score structure for ",
         paste(sprintf("%s (expected 0:%d; observed [%s])", bad,
                       as.integer(maxima[bad]), found), collapse = "; "),
         "; add observations in the missing categories, pool the affected ",
         "group, or leave that structural change unresolved", call. = FALSE)
  invisible(TRUE)
}

.require_fitted_score_structure <- function(fit, maxima,
                                            context = "the structural refit") {
  if (!length(maxima)) return(invisible(TRUE))
  actual <- fit$m[match(names(maxima), fit$items$item)]
  changed <- is.na(actual) | actual != maxima
  if (any(changed))
    stop(context, " cannot preserve the fitted score structure for ",
         paste(names(maxima)[changed], collapse = ", "),
         "; a category or item lacks conditional information after the change. ",
         "Leave the change unresolved or rescore explicitly before refitting",
         call. = FALSE)
  invisible(TRUE)
}

.efrm_source_maxima <- function(fit) {
  vm <- fit$virtual_map
  if (is.null(vm) || !all(c("item", "vkey") %in% names(vm)))
    return(setNames(numeric(0), character(0)))
  mm <- fit$m[match(vm$vkey, colnames(fit$X))]
  out <- tapply(mm, vm$item, function(z) max(z, na.rm = TRUE))
  out[is.finite(out)]
}

# Recover the unexpanded item responses through the explicit virtual map.
# Virtual names are deliberately not parsed: item and level labels may contain
# colons, regular-expression metacharacters, or other valid punctuation.
.efrm_source_matrix <- function(fit, items = names(fit$set_of)) {
  vm <- fit$virtual_map
  if (is.null(vm) || !all(c("item", "vkey") %in% names(vm)))
    stop("the fitted frame-to-item map is not available")
  bad <- setdiff(items, unique(vm$item))
  if (length(bad)) stop("item(s) absent from the frame-to-item map: ",
                        paste(bad, collapse = ", "))
  src <- vapply(items, function(it) {
    cols <- vm$vkey[vm$item == it & vm$vkey %in% colnames(fit$X)]
    v <- rep(NA_real_, nrow(fit$X))
    for (cc in cols) {
      take <- !is.na(fit$X[, cc])
      v[take] <- fit$X[take, cc]
    }
    v
  }, numeric(nrow(fit$X)))
  if (!is.matrix(src)) src <- matrix(src, ncol = length(items))
  colnames(src) <- items
  src
}

# Response columns and person factors occupy separate public namespaces for
# matrix fits, so both may legitimately carry the same name. Structural refits
# temporarily assemble them in one data frame; give only the factor copy a
# private replay name, then restore the public name on the returned fit. This
# keeps role identity without renaming an item or changing any factor values.
.structural_factor_aliases <- function(factors, reserved) {
  if (is.null(factors)) return(list(
    data = NULL, original = character(0), replay = character(0)))
  factors <- as.data.frame(factors, check.names = FALSE,
                           stringsAsFactors = FALSE)
  original <- names(factors)
  replay <- original
  taken <- unique(c(reserved, original))
  for (j in which(replay %in% reserved)) {
    # The trailing delimiter prevents numeric prefix collisions (1 versus 11).
    # Also keep the complete private token out of every public factor name, so
    # later factorial-label restoration cannot rewrite an unchanged label that
    # merely contains the generated alias as a substring.
    candidate <- paste0(".rasch_refit_factor_", j, "__")
    while (candidate %in% taken ||
           any(vapply(original, function(x)
             grepl(candidate, x, fixed = TRUE), logical(1))))
      candidate <- paste0(candidate, ".")
    replay[j] <- candidate
    taken <- c(taken, candidate)
  }
  names(factors) <- replay
  list(data = factors, original = original, replay = replay)
}

.restore_structural_factor_names <- function(fit, aliases,
                                             frame_names = NULL) {
  replay <- aliases$replay
  original <- aliases$original
  if (!is.null(frame_names) && !is.null(fit$frame_group) &&
      length(frame_names) == length(fit$frame_group)) {
    # With crossed frames rasch_efrm() generates the primary cell name from
    # the component selectors. Include that generated name in the restoration
    # map, not only the component aliases supplied directly.
    replay <- c(replay, fit$frame_group)
    original <- c(original, frame_names)
  }
  restore <- function(x) {
    if (is.null(x) || !length(replay)) return(x)
    pos <- match(names(x), replay)
    hit <- !is.na(pos)
    names(x)[hit] <- original[pos[hit]]
    x
  }
  fit$factors <- restore(fit$factors)
  fit$person <- restore(fit$person)
  if (!is.null(frame_names)) {
    # Factorial tables are labelled from the replay component names. Translate
    # any private aliases there as well; estimates and covariance are unchanged.
    changed <- which(aliases$replay != aliases$original)
    display <- vapply(aliases$original[changed], function(x) {
      if (grepl("[:`]", x))
        paste0("`", gsub("`", "``", x, fixed = TRUE), "`") else x
    }, character(1))
    for (table in c("phi_factorial", "phi_factorial_tests")) {
      if (length(changed) && is.data.frame(fit[[table]]) &&
          "term" %in% names(fit[[table]])) {
        # Protect longer aliases first (.rasch_refit_factor_11 contains the
        # prefix .rasch_refit_factor_1), and use placeholders so an original
        # display name cannot itself be rewritten by a later substitution.
        ord <- changed[order(nchar(aliases$replay[changed]), decreasing = TRUE)]
        placeholder <- paste0("\034rasch_factor_", seq_along(ord), "\035")
        term <- fit[[table]]$term
        for (j in seq_along(ord))
          term <- gsub(aliases$replay[ord[j]], placeholder[j], term,
                       fixed = TRUE)
        display_ord <- display[match(ord, changed)]
        for (j in seq_along(ord))
          term <- gsub(placeholder[j], display_ord[j], term, fixed = TRUE)
        fit[[table]]$term <- term
      }
    }
  }
  if (!is.null(frame_names)) fit$frame_group <- frame_names
  fit
}

.efrm_refit <- function(fit, source, set_of, boot_reps = NULL,
                        ids = fit$person$id, factors = fit$factors,
                        se_method = NULL, score_max = NULL) {
  spec <- fit$refit_spec
  if (is.null(spec)) spec <- list()
  if (is.null(score_max)) {
    original_max <- .efrm_source_maxima(fit)
    score_max <- original_max[intersect(names(original_max), colnames(source))]
  }
  .require_score_structure(source, score_max, "the EFRM structural refit")
  group_vars <- spec$groups
  if (is.null(group_vars) || !all(group_vars %in% names(factors)))
    group_vars <- .frame_group_vars(fit)
  extra <- spec$factors
  if (is.null(extra)) extra <- setdiff(names(factors), fit$frame_group)
  extra <- intersect(extra, names(factors))
  need <- unique(c(group_vars, extra))
  taken <- unique(c(colnames(source), names(factors)))
  id_name <- ".rasch_id"
  while (id_name %in% taken) id_name <- paste0(id_name, ".")
  aliases <- .structural_factor_aliases(
    factors[, need, drop = FALSE], c(id_name, colnames(source)))
  replay_name <- stats::setNames(aliases$replay, aliases$original)
  group_replay <- unname(replay_name[group_vars])
  extra_replay <- unname(replay_name[extra])
  d <- data.frame(ids, source, aliases$data,
                  check.names = FALSE, stringsAsFactors = FALSE)
  names(d)[1L] <- id_name
  reps <- if (is.null(boot_reps)) spec$boot_reps else boot_reps
  if (is.null(reps)) {
    # Reproduce what was requested, not the smaller number that happened to
    # survive a bootstrap. The latter is outcome-dependent and can fall below
    # the public uncertainty threshold in a later structural refit.
    reps <- fit$boot_reps_requested %||% fit$boot_reps_used
    if (is.null(reps) || !is.finite(reps))
      reps <- if (any(is.finite(fit$alpha_table$se_log_alpha))) NULL else 0L
  }
  out <- do.call(rasch_efrm, list(
    data = d, item_sets = split(names(set_of), unname(set_of)),
    groups = group_replay, id = id_name,
    factors = if (length(extra_replay)) extra_replay else NULL,
    items = colnames(source),
    n_groups = .refit_n_groups(fit),
    # Source scores are reconstructed from the fitted response matrix; raw
    # missing codes have already been removed, and may now be valid scores.
    na_codes = integer(0),
    maxit = spec$maxit %||% 50, tol = spec$tol %||% 1e-7,
    min_link_persons = spec$min_link_persons %||% 30,
    se_method = se_method %||% spec$se_method %||% fit$se_method,
    boot_reps = reps, workers = spec$workers %||% 1L,
    seed = spec$seed %||% NULL))
  out <- .restore_structural_factor_names(out, aliases, fit$frame_group)
  # The fit correctly records the temporary selectors that produced it.
  # Replace those private names before another structural operation replays the
  # result, including the component selectors of a crossed frame.
  out$refit_spec$groups <- group_vars
  out$refit_spec$factors <- extra
  # A global maximum cannot detect category collapse confined to one frame.
  # Reconstruct the scored source so every observed response must retain its
  # score and frame-specific missingness after refitting.
  actual <- .efrm_source_matrix(out, colnames(source))
  if (!isTRUE(all.equal(unname(actual), unname(as.matrix(source)),
                       check.attributes = FALSE, tolerance = 0)))
    stop("the EFRM structural refit cannot preserve the fitted score structure; ",
         "a category or item lacks conditional information within a frame. ",
         "Leave the change unresolved or rescore explicitly before refitting",
         call. = FALSE)
  out
}

.rasch_refit <- function(fit, source, model = NULL, key_extra = NULL,
                         require_anchor = TRUE,
                         score_max = NULL) {
  .require_refittable_calibration(fit)
  spec <- fit$refit_spec
  if (is.null(spec)) spec <- list()
  source <- as.data.frame(source, check.names = FALSE,
                          stringsAsFactors = FALSE)
  keep <- names(source)
  if (is.null(score_max)) {
    same <- intersect(keep, fit$items$item)
    score_max <- stats::setNames(fit$m[match(same, fit$items$item)], same)
  }
  .require_score_structure(source, score_max, "the structural refit")
  if (inherits(fit, "rasch_explanatory") && all(keep %in% colnames(fit$X))) {
    out <- .explanatory_refit_modified(
      fit, source,
      inherit = stats::setNames(keep, keep))
    .require_fitted_score_structure(out, score_max)
    return(out)
  }
  key <- spec$key
  if (!is.null(fit$mc) && !is.null(key)) {
    raw_items <- intersect(colnames(fit$mc$raw), keep)
    for (it in raw_items) source[[it]] <- fit$mc$raw[, it]
    if (is.data.frame(key)) {
      key <- key[as.character(key$item) %in% keep, , drop = FALSE]
    } else {
      key <- key[intersect(names(key), keep)]
    }
  }
  if (!is.null(key_extra)) {
    if (!is.data.frame(key) && !is.null(key))
      stop("internal refit cannot combine incompatible key formats")
    key <- rbind(key, key_extra)
  }
  if (is.data.frame(key) && !nrow(key)) key <- NULL
  if (!is.data.frame(key) && !is.null(key) && !length(key)) key <- NULL
  anchors <- spec$anchors
  if (!is.null(anchors)) {
    anchor_items <- as.character(anchors$item)
    anchors <- anchors[anchor_items %in% keep, , drop = FALSE]
    if (!nrow(anchors) && require_anchor)
      stop("dropping those items would remove every anchor and lose the ",
           "fitted scale origin; retain an anchor or refit explicitly")
    if (!nrow(anchors)) anchors <- NULL
  }
  aliases <- .structural_factor_aliases(fit$factors, names(source))
  out <- do.call(rasch, list(
    data = source, model = model %||% spec$model %||% fit$model,
    items = names(source),
    id = fit$person$id,
    factors = aliases$data, n_groups = .refit_n_groups(fit),
    anchors = anchors,
    # Both the prepared scores and retained raw keyed answers already use
    # NA for missing responses. Applying the source file's codes again can
    # delete a valid score introduced when its categories were renumbered.
    na_codes = integer(0), key = key,
    pc_components = spec$pc_components,
    maxit = spec$maxit %||% 60, tol = spec$tol %||% 1e-8))
  out <- .restore_structural_factor_names(out, aliases)
  .require_fitted_score_structure(out, score_max)
  out
}

.rasch_refit_after_drop <- function(fit, keep) {
  .rasch_refit(fit, fit$X[, keep, drop = FALSE])
}
# Item screening is a step in an analysis, not a preliminary to it: an item
# is judged by its behaviour in the fit, and judging it means refitting
# without it. For frame models the refit matters more than usual: an item that
# does not follow its set's unit changes both the within-frame calibration and
# the person-side link used to compare that set with another.
# ===========================================================================

#' Drop items and refit
#'
#' Removes named items and refits the analysis with the same model
#' specification.
#'
#' @details
#' The refit retains person identifiers and factors, class-interval settings,
#' optimisation controls, anchors, multiple-choice scoring and PCM component
#' constraints. An EFRM refit also retains the item-set and crossed-frame
#' design, linking controls and uncertainty method. The operation is refused
#' if it would remove an externally anchored item, empty an item set or leave
#' the model unidentified. To change the anchor set, refit explicitly; this
#' prevents an item-removal comparison from silently changing its scale.
#' Older anchored or component-constrained fits without their original refit
#' settings must first be refitted from the source data.
#' A refit is also refused if it drops or rescores a retained item because a
#' response category no longer contributes conditional information.
#'
#' Item removal changes both the item calibration and the person estimates.
#' For an EFRM it can also change the estimated frame units. Compare the
#' original and refitted results as a sensitivity analysis.
#'
#' @param fit A fitted object from \code{\link{rasch}} or
#'   \code{\link{rasch_efrm}}. Many-facet fits are refused: remove the
#'   item's rows from the long-format data and refit
#'   \code{\link{rasch_mfrm}} instead.
#' @param items Item names to remove.
#' @param boot_reps Bootstrap replicates for an EFRM refit. The default retains
#'   the fitted specification; a number overrides it.
#' @return A refitted object of the same class as \code{fit}, carrying a note
#'   recording which items were dropped. Split-item provenance is retained for
#'   the items that remain.
#' @seealso \code{\link{frame_invariance}} and
#'   \code{\link{resolve_frames}} for frame models;
#'   \code{\link{split_items}} and \code{\link{resolve_dif}} for DIF; and
#'   \code{\link{combine_items}} for response dependence.
#' @examples
#' d <- simulate_rasch(300, 8, seed = 1)
#' fit <- rasch(d, id = "id")
#' fit2 <- drop_items(fit, "I03")
#' nrow(fit2$items)
#' @export
drop_items <- function(fit, items, boot_reps = NULL) {
  if (!inherits(fit, "rasch"))
    stop("drop_items needs a fit from rasch() or rasch_efrm()")
  .require_refittable_calibration(fit)
  if (inherits(fit, "rasch_mfrm"))
    stop("remove the item's rows from the long-format data and refit ",
         "rasch_mfrm() instead")
  items <- .structural_item_names(items, "drop")

  if (inherits(fit, "rasch_efrm")) {
    all_items <- names(fit$set_of)
    bad <- setdiff(items, all_items)
    if (length(bad))
      stop("item(s) not in the fit: ", paste(bad, collapse = ", "),
           "; the fit holds: ", paste(utils::head(all_items, 8), collapse = ", "),
           if (length(all_items) > 8) ", ..." else "")
    keep <- setdiff(all_items, items)
    if (!length(keep)) stop("dropping those items would leave no items")
    sets_left <- fit$set_of[keep]
    gone <- setdiff(unique(fit$set_of), unique(sets_left))
    if (length(gone))
      stop("dropping those items would empty set(s): ",
           paste(gone, collapse = ", "),
           "; a set with no items cannot carry a unit")
    if (length(unique(sets_left)) < 2L && length(unique(fit$set_of)) > 1L)
      stop("dropping those items would leave a single item set, so no set ",
           "unit is identified; drop fewer items or refit with rasch()")

    src <- .efrm_source_matrix(fit, keep)
    refit <- .efrm_refit(fit, src, sets_left, boot_reps = boot_reps)
  } else {
    bad <- setdiff(items, colnames(fit$X))
    if (length(bad))
      stop("item(s) not in the fit: ", paste(bad, collapse = ", "))
    keep <- setdiff(colnames(fit$X), items)
    if (length(keep) < 2L)
      stop("dropping those items would leave fewer than two items")
    anchors <- (fit$refit_spec %||% list())$anchors
    dropped_anchors <- if (is.null(anchors)) character(0) else
      intersect(items, as.character(anchors$item))
    if (length(dropped_anchors))
      stop("externally anchored item(s) cannot be dropped: ",
           paste(dropped_anchors, collapse = ", "),
           "; refit explicitly with a revised anchor set so the change of ",
           "scale is deliberate", call. = FALSE)
    refit <- .rasch_refit_after_drop(fit, keep)
  }
  if (!isTRUE(refit$est$converged))
    stop("the reduced calibration did not converge; the dropped-item analysis is unavailable")
  if (!.efrm_link_converged(refit))
    stop("the reduced calibration's set-unit link did not converge; the dropped-item analysis is unavailable")
  if (length(fit$subtest_map)) {
    sk <- intersect(names(fit$subtest_map), refit$items$item)
    refit$subtest_map <- fit$subtest_map[sk]
    refit$subtest_binary <- fit$subtest_binary[sk]
  }
  if (!is.null(fit$split_map)) {
    old_split <- .split_source_map(fit)
    refit$split_map <- old_split[colnames(refit$X)]
  }
  refit$notes <- c(refit$notes,
                   sprintf("dropped item(s): %s", paste(items, collapse = ", ")))
  refit
}
