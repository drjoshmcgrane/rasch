# rasch :: tailored analysis for guessing
# ===========================================================================
# The tailored test-of-fit procedure for detecting and correcting
# guessing on multiple-choice proficiency items (Waller 1989 ARRG; Andrich,
# Marais and Humphry 2012; Andrich & Marais 2019 ch. 17). Where a person's
# probability of success is below the chance level (for example 0.25 with
# four options), a correct response carries more guessing than information;
# tailoring converts those responses to missing and re-estimates. The
# four-step procedure compares the initial analysis with the tailored one
# on a common origin: (3) the initial data re-analysed with the mean
# location of the anchor items fixed at its tailored value (average item
# anchoring shifts the free calibration onto the tailored origin), and (4)
# the initial data with every item anchored at its tailored difficulty,
# re-estimating only the persons. Guessing shows as difficult items
# becoming harder in the tailored analysis.
# ===========================================================================

.tailored_boot_rows <- function(id) {
  # rows with no identifier are not one person: clustering them together
  # would resample them as a single unit and understate the uncertainty
  id_text <- .role_text_values(id)
  unknown <- is.na(id_text) | !nzchar(id_text)
  id_text[unknown] <- NA_character_
  key <- match(id_text, unique(id_text[!unknown]))
  if (any(unknown))
    key[unknown] <- length(unique(id_text[!unknown])) + seq_len(sum(unknown))
  clusters <- split(seq_along(id), key)
  picked <- sample.int(length(clusters), length(clusters), replace = TRUE)
  parts <- clusters[picked]
  list(rows = unlist(parts, use.names = FALSE),
       id = rep(seq_along(parts), lengths(parts)))
}

# Two-sided sign-count probability followed by Holm over m items. The strict
# decision rule is p_adj < .05, so a floor equal to .05 is not resolvable.
.tailored_boot_floor <- function(B, m) 2 * m / (B + 1)

.tailored_nonconvergence <- function(message) {
  stop(structure(list(message = message, call = NULL),
    class = c("rasch_tailored_nonconvergence", "error", "condition")))
}

# No tailoring in a resample is a valid zero change. Keep the complete item
# family before admitting that draw, since a dropped item is a different fit.
.tailored_boot_refit <- function(fit, chance, anchor_items, items, maxima) {
  if (!identical(as.character(fit$items$item), as.character(items)) ||
      !identical(as.integer(fit$m), as.integer(maxima)))
    return(.fit_boot_failure("error"))
  cut <- !is.na(fit$X) & is.finite(fit$moments$E) &
    fit$moments$E < chance
  if (!any(cut)) return(rep(0, length(items)))
  tryCatch({
    result <- tailored_analysis(fit, chance = chance,
      anchor_items = anchor_items, se_method = "none")
    shift <- .tailored_boot_shift(result, items)
    if (is.null(shift)) .fit_boot_failure("error") else shift
  }, rasch_tailored_nonconvergence = function(e)
    .fit_boot_failure("nonconverged"),
    error = function(e) .fit_boot_failure("error"))
}

#' Tailored analysis for guessing
#'
#' Runs the four-step tailored procedure of Andrich, Marais and Humphry
#' (2012) on a dichotomous analysis.
#' Step 1 is the supplied fit. Step 2 (tailored) sets to missing every
#' observed response whose modelled probability of success, at the step-1
#' person and item estimates, is below \code{chance}, and re-estimates
#' items and persons. Step 3 (origin-equated) re-analyses the
#' \emph{original} data with the mean location of the anchor items fixed at
#' its tailored value by average item anchoring (see \code{\link{pcml}}):
#' every item keeps its initial position relative to the others and the
#' calibration as a whole moves onto the tailored origin, so the two
#' calibrations can be compared item by item. Step 4 (all-anchored) fixes
#' every item at its tailored difficulty and re-estimates persons on the
#' original data. Guessing is
#' indicated when difficult items are estimated harder in the tailored
#' analysis than in the origin-equated one; the comparison table and
#' \code{\link{plot_equate}} on the two calibrations show it directly.
#'
#' @param fit An unanchored, unconstrained dichotomous fit from
#'   \code{\link{rasch}}. The procedure estimates its own common origin.
#' @param chance The guessing floor: the probability of success by chance
#'   (1/number of options; default 0.25).
#' @param anchor_items Items whose mean location fixes the common origin in
#'   step 3. The default takes the third of the test (at least two items)
#'   least affected by tailoring -- fewest responses removed, ties broken
#'   towards the easier tailored location -- which are the easy items the
#'   procedure trusts.
#' @param se_method \code{"none"} (default) reports the item shifts
#'   descriptively. \code{"bootstrap"} resamples persons and repeats the
#'   complete four-step procedure, including automatic anchor selection, to
#'   obtain standard errors, percentile intervals, and Holm-adjusted tests.
#'   When a person identifier occurs on several rows, all of that person's
#'   rows are resampled together. A resample requiring no tailoring contributes
#'   zero item shifts.
#' @param boot_reps Person-bootstrap replicates when
#'   \code{se_method = "bootstrap"}; at least 50, default 999. The
#'   sign-count bootstrap p-value has resolution floor \code{2/(boot_reps
#'   + 1)}, so after the Holm adjustment across m items the smallest
#'   achievable adjusted p is \code{2m/(boot_reps + 1)}; a warning fires
#'   when that floor is at or above 0.05 (the procedure declares significance
#'   only below 0.05, so detection would be impossible).
#' @param seed Optional non-negative whole-number seed for the person
#'   bootstrap. The caller's random-number state is restored on exit; see
#'   \code{\link{rasch_rng}} for generator support.
#' @return A list of class \code{"rasch_tailored"}: \code{tailored},
#'   \code{origin_equated}, and \code{anchored} fits, the comparison
#'   \code{table} (initial, tailored, origin-equated locations, the
#'   tailored-minus-equated \code{shift}; bootstrap uncertainty columns when
#'   requested), the number of
#'   responses removed, the anchor items used, \code{se_method}, and bootstrap
#'   accounting: requested, usable, non-converged, other failures, and the
#'   minimum usable count. \code{anchor_items_requested} distinguishes anchors
#'   supplied by the analyst from automatic anchor selection; it is
#'   \code{NULL} for the latter. The algorithm identifier and fitted-model and
#'   result signatures authenticate a saved result against the calibration
#'   and procedure from which it was computed.
#'   The final \code{anchored} component is a fixed-calibration scoring fit.
#'   Its person estimates and observed diagnostics remain available, but
#'   downstream item changes and refit-based bootstraps are not supported.
#'   Returned fits retain keyed scoring and structural records. Raw option
#'   data in the tailored fit exclude the responses removed by tailoring.
#'   For item-shift uncertainty, use this function's person bootstrap on the
#'   original calibration.
#' @references Waller, M. I. (1989). Modeling guessing behavior: A
#'   comparison of two IRT models. Applied Psychological Measurement, 13,
#'   233-243. Andrich, D., Marais, I. and Humphry, S. (2012). Using a
#'   theorem by Andersen and the dichotomous Rasch model to assess the
#'   presence of random guessing in multiple choice items. Journal of
#'   Educational and Behavioral Statistics, 37, 417-442.
#' @examples
#' set.seed(1); N <- 800
#' d <- seq(-2, 2.5, length.out = 10); th <- rnorm(N)
#' P <- plogis(outer(th, d, "-"))
#' P <- 0.25 + 0.75 * P            # uniform guessing floor
#' X <- matrix(rbinom(N * 10, 1, P), N, 10)
#' colnames(X) <- paste0("I", 1:10)
#' ta <- tailored_analysis(rasch(X), chance = 0.25)
#' ta$table
#' @export
tailored_analysis <- function(fit, chance = 0.25, anchor_items = NULL,
                              se_method = c("none", "bootstrap"),
                              boot_reps = 999L, seed = NULL) {
  se_method <- match.arg(se_method)
  if (!inherits(fit, "rasch")) stop("tailored_analysis needs a rasch fit")
  .require_refittable_calibration(fit)
  if (inherits(fit, "rasch_efrm") || inherits(fit, "rasch_mfrm"))
    stop("tailored_analysis currently requires an ordinary dichotomous Rasch fit")
  if (max(fit$m) > 1L)
    stop("tailored analysis applies to dichotomous (multiple-choice) items")
  if (!isTRUE(fit$est$converged))
    stop("the fitted calibration did not converge; tailored analysis is unavailable")
  if (inherits(fit, "rasch_explanatory"))
    stop("tailored_analysis() compares a free calibration with its tailored ",
         "recalibration; an explanatory fit restricts the item locations to ",
         "its design, and the tailored refit would drop that restriction, so ",
         "the two calibrations would not differ only by the tailoring")
  spec <- fit$refit_spec
  if (is.null(spec)) spec <- list()
  if (!is.null(spec$anchors) && nrow(spec$anchors))
    stop("tailored_analysis() requires an unanchored calibration because it estimates its own common origin")
  if (!is.null(spec$pc_components))
    stop("tailored_analysis() is not defined for principal-component constrained thresholds")
  .check_prob(chance, "chance")
  if (!is.null(seed)) seed <- .check_whole(seed, "seed", 0)
  if (!is.null(anchor_items)) {
    if (is.factor(anchor_items)) anchor_items <- as.character(anchor_items)
    if (!is.character(anchor_items) || !is.null(dim(anchor_items)) ||
        !length(anchor_items) || anyNA(anchor_items) ||
        any(!nzchar(trimws(anchor_items))))
      stop("`anchor_items` must be an ordinary vector of non-missing item names",
           call. = FALSE)
    if (anyDuplicated(anchor_items))
      stop("anchor item(s) named more than once: ",
           paste(unique(anchor_items[duplicated(anchor_items)]), collapse = ", "),
           call. = FALSE)
  }

  anchor_requested <- anchor_items
  # step 2: remove responses where the model gives success below chance
  P <- fit$moments$E                     # dichotomous: E = P(x = 1)
  Xt <- fit$X
  cut_cells <- !is.na(Xt) & !is.na(P) & P < chance
  if (!any(cut_cells))
    stop("no responses fall below the chance level; nothing to tailor")
  Xt[cut_cells] <- NA
  tailored <- rasch(Xt, model = fit$model, id = fit$person$id,
                    factors = fit$factors, n_groups = .refit_n_groups(fit),
                    maxit = spec$maxit %||% 60, tol = spec$tol %||% 1e-8)
  if (!isTRUE(tailored$est$converged))
    .tailored_nonconvergence("the tailored calibration did not converge; the comparison is unavailable")
  if (!identical(tailored$items$item, fit$items$item))
    stop("tailoring removed an item entirely; lower 'chance' or drop the item first")

  # step 3: common origin via average anchoring on the easy items -- those
  # least affected by tailoring (fewest responses removed, ties broken by
  # the lower tailored location), at least two and about a third of the test
  removed_per_item <- colSums(cut_cells)
  if (is.null(anchor_items)) {
    n_anchor <- max(2L, ceiling(ncol(fit$X) / 3))
    ord <- order(removed_per_item,
                 tailored$items$location[match(colnames(fit$X),
                                               tailored$items$item)])
    anchor_items <- colnames(fit$X)[ord[seq_len(n_anchor)]]
  }
  unknown <- setdiff(anchor_items, fit$items$item)
  if (length(unknown))
    stop("anchor item(s) not in the fit: ", paste(unknown, collapse = ", "))
  if (length(anchor_items) < 2)
    stop("need at least two anchor items for the common origin; ",
         "nominate easy items via anchor_items")
  ta_loc <- tailored$items$location[match(anchor_items, tailored$items$item)]
  # RUMM's average item anchoring: the initial calibration is estimated free
  # and shifted so the anchor items' mean location equals their tailored
  # mean. Fixing each anchor individually would instead hold the anchors at
  # their tailored values (a zero shift by construction) and estimate the
  # other items against them
  a3 <- data.frame(item = anchor_items, k = NA, tau = ta_loc, average = TRUE)
  origin_equated <- rasch(fit$X, model = fit$model, id = fit$person$id,
                          factors = fit$factors, n_groups = .refit_n_groups(fit),
                          anchors = a3, maxit = spec$maxit %||% 60,
                          tol = spec$tol %||% 1e-8)
  if (!isTRUE(origin_equated$est$converged))
    .tailored_nonconvergence("the common-origin calibration did not converge; the comparison is unavailable")

  # step 4: original data, every item fixed at its tailored value, persons
  # free. With no free item parameter there is nothing for pcml() to do, so
  # the fit is assembled directly on the anchored thresholds.
  thr_t <- tailored$thresholds
  thr4 <- thr_t; thr4$se <- 0; thr4$anchored <- TRUE
  est4 <- list(model = fit$model, thr = thr4,
               cov_tau = matrix(0, nrow(thr4), nrow(thr4)),
               loglik = NA_real_, iterations = 0L, converged = TRUE,
               m = fit$m, anchors = thr4, n_parameters = 0L)
  # step 4 is compared with steps 1-3, so it carries the same reference
  # sample size: hard-coding NA left its item-trait statistics on a
  # different scale from its own siblings
  anchored <- .assemble_fit(fit$model, fit$X, est4, fit$person$id,
                            fit$factors, .refit_n_groups(fit),
                            c(fit$notes,
                              "all item parameters anchored at their tailored values; persons re-estimated"))
  anchored$refit_spec <- list(
    model = fit$model, n_groups = .refit_n_groups(fit),
    anchors = data.frame(item = colnames(fit$X)[thr4$item],
                         k = thr4$k, tau = thr4$tau),
    fixed_calibration = TRUE, calibration_source = "tailored",
    na_codes = -1, key = NULL, pc_components = NULL,
    maxit = spec$maxit %||% 60L, tol = spec$tol %||% 1e-8)

  inherit_records <- function(out) {
    out$split_map <- fit$split_map
    out$subtest_map <- fit$subtest_map
    out$subtest_binary <- fit$subtest_binary
    out$mc <- fit$mc
    out$refit_spec$key <- spec$key
    if (!is.null(out$mc)) {
      keyed_items <- colnames(out$mc$raw)
      if (nrow(out$mc$raw) != nrow(out$X) ||
          !all(keyed_items %in% colnames(out$X)))
        stop("the keyed responses cannot be aligned with the tailored fit")
      # Keep the option identities only where this fit retains the response.
      # Copying the unmasked raw matrix would reinstate censored takers in
      # distractor summaries and later keyed refits.
      out$mc$raw[is.na(out$X[, keyed_items, drop = FALSE])] <- NA_character_
    }
    out
  }
  tailored <- inherit_records(tailored)
  origin_equated <- inherit_records(origin_equated)
  anchored <- inherit_records(anchored)

  idx_t <- match(fit$items$item, tailored$items$item)
  idx_o <- match(fit$items$item, origin_equated$items$item)
  shift <- tailored$items$location[idx_t] - origin_equated$items$location[idx_o]
  tab <- data.frame(item = fit$items$item,
                    initial = fit$items$location,
                    tailored = tailored$items$location[idx_t],
                    origin_equated = origin_equated$items$location[idx_o],
                    removed = removed_per_item,
                    shift = shift, se = NA_real_, ci_low = NA_real_,
                    ci_high = NA_real_, p = NA_real_, p_adj = NA_real_,
                    significant = NA)
  boot_used <- NA_integer_
  if (se_method == "bootstrap") {
    boot_reps <- .check_whole(boot_reps, "boot_reps", 50)
    if (!is.null(seed)) {
      old_seed <- .sim_seed_capture()
      on.exit(.sim_seed_restore(old_seed), add = TRUE)
      set.seed(seed)
    }
    draws <- list()
    boot_nonconverged <- 0L
    boot_errors <- 0L
    # One identifier can occur on several rows in a stacked repeated design.
    # Those rows are one sampling unit: resampling them separately would
    # discard their within-person dependence and understate uncertainty.
    for (bb in seq_len(boot_reps)) {
      bs <- .tailored_boot_rows(fit$person$id)
      take <- bs$rows
      Xb <- fit$X[take, , drop = FALSE]
      fb <- tryCatch(suppressWarnings(rasch(
        Xb, model = fit$model, id = bs$id,
        factors = if (is.null(fit$factors)) NULL else
          fit$factors[take, , drop = FALSE], n_groups = .refit_n_groups(fit),
        maxit = spec$maxit %||% 60, tol = spec$tol %||% 1e-8)),
        error = function(e) NULL)
      if (is.null(fb)) {
        boot_errors <- boot_errors + 1L
        next
      }
      if (!isTRUE(fb$est$converged)) {
        boot_nonconverged <- boot_nonconverged + 1L
        next
      }
      shift_b <- suppressWarnings(.tailored_boot_refit(
        fb, chance, anchor_requested, tab$item, fit$m))
      status <- .fit_boot_status(shift_b)
      if (status == "ok")
        draws[[length(draws) + 1L]] <- shift_b
      else if (status == "nonconverged")
        boot_nonconverged <- boot_nonconverged + 1L
      else boot_errors <- boot_errors + 1L
    }
    minimum_usable <- .fit_min_boot_success(boot_reps)
    if (length(draws) < minimum_usable)
      .fit_boot_refuse(
        "only ", length(draws), " of ", boot_reps,
        " tailored bootstrap replicates were estimable; at least ",
        minimum_usable, " are required for inference (",
        boot_nonconverged, " did not converge and ", boot_errors,
        " otherwise failed)", B = boot_reps, B_used = length(draws),
        B_nonconverged = boot_nonconverged, B_errors = boot_errors)
    B <- do.call(rbind, draws); boot_used <- nrow(B)
    tab$se <- apply(B, 2, stats::sd)
    tab$ci_low <- apply(B, 2, stats::quantile, probs = 0.025, names = FALSE)
    tab$ci_high <- apply(B, 2, stats::quantile, probs = 0.975, names = FALSE)
    tab$p <- vapply(seq_len(ncol(B)), function(j) {
      2 * min((sum(B[, j] <= 0) + 1) / (nrow(B) + 1),
              (sum(B[, j] >= 0) + 1) / (nrow(B) + 1), 0.5)
    }, 0)
    tab$p_adj <- .p_adjust_family(tab$p, method = "holm")
    tab$significant <- tab$p_adj < 0.05
    # the sign-count bootstrap p has resolution floor 2/(B+1); after the
    # Holm step the smallest achievable adjusted p is 2m/(B+1). If that
    # floor sits at or above alpha, NO item could ever be flagged at these
    # settings -- an inert test must say so, not sit quietly
    floor_adj <- .tailored_boot_floor(nrow(B), ncol(B))
    if (floor_adj >= 0.05)
      warning(sprintf(paste0(
        "with %d usable replicates and %d items the smallest achievable ",
        "Holm-adjusted p is %.3f (>= 0.05): no shift can reach significance ",
        "at these settings; at least %d usable draws are needed, so request more if refits fail"),
        nrow(B), ncol(B), floor_adj, 40L * ncol(B)),
        call. = FALSE)
  }
  out <- list(tailored = tailored, origin_equated = origin_equated,
              anchored = anchored, table = tab,
              algorithm = "tailored-four-stage-2",
              n_removed = sum(cut_cells), chance = chance,
              anchor_items = anchor_items,
              anchor_items_requested = anchor_requested,
              se_method = se_method,
              seed = if (se_method == "bootstrap") seed else NULL,
              boot_reps = if (se_method == "bootstrap") boot_reps else NA_integer_,
              boot_reps_used = boot_used,
              boot_reps_nonconverged = if (se_method == "bootstrap")
                boot_nonconverged else NA_integer_,
              boot_reps_errors = if (se_method == "bootstrap")
                boot_errors else NA_integer_,
              boot_minimum_usable = if (se_method == "bootstrap")
                minimum_usable else NA_integer_,
              fit_signature = .fit_boot_signature(fit))
  out <- .tag_tables(out)
  out$result_signature <- .fit_boot_md5(out)
  class(out) <- c("rasch_tailored", "list")
  out
}

# A completed refit is a usable bootstrap replicate only when it yields the
# complete finite shift vector. Letting NA/Inf enter the draw matrix makes its
# standard errors, percentile intervals and sign-count probabilities NA while
# falsely recording the replicate as successful.
.tailored_boot_shift <- function(result, items) {
  if (is.null(result) || !is.list(result) || !is.data.frame(result$table) ||
      !all(c("item", "shift") %in% names(result$table)) ||
      anyDuplicated(result$table$item) ||
      !setequal(as.character(result$table$item), as.character(items)))
    return(NULL)
  shift <- result$table$shift[match(as.character(items),
                                    as.character(result$table$item))]
  if (!is.numeric(shift) || length(shift) != length(items) ||
      any(!is.finite(shift))) return(NULL)
  unname(shift)
}

.validate_tailored_result <- function(result, fit) {
  if (is.null(result)) return(invisible(NULL))
  fail <- function() stop(
    "`tailored` is incomplete or internally inconsistent; recompute it with tailored_analysis()",
    call. = FALSE)
  required <- c(
    "tailored", "origin_equated", "anchored", "table", "n_removed",
    "algorithm", "chance", "anchor_items", "anchor_items_requested",
    "se_method", "boot_reps", "boot_reps_used",
    "boot_reps_nonconverged", "boot_reps_errors", "boot_minimum_usable",
    "fit_signature", "result_signature")
  if (!inherits(result, "rasch_tailored") ||
      !all(required %in% names(result)) ||
      !is.data.frame(result$table) ||
      !all(c("item", "initial", "tailored", "origin_equated", "removed",
             "shift", "se", "ci_low", "ci_high", "p", "p_adj",
             "significant") %in% names(result$table)) ||
      !is.character(result$result_signature) ||
      length(result$result_signature) != 1L || is.na(result$result_signature))
    fail()
  if (!identical(result$algorithm, "tailored-four-stage-2") &&
      !(identical(result$algorithm, "tailored-four-stage-1") &&
        identical(result$se_method, "none"))) fail()
  unsigned <- unclass(result)
  unsigned$result_signature <- NULL
  if (!.fit_boot_hash_matches(result$result_signature, unsigned)) fail()
  if (!inherits(fit, "rasch") ||
      !.fit_boot_signature_matches(result$fit_signature, fit))
    stop("`tailored` was computed from a different fitted model", call. = FALSE)
  if (!identical(result$se_method, "none") &&
      !identical(result$se_method, "bootstrap")) fail()
  numeric_fields <- c("initial", "tailored", "origin_equated", "removed",
                      "shift", "se", "ci_low", "ci_high", "p", "p_adj")
  if (any(!vapply(result$table[numeric_fields], is.numeric, logical(1))) ||
      !inherits(result$tailored, "rasch") ||
      !inherits(result$origin_equated, "rasch") ||
      !inherits(result$anchored, "rasch")) fail()
  if (!is.numeric(result$chance) || length(result$chance) != 1L ||
      !is.finite(result$chance) || result$chance <= 0 || result$chance >= 1 ||
      !is.numeric(result$n_removed) || length(result$n_removed) != 1L ||
      !is.finite(result$n_removed) || result$n_removed < 1 ||
      result$n_removed != floor(result$n_removed) ||
      !is.character(result$anchor_items) || length(result$anchor_items) < 2L ||
      anyNA(result$anchor_items) || any(!nzchar(trimws(result$anchor_items))) ||
      anyDuplicated(result$anchor_items) ||
      any(!result$anchor_items %in% fit$items$item) ||
      !identical(as.character(result$table$item), as.character(fit$items$item)) ||
      any(!is.finite(result$table$removed)) || any(result$table$removed < 0) ||
      any(result$table$removed != floor(result$table$removed)) ||
      sum(result$table$removed) != result$n_removed)
    fail()
  requested <- result$anchor_items_requested
  if (!is.null(requested) &&
      (!is.character(requested) || !length(requested) || anyNA(requested) ||
       any(!nzchar(trimws(requested))) || anyDuplicated(requested) ||
       !identical(requested, result$anchor_items))) fail()
  aligned_location <- function(x)
    x$items$location[match(result$table$item, x$items$item)]
  same <- function(x, y) isTRUE(all.equal(
    as.numeric(x), as.numeric(y), tolerance = 64 * .Machine$double.eps,
    check.attributes = FALSE))
  if (!same(result$table$initial, aligned_location(fit)) ||
      !same(result$table$tailored, aligned_location(result$tailored)) ||
      !same(result$table$origin_equated,
            aligned_location(result$origin_equated)) ||
      !same(result$table$shift,
            result$table$tailored - result$table$origin_equated)) fail()
  if (is.null(requested)) {
    n_anchor <- max(2L, ceiling(nrow(result$table) / 3))
    expected_anchor <- result$table$item[order(
      result$table$removed, result$table$tailored)][seq_len(n_anchor)]
    if (!identical(result$anchor_items, expected_anchor)) fail()
  }
  count_names <- c("boot_reps", "boot_reps_used", "boot_reps_nonconverged",
                   "boot_reps_errors", "boot_minimum_usable")
  whole <- function(x) is.numeric(x) && length(x) == 1L && is.finite(x) &&
    x >= 0 && x == floor(x)
  if (identical(result$se_method, "bootstrap")) {
    if (any(!vapply(result[count_names], whole, logical(1))) ||
        result$boot_reps < 50L ||
        result$boot_reps_used + result$boot_reps_nonconverged +
          result$boot_reps_errors != result$boot_reps ||
        result$boot_minimum_usable != .fit_min_boot_success(result$boot_reps) ||
        result$boot_reps_used < result$boot_minimum_usable ||
        any(!is.finite(unlist(result$table[c(
          "se", "ci_low", "ci_high", "p", "p_adj")], use.names = FALSE))) ||
        any(result$table$p < 0 | result$table$p > 1) ||
        any(result$table$p_adj < 0 | result$table$p_adj > 1) ||
        any(result$table$ci_low > result$table$ci_high) ||
        !same(result$table$p_adj,
              .p_adjust_family(result$table$p, method = "holm")) ||
        anyNA(result$table$significant) ||
        !identical(as.logical(result$table$significant),
                   result$table$p_adj < 0.05)) fail()
  } else {
    if (any(!vapply(result[count_names], function(x)
      is.numeric(x) && length(x) == 1L && is.na(x), logical(1))) ||
        any(!is.na(unlist(result$table[c(
          "se", "ci_low", "ci_high", "p", "p_adj")], use.names = FALSE))) ||
        any(!is.na(result$table$significant))) fail()
  }
  invisible(result)
}

#' @export
print.rasch_tailored <- function(x, ...) {
  cat(sprintf("Tailored analysis: %d response(s) below chance %.2f set to missing\n",
              x$n_removed, x$chance))
  cat(sprintf("Origin from average-anchored items: %s\n",
              paste(x$anchor_items, collapse = ", ")))
  tab <- x$table
  num <- vapply(tab, is.numeric, TRUE)
  tab[num] <- lapply(tab[num], round, 3)
  print(tab, row.names = FALSE)
  if (x$se_method == "bootstrap") {
    up <- sum(x$table$significant %in% TRUE & x$table$shift > 0)
    cat(sprintf(paste0("%d item(s) significantly harder after the ",
                       "Holm-adjusted person bootstrap (%d usable draws).\n"),
                up, x$boot_reps_used))
  } else cat("Item shifts are descriptive; use se_method = 'bootstrap' for inference.\n")
  invisible(x)
}
