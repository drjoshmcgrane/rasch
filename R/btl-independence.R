# rasch :: paired-comparison independence and dimensionality
#
# The Rasch residual-correlation (Q3) and residual-PCA displays operate on the
# persons x items residual matrix, which paired-comparison data do not
# produce (the datum is a comparison, so residuals live on object pairs). The
# analogues that DO fall out of the pair structure are transitivity (are the
# preferences consistent with a single order?) and a decomposition of the
# object-by-object residual preference matrix (is there structured departure
# from the one-dimensional fit, i.e. a second attribute steering contests?).

# weighted "points" object i scored against object j, as a K x K matrix
# (S[i, j]); resp is the per-comparison response (0..m for a, m - x for b).
.btl_scores <- function(ia, ib, resp, w, m, K) {
  S <- matrix(0, K, K)
  fill <- function(rows, cols, val) {
    idx <- (cols - 1L) * K + rows                    # column-major linear index
    ag <- rowsum(val, idx)
    S[as.integer(rownames(ag))] <<- S[as.integer(rownames(ag))] + ag[, 1]
  }
  fill(ia, ib, w * resp)
  fill(ib, ia, w * (m - resp))
  S
}

# count transitive vs circular triads in a tournament: D[i,j] = +1 if i beats
# j, -1 if j beats i, 0 tie; seen[i,j] whether the pair was compared. Returns
# the triple counts and, per object, how many circular triads it sits in.
.btl_triads <- function(D, seen) {
  K <- nrow(D)
  n_tri <- 0L; n_circ <- 0L
  invol <- integer(K)
  if (K >= 3L) for (i in 1:(K - 2L)) for (j in (i + 1L):(K - 1L)) {
    if (!seen[i, j] || D[i, j] == 0) next
    for (k in (j + 1L):K) {
      if (!(seen[i, k] && seen[j, k]) || D[i, k] == 0 || D[j, k] == 0) next
      n_tri <- n_tri + 1L
      # a triple's win-counts are {2,1,0} when transitive, {1,1,1} when circular
      si <- (D[i, j] > 0) + (D[i, k] > 0)
      sj <- (D[i, j] < 0) + (D[j, k] > 0)
      sk <- (D[i, k] < 0) + (D[j, k] < 0)
      if (si == 1L && sj == 1L && sk == 1L) {
        n_circ <- n_circ + 1L
        invol[i] <- invol[i] + 1L; invol[j] <- invol[j] + 1L
        invol[k] <- invol[k] + 1L
      }
    }
  }
  list(n_triples = n_tri, n_circular = n_circ, involvement = invol)
}

#' Transitivity of paired comparisons
#'
#' Summarises circular triads in the observed paired comparisons. A triad is
#' circular when A is preferred to B, B to C, and C to A. For a complete
#' tournament, the function reports Kendall's coefficient of consistency
#' (Kendall and Babington Smith 1940). Judge-specific summaries are returned
#' when judges are available.
#'
#' A circular-triad rate of one quarter is the benchmark for a random
#' tournament. It is not the expected rate under a fitted BTL model with
#' unequal object locations, so this function is a descriptive consistency
#' measure rather than a calibrated goodness-of-fit test. Comparisons involving
#' an undefeated or winless object remain part of this observed-data summary,
#' although that object is set aside from finite maximum-likelihood estimation.
#'
#' @param fit A paired-comparison fit from \code{\link{btl}}.
#' @param min_triples A judge is reported only if this many complete triples
#'   (all three pairs judged) are available.
#' @return A list of class \code{"rasch_btl_transitivity"}: \code{summary} (one
#'   row: objects, pairs compared, complete triples, circular triads, the
#'   circular rate, the chance rate 0.25, the consistency index
#'   \code{1 - rate/0.25}, and Kendall's \code{zeta} when the design is a
#'   complete round-robin with no exactly-tied pair -- \code{NA} otherwise);
#'   \code{objects} (each object's circular-triad
#'   involvement); \code{judges} (per-judge consistency, when judges exist);
#'   and \code{notes}.
#' @references Kendall, M. G., & Babington Smith, B. (1940). On the method of
#'   paired comparisons. \emph{Biometrika}, 31(3/4), 324-345.
#' @examples
#' set.seed(1); objs <- LETTERS[1:6]; beta <- setNames(seq(-1.5, 1.5, len = 6), objs)
#' pr <- t(utils::combn(objs, 2))
#' d <- data.frame(a = rep(pr[, 1], each = 20), b = rep(pr[, 2], each = 20))
#' d$win <- ifelse(runif(nrow(d)) < plogis(beta[d$a] - beta[d$b]), d$a, d$b)
#' btl_transitivity(btl(d, "a", "b", "win"))
#' @export
btl_transitivity <- function(fit, min_triples = 5L) {
  if (!inherits(fit, "rasch_btl")) stop("not a paired-comparison (btl) fit")
  min_triples <- .check_whole(min_triples, "min_triples", 1)
  # Transitivity is a property of the observed tournament, not of the finite
  # calibration. New fits retain the otherwise usable rows that preceded
  # boundary-object removal; saved older fits fall back to their fitted rows.
  cmp <- fit$observed_comparisons %||% fit$comparisons
  objs <- sort(unique(c(as.character(cmp$object_a),
                        as.character(cmp$object_b))))
  K <- length(objs); m <- fit$m
  ia <- match(cmp$object_a, objs); ib <- match(cmp$object_b, objs)
  notes <- paste(
    "the 0.25 chance rate is a random-tournament benchmark, not the fitted",
    "BTL expected circular rate; transitivity is descriptive")

  tournament <- function(rows) {
    S <- .btl_scores(ia[rows], ib[rows], cmp$response[rows], cmp$weight[rows],
                     m, K)
    seen <- (S + t(S)) > 0
    D <- sign(S - t(S)); D[!seen] <- 0
    list(D = D, seen = seen)
  }

  tt <- tournament(seq_len(nrow(cmp)))
  tri <- .btl_triads(tt$D, tt$seen)
  n_pairs <- sum(tt$seen[upper.tri(tt$seen)])
  complete <- n_pairs == choose(K, 2)
  ties <- sum(tt$seen & tt$D == 0) / 2
  rate <- if (tri$n_triples) tri$n_circular / tri$n_triples else NA_real_
  # Kendall's coefficient of consistency, defined for a complete round-robin
  zeta <- NA_real_
  if (complete && ties == 0 && K >= 3) {
    dmax <- if (K %% 2L) (K^3 - K) / 24 else (K^3 - 4 * K) / 24
    zeta <- 1 - tri$n_circular / dmax
  }
  if (!complete)
    notes <- c(notes, sprintf(
      "%d of %d object pairs were compared; triads use complete triples only",
      n_pairs, choose(K, 2)))
  if (ties > 0)
    notes <- c(notes, sprintf("%d pair(s) split exactly evenly and are set aside",
                              ties))

  summary <- data.frame(
    n_objects = K, n_pairs = n_pairs, n_triples = tri$n_triples,
    n_circular = tri$n_circular, circular_rate = rate,
    chance_rate = 0.25, consistency = 1 - rate / 0.25, zeta = zeta)

  objects <- data.frame(object = objs, circular_triads = tri$involvement)
  objects <- objects[order(-objects$circular_triads), ]
  rownames(objects) <- NULL

  judges <- NULL
  if (!all(is.na(cmp$judge))) {
    ju <- sort(unique(cmp$judge[!is.na(cmp$judge)]))
    rows <- lapply(ju, function(j) which(cmp$judge == j))
    jt <- lapply(rows, function(rr) {
      tj <- tournament(rr); .btl_triads(tj$D, tj$seen)
    })
    n_tri <- vapply(jt, `[[`, 0L, "n_triples")
    keep <- n_tri >= min_triples
    if (any(keep)) {
      n_c <- vapply(jt, `[[`, 0L, "n_circular")[keep]
      nt <- n_tri[keep]
      jr <- n_c / nt
      judges <- data.frame(
        judge = ju[keep],
        n_comparisons = vapply(rows[keep], function(rr) sum(cmp$weight[rr]), 0),
        n_triples = nt, n_circular = n_c, circular_rate = jr,
        consistency = 1 - jr / 0.25)
      judges <- judges[order(judges$consistency), ]
      rownames(judges) <- NULL
    } else {
      notes <- c(notes, sprintf(
        "no judge reached %d complete triples; per-judge consistency omitted",
        min_triples))
    }
  }

  out <- list(summary = summary, objects = objects, judges = judges,
              notes = notes)
  out <- .tag_tables(out)
  class(out) <- "rasch_btl_transitivity"
  out
}

#' @export
print.rasch_btl_transitivity <- function(x, ...) {
  s <- x$summary
  cat(sprintf("Paired-comparison transitivity: %d objects, %d complete triples\n",
              s$n_objects, s$n_triples))
  cat(sprintf("Circular triads: %d (%.1f%% of triples; random-tournament benchmark %.0f%%) -> consistency %.2f\n",
              s$n_circular, 100 * s$circular_rate, 100 * s$chance_rate,
              s$consistency))
  if (!is.na(s$zeta))
    cat(sprintf("Kendall coefficient of consistency (complete design): %.3f\n",
                s$zeta))
  if (!is.null(x$judges))
    cat(sprintf("Per-judge consistency reported for %d judge(s); least consistent %.2f\n",
                nrow(x$judges), min(x$judges$consistency)))
  for (n in x$notes) cat("Note:", n, "\n")
  invisible(x)
}

# Pool observed and fitted expected points over the same oriented pairs.
# For ordered responses logit(E[X]/m) is not the object-location gap;
# fitted frame and history effects also change the pooled expectation.
# The same half-point correction on both sides keeps extreme cells finite
# without creating a residual when observed and expected points agree.
.btl_resid_matrix_expected <- function(ia, ib, resp, w, K, expected_a,
                                        m = 1L) {
  S <- .btl_scores(ia, ib, resp, w, m, K)
  E <- .btl_scores(ia, ib, expected_a, w, m, K)
  tot <- S + t(S)
  Pobs <- (S + .5) / (tot + 1)
  Pexp <- (E + .5) / (tot + 1)
  R <- stats::qlogis(Pobs) - stats::qlogis(Pexp)
  R[tot == 0] <- 0
  (R - t(R)) / 2
}

# real skew-symmetric R has eigenvalues in +/- i*lambda pairs; the positive
# lambda are the "bimension" strengths, each a plane of cyclic residual
# structure. Returns strengths (desc) and the leading plane's coordinates.
.btl_bimensions <- function(R) {
  e <- eigen(R)
  lam <- Im(e$values)
  keep <- which(lam > 1e-8)
  if (!length(keep)) return(list(strength = numeric(0), coord = NULL))
  ord <- keep[order(lam[keep], decreasing = TRUE)]
  v <- e$vectors[, ord[1]]
  list(strength = lam[ord], total = sum(R^2),
       coord = cbind(x = Re(v), y = Im(v)))
}

.btl_dimensionality_scope_note <- function(independent_comparisons) {
  if (independent_comparisons) paste(
    "the simulated reference assumes conditionally independent comparison",
    "outcomes given the fitted probabilities and any modeled history;",
    "it is not cluster-robust") else paste(
    "dimensionality inference is withheld because independent comparison",
    "outcomes were not assumed; judge-clustered calibration does not supply",
    "a cluster-robust simulated reference. Set independent_comparisons = TRUE",
    "only for a conditional independent-comparison sensitivity analysis")
}

.btl_dimensionality_efrm <- function(fit, reps, seed = NULL,
                                    independent_comparisons = FALSE) {
  objs <- fit$objects$object; K <- length(objs)
  if (K < 3L) stop("need at least three objects")
  cmp <- fit$comparisons
  ia <- match(cmp$object_a, objs); ib <- match(cmp$object_b, objs)
  w <- cmp$weight
  R <- .btl_resid_matrix_expected(ia, ib, cmp$response, w, K,
                                  cmp$expected)
  bm <- .btl_bimensions(R)
  if (!length(bm$strength)) stop("no residual structure to decompose")
  S_seen <- .btl_scores(ia, ib, cmp$response, w, 1L, K)
  n_seen <- sum(((S_seen + t(S_seen)) > 0)[upper.tri(diag(K))])
  complete_pairs <- n_seen == choose(K, 2)
  lead_ref <- vapply(seq_len(reps), function(r) {
    yr <- stats::rbinom(nrow(cmp), 1L, cmp$expected)
    rr <- .btl_resid_matrix_expected(ia, ib, yr, w, K, cmp$expected)
    s <- .btl_bimensions(rr)$strength
    if (length(s)) s[1] else 0
  }, 0)
  inference <- if (complete_pairs && independent_comparisons)
    .sim_upper_family(bm$strength[1], matrix(lead_ref, ncol = 1L)) else NULL
  ref_mean <- if (is.null(inference)) NA_real_ else inference$mean[1L]
  ref_p95 <- if (is.null(inference)) NA_real_ else inference$critical[1L]
  nb <- length(bm$strength)
  prop <- 2 * bm$strength^2 / bm$total
  lead_flag <- if (is.null(inference)) NA else inference$significant[1L]
  bimensions <- data.frame(
    bimension = seq_len(nb), strength = bm$strength,
    prop_residual = prop,
    ref_mean = c(ref_mean, rep(NA_real_, nb - 1L)),
    ref_p95 = c(ref_p95, rep(NA_real_, nb - 1L)),
    above_reference = c(lead_flag, rep(NA, nb - 1L)))
  coords <- data.frame(object = objs, location = fit$objects$location,
                       x = bm$coord[, "x"], y = bm$coord[, "y"])
  notes <- paste(
    "the simulated reference retains each comparison's fitted panel-by-set",
    "probability and the observed frame allocation")
  notes <- c(notes, .btl_dimensionality_scope_note(independent_comparisons))
  if (!complete_pairs)
    notes <- c(notes, sprintf(
      paste0("%d of %d pairs compared; unseen pairs contribute no residual ",
             "information, so dimensionality inference is withheld"),
      n_seen, choose(K, 2)))
  out <- list(
    bimensions = bimensions, coords = coords,
    leading_structured = lead_flag,
    reference = list(mean = ref_mean, p95 = ref_p95,
                     p = if (is.null(inference)) NA_real_ else inference$p[1L],
                     p_adj = if (is.null(inference)) NA_real_ else
                       inference$p_adjusted[1L], reps = reps,
                     n_used = length(lead_ref),
                     alpha = if (is.null(inference)) 0.05 else inference$alpha,
                     inference_available = !is.null(inference),
                     independent_comparisons = independent_comparisons,
                     draws = lead_ref,
                     method = if (is.null(inference)) NA_character_ else
                       inference$method, seed = seed),
    residual_matrix = R, residual_method = "pooled-expected-score-1",
    notes = notes)
  attr(out, "fit_signature") <- .fit_boot_signature(fit)
  class(out) <- "rasch_btl_dim"
  attr(out, "result_signature") <- .fit_boot_md5(out)
  out
}

.authenticate_btl_dimensionality <- function(result, fit) {
  if (is.null(result)) return(invisible(NULL))
  signature <- attr(result, "result_signature")
  unsigned <- result
  attr(unsigned, "result_signature") <- NULL
  if (!inherits(result, "rasch_btl_dim") || !is.list(result) ||
      !is.data.frame(result$bimensions) || !is.list(result$reference) ||
      !is.character(signature) || length(signature) != 1L || is.na(signature) ||
      !.fit_boot_hash_matches(signature, unsigned) ||
      !.fit_boot_signature_matches(attr(result, "fit_signature"), fit))
    stop("`dimensionality` must be a btl_dimensionality() result from this fitted model")
  invisible(result)
}

.validate_btl_dimensionality <- function(result, fit) {
  if (is.null(result)) return(invisible(NULL))
  .authenticate_btl_dimensionality(result, fit)
  .require_current_btl_residuals(result)
  if (.btl_dimensionality_unsupported_reference(result, fit))
    stop("`dimensionality` retains an inferential reference despite its ",
         "unsupported comparison design; recompute it with btl_dimensionality()")
  invisible(result)
}

.require_current_btl_residuals <- function(result) {
  if (!identical(result$residual_method, "pooled-expected-score-1"))
    stop("this dimensionality result uses an earlier residual definition; ",
         "recompute it with btl_dimensionality()", call. = FALSE)
  invisible(result)
}

.btl_dimensionality_has_inference <- function(x) {
  (is.null(x$reference$inference_available) ||
     isTRUE(x$reference$inference_available)) &&
    length(x$leading_structured) == 1L && !is.na(x$leading_structured)
}

.btl_dimensionality_unsupported_reference <- function(result, fit) {
  finite_reference <- any(is.finite(c(
    result$reference$mean, result$reference$p95, result$reference$p,
    result$reference$p_adj, result$bimensions$ref_mean,
    result$bimensions$ref_p95)))
  if (!finite_reference) return(FALSE)
  if (!.btl_dimensionality_has_inference(result)) return(TRUE)
  if ((isTRUE(fit$clustered) || inherits(fit, "rasch_btl_efrm")) &&
      !isTRUE(result$reference$independent_comparisons)) return(TRUE)
  # Earlier frame results could retain a TRUE/FALSE flag despite unseen
  # pairs. Check the actual retained-object design, not that legacy flag.
  tab <- fit$objects
  if (!inherits(fit, "rasch_btl_efrm") && "extreme" %in% names(tab))
    tab <- tab[!(tab$extreme %in% TRUE), ]
  objs <- tab$object
  cmp <- fit$comparisons
  S <- .btl_scores(match(cmp$object_a, objs), match(cmp$object_b, objs),
                   cmp$response, cmp$weight, fit$m, length(objs))
  seen <- sum(((S + t(S)) > 0)[upper.tri(S)])
  seen != choose(length(objs), 2L)
}

# Diagnose whether comparison position varies across judges. Repeated
# presentations of a pair within one judge are first reduced to that judge's
# mean position; otherwise within-judge repetition can be mistaken for the
# across-judge variation needed to separate order from object structure.
.btl_order_variation <- function(dd, objs) {
  ia <- match(dd$object_a, objs)
  ib <- match(dd$object_b, objs)
  pair <- .factor_keys(data.frame(lo = pmin(ia, ib), hi = pmax(ia, ib)))
  judge <- as.character(dd$judge)
  pos <- numeric(nrow(dd))
  for (ix in split(seq_len(nrow(dd)), judge)) {
    r <- rank(dd$order[ix], ties.method = "first")
    pos[ix] <- if (length(ix) > 1L) (r - 1) / (length(ix) - 1) else 0
  }
  key <- .factor_keys(data.frame(pair = pair, judge = judge))
  cells <- split(seq_len(nrow(dd)), key)
  by_judge_pair <- do.call(rbind, lapply(cells, function(ix) data.frame(
    pair = pair[ix[1L]], judge = judge[ix[1L]], position = mean(pos[ix]),
    repetitions = length(ix), stringsAsFactors = FALSE)))
  rownames(by_judge_pair) <- NULL
  judge_n <- table(judge)
  # Under a random permutation, the variance of the mean of r positions
  # sampled without replacement from a sequence of length n is the discrete
  # uniform variance divided by r, with its finite-population correction.
  # Scaling the observed between-judge variance by that expectation prevents
  # many repeats of one pair from looking like a shared order merely because
  # their mean position is naturally more stable than a single position.
  by_judge_pair$random_variance <- mapply(function(j, r) {
    n <- unname(judge_n[j])
    if (!is.finite(n) || n <= 1L || r >= n) return(0)
    ((n + 1) / (12 * (n - 1))) / r * ((n - r) / (n - 1))
  }, by_judge_pair$judge, by_judge_pair$repetitions)
  pair_var <- tapply(by_judge_pair$position, by_judge_pair$pair, stats::var)
  expected <- tapply(by_judge_pair$random_variance, by_judge_pair$pair, mean)
  relative <- pair_var / expected
  usable <- is.finite(relative) & expected > sqrt(.Machine$double.eps)
  if (!any(usable))
    return(list(shared = TRUE,
                replicated = any(table(by_judge_pair$pair) > 1L),
                variance = pair_var, expected = expected,
                relative_variance = relative))
  list(shared = mean(relative[usable]) < 0.25,
       replicated = TRUE, variance = pair_var, expected = expected,
       relative_variance = relative)
}

#' Experimental residual dimensionality of paired comparisons
#'
#' Decomposes the skew-symmetric matrix of observed-minus-expected pair
#' log-odds into Gower's (1977) rotational planes, or bimensions. A large
#' leading bimension indicates a structured cycle in the residual comparisons.
#' Its strength is compared with simulations from the fitted one-dimensional
#' model using the observed comparison counts.
#'
#' This is an experimental diagnostic. The reference is conditional on the
#' fitted point estimates because the model is not re-estimated in each
#' replicate. Observed and fitted expected points are pooled over each object
#' pair before their proportions are compared on the logit scale. A half-point
#' correction is applied to both totals. This retains category thresholds,
#' frame units and fitted position or history effects in the expectation.
#' The same calculation is used in the simulations; exposure and carry-over
#' expectations follow each draw's generated history. The fitted model must
#' have converged.
#'
#' This reference assumes conditionally independent comparison outcomes given
#' the fitted probabilities and any modeled history. It is not cluster-robust:
#' clustered calibration standard errors do not make this simulated reference
#' valid under general within-judge dependence. Judge-clustered and frame fits
#' therefore return a descriptive decomposition by default. Explicitly set
#' \code{independent_comparisons = TRUE} to request the model-conditional
#' independent-comparison reference as a sensitivity analysis.
#'
#' Inference is withheld if any object pair is unobserved, or if an
#' ordered analysis contains count-weighted rows whose within-row sequence is
#' unavailable. It is also withheld when every judge receives essentially the
#' same comparison sequence and an order effect is fitted, because order and
#' residual structure are then confounded. The result is also withheld when
#' no object pair has an observed position for more than one judge, because
#' across-judge order variation cannot then be assessed. In these cases the
#' observed decomposition remains available, but probabilities, critical
#' values and the reference band are omitted. Any completed simulation draws
#' are retained for descriptive inspection, not as an inferential reference.
#'
#' @param fit A paired-comparison fit from \code{\link{btl}}.
#' @param reps Model-simulated replicates for the noise reference; at least 20.
#'   Larger values give a more stable upper-tail reference.
#' @param seed Optional non-negative whole-number seed. The caller's random-
#'   number state is restored when the calculation finishes; see
#'   \code{\link{rasch_rng}} for generator support.
#' @param independent_comparisons Logical or \code{NULL}. Declare whether
#'   conditionally independent comparison outcomes may be assumed for the
#'   simulated reference. The default \code{NULL} assumes independence only
#'   for an unclustered ordinary BTL fit. \code{FALSE} withholds inference;
#'   \code{TRUE} explicitly enables the reference, subject to design guards.
#' @return A list of class \code{"rasch_btl_dim"}: \code{bimensions} (per
#'   bimension: strength and share of residual size; the reference mean, 5%
#'   upper critical value, and the clears-the-reference flag are reported for the
#'   leading bimension and \code{NA} for the rest);
#'   \code{coords} (each object's position in the leading bimension plane, for
#'   the residual map); \code{leading_structured} (whether bimension 1 clears
#'   its reference); \code{reference} (the simulated mean, finite-simulation
#'   5% critical value and upper-tail probability, calculated as one plus the
#'   exceedance count divided by one plus \code{reps}, together with the
#'   effective \code{independent_comparisons} assumption); \code{residual_matrix}; and
#'   \code{notes}. \code{residual_method} identifies the residual definition.
#' @references Gower, J. C. (1977). The analysis of asymmetry and orthogonality.
#'   In J. R. Barra et al. (Eds.), \emph{Recent Developments in Statistics}
#'   (pp. 109-123). North-Holland.
#' @examples
#' set.seed(1); objs <- LETTERS[1:6]; beta <- setNames(seq(-1.5, 1.5, len = 6), objs)
#' pr <- t(utils::combn(objs, 2))
#' d <- data.frame(a = rep(pr[, 1], each = 30), b = rep(pr[, 2], each = 30))
#' d$win <- ifelse(runif(nrow(d)) < plogis(beta[d$a] - beta[d$b]), d$a, d$b)
#' btl_dimensionality(btl(d, "a", "b", "win"), reps = 20)
#' @export
btl_dimensionality <- function(fit, reps = 200L, seed = NULL,
                               independent_comparisons = NULL) {
  if (!inherits(fit, "rasch_btl")) stop("not a paired-comparison (btl) fit")
  if (!isTRUE(fit$converged))
    stop("the paired-comparison calibration did not converge; dimensionality inference is unavailable")
  reps <- .check_whole(reps, "reps", 20)
  if (!is.null(independent_comparisons) &&
      (!is.logical(independent_comparisons) ||
       length(independent_comparisons) != 1L ||
       is.na(independent_comparisons) || is.object(independent_comparisons) ||
       !is.null(dim(independent_comparisons))))
    stop("`independent_comparisons` must be TRUE, FALSE, or NULL")
  if (is.null(independent_comparisons))
    independent_comparisons <- !isTRUE(fit$clustered) &&
      !inherits(fit, "rasch_btl_efrm")
  if (!is.null(seed)) {
    seed <- .check_whole(seed, "seed", 0)
    old_seed <- .sim_seed_capture()
    on.exit(.sim_seed_restore(old_seed), add = TRUE)
    set.seed(seed)
  }
  if (inherits(fit, "rasch_btl_efrm"))
    return(.btl_dimensionality_efrm(fit, reps, seed = seed,
      independent_comparisons = independent_comparisons))
  tab <- fit$objects
  if ("extreme" %in% names(tab)) tab <- tab[!(tab$extreme %in% TRUE), ]
  objs <- tab$object; K <- length(objs); m <- fit$m
  if (K < 3L) stop("need at least three objects")
  beta <- setNames(tab$location, objs)
  cmp <- fit$comparisons
  ia <- match(cmp$object_a, objs); ib <- match(cmp$object_b, objs)
  w <- cmp$weight
  notes <- .btl_dimensionality_scope_note(independent_comparisons)
  S_seen <- .btl_scores(ia, ib, cmp$response, w, m, K)
  n_seen <- sum(((S_seen + t(S_seen)) > 0)[upper.tri(diag(K))])
  complete_pairs <- n_seen == choose(K, 2)

  expected <- .btl_fitted_moments(fit, cmp)$E
  R <- .btl_resid_matrix_expected(ia, ib, cmp$response, w, K, expected, m)
  bm <- .btl_bimensions(R)
  if (!length(bm$strength)) stop("no residual structure to decompose")

  # parametric-bootstrap reference: simulate unidimensional data from the
  # fitted model at the observed pair counts and recompute the leading
  # strength. When the fit carries within-judge dependence effects the null
  # must carry them too -- order effects push the marginal pair rates away
  # from plogis(beta_i - beta_j) in a structured way, and a reference without
  # them reads that structure as a second attribute (false positives)
  tau <- if (m > 1L) fit$thresholds$tau else NULL
  dep <- fit$dependence
  seq_sim <- NULL
  reference_unavailable <- FALSE
  if (!is.null(dep)) {
    dd <- fit$dependence_data                       # sorted by judge, order
    sa <- match(dd$object_a, objs); sb <- match(dd$object_b, objs)
    sw <- dd$weight; sjd <- dd$judge
    weighted_rows <- any(abs(sw - 1) > sqrt(.Machine$double.eps))
    fractional_rows <- any(abs(sw - round(sw)) > sqrt(.Machine$double.eps))
    history_effects <- any(dep$effect %in% c("exposure", "carry_over"))
    # An aggregated row does not contain the order of the comparisons it
    # represents. One draw multiplied by its count has the wrong sampling
    # variance, while expanding it would invent an unidentified history.
    # Count rows remain usable for a position-only fit: their independent
    # outcomes can be drawn binomially or multinomially without a history.
    reference_unavailable <- weighted_rows &&
      (history_effects || fractional_rows)
    coef_exp <- dep$estimate[match("exposure", dep$effect)]
    coef_cry <- dep$estimate[match("carry_over", dep$effect)]
    coef_pos <- dep$estimate[match("position", dep$effect)]
    if (is.na(coef_exp)) coef_exp <- 0
    if (is.na(coef_cry)) coef_cry <- 0
    if (is.na(coef_pos)) coef_pos <- 0
    # sequential simulation mirroring .btl_exposure's history rules, with
    # the FITTED coefficients: seen-before indicator and running mean verdict
    if (!reference_unavailable) seq_sim <- function() {
      cnt <- new.env(hash = TRUE, parent = emptyenv())
      tot <- new.env(hash = TRUE, parent = emptyenv())
      gets <- function(e, k) if (is.null(v <- e[[k]])) 0 else v
      key_a <- .factor_keys(data.frame(judge = sjd, object = sa,
                                       stringsAsFactors = FALSE))
      key_b <- .factor_keys(data.frame(judge = sjd, object = sb,
                                       stringsAsFactors = FALSE))
      resp <- expected <- numeric(length(sa))
      for (r in seq_along(sa)) {
        ka <- key_a[r]; kb <- key_b[r]
        na_ <- gets(cnt, ka); nb_ <- gets(cnt, kb)
        z_exp <- as.numeric(na_ > 0) - as.numeric(nb_ > 0)
        z_cry <- (if (na_ > 0) gets(tot, ka) / na_ else 0) -
                 (if (nb_ > 0) gets(tot, kb) / nb_ else 0)
        # position is a constant +1 for every (object_a-first) row
        lp <- beta[sa[r]] - beta[sb[r]] + coef_exp * z_exp + coef_cry * z_cry +
              coef_pos
        nr <- as.integer(round(sw[r]))
        x <- if (m == 1L) {
          expected[r] <- stats::plogis(lp)
          stats::rbinom(1L, nr, expected[r]) / nr
        } else {
          moments <- item_moments(lp, tau)
          expected[r] <- moments$E
          cc <- stats::rmultinom(1L, nr, moments$P)
          sum((0:m) * cc) / nr
        }
        resp[r] <- x
        cnt[[ka]] <- na_ + sw[r]; cnt[[kb]] <- nb_ + sw[r]
        tot[[ka]] <- gets(tot, ka) + sw[r] * (2 * x / m - 1)
        tot[[kb]] <- gets(tot, kb) + sw[r] * (2 * (m - x) / m - 1)
      }
      list(response = resp, expected = expected)
    }
  }
  # the model-based reference simulates at the PAIR level: n_ij physical
  # comparisons per unordered pair (count weights summed and rounded, so
  # half-tie rows -- weight 0.5 in each direction -- recombine into whole
  # comparisons), binomial or multinomial sums drawn at the fitted pair
  # probabilities. Row-level simulation either overdispersed count weights
  # (one weighted Bernoulli, variance w^2 p(1-p)) or broke entirely on
  # fractional half-tie weights (as.integer(0.5) = 0).
  if (is.null(seq_sim) && !reference_unavailable) {
    pi_lo <- pmin(ia, ib); pi_hi <- pmax(ia, ib)
    pkey <- paste(pi_lo, pi_hi)
    n_pair <- round(tapply(w, pkey, sum))
    n_pair <- pmax(n_pair, 1L)
    u_keys <- names(n_pair)
    u_lo <- as.integer(sub(" .*", "", u_keys))
    u_hi <- as.integer(sub(".* ", "", u_keys))
    d_pair <- beta[u_lo] - beta[u_hi]
    Pcat_pair <- if (m == 1L) NULL else
      vapply(d_pair, function(dd) item_moments(dd, tau)$P, numeric(m + 1L))
    expected_pair <- if (m == 1L) stats::plogis(d_pair) else
      drop((0:m) %*% Pcat_pair)
  }
  # Shared fixed comparison order confounds the dimensionality test. When
  # every judge is given the same comparison sequence, a real within-judge
  # order effect cannot be separated from the object locations (the fitted
  # order coefficient attenuates, the locations absorb the structured order
  # signal), and the fixed-estimate reference cannot reproduce that
  # bias-induced residual structure -- so the observed matrix clears it in
  # almost every unidimensional draw, a false second dimension. Refitting
  # the reference per replicate does not rescue this: it inflates the
  # reference and destroys power even for randomised orders. The honest
  # course is to detect the shared order and withhold the second-dimension
  # verdict, since the design does not identify it. Varying order removes
  # this particular confounding; it does not establish the conditional
  # independence assumption required separately below.
  shared_order <- FALSE
  replicated_order <- TRUE
  # A static first-position term is not a sequence effect. Shared comparison
  # order matters only when exposure or carry-over uses the preceding
  # judgements; applying this guard to a position-only fit withholds an
  # otherwise identified dimensionality result (and, without judges, there
  # is no sequence to compare in the first place).
  if (!is.null(seq_sim) && history_effects) {
    dd <- fit$dependence_data
    ov <- .btl_order_variation(dd, objs)
    shared_order <- ov$shared
    replicated_order <- ov$replicated
  }
  lead_ref <- if (reference_unavailable) numeric(0) else
    vapply(seq_len(reps), function(r) {
    if (is.null(seq_sim)) {
      resp <- if (m == 1L)
        stats::rbinom(length(d_pair), n_pair, stats::plogis(d_pair)) / n_pair
      else vapply(seq_along(d_pair), function(k) {
        cnt <- stats::rmultinom(1L, n_pair[k], Pcat_pair[, k])
        sum((0:m) * cnt) / n_pair[k]
      }, 0)
      Rr <- .btl_resid_matrix_expected(u_lo, u_hi, resp, n_pair, K,
                                      expected_pair, m)
    } else {
      draw <- seq_sim()
      Rr <- .btl_resid_matrix_expected(sa, sb, draw$response, sw, K,
                                      draw$expected, m)
    }
    s <- .btl_bimensions(Rr)$strength
    if (length(s)) s[1] else 0
  }, 0)
  inference <- if (length(lead_ref) && complete_pairs && !shared_order &&
                    independent_comparisons)
    .sim_upper_family(bm$strength[1], matrix(lead_ref, ncol = 1L)) else NULL
  ref_mean <- if (is.null(inference)) NA_real_ else inference$mean[1L]
  ref_p95 <- if (is.null(inference)) NA_real_ else inference$critical[1L]

  nb <- length(bm$strength)
  prop <- 2 * bm$strength^2 / bm$total
  # under a shared fixed order the verdict is not identifiable: withhold it
  # (NA) rather than report a confounded flag
  lead_flag <- if (is.null(inference)) NA else inference$significant[1L]
  bimensions <- data.frame(
    bimension = seq_len(nb), strength = bm$strength,
    prop_residual = prop,
    ref_mean = c(ref_mean, rep(NA_real_, nb - 1L)),
    ref_p95 = c(ref_p95, rep(NA_real_, nb - 1L)),
    above_reference = c(lead_flag, rep(NA, nb - 1L)))
  if (shared_order) {
    order_note <- if (!replicated_order) paste0(
      "no object pair has an observed position for more than one judge, so ",
      "across-judge order variation cannot be assessed and the order effect ",
      "cannot be separated from object structure") else paste0(
      "every judge shares (nearly) the same comparison order, so the ",
      "within-judge order effect is confounded with the object locations ",
      "and a second dimension cannot be separated from it")
    notes <- c(notes, paste0(
      order_note, ": dimensionality inference is withheld. Randomise ",
      "the comparison order across judges to test dimensionality with an ",
      "order effect present"))
  }
  if (reference_unavailable) {
    why <- if (history_effects)
      paste0("the ordered analysis contains aggregated or fractional rows; ",
             "their within-row comparison sequence is unavailable") else
      paste0("the effect-adjusted analysis contains fractional rows that do ",
             "not identify a whole number of independent outcomes")
    notes <- c(notes, paste0(
      why, ", so a valid noise reference cannot be generated and the ",
      "dimensionality inference is withheld"))
  }

  coords <- data.frame(object = objs, location = unname(beta),
                       x = bm$coord[, "x"], y = bm$coord[, "y"])
  if (!complete_pairs)
    notes <- c(notes, sprintf(
      paste0("%d of %d pairs compared; unseen pairs contribute no residual ",
             "information, so dimensionality inference is withheld"),
      n_seen, choose(K, 2)))

  out <- list(bimensions = bimensions, coords = coords,
              # keep the public verdict consistent with above_reference:
              # TRUE/FALSE when identified, NA when inference is withheld
              leading_structured = lead_flag,
              reference = list(
                mean = ref_mean, p95 = ref_p95,
                p = if (is.null(inference)) NA_real_ else inference$p[1L],
                p_adj = if (is.null(inference)) NA_real_ else
                  inference$p_adjusted[1L],
                reps = if (reference_unavailable) 0L else reps,
                n_used = length(lead_ref),
                alpha = if (is.null(inference)) 0.05 else inference$alpha,
                inference_available = !is.null(inference),
                independent_comparisons = independent_comparisons,
                draws = lead_ref,
                method = if (is.null(inference)) NA_character_ else
                  inference$method, seed = seed),
              residual_matrix = R, residual_method = "pooled-expected-score-1",
              notes = notes)
  attr(out, "fit_signature") <- .fit_boot_signature(fit)
  class(out) <- "rasch_btl_dim"
  attr(out, "result_signature") <- .fit_boot_md5(out)
  out
}

#' @export
print.rasch_btl_dim <- function(x, ...) {
  .require_current_btl_residuals(x)
  b <- x$bimensions
  available <- .btl_dimensionality_has_inference(x)
  verdict <- if (!available)
    "inference withheld"
  else if (x$leading_structured) "above the conditional reference; investigate"
  else "within the conditional reference"
  cat(sprintf("Paired-comparison residual dimensionality: %d bimension(s)\n",
              nrow(b)))
  # The categorical verdict is based on the maximum-statistic familywise
  # probability. Display that same quantity; an older result that did not
  # retain it cannot be repaired from the leading-component draws alone.
  p_ref <- x$reference[["p_adj"]]
  if (available && is.finite(x$reference$p95)) {
    probability <- if (length(p_ref) && is.finite(p_ref))
      sprintf("; adjusted p = %.3f", p_ref) else
        "; adjusted p unavailable"
    cat(sprintf("Leading bimension strength %.3f (%.0f%% of residual; reference 5%% upper limit: %.3f%s) -> %s\n",
                b$strength[1], 100 * b$prop_residual[1], x$reference$p95,
                probability, verdict))
  } else {
    cat(sprintf("Leading bimension strength %.3f (%.0f%% of residual; reference unavailable) -> %s\n",
                b$strength[1], 100 * b$prop_residual[1], verdict))
  }
  for (n in x$notes) cat("Note:", n, "\n")
  invisible(x)
}

#' Consistency plot for paired-comparison transitivity
#'
#' With \code{by = "judge"} (the default when judges exist), plots each
#' judge's consistency -- one minus the circular-triad rate over the chance
#' rate -- as a dot against the chance line at zero: the individual-judge
#' lens, a judge-fit analogue. With \code{by = "object"}, plots each object's
#' circular-triad involvement instead: the structural lens, showing which
#' objects sit in the most contradictions.
#'
#' @param x A \code{"rasch_btl_transitivity"} object.
#' @param by \code{"auto"} (judges if present, else objects), \code{"judge"},
#'   or \code{"object"}.
#' @param ... Unused.
#' @return Called for its plotting side effect.
#' @examples
#' \donttest{
#' d <- simulate_btl(6, 10, reps_per_pair = 20, seed = 1)
#' fit <- btl(d, "object_a", "object_b", winner = "winner", judge = "judge")
#' tr <- btl_transitivity(fit)
#' plot_btl_transitivity(tr)
#' }
#' @export
plot_btl_transitivity <- function(x, by = c("auto", "judge", "object"), ...) {
  if (!inherits(x, "rasch_btl_transitivity"))
    stop("`x` must be a result from btl_transitivity()", call. = FALSE)
  by <- match.arg(by)
  use_judge <- if (by == "auto") !is.null(x$judges) else by == "judge"
  if (use_judge && is.null(x$judges))
    stop("no per-judge consistency (no judges, or too few compared triples)")
  if (use_judge) {
    j <- x$judges[order(x$judges$consistency), ]
    n <- nrow(j)
    op <- .rr_canvas(c(min(0, min(j$consistency)) - 0.02, 1), c(0.5, n + 0.5),
                     "Consistency  (1 = one clean order, 0 = guessing)", "",
                     yaxis = FALSE, grid_x = TRUE, grid_y = FALSE)
    on.exit(par(op))
    abline(v = 0, col = .rr$red, lty = 2, lwd = 1.5)
    abline(v = 1, col = .rr$soft, lty = 3)
    segments(0, seq_len(n), j$consistency, seq_len(n), col = .rr$grid, lwd = 3)
    points(j$consistency, seq_len(n), pch = 21, bg = .rr$blue, col = "white",
           cex = 1.5)
    text(par("usr")[1], seq_len(n), j$judge, pos = 4, cex = 0.8, col = .rr$ink,
         offset = 0.3)
    .rr_legend("bottomright", "chance", lwd = 1.5, lty = 2, col = .rr$red)
  } else {
    o <- x$objects
    n <- nrow(o)
    op <- .rr_canvas(c(0, max(o$circular_triads, 1)), c(0.5, n + 0.5),
                     "Circular triads the object sits in", "",
                     yaxis = FALSE, grid_x = TRUE, grid_y = FALSE)
    on.exit(par(op))
    segments(0, seq_len(n), rev(o$circular_triads), seq_len(n),
             col = .rr$grid, lwd = 3)
    points(rev(o$circular_triads), seq_len(n), pch = 21, bg = .rr$blue,
           col = "white", cex = 1.5)
    text(0, seq_len(n), rev(o$object), pos = 4, cex = 0.8, col = .rr$ink,
         offset = 0.3)
  }
}

#' Scree of paired-comparison residual bimensions
#'
#' Bimension strengths against the model-simulated noise reference (its mean
#' and finite-simulation 5% upper reference band), when that reference is
#' available. A leading bar clearing the band suggests residual structure
#' beyond the fitted model under that reference.
#'
#' @param x A \code{"rasch_btl_dim"} object.
#' @param ... Unused.
#' @return Called for its plotting side effect.
#' @examples
#' \donttest{
#' d <- simulate_btl(7, 12, reps_per_pair = 20, seed = 1)
#' fit <- btl(d, "object_a", "object_b", winner = "winner", judge = "judge")
#' dimensions <- btl_dimensionality(fit, reps = 20)
#' plot_btl_scree(dimensions)
#' }
#' @export
plot_btl_scree <- function(x, ...) {
  if (!inherits(x, "rasch_btl_dim"))
    stop("`x` must be a result from btl_dimensionality()", call. = FALSE)
  .require_current_btl_residuals(x)
  b <- x$bimensions; k <- nrow(b)
  ref_m <- x$reference$mean; ref_p <- x$reference$p95
  has_ref <- .btl_dimensionality_has_inference(x) &&
    is.finite(ref_m) && is.finite(ref_p)
  ymax <- max(c(b$strength, if (has_ref) ref_p), na.rm = TRUE) * 1.15
  op <- .rr_canvas(c(0.5, k + 0.5), c(0, ymax), "Bimension", "Strength",
                   grid_x = FALSE, xaxis = FALSE)
  on.exit(par(op))
  rect(seq_len(k) - 0.32, 0, seq_len(k) + 0.32, b$strength,
       col = ifelse(c(isTRUE(x$leading_structured), rep(FALSE, k - 1)),
                    .rr$blue, .rr$soft), border = NA)
  # Noise reference: mean line with a shaded band up to the finite-simulation
  # 5% upper critical value.
  if (has_ref) {
    rect(0.5, ref_m, k + 0.5, ref_p, col = "#dc262622", border = NA)
    abline(h = ref_m, col = .rr$red, lty = 5, lwd = 1.6)
    abline(h = ref_p, col = .rr$red, lty = 2, lwd = 1.8)
  }
  axis(1, at = seq_len(k), col = .rr$grid, col.ticks = .rr$soft)
  .rr_legend("topright",
             if (has_ref) c("Observed", "Null reference band") else
               "Observed",
             fill = if (has_ref) c(.rr$blue, "#dc262633") else .rr$blue,
             border = NA)
}

#' Residual map of the leading paired-comparison bimension
#'
#' Objects placed in the leading bimension plane to show the pattern of
#' residual comparisons. The coordinates describe the pattern, not its
#' magnitude or significance; use \code{\link{plot_btl_scree}} to compare its
#' strength with the conditional reference. Point size grows with the object's
#' location on the primary scale.
#'
#' @param x A \code{"rasch_btl_dim"} object.
#' @param ... Unused.
#' @return Called for its plotting side effect.
#' @examples
#' \donttest{
#' d <- simulate_btl(7, 12, reps_per_pair = 20, seed = 1)
#' fit <- btl(d, "object_a", "object_b", winner = "winner", judge = "judge")
#' dimensions <- btl_dimensionality(fit, reps = 20)
#' plot_btl_dim_map(dimensions)
#' }
#' @export
plot_btl_dim_map <- function(x, ...) {
  if (!inherits(x, "rasch_btl_dim"))
    stop("`x` must be a result from btl_dimensionality()", call. = FALSE)
  .require_current_btl_residuals(x)
  d <- x$coords
  r <- max(sqrt(d$x^2 + d$y^2), 1e-9)
  lim <- c(-r, r) * 1.25
  op <- .rr_canvas(lim, lim, "Leading bimension", "", grid_y = FALSE)
  on.exit(par(op))
  abline(h = 0, v = 0, col = .rr$grid, lwd = 0.8)
  symbols(0, 0, circles = r, add = TRUE, inches = FALSE, fg = .rr$grid,
          lwd = 0.8)
  loc <- d$location; cex <- 1.2 + 2 * (loc - min(loc)) / (max(loc) - min(loc) + 1e-9)
  points(d$x, d$y, pch = 21, bg = .rr$blue, col = "white", cex = cex)
  text(d$x, d$y, d$object, pos = 3, cex = 0.8, col = .rr$ink, offset = 0.5)
}

# Row-level moments under the complete fitted paired-comparison model. This
# includes position, exposure, and carry-over terms retained in comparisons;
# diagnostics based on beta_a - beta_b alone are mis-centred whenever one of
# those effects is present.
.btl_fitted_moments <- function(fit, cmp) {
  # Frame fits store the row-specific expectation because panel units and
  # object-set units mean it cannot be reconstructed from one location gap.
  if (inherits(fit, "rasch_btl_efrm") && "expected" %in% names(cmp)) {
    E <- pmin(pmax(as.numeric(cmp$expected), 1e-12), 1 - 1e-12)
    return(list(E = E, V = E * (1 - E), lp = stats::qlogis(E)))
  }
  objs <- fit$objects$object
  beta <- setNames(fit$objects$location, objs)
  m <- fit$m
  lp <- unname(beta[cmp$object_a] - beta[cmp$object_b])
  if (!is.null(fit$dependence)) {
    eff <- fit$dependence$effect
    if (!all(eff %in% names(cmp)))
      stop("fitted dependence covariates are unavailable in comparisons")
    Z <- as.matrix(cmp[, eff, drop = FALSE])
    lp <- lp + drop(Z %*% fit$dependence$estimate)
  }
  if (m == 1L) {
    E <- stats::plogis(lp)
    V <- E * (1 - E)
  } else {
    tau <- fit$thresholds$tau
    mo <- lapply(lp, item_moments, tau = tau)
    E <- vapply(mo, `[[`, 0, "E")
    V <- vapply(mo, `[[`, 0, "V")
  }
  list(E = E, V = pmax(V, 1e-12), lp = lp)
}

#' Unexpected judgements of one judge
#'
#' The paired-comparison counterpart of the kidmap. A judge has no ability to
#' condition on, so the reference is the consensus object scale (the pooled
#' locations). For the nominated judge, each object it met is given a
#' standardised residual oriented to the object -- how much more (\code{z > 0},
#' over-rated) or less (\code{z < 0}, under-rated) that judge favoured it than
#' its consensus location predicts. Approximate two-sided normal probabilities
#' are adjusted by Holm across the objects meeting \code{min_n}. A surprise
#' is an eligible object treated against its standing (residual opposite in
#' sign to the location) whose adjusted probability passes the level
#' represented by \code{flag_z}. The fitted model must have converged.
#' An adequately sampled object with unavailable residual inference remains
#' in the adjustment family.
#'
#' @param fit A paired-comparison fit from \code{\link{btl}} with judges.
#' @param judge The judge to profile (a value of the fit's judge column).
#' @param min_n Objects met fewer times are shown but never flagged.
#' @param flag_z Absolute normal-residual threshold defining the familywise
#'   flagging level; 1.96 corresponds to an adjusted two-sided probability
#'   of approximately 0.05.
#' @return A list of class \code{"rasch_btl_judge"}: \code{objects} (per object
#'   met: location, times met \code{n}, residual \code{z}, approximate
#'   \code{p}, Holm-adjusted \code{p_adj}, \code{surprise} flag and its
#'   \code{type}). Strong and weak refer to standing above or below the mean
#'   calibrated-object location, so the classification is unchanged by the
#'   arbitrary scale origin. \code{all_locations} contains every object for
#'   orientation;
#'   the \code{judge} and settings.
#' @examples
#' set.seed(1); objs <- LETTERS[1:6]; beta <- setNames(seq(-1.5, 1.5, len = 6), objs)
#' pr <- t(utils::combn(objs, 2))
#' d <- data.frame(a = rep(pr[, 1], each = 12), b = rep(pr[, 2], each = 12))
#' d$judge <- sample(paste0("J", 1:5), nrow(d), TRUE)
#' d$win <- ifelse(runif(nrow(d)) < plogis(beta[d$a] - beta[d$b]), d$a, d$b)
#' judge_surprise(btl(d, "a", "b", "win", judge = "judge"), "J1")
#' @export
judge_surprise <- function(fit, judge, min_n = 2L, flag_z = 1.96) {
  if (!is.atomic(judge) || !is.null(dim(judge)) || length(judge) != 1L ||
      is.na(judge))
    stop("`judge` must be one judge identifier")
  min_n <- .check_whole(min_n, "min_n", 1)
  if (length(flag_z) != 1L || !is.numeric(flag_z) || is.complex(flag_z) ||
      !is.null(dim(flag_z)) || !is.null(oldClass(flag_z)) ||
      !is.finite(flag_z) || flag_z <= 0)
    stop("`flag_z` must be one positive finite flagging value")
  if (!inherits(fit, "rasch_btl")) stop("not a paired-comparison (btl) fit")
  if (!isTRUE(fit$converged))
    stop("the paired-comparison calibration did not converge; judge residuals are unavailable")
  cmp <- fit$comparisons
  if (all(is.na(cmp$judge))) stop("no judges in this fit")
  judge <- .role_text_values(judge)
  if (!nzchar(judge)) stop("`judge` must be one judge identifier")
  sel <- which(!is.na(cmp$judge) & as.character(cmp$judge) == judge)
  if (!length(sel)) stop("no comparisons for judge ", judge)
  tab <- fit$objects
  if ("extreme" %in% names(tab)) tab <- tab[!(tab$extreme %in% TRUE), ]
  objs <- tab$object; K <- length(objs); m <- fit$m
  beta <- setNames(tab$location, objs)
  d <- cmp[sel, , drop = FALSE]
  mo <- .btl_fitted_moments(fit, d)
  ia <- match(d$object_a, objs); ib <- match(d$object_b, objs)
  obs <- exq <- vr <- nn <- numeric(K)
  for (r in seq_len(nrow(d))) {
    a <- ia[r]; b <- ib[r]; w <- d$weight[r]; x <- d$response[r]
    obs[a] <- obs[a] + w * x;         exq[a] <- exq[a] + w * mo$E[r]
    vr[a]  <- vr[a]  + w * mo$V[r];   nn[a]  <- nn[a]  + w
    obs[b] <- obs[b] + w * (m - x);   exq[b] <- exq[b] + w * (m - mo$E[r])
    vr[b]  <- vr[b]  + w * mo$V[r];   nn[b]  <- nn[b]  + w
  }
  keep <- nn > 0
  z <- (obs - exq) / sqrt(pmax(vr, 1e-9))
  o <- data.frame(object = objs, location = unname(beta), n = nn,
                  z = z)[keep, , drop = FALSE]
  # The normal reference is approximate because the consensus locations are
  # estimated. Holm control remains valid under the dependence among the
  # object residuals.
  o$p <- 2 * stats::pnorm(-abs(o$z))
  # The family is declared by the workload rule, before inspecting whether a
  # particular normal reference is available. Dropping an adequately sampled
  # object because its probability is NA would make the remaining Holm
  # probabilities less conservative merely because one question failed.
  planned <- o$n >= min_n
  eligible <- planned & is.finite(o$p)
  o$p_adj <- NA_real_
  if (any(eligible))
    o$p_adj[eligible] <- stats::p.adjust(o$p[eligible], method = "holm",
                                         n = sum(planned))
  family_alpha <- 2 * stats::pnorm(-flag_z)
  # An anchored BTL scale can be translated without changing a single fitted
  # comparison.  Classify an object from its standing relative to the
  # calibrated-object mean, not from the sign of its arbitrary reported
  # origin; otherwise a pure translation changes the diagnostic verdict.
  standing <- o$location - mean(unname(beta))
  o$surprise <- eligible & o$z * standing < 0 &
    o$p_adj <= family_alpha
  o$type <- ifelse(!o$surprise, "",
                   ifelse(standing > 0, "strong object under-rated",
                          "weak object over-rated"))
  o <- o[order(-abs(o$z)), ]; rownames(o) <- NULL
  structure(list(judge = judge, objects = o, all_locations = beta,
                 n_comparisons = sum(cmp$weight[sel]), flag_z = flag_z,
                 family_alpha = family_alpha, p_adjust = "holm",
                 min_n = min_n),
            class = "rasch_btl_judge")
}

#' @export
print.rasch_btl_judge <- function(x, ...) {
  cat(sprintf("Judge %s: %d comparisons over %d objects\n",
              x$judge, x$n_comparisons, nrow(x$objects)))
  s <- x$objects[x$objects$surprise, , drop = FALSE]
  if (nrow(s)) {
    cat("Unexpected judgements:\n")
    for (i in seq_len(nrow(s)))
      cat(sprintf("  %-6s (loc %+.2f): z = %+.2f, Holm p = %s  [%s]\n",
                  s$object[i], s$location[i], s$z[i], .fmt_p(s$p_adj[i]),
                  s$type[i]))
  } else cat("No object judged against its consensus standing.\n")
  invisible(x)
}

#' Unexpected judgements of one judge, pair by pair
#'
#' The comparison-level companion of \code{\link{judge_surprise}}. Each pair
#' the judge met is oriented to its stronger object (higher consensus
#' location) and given a standardised residual: \code{z < 0} means the
#' stronger object won less than its lead predicts -- the judge backed the
#' underdog. Approximate two-sided normal probabilities are adjusted by Holm
#' across matchups meeting \code{min_n}. A surprise is an eligible matchup
#' with a negative residual whose adjusted probability passes the level
#' represented by \code{flag_z}. The fitted model must have converged.
#' An adequately sampled matchup with unavailable residual inference remains
#' in the adjustment family. Locations within \eqn{10^{-10}} logits are treated
#' as tied: neither object is called stronger, the two-sided residual probability
#' remains descriptive and in the Holm family, and \code{surprise} is false.
#' A tied row is oriented toward its non-negative residual for display only.
#'
#' @param fit A paired-comparison fit from \code{\link{btl}} with judges.
#' @param judge The judge to profile.
#' @param min_n Pairs met fewer times are shown but never flagged.
#' @param flag_z Absolute normal-residual threshold defining the familywise
#'   flagging level; 1.96 corresponds to an adjusted two-sided probability
#'   of approximately 0.05.
#' @return A list of class \code{"rasch_btl_judge_pairs"}: \code{pairs} (per
#'   matchup: the stronger and weaker object and their locations when
#'   \code{tied} is false (at a tie these columns only orient the row), the
#'   location \code{gap}, tie indicator \code{tied}, times met \code{n},
#'   residual \code{z}, approximate \code{p},
#'   Holm-adjusted \code{p_adj}, the \code{net_winner}, and the
#'   \code{surprise} flag); \code{all_locations}; the \code{judge} and settings.
#' @examples
#' set.seed(1); objs <- LETTERS[1:6]; beta <- setNames(seq(-1.5, 1.5, len = 6), objs)
#' pr <- t(utils::combn(objs, 2))
#' d <- data.frame(a = rep(pr[, 1], each = 12), b = rep(pr[, 2], each = 12))
#' d$judge <- sample(paste0("J", 1:5), nrow(d), TRUE)
#' d$win <- ifelse(runif(nrow(d)) < plogis(beta[d$a] - beta[d$b]), d$a, d$b)
#' judge_pair_surprise(btl(d, "a", "b", "win", judge = "judge"), "J1")
#' @export
judge_pair_surprise <- function(fit, judge, min_n = 1L, flag_z = 1.96) {
  if (!is.atomic(judge) || !is.null(dim(judge)) || length(judge) != 1L ||
      is.na(judge))
    stop("`judge` must be one judge identifier")
  min_n <- .check_whole(min_n, "min_n", 1)
  if (length(flag_z) != 1L || !is.numeric(flag_z) || is.complex(flag_z) ||
      !is.null(dim(flag_z)) || !is.null(oldClass(flag_z)) ||
      !is.finite(flag_z) || flag_z <= 0)
    stop("`flag_z` must be one positive finite flagging value")
  if (!inherits(fit, "rasch_btl")) stop("not a paired-comparison (btl) fit")
  if (!isTRUE(fit$converged))
    stop("the paired-comparison calibration did not converge; judge residuals are unavailable")
  cmp <- fit$comparisons
  if (all(is.na(cmp$judge))) stop("no judges in this fit")
  judge <- .role_text_values(judge)
  if (!nzchar(judge)) stop("`judge` must be one judge identifier")
  sel <- which(!is.na(cmp$judge) & as.character(cmp$judge) == judge)
  if (!length(sel)) stop("no comparisons for judge ", judge)
  tab <- fit$objects
  if ("extreme" %in% names(tab)) tab <- tab[!(tab$extreme %in% TRUE), ]
  objs <- tab$object; K <- length(objs); m <- fit$m
  beta <- setNames(tab$location, objs)
  d <- cmp[sel, , drop = FALSE]
  mo <- .btl_fitted_moments(fit, d)
  ia <- match(d$object_a, objs); ib <- match(d$object_b, objs)
  tie_tolerance <- 1e-10
  rows <- list()
  for (i in seq_len(K - 1L)) for (j in (i + 1L):K) {
    take <- (ia == i & ib == j) | (ia == j & ib == i)
    if (!any(take)) next
    # At an estimated location tie there is no favourite or underdog.  Use
    # the direction of the residual only to orient the stored row (and hence
    # report a non-negative descriptive z); object order must not decide
    # whether the same judgement is called a directional surprise.
    i_first <- ia[take] == i
    obs_i <- sum(d$weight[take] * ifelse(
      i_first, d$response[take], m - d$response[take]))
    exp_i <- sum(d$weight[take] * ifelse(
      i_first, mo$E[take], m - mo$E[take]))
    tied <- abs(beta[i] - beta[j]) <= tie_tolerance
    hi <- if (tied) {
      if (is.finite(obs_i - exp_i) && obs_i - exp_i < 0) j else i
    } else if (beta[i] > beta[j]) i else j
    lo <- if (hi == i) j else i
    hi_first <- ia[take] == hi
    obs_r <- ifelse(hi_first, d$response[take], m - d$response[take])
    exp_r <- ifelse(hi_first, mo$E[take], m - mo$E[take])
    ww <- d$weight[take]
    n <- sum(ww); obs <- sum(ww * obs_r); ex <- sum(ww * exp_r)
    vv <- sum(ww * mo$V[take]); dd <- beta[hi] - beta[lo]
    zed <- (obs - ex) / sqrt(max(vv, 1e-9))
    rows[[length(rows) + 1L]] <- data.frame(
      object_hi = objs[hi], object_lo = objs[lo],
      loc_hi = unname(beta[hi]), loc_lo = unname(beta[lo]),
      gap = abs(dd), tied = tied, n = n, z = zed,
      net_winner = if (obs >= n * m / 2) objs[hi] else objs[lo],
      stringsAsFactors = FALSE)
  }
  p <- do.call(rbind, rows)
  p$p <- 2 * stats::pnorm(-abs(p$z))
  # Retain every matchup that met the predeclared workload rule in the Holm
  # denominator, including one whose approximate probability is unavailable.
  planned <- p$n >= min_n
  eligible <- planned & is.finite(p$p)
  p$p_adj <- NA_real_
  if (any(eligible))
    p$p_adj[eligible] <- stats::p.adjust(p$p[eligible], method = "holm",
                                         n = sum(planned))
  family_alpha <- 2 * stats::pnorm(-flag_z)
  p$surprise <- eligible & !p$tied & p$z < 0 & p$p_adj <= family_alpha
  p <- p[order(p$z), ]; rownames(p) <- NULL
  structure(list(judge = judge, pairs = p, all_locations = beta,
                 n_comparisons = sum(cmp$weight[sel]), flag_z = flag_z,
                 family_alpha = family_alpha, p_adjust = "holm",
                 min_n = min_n),
            class = "rasch_btl_judge_pairs")
}

#' @export
print.rasch_btl_judge_pairs <- function(x, ...) {
  cat(sprintf("Judge %s: %d comparisons over %d matchups\n",
              x$judge, x$n_comparisons, nrow(x$pairs)))
  s <- x$pairs[x$pairs$surprise, , drop = FALSE]
  if (nrow(s)) {
    cat("Unexpected judgements (weaker object favoured beyond its lead):\n")
    for (i in seq_len(nrow(s)))
      cat(sprintf("  %s vs %s  (gap %.2f, z = %+.2f, Holm p = %s, %s)\n",
                  s$object_hi[i], s$object_lo[i], s$gap[i], s$z[i],
                  .fmt_p(s$p_adj[i]),
                  if (s$net_winner[i] == s$object_lo[i]) "upset"
                  else "favourite under-performed"))
  } else cat("No matchup went against the consensus beyond noise.\n")
  invisible(x)
}

#' Unexpected-judgement map for one judge (pair level)
#'
#' The judge counterpart of the kidmap, drawn matchup by matchup. Each pair the
#' judge met is a segment on the consensus location axis, spanning its two
#' objects, positioned horizontally by its standardised residual. For a
#' non-tied pair, a negative residual means the stronger object under-performed.
#' A tied pair has no directional interpretation and is oriented to a
#' non-negative residual for display only. A filled dot marks the object the
#' judge's verdict favoured and a hollow dot the other. The rug marks every
#' object's location.
#'
#' @param fit A paired-comparison fit from \code{\link{btl}} with judges.
#' @param judge The judge to map.
#' @param min_n,flag_z Passed to \code{\link{judge_pair_surprise}}.
#' @param ... Unused.
#' @return Called for its plotting side effect; invisibly the
#'   \code{rasch_btl_judge_pairs} object.
#' @examples
#' \donttest{
#' d <- simulate_btl(6, 10, reps_per_pair = 20, seed = 1)
#' fit <- btl(d, "object_a", "object_b", winner = "winner", judge = "judge")
#' plot_btl_judge_map(fit, judge = "J1")
#' }
#' @export
plot_btl_judge_map <- function(fit, judge, min_n = 1L, flag_z = 1.96, ...) {
  jp <- judge_pair_surprise(fit, judge, min_n = min_n, flag_z = flag_z)
  p <- jp$pairs
  if (!nrow(p)) stop("this judge made no usable comparisons")
  xr <- max(abs(p$z), flag_z * 1.3)
  yr <- range(jp$all_locations)
  op <- .rr_canvas(c(-xr, xr) * 1.08,
                   yr + c(-1, 1) * (0.12 * diff(yr) + 0.2),
                   "Matchup residual",
                   "Object location (logits)",
                   main = sprintf("Judge %s  \u00b7  %d matchups",
                                  jp$judge, nrow(p)))
  on.exit(par(op))
  u <- par("usr")
  # The band shows the familiar per-match z scale; red verdicts additionally
  # satisfy the Holm-adjusted familywise rule used by judge_pair_surprise().
  rect(-flag_z, u[3], flag_z, u[4], col = "#94a3b81f", border = NA)
  abline(v = 0, col = .rr$soft, lwd = 1.2, lty = 2)
  rug(jp$all_locations, side = 2, col = .rr$grid, lwd = 1.4)
  # each matchup: a segment spanning its two objects at x = its residual, faint
  # when expected and bold red when the judge backed the underdog
  col <- ifelse(p$surprise, .rr$red, .rr$soft)
  lwd <- ifelse(p$surprise, 2.4, 1)
  segments(p$z, p$loc_lo, p$z, p$loc_hi, col = col, lwd = lwd)
  win_y <- ifelse(p$net_winner == p$object_hi, p$loc_hi, p$loc_lo)
  los_y <- ifelse(p$net_winner == p$object_hi, p$loc_lo, p$loc_hi)
  points(p$z, los_y, pch = 21, bg = "white", col = col, cex = 0.7)  # loser end
  points(p$z, win_y, pch = 21, bg = col, col = "white", cex = 1.2)  # winner end
  if (any(p$surprise)) {
    s <- p[p$surprise, ]
    text(s$z, (s$loc_hi + s$loc_lo) / 2,
         sprintf("%s-%s", s$object_hi, s$object_lo), pos = 2, offset = 0.4,
         cex = 0.7, col = .rr$red)
  }
  .rr_legend("bottomright", c("Directional surprise", "Not flagged"),
             lwd = c(2.4, 1), col = c(.rr$red, .rr$soft))
  invisible(jp)
}
