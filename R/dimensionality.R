# rasch :: dimensionality and local dependence
# ===========================================================================
# Residual correlations for local dependence, principal components of the
# residuals, and the Smith (2002) paired t-test of unidimensionality on the
# first residual component.
# ===========================================================================

#' Residual correlations for local dependence (Yen's Q3)
#'
#' The pairwise correlations of the standardised response residuals are
#' Yen's (1984) Q3 statistics. Under unidimensionality and local
#' independence the off-diagonal values sit near \code{-1/(L-1)}; large
#' positive values flag local dependence between item pairs. Following
#' Christensen, Makransky and Horton (2017), each Q3 is also reported
#' relative to the average off-diagonal value (\code{q3_star}). There is no
#' universal adjusted-Q3 critical value: it depends on sample size, test
#' length, category structure, and missingness. The default therefore reports
#' the statistics without a binary flag. A user-supplied \code{flag} is an
#' explicitly heuristic screening threshold, not a calibrated significance
#' test.
#'
#' @param fit A fitted object from \code{\link{rasch}}.
#' @param flag Optional heuristic excess above the average off-diagonal Q3 at
#'   which a pair is flagged. The default \code{NULL} withholds binary flags.
#' @return A list with the Q3 \code{matrix}, the adjusted-Q3
#'   \code{star_matrix} (each Q3 less the average off-diagonal value, diagonal
#'   empty), the \code{average} off-diagonal value, \code{pairs} (every item
#'   pair with \code{q3}, \code{q3_star} and a \code{flagged} indicator, sorted
#'   by \code{q3}), and the subset of \code{flagged} pairs.
#' @references Yen, W. M. (1984). Effects of local item dependence on the
#'   fit and equating performance of the three-parameter logistic model.
#'   \emph{Applied Psychological Measurement}, 8(2), 125-145.
#'
#'   Christensen, K. B., Makransky, G., & Horton, M. (2017). Critical values
#'   for Yen's Q3: identification of local dependence in the Rasch model
#'   using residual correlations. \emph{Applied Psychological Measurement},
#'   41(3), 178-194.
#' @examples
#' set.seed(1)
#' d <- seq(-2, 2, length.out = 8)
#' X <- matrix(rbinom(300 * 8, 1, plogis(outer(rnorm(300), d, "-"))), 300, 8)
#' colnames(X) <- paste0("I", 1:8)
#' residual_correlations(rasch(X))$average
#' @export
residual_correlations <- function(fit, flag = NULL) {
  if (!inherits(fit, "rasch"))
    stop("residual_correlations needs a rasch fit")
  if (!isTRUE(fit$est$converged))
    stop("the fitted calibration did not converge; residual correlations are unavailable")
  if (!is.null(flag) && (length(flag) != 1L || !is.numeric(flag) ||
                         is.complex(flag) || !is.null(dim(flag)) ||
                         !is.null(oldClass(flag)) || !is.finite(flag) ||
                         flag <= 0))
    stop("flag must be NULL or one positive finite heuristic threshold")
  Z <- fit$residuals
  R <- cor(Z, use = "pairwise.complete.obs")
  off <- R[upper.tri(R)]
  finite_off <- is.finite(off)
  if (!any(finite_off))
    .refuse("residual correlations are unavailable because no item pair has ",
            "enough shared, varying residual observations")
  avg <- mean(off[finite_off])
  idx <- which(upper.tri(R), arr.ind = TRUE)
  pairs <- data.frame(item_a = colnames(Z)[idx[, 1]],
                      item_b = colnames(Z)[idx[, 2]],
                      q3 = R[idx], q3_star = R[idx] - avg,
                      flagged = if (is.null(flag)) NA else
                        (R[idx] - avg) > flag)
  pairs <- pairs[!is.na(pairs$q3), ]
  pairs <- pairs[order(-pairs$q3), ]
  rownames(pairs) <- NULL
  # the adjusted-Q3 (Q3*) matrix: each residual correlation less the average
  # off-diagonal Q3, so 0 marks the local-independence baseline. Self-pairs
  # carry no dependence, so the diagonal is left empty.
  star <- R - avg; diag(star) <- NA
  list(matrix = R, star_matrix = star, average = avg, pairs = pairs,
       flagged = if (is.null(flag)) pairs[FALSE, ] else
         pairs[pairs$flagged %in% TRUE, ],
       flag = flag,
       note = if (is.null(flag))
         "binary Q3 flags withheld: no universal critical value"
       else "binary flags use a user-supplied heuristic, not a calibrated test")
}

#' Principal components of the residual correlations
#'
#' The first residual component (PC1) carries any second dimension; items with
#' opposing loadings define the split used by the unidimensionality t-test.
#' Loadings for the leading components and the eigenvalue table support
#' inspection beyond the first component.
#'
#' @param fit A fitted object from \code{\link{rasch}}.
#' @param n_components Number of leading components to return, capped at the
#'   number of items.
#' @return A list with the residual \code{eigenvalues}, their
#'   \code{prop}ortions, the first-component \code{loadings} (sorted), the
#'   \code{loadings_matrix} for the leading components, the
#'   \code{eigen_table} (component, eigenvalue, proportion, cumulative), and
#'   the \code{first_eigen}value.
#' @examples
#' set.seed(1)
#' d <- seq(-2, 2, length.out = 8)
#' X <- matrix(rbinom(300 * 8, 1, plogis(outer(rnorm(300), d, "-"))), 300, 8)
#' colnames(X) <- paste0("I", 1:8)
#' residual_pca(rasch(X))$first_eigen
#' @export
residual_pca <- function(fit, n_components = 10) {
  if (!inherits(fit, "rasch")) stop("residual_pca needs a rasch fit")
  if (!isTRUE(fit$est$converged))
    stop("the fitted calibration did not converge; residual PCA is unavailable")
  n_components <- .check_whole(n_components, "n_components", 1)
  .residual_pca_matrix(fit$residuals, n_components)
}

# The decomposition itself, on a standardised residual matrix rather than a
# fit object, so that a bootstrap replicate of dimensionality_test() can
# decompose its own residuals by exactly the observed computation
.residual_pca_matrix <- function(Z, n_components) {
  R <- cor(Z, use = "pairwise.complete.obs")
  no_overlap <- is.na(R) & row(R) != col(R)
  if (any(no_overlap)) {
    ij <- which(no_overlap, arr.ind = TRUE)[1L, ]
    .refuse("residual PCA is undefined because some item columns have no ",
         "respondents in common (for example ", colnames(R)[ij[1L]], " and ",
         colnames(R)[ij[2L]], "). This happens for structurally disjoint ",
         "designs (item-by-group columns of an extended-frame fit, facet ",
         "cells of a many-facet fit) and equally for sparse or booklet ",
         "missing data; analyse an observable design block, or persons who ",
         "share items, rather than treating non-overlap as zero correlation")
  }
  if (anyNA(R))
    .refuse("residual PCA is undefined for a constant residual column")
  diag(R) <- 1
  ev0 <- eigen(R, symmetric = TRUE)
  adjusted <- min(ev0$values) < -1e-8
  if (adjusted) {
    # Pairwise-complete correlations need not be positive semidefinite. Use
    # the nearest spectral positive-semidefinite correlation approximation
    # rather than clipping eigenvalues only after the decomposition.
    Rp <- ev0$vectors %*% (pmax(ev0$values, 0) * t(ev0$vectors))
    ds <- sqrt(pmax(diag(Rp), 1e-12))
    R <- Rp / outer(ds, ds); diag(R) <- 1
  }
  ev <- eigen(R, symmetric = TRUE)
  k <- min(n_components, ncol(R))
  loadings <- ev$vectors[, 1] * sqrt(pmax(ev$values[1], 0))
  ld <- data.frame(item = colnames(Z), pc1_loading = loadings)
  lm <- ev$vectors[, seq_len(k), drop = FALSE] %*%
    diag(sqrt(pmax(ev$values[seq_len(k)], 0)), k)
  colnames(lm) <- paste0("PC", seq_len(k))
  # a pairwise-complete correlation matrix need not be positive
  # semi-definite; proportions are taken over the positive eigenvalue mass
  # so the cumulative share cannot exceed one
  tot <- sum(pmax(ev$values, 0))
  list(eigenvalues = ev$values, prop = pmax(ev$values, 0) / tot,
       correlation = R, correlation_adjusted = adjusted,
       note = if (adjusted)
         "pairwise correlation matrix projected to a positive-semidefinite correlation matrix"
       else NULL,
       loadings = ld[order(-ld$pc1_loading), ],
       loadings_matrix = data.frame(item = colnames(Z), lm),
       eigen_table = data.frame(component = seq_len(k),
                                eigenvalue = ev$values[seq_len(k)],
                                proportion = pmax(ev$values[seq_len(k)], 0) / tot,
                                cumulative = cumsum(pmax(ev$values[seq_len(k)], 0)) / tot),
       first_eigen = ev$values[1])
}

# Model-simulated eigenvalue reference for the scree plot. Standardised
# Rasch residuals are not independent noise: each person's estimated
# location is a function of their own responses, which couples the residuals
# within a person (off-diagonal correlations near -1/(L-1)), so a random
# normal reference sits systematically BELOW the null first eigenvalue and
# would call structure on model-true data. The reference therefore draws a
# response pattern from the exact Rasch distribution conditional on each
# person's observed score and missingness pattern, then refits both item and
# person parameters before recomputing the residual eigenvalues. This avoids
# treating estimated person locations as known while carrying calibration
# variability through the same estimation chain as the observed eigenvalues.
.sim_upper_family <- function(observed, draws, alpha = 0.05) {
  .check_prob(alpha, "alpha")
  if (!is.matrix(draws)) draws <- as.matrix(draws)
  B <- nrow(draws); k <- ncol(draws)
  if (length(observed) != k || B < 20L || any(!is.finite(draws)) ||
      any(!is.finite(observed)))
    stop("the simulated upper-tail family is incomplete")

  centre <- colMeans(draws)
  p_raw <- vapply(seq_len(k), function(j)
    (1 + sum(draws[, j] >= observed[j])) / (B + 1), numeric(1))
  permitted <- which((1 + 0:B) / (B + 1) <= alpha) - 1L
  if (!length(permitted))
    stop("too few simulated draws to resolve the requested familywise level")
  critical_index <- B - max(permitted)
  marginal_critical <- apply(draws, 2L, function(x)
    sort(x)[critical_index])

  # A single leading statistic is not a multiple-testing family. Its exact
  # finite-simulation probability and corresponding order statistic are the
  # appropriate decision; maximum-statistic standardisation would add
  # conservatism without controlling anything further.
  if (k == 1L) {
    return(list(
      mean = centre, critical = marginal_critical, p = p_raw,
      p_adjusted = p_raw, significant = p_raw <= alpha,
      max_null = draws[, 1L], alpha = alpha, n_used = B,
      method = "finite simulated upper-tail"))
  }

  spread <- apply(draws, 2L, stats::sd)
  stable <- is.finite(spread) & spread > sqrt(.Machine$double.eps)
  # Each simulated maximum must be externally standardised. If a draw helps
  # estimate its own mean and SD, its departure is shrunk relative to the
  # observed value, which is not part of those estimates, and the familywise
  # test becomes anti-conservative. Leave each row out when standardising that
  # row; the observed vector uses all simulated rows.
  Z <- matrix(0, B, k)
  for (i in seq_len(B)) {
    training <- draws[-i, , drop = FALSE]
    centre_i <- colMeans(training)
    spread_i <- apply(training, 2L, stats::sd)
    stable_i <- is.finite(spread_i) & spread_i > sqrt(.Machine$double.eps)
    Z[i, stable_i] <- (draws[i, stable_i] - centre_i[stable_i]) /
      spread_i[stable_i]
    if (any(!stable_i)) {
      tol_i <- sqrt(.Machine$double.eps) *
        pmax(1, abs(centre_i[!stable_i]))
      delta_i <- draws[i, !stable_i] - centre_i[!stable_i]
      Z[i, !stable_i] <- ifelse(delta_i > tol_i, Inf,
                                ifelse(delta_i < -tol_i, -Inf, 0))
    }
  }
  z_observed <- numeric(k)
  z_observed[stable] <- (observed[stable] - centre[stable]) / spread[stable]
  if (any(!stable)) {
    tol <- sqrt(.Machine$double.eps) * pmax(1, abs(centre[!stable]))
    delta <- observed[!stable] - centre[!stable]
    z_observed[!stable] <- ifelse(delta > tol, Inf,
                                  ifelse(delta < -tol, -Inf, 0))
  }

  max_null <- apply(Z, 1L, max)
  p_adjusted <- vapply(z_observed, function(z)
    (1 + sum(max_null >= z)) / (B + 1), numeric(1))
  p_adjusted <- pmax(p_raw, p_adjusted)

  # Match crossing the plotted curve to p_adjusted <= alpha. With 20 draws
  # the smallest probability is 1/21, so the requested five-percent test is
  # resolvable; larger B makes the critical curve progressively less coarse.
  critical_z <- sort(max_null)[critical_index]
  critical <- pmax(
    marginal_critical,
    centre + critical_z * ifelse(stable, spread, 0))

  list(mean = centre, critical = critical, p = p_raw,
       p_adjusted = p_adjusted, significant = p_adjusted <= alpha,
       max_null = max_null, alpha = alpha, n_used = B,
       method = "single-step leave-one-out maximum-standardised-statistic")
}

.scree_reference <- function(fit, k, reps, seed = NULL) {
  .require_refittable_calibration(fit)
  if (inherits(fit, "rasch_efrm") || inherits(fit, "rasch_mfrm"))
    .refuse("parallel residual reference is not available for mutually exclusive ",
         "EFRM/MFRM virtual designs; fit and analyse an observable design block")
  if (.has_repeated_residual_units(fit))
    .refuse("parallel residual inference is not available when a person ",
            "identifier occurs on several response rows: the current null ",
            "generator treats rows as independent occasions and would not ",
            "preserve within-person dependence. Residual eigenvalues remain ",
            "available with parallel = FALSE")
  reps <- .check_whole(reps, "reps", 20)
  if (!is.null(seed)) {
    seed <- .check_whole(seed, "seed", 0)
    old_seed <- .sim_seed_capture()
    on.exit(.sim_seed_restore(old_seed), add = TRUE)
    set.seed(seed)
  }
  tau_list <- fit$tau_list; L <- length(tau_list)
  disc_v <- if (is.null(fit$disc)) rep(1, L) else fit$disc
  if (any(abs(disc_v - 1) > 1e-12))
    stop("full-refit scree reference currently requires a common Rasch unit")
  keep <- !is.na(fit$person$theta)
  X <- fit$X[keep, , drop = FALSE]
  obs <- !is.na(X)
  spec <- fit$refit_spec
  if (is.null(spec)) spec <- list()
  first_error <- NULL
  status <- rep("used", reps)
  draws <- lapply(seq_len(reps), function(r) {
    Xr <- .fit_gen_conditional(X, tau_list, !obs)
    # the reference distribution must be analysed under the model that was
    # fitted: refitting an explanatory calibration with an unrestricted
    # rasch() would compare observed eigenvalues against a freer model
    fr <- tryCatch(suppressWarnings(
      if (inherits(fit, "rasch_explanatory"))
        .explanatory_refit_modified(fit, Xr, person_rows = which(keep),
                                    inherit_mc = FALSE)
      else
        rasch(Xr, model = fit$model, n_groups = .refit_n_groups(fit),
              anchors = spec$anchors, pc_components = spec$pc_components,
              maxit = spec$maxit %||% 60, tol = spec$tol %||% 1e-8)),
      error = function(e) {
        if (is.null(first_error)) first_error <<- conditionMessage(e)
        status[r] <<- "error"
        NULL
      })
    if (is.null(fr)) return(rep(NA_real_, k))
    if (!isTRUE(fr$est$converged)) {
      status[r] <<- "nonconverged"
      return(rep(NA_real_, k))
    }
    # A sparse polytomous replicate can lose its highest observed category.
    # rasch() then fits a shorter item without necessarily dropping its
    # column, so a column-count check alone admits a different null model.
    # Keep the reference conditional on successful refits of the observed
    # score structure, as the item-fit and dimensionality bootstraps do.
    if (ncol(fr$X) != L ||
        !identical(as.integer(fr$m), as.integer(fit$m))) {
      status[r] <<- "error"
      return(rep(NA_real_, k))
    }
    tryCatch(residual_pca(fr, n_components = k)$eigenvalues[seq_len(k)],
      error = function(e) {
        if (is.null(first_error)) first_error <<- conditionMessage(e)
        status[r] <<- "error"
        rep(NA_real_, k)
      })
  })
  sim <- do.call(rbind, draws)
  # A replicated eigenvalue must be finite. complete.cases() regards Inf as
  # complete, which would contaminate every simulated centre and critical
  # value if a numerically failed refit returned an infinite value.
  complete <- rowSums(is.finite(sim)) == ncol(sim)
  status[!complete & status == "used"] <- "error"
  # The finite-simulation family below needs at least 20 rows to resolve a
  # five-percent upper-tail probability.  The generic bootstrap guard is less
  # demanding for small requested runs, so impose the inferential minimum
  # here as well rather than returning an object that plot_scree() cannot use.
  minimum_usable <- max(20L, .fit_min_boot_success(reps))
  if (sum(complete) < minimum_usable)
    .fit_boot_refuse(
      "only ", sum(complete), " of ", reps,
      " full-refit scree replicates retained the fitted response structure; ",
      "at least ", minimum_usable, " are needed for the 5% reference",
      if (!is.null(first_error))
        paste0("; the first replicate failed with: ", first_error) else "",
      B = reps, B_used = sum(complete),
      B_nonconverged = sum(status == "nonconverged"),
      B_errors = sum(status == "error"))
  sim <- sim[complete, , drop = FALSE]
  # The plotted curve is completed once the observed eigenvalues are known.
  # Retain the full joint draws so plot_scree() can use their maximum
  # standardised departure to control the component family.
  alpha <- 0.05
  reference <- colMeans(sim)
  attr(reference, "mean") <- colMeans(sim)
  attr(reference, "draws") <- sim
  attr(reference, "alpha") <- alpha
  attr(reference, "n_requested") <- reps
  attr(reference, "n_used") <- nrow(sim)
  attr(reference, "n_nonconverged") <- sum(status == "nonconverged")
  attr(reference, "n_errors") <- sum(status == "error")
  attr(reference, "minimum_usable") <- minimum_usable
  attr(reference, "seed") <- seed
  reference
}

#' Scree plot of the residual components with parallel analysis
#'
#' Eigenvalues of the residual correlation matrix for the leading components,
#' with a model-simulated parallel-analysis reference: response patterns are
#' drawn conditional on each person's observed score and missingness pattern,
#' the item calibration and every person are re-estimated, and the residual
#' eigenvalues recomputed. The plotted reference is a finite-simulation 5%
#' familywise upper critical curve, obtained from the maximum standardised
#' departure across the displayed components. Each simulated maximum is
#' standardised against the other simulated draws so that it is comparable
#' with the externally standardised observed value. The returned table also gives
#' the reference mean, marginal upper-tail probability and single-step
#' adjusted probability.
#' Because estimating
#' the person locations couples the residuals within a person, this reference
#' sits above the classical random-normal one and is calibrated under the
#' fitted model (Raiche 2005; Chou & Wang 2010). An observed eigenvalue above
#' the critical reference has a familywise-adjusted simulated upper-tail
#' probability at or below .05 and suggests structure beyond what the fitted model
#' produces.
#'
#' @param fit A fitted object from \code{\link{rasch}}. A simulated reference
#'   requires one response row per person; with repeated identifiers, use
#'   \code{parallel = FALSE} to display the observed eigenvalues alone.
#'   Fully anchored scoring fits also require \code{parallel = FALSE}.
#' @param n_components Number of leading components to display. The familywise
#'   adjustment covers these components.
#' @param parallel Draw the parallel-analysis reference band.
#' @param reps Model-simulated replicates for the reference; at least 20 when
#'   \code{parallel = TRUE}. Larger values give a more stable upper-tail
#'   reference.
#' @param seed Optional non-negative whole-number seed. The caller's random-
#'   number state is restored when the calculation finishes; see
#'   \code{\link{rasch_rng}} for generator support.
#' @param result Optional result returned by an earlier call. Supplying it
#'   redraws that analysis without repeating the simulations.
#' @return Called for its plotting side effect; invisibly the eigen table. With
#'   parallel analysis it also contains \code{reference_mean},
#'   \code{reference_critical}, \code{parallel_p}, \code{parallel_p_adj},
#'   \code{parallel_significant}, and requested, usable, non-converged and
#'   other-failure reference counts. The adjustment is
#'   recorded in the table's \code{parallel_adjustment} attribute.
#' @references Raiche, G. (2005). Critical eigenvalue sizes (variances) in
#'   standardized residual principal components analysis. \emph{Rasch
#'   Measurement Transactions}, 19(1), 1012.
#'
#'   Chou, Y.-T., & Wang, W.-C. (2010). Checking dimensionality in item
#'   response models with principal component analysis on standardized
#'   residuals. \emph{Educational and Psychological Measurement}, 70(5),
#'   717-731.
#'
#'   Westfall, P. H., & Young, S. S. (1993). \emph{Resampling-Based Multiple
#'   Testing: Examples and Methods for p-Value Adjustment}. Wiley.
#' @examples
#' set.seed(1)
#' d <- seq(-2, 2, length.out = 8)
#' X <- matrix(rbinom(300 * 8, 1, plogis(outer(rnorm(300), d, "-"))), 300, 8)
#' colnames(X) <- paste0("I", 1:8)
#' plot_scree(rasch(X), reps = 20)
#' @name plot_scree
NULL

.scree_analysis <- function(fit, n_components = 10, parallel = TRUE, reps = 50,
                            seed = NULL) {
  .check_flag(parallel, "parallel")
  n_components <- .check_whole(n_components, "n_components", 1)
  reps <- .check_whole(reps, "reps", if (parallel) 20 else 1)
  pc <- residual_pca(fit, n_components)
  k <- nrow(pc$eigen_table)
  obs <- pc$eigen_table$eigenvalue
  pa <- if (parallel) .scree_reference(fit, k, reps, seed = seed) else NULL
  if (!is.null(pa)) {
    draws <- attr(pa, "draws")
    n_requested <- attr(pa, "n_requested")
    n_nonconverged <- attr(pa, "n_nonconverged")
    n_errors <- attr(pa, "n_errors")
    reference_seed <- attr(pa, "seed")
    inference <- .sim_upper_family(obs, draws, attr(pa, "alpha"))
    pc$eigen_table$reference_mean <- inference$mean
    pc$eigen_table$reference_critical <- inference$critical
    pc$eigen_table$parallel_p <- inference$p
    pc$eigen_table$parallel_p_adj <- inference$p_adjusted
    pc$eigen_table$parallel_significant <- inference$significant
    pc$eigen_table$n_reference <- inference$n_used
    pc$eigen_table$n_reference_requested <- n_requested
    pc$eigen_table$n_reference_nonconverged <- n_nonconverged
    pc$eigen_table$n_reference_errors <- n_errors
    pc$eigen_table$reference_seed <- reference_seed %||% NA_integer_
    attr(pc$eigen_table, "parallel_adjustment") <- inference$method
  }
  out <- pc$eigen_table
  class(out) <- c("rasch_scree", class(out))
  attr(out, "fit_signature") <- .fit_boot_signature(fit)
  attr(out, "parallel") <- parallel
  attr(out, "result_signature") <- .fit_boot_md5(out)
  out
}

.validate_scree_result <- function(result, fit) {
  if (is.null(result)) return(invisible(NULL))
  signature <- attr(result, "result_signature")
  unsigned <- result
  attr(unsigned, "result_signature") <- NULL
  if (!inherits(result, "rasch_scree") || !is.data.frame(result) ||
      !all(c("component", "eigenvalue") %in% names(result)) ||
      !is.character(signature) || length(signature) != 1L || is.na(signature) ||
      !.fit_boot_hash_matches(signature, unsigned) ||
      !.fit_boot_signature_matches(attr(result, "fit_signature"), fit))
    stop("`result` must be a plot_scree() result from this fitted model")
  if (isTRUE(attr(result, "parallel"))) .require_refittable_calibration(fit)
  invisible(result)
}

#' @rdname plot_scree
#' @export
plot_scree <- function(fit, n_components = 10, parallel = TRUE, reps = 50,
                       seed = NULL, result = NULL) {
  if (is.null(result))
    result <- .scree_analysis(fit, n_components, parallel, reps, seed)
  else .validate_scree_result(result, fit)
  k <- nrow(result)
  obs <- result$eigenvalue
  has_reference <- all(c("reference_mean", "reference_critical") %in%
                         names(result))
  pa_mean <- if (has_reference) result$reference_mean else NULL
  pa <- if (has_reference) result$reference_critical else NULL
  significant <- if ("parallel_significant" %in% names(result))
    result$parallel_significant else rep(FALSE, k)
  ylim <- c(0, max(c(obs, pa, 1), na.rm = TRUE) * 1.12)
  op <- .rr_canvas(c(0.5, k + 0.5), ylim, "Component", "Eigenvalue",
                   grid_x = FALSE, xaxis = FALSE)
  on.exit(par(op))
  if (!is.null(pa))
    polygon(c(seq_len(k), rev(seq_len(k))),
            c(pa_mean, rev(pa)), col = "#dc262622", border = NA)
  abline(h = 1, lty = 3, col = .rr$soft)
  if (!is.null(pa)) {
    lines(seq_len(k), pa_mean, lwd = 1.5, lty = 3, col = .rr$red)
    lines(seq_len(k), pa, lwd = 2, lty = 5, col = .rr$red)
  }
  lines(seq_len(k), obs, lwd = 2.6, col = .rr$blue)
  points(seq_len(k), obs, pch = 21,
         bg = ifelse(significant, .rr$red, .rr$blue),
         col = "white", cex = 1.5)
  axis(1, at = seq_len(k), col = .rr$grid, col.ticks = .rr$soft)
  if (is.null(pa)) {
    .rr_legend("topright", "Observed", lwd = 2.6, lty = 1,
               col = .rr$blue)
  } else {
    .rr_legend(
      "topright",
      c("Observed", "Null reference band", "Adjusted p <= .05"),
      lwd = c(2.6, NA, NA), lty = c(1, NA, NA),
      pch = c(NA, 22, 21), pt.bg = c(NA, "#dc262633", .rr$red),
      pt.cex = c(NA, 1.4, 1.2),
      col = c(.rr$blue, "#dc262666", "white"))
  }
  invisible(result)
}

#' Residual-component test of unidimensionality
#'
#' Estimates each person separately on two item subsets and compares the two
#' estimates with a per-person t-test (Smith 2002). By default the subsets
#' are defined by the sign of a residual-component loading (the first by
#' default; any leading component may be chosen); they can also be nominated
#' manually (for example, by content). Under
#' unidimensionality and local independence the two subset estimates are
#' independent given the person location, so
#' \code{t = (theta_A - theta_B) / sqrt(se_A^2 + se_B^2)} is approximately
#' standard normal and about \code{alpha} of the tests should reach
#' significance. Persons with an extreme score on either subset are excluded
#' (their weighted-likelihood estimates are most biased there). The
#' proportion of significant tests is reported with a Clopper--Pearson
#' binomial confidence interval. For a split fixed in
#' advance, a lower bound above \code{alpha} signals multidimensionality. The
#' test requires a converged calibration and one response row per person.
#'
#' The binomial reading holds for a split fixed in advance. A split chosen
#' from the residuals is chosen to make the two subsets disagree, so its
#' proportion runs above \code{alpha} under unidimensionality. Package
#' simulations confirmed that applying the fixed-split binomial rule after
#' choosing the split from the same residuals is anti-conservative. Without
#' bootstrap calibration the data-driven split therefore has no binary verdict:
#' \code{multidimensional} is \code{NA}, while the interval and uncalibrated
#' binomial reading remain available descriptively. Two inferential routes
#' are available. A content-based split, named through \code{items_positive}
#' and \code{items_negative}, supports the conventional fixed-split binomial
#' rule without the selection induced by the residual-derived split. Otherwise
#' \code{B > 0} calibrates the data-driven split by a parametric bootstrap:
#' each replicate draws responses from the fitted model conditional on every
#' person's raw score and missingness pattern, refits the calibration,
#' repeats the residual-component split on its own residuals and recomputes
#' the proportion, so the bootstrap probability \code{p_boot} carries the
#' same selection the observed proportion carries. With \code{B > 0} the
#' verdict is \code{p_boot <= alpha}; the binomial interval is still reported,
#' as a description of the observed proportion rather than a test of it.
#' A one-sided bootstrap probability cannot be smaller than
#' \code{1/(B_used + 1)}. If that floor exceeds \code{alpha}, a data-driven
#' split has no rejection region and its verdict is withheld. A split fixed
#' in advance retains its binomial verdict in that case.
#'
#' The mean difference between the two subset estimates is reported but not
#' tested. Each subset estimate is a weighted-likelihood estimate on a short
#' test, and the two subsets differ in difficulty, so their estimation bias
#' differs systematically: under a perfectly unidimensional Rasch model the
#' expected difference is non-zero whenever the subsets are not matched in
#' targeting, and it grows relative to its standard error with the number of
#' persons. A t-test of that difference therefore tests the targeting of the
#' split rather than its dimensionality -- it rejects for every sample large
#' enough on unidimensional data -- so the difference is reported as a
#' description of the split and the inference is withheld. The person-level
#' comparisons, whose standard errors carry each person's subset
#' uncertainty, are the test.
#'
#' @param fit A fitted object from \code{\link{rasch}} with one response row
#'   per person. Repeated identifiers are refused because the person-level
#'   comparisons and their binomial count would not be independent.
#'   Fully anchored scoring fits require \code{B = 0}; bootstrap refitting is
#'   not supported for these fits.
#' @param alpha Nominal significance level for the per-person t-tests.
#' @param items_positive,items_negative Optional character vectors naming the
#'   two item subsets; both must be given (disjoint, at least two items
#'   each), otherwise the sign of a residual component defines the split.
#' @param component Which residual principal component's loading sign defines
#'   the default split (ignored when subsets are named). Default the first
#'   component.
#' @param min_score_points Score-point threshold below which the verdict
#'   carries a caution. Andrich and Marais (2019) recommend subtests of
#'   roughly 15 score points for stable subtest estimates; shorter subsets
#'   (the norm for ordinary dichotomous tests) retain the analysis, with a
#'   \code{caution} field noting the reduced stability. A quiet verdict
#'   under caution is inconclusive, not clean: with a four-item subtest the
#'   test lacks power where nonparametric alternatives still flag.
#' @param B Number of parametric-bootstrap replicates that calibrate the
#'   proportion of significant tests under the fitted model (see Details).
#'   The default \code{0} reports the binomial interval and descriptive
#'   reading alone; an automatic split then has no inferential verdict. Each
#'   replicate refits the calibration, so \code{B = 200} costs about two hundred
#'   fits; the bootstrap is available for single-facet fits with a common
#'   unit whose thresholds were estimated directly.
#' @param workers Number of parallel workers for the bootstrap refits.
#' @param seed Optional integer seed for the bootstrap; the replicates are
#'   reproducible for a given seed whatever the worker count. The bootstrap
#'   does not support Box--Muller, including when \code{seed = NULL};
#'   see \code{\link{rasch_rng}}.
#' @return A list with the proportion of significant tests, its
#'   Clopper--Pearson confidence interval, the sample sizes (\code{n} used,
#'   \code{n_excluded_extreme}), the item split and its source, a
#'   \code{multidimensional} verdict, the corresponding uncalibrated
#'   \code{binomial_multidimensional} reading, a \code{caution} note when the
#'   subtests fall short of \code{min_score_points}, and
#'   \code{subset_mean_difference}, the mean and standard deviation of the
#'   person-level differences between the two subset estimates, reported
#'   descriptively with the \code{note} that no test accompanies them (see
#'   Details). With \code{B > 0} the list also carries \code{p_boot}, the
#'   bootstrap probability of a proportion at least as large as the observed
#'   one under the fitted unidimensional model; \code{prop_null}, the mean
#'   replicate proportion (the rate the split produces when nothing is
#'   there); \code{bootstrap_resolution}, the smallest attainable bootstrap
#'   probability; and \code{bootstrap}, the replicate proportions with the
#'   counts requested, used, non-converged and failed. When the comparison
#'   itself is unavailable (undefined split, degenerate subsets, too few
#'   persons) the list carries a \code{note} explaining why and
#'   \code{multidimensional = NA}. Every result carries \code{algorithm},
#'   the stamp of the calculation that produced it; a saved result without
#'   the current stamp reported the superseded paired t-test of the subset
#'   means, so an analysis file carrying one opens with that result dropped
#'   and a warning, and the rest of the analysis intact.
#' @references
#' Smith, E. V. Jr. (2002). Detecting and evaluating the impact of
#' multidimensionality using item fit statistics and principal component
#' analysis of residuals. Journal of Applied Measurement, 3(2), 205--231.
#'
#' Tennant, A., & Pallant, J. F. (2006). Unidimensionality matters! (A tale
#' of two Smiths?). Rasch Measurement Transactions, 20(1), 1048--1051.
#' @examples
#' set.seed(1)
#' d <- seq(-2, 2, length.out = 8)
#' X <- matrix(rbinom(300 * 8, 1, plogis(outer(rnorm(300), d, "-"))), 300, 8)
#' colnames(X) <- paste0("I", 1:8)
#' dimensionality_test(
#'   rasch(X), items_positive = paste0("I", 1:4),
#'   items_negative = paste0("I", 5:8))$multidimensional
#' \donttest{
#' # calibrate the data-driven split under the fitted model
#' dimensionality_test(rasch(X), B = 99, workers = 1, seed = 1)$p_boot
#' }
#' @export
dimensionality_test <- function(fit, alpha = 0.05, items_positive = NULL,
                                items_negative = NULL, component = 1,
                                min_score_points = 15L, B = 0,
                                workers = 4L, seed = NULL) {
  if (!inherits(fit, "rasch")) stop("dimensionality_test needs a rasch fit")
  if (!isTRUE(fit$est$converged))
    stop("the fitted calibration did not converge; the dimensionality test is unavailable")
  if (.has_repeated_residual_units(fit))
    .refuse("the person-subset dimensionality test requires one response ",
            "row per person. Repeated identifiers make the person tests ",
            "dependent and the binomial interval invalid; analyse one ",
            "occasion at a time or use a method defined for repeated ",
            "measurements")
  .check_prob(alpha, "alpha")
  component <- .check_whole(component, "component", 1)
  min_score_points <- .check_whole(min_score_points, "min_score_points", 1)
  B <- .check_whole(B, "B", 0)
  workers <- .check_whole(workers, "workers", 1)
  if (!is.null(seed)) seed <- .check_whole(seed, "seed", 0)
  if (B > 0L) .dim_bootstrap_check(fit)
  X <- fit$X
  disc <- if (is.null(fit$disc)) rep(1, ncol(X)) else fit$disc
  manual <- !is.null(items_positive) || !is.null(items_negative)
  if (manual) {
    if (is.null(items_positive) || is.null(items_negative))
      stop("supply both item subsets, or neither")
    resolve_subset <- function(x, name) {
      if (!(is.character(x) || is.factor(x)) || !is.null(dim(x)) ||
          !length(x) || anyNA(x) || any(!nzchar(trimws(as.character(x)))))
        stop("`", name, "` must be a non-empty ordinary vector of item names",
             call. = FALSE)
      match(as.character(x), colnames(X))
    }
    pos <- resolve_subset(items_positive, "items_positive")
    neg <- resolve_subset(items_negative, "items_negative")
    if (anyNA(pos) || anyNA(neg))
      stop("subset item(s) not in the fit: ",
           paste(c(as.character(items_positive)[is.na(pos)],
                   as.character(items_negative)[is.na(neg)]),
                 collapse = ", "))
    repeated <- c(pos[duplicated(pos)], neg[duplicated(neg)])
    if (length(repeated))
      stop("item(s) named more than once in a subset: ",
           paste(unique(colnames(X)[repeated]), collapse = ", "),
           "; a repeated item would count its score points repeatedly")
    both <- intersect(pos, neg)
    if (length(both))
      stop("item(s) in both subsets: ",
           paste(colnames(X)[both], collapse = ", "),
           "; the subsets must be disjoint")
    split_source <- "manual"
    first_eigen <- NA_real_
  } else {
    # the PCA is needed only to derive the automatic split; when it is
    # undefined (structurally disjoint columns, sparse overlap) that is a
    # reason to report, not an error to crash every downstream consumer
    split <- tryCatch(.dim_split(fit$residuals, component), error = function(e)
      structure(list(msg = conditionMessage(e)), class = "rr_pca_refusal"))
    if (inherits(split, "rr_pca_refusal"))
      return(.dimensionality_test_result(
        list(note = paste0("dimensionality split unavailable: ", split$msg),
             multidimensional = NA), fit))
    if (is.null(split))
      stop("component ", component, " is not available")
    pos <- split$pos; neg <- split$neg
    split_source <- sprintf("residual component %d", as.integer(component))
    first_eigen <- split$first_eigen
  }
  if (length(pos) < 2 || length(neg) < 2)
    return(.dimensionality_test_result(
      list(note = "need >= 2 items in each subset",
           multidimensional = NA, split = split_source,
           items_positive = colnames(X)[pos],
           items_negative = colnames(X)[neg]), fit))
  score_points <- c(positive = sum(fit$m[pos]), negative = sum(fit$m[neg]))
  # Andrich & Marais (2019) recommend subtests of roughly 15 score points
  # for STABLE subtest estimates. Short subtests make the person-level
  # comparison noisier, not undefined -- so the verdict is still computed,
  # with a caution, rather than withheld for every ordinary short test
  caution <- if (any(score_points < min_score_points)) sprintf(
    paste0("subtests carry only %d and %d score points (fewer than the ~%d ",
           "recommended for stable subtest estimates); read the verdict ",
           "cautiously"),
    score_points[1], score_points[2], as.integer(min_score_points)) else NULL
  tt <- .dim_ttest(X, fit$tau_list, disc, pos, neg, alpha)
  n <- tt$n
  if (n < 10) return(.dimensionality_test_result(
    list(note = "fewer than 10 usable persons for the t-test",
         multidimensional = NA, split = split_source,
         items_positive = colnames(X)[pos],
         items_negative = colnames(X)[neg], score_points = score_points,
         n = n, n_excluded_extreme = tt$n_excluded_extreme), fit))
  bt <- stats::binom.test(tt$n_sig, n, p = alpha)
  # the person-level differences behind the comparison. Degenerate subsets
  # (e.g. two-item manual subtests where every usable person has the same
  # difference) carry no comparison at all: report the degeneracy instead
  dd <- tt$difference
  if (!is.finite(stats::sd(dd)) || stats::sd(dd) < 1e-12)
    return(.dimensionality_test_result(list(note = paste0(
        "dimensionality verdict withheld: the subset person estimates are ",
        "degenerate (no variation in the paired differences) -- the subsets ",
        "are too short or too sparsely answered for the comparison"),
        multidimensional = NA, split = split_source,
        items_positive = colnames(X)[pos], items_negative = colnames(X)[neg],
        score_points = score_points), fit))
  binomial_verdict <- bt$conf.int[1] > alpha
  out <- list(prop_significant = tt$n_sig / n, ci = as.numeric(bt$conf.int),
              n = n, n_excluded_extreme = tt$n_excluded_extreme,
              multidimensional = if (manual) binomial_verdict else NA,
              binomial_multidimensional = binomial_verdict,
              verdict_method = if (manual) "fixed-split binomial" else
                "withheld for data-driven split",
              split = split_source,
              score_points = score_points,
              caution = caution,
              items_positive = colnames(X)[pos], items_negative = colnames(X)[neg],
              first_eigenvalue = first_eigen,
              alpha = alpha, component = if (manual) NA_integer_ else component,
              min_score_points = min_score_points, B = B,
              workers = workers, seed = seed,
              # descriptive only: the two subset estimates carry different
              # estimation bias whenever the subsets differ in difficulty, so
              # a test of this mean against zero rejects on unidimensional
              # data as soon as the sample is large enough (see Details)
              subset_mean_difference = list(
                mean_difference = mean(dd), sd_difference = stats::sd(dd),
                n = n, note = paste0(
                  "descriptive: the subsets differ in targeting as well as ",
                  "in content, so no test of this mean against zero is a ",
                  "test of unidimensionality; read the person-level ",
                  "comparisons for that")))
  if (B > 0L) {
    boot <- .dim_bootstrap(fit, pos = pos, neg = neg, manual = manual,
                           component = component, alpha = alpha, B = B,
                           workers = workers, seed = seed)
    out$p_boot <- .boot_p(out$prop_significant, boot$null, "upper")
    out$prop_null <- mean(boot$null)
    out$bootstrap_resolution <- 1 / (boot$B_used + 1)
    if (out$bootstrap_resolution <= alpha) {
      out$multidimensional <- out$p_boot <= alpha
      out$verdict_method <- "parametric bootstrap"
    } else {
      # A p-value floor above alpha gives the bootstrap no rejection region.
      # Do not turn that inability to test into evidence of fit. A manual
      # split still has its pre-specified binomial test; a data-driven split
      # has no fallback inferential verdict.
      warning(sprintf(paste0(
        "with %d usable replicates the smallest attainable bootstrap p is ",
        "%.3f (> alpha %.3f); the bootstrap verdict is unavailable"),
        boot$B_used, out$bootstrap_resolution, alpha), call. = FALSE)
      out$multidimensional <- if (manual) binomial_verdict else NA
      out$verdict_method <- if (manual)
        "fixed-split binomial (bootstrap resolution insufficient)" else
        "withheld: bootstrap resolution insufficient for data-driven split"
    }
    out$bootstrap <- boot
  }
  .dimensionality_test_result(out, fit)
}

.dimensionality_test_algorithm <- "person-subset-comparison-1"

# Computed diagnostic results can be saved in an app project or supplied to a
# report. Bind them to the fitted model, as for scree and fit-bootstrap
# results, so an analysis from an earlier structural fit cannot be presented
# beside a later calibration.
.dimensionality_test_result <- function(x, fit) {
  # Stamp the calculation. Results saved before this stamp reported a paired
  # t-test of the two subset means, which tests the targeting of the split
  # and not its dimensionality (see Details), so they are not restored beside
  # the current verdict: an analysis file drops them on restore and warns.
  x$algorithm <- .dimensionality_test_algorithm
  class(x) <- c("rasch_dimensionality_test", "list")
  attr(x, "fit_signature") <- .fit_boot_signature(fit)
  attr(x, "result_signature") <- .fit_boot_md5(x)
  x
}

.validate_dimensionality_test <- function(result, fit) {
  if (is.null(result)) return(invisible(NULL))
  signature <- attr(result, "result_signature")
  unsigned <- result
  attr(unsigned, "result_signature") <- NULL
  if (!inherits(result, "rasch_dimensionality_test") || !is.list(result) ||
      !is.character(signature) || length(signature) != 1L || is.na(signature) ||
      !.fit_boot_hash_matches(signature, unsigned) ||
      !.fit_boot_signature_matches(attr(result, "fit_signature"), fit))
    stop("`subtest` must be a dimensionality_test() result from this fitted model")
  # A result without the current stamp carries the superseded paired t-test
  # of the subset means, which its displays would present as an inferential
  # reading of the split. Refuse it here; .read_app_project() drops such a
  # result from a saved analysis with a warning before it reaches this point.
  if (!identical(result$algorithm, .dimensionality_test_algorithm))
    stop(paste("`subtest` uses a superseded person-subset comparison;",
               "recompute it with dimensionality_test()"),
         call. = FALSE)
  if (!is.null(result$bootstrap)) .require_refittable_calibration(fit)
  invisible(result)
}

#' @export
print.rasch_dimensionality_test <- function(x, ...) {
  cat("Person-subset unidimensionality t-test\n\n")
  if (!is.null(x$note)) {
    cat(x$note, "\n")
    return(invisible(x))
  }
  cat(sprintf("Split: %s\n", x$split))
  cat(sprintf("Significant comparisons: %.1f%% (95%% CI %.1f%% to %.1f%%; n = %d)\n",
              100 * x$prop_significant, 100 * x$ci[1], 100 * x$ci[2], x$n))
  verdict <- if (isTRUE(x$multidimensional))
    "evidence against unidimensionality" else
    if (identical(x$multidimensional, FALSE))
      "consistent with unidimensionality" else
        if (grepl("resolution insufficient", x$verdict_method %||% "",
                  fixed = TRUE))
          "withheld because the bootstrap resolution is insufficient" else
          "withheld for the data-driven split"
  cat("Verdict:", verdict, "\n")
  if (!is.null(x$p_boot))
    cat(sprintf("Bootstrap p: %s (%d of %d replicates used)\n",
                .fmt_p(x$p_boot), x$bootstrap$B_used, x$bootstrap$B))
  # the description of the split: reported, with no test attached to it
  smd <- x$subset_mean_difference
  if (is.list(smd))
    cat(sprintf(paste0("Mean difference between subset estimates: %.3f ",
                       "(SD %.3f) -- descriptive, not a test of ",
                       "unidimensionality\n"),
                smd$mean_difference, smd$sd_difference))
  if (!is.null(x$caution)) cat("Caution:", x$caution, "\n")
  invisible(x)
}

# The per-person comparison behind dimensionality_test(): each person
# estimated on the two subsets, persons extreme or inestimable on either
# subset set aside, and the standardised differences counted against the
# nominal level. Taking the response matrix and calibration rather than a
# fit object lets a bootstrap replicate compute exactly what the observed
# data did.
.dim_ttest <- function(X, tau_list, disc, pos, neg, alpha) {
  est_sub <- function(cols) {
    if (length(unique(disc[cols])) == 1L)
      .person_estimates(X[, cols, drop = FALSE], tau_list[cols],
                        disc = disc[cols][1])
    else .efrm_person_estimates(X[, cols, drop = FALSE], tau_list[cols],
                                disc[cols])
  }
  a <- est_sub(pos); b <- est_sub(neg)
  # A non-finite estimate or a non-positive uncertainty is not a person-level
  # test. In particular, allowing a zero denominator through would turn a
  # numerical failure into an infinite statistic and a certain rejection.
  usable <- is.finite(a$theta) & is.finite(b$theta) &
    is.finite(a$se) & a$se > 0 & is.finite(b$se) & b$se > 0
  ok <- usable & !a$extreme & !b$extreme
  dd <- a$theta[ok] - b$theta[ok]
  t <- dd / sqrt(a$se[ok]^2 + b$se[ok]^2)
  list(t = t, difference = dd, n = sum(ok),
       n_excluded_extreme = sum(usable) - sum(ok),
       n_sig = sum(abs(t) > qnorm(1 - alpha / 2)))
}

# The data-driven split: items grouped by the sign of their loading on a
# residual component. Computed from a residual matrix so that a bootstrap
# replicate repeats the selection on its own residuals. NULL when the
# component does not exist; a refusal from the decomposition propagates.
.dim_split <- function(Z, component) {
  pca <- .residual_pca_matrix(Z, n_components = component)
  cn <- paste0("PC", as.integer(component))
  if (!cn %in% names(pca$loadings_matrix)) return(NULL)
  ldg <- pca$loadings_matrix[[cn]][match(colnames(Z), pca$loadings_matrix$item)]
  pos <- which(ldg > 0)
  list(pos = pos, neg = setdiff(seq_len(ncol(Z)), pos),
       first_eigen = pca$first_eigen)
}

# The bootstrap generates from a single-facet Rasch model with a common
# unit and thresholds estimated directly, the same conditions as the item
# fit bootstrap, and for the same reasons.
.dim_bootstrap_check <- function(fit) {
  .require_refittable_calibration(fit)
  if (inherits(fit, c("rasch_efrm", "rasch_mfrm", "rasch_explanatory")))
    .refuse("the dimensionality bootstrap generates from a single-facet ",
            "Rasch model; an extended-frame, many-facet or explanatory fit ",
            "has a generating structure this function does not reproduce")
  if (!is.null(fit$disc) && length(unique(fit$disc)) > 1L)
    .refuse("the dimensionality bootstrap generates under equal ",
            "discriminations; this fit carries frame units that differ ",
            "across items")
  if (!is.null((fit$refit_spec %||% list())$pc_components))
    .refuse("the dimensionality bootstrap re-estimates each replicate the ",
            "way the fit was estimated; thresholds estimated through ",
            "principal components are not reproduced")
  invisible(TRUE)
}

# Whether a preparation note reports an anchored item rescored (the same
# "item <name> rescored" match rasch() makes, so anchor I1 does not trip over
# a note about I10)
.dim_anchor_rescored <- function(notes, anchors) {
  if (is.null(anchors) || !length(notes)) return(FALSE)
  a_names <- unique(as.character(anchors$item))
  any(vapply(notes, function(n)
    any(vapply(paste0("item ", a_names, " rescored"), grepl, TRUE,
               x = n, fixed = TRUE)), TRUE))
}

# Null distribution of the proportion of significant person comparisons
# under the fitted model. Every replicate keeps each person's raw score and
# missingness pattern (the score-conditional generator), refits the item
# calibration, and then does what the observed analysis did: a manual split
# is kept, a residual-component split is chosen afresh from the replicate's
# own residuals. The selection that inflates the observed proportion is
# thereby present in every replicate, and the comparison is like with like.
.dim_bootstrap <- function(fit, pos, neg, manual, component, alpha, B,
                           workers, seed) {
  .require_refittable_calibration(fit)
  workers <- min(as.integer(workers), .rasch_available_workers())
  .sim_seed_check()
  if (!is.null(seed)) {
    seed <- .check_whole(seed, "seed", 0)
    old <- .sim_seed_capture()
    on.exit(.sim_seed_restore(old), add = TRUE)
    set.seed(seed)
  }
  X <- fit$X
  na_mask <- if (anyNA(X)) is.na(X) else NULL
  spec <- fit$refit_spec %||% list()
  tau_list <- fit$tau_list; model <- fit$model; anchors <- spec$anchors
  maxit <- spec$maxit %||% 60L; tol <- spec$tol %||% 1e-8
  disc <- rep(1, ncol(X))
  # every random draw a replicate depends on is made here, in the parent,
  # so a worker reproduces its replicate from the seed it is handed and the
  # worker count cannot move the result
  seeds <- sample.int(.Machine$integer.max, B)
  one <- function(b) {
    old_stream <- .sim_seed_capture()
    on.exit(.sim_seed_restore(old_stream), add = TRUE)
    set.seed(seeds[b])
    Xb <- .fit_gen_conditional(X, tau_list, na_mask)
    # A replicate can leave a rare category unvisited. A refit on that shorter
    # scale is not the observed model and must not enter its null reference.
    # An item lost altogether, a shortened scale, or a rescored anchored item
    # therefore fails the replicate and is included in failure accounting.
    prep <- tryCatch(.prepare_X(Xb, model = model, anchors = anchors),
                     error = function(e) NULL)
    if (is.null(prep) || ncol(prep$X) != ncol(X) ||
        .dim_anchor_rescored(prep$notes, anchors))
      return(.fit_boot_failure("error"))
    Xb <- prep$X
    r <- .fit_refit_residuals(Xb, model, anchors, fit$m, maxit, tol)
    if (inherits(r, "rasch_fit_boot_failure")) return(r)
    split <- if (manual) list(pos = pos, neg = neg) else
      tryCatch(.dim_split(r$Z, component), error = function(e) NULL)
    if (is.null(split) || length(split$pos) < 2L || length(split$neg) < 2L)
      return(.fit_boot_failure("error"))
    tt <- .dim_ttest(Xb, r$tau_list, disc, split$pos, split$neg, alpha)
    if (tt$n < 10L) return(.fit_boot_failure("error"))
    tt$n_sig / tt$n
  }
  reps <- .rasch_boot_apply(B, one, workers = workers,
                            label = "dimensionality bootstrap")
  status <- vapply(reps, .fit_boot_status, "")
  keep <- status == "ok"
  B_used <- sum(keep)
  min_success <- .fit_min_boot_success(B)
  if (B_used < min_success)
    .fit_boot_refuse(
      "only ", B_used, " of ", B, " bootstrap replicates were usable (",
      sum(status == "nonconverged"), " did not converge; ",
      sum(status == "error"), " otherwise failed); at least ", min_success,
      " are required for the bootstrap null. The fitted model generates ",
      "data this estimator cannot fit reliably",
      B = B, B_used = B_used,
      B_nonconverged = sum(status == "nonconverged"),
      B_errors = sum(status == "error"))
  if (B_used < 0.9 * B)
    warning(B - B_used, " of ", B, " bootstrap replicates were unusable (",
            sum(status == "nonconverged"), " did not converge; ",
            sum(status == "error"), " otherwise failed); the null is formed ",
            "from the remaining ", B_used,
            "; replicates are not lost at random, so read the result with ",
            "that in mind", call. = FALSE)
  list(null = unlist(reps[keep]), B = B, B_used = B_used,
       B_nonconverged = sum(status == "nonconverged"),
       B_errors = sum(status == "error"), minimum_usable = min_success,
       seed = seed)
}

#' Magnitude of multidimensionality from a subtest analysis
#'
#' Estimates how strongly two or more hypothesised subscales measure
#' distinct traits, by Andrich's (2016) comparison of two reliability
#' calculations: one treating all items as independent (which inflates
#' reliability under multidimensionality) and one on the subtest analysis
#' in which each subscale is combined into a single polytomous super-item
#' (which absorbs the unique subscale variance). Under the bifactor
#' formalisation \eqn{\beta_{ns} = \beta_n + c\,\beta'_{ns}} (Marais and
#' Andrich 2008), with \eqn{S} subscales of \eqn{K} items,
#' \deqn{c^2 = S\,(r_1/r_2 - 1) \frac{SK - 1}{S(K - 1)},}
#' the latent correlation between subscales is \eqn{\rho = 1/(1 + c^2)},
#' and \eqn{A = S/(S + c^2)} is the proportion of common (non-unique,
#' non-error) variance. Both the person separation index and coefficient
#' alpha versions are reported (Andrich and Marais 2019, ch. 24). Both the
#' original and subtest calibrations must converge.
#'
#' Both reliability calculations use rows with responses to every item.
#' PSI further requires a finite person estimate and standard error in both
#' fits. The calibrations are retained; only the reliability sample changes.
#' With missing responses, the result describes this complete-response sample,
#' not necessarily the full population. The table reports rows used and excluded.
#'
#' @param fit A fitted object from \code{\link{rasch}}.
#' @param subtests A list of character vectors assigning \emph{every} item
#'   of the fit to one subscale (at least two subscales of two or more
#'   items). The published magnitude formula requires equal subscale sizes.
#' @return A list of class \code{"rasch_dim_magnitude"}: the comparison
#'   \code{table} (rows PSI and alpha; columns \code{run1}, \code{subtest},
#'   \code{c2}, \code{c}, \code{rho}, \code{A}, \code{n},
#'   \code{n_excluded}), the subtest \code{refit},
#'   and the design constants \code{S} and \code{K}.
#' @references Andrich, D. (2016). Components of variance of scales with a
#'   bifactor subscale structure from two calculations of alpha.
#'   Educational Measurement: Issues and Practice, 35(4), 25-30.
#' @examples
#' set.seed(1); N <- 500
#' common <- rnorm(N); u1 <- rnorm(N); u2 <- rnorm(N)
#' d <- rep(seq(-1, 1, length.out = 5), 2)
#' X <- sapply(1:10, function(i) rbinom(N, 1,
#'   plogis(common + 0.8 * (if (i <= 5) u1 else u2) - d[i])))
#' colnames(X) <- paste0("I", 1:10)
#' fit <- rasch(X)
#' dimensionality_magnitude(fit,
#'   list(paste0("I", 1:5), paste0("I", 6:10)))$table
#' @export
dimensionality_magnitude <- function(fit, subtests) {
  if (!inherits(fit, "rasch")) stop("dimensionality_magnitude needs a rasch fit")
  if (!isTRUE(fit$est$converged))
    stop("the fitted calibration did not converge; dimensionality magnitude is unavailable")
  if (!is.list(subtests) || length(subtests) < 2)
    stop("supply a list of at least two subscales")
  if (!all(vapply(subtests, function(x)
      (is.character(x) || is.factor(x)) && is.null(dim(x)) && length(x) &&
        !anyNA(x) && all(nzchar(trimws(as.character(x)))), logical(1))))
    stop("each subscale must be an ordinary vector of non-missing item names")
  subtests <- lapply(subtests, as.character)
  allit <- unlist(subtests)
  if (!setequal(allit, fit$items$item) || anyDuplicated(allit))
    stop("subtests must assign every item of the fit to exactly one subscale")
  if (length(unique(lengths(subtests))) != 1L)
    stop("the Andrich (2016) magnitude formula requires equal subscale sizes")
  S <- length(subtests)
  K <- lengths(subtests)[1L]
  g <- S * (K - 1) / (S * K - 1)
  # the ratio compares the same calibration with and without the subscale
  # structure. Under an explanatory design the subtest refit relaxes every
  # superitem, so the ratio would measure the design restriction instead
  if (inherits(fit, "rasch_explanatory"))
    stop("the dimensionality magnitude is not defined for an explanatory ",
         "calibration: every subtest is freed in the subtest refit, so the ",
         "reliability ratio would measure the design restriction rather ",
         "than the subscale structure. Use the unrestricted calibration")
  refit <- combine_items(fit, subtests)
  if (!isTRUE(refit$est$converged))
    stop("the subtest calibration did not converge; dimensionality magnitude is unavailable")
  # A partially answered subscale becomes missing when summed. Comparing the
  # full original PSI with that refit's PSI would change both the response
  # coverage and the person distribution, not just the subscale structure.
  p <- fit$person
  q <- refit$person
  if (nrow(fit$X) != nrow(refit$X) || nrow(p) != nrow(fit$X) ||
      nrow(q) != nrow(refit$X) ||
      !identical(as.character(p$id), as.character(q$id)))
    stop("the original and subtest response rows are not aligned")
  complete <- stats::complete.cases(fit$X) & stats::complete.cases(refit$X)
  keep <- complete & is.finite(p$theta) & is.finite(p$se) & p$se >= 0 &
    is.finite(q$theta) & is.finite(q$se) & q$se >= 0
  psi1 <- .psi(p$theta, p$se, keep = keep)
  psi2 <- .psi(q$theta, q$se, keep = keep)
  ratio <- function(r1, r2) {
    # The variance decomposition is a ratio of reliabilities. A zero or
    # negative coefficient is not a boundary estimate of c^2; it means the
    # reliability decomposition itself is unavailable. Truncating the
    # resulting negative ratio to c^2 = 0 would falsely report perfectly
    # correlated subscales from an unusable reliability estimate.
    if (any(!is.finite(c(r1, r2))) || r1 <= 0 || r2 <= 0)
      return(rep(NA_real_, 4))
    c2 <- max(S * (r1 / r2 - 1) / g, 0)
    c(c2, sqrt(c2), 1 / (1 + c2), S / (S + c2))
  }
  psi_row <- ratio(psi1$PSI, psi2$PSI)
  alp_row <- ratio(fit$alpha$alpha, refit$alpha$alpha)
  tab <- data.frame(index = c("PSI", "alpha"),
                    run1 = c(psi1$PSI, fit$alpha$alpha),
                    subtest = c(psi2$PSI, refit$alpha$alpha),
                    c2 = c(psi_row[1], alp_row[1]),
                    c = c(psi_row[2], alp_row[2]),
                    rho = c(psi_row[3], alp_row[3]),
                    A = c(psi_row[4], alp_row[4]),
                    n = c(sum(keep), sum(complete)),
                    n_excluded = c(sum(!keep), sum(!complete)))
  out <- list(table = tab, refit = refit, S = S, K = K,
              algorithm = "complete-panel-1",
              alpha_applicable = fit$alpha$applicable)
  out <- .tag_tables(out)
  class(out) <- "rasch_dim_magnitude"
  out
}

#' @export
print.rasch_dim_magnitude <- function(x, ...) {
  cat(sprintf("Magnitude of multidimensionality (Andrich 2016): %d subscales, %.0f items each\n",
              x$S, x$K))
  tab <- x$table
  num <- vapply(tab, is.numeric, TRUE)
  tab[num] <- lapply(tab[num], round, 3)
  print(tab, row.names = FALSE)
  cat("rho = latent correlation between subscales; A = proportion of common variance\n")
  if (any(x$table$n_excluded > 0))
    cat("Reliabilities use matched complete-response rows; PSI also requires usable estimates in both fits.\n")
  invisible(x)
}
