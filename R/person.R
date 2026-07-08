# rasch :: person estimation
# ===========================================================================
# Warm (1989) weighted likelihood person estimates. The WLE is finite at the
# extreme (zero and maximum) raw scores, so extreme persons receive usable
# locations; they are still flagged so reliability and test-of-fit statistics
# can exclude them. Persons with missing responses are
# estimated on their own observed item subset, grouped by missing-data
# pattern for speed.
# ===========================================================================

#' Category-score moments for a polytomous item
#'
#' @param theta Person location, in logits.
#' @param tau_i Numeric vector of the item's threshold parameters.
#' @param disc Discrimination (frame unit) multiplier on the exponent; 1 for
#'   the ordinary Rasch model.
#' @return A list with category probabilities \code{P}, expected score \code{E},
#'   variance \code{V}, and third and fourth central moments \code{mu3},
#'   \code{mu4}.
#' @examples
#' item_moments(0.5, c(-1, 0, 1))
#' @export
item_moments <- function(theta, tau_i, disc = 1) {
  m <- length(tau_i); x <- 0:m
  num <- exp(disc * (x * theta - c(0, cumsum(tau_i)))); P <- num / sum(num)
  E <- sum(x * P); d <- x - E
  list(P = P, E = E, V = sum(d^2 * P), mu3 = sum(d^3 * P), mu4 = sum(d^4 * P))
}

#' Warm's weighted likelihood estimates by raw score
#'
#' Computes the weighted likelihood estimate (WLE) of person location for every
#' possible raw score on a set of items, with standard errors. WLE estimates
#' are finite at the extreme (zero and maximum) scores, unlike the maximum
#' likelihood estimate.
#'
#' @param tau_list List of per-item threshold vectors.
#' @param disc Common discrimination (frame unit) of the items; with a
#'   constant discrimination the raw score remains sufficient.
#' @return A list with \code{theta} and \code{se}, each named by raw score.
#' @examples
#' person_wle(list(c(-1, 0), c(-0.5, 0.5), c(0, 1)))
#' @export
person_wle <- function(tau_list, disc = 1) {
  Smax <- sum(vapply(tau_list, length, 1L))
  theta <- se <- setNames(rep(NA_real_, Smax + 1L), as.character(0:Smax))
  for (R in 0:Smax) {
    g <- function(th) {
      mo <- lapply(tau_list, item_moments, theta = th, disc = disc)
      E  <- sum(vapply(mo, `[[`, 0, "E"));  V <- sum(vapply(mo, `[[`, 0, "V"))
      m3 <- sum(vapply(mo, `[[`, 0, "mu3"))
      (R - E) + disc * m3 / (2 * V)
    }
    root <- tryCatch(uniroot(g, c(-30, 30), tol = 1e-9)$root, error = function(e) NA_real_)
    theta[as.character(R)] <- root
    if (!is.na(root)) {
      V <- sum(vapply(lapply(tau_list, item_moments, theta = root, disc = disc),
                      `[[`, 0, "V"))
      se[as.character(R)] <- 1 / sqrt(disc^2 * V)      # WLE SE ~ 1/sqrt(information)
    }
  }
  list(theta = theta, se = se)
}

# Person locations for an arbitrary response matrix, grouped by
# missing-data pattern so each WLE score table is solved once. disc is a
# single common discrimination (frame unit); raw scores stay sufficient.
.person_estimates <- function(X, tau_list, disc = 1) {
  N <- nrow(X)
  obs <- !is.na(X)
  m <- vapply(tau_list, length, 1L)
  pat <- apply(obs, 1, function(z) paste(which(z), collapse = ","))
  theta <- se <- rep(NA_real_, N)
  raw <- rowSums(X, na.rm = TRUE); raw[rowSums(obs) == 0L] <- NA
  max_raw <- as.numeric(obs %*% m)
  for (key in unique(pat)) {
    cols <- as.integer(strsplit(key, ",", fixed = TRUE)[[1]])
    if (!length(cols)) next
    sel <- which(pat == key)
    pe <- person_wle(tau_list[cols], disc = disc)
    r <- rowSums(X[sel, cols, drop = FALSE])
    theta[sel] <- pe$theta[as.character(r)]
    se[sel]    <- pe$se[as.character(r)]
  }
  data.frame(n_items = rowSums(obs), raw = raw, max_raw = max_raw,
             theta = theta, se = se,
             extreme = !is.na(raw) & (raw == 0L | raw == max_raw))
}

# Model moments evaluated at each person's location, observed cells only.
# Persons sharing a location (same pattern and raw score) share a row of the
# unique-theta moment tables. disc may be a per-column vector (frame units).
.moment_arrays <- function(theta, tau_list, disc = NULL) {
  L <- length(tau_list)
  if (is.null(disc)) disc <- rep(1, L)
  if (length(disc) == 1L) disc <- rep(disc, L)
  ut <- sort(unique(theta[!is.na(theta)]))
  E <- V <- M3 <- M4 <- matrix(NA_real_, length(ut), L)
  for (u in seq_along(ut)) {
    mo <- lapply(seq_len(L), function(i)
      item_moments(ut[u], tau_list[[i]], disc = disc[i]))
    E[u, ]  <- vapply(mo, `[[`, 0, "E")
    V[u, ]  <- vapply(mo, `[[`, 0, "V")
    M3[u, ] <- vapply(mo, `[[`, 0, "mu3")
    M4[u, ] <- vapply(mo, `[[`, 0, "mu4")
  }
  idx <- match(theta, ut)
  list(E = E[idx, , drop = FALSE], V = V[idx, , drop = FALSE],
       M3 = M3[idx, , drop = FALSE], M4 = M4[idx, , drop = FALSE])
}

#' Raw score to measure conversion table
#'
#' The score-to-logit conversion for complete responders: every possible raw
#' score with its location, standard error, and the frequency and cumulative
#' percentage of complete responders at that score (the complete-data
#' estimates table of Andrich and Marais 2019, ch. 10).
#'
#' Two estimators are available. \code{"wle"} (the default) is Warm's
#' weighted likelihood estimate, finite at the extreme scores. \code{"mle"}
#' is the plain maximum likelihood estimate, infinite at the
#' extremes. \code{extremes = "extrapolated"} replaces the extreme-score
#' entries by the geometric extrapolation described in Andrich and Marais
#' (2019, ch. 10): successive score-to-score
#' differences grow towards the extremes, so the last difference is
#' continued geometrically -- the extrapolated top difference \eqn{d} solves
#' \eqn{b = \sqrt{a d}} where \eqn{a, b} are the two preceding differences
#' (equivalently \eqn{d = b^2/a}), and symmetrically at zero. The standard
#' error at an extrapolated location is \eqn{1/\sqrt{I(\theta)}} evaluated
#' there. With \code{method = "wle"} the extrapolation replaces the finite
#' Warm estimates at the extremes, giving the extrapolated form of the
#' conversion table from a WLE analysis.
#'
#' @param fit A fitted object from \code{\link{rasch}}.
#' @param method \code{"wle"} (Warm, default) or \code{"mle"}.
#' @param extremes \code{"model"} keeps the estimator's own extreme-score
#'   values (\code{NA} for MLE); \code{"extrapolated"} applies the geometric
#'   extrapolation.
#' @return A data frame with \code{score}, \code{theta}, \code{se},
#'   \code{freq}, \code{cum_pct} (omitted when no complete responders
#'   exist), and \code{extrapolated}; \code{NULL} for
#'   fits without a common raw-score metric (EFRM).
#' @examples
#' set.seed(1)
#' d <- seq(-1.5, 1.5, length.out = 6)
#' X <- matrix(rbinom(300 * 6, 1, plogis(outer(rnorm(300), d, "-"))), 300, 6)
#' colnames(X) <- paste0("I", 1:6)
#' score_table(rasch(X), method = "mle", extremes = "extrapolated")
#' @export
score_table <- function(fit, method = c("wle", "mle"),
                        extremes = c("model", "extrapolated")) {
  method <- match.arg(method); extremes <- match.arg(extremes)
  if (is.null(fit$score_table)) return(NULL)
  tab <- fit$score_table[, c("score", "theta", "se")]
  M <- max(tab$score)
  disc <- if (is.null(fit$disc)) 1 else fit$disc[1]
  info <- function(th) sum(vapply(fit$tau_list, function(tt)
    disc^2 * item_moments(th, tt, disc = disc)$V, 0))
  if (method == "mle") {
    for (r in seq_len(M - 1)) {
      tab$theta[tab$score == r] <- uniroot(function(th)
        r - sum(vapply(fit$tau_list, function(tt)
          disc * item_moments(th, tt, disc = disc)$E, 0)),
        c(-30, 30), tol = 1e-9)$root
    }
    tab$theta[c(1, M + 1)] <- NA_real_
    tab$se <- 1 / sqrt(vapply(tab$theta, function(th)
      if (is.na(th)) NA_real_ else info(th), 0))
  }
  tab$extrapolated <- FALSE
  if (extremes == "extrapolated") {
    if (M < 4) stop("extrapolation needs at least three interior scores")
    th <- tab$theta
    lo <- th[2] - (th[3] - th[2])^2 / (th[4] - th[3])
    hi <- th[M] + (th[M] - th[M - 1])^2 / (th[M - 1] - th[M - 2])
    tab$theta[c(1, M + 1)] <- c(lo, hi)
    tab$se[c(1, M + 1)] <- 1 / sqrt(c(info(lo), info(hi)))
    tab$extrapolated[c(1, M + 1)] <- TRUE
  }
  raw <- rowSums(fit$X)
  freq <- as.integer(table(factor(raw[stats::complete.cases(fit$X)],
                                  levels = 0:M)))
  if (sum(freq) > 0) {
    tab$freq <- freq
    tab$cum_pct <- 100 * cumsum(freq) / sum(freq)
  }
  tab
}

#' Person measures with extrapolated extreme scores
#'
#' Extreme persons (zero or maximum raw score on their observed items) are
#' excluded from calibration, but they cannot be left out of group
#' comparisons; Andrich and Marais (2019, ch. 10) therefore describe an
#' extrapolated measure for them,
#' continuing the growth of the score-to-score differences so the last
#' difference is the geometric mean of its neighbours (see
#' \code{\link{score_table}}). This helper applies the same rule to the
#' person table: for each missing-data pattern with extreme persons, the
#' score-to-measure conversion over that pattern's items is extrapolated at
#' its ends, and the extreme persons receive the extrapolated location with
#' the standard error \eqn{1/\sqrt{I(\theta)}} evaluated there. Non-extreme
#' persons keep their estimates unchanged. The extrapolation continues the
#' Warm (weighted likelihood) conversion, matching the package's person
#' estimates.
#'
#' @param fit A fitted object from \code{\link{rasch}} (equal
#'   discriminations; not EFRM).
#' @return The fit's person table with two added columns,
#'   \code{theta_extrapolated} and \code{se_extrapolated}: equal to
#'   \code{theta} and \code{se} for non-extreme persons, extrapolated for
#'   extreme persons. Patterns with fewer than three interior scores cannot
#'   be extrapolated and keep their Warm values.
#' @examples
#' set.seed(1)
#' d <- seq(-2, 2, length.out = 8)
#' X <- matrix(rbinom(300 * 8, 1, plogis(outer(rnorm(300, 0, 2), d, "-"))), 300, 8)
#' colnames(X) <- paste0("I", 1:8)
#' fit <- rasch(X)
#' pe <- person_extrapolated(fit)
#' head(pe[pe$extreme, c("theta", "theta_extrapolated", "se", "se_extrapolated")])
#' @export
person_extrapolated <- function(fit) {
  if (!is.null(fit$disc) && length(unique(fit$disc)) > 1L)
    stop("extrapolation over a common raw score needs equal discriminations; ",
         "EFRM fits report per-group score curves instead")
  p <- fit$person
  p$theta_extrapolated <- p$theta
  p$se_extrapolated <- p$se
  ext <- which(p$extreme)
  if (!length(ext)) return(p)
  X <- fit$X
  pat <- apply(!is.na(X), 1, paste, collapse = "")
  for (pt in unique(pat[ext])) {
    rows <- intersect(which(pat == pt), ext)
    obs <- which(!is.na(X[rows[1], ]))
    if (length(obs) < 2) next
    tl <- fit$tau_list[obs]
    pe <- person_wle(tl)
    th <- unname(pe$theta); M <- length(th) - 1L
    if (M < 4) next                        # too few interior scores
    lo <- th[2] - (th[3] - th[2])^2 / (th[4] - th[3])
    hi <- th[M] + (th[M] - th[M - 1])^2 / (th[M - 1] - th[M - 2])
    info <- function(t) sum(vapply(tl, function(tt)
      item_moments(t, tt)$V, 0))
    for (r in rows) {
      raw <- sum(X[r, obs])
      if (raw == 0L) {
        p$theta_extrapolated[r] <- lo
        p$se_extrapolated[r] <- 1 / sqrt(info(lo))
      } else {
        p$theta_extrapolated[r] <- hi
        p$se_extrapolated[r] <- 1 / sqrt(info(hi))
      }
    }
  }
  p
}
