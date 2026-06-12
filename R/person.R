# RaschR :: person estimation
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
#' score with its WLE location and standard error.
#'
#' @param fit A fitted object from \code{\link{rasch}}.
#' @return A data frame with \code{score}, \code{theta}, and \code{se}.
#' @examples
#' set.seed(1)
#' d <- seq(-1.5, 1.5, length.out = 6)
#' X <- matrix(rbinom(300 * 6, 1, plogis(outer(rnorm(300), d, "-"))), 300, 6)
#' colnames(X) <- paste0("I", 1:6)
#' score_table(rasch(X))
#' @export
score_table <- function(fit) fit$score_table
