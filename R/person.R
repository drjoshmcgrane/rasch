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
  if (length(theta) != 1L || !is.numeric(theta) || is.complex(theta) ||
      !is.null(dim(theta)) || !is.null(oldClass(theta)) || !is.finite(theta))
    stop("`theta` must be one finite location")
  if (!is.numeric(tau_i) || is.complex(tau_i) || !is.null(dim(tau_i)) ||
      !is.null(oldClass(tau_i)) || !length(tau_i) || any(!is.finite(tau_i)))
    stop("`tau_i` must be a non-empty vector of finite thresholds")
  if (length(disc) != 1L || !is.numeric(disc) || is.complex(disc) ||
      !is.null(dim(disc)) || !is.null(oldClass(disc)) ||
      !is.finite(disc) || disc <= 0)
    stop("`disc` must be one positive finite discrimination")
  m <- length(tau_i); x <- 0:m
  cumulative <- c(0, cumsum(tau_i))
  lp <- disc * (x * theta - cumulative)
  if (any(!is.finite(cumulative)) || any(!is.finite(lp)))
    stop("the threshold, location and discrimination combination is outside ",
         "the numerically representable range", call. = FALSE)
  num <- exp(lp - max(lp)); P <- num / sum(num)   # log-sum-exp: no overflow
  E <- sum(x * P); d <- x - E
  list(P = P, E = E, V = sum(d^2 * P), mu3 = sum(d^3 * P), mu4 = sum(d^4 * P))
}

# Root searches for person locations must follow the calibration origin.
# A fixed [-30, 30] interval makes an otherwise valid externally anchored
# scale fail as soon as its origin is translated outside that range. Moving
# thirty inverse-discrimination logits beyond the most extreme threshold
# puts both ends well into the response-function tails while preserving exact
# translation equivariance.
.person_root_interval <- function(tau_list, disc) {
  if (length(disc) == 1L) disc <- rep(disc, length(tau_list))
  if (length(disc) != length(tau_list) || any(!is.finite(disc)) ||
      any(disc <= 0))
    stop("the model discriminations are unavailable for person scoring",
         call. = FALSE)
  lower <- min(vapply(seq_along(tau_list), function(i)
    min(tau_list[[i]]) - 30 / disc[i], 0))
  upper <- max(vapply(seq_along(tau_list), function(i)
    max(tau_list[[i]]) + 30 / disc[i], 0))
  if (!is.finite(lower) || !is.finite(upper) || lower >= upper)
    stop("the threshold scale is too large for a stable person-location root search",
         call. = FALSE)
  c(lower, upper)
}

# Evaluate a common-unit scoring SE without squaring the discrimination.
# The equivalent 1 / sqrt(disc^2 * V) can overflow or underflow even when
# both the response probabilities and the final SE are representable.
.common_person_se <- function(theta, tau_list, disc) {
  V <- sum(vapply(tau_list, function(tt)
    item_moments(theta, tt, disc = disc)$V, 0))
  (1 / sqrt(V)) / disc
}

# Search on the response-function scale. The curve is shared across all
# sufficient scores for an administered item pattern. Evaluating its moments
# in batches keeps a multiple-root search inexpensive on ordinary banks.
.person_wle_curve <- function(tau_list, disc, weights = NULL) {
  disc <- rep_len(disc, length(tau_list))
  scale <- max(disc)
  origin <- mean(range(unlist(tau_list)))
  rp <- disc / scale
  q <- if (is.null(weights)) rep(1, length(disc)) else weights / max(weights)
  interval <- (.person_root_interval(tau_list, disc) - origin) * scale
  if (any(!is.finite(interval)) || any(rp == 0))
    stop("the range of model units is too large for stable person scoring")
  ct <- lapply(seq_along(disc), function(j)
    c(0, cumsum((tau_list[[j]] - origin) * disc[j])))
  evaluate <- function(u) {
    E <- H <- J <- Hprime <- psi <- numeric(length(u))
    for (j in seq_along(disc)) {
      k <- seq_along(ct[[j]]) - 1L
      lp <- outer(u, k * rp[j])
      lp <- sweep(lp, 2L, ct[[j]], "-")
      top <- apply(lp, 1L, max)
      P <- exp(lp - top)
      den <- rowSums(P); P <- P / den
      ej <- drop(P %*% k)
      dev <- matrix(k, length(u), length(k), byrow = TRUE) - ej
      vj <- rowSums(P * dev^2)
      mj <- rowSums(P * dev^3)
      E <- E + q[j] * rp[j] * ej
      H <- H + q[j] * rp[j]^2 * vj
      J <- J + q[j]^2 * rp[j]^2 * vj
      Hprime <- Hprime + q[j] * rp[j]^3 * mj
      psi <- psi + q[j] * (top + log(den))
    }
    cbind(base = E - (J / H) * (Hprime / H) / 2,
          psi = psi, info = H)
  }
  if (diff(interval) <= 160) {
    grid <- seq(interval[1], interval[2], length.out =
      max(3L, ceiling(diff(interval) / .25) + 1L))
  } else {
    # Long gaps need no enormous regular mesh. Include every category-line
    # crossing, not just the adjacent thresholds: disordered thresholds can
    # make non-adjacent categories the dominant pair.
    grid <- interval
    for (j in seq_along(ct)) {
      n <- length(ct[[j]])
      pairs <- utils::combn(seq_len(n), 2L)
      crossings <- (ct[[j]][pairs[2, ]] - ct[[j]][pairs[1, ]]) /
        ((pairs[2, ] - pairs[1, ]) * rp[j])
      offsets <- c(-30, -16, seq(-8, 8, by = .25), 16, 30) / rp[j]
      grid <- c(grid, as.vector(outer(crossings, offsets, "+")))
    }
    grid <- sort(unique(grid[is.finite(grid) & grid >= interval[1] &
                               grid <= interval[2]]))
  }
  list(evaluate = evaluate, grid = grid, values = evaluate(grid),
       origin = origin, scale = scale, ordinary = all(q == 1))
}

.person_wle_maximum <- function(curve, score) {
  g <- score - curve$values[, "base"]
  # A maximum crosses from positive to negative. Include exact grid roots
  # too; comparing objectives also rejects any stationary minimum among them.
  a <- seq_len(length(g) - 1L)
  cross <- a[is.finite(g[a]) & is.finite(g[a + 1L]) &
               g[a] > 0 & g[a + 1L] < 0]
  gf <- function(u) score - curve$evaluate(u)[, "base"]
  roots <- c(curve$grid[is.finite(g) & g == 0],
    vapply(cross, function(i) tryCatch(
      stats::uniroot(gf, curve$grid[c(i, i + 1L)], tol = 1e-10)$root,
      error = function(e) NA_real_), 0))
  roots <- sort(unique(roots[is.finite(roots)]))
  if (!length(roots)) return(NA_real_)
  if (curve$ordinary) {
    v <- curve$evaluate(roots)
    objective <- score * roots - v[, "psi"] + log(v[, "info"]) / 2
  } else {
    # With externally imposed weights the Warm correction is not generally
    # one half log information. Its antiderivative, rather than that ordinary
    # WLE objective, ranks competing roots of the weighted equation.
    objective <- vapply(roots, function(u) {
      if (u == roots[1L]) return(0)
      tryCatch(stats::integrate(gf, roots[1L], u, rel.tol = 1e-9,
        subdivisions = 1000L)$value, error = function(e) NA_real_)
    }, 0)
  }
  if (any(!is.finite(objective))) return(NA_real_)
  best <- max(objective)
  # Equal maxima have no unique WLE. Use the lower maximizer consistently;
  # never average maxima, since their midpoint can be a likelihood minimum.
  take <- which(best - objective <= 1e-10)[1L]
  curve$origin + roots[take] / curve$scale
}

#' Warm's weighted likelihood estimates by raw score
#'
#' Computes the weighted likelihood estimate (WLE) of person location for every
#' possible raw score on a set of items, with standard errors. WLE estimates
#' are finite at the extreme (zero and maximum) scores, unlike the maximum
#' likelihood estimate.
#'
#' @details
#' For raw score \eqn{R}, let \eqn{E(\theta)}, \eqn{V(\theta)}, and
#' \eqn{\mu_3(\theta)} be the sums of the item expected scores, variances, and
#' third central moments. The estimate solves Warm's weighted score equation
#' \deqn{R-E(\theta)+\frac{\mu_3(\theta)}{2V(\theta)}=0.}
#' When the equation has several solutions, competing maxima are compared
#' using log likelihood plus one half log information. Equal maxima use the
#' lower location; their average need not maximize the weighted likelihood.
#' With common discrimination \eqn{d}, its explicit multiplier cancels from
#' this equation, although the moments are evaluated under \eqn{d}. The
#' reported standard error is
#' \deqn{\operatorname{SE}(\hat{\theta})=
#' \{d^2V(\hat{\theta})\}^{-1/2}.}
#'
#' @param tau_list List of per-item threshold vectors.
#' @param disc Common discrimination (frame unit) of the items; with a
#'   constant discrimination the raw score remains sufficient.
#' @return A list with \code{theta} and \code{se}, each named by raw score.
#' @references
#' Warm, T. A. (1989). Weighted likelihood estimation of ability in item
#' response theory. Psychometrika, 54(3), 427--450.
#'
#' Zhang, J. (2005). Bias correction for the maximum likelihood estimate of
#' ability. ETS Research Report RR-05-15.
#' @seealso \code{\link{score_table}} and \code{\link{person_extrapolated}}.
#' @examples
#' person_wle(list(c(-1, 0), c(-0.5, 0.5), c(0, 1)))
#' @export
person_wle <- function(tau_list, disc = 1) {
  if (!is.list(tau_list) || !length(tau_list) ||
      !all(vapply(tau_list, function(t)
        is.numeric(t) && !is.complex(t) && is.null(dim(t)) &&
          is.null(oldClass(t)) && length(t) >= 1L && all(is.finite(t)),
        TRUE)))
    stop("`tau_list` must be a list of non-empty finite threshold vectors")
  if (length(disc) != 1L || !is.numeric(disc) || is.complex(disc) ||
      !is.null(dim(disc)) || !is.null(oldClass(disc)) ||
      !is.finite(disc) || disc <= 0)
    stop("`disc` must be one positive finite discrimination: the raw-score ",
         "WLE requires a common discrimination")
  Smax <- sum(vapply(tau_list, length, 1L))
  theta <- se <- setNames(rep(NA_real_, Smax + 1L), as.character(0:Smax))
  curve <- .person_wle_curve(tau_list, disc)
  for (R in 0:Smax) {
    root <- .person_wle_maximum(curve, R)
    theta[as.character(R)] <- root
    if (!is.na(root)) {
      se[as.character(R)] <- .common_person_se(root, tau_list, disc)
    }
  }
  list(theta = theta, se = se)
}

#' Person estimates with externally imposed weights
#'
#' Calculates a second set of person measures after assigning relative
#' weights to items or item sets. The fitted calibration is not changed.
#' In particular, these estimates do not replace the ordinary Rasch person
#' measures used for fit, reliability, targeting or DIF.
#'
#' @details
#' Let \eqn{q_i} be the external weight and \eqn{a_i} the model unit for
#' response cell \eqn{i}. Write
#' \eqn{H(\theta)=\sum_i q_i a_i^2V_i(\theta)} and
#' \eqn{J(\theta)=\sum_i q_i^2a_i^2V_i(\theta)}. The estimate solves the
#' externally weighted Warm score equation
#' \deqn{\sum_i q_i a_i\{x_i-E_i(\theta)\}+
#' \frac{J(\theta)\sum_i q_i a_i^3\mu_{3i}(\theta)}
#' {2H(\theta)^2}=0.}
#' Competing maxima are ranked by the integral of this estimating score.
#' With equal weights this reduces to the ordinary weighted log likelihood;
#' unequal external weights require their own correction. Equal maxima use
#' the lower location.
#' Its standard error is the sandwich form
#' \deqn{\operatorname{SE}(\hat\theta)=
#' \frac{\{\sum_iq_i^2a_i^2V_i(\hat\theta)\}^{1/2}}
#' {\sum_iq_ia_i^2V_i(\hat\theta)}.}
#' This matters because an external weight changes the estimating equation;
#' it does not create independent replications of an item. With equal weights
#' the equation and standard error reduce to the person estimates in
#' \code{fit$person}. Weights are normalised to mean one over the fitted
#' response cells, so only their relative values matter. A zero weight omits
#' that item or set from the secondary measure. Positive weights must be
#' numerically representable relative to the largest supplied weight.
#'
#' For an MFRM fit, an item weight applies to all response cells belonging to
#' that item. For an EFRM fit, \code{by = "set"} uses the fitted item-set map
#' unless \code{sets} is supplied.
#'
#' @param fit A fitted object from \code{\link{rasch}},
#'   \code{\link{rasch_mfrm}} or \code{\link{rasch_efrm}}.
#' @param weights A named numeric vector of non-negative finite relative
#'   weights. Names must identify every fitted item when \code{by = "item"},
#'   or every item set when \code{by = "set"}.
#' @param by Whether \code{weights} names individual \code{"item"}s or
#'   \code{"set"}s.
#' @param sets For \code{by = "set"}, a named character vector mapping items
#'   to sets, or a named list whose elements contain item names. It can be
#'   omitted for an EFRM fit, which already contains this map.
#' @return A data frame with the person identifiers and factors, observed
#'   response-cell count (item count for an ordinary fit), raw and externally
#'   weighted scores, maximum scores, weighted location, sandwich standard
#'   error and extreme-score flag. The resolved
#'   item weights are retained in the \code{"weighting"} attribute.
#' @references
#' Warm, T. A. (1989). Weighted likelihood estimation of ability in item
#' response theory. Psychometrika, 54(3), 427--450.
#' @examples
#' set.seed(1)
#' d <- simulate_rasch(n_persons = 200, n_items = 6)
#' fit <- rasch(d)
#' weighted_person_estimates(
#'   fit, setNames(c(2, 2, 1, 1, 0.5, 0.5), colnames(fit$X)))
#' @export
weighted_person_estimates <- function(fit, weights,
                                      by = c("item", "set"), sets = NULL) {
  if (!inherits(fit, "rasch") || inherits(fit, "rasch_btl"))
    stop("`fit` must be an ordinary, multiple-ratings or extended-frames Rasch fit",
         call. = FALSE)
  if (!isTRUE(fit$est$converged))
    stop("the fitted calibration did not converge; weighted person estimates are unavailable",
         call. = FALSE)
  if (!.efrm_link_converged(fit))
    stop("the fitted set-unit link did not converge; weighted person estimates are unavailable",
         call. = FALSE)
  by <- match.arg(by)
  # Set names are substantive labels and are canonicalised. Item names are
  # selectors into the fitted response matrix: preserve an exact name first,
  # because a column may legitimately contain leading or trailing spaces.
  if (by == "set" && !is.null(names(weights)))
    names(weights) <- .role_text_values(names(weights))
  if (!is.numeric(weights) || is.complex(weights) || !is.null(dim(weights)) ||
      !is.null(oldClass(weights)) || !length(weights) ||
      is.null(names(weights)) ||
      any(is.na(names(weights))) || any(!nzchar(trimws(names(weights)))) ||
      anyDuplicated(names(weights)) || any(!is.finite(weights)) ||
      any(weights < 0))
    stop("`weights` must be a non-empty named numeric vector of non-negative finite values",
         call. = FALSE)
  if (!any(weights > 0))
    stop("at least one external weight must be positive", call. = FALSE)

  X <- fit$X
  if (!is.matrix(X)) X <- as.matrix(X)
  cells <- colnames(X)
  if (is.null(cells) || !length(cells) || length(fit$tau_list) != ncol(X))
    stop("the fitted response cells and thresholds are incomplete",
         call. = FALSE)
  if (inherits(fit, c("rasch_mfrm", "rasch_efrm"))) {
    vm <- fit$virtual_map
    if (is.null(vm) || !all(c("vkey", "item") %in% names(vm)))
      stop("the fitted response-cell map is unavailable", call. = FALSE)
    ii <- match(cells, as.character(vm$vkey))
    if (anyNA(ii))
      stop("the fitted response-cell map does not cover every response cell",
           call. = FALSE)
    source_item <- as.character(vm$item[ii])
  } else {
    source_item <- cells
  }
  items <- unique(source_item)
  match_items <- function(x) {
    x <- as.character(x)
    vapply(x, function(nm) {
      if (nm %in% items) return(nm)
      hit <- items[trimws(items) == trimws(nm)]
      if (length(hit) == 1L) hit else nm
    }, character(1), USE.NAMES = FALSE)
  }

  set_of <- NULL
  if (by == "item") {
    if (!is.null(sets))
      stop("`sets` is used only when `by = \"set\"`", call. = FALSE)
    names(weights) <- match_items(names(weights))
    if (anyDuplicated(names(weights)))
      stop("item weights must name every fitted item exactly once",
           call. = FALSE)
    missing_w <- setdiff(items, names(weights))
    extra_w <- setdiff(names(weights), items)
    if (length(missing_w) || length(extra_w))
      stop("item weights must name every fitted item exactly",
           if (length(missing_w)) paste0("; missing: ", paste(missing_w, collapse = ", ")) else "",
           if (length(extra_w)) paste0("; unknown: ", paste(extra_w, collapse = ", ")) else "",
           call. = FALSE)
    q <- unname(weights[source_item])
  } else {
    if (is.null(sets) && inherits(fit, "rasch_efrm")) sets <- fit$set_of
    if (is.null(sets))
      stop("`sets` must map every fitted item when `by = \"set\"`",
           call. = FALSE)
    if (is.list(sets)) {
      if (!is.null(names(sets))) names(sets) <- .role_text_values(names(sets))
      if (!length(sets) || is.null(names(sets)) ||
          any(is.na(names(sets))) || any(!nzchar(trimws(names(sets)))) ||
          anyDuplicated(names(sets)) ||
          !all(vapply(sets, function(z)
            (is.character(z) || is.factor(z)) && is.null(dim(z)) &&
              length(z) && !anyNA(z) &&
              all(nzchar(trimws(as.character(z)))), logical(1))))
        stop("a set list needs unique non-empty names and item-name elements",
             call. = FALSE)
      sets <- lapply(sets, match_items)
      set_items <- unlist(sets, use.names = FALSE)
      set_names <- rep(names(sets), lengths(sets))
      if (anyDuplicated(set_items))
        stop("an item cannot belong to more than one externally weighted set",
             call. = FALSE)
      set_of <- stats::setNames(set_names, set_items)
    } else {
      if (!(is.character(sets) || is.factor(sets)) || !is.null(dim(sets)) ||
          !length(sets) || is.null(names(sets)) || any(is.na(names(sets))) ||
          any(!nzchar(trimws(names(sets)))) || anyDuplicated(names(sets)) ||
          anyNA(sets) || any(!nzchar(trimws(as.character(sets)))))
        stop("`sets` must be a named item-to-set vector or a named list",
             call. = FALSE)
      if (is.character(sets) || is.factor(sets)) {
        set_names <- match_items(names(sets))
        if (anyDuplicated(set_names))
          stop("the set map must name every fitted item exactly once",
               call. = FALSE)
        sets <- .role_text_values(sets)
        names(sets) <- set_names
      }
      set_of <- stats::setNames(as.character(sets), names(sets))
    }
    missing_set <- setdiff(items, names(set_of))
    extra_set <- setdiff(names(set_of), items)
    if (length(missing_set) || length(extra_set))
      stop("the set map must name every fitted item exactly",
           if (length(missing_set)) paste0("; missing: ", paste(missing_set, collapse = ", ")) else "",
           if (length(extra_set)) paste0("; unknown: ", paste(extra_set, collapse = ", ")) else "",
           call. = FALSE)
    set_names <- unique(unname(set_of[items]))
    missing_w <- setdiff(set_names, names(weights))
    extra_w <- setdiff(names(weights), set_names)
    if (length(missing_w) || length(extra_w))
      stop("set weights must name every fitted set exactly",
           if (length(missing_w)) paste0("; missing: ", paste(missing_w, collapse = ", ")) else "",
           if (length(extra_w)) paste0("; unknown: ", paste(extra_w, collapse = ", ")) else "",
           call. = FALSE)
    q <- unname(weights[set_of[source_item]])
  }
  # Relative weights are scale free. Rescale by their maximum before setting
  # mean one so a perfectly valid vector near the smallest representable
  # double cannot have mean zero by underflow (and hence become Inf/NaN).
  positive_q <- q > 0
  q <- q / max(q)
  if (any(positive_q & q == 0))
    stop("the relative range of positive `weights` is too large to represent; ",
         "use a less extreme set of relative weights", call. = FALSE)
  q <- q / mean(q)

  disc <- fit$disc
  if (is.null(disc)) disc <- rep(1, ncol(X))
  if (length(disc) == 1L) disc <- rep(disc, ncol(X))
  if (!is.numeric(disc) || length(disc) != ncol(X) ||
      any(!is.finite(disc)) || any(disc <= 0))
    stop("the fitted model units are unavailable", call. = FALSE)

  obs <- !is.na(X) & rep(q > 0, each = nrow(X))
  m <- vapply(fit$tau_list, length, integer(1))
  theta <- se <- rep(NA_real_, nrow(X))
  weighted_score <- rowSums(sweep(X, 2L, q * disc, `*`), na.rm = TRUE)
  weighted_score[rowSums(obs) == 0L] <- NA_real_
  max_weighted_score <- as.numeric(obs %*% (q * disc * m))
  raw <- rowSums(replace(X, !obs, NA_real_), na.rm = TRUE)
  raw[rowSums(obs) == 0L] <- NA_real_
  max_raw <- as.numeric(obs %*% m)
  # Extreme status is a response-pattern property. A numerical tolerance on
  # the weighted total can call a non-extreme response extreme when one item
  # has a very small positive weight or model unit.
  away_from_zero <- obs & X != 0
  away_from_max <- obs & sweep(X, 2L, m, `!=`)
  extreme <- rowSums(obs) > 0L &
    (rowSums(away_from_zero, na.rm = TRUE) == 0L |
       rowSums(away_from_max, na.rm = TRUE) == 0L)
  pat <- apply(obs, 1L, function(z) paste(which(z), collapse = ","))

  for (key in unique(pat)) {
    cols <- as.integer(strsplit(key, ",", fixed = TRUE)[[1L]])
    if (!length(cols)) next
    who_pat <- which(pat == key)
    # Only relative weights among answered items affect the estimate. A
    # pattern may omit every heavily weighted item, leaving tiny weights
    # whose squared moments would underflow on the global scale.
    qp <- q[cols] / max(q[cols])
    # Dividing every unit by the same constant divides the estimating
    # equation by that constant. Keep probabilities on the fitted scale,
    # but form its powers in relative units to avoid overflow/underflow.
    unit_scale <- max(disc[cols])
    rp <- disc[cols] / unit_scale
    pattern_score <- as.numeric(X[who_pat, cols, drop = FALSE] %*%
                                 (qp * rp))
    curve <- .person_wle_curve(fit$tau_list[cols], disc[cols], qp)
    # Equal weighted totals have the same estimating equation within a
    # missingness pattern. Do not round them for caching: distinct totals can
    # be arbitrarily close when external weights are highly unequal.
    for (Wu in unique(pattern_score)) {
      who <- who_pat[pattern_score == Wu]
      root <- .person_wle_maximum(curve, Wu)
      theta[who] <- root
      if (is.finite(root)) {
        V <- vapply(cols, function(j)
          item_moments(root, fit$tau_list[[j]], disc = disc[j])$V, 0)
        H <- sum(qp * rp^2 * V)
        J <- sum(qp^2 * rp^2 * V)
        se[who] <- (sqrt(J) / H) / unit_scale
      }
    }
  }

  identity_cols <- intersect(unique(c("id", names(fit$factors))),
                             names(fit$person))
  out <- cbind(fit$person[, identity_cols, drop = FALSE],
               data.frame(n_items = rowSums(obs), raw = raw,
                          max_raw = max_raw,
                          weighted_score = weighted_score,
                          max_weighted_score = max_weighted_score,
                          theta = theta, se = se, extreme = extreme,
                          check.names = FALSE))
  weighting <- data.frame(response_cell = cells, item = source_item,
                          set = if (is.null(set_of)) NA_character_
                            else unname(set_of[source_item]),
                          supplied_weight = if (by == "item")
                            unname(weights[source_item])
                          else unname(weights[set_of[source_item]]),
                          normalised_weight = q,
                          stringsAsFactors = FALSE)
  attr(out, "weighting") <- weighting
  attr(out, "by") <- by
  attr(out, "algorithm") <- "pattern-max-wle-3"
  out
}

# Validate the signed app result before it is restored or exported. This is
# deliberately internal: weighted_person_estimates() remains the public way to
# calculate the secondary measures, while the app carries enough information
# to prove that a displayed or downloaded table belongs to the active fit.
.authenticate_weighted_person_result <- function(x, fit) {
  if (is.null(x)) return(invisible(NULL))
  scalar_text <- function(z)
    is.character(z) && length(z) == 1L && !is.na(z) && nzchar(z)
  if (inherits(fit, "rasch_btl") || !is.list(x) ||
      !is.data.frame(x$table) || !scalar_text(x$by) ||
      !x$by %in% c("item", "set") || is.null(x$weights) ||
      !is.list(x$fit_signature) || !scalar_text(x$result_signature))
    stop("the externally weighted person result has an invalid structure",
         call. = FALSE)
  unsigned <- x
  unsigned$result_signature <- NULL
  if (!.fit_boot_hash_matches(x$result_signature, unsigned))
    stop("the externally weighted person result has changed since it was calculated",
         call. = FALSE)
  if (!.fit_boot_signature_matches(x$fit_signature, fit))
    stop("the externally weighted person result was calculated from a different fitted model",
         call. = FALSE)
  invisible(x)
}

.validate_weighted_person_result <- function(x, fit) {
  if (is.null(x)) return(invisible(NULL))
  .authenticate_weighted_person_result(x, fit)
  expected <- weighted_person_estimates(
    fit, x$weights, by = x$by, sets = x$sets)
  if (!identical(expected, x$table))
    stop("the externally weighted person table does not reproduce from the fitted model",
         call. = FALSE)
  invisible(x)
}

.report_person_weight_result <- function(fit) {
  x <- attr(fit, "report_person_weights", exact = TRUE)
  .validate_weighted_person_result(x, fit)
  x
}

# Validate a public weighted_person_estimates() table without relying on the
# signed state used by the application. The resolved weighting attribute is
# sufficient to reconstruct the public call; exact reproduction then checks
# both the table and its metadata against the active fit.
.validate_weighted_person_table <- function(x, fit) {
  if (!is.data.frame(x))
    stop("`person_weights` must be a table returned by weighted_person_estimates()",
         call. = FALSE)
  by <- attr(x, "by", exact = TRUE)
  weighting <- attr(x, "weighting", exact = TRUE)
  required <- c("response_cell", "item", "set", "supplied_weight",
                "normalised_weight")
  if (inherits(fit, "rasch_btl") || !is.character(by) || length(by) != 1L ||
      is.na(by) || !by %in% c("item", "set") ||
      !is.data.frame(weighting) || !nrow(weighting) ||
      !all(required %in% names(weighting)) ||
      !is.character(weighting$item) || anyNA(weighting$item) ||
      any(!nzchar(weighting$item)) ||
      !is.numeric(weighting$supplied_weight) ||
      any(!is.finite(weighting$supplied_weight)) ||
      any(weighting$supplied_weight < 0) ||
      !is.numeric(weighting$normalised_weight) ||
      any(!is.finite(weighting$normalised_weight)) ||
      any(weighting$normalised_weight < 0))
    stop("`person_weights` has an invalid weighting structure", call. = FALSE)

  item_rows <- split(seq_len(nrow(weighting)), weighting$item)
  supplied <- lapply(item_rows, function(i)
    unique(weighting$supplied_weight[i]))
  if (any(lengths(supplied) != 1L))
    stop("`person_weights` has inconsistent weights for an item",
         call. = FALSE)

  if (by == "item") {
    weights <- vapply(supplied, `[[`, numeric(1), 1L)
    sets <- NULL
  } else {
    if (!is.character(weighting$set) || anyNA(weighting$set) ||
        any(!nzchar(weighting$set)))
      stop("`person_weights` has an invalid item-set map", call. = FALSE)
    item_sets <- lapply(item_rows, function(i) unique(weighting$set[i]))
    if (any(lengths(item_sets) != 1L))
      stop("`person_weights` maps an item to more than one set",
           call. = FALSE)
    sets <- vapply(item_sets, `[[`, character(1), 1L)
    set_rows <- split(seq_len(nrow(weighting)), weighting$set)
    set_weights <- lapply(set_rows, function(i)
      unique(weighting$supplied_weight[i]))
    if (any(lengths(set_weights) != 1L))
      stop("`person_weights` has inconsistent weights for an item set",
           call. = FALSE)
    weights <- vapply(set_weights, `[[`, numeric(1), 1L)
  }

  expected <- tryCatch(
    weighted_person_estimates(fit, weights, by = by, sets = sets),
    error = function(e) NULL)
  if (is.null(expected) || !identical(expected, x))
    stop("`person_weights` does not reproduce from the fitted model",
         call. = FALSE)
  invisible(x)
}

.resolve_report_person_weights <- function(fit, person_weights = NULL) {
  if (!is.null(person_weights)) {
    .validate_weighted_person_table(person_weights, fit)
    return(list(table = person_weights))
  }
  .report_person_weight_result(fit)
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
#' @param extremes Treatment of the extreme scores. \code{"model"} keeps the
#'   estimator's own values; these are \code{NA} for MLE.
#'   \code{"extrapolated"} applies the geometric extrapolation.
#' @return A data frame with \code{score}, \code{theta}, \code{se},
#'   \code{freq}, \code{cum_pct} (omitted when no complete responders
#'   exist), and \code{extrapolated}; \code{NULL} when the fitted items do not
#'   share one discrimination or an item is represented by several MFRM or
#'   EFRM response cells.
#' @references
#' Andrich, D. and Marais, I. (2019). A Course in Rasch Measurement Theory:
#' Measuring in the Educational, Social and Health Sciences. Springer.
#'
#' Warm, T. A. (1989). Weighted likelihood estimation of ability in item
#' response theory. Psychometrika, 54(3), 427--450.
#' @examples
#' set.seed(1)
#' d <- seq(-1.5, 1.5, length.out = 6)
#' X <- matrix(rbinom(300 * 6, 1, plogis(outer(rnorm(300), d, "-"))), 300, 6)
#' colnames(X) <- paste0("I", 1:6)
#' score_table(rasch(X), method = "mle", extremes = "extrapolated")
#' @export
score_table <- function(fit, method = c("wle", "mle"),
                        extremes = c("model", "extrapolated")) {
  if (!inherits(fit, "rasch") || inherits(fit, "rasch_btl"))
    stop("`fit` must be a response-data Rasch fit", call. = FALSE)
  if (!isTRUE(fit$est$converged))
    stop("the fitted calibration did not converge; the score table is unavailable",
         call. = FALSE)
  if (!.efrm_link_converged(fit))
    stop("the fitted set-unit link did not converge; the score table is unavailable",
         call. = FALSE)
  method <- match.arg(method); extremes <- match.arg(extremes)
  if (is.null(fit$score_table)) return(NULL)
  tab <- fit$score_table[, c("score", "theta", "se")]
  M <- max(tab$score)
  disc <- if (is.null(fit$disc)) 1 else fit$disc[1]
  se_of <- function(th) .common_person_se(th, fit$tau_list, disc)
  if (method == "mle") {
    interval <- .person_root_interval(fit$tau_list, disc)
    for (r in seq_len(M - 1)) {
      # the ML equation is sum(disc*(x_i - E_i)) = 0: the common
      # discrimination cancels, so the root solves r = sum(E), not
      # r = sum(disc*E)
      tab$theta[tab$score == r] <- uniroot(function(th)
        r - sum(vapply(fit$tau_list, function(tt)
          item_moments(th, tt, disc = disc)$E, 0)),
        interval, tol = 1e-9 / disc)$root
    }
    tab$theta[c(1, M + 1)] <- NA_real_
    tab$se <- vapply(tab$theta, function(th)
      if (is.na(th)) NA_real_ else se_of(th), 0)
  }
  tab$extrapolated <- FALSE
  if (extremes == "extrapolated") {
    if (M < 4) stop("extrapolation needs at least three interior scores")
    th <- tab$theta
    ext <- .geometric_score_extremes(th)
    if (is.null(ext))
      stop("extreme-score extrapolation is unavailable because the interior ",
           "score locations do not define stable increasing spacings",
           call. = FALSE)
    lo <- ext[["lower"]]
    hi <- ext[["upper"]]
    ext_se <- c(se_of(lo), se_of(hi))
    if (any(!is.finite(ext_se)) || any(ext_se <= 0))
      stop("extreme-score extrapolation is unavailable because test ",
           "information is not positive and finite at the continued locations",
           call. = FALSE)
    tab$theta[c(1, M + 1)] <- c(lo, hi)
    tab$se[c(1, M + 1)] <- ext_se
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

.geometric_score_extremes <- function(theta) {
  M <- length(theta) - 1L
  if (M < 4L || any(!is.finite(theta[2:M]))) return(NULL)
  spacing <- c(theta[3] - theta[2], theta[4] - theta[3],
               theta[M] - theta[M - 1L],
               theta[M - 1L] - theta[M - 2L])
  if (any(!is.finite(spacing))) return(NULL)
  # The tolerance is relative to the score spacings, not the arbitrary
  # origin of the logit scale. This preserves origin invariance.
  tol <- 100 * .Machine$double.eps * max(abs(spacing))
  if (!is.finite(tol) || any(spacing <= tol)) return(NULL)
  lower <- theta[2] - spacing[1] * (spacing[1] / spacing[2])
  upper <- theta[M] + spacing[3] * (spacing[3] / spacing[4])
  if (!is.finite(lower) || !is.finite(upper) ||
      lower >= theta[2] || upper <= theta[M]) return(NULL)
  c(lower = lower, upper = upper)
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
#' @param fit A fitted object from \code{\link{rasch}} with one common
#'   raw-score conversion and a common discrimination. Expanded many-facet
#'   and frame response-cell designs do not have such a conversion and are
#'   refused.
#' @return The fit's person table with two added columns,
#'   \code{theta_extrapolated} and \code{se_extrapolated}: equal to
#'   \code{theta} and \code{se} for non-extreme persons, extrapolated for
#'   extreme persons. Patterns with fewer than three interior scores cannot
#'   be extrapolated, or whose interior score locations do not have stable
#'   increasing spacings, keep their Warm values.
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
  if (!inherits(fit, "rasch") || inherits(fit, "rasch_btl"))
    stop("`fit` must be a fitted Rasch-family object", call. = FALSE)
  if (!isTRUE(fit$est$converged))
    stop("the fitted calibration did not converge; extrapolated person estimates are unavailable",
         call. = FALSE)
  if (!.efrm_link_converged(fit))
    stop("the fitted set-unit link did not converge; extrapolated person estimates are unavailable",
         call. = FALSE)
  if (is.null(fit$score_table))
    stop("a common raw-score conversion is not defined across expanded ",
         "frame or facet response cells", call. = FALSE)
  disc <- if (is.null(fit$disc)) 1 else unique(fit$disc)
  if (length(disc) != 1L || !is.finite(disc) || disc <= 0)
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
    pe <- person_wle(tl, disc = disc)
    th <- unname(pe$theta); M <- length(th) - 1L
    if (M < 4) next                        # too few interior scores
    endpoints <- .geometric_score_extremes(th)
    if (is.null(endpoints)) next
    lo <- endpoints[["lower"]]
    hi <- endpoints[["upper"]]
    endpoint_se <- c(lower = .common_person_se(lo, tl, disc),
                     upper = .common_person_se(hi, tl, disc))
    for (r in rows) {
      raw <- sum(X[r, obs])
      if (raw == 0L && is.finite(endpoint_se[["lower"]]) &&
          endpoint_se[["lower"]] > 0) {
        p$theta_extrapolated[r] <- lo
        p$se_extrapolated[r] <- endpoint_se[["lower"]]
      } else if (raw != 0L && is.finite(endpoint_se[["upper"]]) &&
                 endpoint_se[["upper"]] > 0) {
        p$theta_extrapolated[r] <- hi
        p$se_extrapolated[r] <- endpoint_se[["upper"]]
      }
    }
  }
  p
}
