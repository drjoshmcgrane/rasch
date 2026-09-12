# rasch :: data simulation
#
# Generate data from the model family with dial-in departures from it, so a
# known pathology can be planted and the matching diagnostic watched as it
# fires. Every simulator returns data ready for its fit function, with the
# true generating parameters attached as attr(x, "truth") and a class that
# prints a summary of what was planted.
#
# Shared truth schema (all simulators populate what applies):
#   list(layout, n_*, theta/locations, difficulty/thresholds, discrimination,
#        guessing, groups, departure_types = <character, model departures>,
#        planted = <character, human-readable generating settings>)

# null-coalescing helper (package-internal; base R gained %||% only in 4.4)
`%||%` <- function(a, b) if (is.null(a)) b else a

# A seeded simulator call should be reproducible without commandeering the
# caller's random number stream: capture the stream, seed, and restore on
# exit, so simulate_*(seed = s) twice gives the same data while code after
# the call draws exactly what it would have drawn anyway. Restoring an
# absent stream means removing the one set.seed() created. Box-Muller's
# cached second normal is not stored in .Random.seed and cannot be restored.
.sim_seed_check <- function() {
  if (identical(RNGkind()[2L], "Box-Muller"))
    stop("this calculation cannot preserve the Box-Muller random-number ",
         "stream; use RNGkind(normal.kind = \"Inversion\") before setting ",
         "your analysis seed (see ?rasch_rng)", call. = FALSE)
  invisible(NULL)
}

.sim_seed_capture <- function() {
  .sim_seed_check()
  if (exists(".Random.seed", envir = globalenv(), inherits = FALSE))
    get(".Random.seed", envir = globalenv(), inherits = FALSE)
  else NULL
}

.sim_seed_restore <- function(old) {
  if (is.null(old)) {
    if (exists(".Random.seed", envir = globalenv(), inherits = FALSE))
      rm(".Random.seed", envir = globalenv())
  } else {
    assign(".Random.seed", old, envir = globalenv())
  }
  invisible(NULL)
}

#' Random-number generation
#'
#' Supplying a \code{seed} makes a simulation or bootstrap reproducible and
#' restores the caller's random-number stream on exit. Bootstrap methods
#' that assign seeds to individual replicates also restore those local streams.
#'
#' These operations do not support R's Box--Muller normal generator: its
#' cached normal value is not part of \code{.Random.seed}, so restoring that
#' vector would change subsequent draws. They refuse before changing the
#' stream. Direct \code{simulate_*} calls with \code{seed = NULL} can still
#' use Box--Muller. \code{\link{sim_replicate}} assigns replicate seeds even
#' when its own \code{seed} is \code{NULL}.
#'
#' The default Inversion generator is supported. To select it explicitly,
#' use \code{RNGkind(normal.kind = "Inversion")} before setting the seed
#' for the analysis. Changing the generator starts a different normal stream;
#' it does not recover a previous Box--Muller stream.
#' @name rasch_rng
#' @seealso \code{\link[base]{Random}}, \code{\link{simulate_rasch}},
#'   \code{\link{fit_bootstrap}}, \code{\link{dif_bootstrap}}.
NULL

.sim_count <- function(x, name, min = 1L) {
  if (length(x) != 1L || !is.numeric(x) || is.complex(x) ||
      !is.null(dim(x)) || !is.null(oldClass(x)) || !is.finite(x) ||
      x != floor(x) || x < min ||
      x > .Machine$integer.max)
    stop(name, " must be one whole number >= ", min)
  as.integer(x)
}

.sim_scalar <- function(x, name, lower = -Inf, upper = Inf,
                        lower_open = FALSE, upper_open = FALSE) {
  ok <- length(x) == 1L && is.numeric(x) && !is.complex(x) &&
    is.null(dim(x)) && is.null(oldClass(x)) && is.finite(x) &&
    if (lower_open) x > lower else x >= lower
  ok <- ok && if (upper_open) x < upper else x <= upper
  if (!ok) {
    left <- if (lower_open) "(" else "["
    right <- if (upper_open) ")" else "]"
    stop(name, " must be one finite value in ", left, lower, ", ", upper, right)
  }
  as.numeric(x)
}

.sim_seed <- function(seed) {
  if (is.null(seed)) return(NULL)
  .sim_count(seed, "seed", 0L)
}

.sim_vector <- function(x, name, lengths, lower = -Inf, upper = Inf,
                        lower_open = FALSE, upper_open = FALSE) {
  ok <- is.numeric(x) && !is.complex(x) && is.null(dim(x)) &&
    is.null(oldClass(x)) && length(x) %in% lengths && all(is.finite(x))
  if (ok) {
    ok <- all(if (lower_open) x > lower else x >= lower) &&
      all(if (upper_open) x < upper else x <= upper)
  }
  if (!ok)
    stop(name, " must contain plain finite numeric values with length ",
         paste(lengths, collapse = " or "))
  as.numeric(x)
}

.sim_structure <- function(x, name, allowed, required = character()) {
  if (is.null(x)) return(NULL)
  if (!is.list(x) || is.data.frame(x))
    stop(name, " must be NULL or a named list")
  if (!length(x)) {
    if (length(required))
      stop(name, " must contain ", paste(required, collapse = ", "))
    return(x)
  }
  nm <- names(x)
  if (is.null(nm) || anyNA(nm) || any(!nzchar(trimws(nm))) || anyDuplicated(nm))
    stop(name, " must be a named list with unique, non-missing names")
  unknown <- setdiff(nm, allowed)
  if (length(unknown))
    stop(name, " has unknown component(s): ", paste(unknown, collapse = ", "))
  missing <- setdiff(required, nm)
  if (length(missing))
    stop(name, " must contain ", paste(missing, collapse = ", "))
  x
}

# Construct a fixed departure that the explanatory design cannot absorb. The
# residual of a unit vector is orthogonal to every fitted design column; scale
# it so the best-supported focal row carries the requested magnitude. A
# saturated design has no such departure and must be refused rather than
# advertised as planted misfit.
.sim_explanatory_departure <- function(B, magnitude) {
  B <- as.matrix(B)
  magnitude <- .sim_scalar(magnitude, "explanatory departure", lower = 0,
                           lower_open = TRUE)
  if (!is.numeric(B) || !nrow(B) || !ncol(B) || any(!is.finite(B)))
    stop("the explanatory simulation design must be a finite numeric matrix")
  q <- qr(B, tol = 1e-10)
  if (q$rank >= nrow(B))
    stop("the explanatory design is saturated; increase the number of items ",
         "or objects before planting a fixed departure")
  Q <- qr.Q(q)[, seq_len(q$rank), drop = FALSE]
  residual_leverage <- pmax(1 - rowSums(Q^2), 0)
  focal <- which.max(residual_leverage)
  if (!length(focal) || residual_leverage[focal] <= 1e-10)
    stop("the explanatory design has no stable row on which to plant a fixed departure")
  direction <- -drop(Q %*% Q[focal, ])
  direction[focal] <- direction[focal] + 1
  direction <- direction * magnitude / direction[focal]
  list(values = direction, index = focal)
}

# Construct a finite object vector with the requested sample correlation to x.
# The random component is centred and projected off x before being rescaled,
# rather than merely combined with x and expected to be orthogonal in a small
# realised object set.
.sim_correlated_fixed <- function(x, rho, target_sd,
                                  draw = function(n) stats::rnorm(n)) {
  if (abs(rho) == 1) return(rho * x)
  xx <- sum(x^2)
  for (attempt in seq_len(20L)) {
    z <- draw(length(x))
    z <- z - mean(z)
    z <- z - x * sum(z * x) / xx
    zsd <- stats::sd(z)
    if (is.finite(zsd) && zsd > sqrt(.Machine$double.eps)) {
      z <- z / zsd * target_sd
      return(rho * x + sqrt(1 - rho^2) * z)
    }
  }
  stop("could not construct the requested finite object correlation")
}

# Allocate a finite set of observations as evenly as possible across the
# requested IDs. Unlike sampling with replacement, this guarantees that every
# declared judge contributes data when the design contains enough rows.
.sim_balanced_ids <- function(ids, n, what = "levels") {
  if (n < length(ids))
    stop("the generated design has ", n, " observations, fewer than the ",
         length(ids), " requested ", what)
  out <- rep(ids, n %/% length(ids))
  rem <- n %% length(ids)
  if (rem) out <- c(out, sample(ids, rem, replace = FALSE))
  sample(out, length(out), replace = FALSE)
}

.sim_planted_count <- function(prop, n) {
  if (prop <= 0) return(0L)
  min(n, max(1L, as.integer(round(prop * n))))
}

.sim_check_sd <- function(x, target, name) {
  if (target == 0) return(invisible(NULL))
  observed <- stats::sd(x)
  # all.equal() switches to an absolute tolerance near zero, so it regards a
  # collapsed SD of zero as equal to a small positive request. Test the
  # relative scale directly, including for subnormal finite values.
  relative_error <- if (is.finite(observed) && observed > 0)
    abs(observed / target - 1) else Inf
  if (!is.finite(relative_error) || relative_error > 1e-10)
    stop("the requested ", name, " cannot be represented at this numerical scale",
         call. = FALSE)
  invisible(NULL)
}

# person locations from one of a few distributions
.sim_theta <- function(n, mean, sd, dist = "normal") {
  n <- .sim_count(n, "n")
  mean <- .sim_scalar(mean, "mean")
  sd <- .sim_scalar(sd, "sd", lower = 0)
  if (sd == 0) return(rep(mean, n))
  z <- switch(dist,
    normal  = stats::rnorm(n),
    uniform = stats::runif(n, -sqrt(3), sqrt(3)),
    skew    = { u <- stats::rgamma(n, 2, 1); (u - 2) / sqrt(2) },
    bimodal = { s <- sample(c(-1, 1), n, TRUE); s * 1.1 + stats::rnorm(n, 0, 0.5) },
    stats::rnorm(n))
  if (n == 1L) return(mean)
  zsd <- stats::sd(z)
  if (!is.finite(zsd) || zsd <= sqrt(.Machine$double.eps))
    z <- seq_len(n)
  out <- mean + sd * as.numeric(scale(z))
  if (any(!is.finite(out)))
    stop("the requested person mean and standard deviation are outside the numerically representable range")
  .sim_check_sd(out, sd, "sample standard deviation")
  out
}

# draw one item's responses for a vector of person locations. tau are the
# item's thresholds (length m; dichotomous m = 1). disc scales the whole
# exponent (a departure when != 1); guess is a lower asymptote (dichotomous).
.sim_item <- function(theta, tau, disc = 1, guess = 0) {
  if (!length(theta) || any(!is.finite(theta)) ||
      !length(tau) || any(!is.finite(tau)))
    stop("theta and tau must contain finite values")
  disc <- .sim_scalar(disc, "disc", lower = 0, lower_open = TRUE)
  guess <- .sim_scalar(guess, "guess", lower = 0, upper = 1,
                       upper_open = TRUE)
  m <- length(tau); xs <- 0:m
  cum <- c(0, cumsum(tau))
  eta <- disc * (outer(theta, xs) - matrix(cum, length(theta), m + 1L, byrow = TRUE))
  if (any(!is.finite(cum)) || any(!is.finite(eta)))
    stop("the generating locations, thresholds and discrimination are ",
         "outside the numerically representable range", call. = FALSE)
  eta <- eta - apply(eta, 1, max)
  P <- exp(eta); P <- P / rowSums(P)
  if (guess > 0 && m == 1L) {
    P[, 2] <- guess + (1 - guess) * P[, 2]; P[, 1] <- 1 - P[, 2]
  }
  cs <- t(apply(P, 1, cumsum))
  as.integer(rowSums(stats::runif(length(theta)) > cs))     # category 0..m
}

# thresholds around an item location: evenly spread (rating-scale pattern),
# or with a supplied relative pattern (partial-credit: pattern varies by
# item); optionally deliberately disordered (one threshold dropped below its
# predecessor, re-centred so the item's mean location is preserved)
.sim_thresholds <- function(delta, m, spread, disordered = FALSE,
                            pattern = NULL) {
  if (m == 1L) return(delta)
  step <- if (is.null(pattern)) seq(-1, 1, length.out = m) * spread
          else pattern
  tau <- delta + step - mean(step)
  if (disordered && m >= 2L) {
    i <- max(2L, ceiling(m / 2))           # always has a predecessor to undercut
    # Set the adjacent reversal directly. Subtracting a fixed amount from the
    # original threshold did not guarantee disorder for PCM patterns whose
    # randomly generated gap happened to be wider than that amount.
    tau[i] <- tau[i - 1L] - 0.2 * spread
    tau <- tau - mean(tau) + delta         # keep the item's location honest
  }
  tau
}

#' Simulate person-by-item Rasch data
#'
#' Generates dichotomous, partial credit, or rating scale data. Optional
#' arguments introduce item misfit, guessing, multidimensionality, local
#' dependence, DIF, response styles, or missingness. Generating values are
#' stored in \code{attr(x, "truth")}. A positive planted proportion selects
#' at least one person or response cell when the requested departures can
#' coexist.
#'
#' @param n_persons,n_items Sample size and test length.
#' @param model \code{"dichotomous"}, \code{"PCM"}, or \code{"RSM"}. Under
#'   \code{"RSM"} every item shares one category-threshold pattern (items
#'   differ by location only); under \code{"PCM"} each item's threshold
#'   spacings and span are drawn afresh, as the partial credit model allows.
#' @param n_categories Response categories for polytomous models (>= 3).
#' @param theta_mean,theta_sd Realised sample mean and standard deviation of
#'   the person distribution.
#' @param theta_dist Shape of the person distribution: \code{"normal"},
#'   \code{"uniform"}, \code{"skew"}, or \code{"bimodal"}.
#' @param difficulty Either the two endpoints of an evenly spaced location
#'   range, or one location per item.
#' @param threshold_spread Half-range of the category thresholds about each
#'   item location (polytomous).
#' @param discrimination The item slope, supplied as one value or one per item.
#'   One common value changes the logit unit without departing from a Rasch
#'   model. Differences between items are a deliberate slope departure: values
#'   above the common pattern produce steeper responses and values below it
#'   produce flatter responses, and require \code{theta_sd > 0}.
#' @param guessing Scalar or length-\code{n_items} lower asymptote
#'   (dichotomous): low-location persons answer correctly by chance. A positive
#'   asymptote requires \code{theta_sd > 0}.
#' @param second_dim \code{NULL}, or \code{list(items=, rho=)}: the named items
#'   load on a second trait whose realised sample correlation with the first
#'   is \code{rho}. It lies in [-1, 1); a correlation of 1 would reproduce the
#'   primary trait rather than plant another dimension. At least three persons
#'   are needed unless \code{rho} is -1. Each item is named once, and at least
#'   one item must remain on
#'   the primary trait; moving every item to the second trait is still a
#'   one-dimensional scale.
#' @param dependence \code{NULL}, or \code{list(pairs=, strength=)}: each pair's
#'   second item responds partly to the first. This departure feeds the
#'   residual-dependence diagnostics. Each directed pair is listed once.
#' @param dif \code{NULL}, or \code{list(items=, uniform=, nonuniform=)}: the
#'   named items function differently for the last person group: a location
#'   shift (\code{uniform}) and/or a slope change (\code{nonuniform}). Needs
#'   \code{n_groups >= 2}; each item is named once and at least one item must
#'   remain invariant. A common shift or slope change on every item changes
#'   the group origin or unit rather than defining item DIF. A non-uniform
#'   effect also requires \code{theta_sd > 0}.
#' @param careless Proportion of persons who answer at random. At least one
#'   person in each generated group must retain model-based responses.
#'   Careless and response-style assignments are disjoint; their requested
#'   counts must fit.
#' @param response_style \code{NULL}, or \code{list(type=, prop=, strength=)}
#'   with \code{type} \code{"extreme"} or \code{"middle"}: a proportion
#'   \code{prop} of persons favour the end (or middle) categories regardless
#'   of the trait, with distortion \code{strength} (default 1.6) on the
#'   log-probability scale (polytomous). At least one person must remain
#'   without the style; applying one category weighting to everyone merely
#'   changes the fitted thresholds.
#' @param speeded Proportion not reached at the last item: a growing tail of
#'   missing responses over the final items. These cells are kept distinct
#'   from any completely-at-random missing cells. Persons are selected
#'   independently of their ability and responses; this does not simulate
#'   non-ignorable missingness or change the response model.
#' @param disordered \code{NULL} or item names/indices given disordered
#'   thresholds (polytomous; feeds the threshold diagnostics).
#' @param n_groups Number of equal person groups (a \code{group} factor column
#'   is added when > 1, for DIF).
#' @param missing Proportion of responses set missing completely at random,
#'   drawn from cells not already missing through speededness. The requested
#'   count must fit among those cells and leave at least one observed response.
#' @param seed Optional non-negative whole-number RNG seed. See
#'   \code{\link{rasch_rng}} for generator support.
#' @return A data frame of class \code{"rasch_sim"} (item columns
#'   \code{I01}..., an \code{id} column, and a \code{group} column when
#'   grouped), with \code{attr(x, "truth")} holding the generating parameters
#'   and the planted departures.
#' @examples
#' # a clean scale with one over-discriminating item and one DIF item
#' d <- simulate_rasch(400, 12, discrimination = c(3, rep(1, 11)),
#'                     dif = list(items = "I06", uniform = 1), n_groups = 2,
#'                     seed = 1)
#' fit <- rasch(d, id = "id", factors = "group")
#' fit$items[c("item", "infit_ms", "outfit_ms")]   # item 1 misfits
#' dif_anova(fit)$summary                           # item 6 flags
#' @export
simulate_rasch <- function(n_persons = 500, n_items = 20,
                           model = c("dichotomous", "PCM", "RSM"),
                           n_categories = 3, theta_mean = 0, theta_sd = 1,
                           theta_dist = "normal", difficulty = c(-2.5, 2.5),
                           threshold_spread = 1.2, discrimination = 1,
                           guessing = 0, second_dim = NULL, dependence = NULL,
                           dif = NULL, careless = 0, response_style = NULL,
                           speeded = 0, disordered = NULL,
                           n_groups = 1, missing = 0, seed = NULL) {
  seed <- .sim_seed(seed)
  if (!is.null(seed)) {
    .old_stream <- .sim_seed_capture()
    on.exit(.sim_seed_restore(.old_stream), add = TRUE)
    set.seed(seed)
  }
  model <- match.arg(model)
  N <- .sim_count(n_persons, "n_persons", 2L)
  I <- .sim_count(n_items, "n_items", 2L)
  n_groups <- .sim_count(n_groups, "n_groups")
  if (n_groups > N) stop("n_groups cannot exceed n_persons")
  theta_mean <- .sim_scalar(theta_mean, "theta_mean")
  theta_sd <- .sim_scalar(theta_sd, "theta_sd", lower = 0)
  threshold_spread <- .sim_scalar(threshold_spread, "threshold_spread",
                                  lower = 0, lower_open = TRUE)
  careless <- .sim_scalar(careless, "careless", lower = 0, upper = 1)
  speeded <- .sim_scalar(speeded, "speeded", lower = 0, upper = 1)
  missing <- .sim_scalar(missing, "missing", lower = 0, upper = 1,
                         upper_open = TRUE)
  theta_dist <- match.arg(theta_dist, c("normal", "uniform", "skew", "bimodal"))
  if (model != "dichotomous")
    n_categories <- .sim_count(n_categories, "n_categories", 3L)
  m <- if (model == "dichotomous") 1L else as.integer(n_categories) - 1L
  second_dim <- .sim_structure(second_dim, "second_dim",
                               c("items", "rho"), "items")
  dependence <- .sim_structure(dependence, "dependence",
                               c("pairs", "strength"), "pairs")
  dif <- .sim_structure(dif, "dif",
                        c("items", "uniform", "nonuniform"), "items")
  response_style <- .sim_structure(response_style, "response_style",
                                   c("type", "prop", "strength"))
  if (!is.null(response_style) && m == 1L) {
    warning("response_style applies to polytomous items only; ignored for dichotomous data")
    response_style <- NULL
  }
  inm <- sprintf("I%02d", seq_len(I))
  as_idx <- function(x, allow_duplicates = FALSE) {
    if (is.null(x) || !length(x)) return(integer(0))
    if (!is.null(dim(x)))
      stop("item selectors must be plain vectors, not matrices or arrays")
    if (is.character(x)) {
      i <- match(x, inm)
    } else if (is.numeric(x) && !is.complex(x) && is.null(oldClass(x))) {
      if (any(!is.finite(x)) || any(x != round(x)))
        stop("item index(es) must be finite whole numbers")
      i <- as.integer(x)
    } else {
      stop("item selectors must be item names or plain numeric indices")
    }
    if (anyNA(i) || any(i < 1L | i > I))
      stop("unknown item name(s)/index(es): ",
           paste(x[is.na(i) | i < 1L | i > I], collapse = ", "))
    if (!allow_duplicates && anyDuplicated(i))
      stop("item selectors must name each generated item at most once")
    i
  }

  # item locations, thresholds, slopes, guessing (with per-item overrides)
  difficulty <- .sim_vector(difficulty, "difficulty", unique(c(2L, I)))
  delta <- setNames(if (length(difficulty) == I) difficulty
                    else seq(difficulty[1], difficulty[2], length.out = I), inm)
  discrimination <- .sim_vector(discrimination, "discrimination",
                                unique(c(1L, I)), lower = 0,
                                lower_open = TRUE)
  guessing <- .sim_vector(guessing, "guessing", unique(c(1L, I)),
                          lower = 0, upper = 1, upper_open = TRUE)
  disc <- if (length(discrimination) == I) discrimination else
    rep(discrimination, I)
  guess <- if (length(guessing) == I) guessing else rep(guessing, I)
  if (theta_sd == 0 && length(unique(disc)) > 1L)
    stop("heterogeneous discrimination requires theta_sd > 0; with every ",
         "person at one trait location, item-specific probabilities are ",
         "absorbed by the item calibration rather than identifying slopes")
  if (theta_sd == 0 && any(guess > 0))
    stop("guessing requires theta_sd > 0; a lower asymptote cannot be ",
         "distinguished from item location at one person-trait value")
  if (m > 1L && any(guess > 0)) {
    warning("guessing applies to dichotomous items only; ignored for ", model)
    guess[] <- 0
  }
  dis_items <- as_idx(disordered)
  # a dichotomous item has one threshold and cannot be disordered: planting
  # it would record a generating feature the data do not carry
  if (m == 1L && length(dis_items)) {
    warning("disordered thresholds apply to polytomous items only; ",
            "ignored for dichotomous data")
    dis_items <- integer(0)
  }
  # RSM: one common threshold pattern (per-item location only). PCM: the
  # pattern itself varies across items, as the model allows -- each item's
  # spacings are drawn afresh (gated so the RNG stream of the other models
  # is untouched)
  patterns <- vector("list", I)
  if (model == "PCM" && m > 1L) patterns <- lapply(seq_len(I), function(i) {
    # ordered spacings with varied gaps and a varied overall span per item
    p <- cumsum(stats::runif(m, 0.5, 1.5))
    p <- (p - mean(p)) / (max(p) - min(p)) * 2 * threshold_spread *
      stats::runif(1, 0.85, 1.15)
    p
  })
  tau <- lapply(seq_len(I), function(i)
    .sim_thresholds(delta[i], m, threshold_spread, i %in% dis_items,
                    pattern = patterns[[i]]))
  names(tau) <- inm
  relative_thresholds <- lapply(tau, function(z) z - mean(z))
  rsm_pattern_departure <- identical(model, "RSM") &&
    length(relative_thresholds) > 1L &&
    any(!vapply(relative_thresholds[-1L], function(z)
      isTRUE(all.equal(z, relative_thresholds[[1L]],
                       tolerance = sqrt(.Machine$double.eps),
                       check.attributes = FALSE)), logical(1)))

  # person locations (primary) and groups
  theta <- .sim_theta(N, theta_mean, theta_sd, theta_dist)
  group <- if (n_groups > 1L)
    factor(sprintf("g%d", (seq_len(N) - 1L) %% n_groups + 1L),
           levels = sprintf("g%d", seq_len(n_groups))) else NULL

  # a second dimension for the nominated items: a correlated latent trait
  theta2 <- NULL; dim_items <- integer(0)
  if (!is.null(second_dim)) {
    dim_items <- as_idx(second_dim$items)
    if (!length(dim_items))
      stop("second_dim$items must name at least one item")
    if (length(dim_items) == I)
      stop("second_dim$items must leave at least one item on the primary ",
           "trait; assigning every item to another trait remains a ",
           "one-dimensional scale")
    rho <- second_dim$rho %||% 0.5
    if (length(rho) != 1L || !is.numeric(rho) || is.complex(rho) ||
        !is.null(dim(rho)) || !is.null(oldClass(rho)) || !is.finite(rho) ||
        rho < -1 || rho >= 1)
      stop("second_dim$rho must be a single correlation in [-1, 1); at 1 ",
           "the nominated items still use the primary trait")
    if (theta_sd == 0)
      stop("second_dim requires theta_sd > 0 because a latent correlation is otherwise undefined")
    if (N < 3L && abs(rho) < 1)
      stop("second_dim needs at least three persons to realise a finite ",
           "correlation strictly between -1 and 1")
    # Construct the orthogonal component in the realised sample. This makes
    # the requested mean, spread and correlation exact rather than relying on
    # an independent random draw to be orthogonal in a small sample.
    theta2 <- theta_mean + .sim_correlated_fixed(
      theta - theta_mean, rho, theta_sd,
      draw = function(n) .sim_theta(n, 0, theta_sd, theta_dist))
  }

  X <- matrix(NA_integer_, N, I, dimnames = list(NULL, inm))
  dif_items <- as_idx(if (is.null(dif)) NULL else dif$items)
  if (!is.null(dif) && !length(dif_items))
    stop("dif$items must name at least one generated item")
  if (!is.null(dif)) {
    du <- .sim_scalar(dif$uniform %||% 0, "dif$uniform")
    dn <- .sim_scalar(dif$nonuniform %||% 0, "dif$nonuniform")
    if (theta_sd == 0 && dn != 0)
      stop("non-uniform DIF requires theta_sd > 0; a group difference in ",
           "slope cannot be identified at one person-trait value")
    if (length(dif_items) && any(disc[dif_items] + dn <= 0))
      stop("dif$nonuniform makes a planted item discrimination non-positive")
    dif$uniform <- du; dif$nonuniform <- dn
    if (du == 0 && dn == 0) {
      dif <- NULL
      dif_items <- integer(0)
    } else if (length(dif_items) == I) {
      stop("dif$items must leave at least one invariant item; a common DIF ",
           "effect on every item changes the group origin or unit and is not ",
           "identified as item DIF")
    }
  }
  if (length(dif_items) && n_groups < 2L)
    stop("dif needs n_groups >= 2 (the last group carries the DIF)")
  dif_grp <- if (n_groups > 1L) levels(group)[n_groups] else NA

  # every regeneration of an item must honour that item's OWN generating
  # structure -- its trait (second dimension), its group-shifted thresholds
  # and slope (DIF), its guessing -- or a later misfit layer would silently
  # erase an earlier planted one
  item_pars <- function(i, p) {
    dif_here <- i %in% dif_items && !is.na(dif_grp) && group[p] == dif_grp
    list(tau = tau[[i]] + if (dif_here) dif$uniform %||% 0 else 0,
         disc = disc[i] + if (dif_here) dif$nonuniform %||% 0 else 0)
  }
  gen_item <- function(i, shift = rep(0, N)) {
    th <- (if (i %in% dim_items) theta2 else theta) + shift
    if (i %in% dif_items && !is.na(dif_grp)) {
      g2 <- group == dif_grp
      out <- integer(N)
      out[!g2] <- .sim_item(th[!g2], tau[[i]], disc[i], guess[i])
      out[g2]  <- .sim_item(th[g2], tau[[i]] + (dif$uniform %||% 0),
                            disc[i] + (dif$nonuniform %||% 0), guess[i])
      out
    } else .sim_item(th, tau[[i]], disc[i], guess[i])
  }
  # the model expectation of item i for every person, under the same
  # generating structure (trait, DIF shift, guessing) the responses used --
  # including any dependence carry-over already applied to that item, or in
  # a chain the residual of the middle item would carry the first pair's
  # shift as a systematic mean and leak it into the second pair
  exp_item <- function(i, shift = NULL) {
    th <- (if (i %in% dim_items) theta2 else theta) + (shift %||% 0)
    E <- vapply(seq_len(N), function(p) {
      pp <- item_pars(i, p)
      sum((0:m) * .p_item(th[p], pp$tau, pp$disc))
    }, 0)
    if (m == 1L && guess[i] > 0) E <- guess[i] + (1 - guess[i]) * E
    E
  }

  for (i in seq_len(I)) X[, i] <- gen_item(i)

  # response dependence: the second item of each pair partly follows the
  # first (adds d*(x1 - E1) to its exponent, inducing residual correlation);
  # the regeneration keeps i2's own DIF / second-dimension structure
  dep_pairs <- list()
  dep_pair_idx <- list()
  dep_shift <- vector("list", I)
  dep_sources <- integer(0)
  dep_keys <- character(0)
  if (!is.null(dependence)) {
    d_str <- .sim_scalar(dependence$strength %||% 1,
                         "dependence$strength")
    if (d_str == 0) dependence <- NULL
    if (!is.null(dependence) && !length(dependence$pairs))
      stop("dependence$pairs must contain at least one directed item pair")
    for (pp in if (is.null(dependence)) list() else dependence$pairs) {
      ij <- as_idx(pp, allow_duplicates = TRUE)
      if (length(ij) != 2L || ij[1L] == ij[2L])
        stop("each dependence pair must name two different items")
      i1 <- ij[1]; i2 <- ij[2]
      pair_key <- paste(i1, i2, sep = "->")
      if (pair_key %in% dep_keys)
        stop("dependence$pairs repeats the directed pair ", inm[i1], " -> ",
             inm[i2], "; list each directed pair once")
      dep_keys <- c(dep_keys, pair_key)
      # regenerating i2 replaces the responses any earlier pair drew its
      # carry-over from, so a source may not become a later target
      if (i2 %in% dep_sources)
        stop("item ", inm[i2], " is the source of an earlier dependence ",
             "pair and the target of a later one: regenerating it would ",
             "erase the earlier dependence -- order the pairs so that each ",
             "item is a target before it is a source")
      dep_sources <- c(dep_sources, i1)
      # the expectation must match X1's actual generating structure
      # (guessing, DIF, second dimension), or the "residual" x1 - E1 has a
      # systematic mean and leaks an unplanted shift into the second item
      # the source's expectation is taken under the structure its responses
      # were actually drawn from, its own carry-over included
      shift <- d_str * (X[, i1] - exp_item(i1, dep_shift[[i1]])) / m
      # an item named as the second of several pairs carries every one of
      # them: regenerating from the new pair alone would erase the earlier
      # dependence and plant only the last
      dep_shift[[i2]] <- (dep_shift[[i2]] %||% 0) + shift
      X[, i2] <- gen_item(i2, shift = dep_shift[[i2]])
      dep_pairs[[length(dep_pairs) + 1L]] <- inm[ij]
      dep_pair_idx[[length(dep_pair_idx) + 1L]] <- ij
    }
  }

  # response styles (polytomous): a proportion of persons distort toward the
  # end or the middle categories regardless of the trait; the base
  # probabilities keep each item's own structure (trait, DIF) per person
  style_idx <- integer(0)
  if (!is.null(response_style) && m >= 2L) {
    stype <- response_style$type %||% "extreme"
    if (!is.character(stype) || length(stype) != 1L || is.na(stype) ||
        !is.null(dim(stype)) || !is.null(oldClass(stype)))
      stop("response_style$type must be one of 'extreme' or 'middle'")
    stype <- match.arg(stype, c("extreme", "middle"))
    response_style$type <- stype
    sprop <- .sim_scalar(response_style$prop %||% 0.15,
                         "response_style$prop", lower = 0, upper = 1)
    ss <- .sim_scalar(response_style$strength %||% 1.6,
                      "response_style$strength", lower = 0)
    # a style of zero strength distorts nothing and no one: drawing the
    # persons anyway would record an effect the data do not carry
    if (ss == 0 || sprop == 0) {
      response_style <- NULL
    } else {
    n_style <- .sim_planted_count(sprop, N)
    if (n_style >= N)
      stop("response_style$prop must leave at least one person without the ",
           "style; a common category weighting is absorbed by the fitted ",
           "thresholds rather than appearing as person-level response style")
    style_idx <- sample(N, n_style)
    mid <- m / 2
    dev2 <- ((0:m - mid) / mid)^2
    # Apply the style on the log-probability scale. Constructing exp(ss *
    # dev2) first overflows at perfectly finite strengths, even though the
    # normalised categorical probabilities remain well defined.
    log_w <- if (stype == "extreme") ss * dev2 else -ss * dev2
    style_prob <- function(i, p, shift = 0) {
      th_p <- if (i %in% dim_items) theta2[p] else theta[p]
      pp <- item_pars(i, p)
      cumulative <- c(0, cumsum(pp$tau))
      eta <- pp$disc * ((0:m) * (th_p + shift) - cumulative) + log_w
      if (any(!is.finite(cumulative)) || any(!is.finite(eta)))
        stop("the response-style strength and generating parameters are ",
             "outside the numerically representable range", call. = FALSE)
      eta <- eta - max(eta)
      pr <- exp(eta)
      pr / sum(pr)
    }
    for (p in style_idx) {
      # Draw the styled marginal responses first. Dependence is then rebuilt
      # in its declared order from these responses. Reusing dep_shift here
      # would condition a target on the source response drawn before the
      # source itself acquired its response style.
      for (i in seq_len(I)) {
        if (is.na(X[p, i])) next
        X[p, i] <- sample.int(m + 1L, 1L,
                              prob = style_prob(i, p)) - 1L
      }
      styled_shift <- numeric(I)
      for (ij in dep_pair_idx) {
        i1 <- ij[1L]; i2 <- ij[2L]
        if (is.na(X[p, i1]) || is.na(X[p, i2])) next
        source_prob <- style_prob(i1, p, styled_shift[i1])
        source_mean <- sum((0:m) * source_prob)
        styled_shift[i2] <- styled_shift[i2] +
          d_str * (X[p, i1] - source_mean) / m
        X[p, i2] <- sample.int(
          m + 1L, 1L, prob = style_prob(i2, p, styled_shift[i2])) - 1L
      }
    }
    }
  }

  # careless responders: answer uniformly at random
  careless_idx <- integer(0)
  if (careless > 0) {
    n_careless <- .sim_planted_count(careless, N)
    if (n_careless >= N)
      stop("careless must leave at least one person with model-based responses")
    # Careless responses replace the whole response vector, so overlap would
    # erase a requested response style and make its realised proportion
    # smaller than the argument states.
    ordinary <- setdiff(seq_len(N), style_idx)
    if (n_careless > length(ordinary))
      stop("the requested careless and response-style proportions cannot ",
           "coexist without overlap; reduce one of them")
    if (n_groups > 1L && n_careless) {
      # A person with a response style still follows the item and DIF model,
      # so only groups with no such protected person need an ordinary keeper.
      protected_groups <- unique(as.character(group[style_idx]))
      need_keeper <- setdiff(levels(group), protected_groups)
      keep <- vapply(need_keeper, function(g) {
        pool <- ordinary[as.character(group[ordinary]) == g]
        if (!length(pool)) return(NA_integer_)
        pool[sample.int(length(pool), 1L)]
      }, integer(1))
      if (anyNA(keep) || n_careless > length(setdiff(ordinary, keep)))
        stop("careless must leave at least one person with model-based ",
             "responses in every generated group")
      eligible <- setdiff(ordinary, keep)
      careless_idx <- eligible[sample.int(length(eligible), n_careless)]
    } else careless_idx <- if (n_careless)
      ordinary[sample.int(length(ordinary), n_careless)] else integer(0)
    X[careless_idx, ] <- matrix(sample(0:m, length(careless_idx) * I, TRUE),
                                length(careless_idx), I)
  }

  # speededness: a contiguous not-reached tail over the last items. `speeded`
  # persons drop out somewhere in the final zone, so the missing rate grows
  # linearly to `speeded` at the last item
  # the not-reached tail spans the last 40% of the test: with fewer than
  # three items there is no tail, so the request cannot be planted and must
  # not be recorded as though it were
  if (speeded > 0 && I < 3L) {
    warning("speededness needs at least three items to place a ",
            "not-reached tail; ignored")
    speeded <- 0
  }
  speeded_idx <- integer(0)
  if (speeded > 0 && I >= 3L) {
    k <- max(1L, round(0.4 * I)); z0 <- I - k
    speeded_idx <- sample(N, .sim_planted_count(speeded, N))
    for (p in speeded_idx) {
      stop_at <- z0 + sample.int(k, 1L) - 1L              # drop point in zone
      if (stop_at < I) X[p, (stop_at + 1L):I] <- NA
    }
  }
  missing_cells <- integer(0)
  if (missing > 0) {
    n_missing <- .sim_planted_count(missing, length(X))
    available <- which(!is.na(X))
    if (n_missing >= length(available))
      stop("the requested speededness and completely-at-random missingness ",
           "would remove every remaining response; reduce one of them")
    missing_cells <- available[sample.int(length(available), n_missing)]
    X[missing_cells] <- NA
  }

  out <- data.frame(id = sprintf("P%04d", seq_len(N)), X,
                    check.names = FALSE, stringsAsFactors = FALSE)
  if (!is.null(group)) out$group <- group

  planted <- character(0)
  common_disc <- length(unique(disc)) == 1L
  if (common_disc && disc[1L] != 1) planted <- c(planted, sprintf(
    "common logit-unit multiplier %.2f", disc[1L]))
  if (!common_disc) planted <- c(planted, sprintf(
    "heterogeneous discrimination: %s",
    paste(sprintf("%s = %.2f", inm, disc), collapse = ", ")))
  if (any(guess > 0)) planted <- c(planted, sprintf("guessing on %s",
    paste(inm[guess > 0], collapse = ", ")))
  if (length(dim_items)) {
    dim_rho <- second_dim$rho %||% 0.5
    planted <- c(planted, if (dim_rho == -1) sprintf(
      "exactly reversed trait on %s (polarity-reversal stress)",
      paste(inm[dim_items], collapse = ", ")) else sprintf(
      "second dimension (rho %.2f) on %s", dim_rho,
      paste(inm[dim_items], collapse = ", ")))
  }
  if (length(dep_pairs)) planted <- c(planted, sprintf(
    "response dependence: %s", paste(vapply(dep_pairs, paste, "",
                                            collapse = "-"), collapse = "; ")))
  if (length(dif_items)) planted <- c(planted, sprintf(
    "DIF (group %s) on %s: uniform %.2f, non-uniform %.2f", dif_grp,
    paste(inm[dif_items], collapse = ", "), dif$uniform %||% 0,
    dif$nonuniform %||% 0))
  if (length(careless_idx)) planted <- c(planted, sprintf(
    "%d careless responder(s)", length(careless_idx)))
  if (length(style_idx)) planted <- c(planted, sprintf(
    "%s response style on %d person(s)",
    response_style$type %||% "extreme", length(style_idx)))
  if (speeded > 0) planted <- c(planted, sprintf(
    "speededness (%.0f%% not-reached at the last item)", 100 * speeded))
  if (length(dis_items)) planted <- c(planted, sprintf(
    "disordered thresholds on %s", paste(inm[dis_items], collapse = ", ")))
  if (length(missing_cells)) planted <- c(planted, sprintf(
    "%d response%s missing completely at random (%.1f%%)",
    length(missing_cells), if (length(missing_cells) == 1L) "" else "s",
    100 * length(missing_cells) / length(X)))

  calibration_truth <- if (common_disc) list(
    unit_multiplier = disc[1L], theta = disc[1L] * theta,
    difficulty = disc[1L] * delta,
    thresholds = lapply(tau, function(z) disc[1L] * z)
  ) else NULL
  departure_types <- c(
    if (!common_disc) "heterogeneous item discriminations",
    if (any(guess > 0)) "guessing",
    if (length(dim_items)) "a second dimension",
    if (length(dep_pairs)) "local response dependence",
    if (length(dif_items)) "differential item functioning",
    if (length(careless_idx)) "careless responding",
    if (length(style_idx)) "response style",
    if (rsm_pattern_departure)
      "item-specific threshold structure")
  attr(out, "truth") <- list(
    layout = "rasch",
    description = sprintf("%s, %d persons x %d items%s", model, N, I,
      if (!is.null(group)) sprintf(", %d groups", nlevels(group)) else ""),
    model = model, n_persons = N, n_items = I,
    person_id = out$id, theta = theta, theta2 = theta2,
    difficulty = delta, thresholds = tau,
    discrimination = disc, guessing = guess,
    calibration_truth = calibration_truth,
    departure_types = departure_types,
    groups = group, dim_items = inm[dim_items], dif_items = inm[dif_items],
    careless_idx = careless_idx, style_idx = style_idx,
    speeded_idx = speeded_idx, missing_cells = missing_cells,
    planted = planted)
  class(out) <- c("rasch_sim", "data.frame")
  out
}

# category probabilities for one location (used by the dependence term)
.p_item <- function(theta, tau, disc = 1) {
  m <- length(tau); cum <- c(0, cumsum(tau))
  e <- disc * ((0:m) * theta - cum)
  if (any(!is.finite(cum)) || any(!is.finite(e)))
    stop("the generating locations, thresholds and discrimination are ",
         "outside the numerically representable range", call. = FALSE)
  e <- e - max(e)
  p <- exp(e); p / sum(p)
}

#' @export
print.rasch_sim <- function(x, ...) {
  tr <- attr(x, "truth")
  cat(sprintf("Simulated %s data: %s\n", tr$layout,
              tr$description %||% sprintf("%d rows", nrow(x))))
  if (length(tr$planted)) {
    cat("Generating values and departures:\n")
    for (p in tr$planted) cat(paste0("  - ", p, "\n"))
  } else cat("Default model-conforming settings (no departures planted).\n")
  invisible(x)
}

#' Simulate paired-comparison data
#'
#' Generates dichotomous or ordered paired comparisons from the
#' Bradley--Terry--Luce model. Optional arguments introduce a second object
#' attribute, erratic judges, or within-judge dependence. Generating values are
#' stored in \code{attr(x, "truth")}.
#'
#' @param n_objects,n_judges Objects to scale and judges comparing them. Every
#'   judge is allocated at least one comparison; the simulator refuses a
#'   design with fewer comparisons than judges.
#' @param reps_per_pair Comparisons made of each object pair.
#' @param model \code{"dichotomous"} (a winner) or \code{"polytomous"} (a rated
#'   margin in \code{n_categories} categories; an earlier development-era value
#'   \code{"graded"} is accepted as an alias).
#' @param n_categories Categories for the polytomous model.
#' @param object_sd Realised sample standard deviation of the object locations
#'   (evenly spaced and sum-zero).
#' @param object_locations Optional numeric vector of generated object
#'   locations. It must have length \code{n_objects}; names, when supplied,
#'   must identify the generated objects. Values are centred to identify the
#'   origin and take precedence over \code{object_sd}.
#' @param second_attribute \code{NULL}, or \code{list(rho=)}: half the judges
#'   rank by a second object attribute whose realised correlation with the
#'   first is \code{rho}. It lies in [-1, 1); at 1 the attributes are identical
#'   and no second attribute is planted.
#'   This introduces residual dimensionality and possible intransitivity.
#' @param erratic_judges Proportion of judges who choose at random. At least
#'   one judge must retain model-based comparisons, in each camp when a
#'   second attribute is generated.
#' @param dependence \code{NULL}, or \code{list(exposure=, carry_over=)}:
#'   within-judge order effects (a seen-before advantage and a pull from the
#'   judge's own earlier verdicts). Adds an \code{order} column. Feeds the
#'   dependence effects fitted by \code{\link{btl}}.
#' @param seed Optional non-negative whole-number RNG seed. See
#'   \code{\link{rasch_rng}} for generator support.
#' @return A data frame of class \code{"rasch_sim"}: \code{object_a},
#'   \code{object_b}, \code{winner} (or \code{response} when polytomous),
#'   \code{judge}, and \code{order} when dependence is planted; with
#'   \code{attr(x, "truth")}.
#' @examples
#' d <- simulate_btl(8, 12, erratic_judges = 0.15, seed = 1)
#' bt <- btl(d, "object_a", "object_b", winner = "winner", judge = "judge")
#' bt$judges          # the erratic judges carry large fit residuals
#' @export
simulate_btl <- function(n_objects = 8, n_judges = 12, reps_per_pair = 25,
                         model = c("dichotomous", "polytomous", "graded"),
                         n_categories = 4,
                         object_sd = 1, second_attribute = NULL,
                         erratic_judges = 0, dependence = NULL, seed = NULL,
                         object_locations = NULL) {
  seed <- .sim_seed(seed)
  if (!is.null(seed)) {
    .old_stream <- .sim_seed_capture()
    on.exit(.sim_seed_restore(.old_stream), add = TRUE)
    set.seed(seed)
  }
  model <- match.arg(model)
  # "graded" is an earlier development name for the polytomous comparison model,
  # kept as a working alias for released user code
  if (model == "graded") model <- "polytomous"
  K <- .sim_count(n_objects, "n_objects", 3L)
  J <- .sim_count(n_judges, "n_judges", 2L)
  reps_per_pair <- .sim_count(reps_per_pair, "reps_per_pair")
  object_sd <- .sim_scalar(object_sd, "object_sd", lower = 0)
  erratic_judges <- .sim_scalar(erratic_judges, "erratic_judges",
                                lower = 0, upper = 1)
  second_attribute <- .sim_structure(second_attribute, "second_attribute",
                                     "rho")
  dependence <- .sim_structure(dependence, "dependence",
                               c("exposure", "carry_over"))
  if (model == "polytomous")
    n_categories <- .sim_count(n_categories, "n_categories", 3L)
  m <- if (model == "polytomous") as.integer(n_categories) - 1L else 1L
  objs <- sprintf("O%d", seq_len(K)); jids <- sprintf("J%d", seq_len(J))
  if (is.null(object_locations)) {
    beta <- setNames(as.numeric(scale(seq_len(K))) * object_sd, objs)
    .sim_check_sd(beta, object_sd, "object-location standard deviation")
  } else {
    location_names <- names(object_locations)
    object_locations <- .sim_vector(object_locations, "object_locations", K)
    if (!is.null(location_names) &&
        (!setequal(location_names, objs) || anyDuplicated(location_names)))
      stop("named object_locations must identify every generated object exactly")
    if (!is.null(location_names))
      object_locations <- stats::setNames(object_locations, location_names)[objs]
    beta <- setNames(object_locations - mean(object_locations), objs)
  }
  tau <- if (m > 1L) .sim_thresholds(0, m, 1.2) else NULL

  # a second object attribute (orthogonal part) for the two-camp design
  beta2 <- NULL; camp <- NULL
  if (!is.null(second_attribute)) {
    beta_sd <- stats::sd(beta)
    if (!is.finite(beta_sd) || beta_sd <= sqrt(.Machine$double.eps))
      stop("second_attribute requires object locations with positive spread because an object correlation is otherwise undefined")
    rho <- second_attribute$rho %||% 0.3
    rho <- .sim_scalar(rho, "second_attribute$rho", lower = -1, upper = 1,
                       upper_open = TRUE)
    beta2 <- setNames(.sim_correlated_fixed(beta, rho, beta_sd), objs)
    camp <- setNames(rep(c("a", "b"), length.out = J), jids)
  }
  pr <- t(utils::combn(objs, 2))
  d <- data.frame(object_a = rep(pr[, 1], each = reps_per_pair),
                  object_b = rep(pr[, 2], each = reps_per_pair),
                  stringsAsFactors = FALSE)
  d$judge <- .sim_balanced_ids(jids, nrow(d), "judges")
  n_erratic <- .sim_planted_count(erratic_judges, J)
  if (n_erratic >= J)
    stop("erratic_judges must leave at least one judge with model-based comparisons")
  if (n_erratic && !is.null(camp)) {
    keep <- unname(vapply(unique(camp), function(g)
      sample(names(camp)[camp == g], 1L), character(1)))
    eligible <- setdiff(jids, keep)
    if (n_erratic > length(eligible))
      stop("erratic_judges must leave at least one model-based judge in ",
           "each second-attribute camp")
    erratic <- sample(eligible, n_erratic, replace = FALSE)
  } else erratic <- if (n_erratic)
    sample(jids, n_erratic, replace = FALSE) else character(0)

  win_prob <- function(a, b, jd) {
    ba <- if (!is.null(camp) && camp[jd] == "b") beta2[a] else beta[a]
    bb <- if (!is.null(camp) && camp[jd] == "b") beta2[b] else beta[b]
    ba - bb
  }
  # dependence needs a per-judge judgment order and running history
  if (!is.null(dependence)) {
    exq <- .sim_scalar(dependence$exposure %||% 0, "dependence$exposure")
    cry <- .sim_scalar(dependence$carry_over %||% 0,
                       "dependence$carry_over")
    if (exq == 0 && cry == 0) dependence <- NULL
  }
  if (!is.null(dependence)) {
    # Give every judge an independently randomised presentation order.  The
    # pair rows are assembled in a common lexicographic order above; merely
    # sorting those rows by judge leaves every judge with essentially the
    # same sequence.  That aliases a genuine order effect with object
    # structure and makes the paired-comparison dimensionality verdict
    # unidentifiable in data generated by the package itself.
    d <- d[order(d$judge, stats::runif(nrow(d))), ]
    d$order <- stats::ave(seq_len(nrow(d)), d$judge, FUN = seq_along)
    seen <- new.env(parent = emptyenv()); hs <- new.env(parent = emptyenv())
    hc <- new.env(parent = emptyenv())
    g0 <- function(e, k) if (is.null(v <- e[[k]])) 0 else v
    resp <- integer(nrow(d))
    for (r in seq_len(nrow(d))) {
      j <- d$judge[r]; a <- d$object_a[r]; b <- d$object_b[r]
      ka <- paste(j, a); kb <- paste(j, b)
      lp <- win_prob(a, b, j) +
        exq * (as.numeric(g0(seen, ka) > 0) - as.numeric(g0(seen, kb) > 0)) +
        cry * ((if (g0(hc, ka) > 0) g0(hs, ka) / g0(hc, ka) else 0) -
               (if (g0(hc, kb) > 0) g0(hs, kb) / g0(hc, kb) else 0))
      x <- if (j %in% erratic) sample(0:m, 1)
           else if (m == 1L) as.integer(stats::runif(1) < stats::plogis(lp))
           else sample(0:m, 1, prob = .p_item(lp, tau))
      resp[r] <- x
      assign(ka, g0(seen, ka) + 1, seen); assign(kb, g0(seen, kb) + 1, seen)
      assign(ka, g0(hc, ka) + 1, hc);     assign(kb, g0(hc, kb) + 1, hc)
      assign(ka, g0(hs, ka) + (2 * x / m - 1), hs)
      assign(kb, g0(hs, kb) + (2 * (m - x) / m - 1), hs)
    }
  } else {
    lp <- vapply(seq_len(nrow(d)), function(r)
      win_prob(d$object_a[r], d$object_b[r], d$judge[r]), 0)
    resp <- integer(nrow(d))
    reg <- !(d$judge %in% erratic)
    resp[reg] <- if (m == 1L) as.integer(stats::runif(sum(reg)) < stats::plogis(lp[reg]))
                 else vapply(which(reg), function(r) sample(0:m, 1, prob = .p_item(lp[r], tau)), 0L)
    if (any(!reg)) resp[!reg] <- sample(0:m, sum(!reg), TRUE)
  }

  if (m == 1L) d$winner <- ifelse(resp == 1L, d$object_a, d$object_b)
  else d$response <- resp
  rownames(d) <- NULL

  planted <- character(0)
  if (length(erratic)) planted <- c(planted,
    sprintf("%d erratic judge(s): %s", length(erratic), paste(erratic, collapse = ", ")))
  if (!is.null(second_attribute)) {
    attribute_rho <- second_attribute$rho %||% 0.3
    planted <- c(planted, if (attribute_rho == -1)
      "exactly reversed object attribute (two-camp reversal stress)" else
      sprintf("second object attribute (rho %.2f), two judge camps",
              attribute_rho))
  }
  if (!is.null(dependence)) planted <- c(planted, sprintf(
    "within-judge dependence: exposure %.2f, carry-over %.2f",
    dependence$exposure %||% 0, dependence$carry_over %||% 0))

  attr(d, "truth") <- list(
    layout = "btl",
    description = sprintf("%s, %d objects, %d judges, %d comparisons",
                          model, K, J, nrow(d)),
    model = model, location = beta, location2 = beta2,
    departure_types = c(
      if (!is.null(beta2)) "a second object attribute",
      if (length(erratic)) "erratic judges",
      if (!is.null(dependence)) "within-judge dependence"),
    dependence = dependence,
    attribute_correlation = if (is.null(beta2)) NULL else
      stats::cor(beta, beta2), camp = camp,
    erratic = erratic, planted = planted)
  class(d) <- c("rasch_sim", "data.frame")
  d
}

#' Simulate many-facet Rasch data
#'
#' Generates fully crossed ratings from a many-facet Rasch model (Linacre
#' 1989), with optional erratic raters, item-by-rater interaction, or halo. A
#' positive rater proportion selects at least one rater when the requested
#' departures can coexist.
#'
#' @param n_persons,n_items,n_raters Facet sizes (fully crossed).
#' @param n_categories Rating categories.
#' @param theta_sd Realised sample standard deviation of person ability.
#' @param item_sd Realised sample standard deviation of the deterministic item
#'   difficulties; zero gives equal item difficulties.
#' @param rater_severity_sd Realised sample standard deviation of rater
#'   severities (the core facet; recovered in \code{facet_effects}).
#' @param erratic_raters Proportion of raters who rate at random (feeds the
#'   rater fit residual). Erratic and halo raters are disjoint, and together
#'   must leave at least one ordinary rater.
#' @param interaction \code{NULL}, or \code{list(rater=, item=, bias=)}: one
#'   rater is unusually harsh (positive) or lenient (negative) on one item.
#'   Feeds the item-by-rater interaction (fit with \code{interaction = }).
#' @param halo Proportion of raters showing a halo effect: they rate by the
#'   person's overall level and barely differentiate items (feeds the rater
#'   fit residual and the item-by-rater interaction). Its requested count
#'   must fit among the non-erratic raters and leave an ordinary rater.
#'   A positive halo proportion requires \code{item_sd > 0}.
#' @param seed Optional non-negative whole-number RNG seed. See
#'   \code{\link{rasch_rng}} for generator support.
#' @return A long data frame of class \code{"rasch_sim"} (\code{person},
#'   \code{item}, \code{rater}, \code{score}) ready for
#'   \code{\link{rasch_mfrm}}, with the truth attached.
#' @examples
#' d <- simulate_mfrm(60, 5, 6, rater_severity_sd = 0.8, seed = 1)
#' mf <- rasch_mfrm(d, person = "person", item = "item", score = "score",
#'                  facets = "rater")
#' cor(mf$facet_effects$rater$severity, attr(d, "truth")$severity)  # recovered
#' @export
simulate_mfrm <- function(n_persons = 80, n_items = 5, n_raters = 6,
                          n_categories = 4, theta_sd = 1.2, item_sd = 1,
                          rater_severity_sd = 0.6, erratic_raters = 0,
                          interaction = NULL, halo = 0, seed = NULL) {
  seed <- .sim_seed(seed)
  if (!is.null(seed)) {
    .old_stream <- .sim_seed_capture()
    on.exit(.sim_seed_restore(.old_stream), add = TRUE)
    set.seed(seed)
  }
  N <- .sim_count(n_persons, "n_persons", 2L)
  I <- .sim_count(n_items, "n_items", 2L)
  R <- .sim_count(n_raters, "n_raters", 2L)
  n_categories <- .sim_count(n_categories, "n_categories", 2L)
  theta_sd <- .sim_scalar(theta_sd, "theta_sd", lower = 0)
  item_sd <- .sim_scalar(item_sd, "item_sd", lower = 0)
  rater_severity_sd <- .sim_scalar(rater_severity_sd, "rater_severity_sd",
                                   lower = 0)
  erratic_raters <- .sim_scalar(erratic_raters, "erratic_raters",
                                lower = 0, upper = 1)
  halo <- .sim_scalar(halo, "halo", lower = 0, upper = 1)
  if (halo > 0 && item_sd == 0)
    stop("halo requires item_sd > 0; when all item difficulties are equal, ",
         "ignoring their differences does not change any generated rating")
  interaction <- .sim_structure(interaction, "interaction",
                                c("rater", "item", "bias"),
                                c("rater", "item", "bias"))
  m <- n_categories - 1L
  pids <- sprintf("P%03d", seq_len(N)); iids <- sprintf("I%d", seq_len(I))
  rids <- sprintf("R%d", seq_len(R))
  theta <- .sim_theta(N, 0, theta_sd)
  delta <- setNames(if (item_sd == 0) rep(0, I) else
    as.numeric(scale(seq_len(I))) * item_sd, iids)
  .sim_check_sd(delta, item_sd, "item-difficulty standard deviation")
  lambda <- setNames(.sim_theta(R, 0, rater_severity_sd), rids)
  base_tau <- .sim_thresholds(0, m, 1.2)
  erratic <- if (erratic_raters > 0)
    rids[seq_len(.sim_planted_count(erratic_raters, R))] else character(0)
  if (length(erratic) >= R)
    stop("erratic_raters must leave at least one rater with model-based ratings")
  # Halo raters are drawn from the end and remain disjoint from erratic
  # raters, whose random scores would erase the halo mechanism.
  halo_r <- if (halo > 0) {
    pool <- setdiff(rev(rids), erratic)
    n_halo <- .sim_planted_count(halo, R)
    if (n_halo > length(pool))
      stop("the requested erratic-rater and halo proportions cannot coexist ",
           "without overlap; reduce one of them")
    pool[seq_len(n_halo)]
  } else character(0)
  if (length(union(erratic, halo_r)) >= R)
    stop("erratic_raters and halo must leave at least one ordinary rater ",
         "whose ratings retain the generated item difficulties")
  int_bias <- matrix(0, I, R, dimnames = list(iids, rids))
  if (!is.null(interaction)) {
    plain_level <- function(x)
      is.atomic(x) && is.null(dim(x)) && is.null(oldClass(x)) &&
        length(x) == 1L && !is.na(x)
    if (!plain_level(interaction$item) || !(interaction$item %in% iids) ||
        !plain_level(interaction$rater) || !(interaction$rater %in% rids))
      stop("interaction$item and interaction$rater must each name one generated level")
    interaction$bias <- .sim_scalar(interaction$bias, "interaction$bias")
    if (interaction$bias == 0) {
      interaction <- NULL
    } else {
      # an erratic rater answers at random, discarding the whole rating model
      # for that rater: a bias planted on one would not be in the data, while
      # the recorded truth would still claim it
      if (interaction$rater %in% erratic)
        stop("interaction$rater '", interaction$rater, "' is one of the ",
             "erratic raters, who answer at random: the planted bias would ",
             "not appear in the data. Nominate another rater, or lower ",
             "erratic_raters")
      int_bias[interaction$item, interaction$rater] <- interaction$bias
    }
  }

  grid <- expand.grid(p = seq_len(N), i = seq_len(I), r = seq_len(R))
  score <- integer(nrow(grid))
  for (i in seq_len(I)) for (r in seq_len(R)) {
    rows <- grid$i == i & grid$r == r
    # item difficulty and rater severity shift the person's thresholds; a halo
    # rater ignores the item's own difficulty (uses the mean instead)
    di <- if (rids[r] %in% halo_r) mean(delta) else delta[i]
    tau_ir <- base_tau + di + lambda[r] + int_bias[i, r]
    score[rows] <- if (rids[r] %in% erratic) sample(0:m, sum(rows), TRUE)
                   else .sim_item(theta[grid$p[rows]], tau_ir)
  }
  d <- data.frame(person = pids[grid$p], item = iids[grid$i],
                  rater = rids[grid$r], score = score, stringsAsFactors = FALSE)

  planted <- character(0)
  if (length(erratic)) planted <- c(planted,
    sprintf("%d erratic rater(s): %s", length(erratic), paste(erratic, collapse = ", ")))
  if (length(halo_r)) planted <- c(planted,
    sprintf("%d halo rater(s): %s", length(halo_r), paste(halo_r, collapse = ", ")))
  if (!is.null(interaction)) planted <- c(planted, sprintf(
    "rater-by-item bias: %s on %s (%.2f)", interaction$rater,
    interaction$item, interaction$bias))
  planted <- c(planted, sprintf("rater severities SD %.2f", stats::sd(lambda)))

  attr(d, "truth") <- list(
    layout = "mfrm",
    description = sprintf("%d persons x %d items x %d raters (%d ratings)",
                          N, I, R, nrow(d)),
    person_id = pids, theta = theta, difficulty = delta, severity = lambda,
    n_categories = n_categories,
    departure_types = c(if (length(erratic)) "erratic raters",
                        if (length(halo_r)) "halo ratings",
                        if (!is.null(interaction)) "item-by-rater interaction"),
    interaction = interaction,
    erratic = erratic, halo = halo_r, planted = planted)
  class(d) <- c("rasch_sim", "data.frame")
  d
}

#' Simulate extended frame-of-reference data with differing units
#'
#' Generates data whose latent unit differs across item-set by person-group
#' frames (Humphry 2005): a person in group g responding to an item in set s
#' does so at the frame unit rho = alpha_set * phi_group scaling the whole
#' exponent. The planted set- and group-unit ratios are recovered by
#' \code{\link{rasch_efrm}}. A positive careless-response or missingness
#' proportion selects at least one person or response cell.
#'
#' @param n_per_group Persons in each group.
#' @param items_per_set Items in each set.
#' @param n_sets,n_groups Numbers of item sets and person groups.
#' @param set_unit_ratio,group_unit_ratio Ratio of the last generated set or
#'   group unit to the first. Intermediate units are geometrically spaced;
#'   1 gives equal units and hence an ordinary Rasch fit.
#' @param n_categories Response categories per item: 2 (the default) gives
#'   dichotomous items; larger values give partial credit items whose
#'   evenly spaced thresholds are centred on the item locations, with the
#'   frame unit scaling the whole exponent as in the dichotomous case.
#' @param theta_sd Realised sample standard deviation of person ability. This
#'   must be positive when more than one item set is generated because relative
#'   set units are identified from person variation shared across sets.
#' @param item_drift Optional \code{list(items=, group=, shift=)}. The named
#'   item or items move by \code{shift} logits in one generated person group,
#'   violating item invariance across frames. At least one item in every
#'   affected set must remain invariant: shifting a set's items together is a
#'   frame-origin difference, not item drift within that set. A non-zero drift
#'   requires at least two person groups.
#' @param careless Proportion of persons whose complete response vectors are
#'   replaced by random category choices. At least one person in every group
#'   must retain model-based responses.
#' @param missing Proportion of response cells set missing completely at
#'   random after the responses are generated. It must leave at least one
#'   observed response.
#' @param seed Optional non-negative whole-number RNG seed. See
#'   \code{\link{rasch_rng}} for generator support.
#' @return A wide data frame of class \code{"rasch_sim"}, containing an ID,
#'   item columns, and group. Its truth attribute contains the item-set map
#'   required by \code{\link{rasch_efrm}}.
#' @examples
#' d <- simulate_efrm(200, 6, set_unit_ratio = 1.3, seed = 1)
#' tr <- attr(d, "truth")
#' ef <- rasch_efrm(d, item_sets = tr$item_sets, groups = "group",
#'                  id = "id", boot_reps = 0)    # point estimates only
#' ef$alpha_table   # planted ratio 1.3, recovered within small-sample noise
#' @export
simulate_efrm <- function(n_per_group = 300, items_per_set = 8, n_sets = 2,
                          n_groups = 2, set_unit_ratio = 1.3,
                          group_unit_ratio = 1, n_categories = 2,
                          theta_sd = 1.3, seed = NULL, item_drift = NULL,
                          careless = 0, missing = 0) {
  # recorded before validation assigns to the formals, which would make
  # missing() report every argument as supplied
  set_ratio_given <- !missing(set_unit_ratio)
  group_ratio_given <- !missing(group_unit_ratio)
  seed <- .sim_seed(seed)
  if (!is.null(seed)) {
    .old_stream <- .sim_seed_capture()
    on.exit(.sim_seed_restore(.old_stream), add = TRUE)
    set.seed(seed)
  }
  S <- .sim_count(n_sets, "n_sets")
  G <- .sim_count(n_groups, "n_groups")
  K <- .sim_count(items_per_set, "items_per_set", 2L)
  npg <- .sim_count(n_per_group, "n_per_group", 2L)
  n_categories <- .sim_count(n_categories, "n_categories", 2L)
  m <- as.integer(n_categories) - 1L
  set_unit_ratio <- .sim_scalar(set_unit_ratio, "set_unit_ratio",
                                lower = 0, lower_open = TRUE)
  group_unit_ratio <- .sim_scalar(group_unit_ratio, "group_unit_ratio",
                                  lower = 0, lower_open = TRUE)
  theta_sd <- .sim_scalar(theta_sd, "theta_sd", lower = 0)
  careless <- .sim_scalar(careless, "careless", lower = 0, upper = 1)
  missing <- .sim_scalar(missing, "missing", lower = 0, upper = 1,
                         upper_open = TRUE)
  item_drift <- .sim_structure(item_drift, "item_drift",
                               c("items", "group", "shift"),
                               c("items", "group", "shift"))
  # a ratio is a comparison between frames: with one frame there is nothing
  # to span, and the span would silently return a unit while the recorded
  # truth claimed the requested ratio
  # only an EXPLICIT request is refused: the defaults describe the ordinary
  # multi-frame design and resolve to the reference when there is one frame
  if (S == 1L) {
    if (set_ratio_given && set_unit_ratio != 1)
      stop("`set_unit_ratio` must be 1 when `n_sets` is 1: a unit ratio is ",
           "a comparison between sets, so none can be planted")
    set_unit_ratio <- 1
  }
  if (S > 1L && theta_sd == 0)
    stop("`theta_sd` must be positive when `n_sets` is greater than 1: ",
         "relative set units require person variation shared across sets")
  if (G == 1L) {
    if (group_ratio_given && group_unit_ratio != 1)
      stop("`group_unit_ratio` must be 1 when `n_groups` is 1: a unit ratio ",
           "is a comparison between groups, so none can be planted")
    group_unit_ratio <- 1
  }
  # Set and group units span the ratio geometrically, normalised to mean one.
  # Centre on the log scale before exponentiating: normalising an already
  # exponentiated sequence can underflow for a valid positive finite ratio.
  gspan <- function(ratio, n) {
    lu <- seq(0, log(ratio), length.out = n)
    u <- exp(lu - mean(lu))
    if (any(!is.finite(u)) || any(u <= 0))
      stop("the requested unit ratio is outside the numerically representable range")
    u
  }
  alpha <- gspan(set_unit_ratio, S)
  phi <- gspan(group_unit_ratio, G)
  set_items <- lapply(seq_len(S), function(s) sprintf("S%dI%02d", s, seq_len(K)))
  inm <- unlist(set_items)
  delta <- setNames(rep(seq(-1.5, 1.5, length.out = K), S), inm)
  set_of <- setNames(rep(seq_len(S), each = K), inm)

  tau_list <- lapply(inm, function(nm) .sim_thresholds(delta[nm], m, 0.8))
  names(tau_list) <- inm

  grp <- factor(rep(sprintf("g%d", seq_len(G)), each = npg))
  gidx <- rep(seq_len(G), each = npg)                  # group NUMBER per person
  names(phi) <- sprintf("g%d", seq_len(G))             # unit by group label
  drift <- rep(0, length(inm)); names(drift) <- inm
  drift_group <- NULL
  if (!is.null(item_drift)) {
    if (!is.atomic(item_drift$items) || !is.null(dim(item_drift$items)) ||
        !is.null(oldClass(item_drift$items)))
      stop("item_drift$items must be a plain vector of generated item names")
    drift_items <- as.character(item_drift$items)
    if (!length(drift_items) || anyNA(drift_items) ||
        any(!nzchar(drift_items)) || anyDuplicated(drift_items) ||
        any(!drift_items %in% inm))
      stop("item_drift$items must name one or more generated items exactly once")
    if (!is.atomic(item_drift$group) || !is.null(dim(item_drift$group)) ||
        !is.null(oldClass(item_drift$group)) ||
        length(item_drift$group) != 1L || is.na(item_drift$group) ||
        !as.character(item_drift$group) %in% levels(grp))
      stop("item_drift$group must name one generated group")
    item_drift$shift <- .sim_scalar(item_drift$shift, "item_drift$shift")
    if (item_drift$shift == 0) {
      item_drift <- NULL
    } else {
      if (G < 2L)
        stop("item_drift requires n_groups >= 2; with one group the shift is ",
             "absorbed by the item's fitted location rather than violating ",
             "invariance across frames")
      affected_sets <- unique(unname(set_of[drift_items]))
      no_reference <- affected_sets[vapply(affected_sets, function(s)
        all(names(set_of)[set_of == s] %in% drift_items), logical(1))]
      if (length(no_reference))
        stop("item_drift$items must leave at least one invariant item in ",
             "every affected set; shifting a set's items together changes ",
             "that frame's origin and is not identified as item drift")
      drift[drift_items] <- item_drift$shift
      drift_group <- as.character(item_drift$group)
    }
  }
  N <- length(grp); theta <- .sim_theta(N, 0, theta_sd)
  X <- matrix(NA_integer_, N, length(inm), dimnames = list(NULL, inm))
  for (col in seq_along(inm)) {
    # phi is indexed by group NUMBER: as.integer() on the factor would order
    # the levels as strings (g1, g10, g2, ...) and attach each unit to the
    # wrong group as soon as there are ten of them
    s <- set_of[inm[col]]
    rho <- alpha[s] * phi[gidx]                        # per-person unit
    shift <- if (is.null(drift_group)) rep(0, N) else
      ifelse(as.character(grp) == drift_group, drift[inm[col]], 0)
    if (m == 1L) {
      eta <- rho * (theta - delta[inm[col]] - shift)
      if (any(!is.finite(eta)))
        stop("the requested frame units and generating parameters are ",
             "outside the numerically representable range", call. = FALSE)
      X[, col] <- as.integer(stats::runif(N) < stats::plogis(eta))
    } else {
      ct <- cumsum(tau_list[[inm[col]]])
      E <- cbind(0, sapply(seq_len(m), function(k)
        rho * (k * theta - ct[k] - k * shift)))
      if (any(!is.finite(ct)) || any(!is.finite(E)))
        stop("the requested frame units and generating parameters are ",
             "outside the numerically representable range", call. = FALSE)
      P <- exp(E - apply(E, 1, max)); P <- P / rowSums(P)
      cum <- P %*% upper.tri(diag(m + 1L), diag = TRUE)
      X[, col] <- as.integer(rowSums(stats::runif(N) > cum))
    }
  }
  careless_idx <- integer(0)
  if (careless > 0) {
    n_careless <- .sim_planted_count(careless, N)
    if (n_careless >= N)
      stop("careless must leave at least one person with model-based responses")
    keep <- vapply(levels(grp), function(g)
      {
        pool <- which(as.character(grp) == g)
        pool[sample.int(length(pool), 1L)]
      }, integer(1))
    eligible <- setdiff(seq_len(N), keep)
    if (n_careless > length(eligible))
      stop("careless must leave at least one person with model-based ",
           "responses in every generated group")
    careless_idx <- eligible[sample.int(length(eligible), n_careless)]
    X[careless_idx, ] <- matrix(
      sample(0:m, length(careless_idx) * ncol(X), replace = TRUE),
      nrow = length(careless_idx), ncol = ncol(X))
  }
  missing_cells <- integer(0)
  if (missing > 0) {
    n_missing <- .sim_planted_count(missing, length(X))
    if (n_missing >= length(X))
      stop("the requested missingness would remove every response; reduce `missing`")
    missing_cells <- sample.int(
      length(X), n_missing)
    X[missing_cells] <- NA_integer_
  }
  out <- data.frame(id = sprintf("P%04d", seq_len(N)), X, group = grp,
                    check.names = FALSE, stringsAsFactors = FALSE)

  planted <- sprintf("set-unit ratio %.2f across %d sets", set_unit_ratio, S)
  if (group_unit_ratio != 1)
    planted <- c(planted, sprintf("group-unit ratio %.2f across %d groups",
                                  group_unit_ratio, G))
  if (!is.null(item_drift) && item_drift$shift != 0)
    planted <- c(planted, sprintf("item drift in %s: %s shifted %.2f logits",
      drift_group, paste(names(drift)[drift != 0], collapse = ", "),
      item_drift$shift))
  if (length(careless_idx)) planted <- c(planted,
    sprintf("%d careless responder(s)", length(careless_idx)))
  if (length(missing_cells)) planted <- c(planted,
    sprintf("%.1f%% responses missing completely at random",
            100 * length(missing_cells) / length(X)))
  attr(out, "truth") <- list(
    layout = "efrm",
    description = sprintf("%d persons, %d sets x %d groups, %d items (%d categories)",
                          N, S, G, length(inm), m + 1L),
    person_id = out$id, theta = theta, difficulty = delta,
    thresholds = if (m > 1L) tau_list else NULL,
    departure_types = c(if (!is.null(item_drift)) "item drift",
                        if (length(careless_idx)) "careless responding"),
    alpha = alpha, phi = phi,
    item_sets = setNames(set_items, sprintf("set%d", seq_len(S))),
    groups = grp, item_drift = if (is.null(item_drift)) NULL else
      list(items = names(drift)[drift != 0], group = drift_group,
           shift = item_drift$shift),
    careless_idx = careless_idx, missing_cells = missing_cells,
    planted = planted)
  class(out) <- c("rasch_sim", "data.frame")
  out
}

#' Replicate a simulation for Monte Carlo studies
#'
#' Calls one of the \code{simulate_*} functions \code{n} times with successive
#' seeds, returning the datasets as a list -- for power, Type-I, or
#' parameter-recovery studies.
#'
#' @param FUN A simulator, e.g. \code{\link{simulate_rasch}}.
#' @param n Number of datasets.
#' @param ... Arguments passed to \code{FUN} (the same each replicate).
#' @param seed Seed of the first replicate (each subsequent one increments it).
#'   See \code{\link{rasch_rng}} for generator support.
#' @return A list of class \code{"rasch_sim_batch"}, one simulated dataset per
#'   element.
#' @examples
#' # 8 datasets with a planted DIF item; how often is it flagged?
#' batch <- sim_replicate(simulate_rasch, 4, n_persons = 300, n_items = 8,
#'                        dif = list(items = "I05", uniform = 0.8), n_groups = 2,
#'                        seed = 1)
#' # sim_apply() is resilient: a replicate the estimator refuses (e.g. a
#' # small or disconnected draw) contributes NA instead of aborting the run
#' flagged <- sim_apply(batch, function(d)
#'   dif_anova(rasch(d, id = "id", factors = "group"))$summary$uniform_DIF[5])
#' mean(flagged, na.rm = TRUE)
#' @export
sim_replicate <- function(FUN, n, ..., seed = NULL) {
  if (!is.function(FUN)) stop("FUN must be a simulation function")
  n <- .sim_count(n, "n")
  .sim_seed_check()
  base <- if (is.null(seed)) sample.int(1e6, 1L)
          else .sim_seed(seed)
  if (base > .Machine$integer.max - n + 1L)
    stop("seed plus the requested replicate count exceeds the integer range")
  # Subtract before adding: at the upper integer boundary, `base + k - 1L`
  # overflows at the intermediate addition even when k is one and the final
  # seed is valid.
  reps <- lapply(seq_len(n), function(k) FUN(..., seed = base + (k - 1L)))
  structure(reps, class = "rasch_sim_batch", n = n,
            layout = attr(reps[[1]], "truth")$layout)
}

#' Apply a statistic across a simulation batch
#'
#' Applies \code{FUN} to each replicate of a \code{\link{sim_replicate}}
#' batch, catching replicates on which \code{FUN} errors -- for example a
#' small or disconnected draw the estimator refuses as unidentified -- so a
#' single failure does not abort the whole Monte-Carlo run. Failed
#' replicates contribute \code{NA}; the number of failures and the distinct
#' error messages are attached as attributes.
#'
#' @param batch A \code{"rasch_sim_batch"} from \code{\link{sim_replicate}}
#'   (or any list of datasets).
#' @param FUN A function of one dataset returning a scalar statistic.
#' @param ... Further arguments passed to \code{FUN}.
#' @return A vector of per-replicate statistics, with \code{NA} where the
#'   function failed. Attribute \code{n_failed} gives the failure count;
#'   \code{failure_messages} contains the distinct messages.
#' @examples
#' batch <- sim_replicate(simulate_rasch, 10, n_persons = 300, n_items = 8,
#'                        seed = 1)
#' psi <- sim_apply(batch, function(d) rasch(d)$psi$PSI)
#' mean(psi, na.rm = TRUE)
#' @export
sim_apply <- function(batch, FUN, ...) {
  if (!is.list(batch) || is.data.frame(batch) || !length(batch))
    stop("batch must be a non-empty list of simulated datasets")
  if (!is.function(FUN)) stop("FUN must be a function")
  res <- lapply(batch, function(d)
    tryCatch(list(ok = TRUE, v = FUN(d, ...)),
             error = function(e) list(ok = FALSE, v = NA, msg = conditionMessage(e))))
  # a one-column data frame and a list holding a vector both have length 1:
  # unlisting them later would expand one dataset into several results while
  # the failure count stayed at zero
  scalar <- function(v) !is.null(v) && is.atomic(v) && !is.list(v) &&
    length(v) == 1L
  valid <- vapply(res, function(r) isTRUE(r$ok) && scalar(r$v), TRUE)
  for (i in which(!valid & vapply(res, `[[`, logical(1), "ok"))) {
    res[[i]]$ok <- FALSE
    res[[i]]$msg <- "FUN must return one atomic scalar value"
  }
  ok <- vapply(res, `[[`, logical(1), "ok")
  vals <- lapply(res, function(r) {
    v <- r$v; if (!scalar(v)) NA else v[[1]]
  })
  out <- tryCatch(unlist(vals, use.names = FALSE),
                  error = function(e) vals)
  msgs <- unique(vapply(res[!ok], function(r) r$msg, ""))
  structure(out, n_failed = sum(!ok), failure_messages = msgs)
}

#' @export
print.rasch_sim_batch <- function(x, ...) {
  cat(sprintf("%d simulated %s datasets (Monte Carlo batch)\n",
              attr(x, "n"), attr(x, "layout")))
  cat("Each element is a simulated dataset; fit and summarise across them, e.g.\n")
  cat("  vapply(batch, function(d) <statistic of fit>, 0)\n")
  invisible(x)
}

# Recovery is meaningful only when the fitted responses are the responses that
# carry the generating truth.  These comparisons deliberately ignore row and
# column presentation, but retain the response values and every fitted design
# allocation.  Lightweight synthetic fits used to exercise the reporting
# helpers do not carry the source data and are left to the parameter-matching
# checks below; every fit returned by a public estimator does carry them.
.recovery_mismatch <- function(detail = NULL) {
  stop("`fit` was not fitted to the responses and design carried by `sim`",
       if (!is.null(detail)) paste0(": ", detail) else "",
       call. = FALSE)
}

.recovery_same_matrix <- function(x, y) {
  is.matrix(x) && is.matrix(y) && identical(dim(x), dim(y)) &&
    isTRUE(all.equal(unname(x), unname(y), tolerance = 0,
                     check.attributes = FALSE))
}

# A fit given no id column labels its persons by row number. Those labels are
# not identifiers: reading them as such would report a different person
# allocation, or an unmatchable group membership, for a fit of these very rows.
.recovery_person_ids <- function(person) {
  if (is.null(person) || !"id" %in% names(person)) return(NULL)
  id <- .role_text_values(person$id)
  if (identical(id, as.character(seq_along(id)))) NULL else id
}

# The fit does carry person labels, but they are row numbers: an allocation
# recorded per person cannot be checked against them either way.
.recovery_rownumber_ids <- function(person) {
  !is.null(person) && "id" %in% names(person) &&
    is.null(.recovery_person_ids(person))
}

# Row names are presentation, not data: a fitted data frame carries the row
# labels of whatever rows were passed and the prepared simulation carries
# none, so compare the response keys alone.
.recovery_row_keys <- function(x) {
  x <- as.matrix(x)
  unname(apply(x, 1L, function(z)
    paste(ifelse(is.na(z), "<NA>", format(z, scientific = FALSE,
                                           trim = TRUE)), collapse = "\034")))
}

.recovery_check_wide <- function(fit, sim, items, fitted) {
  if (is.null(fitted)) return(invisible(NULL))
  if (!is.data.frame(sim) || is.null(colnames(fitted)) ||
      anyDuplicated(colnames(fitted)) ||
      any(!items %in% names(sim)) ||
      any(!colnames(fitted) %in% items))
    .recovery_mismatch("the fitted items do not match the simulation")
  spec <- fit$refit_spec %||% list()
  prep_model <- if (inherits(fit, "rasch_efrm")) "PCM" else
    spec$model %||% fit$model %||% "PCM"
  # Prepare the whole simulated item set, as the estimator did. An item that
  # everyone answered identically is dropped in preparation, so the fit of
  # this simulation legitimately carries fewer items than were generated;
  # any other difference is a different item set.
  expected <- tryCatch(.prepare_X(
    sim[, items, drop = FALSE],
    na_codes = spec$na_codes %||% -1,
    model = prep_model,
    anchors = spec$anchors)$X, error = function(e) NULL)
  if (is.null(expected) || nrow(expected) != nrow(fitted))
    .recovery_mismatch("the fitted score structure differs from the simulation")
  if (!setequal(colnames(expected), colnames(fitted)))
    .recovery_mismatch("the fitted items do not match the simulation")
  expected <- expected[, colnames(fitted), drop = FALSE]

  fid <- .recovery_person_ids(fit$person)
  sid <- if ("id" %in% names(sim)) .role_text_values(sim$id) else NULL
  both_ids <- !is.null(fid) && !is.null(sid)
  aligned <- both_ids && length(fid) == nrow(fitted) &&
    length(sid) == nrow(expected) && !anyNA(fid) && !anyNA(sid) &&
    !anyDuplicated(fid) && !anyDuplicated(sid) && setequal(fid, sid)
  if (both_ids && !aligned)
    .recovery_mismatch("the fitted person allocation differs from the simulation")
  same <- if (aligned) {
    .recovery_same_matrix(fitted,
                          expected[match(fid, sid), , drop = FALSE])
  } else {
    identical(sort(.recovery_row_keys(fitted)),
              sort(.recovery_row_keys(expected)))
  }
  if (!same) .recovery_mismatch("at least one fitted response differs")
  invisible(NULL)
}

# A successful fit after category compression estimates a different response
# scale. Matching it to the same prepared data does not make its locations
# estimates of the original generating threshold means.
.recovery_check_categories <- function(fit, sim, tr) {
  if (is.null(fit$X) || is.null(fit$tau_list)) return(invisible(NULL))
  cells <- colnames(fit$X)
  true_m <- if (identical(tr$layout, "mfrm")) {
    # Older simulations omit the declared count. Their largest observed
    # category still gives a lower bound on the generating response scale.
    rep(if (!is.null(tr$n_categories)) tr$n_categories - 1L else
      max(sim$score, na.rm = TRUE), length(cells))
  } else {
    items <- if (identical(tr$layout, "efrm"))
      as.character(fit$virtual_map$item[
        match(cells, as.character(fit$virtual_map$vkey))]) else cells
    if (is.null(tr$thresholds) && identical(tr$layout, "efrm"))
      rep(1L, length(cells)) else lengths(tr$thresholds)[items]
  }
  fitted_m <- lengths(fit$tau_list)
  if (length(true_m) != length(fitted_m) || anyNA(true_m) ||
      any(fitted_m != true_m))
    stop("the fitted response categories differ from the generating scale; ",
         "parameter recovery is unavailable after category removal or merging",
         call. = FALSE)
  invisible(NULL)
}

.recovery_sort_records <- function(x, fields) {
  x <- x[, fields, drop = FALSE]
  for (nm in setdiff(fields, c("response", "score", "weight", "order")))
    x[[nm]] <- .role_text_values(x[[nm]])
  for (nm in intersect(fields, c("response", "score", "weight", "order")))
    x[[nm]] <- as.numeric(x[[nm]])
  x <- x[do.call(order, c(x, list(na.last = TRUE))), , drop = FALSE]
  rownames(x) <- NULL
  x
}

.recovery_same_records <- function(x, y, fields) {
  x <- .recovery_sort_records(x, fields)
  y <- .recovery_sort_records(y, fields)
  identical(names(x), names(y)) && identical(dim(x), dim(y)) &&
    isTRUE(all.equal(x, y, tolerance = 0, check.attributes = FALSE))
}

.recovery_mfrm_facet <- function(fit, tr) {
  fes <- fit$facet_effects
  overlap <- if (length(fes) && length(names(tr$severity)))
    vapply(fes, function(z) {
      if (!is.data.frame(z) ||
          !all(c("level", "severity") %in% names(z))) return(0L)
      sum(names(tr$severity) %in% as.character(z$level))
    }, integer(1)) else integer(0)
  named_rater <- match("rater", names(fes))
  best <- if (length(overlap) && max(overlap) > 0L)
    which(overlap == max(overlap)) else integer(0)
  fi <- if (!is.na(named_rater) && named_rater %in% best) {
    named_rater
  } else if (length(best) == 1L) best else integer(0)
  if (length(fi) != 1L)
    stop("the planted rater severity cannot be matched to one fitted ",
         "facet; recovery would be ambiguous")
  fi
}

.recovery_check_mfrm <- function(fit, sim, tr, fi) {
  if (is.null(fit$X) || is.null(fit$virtual_map)) return(invisible(NULL))
  if (!is.data.frame(sim) ||
      !all(c("person", "item", "score") %in% names(sim)))
    .recovery_mismatch("the multiple-ratings source columns are unavailable")
  vm <- fit$virtual_map
  facets <- names(fit$facet_effects)
  if (is.null(facets) || anyNA(facets) || any(!nzchar(facets)) ||
      !length(fi) || is.na(fi) ||
      !all(c("vkey", "item", facets) %in% names(vm)) ||
      is.null(fit$person$id))
    .recovery_mismatch("the fitted multiple-ratings design is incomplete")
  source_facets <- facets
  # The simulator calls its generated severity facet `rater`, but users may
  # rename that column before fitting. Prefer the fitted facet's own source
  # column and otherwise match the conventional simulation column by levels.
  rater_facet <- facets[fi]
  if (!rater_facet %in% names(sim)) source_facets[fi] <- "rater"
  if (!all(source_facets %in% names(sim)))
    .recovery_mismatch("a fitted facet allocation is absent from the simulation")
  jj <- match(colnames(fit$X), as.character(vm$vkey))
  if (anyNA(jj)) .recovery_mismatch("the fitted response-cell map is incomplete")
  got <- lapply(seq_len(ncol(fit$X)), function(j) {
    rows <- which(!is.na(fit$X[, j]))
    z <- data.frame(person = fit$person$id[rows], item = vm$item[jj[j]],
                    stringsAsFactors = FALSE)
    for (k in seq_along(facets))
      z[[paste0("facet_", k)]] <- vm[[facets[k]]][jj[j]]
    z$score <- fit$X[rows, j]
    z
  })
  got <- do.call(rbind, got)
  expected <- sim[, c("person", "item"), drop = FALSE]
  for (k in seq_along(source_facets))
    expected[[paste0("facet_", k)]] <- sim[[source_facets[k]]]
  expected$score <- sim$score
  fields <- c("person", "item", paste0("facet_", seq_along(facets)),
              "score")
  if (!.recovery_same_records(got, expected, fields))
    .recovery_mismatch("at least one rating or facet allocation differs")
  invisible(NULL)
}

.recovery_btl_records <- function(a, b, response, weight, judge, m,
                                  presented = FALSE, order = NULL,
                                  panel = NULL) {
  a <- .role_text_values(a); b <- .role_text_values(b)
  judge <- .role_text_values(judge)
  response <- as.numeric(response); weight <- as.numeric(weight)
  if (!presented) {
    swap <- a > b
    lo <- ifelse(swap, b, a); hi <- ifelse(swap, a, b)
    response[swap] <- m - response[swap]
    a <- lo; b <- hi
  }
  out <- data.frame(object_a = a, object_b = b, response = response,
                    weight = weight, judge = judge,
                    stringsAsFactors = FALSE)
  if (!is.null(order)) out$order <- as.numeric(order)
  if (!is.null(panel)) {
    # Panel names are presentation labels, whereas the allocation of judges
    # to panels is part of the fitted design.  Replace each name by a
    # collision-safe signature of its member judges so a consistent relabel
    # is accepted but a changed partition is not.
    panel <- .role_text_values(panel)
    if (length(panel) != length(judge) || anyNA(panel) || anyNA(judge))
      .recovery_mismatch("the judge/frame allocation is incomplete")
    jp <- unique(data.frame(judge = judge, panel = panel,
                            stringsAsFactors = FALSE))
    per_judge <- split(jp$panel, jp$judge)
    if (any(vapply(per_judge, function(z) length(unique(z)) != 1L,
                   logical(1))))
      .recovery_mismatch(
        "the judge/frame allocation assigns a judge to more than one panel")
    members <- split(jp$judge, jp$panel)
    signatures <- vapply(members, function(z) {
      paste(sort(.factor_keys(data.frame(judge = unique(z),
                                          stringsAsFactors = FALSE))),
            collapse = "")
    }, character(1))
    out$panel <- unname(signatures[panel])
  }
  fields <- names(out)
  # Count compression and row expansion are two representations of the same
  # comparison records. Sum the weights over otherwise identical rows.
  key_fields <- setdiff(fields, "weight")
  key <- .factor_keys(out[key_fields])
  uk <- unique(key)
  ans <- out[match(uk, key), key_fields, drop = FALSE]
  ans$weight <- as.numeric(rowsum(out$weight, match(key, uk),
                                  reorder = FALSE))
  # The construction above used sorted keys for the representatives. Re-sort
  # after attaching their summed weights so both sides have one canonical form.
  .recovery_sort_records(ans, c(key_fields, "weight"))
}

.recovery_btl_source <- function(fit) {
  # Current BTL fits retain the usable source rows before free boundary
  # objects are set aside.  Older saved fits have only `comparisons`, so keep
  # that representation as the compatibility fallback.
  fit$observed_comparisons %||% fit$comparisons
}

.recovery_check_btl <- function(fit, sim) {
  source <- .recovery_btl_source(fit)
  if (is.null(source)) return(invisible(NULL))
  if (!is.data.frame(sim) ||
      !all(c("object_a", "object_b", "judge") %in% names(sim)))
    .recovery_mismatch("the paired-comparison source columns are unavailable")
  m <- as.integer(fit$m)
  response <- if ("response" %in% names(sim)) {
    as.numeric(sim$response)
  } else if ("winner" %in% names(sim) && m == 1L) {
    ifelse(.role_text_values(sim$winner) == .role_text_values(sim$object_a),
           1, ifelse(.role_text_values(sim$winner) ==
                       .role_text_values(sim$object_b), 0, NA_real_))
  } else NULL
  if (is.null(response) || anyNA(response))
    .recovery_mismatch("the simulated comparison outcome is unavailable")
  presented <- !is.null(fit$dependence) ||
    isTRUE((fit$refit_spec %||% list())$position)
  ord_fit <- if (presented && "order" %in% names(source))
    source$order else NULL
  ord_sim <- if (!is.null(ord_fit) && "order" %in% names(sim)) sim$order else NULL
  if (!is.null(ord_fit) && is.null(ord_sim))
    .recovery_mismatch("the fitted comparison order is absent from the simulation")
  panel_fit <- if ("panel" %in% names(source))
    source$panel else NULL
  panel_sim <- if (!is.null(panel_fit) && "panel" %in% names(sim))
    sim$panel else NULL
  if (!is.null(panel_fit) && is.null(panel_sim))
    .recovery_mismatch("the fitted panel allocation is absent from the simulation")
  got <- .recovery_btl_records(
    source$object_a, source$object_b,
    source$response, source$weight,
    source$judge, m, presented, ord_fit, panel_fit)
  expected <- .recovery_btl_records(
    sim$object_a, sim$object_b, response, rep(1, nrow(sim)), sim$judge,
    m, presented, ord_sim, panel_sim)
  if (!.recovery_same_records(got, expected, names(got)))
    .recovery_mismatch("at least one comparison, outcome or judge/frame allocation differs")
  invisible(NULL)
}

.recovery_partition_labels <- function(old_membership, old_labels,
                                       fitted_membership, what) {
  old_membership <- .role_text_values(old_membership)
  fitted_membership <- .role_text_values(fitted_membership)
  if (length(old_membership) != length(fitted_membership) ||
      anyNA(old_membership) || anyNA(fitted_membership))
    stop("the planted ", what, " membership cannot be matched to the fit",
         call. = FALSE)
  mapped <- vapply(old_labels, function(lv) {
    z <- unique(fitted_membership[old_membership == lv])
    if (length(z) == 1L) z else NA_character_
  }, character(1))
  reverse_ok <- vapply(unique(fitted_membership), function(lv)
    length(unique(old_membership[fitted_membership == lv])) == 1L, logical(1))
  if (anyNA(mapped) || anyDuplicated(mapped) || !all(reverse_ok))
    stop("the fitted ", what, " partition does not match the simulated ",
         what, " membership; recovery would compare different estimands",
         call. = FALSE)
  unname(mapped)
}

.recovery_efrm_groups <- function(fit, tr) {
  old_labels <- names(tr$phi) %||% sprintf("g%d", seq_along(tr$phi))
  direct <- old_labels
  fid <- .recovery_person_ids(fit$person)
  # Row-number person labels cannot be read as identifiers, so a simulation
  # that records who belongs to which group offers nothing to check the
  # fitted allocation against. Matching the group labels instead would
  # report units for whichever people the fit grouped: withhold.
  if (.recovery_rownumber_ids(fit$person) &&
      !is.null(tr$groups) && !is.null(tr$person_id)) return(NULL)
  if (is.null(tr$groups) || is.null(tr$person_id) ||
      is.null(fid) || is.null(fit$frame_group) ||
      is.null(fit$factors) || !fit$frame_group[1L] %in% names(fit$factors)) {
    if (all(direct %in% fit$phi_table$group)) return(direct)
    stop("the planted group units cannot be matched to the fitted group labels; ",
         "recovery needs the person-to-group membership", call. = FALSE)
  }
  at <- match(fid, .role_text_values(tr$person_id))
  if (anyNA(at) || anyDuplicated(fid) || length(fid) != length(tr$person_id))
    stop("the planted person groups cannot be matched to the fitted people",
         call. = FALSE)
  .recovery_partition_labels(as.character(tr$groups)[at], old_labels,
    fit$factors[[fit$frame_group[1L]]], "person-group")
}

.recovery_btl_panels <- function(fit, tr) {
  old_labels <- names(tr$phi) %||% sprintf("panel%d", seq_along(tr$phi))
  if (is.null(tr$judge_panel)) {
    if (all(old_labels %in% fit$phi_table$panel)) return(old_labels)
    stop("the planted panel units cannot be matched to the fitted panel labels; ",
         "this simulation truth does not record judge-to-panel membership",
         call. = FALSE)
  }
  cmp <- .recovery_btl_source(fit)
  if (is.null(cmp) || !all(c("judge", "panel") %in% names(cmp)))
    stop("the fitted judge-to-panel membership is unavailable", call. = FALSE)
  fp <- split(.role_text_values(cmp$panel), .role_text_values(cmp$judge))
  if (any(vapply(fp, function(z) length(unique(z)) != 1L, logical(1))))
    stop("a fitted judge belongs to more than one panel", call. = FALSE)
  fitted <- vapply(fp, function(z) unique(z)[1L], character(1))
  old <- tr$judge_panel
  at <- match(names(fitted), names(old))
  if (anyNA(at) || length(fitted) != length(old))
    stop("the planted panel membership cannot be matched to the fitted judges",
         call. = FALSE)
  .recovery_partition_labels(unname(old[at]), old_labels, fitted,
                             "judge-panel")
}

#' Compare fitted and generating parameters
#'
#' Compares fitted parameters with the generating values from a
#' \code{simulate_*} function (carried on the data as
#' \code{attr(sim, "truth")}): item difficulties and person abilities for a
#' Rasch fit, object locations for a paired-comparison fit, rater severities
#' (with item and person measures) for a many-facet fit, set and group units
#' for a Rasch frames fit, and common object locations, panel and set units,
#' and set origins for a paired-comparison frames fit. Locations are
#' mean-centred where the model identifies them only up to an origin. An
#' externally anchored Rasch or paired-comparison fit retains its identified
#' origin, so recovery and bias are reported on the anchored scale.
#' The fit must be from these simulated responses and model family, and must
#' have converged. Row and item order do not matter; response values and the
#' person, judge, facet and frame allocations do. A fit given no \code{id}
#' column labels its persons by row number, which are not identifiers: the
#' responses are still verified, but person recovery and the EFRM group units
#' are withheld with a note rather than compared against an unverified
#' person allocation; the EFRM set units are still reported, with the note
#' recording that the fitted group allocation is unchecked.
#' Pass \code{id =} when fitting to recover them.
#' An item the estimator dropped, such as one everyone answered identically,
#' is named in the note and left out of the comparison.
#' Recovery is unavailable when fitting removes or merges generating response
#' categories, because the fitted locations then describe a different scale.
#' For a many-facet simulation, the planted rater facet must be identifiable
#' uniquely by its name or level labels.
#' EFRM set parameters are matched by their item or object membership, not by
#' the spelling of the set labels; a different fitted partition is refused.
#' Paired-comparison frames align generating origins to the fitted reference
#' set. That set must have a generating unit of one: fixing another unit to
#' one imposes a different cross-set restriction, so recovery is refused.
#' A common generating discrimination changes only the logit unit, so ordinary
#' Rasch recovery uses the equivalently rescaled item thresholds and person
#' locations. When the generator includes a departure that the fitted model
#' does not represent, the comparisons are labelled as descriptive rather
#' than as recovery of a single correctly specified target. This includes,
#' for example, heterogeneous item discriminations in an ordinary Rasch fit.
#'
#' @param fit A fit of the simulated data (\code{\link{rasch}},
#'   \code{\link{btl}}, \code{\link{rasch_mfrm}}, \code{\link{rasch_efrm}},
#'   or \code{\link{btl_efrm}}).
#' @param sim The simulated data (from a \code{simulate_*} function).
#' @return A list of class \code{"rasch_recovery"}: \code{summary} (per
#'   parameter type: n, correlation, RMSE, bias) and \code{pieces} (the true
#'   and estimated values behind each). \code{note} identifies unrepresented
#'   generating departures, an unverifiable original response scale,
#'   generated items the fit does not estimate, and comparisons withheld
#'   because the fit carries no person identifiers.
#' @examples
#' d <- simulate_rasch(500, 12, seed = 1)
#' sim_recovery(rasch(d, id = "id"), d)$summary
#' @export
sim_recovery <- function(fit, sim) {
  tr <- attr(sim, "truth")
  if (is.null(tr)) stop("`sim` carries no simulation truth")
  lay <- tr$layout
  if (!is.character(lay) || length(lay) != 1L || is.na(lay) ||
      !lay %in% c("rasch", "btl", "mfrm", "efrm", "btl_efrm"))
    stop("`sim` carries an unsupported or malformed simulation layout")
  family_ok <- switch(lay,
    rasch = inherits(fit, "rasch") &&
      !inherits(fit, c("rasch_btl", "rasch_mfrm", "rasch_efrm")),
    btl = inherits(fit, "rasch_btl") &&
      !inherits(fit, "rasch_btl_efrm"),
    mfrm = inherits(fit, "rasch_mfrm"),
    efrm = inherits(fit, "rasch_efrm") &&
      !inherits(fit, "rasch_btl_efrm"),
    btl_efrm = inherits(fit, "rasch_btl_efrm"))
  if (!isTRUE(family_ok))
    stop("`fit` does not match the simulated ", lay, " model family")
  converged <- if (inherits(fit, "rasch_btl")) fit$converged else
    fit$est$converged
  if (!isTRUE(converged))
    stop("`fit` did not converge; recovery summaries are unavailable")
  rasch_calibration_truth <- if (lay == "rasch") tr$calibration_truth else NULL
  rasch_disc <- if (lay == "rasch")
    as.numeric(tr$discrimination %||% 1) else numeric(0)
  # Older saved simulations do not carry calibration_truth, but do carry the
  # generating discrimination. Reconstruct the equivalent Rasch-unit target
  # so their recovery summaries are not silently compared on the wrong unit.
  if (lay == "rasch" && is.null(rasch_calibration_truth) &&
      length(unique(rasch_disc)) == 1L && is.finite(rasch_disc[1L]) &&
      rasch_disc[1L] > 0) {
    d <- rasch_disc[1L]
    rasch_calibration_truth <- list(
      unit_multiplier = d, theta = d * tr$theta,
      difficulty = d * tr$difficulty,
      thresholds = if (is.null(tr$thresholds)) NULL else
        lapply(tr$thresholds, function(z) d * z))
  }
  heterogeneous_slopes <- lay == "rasch" &&
    length(unique(rasch_disc)) > 1L
  departure_types <- tr$departure_types %||% character(0)
  if (!is.character(departure_types) || !is.null(dim(departure_types)) ||
      !is.null(oldClass(departure_types)) || anyNA(departure_types) ||
      any(!nzchar(departure_types)))
    stop("`sim` carries malformed departure labels", call. = FALSE)
  # A PCM generator's item-specific threshold patterns conform to the PCM and
  # are not labelled as a planted departure in the simulation itself. They do
  # become a fitted-model departure if those responses are instead restricted
  # to one rating-scale threshold pattern.
  if (lay == "rasch" && identical(tr$model, "PCM") &&
      is.list(tr$thresholds) && length(tr$thresholds) > 1L) {
    relative_thresholds <- lapply(tr$thresholds, function(z)
      if (is.numeric(z) && length(z)) z - mean(z) else z)
    if (length(unique(relative_thresholds)) > 1L)
      departure_types <- unique(c(departure_types,
                                  "item-specific threshold structure"))
  }
  # Simulations saved before departure_types was introduced still carry the
  # item discriminations, so retain the earlier misspecification warning.
  if (heterogeneous_slopes)
    departure_types <- unique(c(departure_types,
                                "heterogeneous item discriminations"))
  recovery_theta <- if (!is.null(rasch_calibration_truth))
    rasch_calibration_truth$theta else tr$theta
  mfrm_facet <- if (lay == "mfrm") .recovery_mfrm_facet(fit, tr) else NULL
  mfrm_true_difficulty <- if (lay == "mfrm") tr$difficulty else NULL
  mfrm_true_severity <- if (lay == "mfrm") tr$severity else NULL
  represented <- character(0)
  if (lay == "rasch" && "item-specific threshold structure" %in%
      departure_types && identical(fit$model, "PCM"))
    represented <- c(represented, "item-specific threshold structure")
  if (lay == "btl" && "within-judge dependence" %in% departure_types) {
    requested <- names(tr$dependence)[vapply(
      tr$dependence, function(z) is.numeric(z) && length(z) == 1L &&
        is.finite(z) && z != 0, logical(1))]
    fitted <- if (is.data.frame(fit$dependence))
      as.character(fit$dependence$effect) else character(0)
    if (length(requested) && all(requested %in% fitted))
      represented <- c(represented, "within-judge dependence")
  }
  mfrm_interaction_fitted <- lay == "mfrm" &&
    identical(fit$interaction, names(fit$facet_effects)[mfrm_facet])
  if (mfrm_interaction_fitted) {
    represented <- c(represented,
      intersect(c("item-by-rater interaction", "halo ratings"),
                departure_types))
    # The interaction fit uses sum-zero item, rater and interaction effects.
    # A single planted cell bias, and the item flattening imposed on halo
    # raters, therefore contribute to the fitted main effects as well as to
    # the interaction. Re-express the generating response-cell surface under
    # those same constraints before comparing the main effects.
    surface <- outer(as.numeric(tr$difficulty), as.numeric(tr$severity), `+`)
    rownames(surface) <- names(tr$difficulty)
    colnames(surface) <- names(tr$severity)
    halo <- intersect(as.character(tr$halo %||% character(0)),
                      colnames(surface))
    if (length(halo))
      surface[, halo] <- matrix(
        mean(tr$difficulty) + tr$severity[halo],
        nrow = nrow(surface), ncol = length(halo), byrow = TRUE)
    intr <- tr$interaction
    if (is.list(intr) && length(intr$item) == 1L &&
        length(intr$rater) == 1L && length(intr$bias) == 1L &&
        is.numeric(intr$bias) &&
        intr$item %in% rownames(surface) &&
        intr$rater %in% colnames(surface) && is.finite(intr$bias))
      surface[intr$item, intr$rater] <-
        surface[intr$item, intr$rater] + intr$bias
    origin <- mean(surface)
    mfrm_true_difficulty <- stats::setNames(rowMeans(surface) - origin,
                                             rownames(surface))
    mfrm_true_severity <- stats::setNames(colMeans(surface) - origin,
                                           colnames(surface))
  }
  if (lay == "rasch" && "a fixed explanatory item departure" %in%
      departure_types) {
    relaxed <- if (inherits(fit, "rasch_explanatory") &&
                   is.data.frame(fit$explanatory$relaxations))
      as.character(fit$explanatory$relaxations$item[
        fit$explanatory$relaxations$component == "Item location"]) else
          character(0)
    target <- tr$explanatory_departure$target %||% NA_character_
    target_relaxed <- is.character(target) && length(target) == 1L &&
      !is.na(target) && target %in% relaxed
    if (!inherits(fit, "rasch_explanatory") || target_relaxed)
      represented <- c(represented, "a fixed explanatory item departure")
  }
  if (lay == "btl" && "a fixed explanatory object departure" %in%
      departure_types) {
    relaxed <- if (inherits(fit, "rasch_btl_explanatory") &&
                   is.data.frame(fit$explanatory$relaxations))
      as.character(fit$explanatory$relaxations$object) else character(0)
    target <- tr$explanatory_departure$target %||% NA_character_
    target_relaxed <- is.character(target) && length(target) == 1L &&
      !is.na(target) && target %in% relaxed
    if (!inherits(fit, "rasch_btl_explanatory") || target_relaxed)
      represented <- c(represented, "a fixed explanatory object departure")
  }
  unrepresented <- setdiff(departure_types, represented)
  recovery_note <- if (length(unrepresented)) paste0(
    "The generating process includes ",
    paste(unrepresented, collapse = ", "),
    ", which the fitted model does not represent. These comparisons are ",
    "descriptive and are not unbiased recovery targets.") else NULL
  if (lay == "mfrm" && is.null(tr$n_categories))
    recovery_note <- paste(c(recovery_note,
      "The original category count was not recorded in this simulation, so response-scale equivalence cannot be verified."),
      collapse = " ")
  if (lay == "rasch") {
    .recovery_check_wide(fit, sim, names(tr$difficulty), fit$X)
  } else if (lay == "btl") {
    .recovery_check_btl(fit, sim)
  } else if (lay == "mfrm") {
    .recovery_check_mfrm(fit, sim, tr, mfrm_facet)
  } else if (lay == "efrm" && !is.null(fit$X)) {
    source <- tryCatch(.efrm_source_matrix(fit, names(fit$set_of)),
                       error = function(e) NULL)
    if (is.null(source))
      .recovery_mismatch("the fitted frame-to-item responses are unavailable")
    .recovery_check_wide(fit, sim, names(tr$difficulty), source)
  } else if (lay == "btl_efrm") {
    .recovery_check_btl(fit, sim)
  }
  if (lay %in% c("rasch", "mfrm", "efrm"))
    .recovery_check_categories(fit, sim, tr)
  pieces <- list(); centred <- list()
  add <- function(name, true, est, label = NULL, centre = FALSE) {
    true <- as.numeric(true); est <- as.numeric(est)
    keep <- is.finite(true) & is.finite(est)
    if (!any(keep)) return(invisible())
    if (centre) {                       # location-type: identified up to origin
      true <- true - mean(true[keep]); est <- est - mean(est[keep])
    }
    centred[[name]] <<- isTRUE(centre)
    pieces[[name]] <<- data.frame(
      parameter = name,
      label = if (is.null(label)) NA_character_ else as.character(label)[keep],
      true = true[keep], estimated = est[keep], stringsAsFactors = FALSE)
  }
  add_person <- function(centre = TRUE, name = "person ability") {
    true_id <- tr$person_id %||% names(recovery_theta)
    est_id <- if (!is.null(fit$person) && "id" %in% names(fit$person))
      as.character(fit$person$id) else NULL
    if (!is.null(true_id) && !is.null(est_id)) {
      at <- match(est_id, as.character(true_id))
      keep <- !is.na(at)
      if (any(keep)) {
        add(name, recovery_theta[at[keep]], fit$person$theta[keep],
            est_id[keep], centre = centre)
      } else {
        # A fit given no id column labels its persons by row number, so no
        # fitted person can be paired with the person who generated the row.
        recovery_note <<- paste(c(recovery_note, paste0(
          "No fitted person label matches a simulated person identifier, ",
          "as happens when the fit was given no id= column and labels its ",
          "persons by row number, so person recovery is unavailable.")),
          collapse = " ")
      }
    } else if (!is.null(fit$person) &&
               length(recovery_theta) == length(fit$person$theta)) {
      # Compatibility with simulation objects created before person IDs were
      # recorded in their truth attributes.
      add(name, recovery_theta, fit$person$theta, centre = centre)
    }
  }
  # Set labels are presentation metadata. A caller may give the same item or
  # object partition different names when fitting the simulated data. Recover
  # the fitted label from membership rather than assuming the simulator's
  # conventional set1, set2, ... names survived unchanged.
  fitted_set_labels <- function(truth_sets) {
    fallback <- names(truth_sets)
    if (is.null(fallback)) fallback <- sprintf("set%d", seq_along(truth_sets))
    fitted_map <- fit$set_of
    if ((is.null(fitted_map) || is.null(names(fitted_map))) &&
        is.data.frame(fit$objects) &&
        all(c("object", "set") %in% names(fit$objects)))
      fitted_map <- stats::setNames(as.character(fit$objects$set),
                                    as.character(fit$objects$object))
    if (is.null(fitted_map) || is.null(names(fitted_map))) return(fallback)
    mapped <- vapply(truth_sets, function(members) {
      z <- fitted_map[as.character(members)]
      if (length(z) != length(members) || anyNA(z)) return(NA_character_)
      z <- unique(as.character(z))
      if (length(z) == 1L) z else NA_character_
    }, character(1))
    if (anyNA(mapped) || anyDuplicated(mapped))
      stop("the fitted set partition does not match the simulated set ",
           "membership; set-parameter recovery would compare different ",
           "estimands")
    unname(mapped)
  }
  if (lay == "rasch") {
    anchored <- .has_calibration_anchors(fit)
    ei <- setNames(fit$items$location, fit$items$item)
    true_difficulty <- if (!is.null(rasch_calibration_truth))
      rasch_calibration_truth$difficulty else tr$difficulty
    cm <- intersect(names(true_difficulty), names(ei))
    # An item the estimator dropped (a constant item carries no information
    # about its difficulty) has no estimate to compare: name it rather than
    # reporting a quietly shorter item panel.
    unestimated <- setdiff(names(true_difficulty), cm)
    if (length(unestimated))
      recovery_note <- paste(c(recovery_note, paste0(
        "The fit does not estimate generated item(s) ",
        paste(unestimated, collapse = ", "),
        ", so the comparison omits them.")), collapse = " ")
    item_parameter <- if (heterogeneous_slopes)
      "item location (misspecified slopes)" else "item difficulty"
    person_parameter <- if (heterogeneous_slopes)
      "person ability (misspecified slopes)" else "person ability"
    add(item_parameter, true_difficulty[cm], ei[cm], cm,
        centre = !anchored)
    add_person(centre = !anchored, name = person_parameter)
  } else if (lay == "btl") {
    ot <- fit$objects
    # recovery is judged on calibrated locations; an extrapolated boundary
    # row is a reporting value, not an estimate of the planted location
    if ("extreme" %in% names(ot)) ot <- ot[!(ot$extreme %in% TRUE), ]
    eo <- setNames(ot$location, ot$object)
    cm <- intersect(names(tr$location), names(eo))
    add("object location", tr$location[cm], eo[cm], cm,
        centre = is.null(fit$anchors))
  } else if (lay == "mfrm") {
    # The simulated severity belongs to the rater facet.  Do not assume that
    # this is the first fitted facet: callers may add or reorder facets before
    # passing the fit here.  Prefer the simulator's conventional facet name;
    # otherwise use the single facet whose levels best match the planted
    # rater labels.  An ambiguous match is not a valid recovery comparison.
    fes <- fit$facet_effects
    fe <- fes[[mfrm_facet]]
    es <- setNames(fe$severity, fe$level)
    cm <- intersect(names(mfrm_true_severity), names(es))
    add("rater severity", mfrm_true_severity[cm], es[cm], cm, centre = TRUE)
    # per-item margins live in item_effects (fit$items holds the virtual
    # item-by-facet combinations, whose names never match)
    ie <- fit$item_effects
    if (!is.null(ie)) {
      ei <- setNames(ie$location, ie$item)
      ci <- intersect(names(mfrm_true_difficulty), names(ei))
      add("item difficulty", mfrm_true_difficulty[ci], ei[ci], ci,
          centre = TRUE)
    }
    add_person()
  } else if (lay == "efrm") {
    at <- fit$alpha_table
    set_labels <- fitted_set_labels(tr$item_sets)
    ealpha <- at$alpha[match(set_labels, at$set)]
    if (length(ealpha) != length(tr$alpha) || any(!is.finite(ealpha)))
      stop("the planted set units cannot be matched to finite fitted set units",
           call. = FALSE)
    # units are identified up to a common scale, so compare on the centred
    # log scale (a ratio); the planted alpha is normalised the same way
    add("set unit (log)", log(tr$alpha), log(ealpha), set_labels,
        centre = TRUE)
    # the person-group units phi are a fitted, reported quantity too --
    # recover them, not only the set units
    if (!is.null(tr$phi) && !is.null(fit$phi_table)) {
      glab <- .recovery_efrm_groups(fit, tr)
      if (is.null(glab)) {
        recovery_note <- paste(c(recovery_note, paste0(
          "The fitted persons are labelled by row number, so the simulated ",
          "person-to-group membership cannot be verified (pass id= when ",
          "fitting): the group units are withheld, and the set units are ",
          "those of a fit whose group allocation is unchecked.")),
          collapse = " ")
      } else {
        ephi <- fit$phi_table$phi[match(glab, fit$phi_table$group)]
        if (length(ephi) != length(tr$phi) || any(!is.finite(ephi)))
          stop("the planted group units cannot be matched to finite fitted group units",
               call. = FALSE)
        add("group unit (log)", log(tr$phi), log(ephi), glab, centre = TRUE)
      }
    }
  } else if (lay == "btl_efrm") {
    set_labels <- fitted_set_labels(tr$object_sets)
    reference <- fit$reference_set
    # Older reporting fixtures can omit the reference field or rename only
    # the tables. The first alpha row is the estimator's reference row.
    if (!is.character(reference) || length(reference) != 1L ||
        is.na(reference) || !reference %in% set_labels)
      reference <- as.character(fit$alpha_table$set[1L])
    ref <- match(reference, set_labels)
    if (length(ref) != 1L || is.na(ref) ||
        !is.finite(tr$alpha[ref]) || !is.finite(tr$kappa[ref]))
      stop("the fitted reference set cannot be matched to the generating parameters",
           call. = FALSE)
    # The common origin can move freely, but alpha is not a free scale
    # convention once the within-set calibrations and panel units are fixed.
    # Relabelling sets can change the alphabetically selected reference.
    if (unname(tr$alpha[ref]) != 1)
      stop("the fitted reference set has a generating unit different from one; ",
           "fixing it to one imposes a different cross-set restriction. ",
           "Keep the generating reference set alphabetically first for recovery",
           call. = FALSE)
    true_v <- tr$v - unname(tr$kappa[ref])
    true_kappa <- tr$kappa - unname(tr$kappa[ref])
    ot <- fit$objects
    ev <- setNames(ot$v %||% ot$location, ot$object)
    cm <- intersect(names(tr$v), names(ev))
    add("object location", true_v[cm], ev[cm], cm)

    panel_labels <- .recovery_btl_panels(fit, tr)
    ep <- setNames(fit$phi_table$phi, fit$phi_table$panel)
    if (length(panel_labels) != length(tr$phi) ||
        any(!is.finite(ep[panel_labels])))
      stop("the planted panel units cannot be matched to finite fitted panel units",
           call. = FALSE)
    add("panel unit (log)", log(unname(tr$phi)), log(ep[panel_labels]),
        panel_labels)

    ea <- setNames(fit$alpha_table$alpha, fit$alpha_table$set)
    if (length(set_labels) != length(tr$alpha) ||
        any(!is.finite(ea[set_labels])))
      stop("the planted set units cannot be matched to finite fitted set units",
           call. = FALSE)
    add("set unit (log)", log(unname(tr$alpha)), log(ea[set_labels]),
        set_labels)

    ek <- setNames(fit$kappa_table$kappa, fit$kappa_table$set)
    if (length(set_labels) != length(tr$kappa) ||
        any(!is.finite(ek[set_labels])))
      stop("the planted set origins cannot be matched to finite fitted set origins",
           call. = FALSE)
    add("set origin", unname(true_kappa), ek[set_labels], set_labels)
  } else stop("unsupported layout: ", lay)

  # bias after centring is structurally zero for any parameter identified
  # only up to an origin/scale convention: it is not identifiable, so report
  # NA rather than a misleading ~0. Correlation and RMSE (scatter about the
  # aligned scale) remain meaningful.
  if (!length(pieces))
    stop("the fit and simulation truth have no comparable parameters")
  summ <- do.call(rbind, lapply(pieces, function(d) {
    nm <- d$parameter[1]
    data.frame(
    parameter = nm, n = nrow(d),
    correlation = if (nrow(d) > 2) .safe_cor(d$true, d$estimated) else NA_real_,
    rmse = sqrt(mean((d$estimated - d$true)^2)),
    bias = if (isTRUE(centred[[nm]])) NA_real_ else mean(d$estimated - d$true),
    stringsAsFactors = FALSE)}))
  rownames(summ) <- NULL
  structure(list(summary = summ, pieces = pieces, layout = lay,
                 note = recovery_note),
            class = "rasch_recovery")
}

#' @export
print.rasch_recovery <- function(x, ...) {
  cat("Generating and fitted parameter comparison:\n")
  if (!is.null(x$note)) cat("  ", x$note, "\n", sep = "")
  s <- x$summary
  for (i in seq_len(nrow(s)))
    cat(sprintf("  %-16s n=%-4d r=%.3f  RMSE=%.3f  bias=%+.3f\n",
                s$parameter[i], s$n[i], s$correlation[i], s$rmse[i], s$bias[i]))
  invisible(x)
}

#' Plot fitted against generating parameters
#'
#' One true-versus-estimated panel per parameter type, with the identity line
#' and the correlation and RMSE.
#'
#' @param x A \code{"rasch_recovery"} object.
#' @param ... Unused.
#' @return Called for its plotting side effect.
#' @examples
#' \donttest{
#' d <- simulate_rasch(300, 8, seed = 1)
#' fit <- rasch(d, id = "id")
#' plot_recovery(sim_recovery(fit, d))
#' }
#' @export
plot_recovery <- function(x, ...) {
  if (!inherits(x, "rasch_recovery"))
    stop("`x` must be a recovery result from sim_recovery()", call. = FALSE)
  np <- length(x$pieces)
  op <- par(mfrow = c(1, np), mar = c(4.2, 4.2, 2.4, 1), mgp = c(2.4, 0.7, 0),
            las = 1, col.axis = .rr$ink, col.lab = .rr$ink, col.main = .rr$ink,
            cex.main = 1.0, font.main = 2)
  on.exit(par(op))
  for (i in seq_len(np)) {
    d <- x$pieces[[i]]; s <- x$summary[i, ]
    rng <- range(c(d$true, d$estimated))
    plot(d$true, d$estimated, xlim = rng, ylim = rng, xlab = "generating",
         ylab = "fitted", main = d$parameter[1], axes = FALSE,
         pch = 21, bg = .rr$blue, col = "white", cex = 1.1)
    abline(0, 1, col = .rr$red, lty = 2, lwd = 1.5)
    .rr_axis(1)
    .rr_axis(2)
    mtext(sprintf("r = %.3f   RMSE = %.2f", s$correlation, s$rmse), 3,
          line = 0.2, cex = 0.8, col = .rr$soft)
  }
}

#' Simulate paired-comparison EFRM data with differing frame units
#'
#' Generates dichotomous paired comparisons whose latent unit differs across
#' judge-panel by object-set frames -- the paired-comparison extension of the
#' extended frame of reference model (Humphry 2005) fitted by
#' \code{\link{btl_efrm}}. Objects in set \code{s} have a within-set
#' calibration location \code{beta}; their common-scale value is
#' \code{v = alpha_s beta + kappa_s}. A comparison judged in panel \code{g}
#' carries the panel unit \code{phi_g}: within a set the comparison logit is
#' \code{phi_g (beta_a - beta_b)}, across sets it is \code{phi_g (v_a - v_b)}.
#' The planted panel units, set units and origins are recovered by
#' \code{\link{btl_efrm}}.
#'
#' @param n_objects_per_set,n_sets Objects in each set and number of sets.
#' @param n_judges_per_panel,n_panels Judges in each panel and number of panels.
#'   Comparisons are balanced across panels and judges; a design with fewer
#'   comparisons than judges is refused.
#' @param reps_within Replications of each within-set object pair.
#' @param reps_cross Replications of each cross-set object pair.
#' @param panel_units Panel units \code{phi} (length \code{n_panels}); the
#'   default is all one, and any supplied vector is rescaled to geometric
#'   mean one.
#' @param set_units Set units \code{alpha} (length \code{n_sets}); the default
#'   is all one, and \code{alpha_1} is forced to one (the reference set).
#' @param set_origins Set origins \code{kappa} (length \code{n_sets}); the
#'   default is all zero, and \code{kappa_1} is forced to zero.
#' @param object_sd Realised sample standard deviation of the within-set
#'   calibration locations.
#' @param erratic_judges Proportion of judges who choose between the two
#'   objects at random. At least one judge in every panel must retain
#'   model-based comparisons.
#' @param seed Optional non-negative whole-number RNG seed. See
#'   \code{\link{rasch_rng}} for generator support.
#' @return A data frame of class \code{"rasch_sim"} with columns
#'   \code{object_a}, \code{object_b}, \code{winner}, \code{judge} and
#'   \code{panel}, and \code{attr(x, "truth")} holding the common-scale values
#'   \code{v}, the per-set \code{beta}, the units \code{phi}, \code{alpha},
#'   \code{kappa}, and the \code{object_sets} map to pass to
#'   \code{\link{btl_efrm}}.
#' @examples
#' d <- simulate_btl_efrm(6, 2, set_units = c(1, 1.4), seed = 1)
#' bt <- btl_efrm(d, "object_a", "object_b", winner = "winner",
#'                judge = "judge", panels = "panel",
#'                object_sets = attr(d, "truth")$object_sets,
#'                se_method = "conditional")
#' bt$alpha_table   # recovers the ~1.4 set unit
#' @export
simulate_btl_efrm <- function(n_objects_per_set = 8, n_sets = 2,
                              n_judges_per_panel = 6, n_panels = 2,
                              reps_within = 20, reps_cross = 20,
                              panel_units = NULL, set_units = NULL,
                              set_origins = NULL, object_sd = 1,
                              seed = NULL, erratic_judges = 0) {
  seed <- .sim_seed(seed)
  if (!is.null(seed)) {
    .old_stream <- .sim_seed_capture()
    on.exit(.sim_seed_restore(.old_stream), add = TRUE)
    set.seed(seed)
  }
  S <- .sim_count(n_sets, "n_sets")
  G <- .sim_count(n_panels, "n_panels")
  Kp <- .sim_count(n_objects_per_set, "n_objects_per_set", 2L)
  Jp <- .sim_count(n_judges_per_panel, "n_judges_per_panel")
  reps_within <- .sim_count(reps_within, "reps_within")
  reps_cross <- .sim_count(reps_cross, "reps_cross")
  object_sd <- .sim_scalar(object_sd, "object_sd", lower = 0,
                           lower_open = TRUE)
  erratic_judges <- .sim_scalar(erratic_judges, "erratic_judges",
                                lower = 0, upper = 1)

  phi <- if (is.null(panel_units)) rep(1, G) else panel_units
  if (!is.numeric(phi) || is.complex(phi) || !is.null(dim(phi)) ||
      !is.null(oldClass(phi)) ||
      length(phi) != G || any(!is.finite(phi) | phi <= 0))
    stop("panel_units must contain n_panels positive finite values")
  phi <- as.numeric(phi)
  alpha <- if (is.null(set_units)) rep(1, S) else set_units
  if (!is.numeric(alpha) || is.complex(alpha) || !is.null(dim(alpha)) ||
      !is.null(oldClass(alpha)) ||
      length(alpha) != S || any(!is.finite(alpha) | alpha <= 0))
    stop("set_units must contain n_sets positive finite values")
  alpha <- as.numeric(alpha)
  kappa <- if (is.null(set_origins)) rep(0, S) else set_origins
  if (!is.numeric(kappa) || is.complex(kappa) || !is.null(dim(kappa)) ||
      !is.null(oldClass(kappa)) ||
      length(kappa) != S || any(!is.finite(kappa)))
    stop("set_origins must contain n_sets finite values")
  kappa <- as.numeric(kappa)
  # a unit and an origin are relative to the other frames: with one frame
  # the normalisation below would silently return them to the reference
  # while the recorded truth claimed the requested values
  if (G == 1L && !isTRUE(all.equal(phi, 1)))
    stop("`panel_units` must be 1 when `n_panels` is 1: a panel unit is ",
         "relative to the other panels, so none can be planted")
  if (S == 1L && !isTRUE(all.equal(alpha, 1)))
    stop("`set_units` must be 1 when `n_sets` is 1: a set unit is relative ",
         "to the other sets, so none can be planted")
  if (S == 1L && !isTRUE(all.equal(kappa, 0)))
    stop("`set_origins` must be 0 when `n_sets` is 1: an origin is relative ",
         "to the other sets, so none can be planted")
  # Normalise on the log scale. Direct division can overflow even though all
  # supplied units are positive and finite (for example 1e-300 versus 1e300).
  lphi <- log(phi); phi <- exp(lphi - mean(lphi))       # geometric mean one
  lalpha <- log(alpha) - log(alpha[1])                  # alpha_1 = 1
  alpha <- exp(lalpha)
  kappa <- kappa - kappa[1]                            # kappa_1 = 0
  if (any(!is.finite(phi)) || any(phi <= 0))
    stop("the relative `panel_units` are outside the numerically representable range")
  if (any(!is.finite(alpha)) || any(alpha <= 0))
    stop("the relative `set_units` are outside the numerically representable range")
  if (any(!is.finite(kappa)))
    stop("the relative `set_origins` are outside the numerically representable range")

  set_nm <- sprintf("set%d", seq_len(S))
  panel_nm <- sprintf("panel%d", seq_len(G))
  objs_by_set <- lapply(seq_len(S), function(s) sprintf("S%dO%02d", s, seq_len(Kp)))
  names(objs_by_set) <- set_nm
  beta <- numeric(0)
  for (s in seq_len(S)) {
    bs <- as.numeric(scale(seq_len(Kp))) * object_sd   # sum-zero, spread object_sd
    .sim_check_sd(bs, object_sd, "within-set object-location standard deviation")
    beta <- c(beta, setNames(bs, objs_by_set[[s]]))
  }
  set_of <- setNames(rep(seq_len(S), each = Kp), unlist(objs_by_set))
  v <- alpha[set_of] * beta + kappa[set_of]
  if (any(!is.finite(v)))
    stop("the requested units, origins and object spread produce common-scale values outside the numerically representable range")

  judges <- sprintf("J%03d", seq_len(G * Jp))
  panel_of <- setNames(panel_nm[rep(seq_len(G), each = Jp)], judges)

  # assemble the object pairs (within each set, then across every set pair)
  aa <- bb <- character(0)
  for (s in seq_len(S)) {
    pr <- t(utils::combn(objs_by_set[[s]], 2))
    aa <- c(aa, rep(pr[, 1], reps_within)); bb <- c(bb, rep(pr[, 2], reps_within))
  }
  if (S > 1L) for (i in seq_len(S - 1L)) for (j in (i + 1L):S) {
    grid <- expand.grid(oa = objs_by_set[[i]], ob = objs_by_set[[j]],
                        KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
    aa <- c(aa, rep(grid$oa, reps_cross)); bb <- c(bb, rep(grid$ob, reps_cross))
  }
  R <- length(aa)
  # Balance panels within every set-pair stratum and then balance judges
  # within their panels. Requested panels and judges therefore cannot vanish
  # from a finite generated sample merely by chance.
  stratum <- paste(pmin(set_of[aa], set_of[bb]),
                   pmax(set_of[aa], set_of[bb]), sep = ":")
  if (R < length(judges))
    stop("the generated design has ", R, " comparisons, fewer than the ",
         length(judges), " requested judges")
  pan <- character(R)
  # Continue one randomly ordered panel cycle across strata. Counts differ by
  # at most one both overall and, when a stratum is large enough, within it.
  panel_cycle <- rep(sample(panel_nm, G, replace = FALSE), length.out = R)
  cursor <- 0L
  for (ss in unique(stratum)) {
    rows <- which(stratum == ss)
    take <- cursor + seq_along(rows)
    pan[rows] <- sample(panel_cycle[take], length(rows), replace = FALSE)
    cursor <- cursor + length(rows)
  }
  jd <- character(R)
  for (g in panel_nm) {
    rows <- which(pan == g)
    pool <- names(panel_of)[panel_of == g]
    jd[rows] <- .sim_balanced_ids(pool, length(rows),
                                  paste0("judges in ", g))
  }

  n_erratic <- .sim_planted_count(erratic_judges, length(judges))
  if (n_erratic >= length(judges))
    stop("erratic_judges must leave at least one judge with model-based comparisons")
  if (n_erratic > G * (Jp - 1L))
    stop("erratic_judges must leave at least one model-based judge in every panel")
  if (n_erratic) {
    # Distribute erratic judges across panels as evenly as their number
    # permits, then select them at random within each panel.
    take <- rep(n_erratic %/% G, G)
    extra <- n_erratic %% G
    if (extra) {
      extra_panels <- sample(seq_len(G), extra, replace = FALSE)
      take[extra_panels] <- take[extra_panels] + 1L
    }
    erratic <- unlist(lapply(seq_len(G), function(g) {
      pool <- names(panel_of)[panel_of == panel_nm[g]]
      if (take[g]) sample(pool, take[g], replace = FALSE) else character(0)
    }), use.names = FALSE)
  } else erratic <- character(0)
  # frame-dependent logit: within-set uses beta, cross-set uses the common v
  same <- set_of[aa] == set_of[bb]
  lp <- ifelse(same, phi[match(pan, panel_nm)] * (beta[aa] - beta[bb]),
               phi[match(pan, panel_nm)] * (v[aa] - v[bb]))
  if (any(!is.finite(lp)))
    stop("the requested units and object values produce logits outside the numerically representable range")
  win_a <- stats::runif(R) < stats::plogis(lp)
  bad <- jd %in% erratic
  if (any(bad)) win_a[bad] <- sample(c(FALSE, TRUE), sum(bad), replace = TRUE)
  d <- data.frame(object_a = aa, object_b = bb,
                  winner = ifelse(win_a, aa, bb),
                  judge = jd, panel = unname(pan),
                  stringsAsFactors = FALSE)
  rownames(d) <- NULL

  planted <- character(0)
  if (any(phi != 1)) planted <- c(planted, sprintf(
    "panel units phi = (%s)", paste(sprintf("%.2f", phi), collapse = ", ")))
  if (any(alpha != 1)) planted <- c(planted, sprintf(
    "set units alpha = (%s)", paste(sprintf("%.2f", alpha), collapse = ", ")))
  if (any(kappa != 0)) planted <- c(planted, sprintf(
    "set origins kappa = (%s)", paste(sprintf("%.2f", kappa), collapse = ", ")))
  if (length(erratic)) planted <- c(planted, sprintf(
    "%d erratic judge(s): %s", length(erratic), paste(erratic, collapse = ", ")))
  if (!length(planted)) planted <- "equal units (phi = alpha = 1)"

  attr(d, "truth") <- list(
    layout = "btl_efrm",
    description = sprintf("%d objects (%d sets) x %d panels, %d comparisons",
                          S * Kp, S, G, R),
    v = v, beta = beta, phi = setNames(phi, panel_nm),
    departure_types = if (length(erratic)) "erratic judges" else character(0),
    alpha = setNames(alpha, set_nm), kappa = setNames(kappa, set_nm),
    set_of = set_of, object_sets = objs_by_set, judge_panel = panel_of,
    erratic = erratic,
    planted = planted)
  class(d) <- c("rasch_sim", "data.frame")
  d
}
