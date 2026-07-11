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
#        guessing, groups, planted = <character, human-readable pathologies>)

# null-coalescing helper (package-internal; base R gained %||% only in 4.4)
`%||%` <- function(a, b) if (is.null(a)) b else a

# person locations from one of a few distributions
.sim_theta <- function(n, mean, sd, dist = "normal") {
  z <- switch(dist,
    normal  = stats::rnorm(n),
    uniform = stats::runif(n, -sqrt(3), sqrt(3)),
    skew    = { u <- stats::rgamma(n, 2, 1); (u - 2) / sqrt(2) },
    bimodal = { s <- sample(c(-1, 1), n, TRUE); s * 1.1 + stats::rnorm(n, 0, 0.5) },
    stats::rnorm(n))
  mean + sd * as.numeric(scale(z))
}

# draw one item's responses for a vector of person locations. tau are the
# item's thresholds (length m; dichotomous m = 1). disc scales the whole
# exponent (a departure when != 1); guess is a lower asymptote (dichotomous).
.sim_item <- function(theta, tau, disc = 1, guess = 0) {
  m <- length(tau); xs <- 0:m
  cum <- c(0, cumsum(tau))
  eta <- disc * (outer(theta, xs) - matrix(cum, length(theta), m + 1L, byrow = TRUE))
  eta <- eta - apply(eta, 1, max)
  P <- exp(eta); P <- P / rowSums(P)
  if (guess > 0 && m == 1L) {
    P[, 2] <- guess + (1 - guess) * P[, 2]; P[, 1] <- 1 - P[, 2]
  }
  cs <- t(apply(P, 1, cumsum))
  as.integer(rowSums(stats::runif(length(theta)) > cs))     # category 0..m
}

# even thresholds around an item location with a given spread (ordered), or
# deliberately disordered (a middle threshold moved below the one before it)
.sim_thresholds <- function(delta, m, spread, disordered = FALSE) {
  if (m == 1L) return(delta)
  step <- seq(-1, 1, length.out = m) * spread
  tau <- delta + step - mean(step)
  if (disordered && m >= 2L) { i <- ceiling(m / 2); tau[i] <- tau[i] - 2.2 * spread }
  tau
}

#' Simulate person-by-item Rasch data with dial-in misfit
#'
#' Generates dichotomous or polytomous (partial credit / rating scale) data
#' from the Rasch model, with optional, individually controllable departures
#' from it -- each of which the package's matching diagnostic is built to
#' detect. The result is a data frame ready for \code{\link{rasch}}, with the
#' true parameters attached as \code{attr(x, "truth")}.
#'
#' @param n_persons,n_items Sample size and test length.
#' @param model \code{"dichotomous"}, \code{"PCM"}, or \code{"RSM"}.
#' @param n_categories Response categories for polytomous models (>= 3).
#' @param theta_mean,theta_sd,theta_dist Person distribution: mean, SD, and
#'   shape (\code{"normal"}, \code{"uniform"}, \code{"skew"}, \code{"bimodal"}).
#' @param difficulty Two numbers giving the item-location range (evenly
#'   spaced), or a length-\code{n_items} vector of locations.
#' @param threshold_spread Half-range of the category thresholds about each
#'   item location (polytomous).
#' @param discrimination Scalar or length-\code{n_items}: the slope of each
#'   item. Values above 1 over-discriminate (Guttman-like, negative fit
#'   residual); below 1 under-discriminate (noisy, positive residual). Feeds
#'   infit/outfit and the item-fit F.
#' @param guessing Scalar or length-\code{n_items} lower asymptote
#'   (dichotomous): low-ability persons answer correctly by chance. Feeds
#'   \code{\link{tailored_analysis}}.
#' @param second_dim \code{NULL}, or \code{list(items=, rho=)}: the named items
#'   load on a second trait correlated \code{rho} with the first. Feeds
#'   \code{\link{dimensionality_test}}.
#' @param dependence \code{NULL}, or \code{list(pairs=, strength=)}: each pair's
#'   second item responds partly to the first (response dependence). Feeds
#'   \code{\link{residual_correlations}} / \code{\link{dependence_magnitude}}.
#' @param dif \code{NULL}, or \code{list(items=, uniform=, nonuniform=)}: the
#'   named items function differently for the last person group -- a location
#'   shift (\code{uniform}) and/or a slope change (\code{nonuniform}). Needs
#'   \code{n_groups >= 2}. Feeds \code{\link{dif_anova}} / \code{\link{dif_size}}.
#' @param careless Proportion of persons who answer at random (person misfit;
#'   feeds person infit/outfit).
#' @param disordered \code{NULL} or item names/indices given disordered
#'   thresholds (polytomous; feeds the threshold diagnostics).
#' @param n_groups Number of equal person groups (a \code{group} factor column
#'   is added when > 1, for DIF).
#' @param missing Proportion of responses set missing (completely at random).
#' @param seed Optional RNG seed.
#' @return A data frame of class \code{"rasch_sim"} (item columns
#'   \code{I01}..., an \code{id} column, and a \code{group} column when
#'   grouped), with \code{attr(x, "truth")} holding the generating parameters
#'   and the planted departures.
#' @examples
#' # a clean scale with one over-discriminating item and one DIF item
#' d <- simulate_rasch(400, 12, discrimination = c(3, rep(1, 11)),
#'                     dif = list(items = "I06", uniform = 1), n_groups = 2,
#'                     seed = 1)
#' fit <- rasch(d, factors = "group")
#' fit$items[c("item", "infit_ms", "outfit_ms")]   # item 1 misfits
#' dif_anova(fit)$summary                           # item 6 flags
#' @export
simulate_rasch <- function(n_persons = 500, n_items = 20,
                           model = c("dichotomous", "PCM", "RSM"),
                           n_categories = 3, theta_mean = 0, theta_sd = 1,
                           theta_dist = "normal", difficulty = c(-2.5, 2.5),
                           threshold_spread = 1.2, discrimination = 1,
                           guessing = 0, second_dim = NULL, dependence = NULL,
                           dif = NULL, careless = 0, disordered = NULL,
                           n_groups = 1, missing = 0, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  model <- match.arg(model)
  m <- if (model == "dichotomous") 1L else as.integer(n_categories) - 1L
  I <- as.integer(n_items); N <- as.integer(n_persons)
  inm <- sprintf("I%02d", seq_len(I))
  as_idx <- function(x) if (is.character(x)) match(x, inm) else as.integer(x)

  # item locations, thresholds, slopes, guessing (with per-item overrides)
  delta <- if (length(difficulty) == I) difficulty
           else seq(difficulty[1], difficulty[2], length.out = I)
  disc <- if (length(discrimination) == I) discrimination else rep(discrimination[1], I)
  guess <- if (length(guessing) == I) guessing else rep(guessing[1], I)
  dis_items <- as_idx(disordered)
  tau <- lapply(seq_len(I), function(i)
    .sim_thresholds(delta[i], m, threshold_spread, i %in% dis_items))

  # person locations (primary) and groups
  theta <- .sim_theta(N, theta_mean, theta_sd, theta_dist)
  group <- if (n_groups > 1L)
    factor(sprintf("g%d", (seq_len(N) - 1L) %% n_groups + 1L)) else NULL

  # a second dimension for the nominated items: a correlated latent trait
  theta2 <- NULL; dim_items <- integer(0)
  if (!is.null(second_dim)) {
    dim_items <- as_idx(second_dim$items)
    rho <- second_dim$rho %||% 0.5
    theta2 <- rho * theta + sqrt(1 - rho^2) *
      .sim_theta(N, theta_mean, theta_sd, theta_dist)
  }

  X <- matrix(NA_integer_, N, I, dimnames = list(NULL, inm))
  dif_items <- as_idx(if (is.null(dif)) NULL else dif$items)
  dif_grp <- if (n_groups > 1L) levels(group)[n_groups] else NA

  for (i in seq_len(I)) {
    th <- if (i %in% dim_items) theta2 else theta
    if (i %in% dif_items && !is.na(dif_grp)) {
      # differential functioning for the last group: location shift and/or
      # slope change on this item
      g2 <- group == dif_grp
      X[!g2, i] <- .sim_item(th[!g2], tau[[i]], disc[i], guess[i])
      X[g2, i]  <- .sim_item(th[g2], tau[[i]] + (dif$uniform %||% 0),
                             disc[i] + (dif$nonuniform %||% 0), guess[i])
    } else {
      X[, i] <- .sim_item(th, tau[[i]], disc[i], guess[i])
    }
  }

  # response dependence: the second item of each pair partly follows the
  # first (adds d*(x1 - E1) to its exponent, inducing residual correlation)
  dep_pairs <- list()
  if (!is.null(dependence)) {
    d_str <- dependence$strength %||% 1
    for (pp in dependence$pairs) {
      ij <- as_idx(pp); i1 <- ij[1]; i2 <- ij[2]
      e1 <- vapply(theta, function(t)
        sum((0:m) * .p_item(t, tau[[i1]], disc[i1])), 0)
      shift <- d_str * (X[, i1] - e1) / m           # per-person carry-over
      X[, i2] <- .sim_item(theta + shift, tau[[i2]], disc[i2], guess[i2])
      dep_pairs[[length(dep_pairs) + 1L]] <- inm[ij]
    }
  }

  # careless responders: answer uniformly at random
  careless_idx <- integer(0)
  if (careless > 0) {
    careless_idx <- sample(N, round(careless * N))
    X[careless_idx, ] <- matrix(sample(0:m, length(careless_idx) * I, TRUE),
                                length(careless_idx), I)
  }
  if (missing > 0) X[sample(length(X), round(missing * length(X)))] <- NA

  out <- data.frame(id = sprintf("P%04d", seq_len(N)), X,
                    check.names = FALSE, stringsAsFactors = FALSE)
  if (!is.null(group)) out$group <- group

  planted <- character(0)
  if (any(disc != 1)) planted <- c(planted, sprintf("discrimination != 1 on %s",
    paste(inm[disc != 1], collapse = ", ")))
  if (any(guess > 0)) planted <- c(planted, sprintf("guessing on %s",
    paste(inm[guess > 0], collapse = ", ")))
  if (length(dim_items)) planted <- c(planted, sprintf(
    "second dimension (rho %.2f) on %s", second_dim$rho %||% 0.5,
    paste(inm[dim_items], collapse = ", ")))
  if (length(dep_pairs)) planted <- c(planted, sprintf(
    "response dependence: %s", paste(vapply(dep_pairs, paste, "",
                                            collapse = "-"), collapse = "; ")))
  if (length(dif_items)) planted <- c(planted, sprintf(
    "DIF (group %s) on %s: uniform %.2f, non-uniform %.2f", dif_grp,
    paste(inm[dif_items], collapse = ", "), dif$uniform %||% 0,
    dif$nonuniform %||% 0))
  if (length(careless_idx)) planted <- c(planted, sprintf(
    "%d careless responder(s)", length(careless_idx)))
  if (length(dis_items)) planted <- c(planted, sprintf(
    "disordered thresholds on %s", paste(inm[dis_items], collapse = ", ")))
  if (missing > 0) planted <- c(planted, sprintf("%.0f%% missing", 100 * missing))

  attr(out, "truth") <- list(
    layout = "rasch",
    description = sprintf("%s, %d persons x %d items%s", model, N, I,
      if (!is.null(group)) sprintf(", %d groups", nlevels(group)) else ""),
    model = model, n_persons = N, n_items = I,
    theta = theta, difficulty = delta, thresholds = tau,
    discrimination = disc, guessing = guess,
    groups = group, dim_items = inm[dim_items], dif_items = inm[dif_items],
    careless_idx = careless_idx, planted = planted)
  class(out) <- c("rasch_sim", "data.frame")
  out
}

# category probabilities for one location (used by the dependence term)
.p_item <- function(theta, tau, disc = 1) {
  m <- length(tau); cum <- c(0, cumsum(tau))
  e <- disc * ((0:m) * theta - cum); e <- e - max(e)
  p <- exp(e); p / sum(p)
}

#' @export
print.rasch_sim <- function(x, ...) {
  tr <- attr(x, "truth")
  cat(sprintf("Simulated %s data: %s\n", tr$layout,
              tr$description %||% sprintf("%d rows", nrow(x))))
  if (length(tr$planted)) {
    cat("Planted departures:\n")
    for (p in tr$planted) cat(paste0("  - ", p, "\n"))
  } else cat("Model-conforming (no departures planted).\n")
  invisible(x)
}

#' Simulate paired-comparison (BTL) data with dial-in misfit
#'
#' Generates dichotomous or graded paired comparisons from the
#' Bradley-Terry-Luce model, with optional departures each of which a
#' paired-comparison diagnostic is built to detect. The result is a data frame
#' ready for \code{\link{btl}}, with the truth attached.
#'
#' @param n_objects,n_judges Objects to scale and judges comparing them.
#' @param reps_per_pair Comparisons made of each object pair.
#' @param model \code{"dichotomous"} (a winner) or \code{"graded"} (a rated
#'   margin in \code{n_categories} categories).
#' @param n_categories Categories for the graded model.
#' @param object_sd Spread of the object locations (evenly spaced, sum-zero).
#' @param second_attribute \code{NULL}, or \code{list(rho=)}: half the judges
#'   rank by a second object attribute correlated \code{rho} with the first --
#'   genuine multidimensionality. Feeds \code{\link{btl_dimensionality}} and
#'   \code{\link{btl_transitivity}}.
#' @param erratic_judges Proportion of judges who choose at random. Feeds the
#'   judge fit residual, \code{\link{btl_transitivity}} consistency, and
#'   \code{\link{judge_surprise}}.
#' @param dependence \code{NULL}, or \code{list(exposure=, carry_over=)}:
#'   within-judge order effects (a seen-before advantage and a pull from the
#'   judge's own earlier verdicts). Adds an \code{order} column. Feeds the
#'   dependence effects of \code{\link{btl}}.
#' @param seed Optional RNG seed.
#' @return A data frame of class \code{"rasch_sim"}: \code{object_a},
#'   \code{object_b}, \code{winner} (or \code{response} when graded),
#'   \code{judge}, and \code{order} when dependence is planted; with
#'   \code{attr(x, "truth")}.
#' @examples
#' d <- simulate_btl(8, 12, erratic_judges = 0.15, seed = 1)
#' bt <- btl(d, "object_a", "object_b", winner = "winner", judge = "judge")
#' bt$judges          # the erratic judges carry large fit residuals
#' @export
simulate_btl <- function(n_objects = 8, n_judges = 12, reps_per_pair = 25,
                         model = c("dichotomous", "graded"), n_categories = 4,
                         object_sd = 1, second_attribute = NULL,
                         erratic_judges = 0, dependence = NULL, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  model <- match.arg(model)
  m <- if (model == "graded") as.integer(n_categories) - 1L else 1L
  K <- as.integer(n_objects); J <- as.integer(n_judges)
  objs <- sprintf("O%d", seq_len(K)); jids <- sprintf("J%d", seq_len(J))
  beta <- setNames(as.numeric(scale(seq_len(K))) * object_sd, objs)
  tau <- if (m > 1L) .sim_thresholds(0, m, 1.2) else NULL

  # a second object attribute (orthogonal part) for the two-camp design
  beta2 <- NULL; camp <- NULL
  if (!is.null(second_attribute)) {
    rho <- second_attribute$rho %||% 0.3
    beta2 <- setNames(rho * beta + sqrt(1 - rho^2) *
      as.numeric(scale(stats::rnorm(K))) * object_sd, objs)
    camp <- setNames(rep(c("a", "b"), length.out = J), jids)
  }
  erratic <- if (erratic_judges > 0)
    jids[seq_len(round(erratic_judges * J))] else character(0)

  pr <- t(utils::combn(objs, 2))
  d <- data.frame(object_a = rep(pr[, 1], each = reps_per_pair),
                  object_b = rep(pr[, 2], each = reps_per_pair),
                  stringsAsFactors = FALSE)
  d$judge <- sample(jids, nrow(d), TRUE)

  win_prob <- function(a, b, jd) {
    ba <- if (!is.null(camp) && camp[jd] == "b") beta2[a] else beta[a]
    bb <- if (!is.null(camp) && camp[jd] == "b") beta2[b] else beta[b]
    ba - bb
  }
  # dependence needs a per-judge judgment order and running history
  ord <- NULL
  if (!is.null(dependence)) {
    d <- d[order(d$judge), ]
    d$order <- stats::ave(seq_len(nrow(d)), d$judge, FUN = seq_along)
    seen <- new.env(parent = emptyenv()); hs <- new.env(parent = emptyenv())
    hc <- new.env(parent = emptyenv())
    g0 <- function(e, k) if (is.null(v <- e[[k]])) 0 else v
    exq <- dependence$exposure %||% 0; cry <- dependence$carry_over %||% 0
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
  if (!is.null(second_attribute)) planted <- c(planted,
    sprintf("second object attribute (rho %.2f), two judge camps",
            second_attribute$rho %||% 0.3))
  if (!is.null(dependence)) planted <- c(planted, sprintf(
    "within-judge dependence: exposure %.2f, carry-over %.2f",
    dependence$exposure %||% 0, dependence$carry_over %||% 0))

  attr(d, "truth") <- list(
    layout = "btl",
    description = sprintf("%s, %d objects, %d judges, %d comparisons",
                          model, K, J, nrow(d)),
    model = model, location = beta, location2 = beta2, camp = camp,
    erratic = erratic, planted = planted)
  class(d) <- c("rasch_sim", "data.frame")
  d
}

#' Simulate many-facet (rated) data with dial-in misfit
#'
#' Generates ratings from the many-facet Rasch model (Linacre 1989): every
#' rater rates every person on every item, from person ability, item
#' difficulty, and rater severity. Departures each feed an MFRM diagnostic.
#'
#' @param n_persons,n_items,n_raters Facet sizes (fully crossed).
#' @param n_categories Rating categories.
#' @param theta_sd,item_sd Spread of person ability and item difficulty.
#' @param rater_severity_sd Spread of rater severities (the core facet;
#'   recovered in \code{facet_effects}).
#' @param erratic_raters Proportion of raters who rate at random (feeds the
#'   rater fit residual).
#' @param interaction \code{NULL}, or \code{list(rater=, item=, bias=)}: one
#'   rater is unusually harsh (positive) or lenient (negative) on one item.
#'   Feeds the item-by-rater interaction (fit with \code{interaction = }).
#' @param seed Optional RNG seed.
#' @return A long data frame of class \code{"rasch_sim"} (\code{person},
#'   \code{item}, \code{rater}, \code{score}) ready for
#'   \code{\link{rasch_mfrm}}, with the truth attached.
#' @examples
#' d <- simulate_mfrm(60, 5, 6, rater_severity_sd = 0.8, seed = 1)
#' mf <- rasch_mfrm(d, person = "person", item = "item", score = "score",
#'                  facets = "rater")
#' cor(mf$facet_effects$rater$measure, attr(d, "truth")$severity)  # recovered
#' @export
simulate_mfrm <- function(n_persons = 80, n_items = 5, n_raters = 6,
                          n_categories = 4, theta_sd = 1.2, item_sd = 1,
                          rater_severity_sd = 0.6, erratic_raters = 0,
                          interaction = NULL, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  m <- as.integer(n_categories) - 1L
  N <- as.integer(n_persons); I <- as.integer(n_items); R <- as.integer(n_raters)
  pids <- sprintf("P%03d", seq_len(N)); iids <- sprintf("I%d", seq_len(I))
  rids <- sprintf("R%d", seq_len(R))
  theta <- .sim_theta(N, 0, theta_sd)
  delta <- setNames(seq(-item_sd, item_sd, length.out = I), iids)
  lambda <- setNames(as.numeric(scale(stats::rnorm(R))) * rater_severity_sd, rids)
  base_tau <- .sim_thresholds(0, m, 1.2)
  erratic <- if (erratic_raters > 0) rids[seq_len(round(erratic_raters * R))] else character(0)
  int_bias <- matrix(0, I, R, dimnames = list(iids, rids))
  if (!is.null(interaction))
    int_bias[interaction$item, interaction$rater] <- interaction$bias

  grid <- expand.grid(p = seq_len(N), i = seq_len(I), r = seq_len(R))
  score <- integer(nrow(grid))
  for (i in seq_len(I)) for (r in seq_len(R)) {
    rows <- grid$i == i & grid$r == r
    # item difficulty and rater severity shift the person's thresholds
    tau_ir <- base_tau + delta[i] + lambda[r] + int_bias[i, r]
    score[rows] <- if (rids[r] %in% erratic) sample(0:m, sum(rows), TRUE)
                   else .sim_item(theta[grid$p[rows]], tau_ir)
  }
  d <- data.frame(person = pids[grid$p], item = iids[grid$i],
                  rater = rids[grid$r], score = score, stringsAsFactors = FALSE)

  planted <- character(0)
  if (length(erratic)) planted <- c(planted,
    sprintf("%d erratic rater(s): %s", length(erratic), paste(erratic, collapse = ", ")))
  if (!is.null(interaction)) planted <- c(planted, sprintf(
    "rater-by-item bias: %s on %s (%.2f)", interaction$rater,
    interaction$item, interaction$bias))
  planted <- c(planted, sprintf("rater severities SD %.2f", stats::sd(lambda)))

  attr(d, "truth") <- list(
    layout = "mfrm",
    description = sprintf("%d persons x %d items x %d raters (%d ratings)",
                          N, I, R, nrow(d)),
    theta = theta, difficulty = delta, severity = lambda,
    erratic = erratic, planted = planted)
  class(d) <- c("rasch_sim", "data.frame")
  d
}

#' Simulate extended frame-of-reference data with differing units
#'
#' Generates data whose latent unit differs across item-set by person-group
#' frames (Humphry 2005): a person in group g responding to an item in set s
#' does so at the frame unit rho = alpha_set * phi_group scaling the whole
#' exponent. The planted set- and group-unit ratios are recovered by
#' \code{\link{rasch_efrm}}.
#'
#' @param n_per_group Persons in each group.
#' @param items_per_set Items in each set.
#' @param n_sets,n_groups Numbers of item sets and person groups.
#' @param set_unit_ratio,group_unit_ratio Geometric span of the set and group
#'   units across their levels (1 = equal units, i.e. an ordinary Rasch fit).
#' @param theta_sd Spread of person ability.
#' @param seed Optional RNG seed.
#' @return A wide data frame of class \code{"rasch_sim"} (\code{id}, item
#'   columns, \code{group}) with \code{attr(x, "truth")$item_sets} the set map
#'   to pass to \code{\link{rasch_efrm}}.
#' @examples
#' d <- simulate_efrm(300, 8, set_unit_ratio = 1.3, seed = 1)
#' tr <- attr(d, "truth")
#' ef <- rasch_efrm(d, item_sets = tr$item_sets, groups = "group")
#' ef$alpha_table   # recovers the ~1.3 set-unit ratio
#' @export
simulate_efrm <- function(n_per_group = 300, items_per_set = 8, n_sets = 2,
                          n_groups = 2, set_unit_ratio = 1.3,
                          group_unit_ratio = 1, theta_sd = 1.3, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  S <- as.integer(n_sets); G <- as.integer(n_groups); K <- as.integer(items_per_set)
  npg <- as.integer(n_per_group)
  # set and group units span the ratio geometrically, normalised to mean 1
  gspan <- function(ratio, n) { u <- exp(seq(0, log(ratio), length.out = n)); u / exp(mean(log(u))) }
  alpha <- gspan(set_unit_ratio, S)
  phi <- gspan(group_unit_ratio, G)
  set_items <- lapply(seq_len(S), function(s) sprintf("S%dI%02d", s, seq_len(K)))
  inm <- unlist(set_items)
  delta <- setNames(rep(seq(-1.5, 1.5, length.out = K), S), inm)
  set_of <- setNames(rep(seq_len(S), each = K), inm)

  grp <- factor(rep(sprintf("g%d", seq_len(G)), each = npg))
  N <- length(grp); theta <- .sim_theta(N, 0, theta_sd)
  X <- matrix(NA_integer_, N, length(inm), dimnames = list(NULL, inm))
  for (col in seq_along(inm)) {
    s <- set_of[inm[col]]; rho <- alpha[s] * phi[as.integer(grp)]  # per-person unit
    X[, col] <- as.integer(stats::runif(N) < stats::plogis(rho * (theta - delta[inm[col]])))
  }
  out <- data.frame(id = sprintf("P%04d", seq_len(N)), X, group = grp,
                    check.names = FALSE, stringsAsFactors = FALSE)

  planted <- sprintf("set-unit ratio %.2f across %d sets", set_unit_ratio, S)
  if (group_unit_ratio != 1)
    planted <- c(planted, sprintf("group-unit ratio %.2f across %d groups",
                                  group_unit_ratio, G))
  attr(out, "truth") <- list(
    layout = "efrm",
    description = sprintf("%d persons, %d sets x %d groups, %d items",
                          N, S, G, length(inm)),
    theta = theta, difficulty = delta, alpha = alpha, phi = phi,
    item_sets = setNames(set_items, sprintf("set%d", seq_len(S))),
    groups = grp, planted = planted)
  class(out) <- c("rasch_sim", "data.frame")
  out
}
