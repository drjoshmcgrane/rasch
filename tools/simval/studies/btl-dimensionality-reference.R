#!/usr/bin/env Rscript
# Validate the finite-simulation upper-tail decision in btl_dimensionality().
# This version targets the pooled observed-minus-fitted expected POINT residual
# used for both the observed statistic and every null draw. Cheap ordinary cells
# cover binary, ordered-response, and a fitted first-position effect under a
# fitted null and a directed three-object cycle; BTL-EFRM covers linked sets.

suppressWarnings(pkgload::load_all(".", quiet = TRUE, compile = FALSE))
source("tools/simval/harness.R")

NULL_REPS <- as.integer(Sys.getenv("SV_NULL_REPS", Sys.getenv("SV_REPS", "1000")))
ALT_REPS <- as.integer(Sys.getenv("SV_ALT_REPS", "500"))
EFRM_REPS <- as.integer(Sys.getenv("SV_EFRM_REPS", "250"))
CORES <- as.integer(Sys.getenv("SV_CORES", "4"))
REFS <- as.numeric(strsplit(Sys.getenv("SV_INNER", "20,200"), ",", fixed = TRUE)[[1L]])
PAIR_N <- as.integer(Sys.getenv("SV_PAIR_N", "40"))
if (!is.finite(NULL_REPS) || NULL_REPS < 1000L)
  stop("SV_NULL_REPS (or SV_REPS) must be at least 1000")
if (!is.finite(ALT_REPS) || ALT_REPS < 2L) stop("SV_ALT_REPS must be at least 2")
if (!is.finite(EFRM_REPS) || EFRM_REPS < 2L) stop("SV_EFRM_REPS must be at least 2")
if (!is.finite(CORES) || CORES < 1L) CORES <- 1L
if (!is.finite(PAIR_N) || PAIR_N < 10L || PAIR_N != floor(PAIR_N))
  stop("SV_PAIR_N must be a whole number of at least 10")
if (!length(REFS) || any(!is.finite(REFS) | REFS < 20L | REFS != floor(REFS)))
  stop("SV_INNER must contain comma-separated whole numbers of at least 20")
REFS <- as.integer(REFS)

ordinary_models <- c("binary", "ordered response", "binary + position")
ordinary_grid <- expand.grid(
  model = ordinary_models, cycle_logit = c(0, 0.75, 1.5),
  reference_reps = REFS, stringsAsFactors = FALSE)
ordinary_grid$departure <- ifelse(
  ordinary_grid$cycle_logit == 0, "null", "directed cycle")
scenarios <- rbind(ordinary_grid,
  data.frame(model = "BTL-EFRM", cycle_logit = 0,
             reference_reps = REFS, departure = "null",
             stringsAsFactors = FALSE))
scenarios$n_attempted <- ifelse(
  scenarios$model == "BTL-EFRM", EFRM_REPS,
  ifelse(scenarios$departure == "null", NULL_REPS, ALT_REPS))

# Generate independent outcomes at every unordered pair, then compress to
# integer count rows. For position cells each pair is presented equally often
# in both directions and the true +1.0 first-position term is in the fitted null.
# Alternatives add an antisymmetric +0.75 or +1.50 logit to the directed
# O1>O2, O2>O3, O3>O1 cycle.
simulate_ordinary <- function(model, cycle_logit, seed) {
  set.seed(seed)
  objects <- sprintf("O%d", seq_len(8L))
  beta <- setNames(c(0, 0, 0, -1.2, -0.6, 0.5, 1.0, 1.5), objects)
  beta <- beta - mean(beta)
  pairs <- t(utils::combn(objects, 2L))
  position <- identical(model, "binary + position")
  m <- if (identical(model, "ordered response")) 4L else 1L
  tau <- if (m > 1L) c(-1.2, -0.4, 0.4, 1.2) else NULL
  n_each <- if (position) 2L * ceiling(PAIR_N / 2L) else PAIR_N
  out <- vector("list", nrow(pairs))
  for (k in seq_len(nrow(pairs))) {
    aa <- rep(pairs[k, 1L], n_each); bb <- rep(pairs[k, 2L], n_each)
    if (position) {
      flip <- rep(c(FALSE, TRUE), length.out = n_each)
      a0 <- aa; aa[flip] <- bb[flip]; bb[flip] <- a0[flip]
    }
    cycle <- numeric(n_each)
    if (cycle_logit > 0) {
      forward <- paste(aa, bb, sep = ">") %in% c("O1>O2", "O2>O3", "O3>O1")
      reverse <- paste(bb, aa, sep = ">") %in% c("O1>O2", "O2>O3", "O3>O1")
      cycle <- cycle_logit * (forward - reverse)
    }
    eta <- beta[aa] - beta[bb] + cycle + if (position) 1 else 0
    response <- if (m == 1L) {
      stats::rbinom(n_each, 1L, stats::plogis(eta))
    } else {
      vapply(eta, function(e)
        sample.int(m + 1L, 1L, prob = item_moments(e, tau)$P) - 1L,
        integer(1L))
    }
    out[[k]] <- data.frame(object_a = aa, object_b = bb, response = response,
                           stringsAsFactors = FALSE)
  }
  d <- do.call(rbind, out); d$count <- 1L
  stats::aggregate(count ~ object_a + object_b + response, d, sum)
}

fit_one <- function(model, departure, cycle_logit, seed) {
  if (model == "BTL-EFRM") {
    d <- simulate_btl_efrm(n_objects_per_set = 6, n_sets = 2, n_panels = 2,
      n_judges_per_panel = 10, reps_within = 20, reps_cross = 20, seed = seed)
    btl_efrm(d, "object_a", "object_b", winner = "winner", judge = "judge",
      panels = "panel", object_sets = attr(d, "truth")$object_sets,
      se_method = "conditional")
  } else {
    d <- simulate_ordinary(model, cycle_logit, seed)
    btl(d, "object_a", "object_b", response = "response", count = "count",
        position = identical(model, "binary + position"))
  }
}

one <- function(r, model, departure, cycle_logit, reference_reps, tag) {
  # Data are paired across reference sizes; only the null-reference stream varies.
  model_tag <- match(model, c(ordinary_models, "BTL-EFRM"))
  dep_tag <- match(departure, c("null", "directed cycle"))
  data_seed <- 930000L + model_tag * 100000L + dep_tag * 10000L + r
  fit <- tryCatch(fit_one(model, departure, cycle_logit, data_seed), error = identity)
  if (inherits(fit, "condition"))
    return(c(reject = NA, p = NA, refused = 1, nonconv = 0, error = 0, withheld = 0))
  if (!isTRUE(fit$converged))
    return(c(reject = NA, p = NA, refused = 0, nonconv = 1, error = 0, withheld = 0))
  ans <- tryCatch(btl_dimensionality(fit, reps = reference_reps,
    seed = 1430000L + tag * 10000L + r,
    # EFRM comparisons share judges. This deliberately requests the
    # model-conditional reference, not unconditional evidence of independence.
    independent_comparisons = TRUE), error = identity)
  if (inherits(ans, "condition"))
    return(c(reject = NA, p = NA, refused = 0, nonconv = 0, error = 1, withheld = 0))
  if (is.na(ans$leading_structured))
    return(c(reject = NA, p = NA, refused = 0, nonconv = 0, error = 0, withheld = 1))
  if (!is.finite(ans$reference$p_adj))
    return(c(reject = NA, p = NA, refused = 0, nonconv = 0, error = 1, withheld = 0))
  c(reject = as.numeric(ans$leading_structured), p = ans$reference$p_adj,
    refused = 0, nonconv = 0, error = 0, withheld = 0)
}

run_cell <- function(model, departure, cycle_logit, reference_reps, tag, n_attempted) {
  f <- function(r) one(r, model, departure, cycle_logit, reference_reps, tag)
  z <- if (CORES > 1L && .Platform$OS.type != "windows")
    parallel::mclapply(seq_len(n_attempted), f, mc.cores = CORES,
                       mc.set.seed = FALSE) else lapply(seq_len(n_attempted), f)
  do.call(rbind, z)
}

rows <- vector("list", nrow(scenarios))
for (s in seq_len(nrow(scenarios))) {
  model <- scenarios$model[s]; departure <- scenarios$departure[s]
  cycle_logit <- scenarios$cycle_logit[s]
  reference_reps <- scenarios$reference_reps[s]
  n_attempted <- scenarios$n_attempted[s]
  cat(sprintf("%s, %s%s, B=%d: %d attempts\n", model, departure,
              if (cycle_logit > 0) sprintf(" %.2f", cycle_logit) else "",
              reference_reps, n_attempted))
  z <- run_cell(model, departure, cycle_logit, reference_reps, s, n_attempted)
  ok <- is.finite(z[, "reject"]) & is.finite(z[, "p"])
  rate <- if (any(ok)) mean(z[ok, "reject"]) else NA_real_
  mean_p <- if (any(ok)) mean(z[ok, "p"]) else NA_real_
  design_note <- if (model == "BTL-EFRM") {
    "6 objects/set, 2 sets, 2 panels, conditional link; model-conditional independence opt-in"
  } else if (model == "binary + position") {
    sprintf("8 objects, %d outcomes/pair/direction, fitted true position +1.0; integer count rows",
            ceiling(PAIR_N / 2L))
  } else sprintf("8 objects, %d outcomes/pair; integer count rows", PAIR_N)
  rows[[s]] <- sv_row(
    "BTL dimensionality pooled-expected reference",
    sprintf("%s, %s%s, B = %d", model, departure,
      if (cycle_logit > 0) sprintf(" %.2f", cycle_logit) else "",
      reference_reps),
    "finite simulated upper-tail decision", n_reps = sum(ok),
    n_attempted = n_attempted, n_refused = sum(z[, "refused"]),
    n_nonconv = sum(z[, "nonconv"]), n_error = sum(z[, "error"]),
    n_withheld = sum(z[, "withheld"]),
    effect = if (departure == "directed cycle") cycle_logit else NA_real_,
    type1 = if (departure == "null") rate else NA_real_,
    power = if (departure != "null") rate else NA_real_,
    notes = sprintf(paste0("pooled observed-minus-fitted expected point residual ",
      "with symmetric 0.5 continuity; rejection %.4f, mean p %.4f using ",
      "(1 + exceedances) / (B + 1); %s%s"), rate, mean_p, design_note,
      if (departure == "directed cycle") sprintf(
        "; cycle logit +%.2f on O1>O2>O3>O1", cycle_logit) else ""))
}

out <- do.call(rbind, rows)
sv_write(out, "btl-dimensionality-reference")
print(out[, c("scenario", "n_reps", "n_refused", "n_nonconv", "n_error",
  "n_withheld", "type1", "mc_se_type1", "power", "mc_se_power")],
  row.names = FALSE)
