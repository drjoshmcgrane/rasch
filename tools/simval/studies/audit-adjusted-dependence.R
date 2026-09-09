# STUDY: audit-adjusted-dependence
#
# Validation for the crossed-EFRM and comparative-judgement multiplicity
# decisions introduced in the 1.12.1 audit, together with the corrected
# finite-object correlation used by simulate_btl(). Principal null claims use
# 1,000 replicates; power is examined at two effect sizes. Set SV_CORES on a
# Unix-like system to run independent replicates in parallel.

suppressWarnings(pkgload::load_all(".", quiet = TRUE, compile = FALSE))
source("tools/simval/harness.R")
STUDY <- "audit-adjusted-dependence"
ALPHA <- 0.05
DIM_ONLY <- identical(Sys.getenv("SV_DIM_ONLY", "0"), "1")

env_count <- function(name, default) {
  z <- suppressWarnings(as.integer(Sys.getenv(name, as.character(default))))
  if (!is.finite(z) || z < 1L) stop(name, " must be a positive integer")
  z
}
R_E_NULL <- env_count("SV_EFRM_NULL", 1000L)
R_E_POWER <- env_count("SV_EFRM_POWER", 400L)
R_B_NULL <- env_count("SV_BTL_NULL", 1000L)
R_B_POWER <- env_count("SV_BTL_POWER", 400L)
R_DIM <- env_count("SV_BTL_DIM", 100L)
CORES <- env_count("SV_CORES", 1L)

run_reps <- function(n, fun) {
  z <- if (CORES > 1L && .Platform$OS.type != "windows")
    parallel::mclapply(seq_len(n), fun, mc.cores = CORES,
                       mc.set.seed = FALSE)
  else lapply(seq_len(n), fun)
  bad <- vapply(z, inherits, logical(1), "try-error")
  if (any(bad)) stop("worker failure: ", as.character(z[[which(bad)[1L]]]))
  z
}

check_holm <- function(tab) {
  ok <- is.finite(tab$p)
  expected <- rep(NA_real_, nrow(tab))
  expected[ok] <- stats::p.adjust(tab$p[ok], method = "holm")
  isTRUE(all.equal(tab$p_adj, expected, tolerance = 1e-14,
                   check.attributes = FALSE))
}

# Crossed EFRM: one item set isolates the factorial group-unit inference from
# set linking. Four balanced cells define two main effects and their
# interaction. The planted region log-unit contrast is log(1.25) or log(1.50).
gen_crossed <- function(seed, region_log_ratio = 0) {
  set.seed(seed)
  N <- 400L
  region <- factor(rep(c("North", "South"), each = N / 2L))
  cohort <- factor(rep(rep(c("A", "B"), each = N / 4L), 2L))
  theta <- stats::rnorm(N, sd = 1.3)
  delta <- seq(-1.8, 1.8, length.out = 12L)
  log_unit <- ifelse(region == "South", region_log_ratio / 2,
                     -region_log_ratio / 2)
  X <- vapply(delta, function(b)
    stats::rbinom(N, 1, stats::plogis(exp(log_unit) * (theta - b))),
    integer(N))
  colnames(X) <- sprintf("I%02d", seq_len(ncol(X)))
  data.frame(X, region, cohort, check.names = FALSE)
}

one_efrm <- function(r, effect) {
  blank <- c(raw_any = NA, adjusted_any = NA, target_adjusted = NA,
             refused = 0, nonconv = 0)
  d <- gen_crossed(2100000L + round(effect * 1e5) + r, effect)
  f <- tryCatch(rasch_efrm(
    d, item_sets = list(scale = sprintf("I%02d", 1:12)),
    groups = c("region", "cohort"), boot_reps = 0),
    error = function(e) NULL)
  if (is.null(f)) { blank["refused"] <- 1; return(blank) }
  if (!isTRUE(f$est$converged)) {
    blank["nonconv"] <- 1
    return(blank)
  }
  tab <- f$phi_factorial_tests
  if (!check_holm(tab)) stop("crossed-EFRM p_adj is not the Holm adjustment")
  c(raw_any = any(tab$p < ALPHA, na.rm = TRUE),
    adjusted_any = any(tab$p_adj < ALPHA, na.rm = TRUE),
    target_adjusted = tab$p_adj[tab$term == "region"],
    refused = 0, nonconv = 0)
}

summarise_efrm <- function(effect, n) {
  z <- do.call(rbind, run_reps(n, function(r) one_efrm(r, effect)))
  ok <- is.finite(z[, "adjusted_any"])
  note <- paste(
    "400 persons, 12 dichotomous items in one set, balanced 2 x 2 crossed",
    "frame factors; three factorial terms; Holm decisions")
  if (effect == 0) {
    rbind(
      sv_row(STUDY, "crossed EFRM null", "any raw factorial p below .05",
             sum(ok), familywise = mean(z[ok, "raw_any"]),
             n_attempted = n, n_refused = sum(z[, "refused"]),
             n_nonconv = sum(z[, "nonconv"]), notes = note),
      sv_row(STUDY, "crossed EFRM null", "Holm factorial familywise error",
             sum(ok), familywise = mean(z[ok, "adjusted_any"]),
             n_attempted = n, n_refused = sum(z[, "refused"]),
             n_nonconv = sum(z[, "nonconv"]), notes = note))
  } else {
    sv_row(STUDY, sprintf("crossed EFRM region unit ratio %.2f", exp(effect)),
           "Holm-adjusted power for the region main effect", sum(ok),
           power = mean(z[ok, "target_adjusted"] < ALPHA), effect = effect,
           n_attempted = n, n_refused = sum(z[, "refused"]),
           n_nonconv = sum(z[, "nonconv"]), notes = note)
  }
}

if (!DIM_ONLY) {
  cat(sprintf("Crossed EFRM null: %d replicates on %d core(s)\n",
              R_E_NULL, CORES))
  rows <- list(summarise_efrm(0, R_E_NULL))
  for (effect in log(c(1.25, 1.50))) {
    cat(sprintf("Crossed EFRM power, unit ratio %.2f: %d replicates\n",
                exp(effect), R_E_POWER))
    rows[[length(rows) + 1L]] <- summarise_efrm(effect, R_E_POWER)
  }

# BTL dependence: randomise presentation orientation, fit exposure,
# carry-over and position jointly, and treat their adjusted probabilities as
# one family. Thirty judges puts carry-over inside its validated boundary.
randomise_roles <- function(d, seed) {
  set.seed(seed)
  flip <- stats::runif(nrow(d)) < 0.5
  a <- d$object_a
  d$object_a[flip] <- d$object_b[flip]
  d$object_b[flip] <- a[flip]
  d
}

add_order <- function(d) {
  d <- d[order(d$judge), ]
  d$order <- stats::ave(seq_len(nrow(d)), d$judge, FUN = seq_along)
  d
}

gen_btl_effect <- function(seed, effect, magnitude) {
  dep <- if (effect == "exposure")
    list(exposure = magnitude, carry_over = 0) else if (effect == "carry_over")
    list(exposure = 0, carry_over = magnitude) else NULL
  d <- simulate_btl(8, 30, 12, dependence = dep, seed = seed)
  truth <- attr(d, "truth")$location
  d <- randomise_roles(d, seed + 1e7)
  if (is.null(d$order)) d <- add_order(d)
  if (effect == "position" && magnitude != 0) {
    set.seed(seed + 2e7)
    lp <- truth[d$object_a] - truth[d$object_b] + magnitude
    d$winner <- ifelse(stats::rbinom(nrow(d), 1, stats::plogis(lp)) == 1,
                       d$object_a, d$object_b)
  }
  d
}

one_btl <- function(r, effect, magnitude) {
  blank <- c(raw_any = NA, adjusted_any = NA, target_adjusted = NA,
             refused = 0, nonconv = 0)
  offset <- match(effect, c("null", "position", "exposure", "carry_over"))
  d <- gen_btl_effect(3100000L + offset * 100000L +
                        round(magnitude * 1e4) + r, effect, magnitude)
  f <- tryCatch(btl(d, "object_a", "object_b", winner = "winner",
                    judge = "judge", order = "order", position = TRUE),
                error = function(e) NULL)
  if (is.null(f) || is.null(f$dependence)) {
    blank["refused"] <- 1
    return(blank)
  }
  if (!isTRUE(f$converged)) {
    blank["nonconv"] <- 1
    return(blank)
  }
  tab <- f$dependence
  if (!check_holm(tab)) stop("BTL dependence p_adj is not the Holm adjustment")
  target <- if (effect == "null") NA_real_ else
    tab$p_adj[tab$effect == effect]
  c(raw_any = any(tab$p < ALPHA, na.rm = TRUE),
    adjusted_any = any(tab$p_adj < ALPHA, na.rm = TRUE),
    target_adjusted = target, refused = 0, nonconv = 0)
}

summarise_btl <- function(effect, magnitude, n) {
  z <- do.call(rbind, run_reps(n, function(r) one_btl(r, effect, magnitude)))
  ok <- is.finite(z[, "adjusted_any"])
  note <- paste(
    "8 objects, 30 judges, 12 comparisons per pair; exposure, carry-over",
    "and first-position effects fitted jointly; Holm decisions")
  if (effect == "null") {
    rbind(
      sv_row(STUDY, "BTL dependence null", "any raw dependence p below .05",
             sum(ok), familywise = mean(z[ok, "raw_any"]),
             n_attempted = n, n_refused = sum(z[, "refused"]),
             n_nonconv = sum(z[, "nonconv"]), notes = note),
      sv_row(STUDY, "BTL dependence null", "Holm dependence familywise error",
             sum(ok), familywise = mean(z[ok, "adjusted_any"]),
             n_attempted = n, n_refused = sum(z[, "refused"]),
             n_nonconv = sum(z[, "nonconv"]), notes = note))
  } else {
    sv_row(STUDY, sprintf("BTL %s effect %.1f", effect, magnitude),
           sprintf("Holm-adjusted power for %s", effect), sum(ok),
           power = mean(z[ok, "target_adjusted"] < ALPHA),
           effect = magnitude, n_attempted = n,
           n_refused = sum(z[, "refused"]),
           n_nonconv = sum(z[, "nonconv"]), notes = note)
  }
}

cat(sprintf("BTL dependence null: %d replicates\n", R_B_NULL))
rows[[length(rows) + 1L]] <- summarise_btl("null", 0, R_B_NULL)
for (effect in c("position", "exposure", "carry_over")) {
  for (magnitude in c(0.3, 0.6)) {
    cat(sprintf("BTL %s power %.1f: %d replicates\n",
                effect, magnitude, R_B_POWER))
    rows[[length(rows) + 1L]] <- summarise_btl(effect, magnitude, R_B_POWER)
  }
}

# The simulator property is deterministic for each draw once the random
# component is projected off the primary locations. Exercise short and long
# object sets, negative and positive correlations, and several seeds.
corr_error <- numeric(0)
for (k in c(3L, 5L, 8L, 15L)) for (rho in c(-1, -0.8, -0.3, 0, 0.4, 0.9)) {
  for (seed in seq_len(50L)) {
    d <- simulate_btl(k, 8, 2, second_attribute = list(rho = rho),
                      seed = 4100000L + k * 1000L + seed)
    tr <- attr(d, "truth")
    corr_error <- c(corr_error, abs(stats::cor(tr$location, tr$location2) - rho))
  }
}
rows[[length(rows) + 1L]] <- sv_row(
  STUDY, "finite-object correlation grid", "maximum absolute correlation error",
  length(corr_error), bias = max(corr_error), n_attempted = length(corr_error),
  n_refused = 0L, n_nonconv = 0L,
  notes = "4 object counts x 6 correlations x 50 seeds")
} else {
  # Preserve every historical row verbatim while replacing only the affected
  # dimensionality cell. sv_bind_rows() adds newer accounting columns without
  # altering the recorded provenance of those earlier rows.
  previous <- utils::read.csv(
    "tools/simval/results/audit-adjusted-dependence.csv",
    stringsAsFactors = FALSE, check.names = FALSE)
  previous <- previous[
    previous$scenario != "BTL exact-orthogonal second attribute", , drop = FALSE]
  rows <- list(previous)
}

# Re-run the strong dimensionality-power cell affected by the corrected
# second-attribute generator. The earlier 60-replicate result used a vector
# whose realised correlation only approached the requested value in expectation.
one_dim <- function(r) {
  blank <- c(hit = NA, refused = 0, nonconv = 0, error = 0, withheld = 0)
  d <- simulate_btl(15, 20, 60, object_sd = 1.5,
                    second_attribute = list(rho = 0), seed = 5100000L + r)
  f <- tryCatch(btl(d, "object_a", "object_b", winner = "winner",
                    judge = "judge"), error = identity)
  if (inherits(f, "condition")) { blank["refused"] <- 1; return(blank) }
  if (!isTRUE(f$converged)) { blank["nonconv"] <- 1; return(blank) }
  z <- tryCatch(btl_dimensionality(
    f, reps = 100, seed = 6100000L + r,
    # Comparisons share judges. This is a model-conditional power check and
    # does not claim that the allocation itself proves independence.
    independent_comparisons = TRUE), error = identity)
  if (inherits(z, "condition")) { blank["error"] <- 1; return(blank) }
  if (is.na(z$leading_structured)) {
    blank["withheld"] <- 1
    return(blank)
  }
  blank["hit"] <- as.numeric(z$leading_structured)
  blank
}
cat(sprintf("BTL dimensionality power: %d replicates\n", R_DIM))
dim_result <- do.call(rbind, run_reps(R_DIM, one_dim))
ok_dim <- is.finite(dim_result[, "hit"])
rows[[length(rows) + 1L]] <- sv_row(
  STUDY, "BTL exact-orthogonal second attribute",
  "power of leading-structured dimensionality decision", sum(ok_dim),
  power = mean(dim_result[ok_dim, "hit"]), n_attempted = R_DIM,
  n_refused = sum(dim_result[, "refused"]),
  n_nonconv = sum(dim_result[, "nonconv"]),
  n_error = sum(dim_result[, "error"]),
  n_withheld = sum(dim_result[, "withheld"]),
  notes = paste("15 objects, 20 judges, 60 comparisons per pair, object SD 1.5,",
                "rho 0 exactly, 100 null draws per diagnostic; explicit",
                "model-conditional independence opt-in"))

sv_write(do.call(sv_bind_rows, rows), STUDY)
