#!/usr/bin/env Rscript
# Null calibration and power of dimensionality_test(). The fixed split is
# nominated before seeing the responses. The automatic split is selected from
# each fitted residual matrix; its bootstrap repeats that selection in every
# score-conditional replicate.

suppressWarnings(pkgload::load_all(".", quiet = TRUE))
source("tools/simval/harness.R")

NREP <- as.integer(Sys.getenv("SV_REPS", "100"))
B <- as.integer(Sys.getenv("SV_B", "99"))
CORES <- as.integer(Sys.getenv("SV_CORES", "4"))
if (!is.finite(NREP) || NREP < 20L) stop("SV_REPS must be at least 20")
if (!is.finite(B) || B < 19L) stop("SV_B must be at least 19")
if (!is.finite(CORES) || CORES < 1L) CORES <- 1L

scenarios <- data.frame(
  design = c("dichotomous", "PCM", "PCM"),
  departure = c("null", "null", "second dimension"),
  n_items = c(30L, 20L, 20L), stringsAsFactors = FALSE)

one <- function(r, design, departure, n_items, tag) {
  seed <- 970000L + tag * 10000L + r
  second <- if (departure == "second dimension")
    list(items = sprintf("I%02d", (n_items / 2L + 1L):n_items), rho = .3)
  d <- tryCatch(simulate_rasch(
    n_persons = 400, n_items = n_items, model = design,
    n_categories = if (design == "PCM") 4 else 2,
    difficulty = c(-1.8, 1.8), second_dim = second, seed = seed),
    error = identity)
  if (inherits(d, "condition"))
    return(c(status = "error", prop = NA, binomial = NA, fixed = NA,
             bootstrap = NA, B_used = NA, B_nonconv = NA, B_errors = NA))
  items <- sprintf("I%02d", seq_len(n_items))
  fit <- tryCatch(rasch(d, id = "id", items = items,
                        model = if (design == "PCM") "PCM" else "PCM"),
                  error = identity)
  if (inherits(fit, "condition"))
    return(c(status = "refused", prop = NA, binomial = NA, fixed = NA,
             bootstrap = NA, B_used = NA, B_nonconv = NA, B_errors = NA))
  if (!isTRUE(fit$est$converged))
    return(c(status = "nonconv", prop = NA, binomial = NA, fixed = NA,
             bootstrap = NA, B_used = NA, B_nonconv = NA, B_errors = NA))
  half <- n_items / 2L
  fixed <- tryCatch(dimensionality_test(
    fit, items_positive = items[seq_len(half)],
    items_negative = items[(half + 1L):n_items]), error = identity)
  auto <- tryCatch(suppressWarnings(dimensionality_test(
    fit, B = B, workers = 1L, seed = seed + 500000L)), error = identity)
  if (inherits(fixed, "condition") || inherits(auto, "condition"))
    return(c(status = "refused", prop = NA, binomial = NA, fixed = NA,
             bootstrap = NA,
             B_used = if (inherits(auto, "condition")) auto$B_used %||% NA else
               auto$bootstrap$B_used,
             B_nonconv = if (inherits(auto, "condition"))
               auto$B_nonconverged %||% NA else auto$bootstrap$B_nonconverged,
             B_errors = if (inherits(auto, "condition"))
               auto$B_errors %||% NA else auto$bootstrap$B_errors))
  c(status = "analysed", prop = auto$prop_significant,
    binomial = auto$binomial_multidimensional,
    fixed = fixed$binomial_multidimensional, bootstrap = auto$multidimensional,
    B_used = auto$bootstrap$B_used,
    B_nonconv = auto$bootstrap$B_nonconverged,
    B_errors = auto$bootstrap$B_errors)
}

rows <- list()
for (s in seq_len(nrow(scenarios))) {
  sc <- scenarios[s, ]
  f <- function(r) one(r, sc$design, sc$departure, sc$n_items, s)
  z <- if (.Platform$OS.type != "windows" && CORES > 1L)
    parallel::mclapply(seq_len(NREP), f, mc.cores = min(CORES, NREP),
                       mc.set.seed = FALSE) else lapply(seq_len(NREP), f)
  status <- vapply(z, `[[`, "", "status")
  ok <- status == "analysed"
  val_num <- function(name) as.numeric(vapply(z[ok], `[[`, "", name))
  val_flag <- function(name) as.logical(vapply(z[ok], `[[`, "", name))
  prop <- val_num("prop"); binomial <- val_flag("binomial")
  fixed <- val_flag("fixed"); bootstrap <- val_flag("bootstrap")
  B_used <- val_num("B_used"); B_nonconv <- val_num("B_nonconv")
  B_errors <- val_num("B_errors")
  label <- sprintf("%s, %s, 400 persons x %d items, B = %d",
                   sc$design, sc$departure, sc$n_items, B)
  accounting <- list(n_reps = sum(ok), n_attempted = NREP,
    n_refused = sum(status == "refused"),
    n_nonconv = sum(status == "nonconv"), n_error = sum(status == "error"),
    n_boot_attempted = sum(B_used + B_nonconv + B_errors, na.rm = TRUE),
    n_boot_used = sum(B_used, na.rm = TRUE),
    n_boot_nonconv = sum(B_nonconv, na.rm = TRUE),
    n_boot_errors = sum(B_errors, na.rm = TRUE))
  add <- function(quantity, rate, kind, notes = "") {
    args <- c(list(study = "dimensionality-test bootstrap", scenario = label,
                   quantity = quantity), accounting,
              list(notes = notes))
    args[[kind]] <- rate
    rows[[length(rows) + 1L]] <<- do.call(sv_row, args)
  }
  null <- sc$departure == "null"
  add("mean significant person comparisons, automatic split", mean(prop),
      "effect", "descriptive rate, not a rejection probability")
  add("fixed-split binomial decision", mean(fixed),
      if (null) "type1" else "power",
      "descriptive benchmark; the public inferential verdict requires bootstrap calibration")
  add("automatic-split uncalibrated binomial decision", mean(binomial),
      if (null) "type1" else "power",
      "descriptive benchmark; the public inferential verdict is withheld")
  add("automatic-split parametric-bootstrap decision", mean(bootstrap),
      if (null) "type1" else "power")
}

out <- do.call(rbind, rows)
sv_write(out, "dimensionality-bootstrap")
print(out[, c("scenario", "quantity", "n_reps", "effect", "type1", "power",
              "refusal_rate")], row.names = FALSE)
