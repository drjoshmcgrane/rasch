#!/usr/bin/env Rscript
# Fixed easy/hard subsets under a unidimensional Rasch generator. The
# uncalibrated binomial flag is descriptive, not the package's verdict.
suppressWarnings(pkgload::load_all(".", quiet = TRUE))
source("tools/simval/harness.R")
nrep <- as.integer(Sys.getenv("SV_REPS", "100"))
B <- as.integer(Sys.getenv("SV_B", "39"))
cores <- as.integer(Sys.getenv("SV_CORES", "2"))
stopifnot(nrep >= 20L, B >= 19L, cores >= 1L)
one <- function(r) {
  set.seed(910000L + r)
  n <- 300L; d <- rep(c(-1.5, 1.5), each = 10L)
  X <- matrix(rbinom(n * 20, 1, plogis(outer(rnorm(n), d, "-"))), n, 20,
              dimnames = list(NULL, paste0("I", 1:20)))
  f <- tryCatch(rasch(X), error = identity)
  if (inherits(f, "error")) return(list(status = "refused"))
  if (!isTRUE(f$est$converged)) return(list(status = "nonconv"))
  z <- tryCatch(dimensionality_test(f, items_positive = paste0("I", 1:10),
    items_negative = paste0("I", 11:20), B = B, workers = 1L,
    seed = 1910000L + r), error = identity)
  if (inherits(z, "error") || is.na(z$multidimensional))
    return(list(status = "refused"))
  list(status = "analysed", binomial = z$binomial_multidimensional,
       bootstrap = z$multidimensional, prop = z$prop_significant,
       B_used = z$bootstrap$B_used,
       B_nonconv = z$bootstrap$B_nonconverged,
       B_errors = z$bootstrap$B_errors)
}
z <- if (.Platform$OS.type != "windows" && cores > 1L) {
  parallel::mclapply(seq_len(nrep), one, mc.cores = cores, mc.set.seed = FALSE)
} else lapply(seq_len(nrep), one)
status <- vapply(z, `[[`, "", "status")
ok <- status == "analysed"
stopifnot(any(ok))
val <- function(nm) vapply(z[ok], `[[`, 0, nm)
rows <- lapply(c("binomial", "bootstrap"), function(method) sv_row(
  study = "fixed-split targeting", scenario = sprintf(
    "300 persons, 10+10 dichotomous items at -1.5/+1.5, B=%d", B),
  quantity = method, n_reps = sum(ok), n_attempted = nrep,
  n_refused = sum(status == "refused"), n_nonconv = sum(status == "nonconv"),
  type1 = mean(val(method)),
  n_boot_attempted = sum(val("B_used") + val("B_nonconv") + val("B_errors")),
  n_boot_used = sum(val("B_used")), n_boot_nonconv = sum(val("B_nonconv")),
  n_boot_errors = sum(val("B_errors")),
  notes = paste("Conditional on analysed datasets. Binomial flag descriptive.",
    "This is a targeted calibration screen, not validation across all designs.",
    sprintf("Mean significant-person proportion %.4f.", mean(val("prop"))))))
sv_write(do.call(rbind, rows), "dimensionality-targeting")
print(do.call(rbind, rows))
