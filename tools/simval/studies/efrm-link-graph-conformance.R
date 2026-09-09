# End-to-end checks of EFRM linking and bootstrap accounting on four designs.
# One fixed dataset per design and unit ratio: this checks numerical
# conformance, not Type I error, power, or confidence-interval coverage.
# Run from the package root with Rscript. Each fit uses 60 hybrid draws.
# Cases use fixed seeds and run on up to four processes where supported;
# each model's bootstrap is serial.

pkgload::load_all(".", quiet = TRUE, compile = FALSE)
source("tools/simval/harness.R")

designs <- c("two_sets", "complete_triangle", "pairwise_booklets", "chain")
cases <- expand.grid(design = designs, ratio = c(1, 1.3),
                     stringsAsFactors = FALSE)
run_case <- function(k) {
  design <- cases$design[k]
  ratio <- cases$ratio[k]
  S <- if (design == "two_sets") 2L else 3L
  dat <- simulate_efrm(n_per_group = 600, items_per_set = 5,
                       n_sets = S, n_groups = 1, set_unit_ratio = ratio,
                       group_unit_ratio = 1, seed = 27109)
  sets <- attr(dat, "truth")$item_sets
  if (design == "pairwise_booklets") {
    missing_set <- rep(c(3L, 1L, 2L), each = 200L)
    for (s in seq_len(S)) dat[missing_set == s, sets[[s]]] <- NA
  } else if (design == "chain") {
    dat[1:300, sets[[3]]] <- NA
    dat[301:600, sets[[1]]] <- NA
  }
  cat(design, "ratio", ratio, "\n")
  elapsed <- system.time({
    f <- rasch_efrm(dat, item_sets = sets, groups = "group", id = "id",
                    boot_reps = 60, workers = 1, seed = 9814)
  })[["elapsed"]]
  n_edges <- if (design == "chain") 2L else choose(S, 2L)
  C <- f$unit_cov$cov_log_alpha
  log_ratio <- log(f$alpha_table$alpha[S] / f$alpha_table$alpha[1L])
  ratio_se <- sqrt(C[S, S] + C[1L, 1L] - 2 * C[1L, S])
  stopifnot(isTRUE(f$est$converged), nrow(f$linking$alpha_edges) == n_edges,
            all(f$linking$alpha_edges$converged %in% TRUE),
            f$boot_reps_requested == 60L, f$boot_reps_used >= 31L,
            f$boot_reps_used + f$boot_reps_failed == 60L,
            is.finite(ratio_se), ratio_se > 0,
            nrow(f$thresholds) == S * 5L,
            all(is.finite(f$thresholds$tau)),
            nrow(f$thresholds_arbitrary) == S * 5L,
            all(is.finite(f$thresholds_arbitrary$delta)))
  z <- sv_row(
    study = "efrm-link-graph-conformance",
    scenario = paste(design, ratio, sep = "_"),
    quantity = "log_alpha_last_over_first", n_reps = 1L,
    n_attempted = 1L, n_refused = 0L, n_nonconv = 0L, n_error = 0L,
    n_boot_attempted = 60L, n_boot_used = f$boot_reps_used,
    effect = log(ratio), mean_se = ratio_se,
    notes = paste0("One fixed dataset; numerical and accounting checks only. ",
      "No coverage or rejection claim. B_failed combines numerical, ",
      "non-convergence and support failures; no reason-specific totals."))
  z$estimate <- log_ratio
  z$n_edges <- n_edges
  z$n_boot_failed <- f$boot_reps_failed
  z$elapsed_seconds <- elapsed
  z$pass <- TRUE
  z
}
cores <- if (.Platform$OS.type == "windows") 1L else
  min(4L, rasch:::.rasch_available_workers())
rows <- parallel::mclapply(seq_len(nrow(cases)), run_case, mc.cores = cores,
                           mc.set.seed = FALSE)
if (any(vapply(rows, inherits, logical(1), what = "try-error")))
  stop("A conformance case failed; results were not replaced")
sv_write(do.call(rbind, rows), "efrm-link-graph-conformance")
