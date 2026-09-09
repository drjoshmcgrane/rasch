# Fixed-data checks of full-bootstrap identification, failure accounting and
# covariance reconstruction. Not a coverage, power or Type I error study.
# Run from the package root: Rscript tools/simval/studies/efrm-full-identification.R
# Recording wrappers return the solver and bootstrap results unchanged.

pkgload::load_all(".", quiet = TRUE, compile = FALSE)
source("tools/simval/harness.R")

cases <- data.frame(
  scenario = c("balanced_binary", "balanced_polytomous", "group_12",
               "group_3", "two_sets_two_groups"),
  n = c(150L, 150L, 120L, 120L, 150L),
  items = c(5L, 5L, 3L, 5L, 5L),
  sets = c(1L, 1L, 1L, 1L, 2L),
  categories = c(2L, 4L, 2L, 2L, 2L),
  retained_group_2 = c(150L, 150L, 12L, 3L, 150L),
  data_seed = c(9207L, 731L, 1L, 2L, 291L),
  bootstrap_seed = c(889L, 889L, 719L, 719L, 889L))

run_case <- function(k) {
  ca <- cases[k, ]
  d <- simulate_efrm(ca$n, ca$items, n_sets = ca$sets, n_groups = 2,
                     n_categories = ca$categories, group_unit_ratio = 1,
                     seed = ca$data_seed)
  sets <- attr(d, "truth")$item_sets
  d <- d[c(seq_len(ca$n), ca$n + seq_len(ca$retained_group_2)), ]
  stage <- ""
  refits <- list()
  full_draws <- NULL
  original_solve <- rasch:::.efrm_solve
  original_apply <- rasch:::.efrm_boot_apply
  testthat::local_mocked_bindings(
    .efrm_solve = function(...) {
      z <- tryCatch(original_solve(...), error = identity)
      if (identical(stage, "full person bootstrap"))
        refits[[length(refits) + 1L]] <<- z
      if (inherits(z, "error")) stop(z)
      z
    },
    .efrm_boot_apply = function(...) {
      z <- original_apply(...)
      if (identical(stage, "full person bootstrap")) full_draws <<- z
      z
    }, .package = "rasch")
  cat(ca$scenario, "\n")
  elapsed <- system.time({
    f <- rasch_efrm(d, sets, "group", id = "id", se_method = "bootstrap",
                    boot_reps = 40, workers = 1, seed = ca$bootstrap_seed,
                    progress = function(s, current, total) stage <<- s)
  })[["elapsed"]]
  stopifnot(length(full_draws) == 40L, length(refits) == 40L,
            f$est$converged, identical(f$se_method, "bootstrap"))
  unidentified <- vapply(refits, function(z)
    !inherits(z, "error") && any(z$phi_unident), logical(1))
  nonconv <- vapply(refits, function(z)
    !inherits(z, "error") && !isTRUE(z$converged), logical(1))
  solver_error <- vapply(refits, inherits, logical(1), what = "error")
  keep <- vapply(full_draws, function(z)
    !is.null(z) && all(is.finite(z)), logical(1))
  stopifnot(!any(keep & (unidentified | nonconv | solver_error)),
            sum(keep) == f$full_boot_reps_used,
            sum(!keep) == f$full_boot_reps_failed,
            f$full_boot_reps_attempted == 40L)
  if (ca$scenario %in% c("group_3", "group_12"))
    stopifnot(any(unidentified))
  B <- do.call(rbind, full_draws[keep])
  G <- nrow(f$phi_table)
  S <- nrow(f$alpha_table)
  Md <- nrow(f$unit_cov$cov_dtilde)
  ip <- seq_len(G)
  ia <- G + seq_len(S)
  id <- G + 2L * S + seq_len(Md)
  ic <- G + 2L * S + Md + seq_len(Md)
  direct <- list(
    cov_joint = cov(B[, c(id, ip), drop = FALSE]),
    cov_dtilde = cov(B[, id, drop = FALSE]),
    cov_log_phi = cov(B[, ip, drop = FALSE]),
    cov_delta = cov(B[, ic, drop = FALSE]))
  if (S > 1L) {
    direct$cov_log_alpha <- cov(B[, ia, drop = FALSE])
    direct$cov_log_alpha_phi <-
      cov(B[, ia, drop = FALSE], B[, ip, drop = FALSE])
  }
  differences <- vapply(names(direct), function(nm) {
    equal <- all.equal(unname(f$unit_cov[[nm]]),
                       unname(direct[[nm]]), tolerance = 1e-12)
    if (!isTRUE(equal)) stop(nm, ": ", paste(equal, collapse = "; "))
    max(abs(f$unit_cov[[nm]] - direct[[nm]]))
  }, numeric(1))
  # Remove all failed solves before inspecting their returned estimates.
  # This optional comparison quantifies the rejected draws' contribution;
  # it is not a repeated-sampling assessment of either standard error.
  phi_available <- vapply(refits, function(z)
    !inherits(z, "error") && all(is.finite(log(z$phi))), logical(1))
  all_phi <- vapply(refits[phi_available], function(z) log(z$phi[2L]), 0)
  z <- sv_row("efrm-full-identification", ca$scenario,
    "full_bootstrap_covariance_conformance", n_reps = 1L,
    n_attempted = 1L, n_refused = 0L, n_nonconv = 0L, n_error = 0L,
    n_boot_attempted = 40L, n_boot_used = sum(keep),
    n_boot_nonconv = sum(nonconv), mean_se = f$phi_table$se_log_phi[2L],
    notes = paste("One fixed dataset; recording wrappers do not alter fits.",
      "Covariance reconstructed from accepted draws. No coverage or size claim.",
      "Unidentified and non-converged counts may overlap."))
  z$n_boot_failed <- sum(!keep)
  z$n_boot_unidentified <- sum(unidentified)
  z$n_boot_solver_error <- sum(solver_error)
  z$n_boot_other_failed <- sum(!keep & !unidentified & !nonconv & !solver_error)
  z$max_abs_covariance_difference <- max(differences)
  z$unfiltered_se_log_phi <- sd(all_phi)
  z$phi_inference <- f$unit_support$phi_inference
  z$elapsed_seconds <- elapsed
  z$pass <- TRUE
  z
}
rows <- lapply(seq_len(nrow(cases)), run_case)
sv_write(do.call(rbind, rows), "efrm-full-identification")
