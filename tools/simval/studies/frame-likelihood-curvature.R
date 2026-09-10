# Deterministic likelihood, curvature and bootstrap-accounting regressions.
# Run from the package root. No coverage, Type I error or power claim.
pkgload::load_all(".", quiet = TRUE, compile = FALSE)
source("tools/simval/harness.R")
test_file <- "tests/testthat/test-frame-likelihood-curvature.R"
test_hash <- unname(tools::md5sum(test_file))
results <- testthat::test_local(filter = "^frame-likelihood-curvature$",
  reporter = "summary", stop_on_failure = FALSE)
rows <- lapply(results, function(x) {
  classes <- vapply(x$results, function(y) class(y)[1L], "")
  counts <- table(factor(classes, levels = c("expectation_success",
    "expectation_failure", "expectation_error", "expectation_warning",
    "expectation_skip")))
  z <- sv_row("frame-likelihood-curvature", x$test, "numerical_conformance",
    n_reps = 1L, n_attempted = 1L, n_refused = 0L, n_nonconv = 0L,
    n_error = as.integer(counts[3L]),
    notes = paste("Fixed-data regression block, not a sampling study.",
      "Checks exact curvature, public refusal and accepted-draw covariance."))
  z$test_file <- test_file
  z$test_md5 <- test_hash
  z$passed <- as.integer(counts[1L])
  z$failed <- as.integer(counts[2L])
  z$errors <- as.integer(counts[3L])
  z$warnings <- as.integer(counts[4L])
  z$skipped <- as.integer(counts[5L])
  z$pass <- all(counts[-1L] == 0L) && counts[1L] > 0L
  z
})
stopifnot(identical(unname(tools::md5sum(test_file)), test_hash))
rows <- do.call(rbind, rows)
sv_write(rows, "frame-likelihood-curvature")
stopifnot(all(rows$pass))
