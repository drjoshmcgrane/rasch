# Fixed-data checks of frame-link boundaries and saved-calibration validity.
# Run from the package root. These are not size, power or coverage estimates.
pkgload::load_all(".", quiet = TRUE, compile = FALSE)
source("tools/simval/harness.R")
test_files <- file.path("tests/testthat", c("test-btl-efrm-separation.R",
                                           "test-project-frame-calibration.R"))
test_hashes <- tools::md5sum(test_files)
results <- testthat::test_local(
  filter = "^(btl-efrm-separation|project-frame-calibration)$",
  reporter = "summary", stop_on_failure = FALSE)
rows <- lapply(results, function(x) {
  classes <- vapply(x$results, function(y) class(y)[1L], "")
  counts <- table(factor(classes, levels = c("expectation_success",
    "expectation_failure", "expectation_error", "expectation_warning",
    "expectation_skip")))
  z <- sv_row("frame-link-safety", x$test, "numerical_conformance",
    n_reps = 1L, n_attempted = 1L, n_refused = 0L, n_nonconv = 0L,
    n_error = as.integer(counts[3L]),
    notes = paste("Fixed-data regression block, not a sampling study.",
      "Checks separation certificates, finite overlap, convergence and project refusal."))
  z$test_file <- file.path("tests/testthat", x$file)
  z$test_md5 <- unname(test_hashes[z$test_file])
  z$passed <- as.integer(counts[1L]); z$failed <- as.integer(counts[2L])
  z$errors <- as.integer(counts[3L]); z$warnings <- as.integer(counts[4L])
  z$skipped <- as.integer(counts[5L])
  z$pass <- all(counts[-1L] == 0L) && counts[1L] > 0L
  z
})
stopifnot(identical(tools::md5sum(test_files), test_hashes))
rows <- do.call(rbind, rows)
stopifnot(!anyNA(rows$test_md5))
sv_write(rows, "frame-link-safety")
stopifnot(all(rows$pass))
