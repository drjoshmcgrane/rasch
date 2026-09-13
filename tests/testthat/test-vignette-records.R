.vignette_record_dir <- function() {
  path <- testthat::test_path("..", "..", "vignettes")
  if (!file.exists(file.path(path, "precomputed.R")))
    path <- system.file("doc", package = "rasch")
  path
}

test_that("vignette records retain matching fits and bootstrap results", {
  previous <- setwd(.vignette_record_dir())
  on.exit(setwd(previous), add = TRUE)
  e <- new.env(parent = globalenv())
  sys.source("precomputed.R", e)
  original_hash <- tools::md5sum("precomputed/rasch-workflow.rds")
  z <- e$vignette_result("rasch-workflow")
  expect_identical(z$dimensionality$B, 199L)
  expect_identical(z$dimensionality$bootstrap$B_used, 199L)
  expect_no_error(.validate_dimensionality_test(z$dimensionality, z$fit))
  expect_no_error(.validate_scree_result(z$scree, z$fit))
  expect_identical(tools::md5sum("precomputed/rasch-workflow.rds"), original_hash)
  expect_equal(z$fit$X, as.matrix(z$d[, z$fit$items$item]), ignore_attr = TRUE)
  ef <- e$vignette_result("extended-frame-reference")$fit
  expect_identical(ef$boot_reps_requested, 50L)
  expect_true(ef$boot_reps_used > 0)
  expect_s3_class(e$vignette_result("paired-comparisons")$dimensions,
                  "rasch_btl_dim")
})

test_that("vignette input hashes are independent of source line endings", {
  e <- new.env(parent = globalenv())
  sys.source(file.path(.vignette_record_dir(), "precomputed.R"), e)
  a <- tempfile(); b <- tempfile()
  on.exit(unlink(c(a, b)), add = TRUE)
  writeBin(charToRaw("a\nb\n"), a)
  writeBin(charToRaw("a\r\nb\r\n"), b)
  expect_identical(unname(e$vignette_source_hashes(a)),
                   unname(e$vignette_source_hashes(b)))
})

test_that("vignette readers refuse changed versions and result files", {
  source <- .vignette_record_dir()
  tmp <- tempfile("vignette-records-")
  dir.create(tmp)
  on.exit(unlink(tmp, recursive = TRUE), add = TRUE)
  file.copy(file.path(source, "precomputed"), tmp, recursive = TRUE)
  e <- new.env(parent = globalenv())
  sys.source(file.path(source, "precomputed.R"), e)
  previous <- setwd(tmp)
  # Leave the temporary directory before deleting it, including on Windows.
  on.exit(setwd(previous), add = TRUE, after = FALSE)
  path <- "precomputed/manifest.rds"
  manifest <- readRDS(path)
  changed <- manifest; changed$version <- "0.0.0"
  saveRDS(changed, path)
  expect_error(e$vignette_result("rasch-workflow"), "package version")
  changed <- manifest; changed$verified_fit_pairs <- NULL
  saveRDS(changed, path)
  expect_error(e$vignette_result("rasch-workflow"), "pairing needs verification")
  saveRDS(manifest, path)
  saveRDS(list(), "precomputed/rasch-workflow.rds")
  expect_error(e$vignette_result("rasch-workflow"), "missing or have changed")
})
