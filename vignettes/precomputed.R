# Shared by the three vignettes with precomputed bootstrap calculations.
# The displayed analysis code remains executable with recompute = TRUE.
recompute <- identical(Sys.getenv("RASCH_REBUILD_VIGNETTES"), "true")

vignette_source_hashes <- function(paths) {
  tmp <- tempfile("rasch-vignette-hash-")
  on.exit(unlink(tmp), add = TRUE)
  vapply(paths, function(path) {
    lines <- readLines(path, warn = FALSE, encoding = "UTF-8")
    writeBin(charToRaw(paste0(paste(enc2utf8(lines), collapse = "\n"), "\n")), tmp)
    unname(tools::md5sum(tmp))
  }, "")
}

vignette_result <- function(name) {
  manifest <- readRDS(file.path("precomputed", "manifest.rds"))
  if (!identical(manifest$version, as.character(utils::packageVersion("rasch"))))
    stop("Precomputed vignette results need updating for this package version.")
  # In a source build, also refuse results from older code or vignette inputs.
  # Installed vignette sources do not contain the package's R/ and src/ trees.
  if (file.exists(file.path("..", "R", "rasch.R"))) {
    paths <- file.path("..", names(manifest$source_hashes))
    if (!all(file.exists(paths)) ||
        !identical(unname(vignette_source_hashes(paths)),
                   unname(manifest$source_hashes)))
      stop("Precomputed vignette inputs have changed; run tools/vignettes/precompute.R.")
  }
  path <- file.path("precomputed", paste0(name, ".rds"))
  if (!name %in% names(manifest$result_hashes) || !file.exists(path) ||
      !identical(unname(tools::md5sum(path)), unname(manifest$result_hashes[name])))
    stop("Precomputed vignette results are missing or have changed: ", name)
  result <- readRDS(path)
  if (name == "rasch-workflow") {
    if (!name %in% manifest$verified_fit_pairs)
      stop("Vignette fit/result pairing needs verification; regenerate the records.")
    # Runtime signatures include R's version and native encoding. The file
    # checksum above authenticates the complete fit/result bundle, whose
    # pairing the generator verifies before recording it. Refresh only its
    # in-memory signatures so these fixed records also work in other R builds.
    signature <- getFromNamespace(".fit_boot_signature", "rasch")(result$fit)
    hash <- getFromNamespace(".fit_boot_md5", "rasch")
    for (element in c("dimensionality", "scree")) {
      attr(result[[element]], "fit_signature") <- signature
      attr(result[[element]], "result_signature") <- NULL
      attr(result[[element]], "result_signature") <- hash(result[[element]])
    }
  }
  result
}
