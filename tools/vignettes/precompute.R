# Run from the package root. --check verifies the recorded inputs without
# repeating the calculations. Normal execution refreshes all three records.
args <- commandArgs(trailingOnly = TRUE)
stopifnot(file.exists("DESCRIPTION"), dir.exists("vignettes"))
repo <- normalizePath(".")
helper <- new.env(parent = globalenv())
sys.source("vignettes/precomputed.R", helper)
documents <- c("rasch-workflow", "extended-frame-reference", "paired-comparisons")
sources <- sort(c(
  list.files("R", pattern = "[.]R$", full.names = TRUE),
  list.files("src", pattern = "[.](cpp|h)$", full.names = TRUE),
  "vignettes/precomputed.R", paste0("vignettes/", documents, ".Rmd")))
hashes <- helper$vignette_source_hashes(sources)
version <- read.dcf("DESCRIPTION")[1L, "Version"]
manifest_path <- "vignettes/precomputed/manifest.rds"

if ("--check" %in% args) {
  manifest <- readRDS(manifest_path)
  stopifnot(identical(manifest$version, unname(version)),
            identical(manifest$source_hashes, hashes))
  for (name in names(manifest$result_hashes)) {
    path <- file.path("vignettes/precomputed", paste0(name, ".rds"))
    stopifnot(identical(unname(tools::md5sum(path)),
                        unname(manifest$result_hashes[name])))
  }
  cat("Precomputed vignette inputs and results match the source tree.\n")
} else {
  if (!requireNamespace("pkgload", quietly = TRUE) ||
      !requireNamespace("knitr", quietly = TRUE))
    stop("Install pkgload and knitr to regenerate vignette calculations.")
  pkgload::load_all(repo, quiet = TRUE)
  Sys.setenv(RASCH_REBUILD_VIGNETTES = "true")
  dir.create("vignettes/precomputed", showWarnings = FALSE)
  scratch <- tempfile("rasch-vignettes-")
  dir.create(scratch)
  variables <- list(
    `rasch-workflow` = c("d", "fit", "dimensionality", "scree"),
    `extended-frame-reference` = "fit",
    `paired-comparisons` = "dimensions")
  timings <- numeric(length(documents)); names(timings) <- documents
  for (name in documents) {
    env <- new.env(parent = globalenv())
    knitr::opts_chunk$set(error = FALSE, fig.path = paste0(scratch, "/", name, "-"))
    knitr::opts_knit$set(root.dir = file.path(repo, "vignettes"))
    started <- proc.time()
    knitr::knit(file.path(repo, "vignettes", paste0(name, ".Rmd")),
                output = file.path(scratch, paste0(name, ".md")),
                envir = env, quiet = TRUE)
    timings[name] <- unname((proc.time() - started)["elapsed"])
    result <- mget(variables[[name]], envir = env, inherits = FALSE)
    path <- file.path("vignettes/precomputed", paste0(name, ".rds"))
    # Format 3 preserves the vector representations used by the package's
    # fitted-model fingerprints. Format 2 expands them and changes the hash.
    saveRDS(result, path, compress = "xz", version = 3)
    restored <- readRDS(path)
    if (name == "rasch-workflow") {
      getFromNamespace(".validate_dimensionality_test", "rasch")(
        restored$dimensionality, restored$fit)
      getFromNamespace(".validate_scree_result", "rasch")(
        restored$scree, restored$fit)
      stopifnot(restored$dimensionality$bootstrap$B_used == 199L)
    }
    cat(name, ":", round(timings[name], 1), "seconds\n")
  }
  stopifnot(identical(hashes, helper$vignette_source_hashes(sources)))
  result_hashes <- vapply(documents, function(name)
    unname(tools::md5sum(file.path("vignettes/precomputed", paste0(name, ".rds")))), "")
  manifest <- list(version = unname(version), source_hashes = hashes,
    result_hashes = result_hashes, executed = format(Sys.time(), tz = "UTC"),
    r_version = R.version.string, timings = timings,
    verified_fit_pairs = "rasch-workflow")
  saveRDS(manifest, manifest_path, compress = "xz", version = 2)
  cat("Precomputed results written; figures and rendered source:", scratch, "\n")
}
