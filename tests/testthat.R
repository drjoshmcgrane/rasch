library(testthat)
library(rasch)

# Keep CRAN checks short with representative workflows and focused regression
# tests. The complete suite, including repeated fits and calibration studies,
# runs with NOT_CRAN=true locally and in CI on three operating systems.
if (identical(Sys.getenv("NOT_CRAN"), "true")) {
  test_check("rasch")
} else {
  core <- c(
    "app-project", "btl-efrm", "btl-equating", "btl-targeting", "cran-core",
    "fable-regressions", "format", "identification", "mc-scoring",
    "missing", "pcml-pc", "resolve-frames", "vignette-records", "wright-map"
  )
  test_check("rasch", filter = paste0("^(", paste(core, collapse = "|"), ")$"))
}
