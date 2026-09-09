.app_sim_bundle_path <- function() {
  p <- testthat::test_path("..", "..", "inst", "shiny", "app.R")
  if (!file.exists(p)) p <- system.file("shiny", "app.R", package = "rasch")
  p
}

test_that("a simulation bundle recreates data and truth in a fresh R process", {
  skip_on_cran()
  for (pkg in c("shiny", "bslib", "DT", "bsicons", "rasch"))
    skip_if_not_installed(pkg)

  e <- new.env(parent = globalenv())
  suppressWarnings(sys.source(.app_sim_bundle_path(), envir = e))
  package_root <- normalizePath(testthat::test_path("..", ".."),
                                mustWork = TRUE)
  source_tree <- file.exists(file.path(package_root, "DESCRIPTION")) &&
    file.exists(file.path(package_root, "R", "rasch.R"))
  installed_candidates <- file.path(.libPaths(), "rasch")
  installed_candidates <- installed_candidates[
    file.exists(file.path(installed_candidates, "DESCRIPTION"))]
  if (source_tree)
    installed_candidates <- installed_candidates[
      !vapply(installed_candidates, function(p)
        normalizePath(p, winslash = "/", mustWork = TRUE) ==
          normalizePath(package_root, winslash = "/", mustWork = TRUE),
        logical(1))]
  if (!length(installed_candidates))
    skip("an installed rasch package is required for the fresh child")
  installed_rasch <- installed_candidates[[1L]]
  # Make child processes use the same installed package selected by the
  # coordinator, rather than whichever older copy happens to come first in
  # the child session's default library search path.
  rstr <- function(x) encodeString(x, quote = "\"")
  child_paths <- unique(c(dirname(installed_rasch), .libPaths()))
  path_setup <- sprintf(".libPaths(c(%s))",
    paste(vapply(child_paths, rstr, character(1)), collapse = ", "))
  bundle <- tempfile(fileext = ".zip")
  unpacked <- tempfile("rasch-simulation-bundle-")
  runner <- tempfile(fileext = ".R")
  current_runner <- tempfile(fileext = ".R")
  result <- tempfile(fileext = ".rds")
  current_result <- tempfile(fileext = ".rds")
  on.exit(unlink(c(bundle, unpacked, runner, current_runner,
                  result, current_result), recursive = TRUE),
          add = TRUE)

  shiny::testServer(e$server, {
    session$setInputs(
      sim_layout = "rasch_exp", sim_seed = 17,
      sr_persons = 30, sr_items = 6, sr_model = "dichotomous", sr_cats = 4,
      sr_mean = 0, sr_sd = 1, sr_dist = "normal", sr_diff = c(-2, 2),
      sr_over = 0, sr_under = 0, sr_guess = FALSE, sr_2d = FALSE,
      sr_rho = 0.3, sr_dep = FALSE, sr_dif = FALSE, sr_difmag = 1,
      sr_style = FALSE, sr_styletype = "extreme", sr_speeded = 0,
      sr_careless = 0, sr_missing = 0, sx_cont = 1, sx_cat = 0.5,
      sx_interaction = TRUE, sx_int = 0.4, sx_depart = 0.7, sim_go = 1)
    session$flushReact()
    original_data <- sim_data()
    original_truth <- sim_truth_val()
    write_sim_bundle(bundle)
    utils::unzip(bundle, exdir = unpacked)
    script <- file.path(unpacked, "simulation.R")
    expect_true(file.exists(script))
    expect_match(paste(readLines(script, warn = FALSE), collapse = "\n"),
                 "library\\(rasch\\)", perl = TRUE)

    # Run only the exported script in a fresh child.  In particular, do not
    # load the source tree first: that would mask a missing library(rasch)
    # header by putting simulator functions on the search path.
    writeLines(c(
      path_setup,
      sprintf("source(%s)", rstr(script)),
      sprintf("saveRDS(list(data = data, truth = attr(data, 'truth')), %s)",
              rstr(result))), runner)
    log <- tempfile(fileext = ".log")
    on.exit(unlink(log), add = TRUE)
    status <- system2(file.path(R.home("bin"), "Rscript"),
                      c("--vanilla", shQuote(runner)),
                      stdout = log, stderr = log)
    expect_equal(status, 0L,
                 info = paste(readLines(log, warn = FALSE), collapse = "\n"))
    regenerated <- readRDS(result)
    expect_s3_class(regenerated$data, "rasch_sim")
    # This first child used only the installed package, with no package
    # pre-attached by the test process: it is the standalone-header check.
    expect_equal(lapply(regenerated$data, identity),
                 lapply(original_data, identity))

    # Also replay with the current source tree (without compilation) so the
    # bundle's full truth object is compared against the just-generated app
    # truth, even when the local installed package predates this checkout.
    current_log <- tempfile(fileext = ".log")
    on.exit(unlink(current_log), add = TRUE)
    current_lines <- c(
      path_setup,
      if (source_tree)
        sprintf("pkgload::load_all(%s, quiet = TRUE, compile = FALSE)",
                rstr(package_root)),
      sprintf("source(%s)", rstr(script)),
      sprintf("saveRDS(list(data = data, truth = attr(data, 'truth')), %s)",
              rstr(current_result)))
    writeLines(current_lines, current_runner)
    current_status <- system2(file.path(R.home("bin"), "Rscript"),
                              c("--vanilla", shQuote(current_runner)),
                              stdout = current_log, stderr = current_log)
    expect_equal(current_status, 0L,
                 info = paste(readLines(current_log, warn = FALSE),
                              collapse = "\n"))
    current <- readRDS(current_result)
    expect_equal(lapply(current$data, identity), lapply(original_data, identity))
    expect_equal(current$truth, original_truth)
    exported <- read.csv(file.path(unpacked, "data.csv"),
                         check.names = FALSE, stringsAsFactors = FALSE)
    # data.csv intentionally cannot carry the class/truth/predictor
    # attributes; compare its columns separately from the classed result.
    expect_equal(lapply(current$data, identity),
                 lapply(exported, identity))
    expect_equal(readRDS(file.path(unpacked, "truth.rds")), original_truth)
  })
})
