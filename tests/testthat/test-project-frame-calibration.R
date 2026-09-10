.frame_project <- function(fit, data) .seal_app_project(list(
  format = "rasch-shiny-project", schema = 2L, data = data,
  model_type = .app_fit_family(fit), base_fit = fit,
  rasch_steps = list(), btl_steps = list(), kept_fits = list(),
  kept_fit_code = list(), simulation = list(), results = list(),
  settings = list(), resources = list()))

test_that("frame algorithm records protect every saved-fit location", {
  d <- simulate_efrm(70, 5, n_sets = 1, n_groups = 2, seed = 8651)
  ef <- rasch_efrm(d, attr(d, "truth")$item_sets, "group", id = "id",
                   boot_reps = 0, workers = 1)
  b <- simulate_btl_efrm(n_objects_per_set = 4, n_sets = 1, n_panels = 1,
                         seed = 8621)
  bf <- btl_efrm(b, "object_a", "object_b", "winner", "judge", "panel",
                 attr(b, "truth")$object_sets, se_method = "conditional")
  path <- tempfile(fileext = ".rasch")
  on.exit(unlink(path))
  for (j in 1:2) {
    fit <- list(ef, bf)[[j]]; data <- list(d, b)[[j]]
    expect_identical(fit$calibration_algorithm, "frame-likelihood-1")
    p <- .frame_project(fit, data)
    .save_app_project(p, path)
    expect_identical(.read_app_project(path), p)
    old <- fit; old$calibration_algorithm <- NULL
    expect_error(.validate_app_fit(old), "unverified frame calibration")
    entry <- list(type = "test", label = "Earlier fit", details = list(),
                  code = "fit <- NULL", created = "2026-09-10", fit = old)
    for (where in c("base", "history", "kept", "derived")) {
      legacy <- p
      if (where == "base") legacy$base_fit <- old
      if (where == "history") {
        history <- if (inherits(fit, "rasch_btl")) "btl_steps" else "rasch_steps"
        legacy[[history]] <- list(entry)
      }
      if (where == "kept") legacy$kept_fits <- list(previous = old)
      if (where == "derived") legacy$results$btl_frames <- list(fit = old)
      legacy <- .seal_app_project(legacy)
      saveRDS(legacy, path)
      before <- tools::md5sum(path)
      expect_error(.read_app_project(path), "unverified frame calibration")
      expect_identical(tools::md5sum(path), before)
      expect_identical(readRDS(path)$data, data)
      expect_identical(readRDS(path)$settings, p$settings)
    }
    legacy <- .frame_project(old, data)
    legacy$schema <- 1L; legacy$binding <- NULL
    saveRDS(legacy, path)
    expect_error(.read_app_project(path), "refit this analysis")
    old$calibration_algorithm <- "unknown"
    expect_error(.validate_app_fit(old), "unverified frame calibration")
    # A changed record in a signed file is diagnosed as tampering first.
    p$base_fit <- old
    saveRDS(p, path)
    expect_error(.read_app_project(path), "changed since they were saved")
  }
})

test_that("ordinary saved calibrations do not need frame records", {
  d <- simulate_rasch(100, 6, seed = 8521)
  f <- rasch(d, id = "id")
  expect_null(f$calibration_algorithm)
  expect_no_error(.validate_app_project(.frame_project(f, d)))
})

test_that("the app does not install unverified frame fits over the active analysis", {
  for (pkg in c("shiny", "bslib", "DT", "bsicons")) skip_if_not_installed(pkg)
  app <- testthat::test_path("..", "..", "inst", "shiny", "app.R")
  if (!file.exists(app)) app <- system.file("shiny", "app.R", package = "rasch")
  e <- new.env(parent = globalenv())
  suppressWarnings(sys.source(app, envir = e))
  d <- simulate_efrm(70, 5, n_sets = 1, n_groups = 2, seed = 8651)
  f <- rasch_efrm(d, attr(d, "truth")$item_sets, "group", id = "id",
                  boot_reps = 0, workers = 1)
  f$calibration_algorithm <- NULL
  path <- tempfile(fileext = ".rasch")
  on.exit(unlink(path))
  saveRDS(.frame_project(f, d), path)
  current <- rasch(simulate_rasch(100, 6, seed = 8621))
  shiny::testServer(e$server, {
    fit_val(current)
    session$setInputs(project_file = list(datapath = path, name = "old.rasch",
      size = file.info(path)$size, type = "application/octet-stream"))
    session$flushReact()
    expect_identical(fit_val(), current)
    expect_null(btl_fit())
  })
})
