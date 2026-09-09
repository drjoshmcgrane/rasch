.weighted_migration_fixture <- function() {
  fit <- rasch(simulate_rasch(120, 6, seed = 9038))
  # The active calibration differs from the base fit after an item change.
  active <- drop_items(fit, colnames(fit$X)[6L])
  weights <- setNames(c(2, 1, 1, .5, .5), colnames(active$X))
  weighted <- list(table = weighted_person_estimates(active, weights),
    weights = weights, by = "item", sets = NULL, filename = "weights.csv",
    fit_signature = .fit_boot_signature(active))
  weighted$result_signature <- .fit_boot_md5(weighted)
  project <- .seal_app_project(list(
    format = "rasch-shiny-project", schema = 2L, package_version = "test",
    created = "2026-09-07", data = as.data.frame(fit$X), model_type = "rasch",
    base_fit = fit, rasch_steps = list(list(type = "drop", label = "Drop item",
      details = list(), code = "fit <- drop_items(fit, 6)",
      created = "2026-09-07", fit = active)), btl_steps = list(),
    rcode = "fit <- rasch(data)", kept_fits = list(), kept_fit_code = list(),
    settings = list(model_type = "rasch", item_cols = colnames(fit$X)),
    resources = list(), simulation = list(), results = list(person_weights = weighted)))
  list(project = project, active = active, table = weighted$table)
}

.older_weighted_project <- function(project, algorithm = NULL) {
  weighted <- project$results$person_weights
  attr(weighted$table, "algorithm") <- algorithm
  # Model the final-bit rounding change and an unavailable old-solver score.
  weighted$table$theta[1L] <- weighted$table$theta[1L] + 1e-14
  weighted$table$theta[2L] <- weighted$table$se[2L] <- NA_real_
  weighted$result_signature <- NULL
  weighted$result_signature <- .fit_boot_md5(weighted)
  project$results$person_weights <- weighted
  .seal_app_project(project)
}

test_that("older weighted tables are authenticated and recomputed on the active fit", {
  d <- .weighted_migration_fixture()
  project <- .older_weighted_project(d$project)
  path <- tempfile(fileext = ".rasch")
  on.exit(unlink(path))
  saveRDS(project, path)
  expect_warning(restored <- .read_app_project(path),
                 "weighted person estimates were recomputed")
  expect_identical(restored$results$person_weights$table, d$table)
  expect_identical(restored$base_fit, project$base_fit)
  expect_identical(restored$rasch_steps, project$rasch_steps)
  expect_identical(restored$data, project$data)
  expect_identical(restored$results$person_weights$weights,
                   project$results$person_weights$weights)
  expect_no_error(.validate_app_project(restored))
  .save_app_project(restored, path)
  expect_no_warning(again <- .read_app_project(path))
  expect_identical(again, restored)

  # Schema 1 cannot authenticate an old result; retain its fit but omit it.
  project$schema <- 1L; project$binding <- NULL
  saveRDS(project, path)
  expect_warning(legacy <- .read_app_project(path), "schema-1")
  expect_null(legacy$results$person_weights)
  expect_identical(legacy$rasch_steps, project$rasch_steps)
})

test_that("weighted-table migration does not bypass integrity or fit matching", {
  d <- .weighted_migration_fixture()
  old <- .older_weighted_project(d$project)
  path <- tempfile(fileext = ".rasch")
  on.exit(unlink(path))
  bad <- old
  bad$results$person_weights$table$theta[1L] <- 99
  saveRDS(bad, path)
  expect_error(.read_app_project(path), "changed since they were saved")
  # Resealing the project alone cannot repair the inner result signature.
  saveRDS(.seal_app_project(bad), path)
  expect_error(.read_app_project(path), "changed since it was calculated")
  bad <- old
  bad$results$person_weights$fit_signature <- .fit_boot_signature(bad$base_fit)
  bad$results$person_weights$result_signature <- NULL
  bad$results$person_weights$result_signature <- .fit_boot_md5(bad$results$person_weights)
  saveRDS(.seal_app_project(bad), path)
  expect_error(.read_app_project(path), "different fitted model")

  # Current tables keep their exact-reproduction check, even if resealed.
  bad <- d$project
  bad$results$person_weights$table$theta[1L] <- 99
  bad$results$person_weights$result_signature <- NULL
  bad$results$person_weights$result_signature <- .fit_boot_md5(bad$results$person_weights)
  saveRDS(.seal_app_project(bad), path)
  expect_error(.read_app_project(path), "table does not reproduce")
})

test_that("older weighted algorithm stamps migrate but unknown stamps do not", {
  d <- .weighted_migration_fixture()
  path <- tempfile(fileext = ".rasch")
  on.exit(unlink(path))
  for (stamp in c("pattern-wle-1", "pattern-unit-wle-2")) {
    old <- .older_weighted_project(d$project, stamp)
    saveRDS(old, path)
    expect_warning(restored <- .read_app_project(path),
                   "weighted person estimates were recomputed")
    expect_identical(restored$results$person_weights$table, d$table)
    expect_identical(attr(restored$results$person_weights$table,
                          "algorithm"), "pattern-max-wle-3")
    expect_identical(restored$rasch_steps, old$rasch_steps)
    expect_no_error(.validate_app_project(restored))
    .save_app_project(restored, path)
    expect_no_warning(.read_app_project(path))
  }

  # A recognised old stamp never excuses a changed inner result signature.
  old <- .older_weighted_project(d$project, "pattern-unit-wle-2")
  old$results$person_weights$table$theta[1L] <- 99
  saveRDS(.seal_app_project(old), path)
  expect_error(.read_app_project(path), "changed since it was calculated")
  unknown <- .older_weighted_project(d$project, "unrecognised-solver")
  saveRDS(unknown, path)
  expect_error(.read_app_project(path), "table does not reproduce")
})
