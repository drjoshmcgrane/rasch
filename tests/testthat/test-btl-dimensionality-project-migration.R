legacy_btl_dimension_reference <- function(result, flag = result$leading_structured) {
  # Keep this fixture current with the residual definition while making its
  # reference metadata deliberately legacy/unsupported.
  result$residual_method <- "pooled-expected-score-1"
  result$reference$inference_available <- NULL
  result$reference$mean <- mean(result$reference$draws)
  result$reference$p95 <- as.numeric(quantile(result$reference$draws, .95))
  result$reference$p <- result$reference$p_adj <- .01
  result$leading_structured <- flag
  result$bimensions$above_reference[1L] <- flag
  result$bimensions$ref_mean[1L] <- result$reference$mean
  result$bimensions$ref_p95[1L] <- result$reference$p95
  attr(result, "result_signature") <- NULL
  attr(result, "result_signature") <- .fit_boot_md5(result)
  result
}

btl_dimension_project <- function(data, base_fit, result, active_fit = NULL) {
  history <- if (is.null(active_fit)) list() else list(list(
    type = "subset", label = "Retained comparisons", details = list(),
    created = "2026-09-07", code = "btl(retained)", fit = active_fit))
  .seal_app_project(list(
    format = "rasch-shiny-project", schema = 2L,
    data = as.data.frame(data), model_type = "btl", base_fit = base_fit,
    rasch_steps = list(), btl_steps = history,
    kept_fits = list(original = base_fit),
    kept_fit_code = list(original = "btl(source_data)"),
    settings = list(display = "objects"),
    results = list(dimensionality = result,
                   display = list(selected_object = base_fit$objects$object[1L]))))
}

test_that("authenticated projects omit only unsupported legacy CJ references", {
  d <- simulate_btl(5, 15, reps_per_pair = 30, seed = 9681)
  pair <- sort(c(d$object_a[1L], d$object_b[1L]))
  removed <- (d$object_a == pair[1L] & d$object_b == pair[2L]) |
    (d$object_a == pair[2L] & d$object_b == pair[1L])
  base <- btl(d, "object_a", "object_b", winner = "winner", judge = "judge")
  active <- btl(d[!removed, ], "object_a", "object_b", winner = "winner",
                  judge = "judge")
  current <- btl_dimensionality(active, reps = 20L, seed = 9691)
  old <- legacy_btl_dimension_reference(current)
  project <- btl_dimension_project(d, base, old, active)
  path <- tempfile(fileext = ".rasch")
  on.exit(unlink(path), add = TRUE)
  saveRDS(project, path)

  expect_warning(restored <- .read_app_project(path),
                 "dimensionality reference.*unsupported comparison design")
  expect_null(restored$results$dimensionality)
  for (field in c("data", "base_fit", "btl_steps", "rasch_steps", "kept_fits",
                  "kept_fit_code", "settings"))
    expect_identical(restored[[field]], project[[field]])
  expect_identical(restored$results$display, project$results$display)
  expect_match(attr(restored, "rasch_project_legacy_dropped"), "dimensionality")
  expect_no_error(.validate_app_project(restored))
  expect_no_error(.save_app_project(restored, path))
  expect_no_warning(.read_app_project(path))

  # A changed bundle must fail before migration can replace its seal.
  bad <- project
  bad$settings$display <- "changed after saving"
  saveRDS(bad, path)
  expect_error(.read_app_project(path), "changed since they were saved")

  # A valid enclosing seal does not excuse a stale result signature or a
  # dimensionality result from the base fit rather than the active history.
  bad <- project
  bad$results$dimensionality$reference$p_adj <- .02
  saveRDS(.seal_app_project(bad), path)
  expect_error(.read_app_project(path), "result from this fitted model")
  bad <- project
  bad$btl_steps <- list()
  saveRDS(.seal_app_project(bad), path)
  expect_error(.read_app_project(path), "result from this fitted model")

  # Current withheld analyses retain their descriptive results unchanged.
  project$results$dimensionality <- current
  saveRDS(.seal_app_project(project), path)
  expect_no_warning(restored <- .read_app_project(path))
  expect_identical(restored$results$dimensionality, current)

  project$results$dimensionality <- old
  project$schema <- 1L
  project$binding <- NULL
  saveRDS(project, path)
  expect_warning(restored <- .read_app_project(path), "schema-1")
  expect_null(restored$results$dimensionality)
  expect_identical(restored$btl_steps, project$btl_steps)
  expect_identical(restored$base_fit, project$base_fit)
  expect_no_error(.validate_app_project(restored))

  attr(project$results$dimensionality, "result_signature") <- NULL
  saveRDS(project, path)
  expect_warning(restored <- .read_app_project(path), "schema-1")
  expect_null(restored$results$dimensionality)
})

test_that("supported legacy CJ references survive project migration", {
  d <- simulate_btl(4, 12, reps_per_pair = 20, seed = 9711)
  fit <- btl(d, "object_a", "object_b", winner = "winner", judge = "judge")
  old <- btl_dimensionality(fit, reps = 20L, seed = 9712,
                            independent_comparisons = TRUE)
  old$reference$inference_available <- NULL
  attr(old, "result_signature") <- NULL
  attr(old, "result_signature") <- .fit_boot_md5(old)
  project <- btl_dimension_project(d, fit, old)
  path <- tempfile(fileext = ".rasch")
  on.exit(unlink(path), add = TRUE)
  saveRDS(project, path)
  expect_no_error(.validate_btl_dimensionality(old, fit))
  expect_no_warning(restored <- .read_app_project(path))
  expect_identical(restored$results$dimensionality, old)
  expect_identical(restored$binding, project$binding)
  expect_null(attr(restored, "rasch_project_legacy_dropped"))
})

test_that("clustered legacy CJ references without explicit independence are dropped", {
  d <- simulate_btl(4, 12, reps_per_pair = 20, seed = 9721)
  fit <- btl(d, "object_a", "object_b", winner = "winner", judge = "judge")
  expect_true(isTRUE(fit$clustered))
  old <- btl_dimensionality(fit, reps = 20L, seed = 9722,
                            independent_comparisons = TRUE)
  old$reference$independent_comparisons <- NULL
  attr(old, "result_signature") <- NULL
  attr(old, "result_signature") <- .fit_boot_md5(old)
  project <- btl_dimension_project(d, fit, old)
  path <- tempfile(fileext = ".rasch")
  on.exit(unlink(path), add = TRUE)
  saveRDS(project, path)

  expect_warning(restored <- .read_app_project(path),
                 "dimensionality reference.*unsupported comparison design")
  expect_null(restored$results$dimensionality)
  expect_identical(restored$base_fit, fit)
  expect_no_error(.validate_app_project(restored))
})

test_that("old frame references are checked against actual pair coverage", {
  d <- simulate_btl_efrm(4, 2, 5, 2, 10, 10, seed = 71)
  sets <- attr(d, "truth")$object_sets
  pair <- sets[[1L]][1:2]
  removed <- (d$object_a == pair[1L] & d$object_b == pair[2L]) |
    (d$object_a == pair[2L] & d$object_b == pair[1L])
  d <- d[!removed, ]
  fit <- btl_efrm(d, "object_a", "object_b", winner = "winner", judge = "judge",
    panels = "panel", object_sets = sets, se_method = "conditional", boot_reps = 0)
  old <- legacy_btl_dimension_reference(
    btl_dimensionality(fit, reps = 20L, seed = 9691), flag = TRUE)
  expect_true(old$leading_structured)
  expect_true(is.finite(old$reference$p_adj))
  expect_error(.validate_btl_dimensionality(old, fit), "recompute")
  project <- btl_dimension_project(d, fit, old)
  path <- tempfile(fileext = ".rasch")
  on.exit(unlink(path), add = TRUE)
  saveRDS(project, path)
  expect_warning(restored <- .read_app_project(path),
                 "dimensionality reference.*unsupported comparison design")
  expect_null(restored$results$dimensionality)
  expect_identical(restored$base_fit, fit)
  expect_no_error(.validate_app_project(restored))
})

test_that("CJ dimensionality with an obsolete residual definition is omitted", {
  d <- simulate_btl(4, 12, reps_per_pair = 20, seed = 9731)
  fit <- btl(d, "object_a", "object_b", winner = "winner", judge = "judge")
  current <- btl_dimensionality(fit, reps = 20L, seed = 9732)
  # The current result is descriptive for clustered judges; no finite
  # reference is needed to trigger this migration.
  current$residual_method <- "pooled-expected-score-1"
  for (method in list(NULL, "row-residuals-0")) {
    stale <- current
    stale$residual_method <- method
    stale$reference$mean <- stale$reference$p95 <- NA_real_
    stale$reference$p <- stale$reference$p_adj <- NA_real_
    stale$bimensions$ref_mean <- stale$bimensions$ref_p95 <- NA_real_
    attr(stale, "result_signature") <- NULL
    attr(stale, "result_signature") <- .fit_boot_md5(stale)
    project <- btl_dimension_project(d, fit, stale)
    path <- tempfile(fileext = ".rasch")
    on.exit(unlink(path), add = TRUE)
    saveRDS(project, path)

    expect_warning(restored <- .read_app_project(path),
                   "earlier residual definition.*rerun dimensionality")
    expect_null(restored$results$dimensionality)
    expect_identical(restored$data, project$data)
    expect_identical(restored$base_fit, project$base_fit)
    expect_identical(restored$btl_steps, project$btl_steps)
    expect_no_error(.validate_app_project(restored))

    # A signed result that no longer authenticates must fail rather than be
    # silently accepted as an obsolete result.
    bad <- project
    bad$results$dimensionality$residual_matrix[1L, 1L] <-
      bad$results$dimensionality$residual_matrix[1L, 1L] + 0.01
    saveRDS(.seal_app_project(bad), path)
    expect_error(.read_app_project(path), "result from this fitted model")

    # Schema-1 files have no enclosing binding, but a valid result signature
    # still authenticates the old result before it is omitted.
    project$schema <- 1L
    project$binding <- NULL
    saveRDS(project, path)
    expect_warning(restored <- .read_app_project(path),
                   "schema-1.*rerun dimensionality")
    expect_null(restored$results$dimensionality)
    expect_identical(restored$base_fit, project$base_fit)
    expect_no_error(.validate_app_project(restored))
  }
})
