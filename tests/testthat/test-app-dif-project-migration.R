test_that("saved DIF follow-ups without the normalized-data stamp are omitted", {
  set.seed(227)
  n <- 240L
  group <- factor(rep(c("A", "B"), each = n / 2L))
  theta <- rnorm(n)
  difficulty <- seq(-1.2, 1.2, length.out = 5L)
  shift <- matrix(0, n, length(difficulty))
  shift[group == "B", 3L] <- 1.5
  X <- matrix(rbinom(n * ncol(shift), 1,
                     plogis(outer(theta, difficulty, "-") - shift)),
              n, ncol(shift))
  colnames(X) <- paste0("I", seq_len(ncol(X)))
  fit <- rasch(X)
  current <- dif_anova(fit, factors = data.frame(group = group),
                       n_groups = 3L, sizes = TRUE)
  expect_true(nrow(current$posthoc) > 0L)

  legacy <- current
  legacy$followup_algorithm <- NULL
  legacy$result_signature <- NULL
  legacy$result_signature <- .fit_boot_md5(unclass(legacy))
  expect_false(.dif_followups_current(legacy))

  project <- .seal_app_project(list(
    format = "rasch-shiny-project", schema = 2L,
    data = data.frame(X, group = group), model_type = "rasch", base_fit = fit,
    rasch_steps = list(), btl_steps = list(), kept_fits = list(),
    kept_fit_code = list(), simulation = list(),
    settings = list(model_type = "rasch", dif_effects = "main",
                    dif_alpha = 0.05), resources = list(),
    results = list(dif = legacy,
                   dif_bootstrap = list(db = list(stale = TRUE)),
                   resolve = structure(list(
                     algorithm = "factor-design-resolution-2", effects = "main"),
                     class = "rasch_resolve_dif"))))
  path <- tempfile(fileext = ".rasch")
  on.exit(unlink(path), add = TRUE)
  saveRDS(project, path)

  expect_warning(restored <- .read_app_project(path),
                 "earlier factor values")
  expect_null(restored$results$dif)
  expect_null(restored$results$dif_bootstrap)
  expect_null(restored$results$resolve)
  expect_identical(restored$data, project$data)
  expect_identical(restored$base_fit, project$base_fit)
  expect_no_error(.validate_app_project(restored))
  expect_error(.save_app_project(project, path), "earlier factor values")

  other <- rasch(simulate_rasch(n, 5, seed = 1127))
  for (schema in c(1L, 2L)) {
    old <- project
    old$schema <- schema
    old <- if (schema == 2L) .seal_app_project(old) else {
      old$binding <- NULL
      old
    }
    saveRDS(old, path)
    before <- tools::md5sum(path)
    expect_warning(read <- .read_app_project(path),
                   if (schema == 1L) "schema-1" else "earlier factor values")
    expect_null(read$results$dif)
    expect_identical(read$base_fit, old$base_fit)
    expect_identical(tools::md5sum(path), before)

    # An enclosing seal cannot make an altered inner result authentic.
    bad <- old
    bad$results$dif$posthoc$estimate[1L] <-
      bad$results$dif$posthoc$estimate[1L] + 1
    if (schema == 2L) bad <- .seal_app_project(bad)
    saveRDS(bad, path)
    expect_error(.read_app_project(path), "incomplete or internally inconsistent")

    # A valid result signature for another fit must not be hidden by migration.
    bad <- old
    bad$results$dif$fit_signature <- .fit_boot_signature(other)
    bad$results$dif$result_signature <- NULL
    bad$results$dif$result_signature <- .fit_boot_md5(unclass(bad$results$dif))
    if (schema == 2L) bad <- .seal_app_project(bad)
    saveRDS(bad, path)
    expect_error(.read_app_project(path), "different fitted model")
  }

  changed <- project
  changed$data[1L, 1L] <- 1 - changed$data[1L, 1L]
  saveRDS(changed, path)
  expect_error(.read_app_project(path), "changed since they were saved")

  tampered <- project
  tampered$schema <- 1L
  tampered$binding <- NULL
  tampered$results$dif <- legacy
  tampered$results$dif$fit_signature <-
    strrep("0", nchar(legacy$fit_signature))
  tampered$results$dif$result_signature <- NULL
  tampered$results$dif$result_signature <-
    .fit_boot_md5(unclass(tampered$results$dif))
  tampered_path <- tempfile(fileext = ".rasch")
  on.exit(unlink(tampered_path), add = TRUE)
  saveRDS(tampered, tampered_path)
  expect_error(.read_app_project(tampered_path), "cannot be authenticated")
})

test_that("obsolete CJ role migration authenticates results in both schemas", {
  d <- as.data.frame(simulate_btl(4, 24, 30, seed = 816))
  j <- as.integer(sub("^J", "", d$judge))
  d$group <- ifelse(j <= 12, "A", "B")
  d$judge_perm <- paste0("J", j %% 24L + 1L)
  fit <- btl(d, "object_a", "object_b", winner = "winner", judge = "judge")
  settings <- list(model_type = "btl", bt_a = "object_a", bt_b = "object_b",
                   bt_win = "winner", bt_judge = "judge")
  attr(fit, "rasch_app_source") <- list(data = d, settings = settings,
                                       resources = list(), simulation = list())
  ids <- unique(d$judge)
  maps <- setNames(vapply(ids, function(id)
    unique(d$group[d$judge_perm == id]), character(1)), ids)
  result <- btl_dif(fit, factors = list(group = maps))
  expect_identical(unname(maps["J13"]), "A")
  expect_identical(unique(d$group[d$judge == "J13"]), "B")
  project <- .seal_app_project(list(
    format = "rasch-shiny-project", schema = 2L, data = d, model_type = "btl",
    base_fit = fit, rasch_steps = list(), btl_steps = list(), kept_fits = list(),
    kept_fit_code = list(), settings = settings, resources = list(),
    simulation = list(), results = list(btl_dif = result,
      btl_dif_meta = list(judge_col = "judge_perm"),
      dif_bootstrap = list(db = list(stale = TRUE)))))
  path <- tempfile(fileext = ".rasch")
  on.exit(unlink(path), add = TRUE)
  expect_error(.save_app_project(project, path), "judge-role provenance")
  for (schema in c(1L, 2L)) {
    old <- project
    old$schema <- schema
    if (schema == 2L) old <- .seal_app_project(old) else old$binding <- NULL
    saveRDS(old, path)
    expect_warning(read <- .read_app_project(path),
                   if (schema == 1L) "schema-1" else "Comparative Judgement DIF")
    expect_null(read$results$btl_dif)
    expect_null(read$results$btl_dif_meta)
    expect_null(read$results$dif_bootstrap)
    expect_identical(read$data, old$data)
    expect_identical(read$base_fit, old$base_fit)
    expect_no_error(.validate_app_project(read))

    bad <- old
    bad$results$btl_dif$summary$F_uniform[1L] <- 999
    if (schema == 2L) bad <- .seal_app_project(bad)
    saveRDS(bad, path)
    expect_error(.read_app_project(path), "incomplete or internally inconsistent")

    bad <- old
    bad$results$btl_dif$fit_signature <- strrep("0", nchar(result$fit_signature))
    bad$results$btl_dif$result_signature <- NULL
    bad$results$btl_dif$result_signature <-
      .fit_boot_md5(unclass(bad$results$btl_dif))
    if (schema == 2L) bad <- .seal_app_project(bad)
    saveRDS(bad, path)
    expect_error(.read_app_project(path), "different fitted model")
  }
})
