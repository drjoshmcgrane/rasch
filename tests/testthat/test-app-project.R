.app_test_path <- function() {
  app <- testthat::test_path("..", "..", "inst", "shiny", "app.R")
  if (!file.exists(app)) app <- system.file("shiny", "app.R", package = "rasch")
  app
}

test_that("the Shiny application resolves in source and installed layouts", {
  expect_true(file.exists(.app_test_path()))
})

test_that("app formula and code symbols preserve exact source names", {
  skip_if_not_installed("shiny")
  skip_if_not_installed("bslib")
  skip_if_not_installed("DT")
  skip_if_not_installed("bsicons")
  e <- new.env(parent = globalenv())
  suppressWarnings(sys.source(.app_test_path(), envir = e))

  nms <- c("if", ".2x", "a b", "a:b", "x`y")
  quoted <- e$.app_quote_name(nms)
  expect_no_error(lapply(quoted, function(z)
    parse(text = paste0("dat$", z))))

  spec <- e$.app_explanatory_interactions(c("a:b", "if", "c"))
  expect_length(unique(unname(spec$choices)), 3L)
  expect_identical(spec$map[["interaction_0001"]], c("a:b", "if"))
  f <- e$.app_explanatory_formula(
    c("a:b", "if"), list(c("a:b", "if")))
  d <- data.frame(row = seq_len(6L)); d$row <- NULL
  d[["a:b"]] <- seq_len(6L)
  d[["if"]] <- rep(0:1, 3L)
  expect_no_error(mm <- model.matrix(f, d))
  expect_setequal(colnames(mm),
                  c("(Intercept)", "`a:b`", "`if`", "`a:b`:`if`"))

  mf <- structure(list(
    items = data.frame(item = c("I1:r1", "I1:r2", "I2:r1")),
    virtual_map = data.frame(item = c("I1", "I1", "I2"))),
    class = c("rasch_mfrm", "rasch"))
  expect_identical(e$.app_dif_item_choices(mf), c("I1", "I2"))
  ordinary <- structure(list(items = data.frame(item = c("I1", "I2"))),
                        class = "rasch")
  expect_identical(e$.app_dif_item_choices(ordinary), c("I1", "I2"))

  expect_identical(e$.app_selected_row(NULL, 4L), 1L)
  expect_identical(e$.app_selected_row(3L, 4L), 3L)
  expect_identical(e$.app_selected_row(8L, 4L), 1L)
  expect_identical(e$.app_selected_row(c(2L, 3L), 4L), 1L)
  expect_error(e$.app_selected_row(1L, 0L), "positive whole row count")

  current_wald <- e$.app_wald_reference(list(t = 2.25, df = 17, z = 99))
  expect_identical(current_wald$statistic, "t")
  expect_identical(current_wald$label, "t(17)")
  expect_equal(current_wald$value, 2.25)
  expect_identical(e$.app_wald_columns(
    data.frame(t = 2.25, df = 17, z = 99)), c("t", "df"))
  legacy_wald <- e$.app_wald_reference(list(z = -1.75))
  expect_identical(legacy_wald$statistic, "z")
  expect_identical(legacy_wald$label, "z")
  expect_equal(legacy_wald$value, -1.75)
  expect_identical(legacy_wald$df, Inf)
  expect_identical(e$.app_wald_columns(data.frame(z = -1.75)), "z")

  expect_true(e$.app_equating_fit(ordinary))
  mfrm <- ordinary
  class(mfrm) <- c("rasch_mfrm", class(mfrm))
  expect_false(e$.app_equating_fit(mfrm))
  btl <- structure(list(), class = "rasch_btl")
  expect_false(e$.app_equating_fit(btl))
  explanatory <- ordinary
  class(explanatory) <- c("rasch_explanatory", class(explanatory))
  expect_false(e$.app_equating_fit(explanatory))

  position_only <- structure(list(
    dependence_data = data.frame(position = rep(1, 3))),
    class = "rasch_btl")
  ordered <- structure(list(
    dependence_data = data.frame(exposure = c(0, 1), carry_over = c(0, .5))),
    class = "rasch_btl")
  expect_false(e$.app_has_btl_history(position_only))
  expect_true(e$.app_has_btl_history(ordered))
})

test_that("sourcing the in-tree app reuses the active source namespace", {
  skip_if_not_installed("shiny")
  skip_if_not_installed("bslib")
  skip_if_not_installed("DT")
  skip_if_not_installed("bsicons")
  before <- get(".rasch_boot_apply", envir = asNamespace("rasch"),
                inherits = FALSE)
  e <- new.env(parent = globalenv())
  suppressWarnings(sys.source(.app_test_path(), envir = e))
  after <- get(".rasch_boot_apply", envir = asNamespace("rasch"),
               inherits = FALSE)
  expect_identical(after, before)
})

test_that("saved app analyses make a validated round trip", {
  set.seed(81)
  theta <- rnorm(120)
  X <- sapply(seq(-1, 1, length.out = 5), function(d)
    rbinom(length(theta), 1, plogis(theta - d)))
  colnames(X) <- paste0("I", seq_len(ncol(X)))
  fit <- rasch(X)
  project <- .seal_app_project(list(
    format = "rasch-shiny-project",
    schema = 2L,
    package_version = "test",
    created = "2026-08-16T00:00:00+1000",
    data = as.data.frame(X),
    model_type = "rasch",
    base_fit = fit,
    rasch_steps = list(),
    btl_steps = list(),
    rcode = "fit <- rasch(X)",
    kept_fits = list(reference = fit),
    kept_fit_code = list(reference = list(code = "fit <- rasch(X)",
                                          value = "fit")),
    settings = list(model_type = "rasch", item_cols = colnames(X),
                    ng_auto = TRUE, maxit = 75),
    resources = list(key = NULL),
    simulation = list(),
    results = list()
  ))
  path <- tempfile(fileext = ".rasch")
  on.exit(unlink(path))

  expect_identical(.save_app_project(project, path), path)
  restored <- .read_app_project(path)
  expect_identical(restored$format, "rasch-shiny-project")
  expect_equal(restored$data, project$data)
  expect_equal(restored$base_fit$items, fit$items)
  expect_identical(restored$rcode, project$rcode)
  expect_identical(restored$kept_fit_code, project$kept_fit_code)
  expect_identical(restored$settings, project$settings)
  expect_identical(restored$resources, project$resources)

  weights <- stats::setNames(c(2, 1, 1, 0.5, 0.5), colnames(fit$X))
  weighted <- list(
    table = weighted_person_estimates(fit, weights), weights = weights,
    by = "item", sets = NULL, filename = "weights.csv",
    fit_signature = .fit_boot_signature(fit))
  weighted$result_signature <- .fit_boot_md5(weighted)
  with_weights <- project
  with_weights$results$person_weights <- weighted
  with_weights <- .seal_app_project(with_weights)
  expect_no_error(.validate_app_project(with_weights))

  stale_weights <- with_weights
  stale_weights$results$person_weights$table$theta[1L] <-
    stale_weights$results$person_weights$table$theta[1L] + 0.1
  stale_weights$results$person_weights$result_signature <- NULL
  stale_weights$results$person_weights$result_signature <- .fit_boot_md5(
    stale_weights$results$person_weights)
  stale_weights <- .seal_app_project(stale_weights)
  expect_error(.validate_app_project(stale_weights),
               "table does not reproduce")

  # Comparison fits are active analysis objects after reopening and therefore
  # receive the same structural validation as the base and history fits.
  malformed_kept <- project
  malformed_kept$kept_fits$reference$items <- NULL
  malformed_kept <- .seal_app_project(malformed_kept)
  expect_error(.validate_app_project(malformed_kept),
               "kept fit 'reference'.*item calibration")
  unnamed_kept <- project
  names(unnamed_kept$kept_fits) <- NULL
  unnamed_kept <- .seal_app_project(unnamed_kept)
  expect_error(.validate_app_project(unnamed_kept), "kept-fit names")

  # Schema-2 projects written with the earlier text signature remain valid
  # across the LF/CRLF boundary.
  legacy <- project
  unsigned <- legacy
  unsigned$binding <- NULL
  legacy$binding <- .fit_boot_md5_legacy_candidates(unsigned)[2L]
  expect_no_error(.validate_app_project(legacy))
})

test_that("saved dimensionality results are authenticated against the active fit", {
  set.seed(811)
  X <- matrix(rbinom(600, 1, .5), 120, 5,
              dimnames = list(NULL, paste0("I", 1:5)))
  fit <- rasch(X)
  dimensionality <- .scree_analysis(fit, n_components = 3,
                                    parallel = FALSE)
  project <- .seal_app_project(list(
    format = "rasch-shiny-project", schema = 2L,
    data = as.data.frame(X), model_type = "rasch", base_fit = fit,
    rasch_steps = list(), btl_steps = list(), settings = list(),
    results = list(dimensionality = dimensionality)))
  expect_no_error(.validate_app_project(project))

  project$results$dimensionality$eigenvalue[1] <-
    project$results$dimensionality$eigenvalue[1] + 1
  expect_error(.validate_app_project(project),
               "saved dimensionality analysis")
})

test_that("saved person-subset tests are authenticated against the active fit", {
  set.seed(8111)
  X <- matrix(rbinom(720, 1, .5), 120, 6,
              dimnames = list(NULL, paste0("I", 1:6)))
  fit <- rasch(X)
  subtest <- dimensionality_test(
    fit, items_positive = colnames(X)[1:3],
    items_negative = colnames(X)[4:6], min_score_points = 2)
  project <- .seal_app_project(list(
    format = "rasch-shiny-project", schema = 2L,
    data = as.data.frame(X), model_type = "rasch", base_fit = fit,
    rasch_steps = list(), btl_steps = list(), settings = list(),
    results = list(subtest = subtest)))
  expect_no_error(.validate_app_project(project))

  project$results$subtest$prop_significant <-
    project$results$subtest$prop_significant + .01
  expect_error(.validate_app_project(project),
               "saved person-subset dimensionality test")
})

test_that("schema-2 projects drop results made with the earlier maxT adjustment", {
  set.seed(812)
  X <- matrix(rbinom(600, 1, .5), 120, 5,
              dimnames = list(NULL, paste0("I", 1:5)))
  fit <- rasch(X)
  old_bs <- suppressWarnings(fit_bootstrap(fit, B = 3, workers = 1,
                                            seed = 812))
  old_bs$algorithm <- NULL
  unsigned <- old_bs
  unsigned$result_signature <- NULL
  old_bs$result_signature <- .fit_boot_md5(unsigned)
  project <- .seal_app_project(list(
    format = "rasch-shiny-project", schema = 2L,
    data = as.data.frame(X), model_type = "rasch", base_fit = fit,
    rasch_steps = list(), btl_steps = list(), settings = list(),
    results = list(bootstrap = list(bs = old_bs, B = 3L, seed = 812L,
                                    kind = "rasch"))))
  path <- tempfile(fileext = ".rasch")
  on.exit(unlink(path), add = TRUE)
  saveRDS(project, path)
  expect_warning(restored <- .read_app_project(path), "earlier maxT")
  expect_null(restored$results$bootstrap)
  expect_match(attr(restored, "rasch_project_legacy_dropped"),
               "superseded maxT")
  expect_no_error(.validate_app_project(restored))

  # Early schema-2 files used the legacy text signature for the enclosing
  # project as well as for the bootstrap result. They must reach the same
  # migration rather than being refused as altered files.
  legacy <- project
  legacy$binding <- NULL
  legacy$binding <- .fit_boot_md5_legacy_candidates(legacy)[2L]
  saveRDS(legacy, path)
  expect_warning(restored_legacy <- .read_app_project(path), "earlier maxT")
  expect_null(restored_legacy$results$bootstrap)
  expect_match(attr(restored_legacy, "rasch_project_legacy_dropped"),
               "superseded maxT")
  expect_no_error(.validate_app_project(restored_legacy))
})

test_that("schema-2 projects omit superseded frame-invariance inference", {
  set.seed(813)
  X <- matrix(rbinom(600, 1, .5), 120, 5,
              dimnames = list(NULL, paste0("I", 1:5)))
  fit <- rasch(X)
  # This result has complete accounting and resampling provenance but predates
  # the exact-comparison-family algorithm identifier. The project seal is
  # valid, so the reader should migrate the derived result rather than reject
  # the whole analysis file.
  old_invariance <- structure(list(
    boot_reps = 100L, boot_reps_used = 85L,
    boot_reps_nonconverged = 10L, boot_reps_errors = 5L,
    boot_minimum_usable = 51L, family_n = 1L,
    bootstrap_stratified = TRUE,
    fit_signature = .fit_boot_signature(fit),
    result_signature = "legacy-result"),
    class = c("rasch_frame_invariance", "list"))
  project <- .seal_app_project(list(
    format = "rasch-shiny-project", schema = 2L,
    data = as.data.frame(X), model_type = "rasch", base_fit = fit,
    rasch_steps = list(), btl_steps = list(), settings = list(),
    results = list(frame_invariance = old_invariance)))
  path <- tempfile(fileext = ".rasch")
  on.exit(unlink(path), add = TRUE)
  saveRDS(project, path)

  expect_warning(restored <- .read_app_project(path),
                 "frame-invariance analysis used earlier")
  expect_null(restored$results$frame_invariance)
  expect_match(attr(restored, "rasch_project_legacy_dropped"),
               "frame-invariance")
  expect_no_error(.validate_app_project(restored))
})

test_that("invalid app analysis files are refused", {
  bad <- tempfile(fileext = ".rasch")
  on.exit(unlink(bad))
  saveRDS(list(format = "other", schema = 1L), bad)
  expect_error(.read_app_project(bad), "not a rasch analysis file")

  set.seed(2)
  X <- matrix(rbinom(300, 1, .5), 60, 5,
              dimnames = list(NULL, paste0("I", 1:5)))
  fit <- rasch(X)
  future <- list(format = "rasch-shiny-project", schema = 3L,
                 data = as.data.frame(X), base_fit = fit)
  saveRDS(future, bad)
  expect_error(.read_app_project(bad), "unsupported")

  shaped_schema <- .seal_app_project(list(
    format = "rasch-shiny-project", schema = 2L,
    data = as.data.frame(X), model_type = "rasch", base_fit = fit,
    settings = list(), results = list()))
  shaped_schema$schema <- matrix(2L)
  expect_error(.validate_app_project(shaped_schema), "unsupported")

  malformed <- .seal_app_project(list(
    format = "rasch-shiny-project", schema = 2L,
    data = as.data.frame(X), base_fit = fit, settings = "not a list"))
  saveRDS(malformed, bad)
  expect_error(.read_app_project(bad), "invalid settings")

  malformed_results <- .seal_app_project(list(
    format = "rasch-shiny-project", schema = 2L,
    data = as.data.frame(X), base_fit = fit, results = "not a list"))
  saveRDS(malformed_results, bad)
  expect_error(.read_app_project(bad), "invalid results")

  legacy <- list(format = "rasch-shiny-project", schema = 1L,
                 data = as.data.frame(X), model_type = "rasch",
                 base_fit = fit, rasch_steps = list(), btl_steps = list(),
                 results = list())
  saveRDS(legacy, bad)
  expect_warning(upgraded <- .read_app_project(bad), "schema-1")
  expect_identical(upgraded$schema, 2L)
  expect_true(isTRUE(attr(upgraded, "rasch_project_legacy")))
  expect_no_error(.validate_app_project(upgraded))

  # Schema 1 could save a bootstrap before result fingerprints existed. The
  # analysis remains reopenable, but that unverifiable result is not restored.
  old_bs <- suppressWarnings(fit_bootstrap(fit, B = 3, seed = 17))
  old_bs$result_signature <- NULL
  legacy$results <- list(bootstrap = list(bs = old_bs, B = 3L,
                                          seed = 17L, kind = "rasch"))
  saveRDS(legacy, bad)
  expect_warning(upgraded_bs <- .read_app_project(bad),
                 "omitted: fit bootstrap")
  expect_null(upgraded_bs$results$bootstrap)
  expect_identical(attr(upgraded_bs, "rasch_project_legacy_dropped"),
                   "fit bootstrap")
  expect_s3_class(upgraded_bs$base_fit, "rasch")
  expect_no_error(.validate_app_project(upgraded_bs))

  empty_fit <- structure(list(), class = "rasch")
  malformed_fit <- .seal_app_project(list(
    format = "rasch-shiny-project", schema = 2L,
    data = as.data.frame(X), model_type = "rasch", base_fit = empty_fit,
    rasch_steps = list(), btl_steps = list(), results = list()))
  expect_error(.validate_app_project(malformed_fit), "saved base fit")

  bound <- .seal_app_project(list(
    format = "rasch-shiny-project", schema = 2L,
    data = as.data.frame(X), model_type = "rasch", base_fit = fit,
    rasch_steps = list(), btl_steps = list(), results = list()))
  expect_no_error(.validate_app_project(bound))
  matrix_bound <- bound
  matrix_bound$data <- X
  matrix_bound <- .seal_app_project(matrix_bound)
  expect_no_error(.validate_app_project(matrix_bound))
  blank_named <- bound
  names(blank_named$data)[1] <- "   "
  blank_named <- .seal_app_project(blank_named)
  expect_error(.validate_app_project(blank_named), "column names")
  bound$data <- data.frame(unrelated = seq_len(nrow(X)))
  expect_error(.validate_app_project(bound), "changed since they were saved")
})

test_that("saved bootstrap results must belong to the active saved fit", {
  set.seed(83)
  X1 <- matrix(rbinom(600, 1, 0.5), 120, 5,
               dimnames = list(NULL, paste0("I", 1:5)))
  X2 <- X1
  X2[1:10, 1] <- 1L - X2[1:10, 1]
  fit1 <- rasch(X1)
  fit2 <- rasch(X2)
  bs <- suppressWarnings(fit_bootstrap(
    fit1, B = 1L, workers = 1L, seed = 83L))
  project <- .seal_app_project(list(
    format = "rasch-shiny-project", schema = 2L,
    data = as.data.frame(X1), model_type = "rasch", base_fit = fit1,
    rasch_steps = list(), btl_steps = list(), results = list(
      bootstrap = list(bs = bs, B = 1L, seed = 1L, kind = "rasch"))))
  expect_no_error(.validate_app_project(project))

  history_entry <- function(type, fit) list(
    type = type, label = type, fit = fit, details = list(), code = NULL,
    created = "2026-08-16T00:00:00+1000")
  project$rasch_steps <- list(history_entry("test", fit2))
  expect_error(.validate_app_project(project),
               "does not belong to the active fit")
  project$rasch_steps <- list(
    history_entry("corrupt", "not a fitted model"),
    history_entry("test", fit1))
  expect_error(.validate_app_project(project), "fitted-model history")
  wrong_family <- fit1
  class(wrong_family) <- c("rasch_efrm", class(wrong_family))
  project$rasch_steps <- list(history_entry("wrong", wrong_family))
  expect_error(.validate_app_project(project), "fitted-model history")
  project$rasch_steps <- list()
  project$results$bootstrap <- list(bs = list())
  expect_error(.validate_app_project(project), "saved bootstrap result")

  project$results <- list()
  project$base_fit <- wrong_family
  expect_error(.validate_app_project(project), "declares model type")
})

test_that("opening a project retains results tied to its active fit", {
  skip_on_cran()
  skip_if_not_installed("shiny")
  skip_if_not_installed("bslib")
  skip_if_not_installed("DT")
  skip_if_not_installed("bsicons")

  app <- .app_test_path()
  e <- new.env(parent = globalenv())
  suppressWarnings(sys.source(app, envir = e))
  set.seed(82)
  X <- matrix(rbinom(500, 1, .5), 100, 5,
              dimnames = list(NULL, paste0("I", 1:5)))
  fit <- rasch(X)
  project <- .seal_app_project(list(
    format = "rasch-shiny-project", schema = 2L,
    data = as.data.frame(X), model_type = "rasch", base_fit = fit,
    rasch_steps = list(), btl_steps = list(), rcode = "fit <- rasch(dat)",
    kept_fits = list(), kept_fit_code = list(),
    settings = list(model_type = "rasch", id_col = "(none)",
                    factor_cols = character(0), item_cols = colnames(X),
                    thr_structure = "pcm", ng_auto = FALSE, ng = "6",
                    maxit = 75, tol = 1e-7,
                    eq_source = "csv", eq_shift = "mean",
                    eq_csv_independent = TRUE,
                    exp_type_1 = "ordinal", exp_order_1 = "low,high",
                    exp_type_2 = "categorical", exp_ref_2 = "B"),
    resources = list(
      anchors = NULL, key = NULL,
      eq_reference = fit$items[, c("item", "location", "se", "max")],
      predictors = data.frame(item = colnames(X),
        stage = rep(c("low", "high"), length.out = ncol(X)),
        family = rep(c("A", "B"), length.out = ncol(X)))),
    simulation = list(),
    results = list(lr = list(marker = 17L))))
  path <- tempfile(fileext = ".rasch")
  on.exit(unlink(path), add = TRUE)
  .save_app_project(project, path)

  shiny::testServer(e$server, {
    session$setInputs(project_file = list(
      datapath = path, name = "test.rasch", size = file.info(path)$size,
      type = "application/octet-stream"))
    session$flushReact()
    expect_identical(lr_res()$marker, 17L)
    expect_s3_class(fit(), "rasch")
    expect_equal(raw_data(), as.data.frame(X))
    expect_no_error(src <- data_source_code())
    expect_match(src, 'readRDS\\("test[.]rasch"\\)')
    expect_match(src, "dat <- project[$]data")
    expect_identical(restored_project_settings(), project$settings)
    expect_identical(restored_project_resources(), project$resources)
    expect_identical(eq_ref(), project$resources$eq_reference)
    expect_s3_class(eq_res(), "rasch_equate")
    p <- exp_predictors_raw()
    expect_identical(exp_predictor_type(p, "stage"), "ordinal")
    expect_identical(exp_level_order(p, "stage"), c("low", "high"))
    expect_identical(exp_predictor_type(p, "family"), "categorical")
    expect_identical(exp_category_levels(p, "family"), c("B", "A"))
    html <- paste(as.character(output$data_main), collapse = " ")
    expect_false(grepl("Welcome to rasch", html, fixed = TRUE))
    expect_match(html, "Data preview", fixed = TRUE)

    # Saving again captures the current run-defining controls rather than
    # relying on column-name guesses at the next opening.
    session$setInputs(id_col = "(none)", item_cols = colnames(X),
                      factor_cols = character(0), thr_structure = "pcm",
                      ng_auto = FALSE, ng = "6", maxit = 75, tol = 1e-7)
    saved <- project_state()
    expect_identical(saved$settings$item_cols, colnames(X))
    expect_identical(saved$settings$ng, "6")
    expect_identical(saved$settings$maxit, 75)
  })
})

test_that("CJ results and background work remain tied to their launching fit", {
  skip_on_cran()
  skip_if_not_installed("shiny")
  skip_if_not_installed("bslib")
  skip_if_not_installed("DT")
  skip_if_not_installed("bsicons")

  e <- new.env(parent = globalenv())
  suppressWarnings(sys.source(.app_test_path(), envir = e))
  d1 <- simulate_btl(5, 12, 8, seed = 91)
  d2 <- simulate_btl(5, 12, 8, seed = 92)
  bt1 <- btl(d1, "object_a", "object_b", "winner", judge = "judge")
  attr(bt1, "rasch_app_source") <- list(
    data = as.data.frame(d1),
    settings = list(model_type = "btl", bt_a = "object_a",
                    bt_b = "object_b", bt_win = "winner",
                    bt_judge = "judge"),
    resources = list(), simulation = list())
  bt2 <- btl(d2, "object_a", "object_b", "winner", judge = "judge")
  judges <- unique(bt1$comparisons$judge)
  group <- setNames(rep(c("A", "B"), length.out = length(judges)), judges)
  bdif1 <- btl_dif(bt1, group, objects = "O3", min_n = 2)
  project <- .seal_app_project(list(
    format = "rasch-shiny-project", schema = 2L,
    data = as.data.frame(d1), model_type = "btl", base_fit = bt1,
    rasch_steps = list(), btl_steps = list(), rcode = "bt <- btl(dat)",
    kept_fits = list(), kept_fit_code = list(),
    settings = list(model_type = "btl", bt_a = "object_a",
                    bt_b = "object_b", bt_win = "winner",
                    bt_judge = "judge"),
    resources = list(), simulation = list(),
    results = list(btl_dif = bdif1,
                   btl_dif_meta = list(judge_col = "judge",
                                       fitted_judge_col = "judge"))))
  no_source <- project
  no_source$results$btl_dif_meta$fitted_judge_col <- NULL
  no_source <- .seal_app_project(no_source)
  expect_error(.validate_app_project(no_source),
               "lacks authenticated fitted judge-role provenance")
  path <- tempfile(fileext = ".rasch")
  on.exit(unlink(path), add = TRUE)
  .save_app_project(project, path)

  legacy_path <- tempfile(fileext = ".rasch")
  on.exit(unlink(legacy_path), add = TRUE)
  legacy <- project
  legacy$schema <- 1L
  legacy$binding <- NULL
  legacy$results$btl_dif$result_signature <- NULL
  saveRDS(legacy, legacy_path)
  expect_warning(upgraded <- .read_app_project(legacy_path),
                 "Comparative Judgement DIF")
  expect_null(upgraded$results$btl_dif)
  expect_null(upgraded$results$btl_dif_meta)
  expect_s3_class(upgraded$base_fit, "rasch_btl")

  shiny::testServer(e$server, {
    # A successful replacement owns none of the preceding fit's requested
    # results, even when the replacement is another CJ calibration.
    bdif_res(list(marker = 12L))
    btlef_res(list(marker = 11L))
    person_weight_state(list(marker = 10L))
    complete_fit(bt2, "bt <- btl(dat)", character(0), "dat <- replacement")
    session$flushReact()
    expect_null(bdif_res())
    expect_null(btlef_res())
    expect_null(person_weight_state())

    # The fit observer must not erase results serialised with a reopened fit.
    session$setInputs(project_file = list(
      datapath = path, name = "cj.rasch", size = file.info(path)$size,
      type = "application/octet-stream"))
    session$flushReact()
    expect_identical(bdif_res()$result_signature, bdif1$result_signature)
    expect_identical(bdif_meta()$judge_col, "judge")
    expect_no_error(rasch:::.validate_btl_dif_result(bdif_res(), bfit()))
    expect_no_error(.validate_app_project(project_state()))
    expect_identical(btl_fit()$objects$object, bt1$objects$object)

    # Opening a project cancels either kind of worker before restoring state.
    fake_process <- function() {
      p <- new.env(parent = emptyenv())
      p$alive <- TRUE
      p$is_alive <- function() p$alive
      p$kill_tree <- function() p$alive <- FALSE
      p$kill <- function() p$alive <- FALSE
      p
    }
    ep <- fake_process(); bp <- fake_process()
    fake_job <- function(process) list(
      process = process, progress = NULL,
      progress_file = tempfile(), log_file = tempfile())
    efrm_job(fake_job(ep)); btlef_job(fake_job(bp))
    session$setInputs(project_file = list(
      datapath = path, name = "cj-again.rasch", size = file.info(path)$size,
      type = "application/octet-stream"))
    session$flushReact()
    expect_false(ep$alive)
    expect_false(bp$alive)
    expect_null(efrm_job())
    expect_null(btlef_job())
    expect_identical(bdif_res()$result_signature, bdif1$result_signature)
    expect_identical(bdif_meta()$judge_col, "judge")

    # A worker result is current only while its context, data and CJ base fit
    # are the same as at launch.
    st <- list(context = analysis_context(), data = raw_data(),
               check_btl_fit = TRUE, base_fit = btl_fit())
    expect_true(background_job_is_current(st))
    advance_analysis_context()
    expect_false(background_job_is_current(st))

    st$context <- analysis_context()
    st$data <- raw_data()
    expect_true(background_job_is_current(st))
    btl_fit(bt2)
    expect_false(background_job_is_current(st))

    btl_fit(bt1)
    st$base_fit <- bt1
    sim_data(as.data.frame(d2))
    expect_false(background_job_is_current(st))

    # Both completion observers apply the guard before storing their result.
    completed_process <- function(value) {
      p <- new.env(parent = emptyenv())
      p$is_alive <- function() FALSE
      p$get_result <- function() list(value = value, warnings = character(0))
      p
    }
    common <- list(progress = NULL, progress_file = tempfile(),
                   log_file = tempfile(), context = analysis_context(),
                   data = raw_data())
    ej <- c(list(process = completed_process(bt2), code_call = NULL,
                 src_line = "dat <- old", se_method = "hybrid",
                 simulation_stamp = NULL), common)
    advance_analysis_context()
    efrm_job(ej)
    session$flushReact()
    expect_null(efrm_job())
    expect_identical(btl_fit()$objects$location, bt1$objects$location)

    btlef_res(NULL); clear_btl_analysis_steps()
    bj <- c(list(process = completed_process(list(marker = 99L)),
                 check_btl_fit = TRUE, base_fit = bt1), common)
    btlef_job(bj)
    session$flushReact()
    expect_null(btlef_job())
    expect_null(btlef_res())
    expect_length(btl_analysis_steps(), 0L)
  })
})

test_that("schema-2 CJ DIF without fitted judge provenance is omitted", {
  d <- as.data.frame(simulate_btl(5, 16, 20, seed = 919))
  j <- as.integer(sub("^J", "", d$judge))
  d$judge_perm <- paste0("J", ifelse(j %% 16 == 0, 1, j %% 16 + 1))
  d$group <- ifelse(j <= 8, "A", "B")
  bt <- btl(d, "object_a", "object_b", winner = "winner", judge = "judge")
  ids <- unique(as.character(d$judge))
  perm_group <- setNames(vapply(ids, function(id)
    unique(as.character(d$group[d$judge_perm == id]))[1L], character(1)), ids)
  old_dif <- btl_dif(bt, factors = list(group = perm_group), min_n = 2)
  project <- .seal_app_project(list(
    format = "rasch-shiny-project", schema = 2L,
    data = d, model_type = "btl", base_fit = bt,
    rasch_steps = list(), btl_steps = list(), kept_fits = list(),
    kept_fit_code = list(), settings = list(model_type = "btl",
      bt_a = "object_a", bt_b = "object_b", bt_win = "winner",
      bt_judge = "judge_perm"), resources = list(), simulation = list(),
    results = list(btl_dif = old_dif,
      btl_dif_meta = list(judge_col = "judge_perm"),
      dif_bootstrap = list(db = list(stale = TRUE)))))
  path <- tempfile(fileext = ".rasch")
  on.exit(unlink(path), add = TRUE)
  saveRDS(project, path)
  expect_warning(restored <- .read_app_project(path),
                 "Comparative Judgement DIF")
  expect_null(restored$results$btl_dif)
  expect_null(restored$results$btl_dif_meta)
  expect_null(restored$results$dif_bootstrap)
  expect_identical(restored$data, project$data)
  expect_identical(restored$base_fit, project$base_fit)
  expect_no_error(.validate_app_project(restored))
})

test_that("a failed replacement leaves the current app analysis intact", {
  skip_on_cran()
  skip_if_not_installed("shiny")
  skip_if_not_installed("bslib")
  skip_if_not_installed("DT")
  skip_if_not_installed("bsicons")

  app <- .app_test_path()
  e <- new.env(parent = globalenv())
  suppressWarnings(sys.source(app, envir = e))
  set.seed(83)
  X <- matrix(rbinom(500, 1, .5), 100, 5,
              dimnames = list(NULL, paste0("I", 1:5)))
  current <- rasch(X)

  shiny::testServer(e$server, {
    fit_val(current)
    push_analysis_step("test", "Existing change", current)
    complete_fit(simpleError("replacement failed"), NULL, character(0),
                 "dat <- existing")
    expect_identical(fit_val(), current)
    expect_length(analysis_steps(), 1L)
  })
})

test_that("the app calculates and restores external person weights", {
  skip_on_cran()
  skip_if_not_installed("shiny")
  skip_if_not_installed("bslib")
  skip_if_not_installed("DT")
  skip_if_not_installed("bsicons")

  e <- new.env(parent = globalenv())
  suppressWarnings(sys.source(.app_test_path(), envir = e))
  current <- rasch(simulate_rasch(80, 5, seed = 921))
  w <- data.frame(item = colnames(current$X),
                  weight = c(2, 2, 1, 1, 0.5))
  path <- tempfile(fileext = ".csv")
  on.exit(unlink(path), add = TRUE)
  write.csv(w, path, row.names = FALSE)

  shiny::testServer(e$server, {
    fit_val(current)
    session$setInputs(
      person_weight_level = "item",
      person_weights_file = list(datapath = path, name = "weights.csv",
        size = file.info(path)$size, type = "text/csv"),
      person_weights_go = 1)
    session$flushReact()
    z <- person_weight_state()
    expect_equal(nrow(z$table), nrow(current$X))
    expect_identical(z$by, "item")
    expect_true(.fit_boot_signature_matches(z$fit_signature, current))
    unsigned <- z; unsigned$result_signature <- NULL
    expect_true(.fit_boot_hash_matches(z$result_signature, unsigned))
    expect_match(person_weight_code(), "weighted_person_estimates")
  })
})

test_that("the app simulator preserves explanatory metadata and frame calls", {
  skip_on_cran()
  skip_if_not_installed("shiny")
  skip_if_not_installed("bslib")
  skip_if_not_installed("DT")
  skip_if_not_installed("bsicons")

  e <- new.env(parent = globalenv())
  suppressWarnings(sys.source(.app_test_path(), envir = e))
  shiny::testServer(e$server, {
    session$setInputs(sim_layout = "rasch_exp", sim_seed = 17.5,
                      sim_go = 1)
    session$flushReact()
    expect_null(sim_data())
    expect_null(sim_code_val())

    session$setInputs(
      sim_layout = "rasch_exp", sim_seed = 17,
      sr_persons = 100, sr_items = 6, sr_model = "dichotomous", sr_cats = 4,
      sr_mean = 0, sr_sd = 1, sr_dist = "normal", sr_diff = c(-2, 2),
      sr_over = 0, sr_under = 0, sr_guess = FALSE, sr_2d = FALSE,
      sr_rho = 0.3, sr_dep = FALSE, sr_dif = FALSE, sr_difmag = 1,
      sr_style = FALSE, sr_styletype = "extreme", sr_speeded = 0,
      sr_careless = 0, sr_missing = 0, sx_cont = 1, sx_cat = 0.5,
      sx_interaction = TRUE, sx_int = 0.4, sx_depart = 0.7, sim_go = 2)
    session$flushReact()
    expect_equal(nrow(sim_predictors_val()), 6L)
    expect_identical(sim_interactions_val(), "exposure:type")
    regenerated <- eval(parse(text = sim_code_val()))
    expect_s3_class(regenerated, "rasch_sim")
    expect_equal(nrow(attr(regenerated, "predictors")), 6L)
    expect_identical(attr(regenerated, "truth")$explanatory_formula,
                     "~ exposure + type + exposure:type")
    expect_true(any(grepl("fixed explanatory departure",
                          attr(regenerated, "truth")$planted)))
    expect_identical(attr(regenerated, "truth")$departure_types,
                     "a fixed explanatory item departure")
    expect_true(attr(regenerated, "truth")$explanatory_departure$target %in%
                  attr(regenerated, "predictors")$item)

    # Item-specific slope departures need at least one reference item. If all
    # items carry the same non-unit slope, that is only a change of logit unit
    # and must not be presented as planted item misfit.
    previous_generation <- sim_gen()
    previous_code <- sim_code_val()
    session$setInputs(sr_over = 5, sr_under = 1, sr_items = 6,
                      sim_go = 2.5)
    session$flushReact()
    expect_identical(sim_gen(), previous_generation)
    expect_identical(sim_code_val(), previous_code)
    session$setInputs(sr_over = 0, sr_under = 0)

    bundle <- tempfile(fileext = ".zip")
    unpacked <- tempfile("rasch-simulation-test-")
    dir.create(unpacked)
    on.exit(unlink(c(bundle, unpacked), recursive = TRUE, force = TRUE),
            add = TRUE)
    write_sim_bundle(bundle)
    bundle_files <- utils::unzip(bundle, list = TRUE)$Name
    expect_setequal(bundle_files,
      c("data.csv", "simulation.R", "truth.rds", "predictors.csv",
        "README.txt"))
    utils::unzip(bundle, exdir = unpacked)
    expect_equal(nrow(read.csv(file.path(unpacked, "data.csv"))), 100L)
    expect_identical(readRDS(file.path(unpacked, "truth.rds"))$layout,
                     "rasch")
    expect_equal(nrow(read.csv(file.path(unpacked, "predictors.csv"))), 6L)

    # Erratic raters are generated from R1 upwards. The app must place a
    # simultaneous item-by-rater bias on a rater whose model-based responses
    # have not been replaced by random ratings.
    session$setInputs(
      sim_layout = "mfrm", sim_seed = 18,
      sm_persons = 30, sm_items = 4, sm_raters = 6, sm_cats = 4,
      sm_thsd = 1.2, sm_itemsd = 1, sm_sev = 0.6,
      sm_erratic = 0.4, sm_halo = 0, sm_int = TRUE, sim_go = 3)
    session$flushReact()
    regenerated_mfrm <- eval(parse(text = sim_code_val()))
    expect_identical(attr(regenerated_mfrm, "truth")$erratic,
                     c("R1", "R2"))
    expect_match(sim_code_val(), 'rater = "R3"', fixed = TRUE)
    expect_true(any(grepl("R3 on I2", attr(regenerated_mfrm, "truth")$planted,
                          fixed = TRUE)))

    # Fail closed if a future control range permits every rater to be made
    # erratic. There is then no model-based rater on which an interaction can
    # be planted, and the successful simulation already loaded must survive.
    previous_generation <- sim_gen()
    previous_code <- sim_code_val()
    session$setInputs(sm_erratic = 1, sm_int = TRUE, sim_go = 4)
    session$flushReact()
    expect_identical(sim_gen(), previous_generation)
    expect_identical(sim_code_val(), previous_code)

    session$setInputs(
      sim_layout = "btl_efrm", sim_seed = 18,
      sbf_objects = 4, sbf_sets = 2, sbf_judges = 4, sbf_panels = 2,
      sbf_within = 5, sbf_cross = 5, sbf_objsd = 1,
      sbf_setratio = 1.3, sbf_panelratio = 1.2, sbf_origin = 0.5,
      sbf_erratic = 0.25, sim_go = 5)
    session$flushReact()
    regenerated <- eval(parse(text = sim_code_val()))
    tr <- attr(regenerated, "truth")
    expect_identical(tr$layout, "btl_efrm")
    expect_length(tr$erratic, 2L)

    # Recovery must use the active frame-adjusted comparison fit, not the
    # equal-unit base fit. A truth-valued stand-in isolates that app routing
    # from the estimator, which is covered by test-simulate.R. Retain the
    # exact comparison and panel allocation so sim_recovery() can enforce its
    # fitted-data identity check rather than treating this as an unverified
    # reporting object.
    comparisons <- data.frame(
      object_a = regenerated$object_a,
      object_b = regenerated$object_b,
      response = as.numeric(regenerated$winner == regenerated$object_a),
      weight = 1,
      judge = regenerated$judge,
      panel = regenerated$panel,
      stringsAsFactors = FALSE)
    framed <- structure(list(
      converged = TRUE,
      m = 1L,
      comparisons = comparisons,
      objects = data.frame(object = names(tr$v), v = unname(tr$v)),
      phi_table = data.frame(panel = names(tr$phi), phi = unname(tr$phi)),
      alpha_table = data.frame(set = names(tr$alpha), alpha = unname(tr$alpha)),
      kappa_table = data.frame(set = names(tr$kappa), kappa = unname(tr$kappa)),
      set_of = stats::setNames(paste0("set", tr$set_of), names(tr$set_of))),
      class = c("rasch_btl_efrm", "rasch_btl"))
    btl_fit(framed)
    fitted_sim_gen(sim_gen())
    rr <- sim_recovery_val()
    expect_s3_class(rr, "rasch_recovery")
    expect_setequal(rr$summary$parameter,
      c("object location", "panel unit (log)", "set unit (log)", "set origin"))
  })
})

test_that("the app exposes reproducible WrightMap panel controls", {
  skip_on_cran()
  skip_if_not_installed("shiny")
  skip_if_not_installed("bslib")
  skip_if_not_installed("DT")
  skip_if_not_installed("bsicons")
  skip_if_not_installed("WrightMap")
  skip_if_not("item.groups" %in% names(formals(WrightMap::wrightMap)),
              "WrightMap item panels require version 1.5")

  app <- .app_test_path()
  e <- new.env(parent = globalenv())
  suppressWarnings(sys.source(app, envir = e))

  d <- simulate_efrm(n_per_group = 55, items_per_set = 4, n_sets = 2,
                     n_groups = 2, n_categories = 2, seed = 73)
  truth <- attr(d, "truth")
  current <- rasch_efrm(d, item_sets = truth$item_sets, groups = "group",
                        id = "id", boot_reps = 0)

  shiny::testServer(e$server, {
    fit_val(current)
    session$flushReact()
    session$setInputs(wright_renderer = "wrightmap",
                      wright_type = "thresholds",
                      wright_person_panels = "groups",
                      wright_item_panels = "sets_groups",
                      wright_person_style = "histogram",
                      tg_bins = 35,
                      tg_axis_mode = "standard")
    session$flushReact()
    expect_match(output$wright_code, "wright_map\\(fit")
    expect_match(output$wright_code, 'person_panels = "groups"', fixed = TRUE)
    expect_match(output$wright_code,
                 'item_panels = c("sets", "groups")', fixed = TRUE)
    expect_type(output$wright, "list")
  })
})

test_that("uploaded frame maps refuse conflicting or unknown assignments", {
  skip_on_cran()
  skip_if_not_installed("shiny")
  skip_if_not_installed("bslib")
  skip_if_not_installed("DT")
  skip_if_not_installed("bsicons")

  app <- .app_test_path()
  e <- new.env(parent = globalenv())
  suppressWarnings(sys.source(app, envir = e))
  duplicate <- tempfile(fileext = ".csv")
  unknown <- tempfile(fileext = ".csv")
  blank <- tempfile(fileext = ".csv")
  spaced <- tempfile(fileext = ".csv")
  on.exit(unlink(c(duplicate, unknown, blank, spaced)), add = TRUE)
  write.csv(data.frame(item = c("I1", "I1"), set = c("A", "B")),
            duplicate, row.names = FALSE)
  write.csv(data.frame(item = c("I1", "wrong"), set = c("A", "B")),
            unknown, row.names = FALSE)
  write.csv(data.frame(item = c("I1", "I2"), set = c("A", " ")),
            blank, row.names = FALSE)
  write.csv(data.frame(item = c(" I1 ", "I2"), set = c(" A ", "B")),
            spaced, row.names = FALSE)

  shiny::testServer(e$server, {
    expect_error(read_frame_map(list(datapath = duplicate), "item",
                                c("I1", "I2"), "item"),
                 "assigns item values more than once")
    expect_error(read_frame_map(list(datapath = unknown), "item",
                                c("I1", "I2"), "item"),
                 "unknown item values")
    expect_error(read_frame_map(list(datapath = blank), "item",
                                c("I1", "I2"), "item"),
                 "non-blank set")
    expect_identical(
      read_frame_map(list(datapath = spaced), "item", c("I1", "I2"),
                     "item"),
      c(I1 = "A", I2 = "B"))
  })
})

test_that("supplied app anchors fail closed", {
  skip_on_cran()
  skip_if_not_installed("shiny")
  skip_if_not_installed("bslib")
  skip_if_not_installed("DT")
  skip_if_not_installed("bsicons")

  app <- .app_test_path()
  e <- new.env(parent = globalenv())
  suppressWarnings(sys.source(app, envir = e))
  bad_rasch <- tempfile(fileext = ".csv")
  bad_btl <- tempfile(fileext = ".csv")
  on.exit(unlink(c(bad_rasch, bad_btl)), add = TRUE)
  write.csv(data.frame(item = "I1", value = 0), bad_rasch,
            row.names = FALSE)
  write.csv(data.frame(object = "A", value = 0), bad_btl,
            row.names = FALSE)

  shiny::testServer(e$server, {
    session$setInputs(anchor_file = list(
      datapath = bad_rasch, name = "bad.csv", size = file.info(bad_rasch)$size,
      type = "text/csv"))
    expect_error(anchors_in(), "needs columns item, k, tau")
    session$setInputs(bt_anchor_file = list(
      datapath = bad_btl, name = "bad.csv", size = file.info(bad_btl)$size,
      type = "text/csv"))
    expect_error(bt_anchors_in(), "needs columns object, location")
  })
})

test_that("app BTL DIF refuses a factor that varies within judge", {
  skip_on_cran()
  skip_if_not_installed("shiny")
  skip_if_not_installed("bslib")
  skip_if_not_installed("DT")
  skip_if_not_installed("bsicons")

  app <- .app_test_path()
  e <- new.env(parent = globalenv())
  suppressWarnings(sys.source(app, envir = e))
  pair <- as.data.frame(t(utils::combn(LETTERS[1:5], 2)),
                        stringsAsFactors = FALSE)
  names(pair) <- c("object_a", "object_b")
  beta <- setNames(seq(-1, 1, length.out = 5), LETTERS[1:5])
  set.seed(902)
  d <- do.call(rbind, lapply(sprintf("J%02d", 1:12), function(j) {
    z <- pair
    z$judge <- j
    p <- plogis(beta[z$object_a] - beta[z$object_b])
    z$winner <- ifelse(runif(nrow(z)) < p, z$object_a, z$object_b)
    z$group <- if (j <= "J06") "G1" else "G2"
    z
  }))
  d$group[d$judge == "J01" & seq_len(nrow(d)) == 2L] <- "changed"
  bt_current <- btl(d, "object_a", "object_b", winner = "winner",
                    judge = "judge")
  attr(bt_current, "rasch_app_source") <- list(
    data = d,
    settings = list(model_type = "btl", bt_a = "object_a",
                    bt_b = "object_b", bt_win = "winner",
                    bt_judge = "judge"),
    resources = list(), simulation = list())
  path <- tempfile(fileext = ".csv")
  stable_path <- tempfile(fileext = ".csv")
  on.exit(unlink(c(path, stable_path)), add = TRUE)
  write.csv(d, path, row.names = FALSE)
  stable <- d
  stable$group[stable$judge == "J01"] <- "G1"
  bt_stable <- btl(stable, "object_a", "object_b", winner = "winner",
                   judge = "judge")
  attr(bt_stable, "rasch_app_source") <- list(
    data = stable,
    settings = list(model_type = "btl", bt_a = "object_a",
                    bt_b = "object_b", bt_win = "winner",
                    bt_judge = "judge"),
    resources = list(), simulation = list())
  write.csv(stable, stable_path, row.names = FALSE)

  shiny::testServer(e$server, {
    session$setInputs(
      file = list(datapath = path, name = "comparisons.csv",
                  size = file.info(path)$size, type = "text/csv"),
      bt_judge = "judge", bdif_factors = "group")
    session$flushReact()
    btl_fit(bt_current)
    expect_error(bdif_factor_maps(), "varies within judge")
    expect_match(bdif_code_grp(), "judge_factor <- function", fixed = TRUE)
    expect_silent(parse(text = bdif_code_grp()))
    session$setInputs(file = list(
      datapath = stable_path, name = "stable.csv",
      size = file.info(stable_path)$size, type = "text/csv"))
    session$flushReact()
    btl_fit(bt_stable)
    maps <- bdif_factor_maps()
    expect_identical(unname(maps$group["J01"]), "G1")
  })
})

test_that("app captions escape item names and table decisions use strict cuts", {
  skip_on_cran()
  skip_if_not_installed("shiny")
  skip_if_not_installed("bslib")
  skip_if_not_installed("DT")
  skip_if_not_installed("bsicons")

  e <- new.env(parent = globalenv())
  suppressWarnings(sys.source(.app_test_path(), envir = e))
  set.seed(903)
  X <- matrix(rbinom(600, 1, .5), 120, 5,
              dimnames = list(NULL, c("<em>Injected</em>",
                                      paste0("I", 2:5))))
  f <- rasch(X)

  shiny::testServer(e$server, {
    fit_val(f)
    session$flushReact()
    html <- paste(as.character(output$chisq_caption), collapse = "")
    expect_false(grepl("<em>Injected</em>", html, fixed = TRUE))
    expect_match(html, "&lt;em&gt;Injected&lt;/em&gt;", fixed = TRUE)
    expect_match(as.character(colour_if("x<0.05")), "x<0.05", fixed = TRUE)
    expect_match(as.character(weight_if("x<0.05")), "x<0.05", fixed = TRUE)
    expect_match(as.character(colour_if("x<0.05")), "value===null",
                 fixed = TRUE)
    expect_match(as.character(colour_if("Math.abs(x)>=0.5")),
                 "Math.abs(x)>=0.5", fixed = TRUE)
  })
})

test_that("uploaded equating and Wright-map code reproduces input normalisation", {
  skip_on_cran()
  skip_if_not_installed("shiny")
  skip_if_not_installed("bslib")
  skip_if_not_installed("DT")
  skip_if_not_installed("bsicons")

  e <- new.env(parent = globalenv())
  suppressWarnings(sys.source(.app_test_path(), envir = e))
  bank_path <- tempfile(fileext = ".csv")
  map_path <- tempfile(fileext = ".csv")
  on.exit(unlink(c(bank_path, map_path)), add = TRUE)
  write.csv(data.frame(object = c(" O1 ", "O2"), location = c(-.2, .2),
                       se = c(.1, .1), m = c("2", "2")), bank_path,
            row.names = FALSE)
  write.csv(data.frame(item = c(" I1 ", "I2"), panel = c(" A ", "B")),
            map_path, row.names = FALSE)
  f <- rasch(matrix(rbinom(400, 1, .5), 100, 4,
                    dimnames = list(NULL, paste0("I", 1:4))))

  shiny::testServer(e$server, {
    fit_val(f)
    session$setInputs(
      bt_eq_file = list(datapath = bank_path, name = bank_path,
                        size = file.info(bank_path)$size, type = "text/csv"),
      wright_renderer = "wrightmap", wright_type = "thresholds",
      wright_person_panels = "", wright_item_panels = "uploaded",
      wright_item_map = list(datapath = map_path, name = map_path,
                             size = file.info(map_path)$size, type = "text/csv"),
      wright_person_style = "histogram", tg_bins = 35,
      tg_axis_mode = "standard")
    session$flushReact()

    parsed_bank <- bt_eq_bank()
    expect_identical(parsed_bank$object, c("O1", "O2"))
    expect_identical(attr(parsed_bank, "m"), 2L)
    code_env <- new.env(parent = globalenv())
    eval(parse(text = bt_eq_reference_code()), envir = code_env)
    expect_identical(code_env$bank$object, parsed_bank$object)
    expect_identical(attr(code_env$bank, "m"), 2L)
    expect_match(bt_eq_reference_code(), "[.]Machine[$]integer[.]max")

    parsed_map <- wright_item_arg()
    expect_identical(parsed_map, c(I1 = "A", I2 = "B"))
    code_env$fit <- f
    code_env$wright_map <- function(..., item_panels) item_panels
    generated_map <- eval(parse(text = wright_code()), envir = code_env)
    expect_identical(generated_map, parsed_map)
  })
})

test_that("changed dimensionality controls invalidate a restored t-test", {
  skip_on_cran()
  skip_if_not_installed("shiny")
  skip_if_not_installed("bslib")
  skip_if_not_installed("DT")
  skip_if_not_installed("bsicons")

  e <- new.env(parent = globalenv())
  suppressWarnings(sys.source(.app_test_path(), envir = e))
  set.seed(904)
  X <- matrix(rbinom(720, 1, .5), 120, 6,
              dimnames = list(NULL, paste0("I", 1:6)))
  f <- rasch(X)
  saved <- dimensionality_test(
    f, items_positive = colnames(X)[1:3],
    items_negative = colnames(X)[4:6], min_score_points = 2)

  shiny::testServer(e$server, {
    fit_val(f)
    session$flushReact()
    restored_subtest(saved)
    expect_identical(dim_res()$result_signature, saved$result_signature)
    session$setInputs(pca_component = 2)
    session$flushReact()
    expect_null(restored_subtest())
    expect_error(dim_res(), "press Run t-test")
  })
})

test_that("opening a project retains its dimensionality specification", {
  skip_on_cran()
  skip_if_not_installed("shiny")
  skip_if_not_installed("bslib")
  skip_if_not_installed("DT")
  skip_if_not_installed("bsicons")

  e <- new.env(parent = globalenv())
  suppressWarnings(sys.source(.app_test_path(), envir = e))
  set.seed(9041)
  X <- matrix(rbinom(840, 1, .5), 140, 6,
              dimnames = list(NULL, paste0("I", 1:6)))
  f <- rasch(X)
  saved <- dimensionality_test(f, component = 2, min_score_points = 2)
  project <- .seal_app_project(list(
    format = "rasch-shiny-project", schema = 2L,
    data = as.data.frame(X), model_type = "rasch", base_fit = f,
    rasch_steps = list(), btl_steps = list(), kept_fits = list(),
    kept_fit_code = list(), simulation = list(),
    settings = list(model_type = "rasch", pca_component = "2",
                    dim_boot_B = 0, dim_workers = "1", dim_boot_seed = 1),
    resources = list(), results = list(subtest = saved)))
  path <- tempfile(fileext = ".rasch")
  on.exit(unlink(path), add = TRUE)
  .save_app_project(project, path)

  shiny::testServer(e$server, {
    session$setInputs(project_file = list(
      datapath = path, name = "dimensionality.rasch",
      size = file.info(path)$size, type = "application/octet-stream"))
    session$flushReact()
    session$flushReact()
    expect_false(is.null(restored_subtest()))
    expect_identical(dim_res()$result_signature, saved$result_signature)
  })
})

test_that("saved analyses retain uploaded equating and item-panel resources", {
  set.seed(905)
  X <- matrix(rbinom(600, 1, .5), 120, 5,
              dimnames = list(NULL, paste0("I", 1:5)))
  f <- rasch(X)
  reference <- f$items[, c("item", "location", "se", "max")]
  panels <- data.frame(item = c("I1", "I2"), panel = c("A", "B"))
  project <- .seal_app_project(list(
    format = "rasch-shiny-project", schema = 2L,
    data = as.data.frame(X), model_type = "rasch", base_fit = f,
    rasch_steps = list(), btl_steps = list(), kept_fits = list(),
    kept_fit_code = list(), simulation = list(), results = list(),
    settings = list(eq_source = "csv", eq_shift = "mean",
                    eq_csv_independent = TRUE,
                    wright_item_panels = "uploaded"),
    resources = list(eq_reference = reference, wright_item_map = panels)))
  path <- tempfile(fileext = ".rasch")
  on.exit(unlink(path), add = TRUE)
  expect_no_error(.save_app_project(project, path))
  restored <- .read_app_project(path)
  expect_identical(restored$resources$eq_reference, reference)
  expect_identical(restored$resources$wright_item_map, panels)

  bad <- project
  bad$resources$wright_item_map$item[2] <- " I1 "
  bad <- .seal_app_project(bad)
  expect_error(.validate_app_project(bad), "invalid Wright-map item panel map")
})

test_that("saved CJ analyses retain a polytomous equating bank", {
  d <- simulate_btl(4, 8, 3, model = "polytomous", n_categories = 3,
                    seed = 906)
  fit <- btl(d, "object_a", "object_b", response = "response",
             judge = "judge")
  bank <- fit$objects[, c("object", "location", "se")]
  attr(bank, "m") <- fit$m
  project <- .seal_app_project(list(
    format = "rasch-shiny-project", schema = 2L,
    data = as.data.frame(d), model_type = "btl", base_fit = fit,
    rasch_steps = list(), btl_steps = list(), kept_fits = list(),
    kept_fit_code = list(), simulation = list(), results = list(),
    settings = list(bt_eq_independent = TRUE),
    resources = list(bt_eq_bank = bank)))
  path <- tempfile(fileext = ".rasch")
  on.exit(unlink(path), add = TRUE)

  expect_no_error(.save_app_project(project, path))
  restored <- .read_app_project(path)
  expect_identical(restored$resources$bt_eq_bank, bank)
  expect_identical(attr(restored$resources$bt_eq_bank, "m"), fit$m)

  wrong_scale <- project
  attr(wrong_scale$resources$bt_eq_bank, "m") <- fit$m + 1L
  wrong_scale <- .seal_app_project(wrong_scale)
  expect_error(.validate_app_project(wrong_scale),
               "response-scale metadata do not match")
})
