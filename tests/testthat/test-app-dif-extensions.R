.dif_ext_app <- function() {
  for (package in c("shiny", "bslib", "DT", "bsicons"))
    skip_if_not_installed(package)
  env <- new.env(parent = globalenv())
  app_path <- testthat::test_path("..", "..", "inst", "shiny", "app.R")
  if (!file.exists(app_path))
    app_path <- system.file("shiny", "app.R", package = "rasch")
  suppressWarnings(sys.source(app_path, envir = env))
  env
}

.dif_ext_data <- function() {
  d <- simulate_rasch(n_persons = 300, n_items = 8,
                      dif = list(items = c("I02", "I05"), uniform = 1.5),
                      n_groups = 2, seed = 12)
  X <- as.data.frame(d)
  X[, c(grep("^I", names(X), value = TRUE), "group")]
}

test_that("the bundle text parses to the list dif_anova() takes", {
  expect_null(.parse_dif_bundles(NULL))
  expect_null(.parse_dif_bundles(""))
  expect_null(.parse_dif_bundles("  \n\n"))
  expect_identical(
    .parse_dif_bundles("reading: I1, I2, I3\r\nmaths: I4,I5\n"),
    list(reading = c("I1", "I2", "I3"), maths = c("I4", "I5")))
  expect_identical(.parse_dif_bundles("a: I1, I1, I2"),
                   list(a = c("I1", "I2")))
  expect_error(.parse_dif_bundles("I1, I2"),
               "bundle line 1 needs the form name: item, item, ...")
  expect_error(.parse_dif_bundles("a: I1, I2\n: I3, I4"),
               "bundle line 2 needs the form")
  expect_error(.parse_dif_bundles("a: I1"),
               "bundle line 1 needs a name and at least two items")
  expect_true(.same_dif_bundles(NULL, list()))
  expect_true(.same_dif_bundles(list(a = c("I1", "I2")),
                                list(a = c("I1", "I2"))))
  expect_false(.same_dif_bundles(list(a = c("I1", "I2")),
                                 list(a = c("I2", "I1"))))
  expect_false(.same_dif_bundles(list(a = c("I1", "I2")), NULL))
  expect_false(.same_dif_bundles(list(a = c("I1", "I2")),
                                 list(b = c("I1", "I2"))))
})

test_that("a split item in a bundle stands for its group copies in the ANOVA", {
  dat <- .dif_ext_data()
  items <- grep("^I", names(dat), value = TRUE)
  f <- rasch(dat[, items], factors = data.frame(group = dat$group))
  b <- list(pair = c("I02", "I03"))
  expect_identical(.app_dif_bundles(b, f), b)
  s <- split_items(f, "I02", by = "group")
  expanded <- .app_dif_bundles(b, s)
  expect_identical(names(expanded), "pair")
  expect_setequal(expanded$pair, c("I02 (g1)", "I02 (g2)", "I03"))
  expect_no_error(dif_anova(s, bundles = expanded))
})

test_that("the DIF panel reads bundles, the Wald tests and test functioning", {
  skip_on_cran()
  e <- .dif_ext_app()
  dat <- .dif_ext_data()
  items <- grep("^I", names(dat), value = TRUE)
  f <- rasch(dat[, items], factors = data.frame(group = dat$group))
  shiny::testServer(e$server, {
    session$flushReact()
    sim_data(dat); fit_val(f); session$flushReact()
    session$setInputs(dif_effects = "main", dif_alpha = 0.05,
                      dif_bundles = "", dif_criterion = "anova",
                      dif_full = FALSE, dif_wald_full = FALSE)
    session$flushReact()
    expect_identical(output$dif_tbl_code,
                     'dif_anova(fit, p_adjust = "holm", alpha = 0.05)$summary')
    expect_null(dif_bundles())
    expect_null(dif_res()$bundles)

    # bundles reach the analysis, its footer and its frozen code
    session$setInputs(dif_bundles = "pair: I01, I03\nother: I04, I06, I07")
    session$flushReact()
    expect_match(output$dif_bundles_note$html,
                 "2 bundle(s): pair (2 items), other (3 items).", fixed = TRUE)
    expect_true(all(c("pair", "other") %in% dif_res()$summary$item))
    code <- output$dif_tbl_code
    expect_match(code, 'bundles = list(pair = c("I01", "I03"), other = c("I04", "I06", "I07"))',
                 fixed = TRUE)
    env <- new.env(parent = globalenv()); env$fit <- f
    expect_equal(eval(parse(text = code), envir = env), dif_res()$summary)
    expect_match(output$dif_full_tbl_code, "bundles = list(", fixed = TRUE)

    # a bundle row has no curve, no post-hoc comparison and no split
    session$setInputs(dif_tbl_rows_selected = which(dif_tbl_items() == "pair"))
    session$flushReact()
    expect_identical(dif_sel_item(), "pair")
    expect_error(output$dif_icc, "no single characteristic curve")
    expect_s3_class(dif_posthoc_res(), "error")
    expect_match(conditionMessage(dif_posthoc_res()), "no single resolved location")
    session$setInputs(make_split = 1L); session$flushReact()
    expect_null(active_step())

    # a problem in the text withholds the analysis rather than dropping it
    session$setInputs(dif_bundles = "pair: I01, ZZZ"); session$flushReact()
    expect_match(output$dif_bundles_note$html, "text-danger", fixed = TRUE)
    expect_match(output$dif_bundles_note$html, "ZZZ", fixed = TRUE)
    expect_error(output$dif_tbl, "Item bundles: ")
    expect_null(app_dif_res())
    session$setInputs(dif_bundles = "pair: I01, I03"); session$flushReact()

    # conditional Wald tests with the selected row's level locations
    w <- dif_wald_res()
    expect_s3_class(w, "rasch_dif_wald")
    expect_identical(output$dif_wald_tbl_code,
                     'dif_wald(fit, p_adjust = "holm", alpha = 0.05)$summary')
    expect_match(output$dif_wald_note$html, "item-by-factor tests significant")
    expect_match(output$dif_wald_note$html, "second level minus first")
    session$setInputs(dif_wald_tbl_rows_selected = 2L); session$flushReact()
    expect_identical(dif_wald_sel()$item, w$summary$item[2L])
    lv <- output$dif_wald_levels_tbl_csv
    expect_true(!is.null(output$dif_wald_levels_tbl))

    # no test functioning before a split
    expect_null(fit()$split_map)
    expect_error(dtf_args(), class = "shiny.silent.error")

    # the automatic run under the Wald criterion, recorded and reproducible
    session$setInputs(dif_criterion = "wald"); session$flushReact()
    session$setInputs(resolve_all = 1L); session$flushReact()
    rr <- resolve_res()
    expect_identical(rr$criterion, "wald")
    expect_identical(rr$run_criterion, "wald")
    expect_identical(rr$run_effects, "main")
    expect_equal(rr$n_splits, 2)
    expect_match(output$resolve_summary$html, "by the conditional Wald criterion")
    code <- active_step()$code
    expect_match(code, 'criterion = "wald"', fixed = TRUE)
    env <- new.env(parent = globalenv()); env$fit <- f
    eval(parse(text = code), envir = env)
    expect_equal(env$dif_resolution$splits, rr$splits)

    # test functioning defaults to the split term and the first level
    expect_identical(dtf_split_term(), "group")
    expect_identical(dtf_by_choices(), "group")
    session$setInputs(dtf_by = "group"); session$flushReact()
    expect_identical(dtf_levels(), c("g1", "g2"))
    session$setInputs(dtf_reference = "g1", dtf_group = "g2",
                      dtf_test_full = FALSE, dtf_items_full = FALSE,
                      dtf_bundles_full = FALSE)
    session$flushReact()
    r <- dtf_res()
    expect_s3_class(r, "rasch_dtf")
    expect_match(output$dtf_note$html,
                 "8 item(s) compared by group against the reference level g1: 2 split, 6 anchoring",
                 fixed = TRUE)
    expect_identical(
      output$dtf_test_tbl_code,
      'dtf(fit, by = "group", reference = "g1", bundles = list(pair = c("I01", "I03")), p_adjust = "holm", alpha = 0.05)$test')
    expect_match(output$dtf_plot_code, 'plot_dtf(dtf_result, group = "g2")',
                 fixed = TRUE)
    env <- new.env(parent = globalenv()); env$fit <- fit()
    eval(parse(text = output$dtf_plot_code), envir = env)
    expect_equal(env$dtf_result$test, r$test)
    for (id in c("dtf_test_tbl", "dtf_items_tbl", "dtf_scores_tbl",
                 "dtf_bundles_tbl"))
      expect_no_error(output[[id]])
    expect_no_error(output$dtf_plot)

    # a bundle over split items is typed over the sources: dtf() takes it
    # as typed and dif_anova() takes each split item as its copies
    session$setInputs(dif_bundles = "pair: I02, I05"); session$flushReact()
    expect_identical(dif_bundles(), list(pair = c("I02", "I05")))
    expect_setequal(dif_anova_bundles()$pair,
                    c("I02 (g1)", "I02 (g2)", "I05 (g1)", "I05 (g2)"))
    expect_true("pair" %in% dif_res()$summary$item)
    expect_match(output$dif_tbl_code, '"I02 (g1)"', fixed = TRUE)
    expect_match(output$dtf_test_tbl_code,
                 'bundles = list(pair = c("I02", "I05"))', fixed = TRUE)
    expect_identical(dtf_res()$bundles$bundle, "pair")
  })
})

test_that("bundles and the resolution criterion survive a project round trip", {
  skip_on_cran()
  e <- .dif_ext_app()
  dat <- .dif_ext_data()
  items <- grep("^I", names(dat), value = TRUE)
  f <- rasch(dat[, items], factors = data.frame(group = dat$group))
  path_a <- tempfile(fileext = ".rasch")
  path_b <- tempfile(fileext = ".rasch")
  path_c <- tempfile(fileext = ".rasch")
  on.exit(unlink(c(path_a, path_b, path_c)), add = TRUE)
  shiny::testServer(e$server, {
    session$flushReact()
    sim_data(dat); fit_val(f); session$flushReact()
    session$setInputs(dif_effects = "main", dif_alpha = 0.05,
                      dif_bundles = "pair: I01, I03", dif_criterion = "wald")
    session$flushReact()
    bs <- suppressWarnings(dif_bootstrap(fit(), dif_res(), B = 19,
                                         workers = 1, seed = 3))
    dif_boot_val(list(db = bs, B = 19L, seed = 3L))
    session$flushReact()
    saved <- project_state()
    expect_identical(saved$settings$dif_bundles, "pair: I01, I03")
    expect_identical(saved$settings$dif_criterion, "wald")
    expect_identical(saved$results$dif$bundles, list(pair = c("I01", "I03")))
    expect_no_error(.save_app_project(saved, path_a))
    session$setInputs(resolve_all = 1L); session$flushReact()
    expect_null(dif_boot_val())
    saved <- project_state()
    expect_identical(saved$results$resolve$criterion, "wald")
    expect_no_error(.save_app_project(saved, path_b))
  })

  # the bootstrap survives the restore: the text is recorded as restored,
  # so the change observer does not read it as new input
  shiny::testServer(e$server, {
    session$flushReact()
    session$setInputs(project_file = list(
      datapath = path_a, name = "a.rasch", size = file.info(path_a)$size,
      type = "application/octet-stream"))
    session$flushReact()
    expect_false(is.null(dif_boot_val()))
    session$setInputs(dif_bundles = "pair: I01, I03", dif_criterion = "wald",
                      dif_effects = "main", dif_alpha = 0.05)
    session$flushReact()
    expect_false(is.null(dif_boot_val()))
    expect_identical(dif_res()$bundles, list(pair = c("I01", "I03")))
    session$setInputs(dif_bundles = "pair: I01, I04"); session$flushReact()
    expect_null(dif_boot_val())
  })
  shiny::testServer(e$server, {
    session$flushReact()
    session$setInputs(project_file = list(
      datapath = path_b, name = "b.rasch", size = file.info(path_b)$size,
      type = "application/octet-stream"))
    session$flushReact()
    expect_identical(resolve_res()$criterion, "wald")
    session$setInputs(dif_bundles = "pair: I01, I03", dif_criterion = "wald",
                      dif_effects = "main", dif_alpha = 0.05)
    session$flushReact()
    expect_match(output$resolve_summary$html, "conditional Wald criterion")
    expect_no_error(.validate_app_project(project_state()))
  })

  # a project saved before the control existed restores blank text, so
  # text left in the session never reaches the older analysis
  p <- readRDS(path_a)
  p$settings$dif_bundles <- NULL
  p$results$dif <- NULL
  p$results$dif_bootstrap <- NULL
  saveRDS(.seal_app_project(unclass(p)), path_c, version = 3)
  shiny::testServer(e$server, {
    session$flushReact()
    session$setInputs(dif_bundles = "stale: I01, I02")
    session$setInputs(project_file = list(
      datapath = path_c, name = "c.rasch", size = file.info(path_c)$size,
      type = "application/octet-stream"))
    session$flushReact()
    restored <- e$.restored_input_values(restored_project_settings())
    expect_identical(restored$dif_bundles, "")
  })

  # the validator pins the saved analysis to the restored text and the
  # resolution to a known criterion
  p <- readRDS(path_a)
  p$settings$dif_bundles <- "pair: I01, I04"
  expect_error(.validate_app_project(.seal_app_project(unclass(p))),
               "does not match the restored item bundles")
  p <- readRDS(path_a)
  p$settings$dif_bundles <- "not a bundle"
  expect_error(.validate_app_project(.seal_app_project(unclass(p))),
               "does not match the restored item bundles")
  p <- readRDS(path_b)
  p$results$resolve$criterion <- "bogus"
  expect_error(.validate_app_project(.seal_app_project(unclass(p))),
               "superseded")
  p <- readRDS(path_b)
  p$results$resolve$criterion <- NULL
  expect_no_error(.validate_app_project(.seal_app_project(unclass(p))))
})
