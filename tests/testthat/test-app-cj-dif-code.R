test_that("CJ DIF code reproduces the displayed judge maps after metadata edits", {
  skip_on_cran()
  for (pkg in c("shiny", "bslib", "DT", "bsicons"))
    skip_if_not_installed(pkg)

  e <- new.env(parent = globalenv())
  app_path <- test_path("..", "..", "inst", "shiny", "app.R")
  if (!file.exists(app_path))
    app_path <- system.file("shiny", "app.R", package = "rasch")
  suppressWarnings(sys.source(
    app_path, envir = e))
  d <- as.data.frame(simulate_btl(4, 40, 160, seed = 782))
  judge_number <- as.integer(sub("J", "", d$judge))
  d$judge_alias <- d$judge
  d[["judge group"]] <- ifelse(judge_number <= 20, "A", "B")
  bt <- btl(d, "object_a", "object_b", winner = "winner", judge = "judge")
  source <- list(data = d, settings = list(
    model_type = "btl", bt_a = "object_a", bt_b = "object_b",
    bt_win = "winner", bt_judge = "judge"),
    resources = list(), simulation = list())
  attr(bt, "rasch_app_source") <- source
  revised <- d
  revised[["judge group"]] <- ifelse(judge_number %% 2, "A", "B")
  revised$training <- ifelse(judge_number <= 20, "trained \"yes\"", "untrained")
  revised$training[judge_number == 40] <- NA_character_
  project_path <- tempfile(fileext = ".rasch")
  data_path <- tempfile(fileext = ".csv")
  on.exit(unlink(c(project_path, data_path)), add = TRUE)

  shiny::testServer(e$server, {
    btl_fit(bt)
    sim_data(revised)
    session$setInputs(bt_judge = "judge_alias",
                      bt_jfactors = c("judge group", "training"),
                      bdif_factors = c("judge group", "training"),
                      bdif_effects = "main", bdif_alpha = .05)
    session$flushReact()
    session$setInputs(bdif_run = 1)
    session$flushReact()
    result <- bdif_res()
    expect_s3_class(result, "rasch_btl_dif")
    maps <- result$bootstrap_design$factors
    expect_true(is.na(maps$training["J40"]))
    old_group <- setNames(d[["judge group"]], d$judge)[names(maps[[1]])]
    expect_equal(sum(maps[[1]] != old_group), 20)

    # The comparison data in a saved project belong to its calibration.
    # Reconstructing the DIF maps from those data used the old groups.
    saved <- project_state()
    expect_identical(saved$data, d)
    expect_identical(saved$settings$bt_judge, "judge")
    expect_identical(saved$results$btl_dif_meta$judge_col, "judge_alias")
    env <- new.env(parent = globalenv())
    env$dat <- saved$data
    env$bt <- bt
    eval(parse(text = bdif_code_grp()), envir = env)
    expect_identical(env$factors, maps)
    expect_equal(eval(parse(text = output$bdif_anova_tbl_code), envir = env),
                 result$summary)

    # Neither changing the controls nor changing the source metadata after
    # Run is allowed to rewrite code for the result still on screen.
    code <- bdif_code_grp()
    session$setInputs(bt_judge = "missing column", bdif_factors = "missing factor")
    sim_data(d)
    session$flushReact()
    expect_identical(bdif_code_grp(), code)
    expect_equal(eval(parse(text = output$bdif_anova_tbl_code), envir = env),
                 result$summary)

    # Reopening retains the displayed design even though the stored base
    # data predate the metadata correction.
    .save_app_project(saved, project_path)
    session$setInputs(project_file = list(
      datapath = project_path, name = "cj.rasch",
      size = file.info(project_path)$size, type = "application/octet-stream"))
    session$flushReact()
    expect_identical(bdif_res()$bootstrap_design$factors, maps)
    expect_identical(bdif_code_grp(), code)
    expect_equal(eval(parse(text = output$bdif_anova_tbl_code), envir = env),
                 result$summary)

    # Run uses the restored assignments, including missing judge metadata,
    # without changing any of the calibration's raw source columns.
    session$setInputs(bt_judge = "judge", bdif_factors = c("judge group", "training"),
                      bdif_run = 2)
    session$flushReact()
    expect_identical(bdif_res()$bootstrap_design$factors, maps)
    expect_equal(bdif_res()$summary, result$summary)
    expect_identical(raw_data(), d)

    # A genuine new upload is an explicit metadata replacement, even if its
    # bytes happen to match the calibration data stored in the project.
    write.csv(d, data_path, row.names = FALSE)
    session$setInputs(file = list(
      datapath = data_path, name = "new-judge-metadata.csv",
      size = file.info(data_path)$size, type = "text/csv"))
    session$flushReact()
    # A genuine upload replaces the source; the previous calibration and its
    # downstream DIF result are invalidated before the new columns are used.
    expect_null(btl_fit())
    expect_null(bdif_res())
    expect_error(bdif_factor_maps(), "run a Comparative Judgement analysis first")
    session$setInputs(bdif_factors = "judge group")
    session$flushReact()
    expect_error(bdif_factor_maps(), "run a Comparative Judgement analysis first")
  })
})
