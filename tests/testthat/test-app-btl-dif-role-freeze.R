.app_btl_role_test_path <- function() {
  p <- testthat::test_path("..", "..", "inst", "shiny", "app.R")
  if (!file.exists(p)) p <- system.file("shiny", "app.R", package = "rasch")
  p
}

test_that("CJ DIF refuses a permuted judge role and emits fitted-role code", {
  skip_on_cran()
  for (pkg in c("shiny", "bslib", "DT", "bsicons"))
    skip_if_not_installed(pkg)

  e <- new.env(parent = globalenv())
  suppressWarnings(sys.source(.app_btl_role_test_path(), envir = e))
  d <- as.data.frame(simulate_btl(8, 12, 25, seed = 912))
  j <- as.integer(sub("^J", "", d$judge))
  d$judge_group <- ifelse(j <= 6, "A", "B")
  d$judge_perm <- paste0("J", ifelse(j %% 12 == 0, 1, j %% 12 + 1))
  bt <- btl(d, "object_a", "object_b", winner = "winner", judge = "judge")
  attr(bt, "rasch_app_source") <- list(
    data = d,
    settings = list(model_type = "btl", bt_a = "object_a",
                    bt_b = "object_b", bt_win = "winner",
                    bt_judge = "judge"),
    resources = list(), simulation = list())
  shiny::testServer(e$server, {
    btl_fit(bt); sim_data(d)
    session$setInputs(bt_judge = "judge_perm",
                      bt_jfactors = "judge_group",
                      bdif_factors = "judge_group")
    session$flushReact()
    expect_error(bdif_factor_maps(), "does not reproduce the fitted judge IDs")
    expect_null(bdif_res())
    # Before a run, reproducible code names the fitted role rather than the
    # edited sidebar column; this prevents a misleading runnable snippet.
    expect_match(bdif_code_grp(), "dat\\$judge[^_a-zA-Z]", perl = TRUE)
    expect_false(grepl("dat\\$judge_perm", bdif_code_grp(), fixed = FALSE))
    session$setInputs(bdif_run = 1)
    session$flushReact()
    expect_null(bdif_res())
  })
})

test_that("CJ DIF accepts equivalent judge aliases and records provenance", {
  skip_on_cran()
  for (pkg in c("shiny", "bslib", "DT", "bsicons"))
    skip_if_not_installed(pkg)

  e <- new.env(parent = globalenv())
  suppressWarnings(sys.source(.app_btl_role_test_path(), envir = e))
  d <- as.data.frame(simulate_btl(5, 40, 25, seed = 913))
  d$judge_alias <- d$judge
  j <- as.integer(sub("^J", "", d$judge))
  d$judge_group <- ifelse(j <= 20, "A", "B")
  d$judge_perm <- paste0("J", ifelse(j %% 40 == 0, 1, j %% 40 + 1))
  bt <- btl(d, "object_a", "object_b", winner = "winner", judge = "judge")
  attr(bt, "rasch_app_source") <- list(
    data = d,
    settings = list(model_type = "btl", bt_a = "object_a",
                    bt_b = "object_b", bt_win = "winner",
    bt_judge = "judge"),
    resources = list(), simulation = list())
  project_path <- tempfile(fileext = ".rasch")
  on.exit(unlink(project_path), add = TRUE)

  shiny::testServer(e$server, {
    btl_fit(bt); sim_data(d)
    session$setInputs(bt_judge = "judge_alias",
                      bt_jfactors = "judge_group",
                      bdif_factors = "judge_group",
                      bdif_effects = "main", bdif_alpha = .05)
    session$flushReact()
    session$setInputs(bdif_run = 1)
    session$flushReact()
    expect_s3_class(bdif_res(), "rasch_btl_dif")
    expect_identical(bdif_meta()$judge_col, "judge_alias")
    expect_identical(bdif_meta()$fitted_judge_col, "judge")
    env <- new.env(parent = globalenv())
    env$dat <- d; env$bt <- bt
    eval(parse(text = bdif_code_grp()), envir = env)
    expect_identical(env$factors, bdif_res()$bootstrap_design$factors)

    saved <- project_state()
    .save_app_project(saved, project_path)
    session$setInputs(project_file = list(
      datapath = project_path, name = "cj.rasch",
      size = file.info(project_path)$size,
      type = "application/octet-stream"))
    session$flushReact()
    expect_s3_class(bdif_res(), "rasch_btl_dif")
    session$setInputs(bt_judge = "judge_perm", bdif_run = 2)
    session$flushReact()
    expect_null(bdif_res())
  })
})

test_that("CJ DIF code refuses a legacy fit without source provenance", {
  skip_on_cran()
  for (pkg in c("shiny", "bslib", "DT", "bsicons"))
    skip_if_not_installed(pkg)

  e <- new.env(parent = globalenv())
  suppressWarnings(sys.source(.app_btl_role_test_path(), envir = e))
  d <- as.data.frame(simulate_btl(4, 12, 25, seed = 914))
  bt <- btl(d, "object_a", "object_b", winner = "winner", judge = "judge")
  shiny::testServer(e$server, {
    btl_fit(bt); sim_data(d)
    session$setInputs(bt_judge = "judge", bdif_factors = "judge_group")
    session$flushReact()
    expect_error(bdif_factor_maps(), "predates saved run metadata")
    expect_error(bdif_code_grp(), "has no saved source metadata")
  })
})
