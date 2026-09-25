.rank_app <- function() {
  for (package in c("shiny", "bslib", "DT", "bsicons"))
    skip_if_not_installed(package)
  env <- new.env(parent = globalenv())
  app_path <- testthat::test_path("..", "..", "inst", "shiny", "app.R")
  if (!file.exists(app_path))
    app_path <- system.file("shiny", "app.R", package = "rasch")
  suppressWarnings(sys.source(app_path, envir = env))
  env
}

.rank_demo <- function() {
  path <- testthat::test_path("..", "..", "inst", "shiny", "examples.R")
  if (!file.exists(path))
    path <- system.file("shiny", "examples.R", package = "rasch")
  e <- new.env(parent = asNamespace("rasch"))
  sys.source(path, envir = e)
  e$.demo_pl()
}

test_that("the bundled rankings example is one ranked object per row", {
  d <- .rank_demo()
  expect_identical(names(d), c("judge", "ranking", "object", "rank"))
  expect_equal(nrow(d), 24 * 3 * 5)
  expect_true(all(table(d$ranking) == 5))
  expect_identical(sort(unique(d$object)), sprintf("O%d", 1:8))
  k <- pl(d, judge = "judge")
  expect_s3_class(k, "rasch_pl")
  expect_true(k$converged)
  expect_equal(k$n_rankings, 72L)
  expect_equal(nrow(k$judges), 24L)
  # the erratic judge stands out as surprising
  j7 <- k$judges[k$judges$judge == "J07", ]
  expect_gt(j7$mean_surprise_z, 1)
  expect_gt(j7$infit_ms, 1.3)
  # the withheld invariance probability prints as withheld, not blank
  expect_output(print(k), "p = withheld; objects moving")
})

test_that("a rank analysis runs from the sidebar and is reproduced by its code", {
  skip_on_cran()
  e <- .rank_app()
  d <- .rank_demo()
  shiny::testServer(e$server, {
    session$flushReact()
    session$setInputs(demo_choice = "pl", model_type = "btl",
                      bt_layout = "ranks", pl_ranking = "ranking",
                      pl_object = "object", pl_rank = "rank",
                      pl_judge = "judge", pl_se = "sandwich",
                      pl_split = "first", pl_ties = "drop",
                      maxit = 100, tol = 1e-8, ng_auto = TRUE,
                      pl_objects_full = FALSE)
    session$flushReact()
    expect_identical(raw_data(), d)
    expect_null(pl_anchors_in())
    session$setInputs(run = 1L); session$flushReact()
    k <- pl_fit()
    expect_s3_class(k, "rasch_pl")
    expect_null(fit_or_null())
    expect_null(btl_fit())

    # the disclosed code reproduces the fitted object
    code <- current_rcode()
    expect_match(code, 'dat <- rasch:::.app_example_data("pl")', fixed = TRUE)
    expect_match(code, 'rk <- pl(dat,', fixed = TRUE)
    expect_match(code, 'judge = "judge"', fixed = TRUE)
    expect_false(grepl("anchors", code, fixed = TRUE))
    env <- new.env(parent = globalenv())
    eval(parse(text = code), envir = env)
    expect_equal(env$rk$objects, k$objects)
    expect_equal(env$rk$loglik, k$loglik)

    # every card renders and names the object its code reads
    for (id in c("pl_boxes", "pl_fitsum_tbl", "pl_reversal_tbl", "pl_map",
                 "pl_objects_tbl", "pl_judges_tbl", "pl_rankings_tbl",
                 "pl_invariance_tbl"))
      expect_no_error(output[[id]])
    expect_identical(output$pl_objects_tbl_code, "rk$objects")
    expect_identical(output$pl_map_code, "plot_pl(rk)")
    expect_identical(output$pl_invariance_tbl_code, "rk$invariance$objects")
    expect_match(output$pl_fitsum_notes$html, "invariance likelihood-ratio p value is withheld")
    expect_match(output$pl_invariance_note$html, "first choice", fixed = TRUE)
    expect_match(output$pl_invariance_note$html, "none", fixed = TRUE)
    csv <- output$pl_fitsum_tbl_csv
    expect_true(!is.null(csv))
    expect_match(output$nav_status$html, "Rankings", fixed = TRUE)
    expect_match(output$nav_status$html, "72 rankings", fixed = TRUE)

    # a rank analysis is not stored, reported or archived
    expect_error(project_state(), "no analysis to save")
    expect_error(report_content(), "no fit to report")

    # a missing role stops the run and leaves the analysis in place
    session$setInputs(pl_rank = e$NONE); session$flushReact()
    session$setInputs(run = 2L); session$flushReact()
    expect_identical(pl_fit(), k)

    # a paired-comparison run replaces the rank analysis
    session$setInputs(demo_choice = "btl", bt_layout = "pairs",
                      bt_a = "object_a", bt_b = "object_b", bt_win = "winner",
                      bt_judge = "judge", bt_count = e$NONE, bt_ties = "drop",
                      bt_response = "", bt_margin = "", bt_order = "",
                      bt_position = FALSE, rasch_calibration = "free")
    session$flushReact()
    session$setInputs(run = 3L); session$flushReact()
    expect_null(pl_fit())
    expect_s3_class(btl_fit(), "rasch_btl")
  })
})

test_that("the anchor CSV for a rank analysis reaches pl() and its code", {
  skip_on_cran()
  e <- .rank_app()
  d <- .rank_demo()
  anchors <- tempfile("anchors", fileext = ".csv")
  on.exit(unlink(anchors), add = TRUE)
  write.csv(data.frame(object = c("O1", "O8"), location = c(-1.5, 1.5)),
            anchors, row.names = FALSE)
  shiny::testServer(e$server, {
    session$flushReact()
    session$setInputs(demo_choice = "pl", model_type = "btl",
                      bt_layout = "ranks", pl_ranking = "ranking",
                      pl_object = "object", pl_rank = "rank",
                      pl_judge = e$NONE, pl_se = "model",
                      pl_split = "half", pl_ties = "drop",
                      maxit = 100, tol = 1e-8, ng_auto = TRUE,
                      pl_anchor_file = list(name = "anchors.csv",
                                            datapath = anchors))
    session$flushReact()
    expect_identical(pl_anchors_in(), c(O1 = -1.5, O8 = 1.5))
    session$setInputs(run = 1L); session$flushReact()
    k <- pl_fit()
    expect_s3_class(k, "rasch_pl")
    expect_identical(k$anchors, c(O1 = -1.5, O8 = 1.5))
    expect_identical(k$se_type, "model")
    expect_null(k$judges)
    code <- current_rcode()
    expect_match(code, 'rk_anchors <- read.csv("anchors.csv"', fixed = TRUE)
    expect_match(code, "anchors = rk_anchor_values", fixed = TRUE)
    expect_match(code, 'se = "model"', fixed = TRUE)
    expect_match(code, 'split = "half"', fixed = TRUE)
    expect_false(grepl("judge =", code, fixed = TRUE))
    env <- new.env(parent = globalenv())
    code <- sub('read.csv("anchors.csv"', paste0('read.csv("', anchors, '"'),
                code, fixed = TRUE)
    eval(parse(text = code), envir = env)
    expect_equal(env$rk$objects$location, k$objects$location)
    expect_match(output$pl_fitsum_tbl$html, "Anchored objects", fixed = TRUE)
  })
})

test_that("a rank anchor CSV is checked before it is used", {
  skip_on_cran()
  e <- .rank_app()
  bad <- tempfile("anchors", fileext = ".csv")
  on.exit(unlink(bad), add = TRUE)
  shiny::testServer(e$server, {
    session$flushReact()
    write.csv(data.frame(name = "O1", location = 1), bad, row.names = FALSE)
    session$setInputs(pl_anchor_file = list(name = "a.csv", datapath = bad))
    expect_error(pl_anchors_in(), "needs columns object, location")
    write.csv(data.frame(object = c("O1", ""), location = c(1, 2)), bad,
              row.names = FALSE)
    session$setInputs(pl_anchor_file = list(name = "b.csv", datapath = bad))
    expect_error(pl_anchors_in(), "non-blank object name")
    write.csv(data.frame(object = "O1", location = "x"), bad,
              row.names = FALSE)
    session$setInputs(pl_anchor_file = list(name = "c.csv", datapath = bad))
    expect_error(pl_anchors_in(), "finite location")
  })
})

test_that("the comparison layout defaults to pairs for saved analyses", {
  e <- .rank_app()
  vals <- e$.restored_input_values(list(model_type = "btl"))
  expect_identical(vals$bt_layout, "pairs")
  vals <- e$.restored_input_values(list(model_type = "btl",
                                        bt_layout = "ranks"))
  expect_identical(vals$bt_layout, "ranks")
})
