.selected_dif_app <- function() {
  for (package in c("shiny", "bslib", "DT", "bsicons"))
    skip_if_not_installed(package)
  env <- new.env(parent = globalenv())
  app_path <- testthat::test_path("..", "..", "inst", "shiny", "app.R")
  if (!file.exists(app_path))
    app_path <- system.file("shiny", "app.R", package = "rasch")
  suppressWarnings(sys.source(app_path, envir = env))
  env
}

test_that("app exports do not replace an unavailable factorial DIF analysis", {
  skip_on_cran()
  e <- .selected_dif_app()
  set.seed(9721)
  n <- 80L
  X <- matrix(rbinom(n * 6, 1, .5), n, 6,
               dimnames = list(NULL, paste0("I", 1:6)))
  factors <- data.frame(a = factor(rep(1:10, length.out = n)),
    b = factor(rep(1:10, each = 10, length.out = n)))
  f <- rasch(X, factors = factors)
  expect_error(dif_anova(f, effects = "factorial"),
               "no item yielded an estimable factorial ANOVA")
  main <- dif_anova(f, effects = "main")
  expect_true(any(is.finite(main$terms$p)))
  html <- tempfile(fileext = ".html")
  docx <- tempfile(fileext = ".docx")
  project <- tempfile(fileext = ".rasch")
  on.exit(unlink(c(html, docx, project)), add = TRUE)
  shiny::testServer(e$server, {
    sim_data(data.frame(X, factors))
    fit_val(f)
    session$setInputs(dif_effects = "factorial", dif_alpha = .01)
    session$flushReact()
    expect_error(dif_res(), "no item yielded an estimable factorial ANOVA")
    expect_null(app_dif_res())
    expect_error(app_dif_res(strict = TRUE), "selected DIF analysis is unavailable")
    expect_error(report_content(html, "html"), "cannot substitute a different analysis")
    expect_false(file.exists(html))
    expect_error(report_content(docx, "docx"), "cannot substitute a different analysis")
    expect_false(file.exists(docx))
    expect_error(output$dl_zip, "cannot substitute a different analysis")
    saved <- project_state()
    expect_identical(saved$settings$dif_effects, "factorial")
    expect_null(saved$results$dif)
    expect_no_error(.save_app_project(saved, project))
    session$setInputs(dif_effects = "main")
    session$flushReact()
    selected <- app_dif_res(strict = TRUE)
    expect_identical(selected$effects, "main")
    expect_identical(selected$alpha, .01)
    expect_identical(selected$summary, dif_res()$summary)
  })
})

test_that("strict app export retrieval permits fits without eligible DIF factors", {
  skip_on_cran()
  e <- .selected_dif_app()
  data <- simulate_efrm(n_per_group = 150, items_per_set = 5,
    n_sets = 1, n_groups = 2, seed = 9722)
  framed <- rasch_efrm(data, item_sets = attr(data, "truth")$item_sets,
    groups = "group", id = "id", boot_reps = 0)
  ordinary <- rasch(as.matrix(data[unlist(attr(data, "truth")$item_sets)]))
  expect_length(setdiff(names(framed$factors), framed$frame_group), 0L)
  shiny::testServer(e$server, {
    fit_val(ordinary)
    session$flushReact()
    expect_null(app_dif_res(strict = TRUE))
    fit_val(framed)
    session$flushReact()
    expect_null(app_dif_res(strict = TRUE))
  })
})

test_that("app exports retain the selected person-subset dimensionality analysis", {
  skip_on_cran()
  e <- .selected_dif_app()
  set.seed(9731)
  X <- matrix(rbinom(240 * 8, 1, .5), 240, 8,
               dimnames = list(NULL, paste0("I", 1:8)))
  f <- rasch_explanatory(X,
    predictors = data.frame(item = colnames(X), z = rep(0:3, each = 2)),
    formula = ~ z)
  html <- tempfile(fileext = ".html")
  docx <- tempfile(fileext = ".docx")
  on.exit(unlink(c(html, docx)), add = TRUE)
  shiny::testServer(e$server, {
    fit_val(f)
    session$flushReact()
    expect_null(app_subtest_res(strict = TRUE))
    session$setInputs(dim_pos = colnames(X)[1:4], dim_neg = colnames(X)[5:8],
                       dim_boot_B = 0, dim_boot_seed = 1, dim_workers = "1")
    session$setInputs(dim_apply = 1)
    chosen <- app_subtest_res(strict = TRUE)
    expect_identical(chosen$split, "manual")
    expect_identical(chosen$items_positive, colnames(X)[1:4])
    expect_identical(chosen$items_negative, colnames(X)[5:8])

    # Editing unapplied controls leaves the displayed run intact. Once Run
    # is pressed, rejecting that split must not retain the last result.
    session$setInputs(dim_neg = colnames(X)[4:8])
    expect_identical(app_subtest_res(strict = TRUE), chosen)
    session$setInputs(dim_apply = 2)
    expect_null(dim_computed())
    expect_null(dim_subsets())
    expect_error(app_subtest_res(strict = TRUE), "selected dimensionality t-test")

    # Bootstrap refits are deliberately unsupported for explanatory fits;
    # an export must not turn the requested manual/bootstrap analysis into
    # the default residual split without bootstrap calibration.
    session$setInputs(dim_neg = colnames(X)[5:8], dim_boot_B = 99)
    session$setInputs(dim_apply = 3)
    expect_null(dim_computed())
    expect_null(app_subtest_res())
    expect_error(report_content(html, "html"), "cannot substitute its default item split")
    expect_false(file.exists(html))
    expect_error(report_content(docx, "docx"), "cannot substitute its default item split")
    expect_false(file.exists(docx))
    expect_error(output$dl_zip, "cannot substitute its default item split")
  })
})

test_that("universal person-subset refusals do not block unrelated app exports", {
  skip_on_cran()
  e <- .selected_dif_app()
  set.seed(9732)
  X <- matrix(rbinom(80 * 6, 1, .5), 80, 6,
               dimnames = list(NULL, paste0("I", 1:6)))
  f <- rasch(X[rep(1:80, each = 2), ], id = rep(1:80, each = 2))
  shiny::testServer(e$server, {
    fit_val(f)
    session$setInputs(dim_pos = colnames(X)[1:3], dim_neg = colnames(X)[4:6],
                       dim_boot_B = 99)
    expect_null(app_subtest_res(strict = TRUE))
  })
})
