.report_app <- function() {
  for (package in c("shiny", "bslib", "DT", "bsicons", "rmarkdown"))
    skip_if_not_installed(package)
  skip_if_not(rmarkdown::pandoc_available("2.8"))
  env <- new.env(parent = globalenv())
  app_path <- testthat::test_path("..", "..", "inst", "shiny", "app.R")
  if (!file.exists(app_path))
    app_path <- system.file("shiny", "app.R", package = "rasch")
  suppressWarnings(sys.source(app_path, envir = env))
  env
}

test_that("the app's report export offers Markdown in place of HTML", {
  skip_on_cran()
  e <- .report_app()
  ui <- as.character(e$ui)
  expect_match(ui, 'value="md"', fixed = TRUE)
  expect_false(grepl('value="html"', ui, fixed = TRUE))
  expect_match(ui, 'title="Analysis report (Markdown)"', fixed = TRUE)
})

test_that("the app writes its Markdown report from the Export page and the navbar", {
  skip_on_cran()
  e <- .report_app()
  set.seed(9731)
  X <- matrix(rbinom(90 * 6, 1, .5), 90, 6,
              dimnames = list(NULL, paste0("I", 1:6)))
  f <- rasch(X)
  shiny::testServer(e$server, {
    sim_data(data.frame(X))
    fit_val(f)
    session$flushReact()
    md <- output$dl_report
    expect_true(file.exists(md))
    lines <- readLines(md, warn = FALSE)
    expect_identical(lines[1L], "# Rasch measurement analysis")
    expect_true("## Analysis summary" %in% lines)
    expect_true("## Diagnostic figures" %in% lines)
    expect_false(any(grepl("<style", lines, fixed = TRUE)))
    nav <- output$dl_report_nav
    expect_true(file.exists(nav))
    expect_identical(readLines(nav, warn = FALSE)[-2L], lines[-2L])
  })
})
