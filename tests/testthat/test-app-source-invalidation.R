.app_source_invalidation_path <- function() {
  p <- testthat::test_path("..", "..", "inst", "shiny", "app.R")
  if (!file.exists(p)) p <- system.file("shiny", "app.R", package = "rasch")
  p
}

test_that("changing upload, example, or simulation clears the active fit", {
  skip_if_not_installed("shiny")
  skip_if_not_installed("bslib")
  skip_if_not_installed("DT")
  skip_if_not_installed("bsicons")

  e <- new.env(parent = globalenv())
  suppressWarnings(sys.source(.app_source_invalidation_path(), envir = e))
  set.seed(91)
  X <- matrix(rbinom(100 * 5, 1, .5), 100, 5,
              dimnames = list(NULL, paste0("I", 1:5)))
  current <- rasch(X)
  path <- tempfile(fileext = ".csv")
  on.exit(unlink(path), add = TRUE)
  write.csv(data.frame(I1 = 0:2, I2 = 1:3, I3 = 0:2,
                       I4 = 1:3, I5 = 0:2), path, row.names = FALSE)

  shiny::testServer(e$server, {
    kept_fits(list(previous = current))
    kept_fit_code(list(previous = "fit <- rasch(old_data)"))
    fit_val(current)
    rcode_str("library(rasch)\ndat <- read.csv(\"old.csv\")\nfit <- rasch(dat)")
    session$setInputs(file = list(datapath = path, name = "new.csv",
                                  size = file.info(path)$size,
                                  type = "text/csv"))
    session$flushReact()
    expect_null(fit_val())
    expect_null(btl_fit())
    expect_length(analysis_steps(), 0L)
    expect_length(btl_analysis_steps(), 0L)
    expect_null(current_rcode())
    expect_equal(nrow(raw_data()), 3L)
    expect_identical(kept_fits(), list(previous = current))
    expect_identical(kept_fit_code(), list(previous = "fit <- rasch(old_data)"))

    fit_val(current)
    session$setInputs(demo_choice = "dich")
    session$flushReact()
    expect_null(fit_val())

    fit_val(current)
    session$setInputs(demo_choice = "none")
    session$flushReact()
    expect_null(fit_val())

    fit_val(current)
    session$setInputs(
      sim_layout = "rasch", sim_seed = 1, sr_persons = 30, sr_items = 5,
      sr_model = "dichotomous", sr_over = 0, sr_under = 0,
      sr_guess = FALSE, sr_2d = FALSE, sr_dep = FALSE, sr_dif = FALSE,
      sr_style = FALSE, sr_speeded = 0, sr_careless = 0, sr_missing = 0,
      sr_mean = 0, sr_sd = 1, sr_dist = "normal",
      sr_diff = c(-2.5, 2.5), sim_go = 1)
    session$flushReact()
    expect_null(fit_val())
    expect_equal(nrow(raw_data()), 30L)
  })
})

test_that("HTML reports label withheld total-fit probabilities", {
  set.seed(212)
  d <- simulate_rasch(120, 8, seed = 212)
  repeated <- d[rep(seq_len(nrow(d)), each = 2L), ]
  fit <- rasch(repeated, id = "id", items = sprintf("I%02d", 1:8))
  expect_true(isTRUE(.has_repeated_residual_units(fit)))
  expect_true(is.na(fit$total_chisq_p))

  out <- tempfile(fileext = ".html")
  on.exit(unlink(out), add = TRUE)
  suppressWarnings(report_html(fit, out, dpi = 40))
  html <- paste(readLines(out, warn = FALSE), collapse = "\n")
  expect_match(html, "p = unavailable", fixed = TRUE)
  expect_false(grepl("p = ).", html, fixed = TRUE))
})
