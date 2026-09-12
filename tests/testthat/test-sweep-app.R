.sweep_app_path <- function() {
  app <- testthat::test_path("..", "..", "inst", "shiny", "app.R")
  if (!file.exists(app)) app <- system.file("shiny", "app.R", package = "rasch")
  app
}

test_that("DT cell flags decide from the cell value, not the row index", {
  skip_on_cran()
  skip_if_not_installed("shiny")
  skip_if_not_installed("bslib")
  skip_if_not_installed("DT")
  skip_if_not_installed("bsicons")

  e <- new.env(parent = globalenv())
  suppressWarnings(sys.source(.sweep_app_path(), envir = e))

  shiny::testServer(e$server, {
    # formatStyle() interpolates the style as a jQuery .css() property. A
    # function literal there is called with (index, current value), so the
    # style must be an expression over the `value` variable DT declares.
    for (style in list(colour_if("x<0.05"), weight_if("x<0.05"))) {
      js <- as.character(style)
      expect_false(startsWith(js, "function"))
      expect_true(endsWith(js, "(value)"))
    }
    d <- data.frame(item = c("I1", "I2"), infit_ms = c(0.5, 1.0),
                    fit_resid = c(-3, 0.2), p = c(0.001, 0.9),
                    stringsAsFactors = FALSE)
    callback <- num_dt(d, p_bold = "p")$x$options$rowCallback
    expect_false(grepl("css({'color':function", callback, fixed = TRUE))
    expect_false(grepl("css({'font-weight':function", callback, fixed = TRUE))
    expect_true(grepl("var value=data[1];", callback, fixed = TRUE))
    expect_equal(lengths(regmatches(callback,
                                    gregexpr("})(value)", callback,
                                             fixed = TRUE))), 3L)
  })
})

test_that("the app explains a withheld dependence or spread inference", {
  skip_on_cran()
  skip_if_not_installed("shiny")
  skip_if_not_installed("bslib")
  skip_if_not_installed("DT")
  skip_if_not_installed("bsicons")

  e <- new.env(parent = globalenv())
  suppressWarnings(sys.source(.sweep_app_path(), envir = e))

  shiny::testServer(e$server, {
    dep_res(list(dependent = "I02", independent = "I01", d = 0.757,
                 se = NA_real_, t = NA_real_, df = NA_real_, p = NA_real_,
                 note = "the resolved-threshold covariance is unavailable",
                 thresholds = data.frame(k = 1L, d = 0.757)))
    session$flushReact()
    expect_match(output$dep_txt, "descriptive; inference withheld",
                 fixed = TRUE)
    expect_match(output$dep_txt,
                 "Note: the resolved-threshold covariance is unavailable",
                 fixed = TRUE)

    dep_res(list(dependent = "I02", independent = "I01", d = 0.5,
                 se = 0.2, t = 2.5, df = 100, p = 0.014,
                 thresholds = data.frame(k = 1L, d = 0.5)))
    session$flushReact()
    expect_match(output$dep_txt, "(se 0.200), t(100) = 2.50", fixed = TRUE)
    expect_false(grepl("Note:", output$dep_txt, fixed = TRUE))

    # a result carrying no standard error at all: the panel reports the
    # magnitude as descriptive rather than failing on a zero-length test
    dep_res(list(dependent = "I02", independent = "I01", d = 0.5, p = 0.012,
                 thresholds = data.frame(k = 1L, d = 0.5)))
    session$flushReact()
    expect_match(output$dep_txt, "descriptive; inference withheld",
                 fixed = TRUE)

    spread <- data.frame(subtest = "S1", lambda = 0.6, lub = 0.55,
                         se = NA_real_, t = NA_real_, p = NA_real_,
                         p_adj = NA_real_, below_bound = NA, dependent = NA)
    spread_res(structure(spread, class = c("rasch_spread", "data.frame"),
                         alpha = 0.05, p_adjust = "holm",
                         note = "the reference refit is unavailable"))
    session$flushReact()
    expect_match(as.character(output$spread_note$html),
                 "Note: the reference refit is unavailable", fixed = TRUE)

    spread_res(structure(spread, class = c("rasch_spread", "data.frame"),
                         alpha = 0.05, p_adjust = "holm"))
    session$flushReact()
    expect_null(output$spread_note)
  })
})

test_that("the resolution summary withholds an unknown remaining-DIF count", {
  skip_on_cran()
  for (pkg in c("shiny", "bslib", "DT", "bsicons"))
    skip_if_not_installed(pkg)

  e <- new.env(parent = globalenv())
  suppressWarnings(sys.source(.sweep_app_path(), envir = e))

  empty_splits <- data.frame(
    order = integer(), item = character(), factor = character(),
    base_item = character(), eta2 = numeric(), magnitude = numeric(),
    stringsAsFactors = FALSE)

  shiny::testServer(e$server, {
    # resolve_dif() returns NA for n_remaining_dif when no item-term test in
    # the final assessment was estimable. "NA item(s) still flag DIF" would
    # print that withheld verdict as though it were a count.
    resolve_res(list(
      n_splits = 2L, n_remaining_dif = NA_integer_, splits = empty_splits,
      stopped = paste("no further item flagged DIF; no item-term test in the",
                      "final DIF assessment was estimable; remaining DIF is",
                      "unknown")))
    session$flushReact()
    txt <- as.character(output$resolve_summary$html)
    expect_match(txt, "the number of items still flagging DIF is unknown",
                 fixed = TRUE)
    expect_false(grepl("NA item(s)", txt, fixed = TRUE))
    # the qualification resolve_dif() carries on `stopped` is passed through
    expect_match(txt, "remaining DIF is unknown", fixed = TRUE)

    resolve_res(list(n_splits = 1L, n_remaining_dif = 3L,
                     splits = empty_splits,
                     stopped = "no further item flagged DIF"))
    session$flushReact()
    expect_match(
      as.character(output$resolve_summary$html),
      "1 split(s); no further item flagged DIF; 3 item(s) still flag DIF.",
      fixed = TRUE)
  })
})

test_that("the dimensionality panel reads the subset means descriptively", {
  skip_on_cran()
  for (pkg in c("shiny", "bslib", "DT", "bsicons"))
    skip_if_not_installed(pkg)

  f <- rasch(simulate_rasch(60, 6, seed = 11))
  e <- new.env(parent = globalenv())
  suppressWarnings(sys.source(.sweep_app_path(), envir = e))

  # dimensionality_test() no longer returns a paired t-test of the subset
  # means: the two estimates come from the same persons under one model, so
  # the panel reports the difference and says it tests nothing.
  dt <- list(
    split = "manual", first_eigenvalue = 1.8, prop_significant = 0.07,
    ci = c(0.03, 0.14), n = 60L, n_excluded_extreme = 2L,
    multidimensional = FALSE, verdict_method = "fixed-split binomial",
    alpha = 0.05,
    subset_mean_difference = list(mean_difference = 0.123,
                                  sd_difference = 0.456),
    items_positive = c("I1", "I2", "I3"),
    items_negative = c("I4", "I5", "I6"))
  sentence <- paste("Mean difference between subset estimates: 0.123",
                    "(SD 0.456) -- descriptive, not a test of",
                    "unidimensionality")

  shiny::testServer(e$server, {
    fit_val(f)
    session$flushReact()
    # the control observers clear a stale result on the flush that follows a
    # fit, so the computed result is installed after that flush
    dim_computed(dt)
    session$flushReact()
    txt <- paste(output$dim_txt, collapse = "\n")
    expect_match(txt, sentence, fixed = TRUE)
    expect_false(grepl("Paired t-test", txt, fixed = TRUE))
  })

  # the console print method emits the same sentence, so panel and console
  # agree on what the subset means are worth
  console <- paste(utils::capture.output(print(
    structure(dt, class = "rasch_dimensionality_test"))), collapse = "\n")
  expect_match(console, sentence, fixed = TRUE)
})
