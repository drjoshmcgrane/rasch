.round4_app_path <- function(file = "app.R") {
  path <- testthat::test_path("..", "..", "inst", "shiny", file)
  if (!file.exists(path)) path <- system.file("shiny", file, package = "rasch")
  path
}

.round4_app_env <- function(notes = NULL) {
  e <- new.env(parent = globalenv())
  suppressWarnings(sys.source(.round4_app_path(), envir = e))
  if (!is.null(notes))
    e$showNotification <- function(ui, ..., type = "default", duration = NULL) {
      notes$msgs <- c(notes$msgs, paste0("[", type, "] ",
                                         paste(as.character(ui),
                                               collapse = " ")))
      invisible("id")
    }
  e
}

.round4_notes <- function() {
  n <- new.env(parent = emptyenv())
  n$msgs <- character(0)
  n
}

.round4_csv <- function(txt) {
  p <- tempfile(fileext = ".csv")
  writeLines(txt, p)
  list(datapath = p, name = "map.csv", size = nchar(txt), type = "text/csv")
}

# The observeEvent() whose event expression is input$<id>, as the unevaluated
# call, so a test can read the order in which its body reads its reactives.
.round4_observer <- function(id, path = .round4_app_path()) {
  want <- call("$", as.name("input"), as.name(id))
  found <- list()
  walk <- function(x) {
    if (is.call(x) && identical(x[[1L]], as.name("observeEvent")) &&
        identical(x[[2L]], want))
      found[[length(found) + 1L]] <<- x
    if (is.call(x) || is.expression(x) || is.list(x)) {
      els <- as.list(x)
      for (i in seq_along(els))
        if (!identical(els[[i]], quote(expr = ))) walk(els[[i]])
    }
    invisible(NULL)
  }
  walk(parse(path, keep.source = FALSE))
  found
}

# How many times `fname` is called outside tryCatch()'s protected expression.
.round4_uncaught <- function(expr, fname) {
  n <- 0L
  rec <- function(x, caught) {
    if (!is.call(x)) return(invisible(NULL))
    if (is.name(x[[1L]]) && identical(as.character(x[[1L]]), fname) && !caught)
      n <<- n + 1L
    protects <- is.name(x[[1L]]) && identical(as.character(x[[1L]]), "tryCatch")
    els <- as.list(x)[-1L]
    for (i in seq_along(els))
      if (!identical(els[[i]], quote(expr = )))
        rec(els[[i]], caught || (protects && i == 1L))
    invisible(NULL)
  }
  rec(expr, FALSE)
  n
}

test_that("a restored project survives the echo of its own controls", {
  skip_on_cran()
  skip_if_not_installed("shiny")
  skip_if_not_installed("bslib")
  skip_if_not_installed("DT")
  skip_if_not_installed("bsicons")

  e <- .round4_app_env()

  set.seed(4041)
  d <- simulate_rasch(n_persons = 800, n_items = 12, seed = 77)
  X <- as.data.frame(d[, grep("^I", names(d))])
  f <- rasch(X)
  pos <- names(X)[1:6]
  neg <- names(X)[7:12]
  # few replicates keep the test quick; the resolution warning that comes
  # with them is not what is pinned here
  saved <- suppressWarnings(
    dimensionality_test(f, items_positive = pos, items_negative = neg,
                        B = 9, seed = 1, workers = 1))
  magnitude <- dimensionality_magnitude(f, list(pos, neg))
  project <- .seal_app_project(list(
    format = "rasch-shiny-project", schema = 2L, data = X,
    model_type = "rasch", base_fit = f, rasch_steps = list(),
    btl_steps = list(), kept_fits = list(), kept_fit_code = list(),
    simulation = list(),
    settings = list(model_type = "rasch", dim_pos = pos, dim_neg = neg,
                    dim_boot_B = 9, dim_workers = "1", dim_boot_seed = 1,
                    dif_effects = "main", dif_alpha = 0.05),
    resources = list(),
    results = list(subtest = saved, dimension_magnitude = magnitude,
                   dimension_subsets = list(pos = pos, neg = neg))))
  pf <- tempfile(fileext = ".rasch")
  .save_app_project(project, pf)

  shiny::testServer(e$server, {
    session$flushReact()
    session$setInputs(project_file = list(datapath = pf,
                                          name = "analysis.rasch"))
    session$flushReact()
    # the restore ran to completion and reinstated the saved results
    expect_false(isTRUE(restoring_project()))
    expect_false(is.null(restored_subtest()))
    expect_false(is.null(dm_res()))

    # the browser now echoes the controls the restore updated. This project
    # saved a value for each of them that differs from the one the control
    # held, which is the only shape the browser does echo; the shapes where
    # it sends nothing back are in test-fix-app-round5.R. Without the record
    # of what the restore sent, every guarded observer reads this echo as a
    # user change and clears exactly what the project just reinstated.
    session$setInputs(dim_boot_B = 9, dim_workers = "1", dim_boot_seed = 1,
                      pca_component = "1", dim_pos = pos, dim_neg = neg)
    session$flushReact()
    expect_false(is.null(restored_subtest()))
    expect_false(is.null(dm_res()))
    expect_match(paste(output$dim_txt, collapse = "\n"), "Item split: manual",
                 fixed = TRUE)

    # the DIF panel is guarded the same way: a restored DIF bootstrap keeps
    # its place through the echo of the controls it was saved with
    dif_boot_val("a restored DIF bootstrap")
    session$setInputs(dif_effects = "main", dif_alpha = 0.05)
    session$flushReact()
    expect_identical(dif_boot_val(), "a restored DIF bootstrap")
    session$setInputs(dif_effects = "factorial")
    session$flushReact()
    expect_null(dif_boot_val())

    # a change of the user's own, after the echo, still clears the result
    session$setInputs(dim_boot_B = 25)
    session$flushReact()
    expect_null(restored_subtest())
    # with nothing to print, the panel asks for a run again
    expect_match(tryCatch(output$dim_txt, error = function(e)
      conditionMessage(e)), "press Run t-test", fixed = TRUE)
    session$setInputs(pca_component = "2")
    session$flushReact()
    expect_null(dm_res())
  })
})

test_that("an echoed control is recognised whatever mode it comes back in", {
  skip_on_cran()
  skip_if_not_installed("shiny")

  e <- .round4_app_env()
  # the browser returns a number as a number, not necessarily as the mode it
  # was saved in, and a selectize choice as the character vector it sent
  expect_true(e$.same_input_value(19L, 19))
  expect_true(e$.same_input_value(c("I01", "I02"), c("I01", "I02")))
  expect_true(e$.same_input_value("1", 1))
  expect_false(e$.same_input_value(19, 25))
  expect_false(e$.same_input_value(NULL, 19))
  expect_false(e$.same_input_value("I01", c("I01", "I02")))

  # only the controls a restore actually updates are recorded
  restored <- e$.restored_input_values(
    list(dim_boot_B = 9, exp_type_1 = "numeric", not_a_control = 1))
  expect_identical(sort(names(restored)), c("dim_boot_B", "exp_type_1"))
  expect_identical(e$.restored_input_values(list()), list())
})

test_that("the magnitude button reports a t-test the fit refuses", {
  skip_on_cran()
  skip_if_not_installed("shiny")
  skip_if_not_installed("bslib")
  skip_if_not_installed("DT")
  skip_if_not_installed("bsicons")

  notes <- .round4_notes()
  e <- .round4_app_env(notes)

  set.seed(11)
  n <- 300L
  L <- 10L
  th <- rnorm(n)
  d <- seq(-1.5, 1.5, length.out = L)
  X <- matrix(rbinom(n * L, 1, plogis(outer(th, d, "-"))), n, L,
              dimnames = list(NULL, sprintf("I%02d", 1:L)))
  f <- rasch(as.data.frame(X))
  # unequal discriminations: the bootstrap refuses to generate from this fit
  f$disc <- c(rep(1, L - 1L), 1.8)

  shiny::testServer(e$server, {
    session$flushReact()
    fit_val(f); session$flushReact()
    session$setInputs(dim_pos = character(0), dim_neg = character(0),
                      pca_component = "1", dim_boot_B = 19, dim_workers = 1,
                      dim_boot_seed = 1)
    session$flushReact()
    session$setInputs(dim_apply = 1L); session$flushReact()
    expect_false(is.null(dim_refusal()))
    expect_null(dim_subsets())

    # the magnitude run reads the same refused t-test for its split; an
    # observer discards a validation condition, so it must catch it itself
    notes$msgs <- character(0)
    session$setInputs(dm_run = 1L); session$flushReact()
    expect_length(notes$msgs, 1L)
    expect_match(notes$msgs, "^\\[warning\\] No usable subsets: ")
    expect_match(notes$msgs, "equal discriminations", fixed = TRUE)
    expect_null(dm_res())
  })
})

test_that("the split button reports a DIF analysis that refuses", {
  skip_on_cran()
  skip_if_not_installed("shiny")
  skip_if_not_installed("bslib")
  skip_if_not_installed("DT")
  skip_if_not_installed("bsicons")

  notes <- .round4_notes()
  e <- .round4_app_env(notes)

  set.seed(9721)
  n <- 80L
  X <- matrix(rbinom(n * 6, 1, .5), n, 6,
              dimnames = list(NULL, paste0("I", 1:6)))
  factors <- data.frame(a = factor(rep(1:10, length.out = n)),
                        b = factor(rep(1:10, each = 10, length.out = n)))
  f <- rasch(as.data.frame(X), factors = factors)

  shiny::testServer(e$server, {
    session$flushReact()
    fit_val(f); session$flushReact()
    # no item yields an estimable factorial ANOVA, so dif_res() refuses
    session$setInputs(dif_effects = "factorial", dif_alpha = 0.05)
    session$flushReact()
    expect_error(dif_res())
    notes$msgs <- character(0)
    session$setInputs(make_split = 1L); session$flushReact()
    expect_length(notes$msgs, 1L)
    expect_match(notes$msgs, "^\\[error\\] DIF analysis failed: ")
    expect_match(notes$msgs, "estimable factorial ANOVA", fixed = TRUE)
  })
})

test_that("Run Analysis says why an unusable set map stopped it", {
  skip_on_cran()
  skip_if_not_installed("shiny")
  skip_if_not_installed("bslib")
  skip_if_not_installed("DT")
  skip_if_not_installed("bsicons")

  notes <- .round4_notes()
  e <- .round4_app_env(notes)
  d <- as.data.frame(simulate_rasch(n_persons = 120, n_items = 6, seed = 81))
  its <- sprintf("I%02d", 1:6)

  shiny::testServer(e$server, {
    sim_data(d)
    sim_code_val("simulate_rasch(n_persons = 120, n_items = 6, seed = 81)")
    session$flushReact()
    session$setInputs(model_type = "efrm", ef_items = its,
                      ef_group = "group", ef_id = "(none)",
                      ef_sets = .round4_csv("item,set"))
    session$flushReact()
    notes$msgs <- character(0)
    session$setInputs(run = 1L); session$flushReact()
    expect_length(notes$msgs, 1L)
    expect_match(notes$msgs, "^\\[error\\] Analysis failed: ")
    expect_match(notes$msgs, "item-set CSV has no rows", fixed = TRUE)
    expect_null(fit_val())

    # the older refusals in the same reader reach the user too
    session$setInputs(ef_sets = .round4_csv(c("item,set", "I01,A", "I01,B")))
    session$flushReact()
    notes$msgs <- character(0)
    session$setInputs(run = 2L); session$flushReact()
    expect_length(notes$msgs, 1L)
    expect_match(notes$msgs, "more than once", fixed = TRUE)
  })

  # the paired-comparison frames button reads the same map through
  # btlef_build_sets(); its refusal must not be read bare either
  obs <- .round4_observer("btlef_run")
  expect_length(obs, 1L)
  expect_identical(.round4_uncaught(obs[[1L]][[3L]], "btlef_build_sets"), 0L)
  run_obs <- .round4_observer("run")
  expect_length(run_obs, 1L)
  expect_identical(.round4_uncaught(run_obs[[1L]][[3L]], "ef_setmap"), 0L)
})

test_that("the strict export refusal reads as two sentences", {
  skip_on_cran()
  skip_if_not_installed("shiny")
  skip_if_not_installed("bslib")
  skip_if_not_installed("DT")
  skip_if_not_installed("bsicons")

  e <- .round4_app_env(.round4_notes())
  set.seed(5)
  n <- 150L
  X <- matrix(rbinom(n * 8, 1, .5), n, 8,
              dimnames = list(NULL, paste0("I", 1:8)))
  scoring <- rasch(as.data.frame(X))
  # a scoring fit: the package refuses the bootstrap, with no full stop
  scoring$refit_spec$fixed_calibration <- TRUE

  shiny::testServer(e$server, {
    session$flushReact()
    fit_val(scoring); session$flushReact()
    session$setInputs(dim_boot_B = 19, dim_workers = 1, dim_boot_seed = 1,
                      pca_component = "1", dim_pos = character(0),
                      dim_neg = character(0))
    session$flushReact()
    session$setInputs(dim_apply = 1L); session$flushReact()
    err <- tryCatch(app_subtest_res(strict = TRUE),
                    error = function(e) conditionMessage(e))
    expect_match(err, "for a new analysis. Run a supported t-test",
                 fixed = TRUE)
    expect_false(grepl("analysis Run a supported", err, fixed = TRUE))

    # an app validation string already ends in a full stop; it keeps one
    session$setInputs(dim_boot_B = 0); session$flushReact()
    session$setInputs(dim_boot_B = 19); session$flushReact()
    err2 <- tryCatch(app_subtest_res(strict = TRUE),
                     error = function(e) conditionMessage(e))
    expect_match(err2, "press Run t-test. Run a supported t-test",
                 fixed = TRUE)
    expect_false(grepl("t-test.. Run a supported", err2))
  })
})

test_that("the DIF post-hoc caption lists its notes", {
  skip_on_cran()
  skip_if_not_installed("shiny")
  skip_if_not_installed("bslib")
  skip_if_not_installed("DT")
  skip_if_not_installed("bsicons")

  e <- .round4_app_env(.round4_notes())
  set.seed(21)
  n <- 420L
  X <- as.data.frame(simulate_rasch(n_persons = n, n_items = 8,
                                    seed = 21))[, sprintf("I%02d", 1:8)]
  grp <- c(rep("g1", 200), rep("g2", 200), rep("g3", 20))
  # g3 is observed at one site only, so its contrasts are withheld
  site <- c(rep(c("s1", "s2"), 100), rep(c("s1", "s2"), 100), rep("s1", 20))
  X$I02[grp == "g2"] <- pmin(1, X$I02[grp == "g2"] +
                               rbinom(sum(grp == "g2"), 1, .35))
  f <- rasch(X, factors = data.frame(grp = grp, site = site,
                                     stringsAsFactors = FALSE))

  shiny::testServer(e$server, {
    session$flushReact()
    fit_val(f); session$flushReact()
    session$setInputs(dif_effects = "main", dif_alpha = 0.05,
                      dif_size_minn = 30, dif_size_flag = 0.5)
    session$flushReact()
    row <- which(dif_res()$summary$item == "I02" &
                   dif_res()$summary$term == "grp")
    expect_length(row, 1L)
    session$setInputs(dif_tbl_rows_selected = row); session$flushReact()
    ph <- dif_posthoc_res()
    expect_true(length(ph$notes) > 1L)
    html <- output$dif_posthoc_note$html
    expect_match(html, "<li>", fixed = TRUE)
    # the boundary between two withheld contrasts is visible
    expect_false(grepl(paste(ph$notes, collapse = " "), html, fixed = TRUE))
    for (note in ph$notes) expect_match(html, note, fixed = TRUE)
  })
})
