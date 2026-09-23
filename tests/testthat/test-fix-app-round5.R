.round5_app_path <- function(file = "app.R") {
  path <- testthat::test_path("..", "..", "inst", "shiny", file)
  if (!file.exists(path)) path <- system.file("shiny", file, package = "rasch")
  path
}

.round5_notes <- function() {
  n <- new.env(parent = emptyenv())
  n$msgs <- character(0)
  n
}

.round5_app_env <- function(notes = NULL) {
  e <- new.env(parent = globalenv())
  suppressWarnings(sys.source(.round5_app_path(), envir = e))
  if (!is.null(notes))
    e$showNotification <- function(ui, ..., type = "default", duration = NULL) {
      notes$msgs <- c(notes$msgs, paste0("[", type, "] ",
                                         paste(as.character(ui),
                                               collapse = " ")))
      invisible("id")
    }
  e
}

# A sealed project carrying a descriptive t-test on a nominated split, the
# magnitude table computed from it, and the settings it was saved under.
.round5_project <- function(settings = list()) {
  set.seed(4041)
  d <- simulate_rasch(n_persons = 800, n_items = 12, seed = 77)
  X <- as.data.frame(d[, grep("^I", names(d))])
  f <- rasch(X)
  pos <- names(X)[1:6]
  neg <- names(X)[7:12]
  saved <- suppressWarnings(
    dimensionality_test(f, items_positive = pos, items_negative = neg, B = 0))
  magnitude <- dimensionality_magnitude(f, list(pos, neg))
  base <- list(model_type = "rasch", dim_pos = pos, dim_neg = neg,
               dim_boot_B = 0, dim_workers = "1", dim_boot_seed = 1,
               dif_effects = "main", dif_alpha = 0.05)
  base[names(settings)] <- settings
  project <- .seal_app_project(list(
    format = "rasch-shiny-project", schema = 2L, data = X,
    model_type = "rasch", base_fit = f, rasch_steps = list(),
    btl_steps = list(), kept_fits = list(), kept_fit_code = list(),
    simulation = list(), settings = base, resources = list(),
    results = list(subtest = saved, dimension_magnitude = magnitude,
                   dimension_subsets = list(pos = pos, neg = neg))))
  pf <- tempfile(fileext = ".rasch")
  .save_app_project(project, pf)
  list(file = pf, data = X, fit = f, pos = pos, neg = neg,
       subtest = saved, magnitude = magnitude)
}

test_that("a restore that changes no control still clears on the next change", {
  skip_on_cran()
  skip_if_not_installed("shiny")
  skip_if_not_installed("bslib")
  skip_if_not_installed("DT")
  skip_if_not_installed("bsicons")

  e <- .round5_app_env()
  p <- .round5_project()

  shiny::testServer(e$server, {
    session$flushReact()
    # the controls as the browser already shows them: every one the guarded
    # observers watch is on the value this project saved, so the restore's
    # updates change nothing and the browser echoes nothing back for them
    session$setInputs(dim_boot_B = 0, dim_workers = "1", dim_boot_seed = 1,
                      pca_component = "1", dim_pos = character(0),
                      dim_neg = character(0), dif_effects = "main",
                      dif_alpha = 0.05)
    session$flushReact()
    session$setInputs(project_file = list(datapath = p$file,
                                          name = "saved.rasch"))
    session$flushReact()
    expect_false(is.null(restored_subtest()))
    expect_false(is.null(dm_res()))
    dif_boot_val("a restored DIF bootstrap")

    # the user raises the replicate count by hand. The saved t-test is a
    # descriptive one, computed at B = 0, so it must go -- and with no echo
    # to consume, the round-4 record kept it instead.
    session$setInputs(dim_boot_B = 99)
    session$flushReact()
    expect_null(restored_subtest())
    expect_null(dim_computed())
    expect_match(tryCatch(output$dim_txt, error = function(e)
      conditionMessage(e)), "press Run t-test", fixed = TRUE)
    # and the export gate refuses rather than handing the report a result
    # computed under settings the panel no longer shows
    ex <- tryCatch(app_subtest_res(strict = TRUE), error = function(e) e)
    expect_s3_class(ex, "error")

    # the same for a group whose other control is still on its restored
    # value: dif_alpha holds 0.05, and does not mask the change of model
    session$setInputs(dif_effects = "factorial")
    session$flushReact()
    expect_null(dif_boot_val())
  })
})

test_that("a value the restore sent has no privilege later in the session", {
  skip_on_cran()
  skip_if_not_installed("shiny")
  skip_if_not_installed("bslib")
  skip_if_not_installed("DT")
  skip_if_not_installed("bsicons")

  e <- .round5_app_env()
  p <- .round5_project()

  shiny::testServer(e$server, {
    session$flushReact()
    session$setInputs(dim_boot_B = 0, dim_workers = "1", dim_boot_seed = 1,
                      pca_component = "1", dim_pos = character(0),
                      dim_neg = character(0))
    session$flushReact()
    session$setInputs(project_file = list(datapath = p$file,
                                          name = "saved.rasch"))
    session$flushReact()
    expect_false(is.null(restored_subtest()))

    # the user looks at the second residual component: the restored result
    # belongs to the split the project saved, so it goes
    session$setInputs(pca_component = "2")
    session$flushReact()
    expect_null(restored_subtest())
    expect_null(dm_res())
    # and runs the t-test there
    session$setInputs(dim_pos = character(0), dim_neg = character(0))
    session$flushReact()
    session$setInputs(dim_apply = 1L)
    session$flushReact()
    expect_false(is.null(dim_computed()))
    expect_identical(dim_computed()$split, "residual component 2")
    dm_res(p$magnitude)

    # back to the component the restore had sent. That is a change from the
    # component this result was run on, so the result on screen goes with it
    session$setInputs(pca_component = "1")
    session$flushReact()
    expect_null(dim_computed())
    expect_null(dm_res())
    expect_match(tryCatch(output$dim_txt, error = function(e)
      conditionMessage(e)), "press Run t-test", fixed = TRUE)

    # the same for the bootstrap group: a return to the replicate count the
    # restore sent clears the result computed at the count in between
    session$setInputs(dim_boot_B = 25)
    session$flushReact()
    session$setInputs(dim_apply = 1L)
    session$flushReact()
    expect_false(is.null(dim_computed()))
    session$setInputs(dim_boot_B = 0)
    session$flushReact()
    expect_null(dim_computed())
  })
})

test_that("an example dataset selected before a restore survives it", {
  skip_on_cran()
  skip_if_not_installed("shiny")
  skip_if_not_installed("bslib")
  skip_if_not_installed("DT")
  skip_if_not_installed("bsicons")

  notes <- .round5_notes()
  e <- .round5_app_env(notes)
  p <- .round5_project()

  shiny::testServer(e$server, {
    session$flushReact()
    # the user was looking at a built-in example dataset
    session$setInputs(demo_choice = "dich")
    session$flushReact()
    notes$msgs <- character(0)
    session$setInputs(project_file = list(datapath = p$file,
                                          name = "saved.rasch"))
    session$flushReact()
    expect_false(is.null(fit_val()))
    expect_false(is.null(restored_subtest()))

    # the restore sent demo_choice = "none", and the browser echoes it back
    # a round trip later. It is not a change of dataset: nothing is cleared,
    # and nothing is announced as cleared.
    session$setInputs(demo_choice = "none")
    session$flushReact()
    expect_false(is.null(fit_val()))
    expect_false(is.null(restored_subtest()))
    expect_false(is.null(dim_subsets()))
    expect_false(is.null(dm_res()))
    expect_length(grep("example dataset selection changed", notes$msgs), 0L)

    # a dataset the user does select afterwards clears the analysis and says
    # so, as it always did
    session$setInputs(demo_choice = "pcm")
    session$flushReact()
    expect_null(fit_val())
    expect_null(restored_subtest())
    expect_length(grep("example dataset selection changed", notes$msgs), 1L)

    # and so does a return to "none", the value the restore had sent
    fit_val(p$fit)
    session$setInputs(demo_choice = "none")
    session$flushReact()
    expect_null(fit_val())
    expect_length(grep("example dataset selection changed", notes$msgs), 2L)
  })
})

test_that("the change record answers for each control on its own", {
  skip_on_cran()
  skip_if_not_installed("shiny")
  skip_if_not_installed("bslib")
  skip_if_not_installed("DT")
  skip_if_not_installed("bsicons")

  e <- .round5_app_env()

  shiny::testServer(e$server, {
    session$flushReact()
    # a control with no record has changed, whatever it holds
    expect_true(inputs_changed("dif_alpha"))
    session$setInputs(dim_boot_B = 5, dim_workers = "1")
    session$flushReact()
    # the observer that acted on the change recorded what it acted on
    expect_false(inputs_changed(c("dim_boot_B", "dim_workers",
                                  "dim_boot_seed")))
    # an integer that comes back from the browser as a double still matches
    record_input_values(list(dim_boot_B = 5L))
    expect_false(inputs_changed("dim_boot_B"))
    # each control is compared with its own record: one off its recorded
    # value is a change even when its siblings are on theirs
    record_input_values(list(dim_boot_B = 5, dim_workers = "4"))
    expect_true(inputs_changed(c("dim_boot_B", "dim_workers")))
    # and the record has advanced to the values that change was read at
    expect_false(inputs_changed(c("dim_boot_B", "dim_workers")))
  })
})
