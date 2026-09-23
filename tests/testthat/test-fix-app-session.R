.fix_session_app_path <- function(file = "app.R") {
  path <- testthat::test_path("..", "..", "inst", "shiny", file)
  if (!file.exists(path)) path <- system.file("shiny", file, package = "rasch")
  path
}

# Every `session$onFlush()` / `session$onFlushed()` registration in the app,
# as the unevaluated call, so a test can both inspect and invoke the callback
# outside the reactive context a real ShinySession does not supply.
.fix_session_flush_calls <- function(path = .fix_session_app_path()) {
  found <- list()
  walk <- function(x) {
    if (is.call(x)) {
      fn <- x[[1L]]
      if (is.call(fn) && identical(fn[[1L]], as.name("$")) &&
          identical(fn[[2L]], as.name("session")) &&
          as.character(fn[[3L]]) %in% c("onFlush", "onFlushed"))
        found[[length(found) + 1L]] <<- x
    }
    if (is.call(x) || is.pairlist(x) || is.expression(x) || is.list(x)) {
      els <- as.list(x)
      # a call may carry empty arguments; those are not expressions to walk
      for (i in seq_along(els))
        if (!identical(els[[i]], quote(expr = ))) walk(els[[i]])
    }
    invisible(NULL)
  }
  walk(parse(path, keep.source = FALSE))
  found
}

test_that("a flushed callback isolates the reactives it reads", {
  skip_on_cran()
  skip_if_not_installed("shiny")

  calls <- .fix_session_flush_calls()
  # one registration today: the project-open observer's settings restore
  expect_length(calls, 1L)

  for (cl in calls) {
    # the registered function must do nothing outside isolate(): a flushed
    # callback runs in the session's reactive domain but in no reactive
    # context, so an unisolated reactiveVal read aborts the R process
    f <- cl[[2L]]
    expect_true(is.call(f) && identical(f[[1L]], as.name("function")))
    body_expr <- f[[3L]]
    if (identical(body_expr[[1L]], as.name("{"))) {
      exprs <- as.list(body_expr)[-1L]
      expect_length(exprs, 1L)
      body_expr <- exprs[[1L]]
    }
    expect_identical(body_expr[[1L]], as.name("isolate"))
  }

  # and the property itself: invoke the callback with real reactiveVals, from
  # no reactive context, the way ShinySession$flushOutput() invokes it
  restored <- NULL
  env <- new.env(parent = asNamespace("rasch"))
  env$session <- structure(list(), class = "rasch_test_session")
  env$isolate <- shiny::isolate
  env$updateSelectizeInput <- function(...) invisible(NULL)
  env$.restore_app_settings <- function(session, settings) {
    restored <<- settings
    invisible(NULL)
  }
  env$restored_project_settings <- shiny::reactiveVal(list(dim_boot_B = 19))
  env$dim_subsets <- shiny::reactiveVal(list(pos = "I1", neg = "I2"))
  env$restoring_project <- shiny::reactiveVal(TRUE)
  recorded <- NULL
  env$record_input_values <- function(values) {
    recorded <<- values
    invisible(NULL)
  }
  env$.restored_input_values <- function(settings) settings

  cb <- eval(calls[[1L]][[2L]], envir = env)
  err <- tryCatch({ cb(); NULL }, error = function(e) conditionMessage(e))
  expect_null(err)
  # the saved settings reached the restore, and the open flag was cleared
  expect_identical(restored, list(dim_boot_B = 19))
  expect_false(shiny::isolate(env$restoring_project()))
  # the updated controls echo back later, so the callback records the value
  # it sent to each of them: the saved settings, the example selection the
  # restore resets, and the subsets it reflects in the two selectors
  expect_identical(recorded,
                   list(demo_choice = "none", dim_boot_B = 19,
                        dim_pos = "I1", dim_neg = "I2"))
})

test_that("Run t-test reports a bootstrap the fit refuses", {
  skip_on_cran()
  skip_if_not_installed("shiny")
  skip_if_not_installed("bslib")
  skip_if_not_installed("DT")
  skip_if_not_installed("bsicons")

  e <- new.env(parent = globalenv())
  suppressWarnings(sys.source(.fix_session_app_path(), envir = e))
  notes <- character(0)
  e$showNotification <- function(ui, ..., type = "default") {
    notes <<- c(notes, paste0("[", type, "] ",
                             paste(as.character(ui), collapse = " ")))
    invisible(NULL)
  }

  set.seed(4471)
  n <- 200L
  predictors <- data.frame(item = paste0("I", 1:10),
                           operation = rep(0:1, each = 5),
                           format = rep(c("A", "B"), 5),
                           stringsAsFactors = FALSE)
  theta <- rnorm(n)
  delta <- -0.5 + 0.7 * predictors$operation + 0.35 * (predictors$format == "B")
  X <- sapply(delta, function(d) rbinom(n, 1, plogis(theta - d)))
  colnames(X) <- predictors$item
  f_expl <- rasch_explanatory(X, predictors, ~ operation + format)
  f_ord <- rasch(X)

  shiny::testServer(e$server, {
    session$flushReact()
    fit_val(f_expl); session$flushReact()
    session$setInputs(dim_pos = paste0("I", 1:5), dim_neg = paste0("I", 6:10),
                      dim_boot_B = 19, dim_workers = 1, dim_boot_seed = 1)
    session$flushReact()
    notes <<- character(0)
    session$setInputs(dim_apply = 1L); session$flushReact()
    # the refusal arrives as a notification, not as silence
    expect_length(notes, 1L)
    expect_match(notes, "^\\[warning\\] The t-test was not run: ")
    expect_match(notes, "single-facet Rasch model", fixed = TRUE)
    # and the panel prints its reason instead of asking for another press
    txt <- tryCatch(output$dim_txt, error = function(e) conditionMessage(e))
    expect_match(txt, "single-facet Rasch model", fixed = TRUE)
    expect_false(grepl("press Run t-test", txt, fixed = TRUE))
    # the Export path names it too
    err <- tryCatch(app_subtest_res(strict = TRUE),
                    error = function(e) conditionMessage(e))
    expect_match(err, "dimensionality t-test is unavailable", fixed = TRUE)
    expect_match(err, "single-facet Rasch model", fixed = TRUE)
    # the refusal is held beside the result, never as the result
    expect_null(dim_computed())
    expect_match(dim_refusal(), "single-facet Rasch model", fixed = TRUE)
    # and it clears with the settings that provoked it
    session$setInputs(dim_boot_B = 0); session$flushReact()
    expect_null(dim_refusal())

    # a fit the bootstrap does generate from is unaffected
    session$setInputs(dim_boot_B = 19); session$flushReact()
    fit_val(f_ord); session$flushReact()
    notes <<- character(0)
    session$setInputs(dim_apply = 2L); session$flushReact()
    expect_match(notes, "^\\[message\\] Ran the t-test on the nominated item subsets\\.")
    expect_null(dim_refusal())
    expect_false(is.null(dim_computed()))
  })
})

test_that("the t-test panel names the fits the bootstrap refuses", {
  # the blurb and the replicates label send every fit to B > 0; a fit whose
  # class has no bootstrap reference must not be promised a verdict
  src <- readLines(.fix_session_app_path(), warn = FALSE)
  # the panel text is a paste() of one fragment per source line; rejoin them
  # so an assertion can read the sentence the user reads
  rejoin <- function(x) gsub('",\\s*"', " ", paste(x, collapse = " "))
  from <- grep('title = "Unidimensionality t-test"', src, fixed = TRUE)
  expect_length(from, 1L)
  blurb <- rejoin(src[from:(from + 12L)])
  expect_match(blurb, "The bootstrap generates from a single-facet Rasch model",
               fixed = TRUE)
  expect_match(blurb, "explanatory, many-facet or extended-frame fit",
               fixed = TRUE)
  expect_match(blurb, "the comparison stays descriptive", fixed = TRUE)

  at <- grep('numericInput("dim_boot_B", info_label(', src, fixed = TRUE)
  expect_length(at, 1L)
  label <- rejoin(src[at:(at + 6L)])
  expect_match(label, "Bootstrap replicates", fixed = TRUE)
  expect_match(label, "refuses a positive value and says why", fixed = TRUE)
})

test_that("a withheld DIF test is not reported as a null result", {
  skip_on_cran()
  skip_if_not_installed("shiny")
  skip_if_not_installed("bslib")
  skip_if_not_installed("DT")
  skip_if_not_installed("bsicons")

  e <- new.env(parent = globalenv())
  suppressWarnings(sys.source(.fix_session_app_path(), envir = e))
  notes <- character(0)
  e$showNotification <- function(ui, ..., type = "default") {
    notes <<- c(notes, paste0("[", type, "] ",
                             paste(as.character(ui), collapse = " ")))
    invisible(NULL)
  }
  strip <- function(h) gsub("\\s+", " ", gsub("<[^>]+>", " ", as.character(h)))

  # I2 is answered by group A only below the median score on the other items
  # and by group B only above it, so grp is aliased with the class interval
  # for that item and neither of its DIF tests is estimable
  set.seed(3); n <- 600L; L <- 8L
  d <- seq(-2, 2, length.out = L); th <- rnorm(n)
  g <- rep(c("A", "B"), each = n / 2)
  X <- matrix(rbinom(n * L, 1, plogis(outer(th, d, "-"))), n, L)
  colnames(X) <- paste0("I", 1:L)
  raw <- rowSums(X[, -2]); med <- median(raw)
  X[g == "A" & raw > med, "I2"] <- NA
  X[g == "B" & raw <= med, "I2"] <- NA
  f <- rasch(as.data.frame(X),
             factors = data.frame(grp = g, stringsAsFactors = FALSE))
  da <- dif_anova(f)
  row2 <- which(da$summary$item == "I2")
  expect_length(row2, 1L)
  expect_true(is.na(da$summary$p_uniform_adj[row2]))
  expect_true(is.na(da$summary$p_nonuniform_adj[row2]))
  # the flags read FALSE, which is the trap the caption used to fall into
  expect_false(da$summary$uniform_DIF[row2])
  expect_false(da$summary$nonuniform_DIF[row2])

  shiny::testServer(e$server, {
    session$flushReact()
    fit_val(f); session$flushReact()
    session$setInputs(dif_tbl_rows_selected = row2); session$flushReact()
    expect_identical(dif_sel_item(), "I2")
    note <- strip(output$dif_posthoc_note$html)
    expect_match(note, "untested for this item, not non-significant",
                 fixed = TRUE)
    expect_false(grepl("The omnibus term is not significant", note,
                       fixed = TRUE))

    # a manual split on the same item warns that discrimination is untested
    notes <<- character(0)
    session$setInputs(make_split = 1L); session$flushReact()
    expect_true(any(grepl("non-uniform term is untested for this item", notes,
                          fixed = TRUE)))

    # an item whose tests did run keeps the plain reading
    ok <- which(da$summary$item == "I5")
    session$setInputs(dif_tbl_rows_selected = ok); session$flushReact()
    note_ok <- strip(output$dif_posthoc_note$html)
    expect_false(grepl("untested", note_ok, fixed = TRUE))
    notes <<- character(0)
    session$setInputs(make_split = 2L); session$flushReact()
    expect_false(any(grepl("untested", notes, fixed = TRUE)))
  })
})
