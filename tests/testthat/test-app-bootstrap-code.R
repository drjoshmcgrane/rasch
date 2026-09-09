test_that("bootstrap table code reproduces the app's complete exported tables", {
  skip_on_cran()
  for (pkg in c("shiny", "bslib", "DT", "bsicons"))
    skip_if_not_installed(pkg)
  f <- rasch(simulate_rasch(80, 6, seed = 779))
  bs <- suppressWarnings(fit_bootstrap(f, B = 99, workers = 1, seed = 72))
  e <- new.env(parent = globalenv())
  app_path <- test_path("..", "..", "inst", "shiny", "app.R")
  if (!file.exists(app_path))
    app_path <- system.file("shiny", "app.R", package = "rasch")
  suppressWarnings(sys.source(
    app_path, envir = e))
  shiny::testServer(e$server, {
    fit_val(f)
    session$flushReact()
    boot_val(list(bs = bs, B = 99, seed = 72, kind = "rasch"))
    session$flushReact()
    env <- new.env(parent = globalenv())
    env$fit <- f
    persons <- suppressWarnings(eval(parse(text = output$person_tbl_code),
                                      envir = env))
    expect_equal(anyDuplicated(names(persons)), 0L)
    expect_identical(names(persons), names(persons_with_boot()))
    expect_equal(as.data.frame(persons), as.data.frame(persons_with_boot()))
    items <- suppressWarnings(eval(parse(text = output$items_tbl_code),
                                    envir = env))
    expect_true("fit_resid_p_boot" %in% names(items))
    expect_identical(names(items), names(items_with_boot()))
    expect_equal(as.data.frame(items), as.data.frame(items_with_boot()))
  })
})

test_that("installed app binds package internals used by reactive paths", {
  for (pkg in c("shiny", "bslib", "DT", "bsicons"))
    skip_if_not_installed(pkg)
  app_path <- test_path("..", "..", "inst", "shiny", "app.R")
  if (!file.exists(app_path))
    app_path <- system.file("shiny", "app.R", package = "rasch")
  e <- new.env(parent = globalenv())
  suppressWarnings(sys.source(app_path, envir = e))
  for (nm in c(".has_repeated_residual_units", ".fit_boot_signature",
               ".fit_boot_signature_matches", ".fit_boot_md5",
               ".sim_planted_count"))
    expect_true(is.function(get(nm, envir = e, inherits = FALSE)), nm)
})
