.resolution_interaction_fit <- function(nonuniform = FALSE) {
  set.seed(if (nonuniform) 88024 else 88023)
  n <- if (nonuniform) 2000L else 1600L
  A <- rep(c("a", "b"), each = n / 2)
  B <- rep(rep(c("c", "d"), each = n / 4), 2)
  cell <- ifelse((A == "a") == (B == "c"), 1, -1)
  theta <- rnorm(n)
  eta <- outer(theta, seq(-1, 1, length.out = 8), "-")
  eta[, 1] <- if (nonuniform) exp(.9 * cell) * theta + 1 else
    eta[, 1] - 1.4 * cell
  X <- matrix(rbinom(n * 8, 1, plogis(eta)), n, 8)
  colnames(X) <- paste0("I", seq_len(8))
  rasch(data.frame(X, A, B), factors = c("A", "B"))
}

test_that("automatic resolution uses the requested factorial model throughout", {
  f <- .resolution_interaction_fit()
  before <- dif_anova(f, effects = "factorial")$summary
  planted <- before$item == "I1" & before$term == "A:B"
  expect_true(before$uniform_DIF[planted])
  expect_false(before$nonuniform_DIF[planted])
  additive <- resolve_dif(f)
  expect_identical(additive$effects, "main")
  expect_equal(additive$n_splits, 0)

  r <- resolve_dif(f, effects = "factorial")
  expect_identical(r$algorithm, "factor-design-resolution-1")
  expect_identical(r$effects, "factorial")
  expect_equal(r$n_splits, 1)
  expect_identical(r$splits$item, "I1")
  expect_identical(r$splits$factor, "A:B")
  expect_equal(sum(.split_source_map(r$fit) == "I1"), 4)
  expect_equal(r$n_remaining_dif, 0)
  after <- dif_anova(r$fit, effects = "factorial")$summary
  expect_false(any(after$uniform_DIF | after$nonuniform_DIF))

  # The final assessment must retain the selected model even when no split
  # is permitted by the caller.
  capped <- resolve_dif(f, effects = "factorial", max_splits = 0)
  expect_equal(capped$n_splits, 0)
  expect_equal(capped$n_remaining_dif, length(unique(before$item[
    (before$uniform_DIF | before$nonuniform_DIF) & !before$superseded])))
  expect_identical(unique(capped$dif$factor), "A:B")
  expect_error(resolve_dif(f, effects = "other"), "arg")
})

test_that("factorial automatic resolution leaves non-uniform DIF for review", {
  f <- .resolution_interaction_fit(nonuniform = TRUE)
  before <- dif_anova(f, effects = "factorial")$summary
  planted <- before$item == "I1" & before$term == "A:B"
  expect_true(before$nonuniform_DIF[planted])
  r <- resolve_dif(f, effects = "factorial", max_splits = 1)
  expect_false("I1" %in% r$splits$item)
  expect_true(any(r$dif$item == "I1" & r$dif$nonuniform))
  expect_gt(r$n_remaining_dif, 0)
  expect_match(r$stopped, "non-uniform DIF requires item review")
})

test_that("automatic resolution retains external factor values", {
  d <- simulate_rasch(1000, 8, n_groups = 2,
    dif = list(items = "I04", uniform = 1.4), seed = 1284)
  f <- rasch(d, id = "id")
  factors <- data.frame(group = d$group)
  expect_null(f$factors)
  before <- dif_anova(f, factors = factors)$summary
  expect_true(before$uniform_DIF[before$item == "I04"])
  r <- resolve_dif(f, factors = factors)
  expect_identical(r$splits$item, "I04")
  expect_identical(r$splits$factor, "group")
  expect_equal(r$n_remaining_dif, 0)
  expect_identical(r$fit$person$id, f$person$id)
  expect_equal(nrow(r$fit$X), nrow(f$X))
  vector <- resolve_dif(f, factors = d$group)
  expect_equal(vector$splits, r$splits)

  # An external column with the same name also replaces stored values;
  # choosing by name again would silently test a different grouping.
  d$group <- rev(d$group)
  stored <- rasch(d, id = "id", factors = "group")
  overridden <- resolve_dif(stored, factors = factors)
  expect_equal(overridden$splits, r$splits)
  expect_identical(overridden$fit$X, r$fit$X)
  expect_error(resolve_dif(f, factors = factors[, FALSE, drop = FALSE]),
               "at least one person factor")
})

test_that("the default additive resolver is unchanged on a null design", {
  d <- simulate_rasch(500, 8, n_groups = 2, seed = 671)
  f <- rasch(d, id = "id", factors = "group")
  old <- resolve_dif(f)
  explicit <- resolve_dif(f, effects = "main")
  expect_identical(old, explicit)
  expect_equal(old$n_splits, 0)
  expect_equal(old$n_remaining_dif, 0)
})

test_that("app automatic DIF resolution and its frozen code keep factorial effects", {
  skip_on_cran()
  for (pkg in c("shiny", "bslib", "DT", "bsicons"))
    skip_if_not_installed(pkg)
  f <- .resolution_interaction_fit()
  e <- new.env(parent = globalenv())
  app_path <- test_path("..", "..", "inst", "shiny", "app.R")
  if (!file.exists(app_path))
    app_path <- system.file("shiny", "app.R", package = "rasch")
  suppressWarnings(sys.source(
    app_path, envir = e))
  shiny::testServer(e$server, {
    fit_val(f)
    session$setInputs(dif_effects = "factorial", dif_alpha = .05)
    session$flushReact()
    session$setInputs(resolve_all = 1)
    session$flushReact()
    expect_equal(resolve_res()$n_splits, 1)
    expect_identical(resolve_res()$run_effects, "factorial")
    expect_identical(active_step()$type, "dif_auto")
    code <- active_step()$code
    expect_match(code, 'effects = "factorial"', fixed = TRUE)
    env <- new.env(parent = globalenv())
    env$fit <- f
    eval(parse(text = code), envir = env)
    expect_identical(env$fit$X, fit()$X)
    expect_equal(env$dif_resolution$splits, resolve_res()$splits)
    session$setInputs(dif_effects = "main")
    session$flushReact()
    expect_identical(active_step()$code, code)
  })
})
