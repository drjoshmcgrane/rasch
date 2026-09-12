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
  expect_identical(r$algorithm, "factor-design-resolution-2")
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

test_that("automatic resolution withholds its verdict when no test was run", {
  # Only four persons hold a complete three-wave panel, so every reported
  # term is a within-person term with no estimable test. dif_anova reports
  # those as NA rather than as an error, and NA is not "no DIF".
  set.seed(4106)
  n <- 120L; L <- 6L
  d <- seq(-1.5, 1.5, length.out = L); th <- rnorm(n)
  draw <- function(theta, shift) {
    sh <- matrix(0, length(theta), L); sh[, 3] <- shift
    matrix(rbinom(length(theta) * L, 1,
                  plogis(outer(theta, d, "-") - sh)), length(theta), L)
  }
  few <- 1:4
  X <- rbind(draw(th, 0), draw(th + .2, 1.6), draw(th[few] + .4, 1.6))
  colnames(X) <- paste0("I", seq_len(L))
  fit <- rasch(data.frame(X, occasion = c(rep("T1", n), rep("T2", n),
                                          rep("T3", length(few)))),
               id = c(seq_len(n), seq_len(n), few), factors = "occasion")
  expect_true(all(is.na(dif_anova(fit)$summary$F_uniform)))

  r <- resolve_dif(fit)
  expect_equal(r$n_splits, 0)
  expect_match(r$stopped, "no item-term test in the final DIF assessment")
  expect_true(is.na(r$n_remaining_dif))
  expect_true(is.na(r$n_nonuniform))
  expect_match(paste(r$notes, collapse = " "),
               "test\\(s\\) were not estimable")
  expect_equal(r$n_untested, nrow(dif_anova(r$fit)$summary))
  expect_output(print(r), "Remaining items with significant DIF: NA")
  # the assessment's own notes reach the printed output
  expect_output(print(r), "final DIF assessment: ")
})

test_that("a split copy's structurally absent term is not counted untested", {
  # Complete three-wave panels, so every pre-split term is estimable. The
  # split copies of I3 each live in one wave only, so their occasion terms
  # are structurally absent rather than lost: they say nothing about the DIF
  # that remains and must not qualify a clean verdict.
  set.seed(4106)
  n <- 120L; L <- 6L
  d <- seq(-1.5, 1.5, length.out = L); th <- rnorm(n)
  draw <- function(theta, shift) {
    sh <- matrix(0, length(theta), L); sh[, 3] <- shift
    matrix(rbinom(length(theta) * L, 1,
                  plogis(outer(theta, d, "-") - sh)), length(theta), L)
  }
  X <- rbind(draw(th, 0), draw(th + .2, 1.6), draw(th + .4, 1.6))
  colnames(X) <- paste0("I", seq_len(L))
  fit <- rasch(data.frame(X, occasion = rep(c("T1", "T2", "T3"), each = n)),
               id = rep(seq_len(n), 3), factors = "occasion")

  r <- resolve_dif(fit)
  expect_equal(r$n_splits, 1)
  expect_equal(r$n_remaining_dif, 0)
  # the only terms the final assessment could not estimate are the three
  # split copies of I3
  s <- dif_anova(r$fit)$summary
  absent <- s$item[is.na(s$p_uniform_adj) & is.na(s$p_nonuniform_adj)]
  expect_identical(sort(absent), sort(paste0("I3 (T", 1:3, ")")))
  expect_equal(r$n_untested, 0L)
  expect_identical(r$stopped, "no significant DIF remains")
  expect_false(any(grepl("Item-term tests not estimable",
                         utils::capture.output(print(r)))))

  # a stop reason of the loop's own is a fact about the loop, not a verdict,
  # and survives untouched when every term was tested
  capped <- resolve_dif(fit, max_splits = 0)
  expect_identical(capped$stopped, "reached the split cap")
  expect_equal(capped$n_untested, 0L)
})

test_that("automatic resolution qualifies a verdict reached on partial tests", {
  # I1 is answered by the lowest-scoring group a and the highest-scoring
  # group b only, so the group is aliased with the class interval for that
  # item alone and its term is not estimable. The item does span both
  # groups, so the question was answerable: silence there is not no DIF.
  set.seed(5)
  n <- 400L; L <- 7L
  d <- seq(-1.5, 1.5, length.out = L); th <- rnorm(n)
  X <- matrix(rbinom(n * L, 1, plogis(outer(th, d, "-"))), n, L)
  colnames(X) <- paste0("I", seq_len(L))
  g <- rep(c("a", "b"), each = n / 2)
  tot <- rowSums(X[, -1])
  keep <- rep(FALSE, n)
  keep[g == "a"][order(tot[g == "a"])[1:60]] <- TRUE
  keep[g == "b"][order(-tot[g == "b"])[1:60]] <- TRUE
  X[!keep, "I1"] <- NA
  fit <- rasch(data.frame(X, grp = g), factors = "grp")
  expect_equal(nlevels(droplevels(factor(g[!is.na(fit$X[, "I1"])]))), 2L)
  expect_true(is.na(dif_anova(fit)$summary$p_uniform_adj[1]))

  r <- resolve_dif(fit)
  expect_equal(r$n_splits, 0)
  expect_equal(r$n_remaining_dif, 0)
  expect_equal(r$n_untested, 1L)
  # the verdict the loop reached is kept, and qualified rather than replaced
  expect_match(r$stopped, "^no significant DIF remains; ")
  expect_match(r$stopped, "1 of 7 item-term test\\(s\\)")
  expect_output(print(r), "Item-term tests not estimable: 1")
})

test_that("a resolution saved under the superseded tag is not restored", {
  d <- simulate_rasch(200, 6, n_groups = 2, seed = 814)
  fit <- rasch(d, factors = "group")
  current <- resolve_dif(fit)
  old <- current
  old$algorithm <- "factor-design-resolution-1"
  project <- .seal_app_project(list(
    format = "rasch-shiny-project", schema = 2L,
    model_type = "rasch", data = as.data.frame(d), base_fit = fit,
    rasch_steps = list(), btl_steps = list(), kept_fits = list(base = fit),
    settings = list(), results = list(resolve = old)))
  path <- tempfile(fileext = ".rasch")
  on.exit(unlink(path), add = TRUE)
  saveRDS(project, path)
  expect_error(.save_app_project(project, path), "superseded calculation")
  expect_warning(restored <- .read_app_project(path),
                 "factor model and reporting rules")
  expect_null(restored$results$resolve)
  expect_identical(restored$base_fit, fit)
  expect_no_error(.validate_app_project(restored))

  project$results$resolve <- current
  project <- .seal_app_project(project)
  expect_no_error(.save_app_project(project, path))
  expect_no_warning(restored <- .read_app_project(path))
  expect_identical(restored$results$resolve, current)
})

test_that("a second factor's test lost to a split is still counted untested", {
  # Splitting I4 by grp leaves grp constant within each copy, which makes
  # the copy's whole joint ANOVA rank deficient and takes its site test
  # with it. The grp rows are structurally absent and are not counted; the
  # site rows are a real question the split destroyed, and are.
  set.seed(3); n <- 500L; L <- 6L
  d <- seq(-1.8, 1.8, length.out = L); th <- rnorm(n)
  g <- rep(c("a", "b"), each = n / 2)
  sh <- matrix(0, n, L); sh[g == "b", 4] <- 1.4
  X <- matrix(rbinom(n * L, 1, plogis(outer(th, d, "-") - sh)), n, L)
  colnames(X) <- paste0("I", seq_len(L))
  fit <- rasch(data.frame(X, grp = g, site = rep(c("x", "y"), length.out = n)),
               factors = c("grp", "site"))

  r <- resolve_dif(fit)
  expect_equal(r$n_splits, 1)
  expect_equal(r$n_remaining_dif, 0)
  s <- dif_anova(r$fit, factors = c("grp", "site"))$summary
  absent <- s[is.na(s$p_uniform_adj) & is.na(s$p_nonuniform_adj), ]
  expect_identical(sort(paste(absent$item, absent$term)),
                   sort(paste(rep(paste0("I4 (", c("a", "b"), ")"), each = 2),
                              c("grp", "site"))))
  expect_equal(r$n_untested, 2L)
  expect_match(r$stopped, "2 of 12 item-term test\\(s\\)")
})
