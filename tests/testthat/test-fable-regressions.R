test_that("uniform and non-uniform availability are counted separately", {
  set.seed(1030)
  L <- sample(5:8, 1); K <- sample(2:4, 1)
  ns <- sample(c(15, 25, 40, 80, 150), K, replace = TRUE); n <- sum(ns)
  g <- factor(rep(letters[1:K], ns)); mu <- rnorm(K, 0, 1.5)
  th <- rnorm(n, rep(mu, ns), .8); d <- seq(-2, 2, length.out = L)
  sh <- matrix(0, n, L)
  sh[g == levels(g)[K], sample(L, 1)] <- runif(1, -2, 2)
  X <- matrix(rbinom(n * L, 1, plogis(outer(th, d, "-") - sh)), n, L)
  colnames(X) <- paste0("I", 1:L)
  it <- sample(L, 1); who <- sample(n, floor(n * runif(1, .3, .8)))
  X[who, it] <- NA; ng <- sample(2:5, 1)
  fit <- rasch(data.frame(X, grp = g), factors = "grp", n_groups = ng)
  s <- dif_anova(fit)$summary
  expect_true(any(is.finite(s$p_uniform_adj) & is.na(s$p_nonuniform_adj)))
  r <- resolve_dif(fit, max_splits = 0)
  expect_equal(r$n_untested,
    sum(!is.finite(s$p_uniform_adj)) + sum(!is.finite(s$p_nonuniform_adj)))
  expect_true(is.na(r$n_nonuniform))
  expect_output(print(r), "Non-uniform DIF count unavailable")
})

test_that("a structurally absent factorial contrast is not a lost DIF test", {
  set.seed(12); n <- 800; L <- 7
  g <- rep(c("a", "b"), each = n / 2); s <- rep(c("x", "y"), length.out = n)
  X <- matrix(rbinom(n * L, 1, plogis(outer(rnorm(n),
    seq(-1.6, 1.6, length.out = L), "-"))), n, L)
  colnames(X) <- paste0("I", 1:L)
  X[g == "b" & s == "y", "I1"] <- NA
  fit <- rasch(data.frame(X, grp = g, site = s), factors = c("grp", "site"))
  r <- resolve_dif(fit, effects = "factorial", max_splits = 0)
  a <- dif_anova(fit, effects = "factorial")
  row <- a$summary$item == "I1" & a$summary$term == "grp:site"
  expect_equal(sum(row), 1)
  expect_true(all(is.na(a$summary$p_uniform_adj[row])))
  expect_match(paste(a$notes, collapse = " "), "I1 \\[grp:site\\]: unavailable")
  expect_match(paste(a$notes, collapse = " "), "retained in the holm adjustment family")
  # The two main terms remain questions; the interaction is structurally absent.
  expected <- sum(!is.finite(a$summary$p_uniform_adj[!row])) +
    sum(!is.finite(a$summary$p_nonuniform_adj[!row]))
  expect_equal(r$n_untested, expected)
})

test_that("automatic resolution retains an unexpected split refusal", {
  set.seed(56); n <- 700L
  g <- rep(c("a", "b"), each = n / 2)
  eta <- outer(rnorm(n), seq(-1.5, 1.5, length.out = 8), "-")
  eta[g == "b", 3] <- eta[g == "b", 3] - 2
  X <- matrix(rbinom(n * 8, 1, plogis(eta)), n, 8,
              dimnames = list(NULL, paste0("I", 1:8)))
  fit <- rasch(data.frame(X, grp = g), factors = "grp")
  expect_true(any(dif_anova(fit)$summary$uniform_DIF))
  z <- with_mocked_bindings(resolve_dif(fit),
    .dif_resolve = function(...) stop("calibration link disconnected"),
    .package = "rasch")
  expect_match(paste(z$notes, collapse = " "), "calibration link disconnected")
  expect_false(any(grepl("weak boundary", z$notes, fixed = TRUE)))
})

test_that("DIF withholding retains thin-cell and refit explanations", {
  d <- simulate_rasch(200, 6, n_groups = 2, seed = 918)
  f <- rasch(d, id = "id", factors = "group")
  z <- dif_contrasts(f, items = "I01", min_n = 1000)
  expect_true(all(is.na(z$table$estimate)))
  expect_match(paste(z$notes, collapse = " "), "min_n|fewer than")
  expect_error(dif_posthoc(f, "I01", term = "group", min_n = 1000),
               "fewer than")
  z <- with_mocked_bindings(
    dif_contrasts(f, items = "I01"),
    split_items = function(...) stop("specific refit failure"), .package = "rasch")
  expect_match(paste(z$notes, collapse = " "), "specific refit failure")
  expect_error(with_mocked_bindings(dif_posthoc(f, "I01", term = "group"),
    split_items = function(...) stop("specific refit failure"), .package = "rasch"),
    "specific refit failure")
  m <- simulate_mfrm(60, 4, 3, n_categories = 3, seed = 17)
  fm <- rasch_mfrm(m, person = "person", item = "item", score = "score", facets = "rater")
  grp <- factor(rep(c("a", "b"), length.out = nrow(fm$X)))
  item <- unique(fm$virtual_map$item)[1]
  z <- with_mocked_bindings(.dif_resolve(fm, item, grp, 2),
    split_items = function(...) stop("specific MFRM failure"), .package = "rasch")
  expect_identical(z$refit_error, "specific MFRM failure")
  expect_match(paste(z$notes, collapse = " "), "virtual-item level")
  expect_match(paste(z$notes, collapse = " "), "specific MFRM failure")
})

test_that("colliding numeric factor labels have an actionable contrast error", {
  f <- factor(c("1", "01", "2"), levels = c("1", "01", "2"))
  expect_error(.dif_factor_contrasts(f, "visit"), "colliding numeric scores.*1, 01")
  expect_error(.dif_leading(f), "supply explicit contrasts")
  expect_length(.dif_factor_contrasts(factor(c("1", "2", "3")), "visit"), 2)
})

test_that("fixed subset verdicts need a usable bootstrap reference", {
  f <- rasch(simulate_rasch(300, 12, seed = 926))
  args <- list(fit = f, items_positive = sprintf("I%02d", 1:6),
               items_negative = sprintf("I%02d", 7:12))
  z <- do.call(dimensionality_test, args)
  expect_true(is.na(z$multidimensional))
  expect_identical(z$verdict_method,
                   "withheld for fixed split without bootstrap reference")
  expect_match(z$verdict_note, "targeting")
  expect_true(is.logical(z$binomial_multidimensional))
  expect_output(print(z), "withheld without a bootstrap reference")
  expect_warning(b <- do.call(dimensionality_test,
    c(args, list(B = 9, seed = 927, workers = 1))), "smallest attainable")
  expect_true(is.na(b$multidimensional))
  expect_match(b$verdict_method, "resolution insufficient")
})

test_that("legacy truth follows unique response keys, never unchecked position", {
  d <- simulate_rasch(300, 30, model = "PCM", n_categories = 3,
                       difficulty = c(-1, 1), seed = 14)
  attr(d, "truth")$person_id <- NULL
  expect_null(names(attr(d, "truth")$theta))
  items <- names(attr(d, "truth")$difficulty)
  ord <- rev(seq_len(nrow(d)))
  X <- as.matrix(d[ord, items]); rownames(X) <- NULL
  f <- rasch(X)
  expect_identical(as.character(f$person$id), as.character(seq_len(nrow(d))))
  z <- sim_recovery(f, d)$pieces[["person ability"]]
  expect_equal(z$true, unname(attr(d, "truth")$theta[ord] - mean(attr(d, "truth")$theta)))
  d2 <- simulate_rasch(200, 5, seed = 15)
  attr(d2, "truth")$person_id <- NULL
  X2 <- as.matrix(d2[rev(seq_len(nrow(d2))), names(attr(d2, "truth")$difficulty)])
  rownames(X2) <- NULL
  z2 <- sim_recovery(rasch(X2), d2)
  expect_null(z2$pieces[["person ability"]])
  expect_match(z2$note, "do not uniquely pair")
})

test_that("app restoration clears computed dimensionality and displays DIF notes", {
  skip_if_not_installed("shiny"); skip_if_not_installed("bslib")
  skip_if_not_installed("DT"); skip_if_not_installed("bsicons")
  app <- testthat::test_path("..", "..", "inst", "shiny", "app.R")
  if (!file.exists(app)) app <- system.file("shiny", "app.R", package = "rasch")
  e <- new.env(parent = globalenv()); suppressWarnings(sys.source(app, e))
  d <- simulate_rasch(120, 6, seed = 81); f <- rasch(d, id = "id")
  project <- .seal_app_project(list(format = "rasch-shiny-project", schema = 2L,
    data = as.data.frame(d), model_type = "rasch", base_fit = f,
    rasch_steps = list(), btl_steps = list(), settings = list(), results = list()))
  path <- tempfile(fileext = ".rasch"); bank <- tempfile(fileext = ".csv")
  on.exit(unlink(c(path, bank)), add = TRUE)
  .save_app_project(project, path)
  write.csv(data.frame(object = character(), location = numeric()), bank, row.names = FALSE)
  shiny::testServer(e$server, {
    session$flushReact()
    resolve_res(list(notes = "specific split failure", n_remaining_dif = NA,
                     n_splits = 0, stopped = "withheld", splits = data.frame()))
    session$flushReact()
    expect_match(output$resolve_notes$html, "specific split failure")
    dim_computed(list(stale = TRUE))
    session$setInputs(project_file = list(datapath = path, name = "test.rasch"))
    session$flushReact()
    expect_null(dim_computed())
    session$setInputs(bt_eq_file = list(datapath = bank, name = "bank.csv"))
    session$flushReact()
    expect_null(bt_eq_bank())
    no_p <- f; no_p$total_chisq_p <- NA_real_
    fit_val(no_p)
    session$flushReact()
    expect_match(output$fitsum_tbl$html, "probability unavailable")
    expect_false(grepl("p = NA", output$fitsum_tbl$html, fixed = TRUE))
  })
})
