.fix_app_path <- function(file = "app.R") {
  path <- testthat::test_path("..", "..", "inst", "shiny", file)
  if (!file.exists(path)) path <- system.file("shiny", file, package = "rasch")
  path
}

test_that("the item chi-square caption withholds a probability the fit withholds", {
  skip_on_cran()
  skip_if_not_installed("shiny")
  skip_if_not_installed("bslib")
  skip_if_not_installed("DT")
  skip_if_not_installed("bsicons")

  e <- new.env(parent = globalenv())
  suppressWarnings(sys.source(.fix_app_path(), envir = e))
  strip <- function(h) gsub("\\s+", " ", gsub("<[^>]+>", " ", as.character(h)))

  d <- simulate_rasch(n_persons = 200, n_items = 8, seed = 5)
  X <- d[, grep("^I", names(d))]
  # every person answers twice: chisq_detail() withholds the per-item
  # probability because its reference assumes independent response rows
  f_rep <- rasch(data.frame(id = rep(sprintf("P%03d", 1:100), 2), X),
                 id = "id", items = names(X))
  expect_true(.has_repeated_residual_units(f_rep))
  expect_true(is.na(chisq_detail(f_rep, f_rep$items$item[1])$p))
  f_one <- rasch(X)
  expect_false(.has_repeated_residual_units(f_one))

  shiny::testServer(e$server, {
    session$flushReact()
    fit_val(f_rep); session$flushReact()
    session$setInputs(items_tbl_rows_selected = 1L); session$flushReact()
    cap <- strip(output$chisq_caption$html)
    expect_match(cap, "on [0-9]+ df, probability unavailable;")
    expect_false(grepl("p = NA", cap, fixed = TRUE))
    expect_match(cap, "withheld because person IDs repeat")
    # the caption agrees with the stat box for the same fit
    expect_match(output$fitsum_tbl$html, "probability unavailable")
    expect_false(grepl("p = NA", output$fitsum_tbl$html, fixed = TRUE))

    # an ordinary fit keeps its probability and carries no withholding note
    fit_val(f_one); session$flushReact()
    cap <- strip(output$chisq_caption$html)
    expect_match(cap, "on [0-9]+ df, p (= 0\\.[0-9]{3}|< 0\\.001);")
    expect_false(grepl("probability unavailable", cap, fixed = TRUE))
    expect_false(grepl("withheld", cap, fixed = TRUE))
  })
})

test_that("the predictor-effects explainer states the df rule of each family", {
  h <- new.env(parent = baseenv())
  sys.source(.fix_app_path("help.R"), envir = h)
  txt <- h$APP_HELP[["expl_coef_tbl"]]
  # the one card serves rasch_explanatory and rasch_btl_explanatory alike,
  # so it must not claim a finite reference for every supported fit
  expect_false(grepl("All supported fits", txt, fixed = TRUE))
  expect_false(grepl("repeated-person|judge-clustered", txt))
  expect_match(txt, "Rasch fits use a t reference on independent person units minus one, whether or not identifiers repeat",
               fixed = TRUE)
  expect_match(txt, "comparative judgement fits use judge-cluster degrees of freedom when judges are identified and a normal reference otherwise",
               fixed = TRUE)
})

test_that("the explainer's df rule is the one the explanatory fits apply", {
  skip_on_cran()
  set.seed(4471)
  n <- 60L
  predictors <- data.frame(item = paste0("I", 1:10),
                           operation = rep(0:1, each = 5),
                           format = rep(c("A", "B"), 5),
                           stringsAsFactors = FALSE)
  theta <- rnorm(n)
  delta <- -0.5 + 0.7 * predictors$operation + 0.35 * (predictors$format == "B")
  X <- sapply(delta, function(d) rbinom(n, 1, plogis(theta - d)))
  colnames(X) <- predictors$item
  f0 <- rasch_explanatory(X, predictors, ~ operation + format)
  # the units are the persons carrying information, not the rows of X
  n_units <- f0$est$cluster_support$n
  expect_true(is.finite(n_units) && n_units <= n)
  expect_equal(unique(f0$est$coefficients$df), n_units - 1L)
  d2 <- as.data.frame(rbind(X, X)); d2$id <- rep(sprintf("P%02d", 1:n), 2)
  f2 <- rasch_explanatory(d2, predictors, ~ operation + format, id = "id",
                          items = predictors$item)
  expect_true(isTRUE(f2$est$cluster_support$repeated))
  expect_equal(unique(f2$est$coefficients$df), f2$est$cluster_support$n - 1L)
  expect_true(is.finite(unique(f2$est$coefficients$df)))

  set.seed(1)
  q <- data.frame(object = LETTERS[1:6], domain = rep(0:1, each = 3))
  beta <- setNames(0.8 * q$domain, q$object)
  pr <- t(combn(q$object, 2))
  d <- data.frame(a = rep(pr[, 1], each = 20), b = rep(pr[, 2], each = 20))
  p <- plogis(beta[d$a] - beta[d$b])
  d$winner <- ifelse(runif(nrow(d)) < p, d$a, d$b)
  cj <- btl_explanatory(d, q, ~ domain, "a", "b", winner = "winner")
  expect_true(all(is.infinite(cj$object_coefficients$df)))
  d$judge <- rep(sprintf("J%02d", 1:10), length.out = nrow(d))
  cj_j <- btl_explanatory(d, q, ~ domain, "a", "b", winner = "winner",
                          judge = "judge")
  expect_equal(unique(cj_j$object_coefficients$df), 9L)
})

test_that("the t-test panel blurb routes both splits to the bootstrap for a verdict", {
  # dimensionality_test() withholds the verdict of a fixed split without a
  # bootstrap reference, so the panel must not present the binomial interval
  # as that split's inferential route
  src <- readLines(.fix_app_path(), warn = FALSE)
  from <- grep('title = "Unidimensionality t-test"', src, fixed = TRUE)
  expect_length(from, 1L)
  blurb <- paste(src[from:(from + 12L)], collapse = " ")
  expect_false(grepl("fixed in advance uses the binomial interval", blurb, fixed = TRUE))
  expect_match(blurb, "binomial interval describes", fixed = TRUE)
  expect_match(blurb, "verdict needs bootstrap replicates (B > 0)", fixed = TRUE)
  expect_match(blurb, "whether the split is fixed in advance or residual-derived", fixed = TRUE)
})
