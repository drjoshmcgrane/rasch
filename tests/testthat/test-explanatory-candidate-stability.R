test_that("comparative-judgement candidate failures stay in the family", {
  d <- simulate_btl(n_objects = 6, n_judges = 30,
                    reps_per_pair = 80, seed = 51)
  objects <- sort(unique(c(d$object_a, d$object_b)))
  predictors <- data.frame(
    object = objects,
    x = c(1, 0, 0, 0, 0, 0) +
      1e-5 * c(0, -1, -.5, 0, .5, 1))
  fit <- btl_explanatory(d, predictors, ~ x,
                         "object_a", "object_b", winner = "winner",
                         judge = "judge")
  real_refit <- rasch:::.btl_explanatory_refit
  testthat::local_mocked_bindings(
    .btl_explanatory_refit = function(fit, B, relaxations) {
      if (any(grepl("departure[O1]", colnames(B), fixed = TRUE)))
        stop("forced candidate failure")
      real_refit(fit, B, relaxations)
    },
    .package = "rasch")
  dg <- explanatory_diagnostics(fit)
  expect_equal(nrow(dg), length(objects))
  failed <- dg$object == "O1"
  expect_true(!dg$converged[failed])
  expect_true(is.na(dg$p[failed]))
  expect_true(all(is.finite(dg$p_adj[!failed])))
  expect_equal(dg$p_adj[!failed],
               stats::p.adjust(dg$p[!failed], "holm", n = nrow(dg)))
  expect_match(attr(dg, "note"), "failed or")
})

test_that("near-collinear Rasch candidate failures stay in the family", {
  set.seed(4101)
  predictors <- data.frame(
    item = paste0("I", 1:6),
    operation = rep(0:1, each = 3),
    format = rep(c("A", "B"), 3))
  theta <- rnorm(300)
  delta <- 0.75 * predictors$operation +
    0.35 * (predictors$format == "B")
  X <- sapply(delta, function(d) rbinom(length(theta), 1,
                                        plogis(theta - d)))
  colnames(X) <- predictors$item
  fit <- rasch_explanatory(X, predictors, ~ operation + format)
  real_design <- rasch:::.pcml_design
  testthat::local_mocked_bindings(
    .pcml_design = function(X, B, parameter_names = colnames(B), ...) {
      if (any(grepl("departure_location[I1]", colnames(B), fixed = TRUE)))
        stop("forced candidate failure")
      real_design(X, B, parameter_names, ...)
    },
    .package = "rasch")
  dg <- explanatory_diagnostics(fit)
  expect_equal(nrow(dg), ncol(X))
  failed <- dg$item == "I1"
  expect_true(!dg$converged[failed])
  expect_true(is.na(dg$p[failed]))
  expect_true(all(is.finite(dg$p_adj[!failed])))
  expect_equal(dg$p_adj[!failed],
               stats::p.adjust(dg$p[!failed], "holm", n = nrow(dg)))
  expect_match(attr(dg, "note"), "failed or")
})

test_that("residualized comparative-judgement departure is estimable", {
  d <- simulate_btl(n_objects = 6, n_judges = 30,
                    reps_per_pair = 80, seed = 51)
  objects <- sort(unique(c(d$object_a, d$object_b)))
  predictors <- data.frame(
    object = objects,
    x = c(1, 0, 0, 0, 0, 0) +
      1e-5 * c(0, -1, -.5, 0, .5, 1))
  fit <- btl_explanatory(d, predictors, ~ x,
                         "object_a", "object_b", winner = "winner",
                         judge = "judge")
  D <- matrix(-1 / 6, 6, 1,
              dimnames = list(objects, "departure[O1]"))
  D[1, 1] <- 5 / 6
  stable <- rasch:::.explanatory_addition(fit$explanatory$active_B, D)
  expect_equal(qr(cbind(fit$explanatory$active_B, D), tol = 1e-10)$rank,
               qr(cbind(fit$explanatory$active_B, stable), tol = 1e-10)$rank)
  candidate <- rasch:::.btl_explanatory_refit(
    fit, cbind(fit$explanatory$active_B, stable),
    fit$explanatory$relaxations)
  expect_true(candidate$converged)
  expect_true(is.finite(candidate$loglik))
})
