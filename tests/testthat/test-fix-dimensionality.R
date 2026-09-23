# The B = 0 verdict note points a fit at the bootstrap only when the
# bootstrap will accept that fit; otherwise it says the comparison stays
# descriptive and why, so the reader is not sent down a route that errors.

test_that("the B = 0 note points a refittable fit at the bootstrap", {
  f <- rasch(simulate_rasch(300, 12, seed = 926))
  z <- dimensionality_test(f, items_positive = sprintf("I%02d", 1:6),
                           items_negative = sprintf("I%02d", 7:12))
  expect_true(is.na(z$multidimensional))
  expect_identical(z$verdict_method,
                   "withheld for fixed split without bootstrap reference")
  expect_match(z$verdict_note, "targeting")
  expect_match(z$verdict_note, "Use B > 0 for a model-based reference\\.$")
  expect_false(grepl("stays descriptive", z$verdict_note, fixed = TRUE))
  expect_output(print(z), "withheld without a bootstrap reference")
  expect_output(print(z), "Use B > 0", fixed = TRUE)
})

test_that("the B = 0 note names the refusal for a fit the bootstrap rejects", {
  set.seed(1)
  q <- data.frame(item = paste0("I", 1:8), operation = rep(0:1, each = 4),
                  format = rep(c("A", "B"), 4))
  difficulty <- -1 + 0.7 * q$operation + 0.4 * (q$format == "B")
  X <- matrix(rbinom(500 * 8, 1, plogis(outer(rnorm(500), difficulty, "-"))),
              500, 8)
  colnames(X) <- q$item
  fit <- rasch_explanatory(X, predictors = q, formula = ~ operation + format)
  args <- list(fit = fit, items_positive = paste0("I", 1:4),
               items_negative = paste0("I", 5:8))
  z <- do.call(dimensionality_test, args)
  expect_true(is.na(z$multidimensional))
  expect_identical(z$verdict_method,
                   "withheld for fixed split without bootstrap reference")
  expect_match(z$verdict_note, "targeting")
  expect_false(grepl("Use B > 0", z$verdict_note, fixed = TRUE))
  expect_match(z$verdict_note, "No bootstrap reference is available for this fit")
  expect_match(z$verdict_note, "stays descriptive\\.$")
  # the reason quoted in the note is the refusal B > 0 raises
  refusal <- tryCatch(do.call(dimensionality_test,
                              c(args, list(B = 19, workers = 1))),
                      rasch_refusal = function(e) conditionMessage(e))
  expect_type(refusal, "character")
  expect_true(grepl(refusal, z$verdict_note, fixed = TRUE))
  expect_output(print(z), "stays descriptive", fixed = TRUE)
})

test_that("the note advice follows the bootstrap check for other refusals", {
  f <- rasch(simulate_rasch(300, 12, seed = 926))
  expect_identical(.dim_reference_advice(f),
                   "Use B > 0 for a model-based reference.")
  # unequal frame units are refused by the bootstrap and named in the note
  g <- f
  g$disc <- rep(c(1, 1.5), length.out = ncol(f$X))
  expect_match(.dim_reference_advice(g), "equal discriminations")
  expect_match(.dim_reference_advice(g), "^No bootstrap reference")
})
