.joint_dif_fixture <- function(n = 240L, seed = 9151) {
  set.seed(seed)
  A <- factor(rep(c("a0", "a1"), each = n / 2))
  B <- factor(ifelse(runif(n) < ifelse(A == "a0", .15, .85), "b1", "b0"))
  ci <- factor(rep(1:3, length.out = n))
  ids <- c(seq_len(n), which(A == "a1"))
  occasion <- factor(c(rep("pre", n), rep("post", sum(A == "a1"))))
  random_person <- rnorm(n, sd = .4)
  data.frame(pid = ids, A = A[ids], B = B[ids], ci = ci[ids], occasion,
    z = 2 * (B[ids] == "b1") + .7 * (occasion == "post") +
      random_person[ids] + rnorm(length(ids), sd = .3))
}

test_that("incomplete DIF CR3 agrees with independent delete-person fits", {
  d <- .joint_dif_fixture()
  terms <- c("A", "B", "ci", "occasion", "A:ci", "B:ci", "occasion:ci")
  w <- 1 / as.numeric(table(d$pid)[as.character(d$pid)])
  got <- .dif_type2(d, terms, variance = "cr3", cluster = d$pid,
                    weights = w, report_terms = c("A", "A:ci"))
  fit <- lm(z ~ B * ci + occasion * ci + A, d, weights = w)
  beta <- coef(fit)
  deltas <- vapply(unique(d$pid), function(id) {
    keep <- d$pid != id
    coef(lm(z ~ B * ci + occasion * ci + A, d[keep, ], weights = w[keep])) - beta
  }, beta)
  V <- tcrossprod(deltas)
  F_A <- beta["Aa1"]^2 / V["Aa1", "Aa1"]
  full <- lm(z ~ (A + B + occasion) * ci, d, weights = w)
  row <- got[got$term == "A", ]
  expect_equal(row$F_value, unname(F_A), tolerance = 1e-9)
  expect_equal(row$df_denom, length(unique(d$pid)) - full$rank)
  expect_equal(row$p, unname(pf(F_A, 1, row$df_denom, lower.tail = FALSE)))
  expect_gt(row$p, .05)
  reordered <- .dif_type2(d, rev(terms), variance = "cr3", cluster = d$pid,
    weights = w, report_terms = c("A", "A:ci"))
  expect_equal(got$F_value, reordered$F_value, tolerance = 1e-9)
})

test_that("public repeated DIF jointly adjusts correlated factors and occasions", {
  d <- .joint_dif_fixture(n = 1200L)
  set.seed(9152)
  X <- matrix(rbinom(nrow(d) * 5, 1, .5), nrow(d), 5,
              dimnames = list(NULL, paste0("I", 1:5)))
  fit <- rasch(X, id = d$pid, factors = d[c("A", "B", "occasion")])
  # Fix known residuals to isolate the ANOVA adjustment from calibration.
  fit$residuals[, 3] <- d$z
  fit$person$theta <- as.numeric(d$ci)
  fit$person$extreme <- FALSE
  out <- dif_anova(fit, within = "occasion", n_groups = 3)
  A <- out$terms[out$terms$item == "I3" & out$terms$term == "A", ]
  B <- out$terms[out$terms$item == "I3" & out$terms$term == "B", ]
  expect_gt(A$p_adj, .05)
  expect_lt(B$p_adj, .001)
  expect_match(out$between_covariance, "CR3")
  expect_identical(out$algorithm, "joint-between-1")
  expect_no_error(.validate_dif_result(out, fit))
  old <- out; old$algorithm <- NULL; old$result_signature <- NULL
  old$result_signature <- .fit_boot_md5(unclass(old))
  expect_error(.validate_dif_result(old, fit), "predates joint adjustment")
})

test_that("project migration drops authenticated obsolete mixed DIF dependencies", {
  d <- .joint_dif_fixture(n = 240L, seed = 9161)
  X <- matrix(rbinom(nrow(d) * 5, 1, .5), nrow(d), 5,
              dimnames = list(NULL, paste0("I", 1:5)))
  fit <- rasch(X, id = d$pid, factors = d[c("A", "B", "occasion")])
  fit$residuals[, 3] <- d$z
  fit$person$theta <- as.numeric(d$ci)
  fit$person$extreme <- FALSE
  current <- dif_anova(fit, within = "occasion", n_groups = 3)
  old <- current
  old$algorithm <- "joint-between-legacy"
  old$result_signature <- NULL
  old$result_signature <- .fit_boot_md5(unclass(old))
  project <- .seal_app_project(list(
    format = "rasch-shiny-project", schema = 2L,
    data = d, model_type = "rasch", base_fit = fit,
    rasch_steps = list(), btl_steps = list(), kept_fits = list(),
    kept_fit_code = list(), settings = list(), resources = list(),
    simulation = list(), results = list(dif = old,
      dif_bootstrap = list(db = "obsolete"),
      resolve = list(algorithm = "factor-design-resolution-1", effects = "main"))))
  path <- tempfile(fileext = ".rasch")
  on.exit(unlink(path), add = TRUE)
  saveRDS(project, path)

  expect_warning(restored <- .read_app_project(path),
                 "mixed-panel DIF analysis.*omitted")
  expect_null(restored$results$dif)
  expect_null(restored$results$dif_bootstrap)
  expect_null(restored$results$resolve)
  expect_identical(restored$base_fit, fit)
  expect_no_error(.validate_app_project(restored))

  # A valid project seal cannot bypass the result's active-fit binding.
  bad <- project
  bad$base_fit <- rasch(X, id = d$pid)
  saveRDS(.seal_app_project(bad), path)
  expect_error(.read_app_project(path), "computed from a different fitted model")
})
