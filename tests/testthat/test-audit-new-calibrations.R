# Independent edge cases for the September calibration additions.
audit_rank_data <- function(n = 100L, judges = 20L) {
  set.seed(914)
  beta <- c(A = -1, B = -0.4, C = 0.2, D = 1.2)
  do.call(rbind, lapply(seq_len(n), function(i) {
    ord <- sample(names(beta), prob = exp(beta))
    data.frame(ranking = i, object = ord, rank = seq_along(ord),
               judge = paste0("J", (i - 1L) %% judges + 1L))
  }))
}

test_that("joint Newton convergence requires an identified maximum", {
  z <- rasch:::.cj_newton(0, function(x) x^2 - x^4,
    function(x) 2*x - 4*x^3, function(x) matrix(2 - 12*x^2), maxit = 10)
  expect_false(z$converged)
  flat <- rasch:::.cj_newton(c(0, 0), function(x) -x[1]^2,
    function(x) c(-2*x[1], 0), function(x) diag(c(-2, 0)))
  expect_false(flat$converged)
  peak <- rasch:::.cj_newton(0, function(x) -x^2,
    function(x) -2*x, function(x) matrix(-2))
  expect_true(peak$converged)
})

test_that("person-mode responses and anchor numbers cannot be truncated", {
  a <- data.frame(item = c("I1", "I2"), k = 1, tau = c(-1, 1))
  d <- data.frame(I1 = c(0, 1, 1), I2 = c(0, 0, 1))
  prepare <- function(d, anchors = a) rasch:::.cj_person_test(
    d, anchors, NULL, NULL, -1, "responses", FALSE)
  for (bad in c(0.9, Inf, NaN, "wrong")) {
    x <- d; x$I1[2] <- bad
    expect_error(prepare(x), "score|response", info = as.character(bad))
  }
  b <- a; b$k[1] <- 1.9
  expect_error(prepare(d, b), "threshold")
  b <- a; b$tau <- factor(c("-1", "1"), levels = c("1", "-1"))
  expect_equal(unlist(prepare(d, b)$tau_list, use.names = FALSE), c(-1, 1))
})

test_that("judgement numeric labels are read as labels, not factor codes", {
  key <- rasch:::.cj_object_keys
  expect_error(key("I", "nonsense", "I", 3L, "comparisons"), "threshold")
  expect_equal(key("I", factor("2", levels = c("2", "1")), "I", 3L,
                   "comparisons"), "I:2")
  rk <- data.frame(ranking = 1, item = c("A", "B", "C"),
                   rank = factor(c("1", "2", "3"), levels = c("3", "1", "2")))
  z <- rasch:::.cj_ranking_keys(rk, "ranking", "item", "rank", "threshold",
                               c("A", "B", "C"), rep(1L, 3), character())
  expect_identical(z$rk[[1]], c("A", "B", "C"))
  rk$rank <- c(1, 2, Inf)
  expect_error(rasch:::.cj_ranking_keys(rk, "ranking", "item", "rank", "threshold",
                 c("A", "B", "C"), rep(1L, 3), character()), "rank")
})

test_that("ranking inference never bypasses the judge support guard", {
  f <- pl(audit_rank_data(judges = 5), judge = "judge")
  expect_false(f$se_available)
  expect_true(all(is.na(f$invariance$objects$se)))
  expect_true(all(is.na(f$invariance$objects$p)))
  expect_true(is.na(f$invariance$p))
  expect_true(is.na(f$reversal$p))
})

test_that("ranking role columns cannot be NULL", {
  for (nm in c("ranking", "object", "rank")) {
    args <- list(data = audit_rank_data())
    args[nm] <- list(NULL)
    expect_error(do.call(pl, args), "exactly one column")
  }
})

test_that("ranking estimates are stable after a large origin shift", {
  d <- audit_rank_data()
  a <- pl(d, anchors = c(A = -1), se = "model")
  b <- pl(d, anchors = c(A = 999), se = "model")
  expect_true(b$converged)
  expect_equal(b$objects$location - 1000, a$objects$location, tolerance = 1e-6)
  expect_equal(b$objects$se, a$objects$se, tolerance = 1e-6)
  expect_equal(b$loglik, a$loglik, tolerance = 1e-7)
  expect_equal(b$reversal$loglik_reversed, a$reversal$loglik_reversed,
               tolerance = 1e-7)
})

test_that("DIF Wald retains an unavailable row when no levels meet min_n", {
  set.seed(192)
  d <- as.data.frame(matrix(rbinom(400, 1, .5), 100, 4))
  d$grp <- rep(c("a", "b"), each = 50)
  f <- rasch(d, factors = "grp")
  w <- dif_wald(f, min_n = 101)
  expect_equal(nrow(w$summary), 4L)
  expect_equal(nrow(w$levels), 0L)
  expect_true(all(is.na(w$summary$p_adj)))
  expect_true(all(w$summary$df == 0))
})

test_that("DTF does not infer from an indefinite calibration covariance", {
  set.seed(81)
  d <- as.data.frame(matrix(rbinom(600, 1, .5), 100, 6))
  d$grp <- rep(c("a", "b"), each = 50)
  f <- split_items(rasch(d, factors = "grp"), "V1", by = "grp")
  f$est$cov_tau <- -f$est$cov_tau
  z <- dtf(f, "grp", grid = 5)
  expect_true(all(is.na(z$items$se)))
  expect_true(all(is.na(z$items$significant)))
  expect_true(all(is.na(z$test$p)))
})

test_that("ranking invariance uses within-cluster contrast influences", {
  d <- audit_rank_data(160, 20)
  f <- pl(d, judge = "judge")
  k <- nrow(f$objects)
  B <- rbind(diag(k - 1L), -1)
  cells <- rasch:::.pl_cells(match(d$object, f$objects$object), d$rank,
                             d$ranking, k)
  cl <- (cells$ranking - 1L) %% 20L + 1L
  influences <- lapply(1:2, function(j) {
    stages <- if (j == 1L) cells$position == 1 else cells$position > 1
    selected <- stages[cells$stage]
    keep <- which(stages)
    z <- list(stage = match(cells$stage[selected], keep),
              obj = cells$obj[selected], y = cells$y[selected],
              n_stages = length(keep))
    beta <- if (j == 1L) f$invariance$objects$first else f$invariance$objects$later
    H <- -crossprod(B, rasch:::.pl_hess(beta, z, k) %*% B)
    G <- rasch:::.pl_scores(beta, z, k, cl[keep]) %*% B
    sqrt(20 / 19) * G %*% solve(H) %*% t(B)
  })
  expected <- sqrt(colSums((influences[[1]] - influences[[2]])^2))
  marginal_only <- sqrt(colSums(influences[[1]]^2 + influences[[2]]^2))
  expect_gt(max(abs(expected - marginal_only)), .001)
  expect_equal(f$invariance$objects$se, expected, tolerance = 1e-10)
  expect_equal(f$invariance$objects$p,
    2 * pt(-abs(f$invariance$objects$z), df = 19), tolerance = 1e-12)
})

test_that("an extreme fixed ranking anchor is never silently discarded", {
  d <- audit_rank_data()
  ext <- do.call(rbind, lapply(1:30, function(i)
    data.frame(ranking = 100L + i, object = c("Z", "B"), rank = 1:2,
               judge = "J1")))
  f <- pl(rbind(d, ext), anchors = c(A = -1, Z = 3), se = "model")
  expect_equal(f$objects$location[f$objects$object == "Z"], 3)
  expect_equal(f$objects$se[f$objects$object == "Z"], 0)
  expect_false(f$objects$extreme[f$objects$object == "Z"])
})

test_that("linked covariance includes direct and link-estimation cross terms", {
  I1 <- diag(c(4, 5)); I2 <- diag(c(3, 6))
  J1 <- cbind(diag(2), 0)
  J2 <- cbind(diag(2), c(1, -1))
  info <- list(a = list(I = I1, J = J1), b = list(I = I2, J = J2))
  V <- solve(t(J1) %*% I1 %*% J1 + t(J2) %*% I2 %*% J2)
  direct <- list(a = solve(I1), b = -solve(I2))
  link <- cbind(matrix(0, 2, 2), c(1, -1))
  actual <- rasch:::.cj_linked_cov(direct, link, V, info)
  A <- cbind(direct$a, direct$b) + link %*% V %*% cbind(t(J1), t(J2))
  expected <- A %*% diag(c(diag(I1), diag(I2))) %*% t(A)
  expect_equal(actual, expected, tolerance = 1e-12)
  naive <- solve(I1) + solve(I2)
  expect_gt(max(abs(actual - naive)), .1)
  expect_true(rasch:::.covariance_supports_wald(actual))
})

test_that("shared response items retain the same category mapping across tests", {
  set.seed(902)
  x <- matrix(rbinom(400, 1, .5), 100, 4,
               dimnames = list(NULL, c("A", "B", "C", "D")))
  y <- x; y[, "A"] <- y[, "A"] + 1L
  cj <- data.frame(object_a = c("A", "B"), object_b = c("B", "A"),
                    winner = c("A", "B"))
  expect_error(rasch_cj(list(first = x, second = y), cj), "category mappings")
})

test_that("joint object names with a colon do not become threshold selectors", {
  set.seed(92)
  x <- matrix(rbinom(480, 1, .5), 120, 4,
               dimnames = list(NULL, c("A:1", "B", "C", "D")))
  cj <- data.frame(object_a = rep("A:1", 20), object_b = rep("B", 20),
                    winner = rep(c("A:1", "B"), 10))
  f <- rasch_cj(x, cj, units = c(comparisons = 1))
  expect_true(f$converged)
  expect_true("A:1" %in% f$objects$item)
  expect_true(all(is.na(f$objects$threshold)))
  expect_error(rasch:::.cj_object_keys("A", 1, c("A", "A:1"), c(1, 1),
                                     "comparisons"), "collide")
})

test_that("DTF bundle homogeneity counts independent contrasts", {
  set.seed(925)
  d <- as.data.frame(matrix(rbinom(2400, 1, .5), 400, 6))
  d$grp <- rep(c("a", "b"), each = 200)
  f <- split_items(rasch(d, factors = "grp"), "V1", by = "grp")
  z <- dtf(f, "grp", bundles = list(one_free = c("V1", "V2", "V3")), grid = 5)
  expect_equal(z$bundles$df_hom, 1)
  expect_true(is.finite(z$bundles$p_hom))
  p_item <- z$items$p[z$items$item == "V1"]
  expect_equal(z$bundles$p_hom, p_item, tolerance = 1e-10)
  # A retained repeated-person covariance uses the supported cluster df.
  f$est$cluster_support <- list(repeated = TRUE, n = 12L)
  f$est$cluster_inference <- TRUE
  rep <- dtf(f, "grp", grid = 5)
  expect_equal(rep$ref_df, 11)
  expect_equal(rep$test$p, 2 * pt(-abs(rep$test$z), 11))
})

test_that("DTF plot defaults can be overridden", {
  x <- structure(list(groups = "b", reference = "a", max_score = 3,
    curves = data.frame(group = "b", theta = 1:3, expected_reference = 1:3,
      expected_group = 1:3, shift_score = 0, se_score = .1)), class = "rasch_dtf")
  pdf(NULL)
  on.exit(dev.off())
  expect_silent(plot_dtf(x, main = "Custom title", ylim = c(0, 4)))
})

test_that("separate response booklets have separate origin constraints", {
  set.seed(924)
  X <- matrix(rbinom(1600, 1, .5), 400, 4,
              dimnames = list(NULL, c("A", "B", "C", "D")))
  X[1:200, 3:4] <- NA; X[201:400, 1:2] <- NA
  pairs <- t(combn(colnames(X), 2))
  cj <- data.frame(object_a = rep(pairs[, 1], each = 20),
                    object_b = rep(pairs[, 2], each = 20))
  cj$winner <- ifelse(rep(c(TRUE, FALSE), 60), cj$object_a, cj$object_b)
  f <- rasch_cj(X, cj, units = c(comparisons = 1))
  expect_true(f$converged)
  # Separate response booklets: 1 + 1 free parameters; comparisons: 3.
  # Joint model: 3 free locations. Thus the LR has 2, not 3, df.
  expect_equal(f$invariance$lr$df, 2)
  expect_true(is.finite(f$invariance$lr$p))
})

test_that("saved DIF cannot reuse pre-common-interval split results", {
  set.seed(930)
  d <- as.data.frame(matrix(rbinom(1200, 1, .5), 200, 6))
  d$grp <- rep(c("a", "b"), each = 100)
  f <- split_items(rasch(d, factors = "grp"), "V1", by = "grp")
  da <- dif_anova(f)
  expect_silent(rasch:::.validate_dif_result(da, f))
  old <- da; old$interval_algorithm <- NULL
  old$result_signature <- NULL
  old$result_signature <- rasch:::.fit_boot_md5(unclass(old))
  expect_error(rasch:::.validate_dif_result(old, f), "common split-item")
  project <- rasch:::.seal_app_project(list(
    format = "rasch-shiny-project", schema = 2L, package_version = "test",
    created = "2026-09-24", data = data.frame(f$X, grp = d$grp, check.names = FALSE),
    model_type = "rasch", base_fit = f, rasch_steps = list(), btl_steps = list(),
    rcode = "fit <- rasch(data)", kept_fits = list(),
    settings = list(model_type = "rasch", item_cols = colnames(f$X),
                    dif_effects = "main", dif_alpha = .05),
    resources = list(), simulation = list(), results = list(dif = old)))
  path <- tempfile(fileext = ".rasch")
  on.exit(unlink(path))
  saveRDS(project, path)
  expect_warning(restored <- rasch:::.read_app_project(path), "split-item class intervals")
  expect_null(restored$results[["dif"]])
  expect_equal(restored$base_fit$items, f$items)
})
