simP <- function(theta, tau) { x <- 0:length(tau); p <- exp(x * theta - c(0, cumsum(tau))); p / sum(p) }

test_that("dimensionality test separates 1D from 2D data", {
  set.seed(1); Np <- 1500; L <- 20
  d <- scale(seq(-2, 2, length.out = L), scale = FALSE)[, 1]; th <- rnorm(Np, 0, 1.4)
  X1 <- matrix(rbinom(Np * L, 1, plogis(outer(th, d, "-"))), Np, L); colnames(X1) <- sprintf("U%02d", 1:L)
  dt1 <- dimensionality_test(
    rasch(X1, model = "PCM"), items_positive = colnames(X1)[1:10],
    items_negative = colnames(X1)[11:20], min_score_points = 2)
  expect_false(dt1$multidimensional)

  set.seed(2)
  thA <- rnorm(Np, 0, 1.4); thB <- 0.3 * thA + sqrt(1 - 0.3^2) * rnorm(Np, 0, 1.4)
  XA <- matrix(rbinom(Np * 10, 1, plogis(outer(thA, d[1:10], "-"))), Np, 10)
  XB <- matrix(rbinom(Np * 10, 1, plogis(outer(thB, d[11:20], "-"))), Np, 10)
  X2 <- cbind(XA, XB); colnames(X2) <- sprintf("D%02d", 1:20)
  dt2 <- dimensionality_test(
    rasch(X2, model = "PCM"), items_positive = colnames(X2)[1:10],
    items_negative = colnames(X2)[11:20], min_score_points = 2)
  expect_true(dt2$multidimensional)
  expect_gt(dt2$prop_significant, dt1$prop_significant)
})

test_that("residual diagnostics refuse an empty pairwise reference", {
  set.seed(101)
  X <- matrix(rbinom(300 * 4, 1, 0.5), 300, 4,
              dimnames = list(NULL, paste0("I", 1:4)))
  fit <- rasch(X)
  # Preserve a valid fitted object while representing a structurally empty
  # residual-overlap graph. Returning an average of NaN here is not a result.
  fit$residuals[,] <- NA_real_
  expect_error(residual_correlations(fit), "no item pair",
               class = "rasch_refusal")
})

test_that("a short automatic dimensionality split is descriptive with a caution", {
  set.seed(11)
  X <- matrix(rbinom(500 * 12, 1, 0.5), 500, 12,
              dimnames = list(NULL, paste0("S", 1:12)))
  dt <- dimensionality_test(rasch(X))
  # The short-subtest caution and the data-driven-split withholding are
  # separate: the descriptive binomial reading remains available.
  expect_true(is.na(dt$multidimensional))
  expect_true(dt$binomial_multidimensional %in% c(TRUE, FALSE))
  expect_match(dt$caution, "score points")
  expect_true(all(dt$score_points < 15))
  printed <- capture.output(print(dt))
  expect_true(any(grepl("withheld for the data-driven split", printed,
                        fixed = TRUE)))
  expect_false(any(grepl("fit_signature", printed, fixed = TRUE)))
})

test_that("the dimensionality bootstrap calibrates the split it is given", {
  set.seed(21); N <- 300; L <- 12
  d <- seq(-1.5, 1.5, length.out = L)
  X <- matrix(rbinom(N * L, 1, plogis(outer(rnorm(N), d, "-"))), N, L,
              dimnames = list(NULL, sprintf("I%02d", 1:L)))
  fit <- rasch(X)
  dt <- dimensionality_test(fit, min_score_points = 2)
  expect_s3_class(dt, "rasch_dimensionality_test")
  expect_no_error(.validate_dimensionality_test(dt, fit))
  expect_null(dt$p_boot)
  expect_null(dt$bootstrap)
  b1 <- dimensionality_test(fit, min_score_points = 2, B = 19, seed = 3)
  # the observed reading is what it was; the bootstrap adds its own
  expect_equal(b1$prop_significant, dt$prop_significant)
  expect_equal(b1$ci, dt$ci)
  expect_equal(b1$items_positive, dt$items_positive)
  expect_true(b1$p_boot > 0 && b1$p_boot <= 1)
  expect_identical(b1$multidimensional, b1$p_boot <= 0.05)
  expect_equal(b1$bootstrap$B, 19L)
  expect_length(b1$bootstrap$null, b1$bootstrap$B_used)
  expect_true(all(b1$bootstrap$null >= 0 & b1$bootstrap$null <= 1))
  expect_equal(b1$prop_null, mean(b1$bootstrap$null))
  expect_equal(b1$bootstrap$minimum_usable, 10L)
  b_alpha <- dimensionality_test(fit, min_score_points = 2, B = 19,
                                 seed = 3, alpha = .10)
  expect_identical(b_alpha$multidimensional, b_alpha$p_boot <= .10)
  # a seed reproduces the replicates
  b2 <- dimensionality_test(fit, min_score_points = 2, B = 19, seed = 3)
  expect_equal(b1$bootstrap$null, b2$bootstrap$null)
  # a split fixed in advance is kept in every replicate
  m <- dimensionality_test(fit, items_positive = colnames(X)[1:6],
                           items_negative = colnames(X)[7:12],
                           min_score_points = 2, B = 19, seed = 3)
  expect_identical(m$split, "manual")
  expect_true(m$p_boot > 0 && m$p_boot <= 1)
  expect_length(m$bootstrap$null, m$bootstrap$B_used)

  expect_error(dimensionality_test(fit, B = -1), "`B`")
  expect_error(dimensionality_test(fit, B = 2.5), "`B`")
  expect_error(dimensionality_test(fit, B = 5, workers = 0), "`workers`")
  expect_error(dimensionality_test(fit, B = 5, seed = -1), "`seed`")
  expect_error(dimensionality_test(fit, seed = -1), "`seed`")
  for (cl in c("rasch_efrm", "rasch_mfrm", "rasch_explanatory")) {
    f2 <- fit; class(f2) <- c(cl, "rasch")
    expect_error(dimensionality_test(f2, B = 5), "generating structure",
                 class = "rasch_refusal")
  }
  f3 <- fit; f3$disc <- c(rep(1, L - 1), 2)
  expect_error(dimensionality_test(f3, B = 5), "frame units",
               class = "rasch_refusal")
  f4 <- fit; f4$refit_spec$pc_components <- 2L
  expect_error(dimensionality_test(f4, B = 5), "principal components",
               class = "rasch_refusal")
})

test_that("person-subset and simulated residual inference refuse repeated IDs", {
  d <- simulate_rasch(120, 8, seed = 212)
  items <- sprintf("I%02d", 1:8)
  repeated <- d[rep(seq_len(nrow(d)), each = 2L), ]
  fit <- rasch(repeated, id = "id", items = items)

  expect_error(
    dimensionality_test(
      fit, items_positive = items[1:4], items_negative = items[5:8]),
    "requires one response row per person", class = "rasch_refusal")
  expect_no_error(.scree_analysis(fit, parallel = FALSE))
  expect_error(plot_scree(fit, parallel = TRUE, reps = 20),
               "person identifier occurs on several response rows",
               class = "rasch_refusal")

  # Missing identifiers denote unknown people, not a repeated person. They
  # must not trigger the repeated-row refusal merely because NA recurs.
  partly_unknown <- d
  partly_unknown$id <- as.character(partly_unknown$id)
  partly_unknown$id[1:4] <- c(NA, NA, "", "  ")
  fit_unknown <- rasch(partly_unknown, id = "id", items = items)
  expect_no_error(dimensionality_test(
    fit_unknown, items_positive = items[1:4], items_negative = items[5:8]))
  expect_no_error(plot_scree(fit_unknown, parallel = TRUE, reps = 20,
                             seed = 212))
})

test_that("unavailable dimensionality comparisons keep a stable result shape", {
  set.seed(211)
  X <- matrix(rbinom(120 * 6, 1, .5), 120, 6,
              dimnames = list(NULL, paste0("I", 1:6)))
  fit <- rasch(X)
  thin <- fit
  thin$X <- thin$X[seq_len(9), , drop = FALSE]
  out <- dimensionality_test(
    thin, items_positive = colnames(X)[1:3],
    items_negative = colnames(X)[4:6], min_score_points = 2)
  expect_s3_class(out, "rasch_dimensionality_test")
  expect_true(is.na(out$multidimensional))
  expect_match(out$note, "fewer than 10 usable persons")
  expect_no_error(.validate_dimensionality_test(out, thin))
})

test_that("a residual-component split runs above a fixed split under the model", {
  skip_on_cran()
  # well-populated four-category items, so the replicates lose no category
  set.seed(41); N <- 300; L <- 12
  th <- rnorm(N, 0, 1.5); loc <- seq(-1, 1, length.out = L)
  X <- sapply(seq_len(L), function(i) {
    tau <- loc[i] + c(-1, 0, 1)
    vapply(th, function(t) sample(0:3, 1, prob = simP(t, tau)), 0L)
  })
  colnames(X) <- sprintf("P%02d", seq_len(L))
  fit <- rasch(X, model = "PCM")
  auto <- dimensionality_test(fit, B = 19, seed = 3)
  fixed <- dimensionality_test(fit, items_positive = colnames(X)[1:6],
                               items_negative = colnames(X)[7:12],
                               B = 19, seed = 3)
  expect_equal(auto$bootstrap$B_used, 19L)
  expect_equal(fixed$bootstrap$B_used, 19L)
  # the same seed generates the same replicates, so the two nulls differ
  # only in how each replicate is split: the split chosen from the
  # replicate's own residuals is chosen to disagree, and its proportion
  # runs higher than the fixed split's on the same data
  expect_gt(auto$prop_null, fixed$prop_null)
  expect_true(auto$p_boot > 0 && auto$p_boot <= 1)
})

test_that("the dimensionality bootstrap withholds a depleted null", {
  set.seed(22); N <- 200; L <- 10
  X <- matrix(rbinom(N * L, 1, plogis(outer(rnorm(N), seq(-1, 1, length.out = L), "-"))),
              N, L, dimnames = list(NULL, sprintf("I%02d", 1:L)))
  fit <- rasch(X)
  fake <- function(n_ok, n_nonconv, n_err) function(n, fun, ...) c(
    as.list(rep(0.05, n_ok)),
    rep(list(.fit_boot_failure("nonconverged")), n_nonconv),
    rep(list(.fit_boot_failure("error")), n_err))
  refused <- tryCatch(with_mocked_bindings(
    dimensionality_test(fit, min_score_points = 2, B = 40, seed = 1),
    .rasch_boot_apply = fake(29L, 6L, 5L), .package = "rasch"),
    error = identity)
  expect_s3_class(refused, "rasch_fit_bootstrap_refusal")
  expect_match(conditionMessage(refused), "at least 36")
  expect_equal(refused$B_used, 29L)
  expect_equal(refused$B_nonconverged, 6L)
  expect_equal(refused$B_errors, 5L)
  # enough to form the null, but thinned: said out loud
  seen <- character()
  thinned <- withCallingHandlers(
    with_mocked_bindings(
      dimensionality_test(fit, min_score_points = 2, B = 20, seed = 1),
      .rasch_boot_apply = fake(12L, 5L, 3L), .package = "rasch"),
    warning = function(w) {
      seen <<- c(seen, conditionMessage(w))
      invokeRestart("muffleWarning")
    })
  expect_true(any(grepl("not lost at random", seen, fixed = TRUE)))
  expect_true(any(grepl("smallest attainable", seen, fixed = TRUE)))
  expect_equal(thinned$bootstrap$B_used, 12L)
  expect_equal(thinned$bootstrap$B_nonconverged, 5L)
  expect_equal(thinned$bootstrap$B_errors, 3L)
  expect_equal(thinned$prop_null, 0.05)
  expect_equal(thinned$p_boot,
               .boot_p(thinned$prop_significant, rep(0.05, 12), "upper"))
  expect_equal(thinned$bootstrap_resolution, 1 / 13)
  expect_true(is.na(thinned$multidimensional))
  expect_match(thinned$verdict_method, "resolution insufficient")
  expect_match(paste(capture.output(print(thinned)), collapse = "\n"),
               "bootstrap resolution is insufficient")
})

test_that("Q3 binary flags require an explicit heuristic", {
  set.seed(12)
  X <- matrix(rbinom(500 * 8, 1, 0.5), 500, 8,
              dimnames = list(NULL, paste0("Q", 1:8)))
  rc <- residual_correlations(rasch(X))
  expect_true(all(is.na(rc$pairs$flagged)))
  expect_equal(nrow(rc$flagged), 0L)
  expect_match(rc$note, "withheld")
})

test_that("local dependence is flagged for a near-duplicated item", {
  set.seed(1); Np <- 1500; L <- 20
  d <- scale(seq(-2, 2, length.out = L), scale = FALSE)[, 1]; th <- rnorm(Np, 0, 1.4)
  X <- matrix(rbinom(Np * L, 1, plogis(outer(th, d, "-"))), Np, L); colnames(X) <- sprintf("U%02d", 1:L)
  set.seed(3); X[, 5] <- ifelse(runif(Np) < 0.85, X[, 4], X[, 5])
  fl <- residual_correlations(rasch(X, model = "PCM"), flag = 0.2)$flagged
  expect_true(any((fl$item_a == "U04" & fl$item_b == "U05") |
                  (fl$item_a == "U05" & fl$item_b == "U04")))
})

test_that("uniform DIF is detected on planted items only, across two factors", {
  set.seed(4); Np <- 2000; L <- 15
  d <- scale(seq(-2, 2, length.out = L), scale = FALSE)[, 1]
  grp <- rep(c("ref", "foc"), each = Np / 2)
  sex <- sample(c("f", "m"), Np, replace = TRUE)
  th <- rnorm(Np, 0, 1.4)
  shift <- matrix(0, Np, L); shift[grp == "foc", 3] <- 1.0; shift[grp == "foc", 10] <- -1.0
  X <- matrix(rbinom(Np * L, 1, plogis(outer(th, d, "-") - shift)), Np, L)
  colnames(X) <- sprintf("G%02d", 1:L)

  fit <- rasch(data.frame(X, group = grp, sex = sex), factors = c("group", "sex"),
               n_groups = 6)
  # the joint main-effects model gives one row per item and factor term
  s <- dif_anova(fit)$summary
  expect_setequal(unique(s$term), c("group", "sex"))
  dg <- s[s$term == "group", ]
  expect_true(dg$uniform_DIF[dg$item == "G03"])
  expect_true(dg$uniform_DIF[dg$item == "G10"])
  expect_equal(sum(dg$uniform_DIF[!dg$item %in% c("G03", "G10")]), 0)
  expect_equal(sum(s$uniform_DIF[s$term == "sex"]), 0)
})

test_that("threshold disordering is detected", {
  set.seed(5); Np <- 2000; th <- rnorm(Np, 0, 1.4)
  tau_dis <- c(1.5, -1.5, 0.8); tau_ok <- c(-1.0, 0.0, 1.0)
  X <- cbind(sapply(th, function(t) sample(0:3, 1, prob = simP(t, tau_dis))),
             sapply(th, function(t) sample(0:3, 1, prob = simP(t, tau_ok))),
             sapply(th, function(t) sample(0:3, 1, prob = simP(t, tau_ok + 0.4))))
  colnames(X) <- c("DIS", "ok1", "ok2")
  td <- rasch(X, model = "PCM", n_groups = 6)$thresholds_diag
  expect_false(td[["DIS"]]$ordered)
  expect_true(td[["ok1"]]$ordered)
})

test_that("ID and factors carry through to the person table", {
  set.seed(6); Np <- 400; L <- 8
  d <- scale(seq(-1.5, 1.5, length.out = L), scale = FALSE)[, 1]
  X <- matrix(rbinom(Np * L, 1, plogis(outer(rnorm(Np), d, "-"))), Np, L)
  colnames(X) <- paste("My item", 1:L)   # names with spaces survive
  df <- data.frame(sid = sprintf("S%03d", 1:Np), X, grp = rep(c("a", "b"), Np / 2),
                   check.names = FALSE)
  fit <- rasch(df, id = "sid", factors = "grp")
  expect_identical(fit$person$id, df$sid)
  expect_identical(as.character(fit$person$grp), df$grp)
  expect_identical(fit$items$item, paste("My item", 1:L))
  expect_identical(colnames(fit$residuals), paste("My item", 1:L))
})

test_that("data preparation collapses gaps and drops constants with notes", {
  set.seed(8); Np <- 500
  th <- rnorm(Np)
  x1 <- rbinom(Np, 1, plogis(th)) * 2L              # categories 0, 2 -> collapse
  x2 <- rbinom(Np, 1, plogis(th - 0.5))
  x3 <- rbinom(Np, 1, plogis(th + 0.5))
  x4 <- rep(1L, Np)                                  # constant -> dropped
  fit <- rasch(cbind(a = x1, b = x2, c = x3, d = x4))
  expect_equal(ncol(fit$X), 3)
  expect_equal(max(fit$X[, "a"]), 1)
  expect_true(any(grepl("rescored", fit$notes)))
  expect_true(any(grepl("constant", fit$notes)))
})

test_that("a category seen only in extreme response patterns is merged", {
  # a PCM category that only all-minimum (or all-maximum) persons use enters
  # no pairwise conditional contribution: its threshold has no finite
  # estimate, and before the merge the solver ran it to -34 and then
  # reported a singular information matrix
  sim <- simulate_rasch(n_persons = 400, n_items = 20, model = "PCM",
                        n_categories = 4, seed = 5)
  X <- as.matrix(sim[, grep("^I", names(sim))])
  who <- which(X[, "I01"] == 0)
  expect_length(who, 1L)
  expect_equal(sum(X[who, ]), 0)                    # an all-zero person
  fit <- rasch(X)
  expect_true(fit$est$converged)
  expect_equal(unname(fit$m[1]), 2L)
  expect_true(any(grepl("item I01 rescored: category 0 observed only in extreme",
                        fit$notes, fixed = TRUE)))
  # identical to collapsing the category by hand
  XA <- X; XA[, "I01"] <- pmax(XA[, "I01"] - 1L, 0L)
  expect_equal(fit$thresholds$tau, rasch(XA)$thresholds$tau)
  # the mirror image at the top of the scale
  ft <- rasch(3L - X)
  expect_equal(unname(ft$m[1]), 2L)
  expect_true(any(grepl("category 3 observed only in extreme", ft$notes)))
  # the RSM shares its thresholds across items, so nothing is merged
  fr <- rasch(X, model = "RSM")
  expect_true(fr$est$converged)
  expect_false(any(grepl("extreme response patterns", fr$notes)))
  # pcml() names the category instead of reporting a singular matrix
  expect_error(pcml(X), "I01 \\(category 0\\)")
  # an item with no conditional information at all is dropped
  Xd <- X
  Xd[, "I01"] <- 0L
  Xd[who, ] <- 3L; Xd[who, "I01"] <- 1L            # its one 1 is all-maximum
  fd <- rasch(Xd)
  expect_false("I01" %in% fd$items$item)
  expect_true(any(grepl("dropped item(s) with no conditional information",
                        fd$notes, fixed = TRUE)))
  # an anchored item cannot be merged silently
  expect_error(rasch(X, anchors = data.frame(item = "I01", k = 1, tau = -2)),
               "rescored during data preparation")
  # an item whose every threshold is fixed estimates nothing, so it keeps
  # its coding: the merge exempts it exactly as pcml() exempts it from the
  # uninformative-category check. A numeric index resolves to the same item
  full <- data.frame(item = "I01", k = 1:3, tau = c(-2, 0, 2))
  ff <- rasch(X, anchors = full)
  expect_equal(unname(ff$m[1]), 3L)
  expect_false(any(grepl("I01 rescored", ff$notes)))
  expect_true(all(ff$thresholds$anchored[ff$thresholds$item == 1]))
  full$item <- 1
  expect_equal(rasch(X, anchors = full)$thresholds$tau, ff$thresholds$tau)
  # a location anchor on a polytomous item leaves its thresholds free, so
  # the category is merged and the anchor refused as before
  expect_error(rasch(X, anchors = data.frame(item = "I01", k = NA, tau = 0)),
               "rescored during data preparation")
})

test_that("reliability and fit summaries are coherent", {
  set.seed(9); Np <- 1000; L <- 15
  d <- scale(seq(-2, 2, length.out = L), scale = FALSE)[, 1]
  X <- matrix(rbinom(Np * L, 1, plogis(outer(rnorm(Np, 0, 1.5), d, "-"))), Np, L)
  colnames(X) <- sprintf("I%02d", 1:L)
  fit <- rasch(X)
  expect_gt(fit$psi$PSI, 0.5); expect_lt(fit$psi$PSI, 1)
  expect_gt(fit$alpha$alpha, 0.5); expect_lt(fit$alpha$alpha, 1)
  expect_gt(fit$total_chisq_p, 1e-6)   # well-fitting data should not collapse
  expect_lt(abs(fit$item_fit_summary$mean), 1)
  expect_false(any(fit$items$p_adj < 0.05, na.rm = TRUE))
  expect_true(fit$separation_quality %in%
                c("reasonable", "good", "excellent"))
  expect_identical(fit$power_of_fit, fit$separation_quality)
})

test_that("save_outputs writes the full set of tables and plots", {
  set.seed(10); Np <- 300; L <- 6
  d <- scale(seq(-1.5, 1.5, length.out = L), scale = FALSE)[, 1]
  X <- matrix(rbinom(Np * L, 1, plogis(outer(rnorm(Np), d, "-"))), Np, L)
  colnames(X) <- sprintf("I%02d", 1:L)
  fit <- rasch(data.frame(X, g = rep(c("x", "y"), Np / 2)), factors = "g")
  bs <- suppressWarnings(fit_bootstrap(fit, B = 19, seed = 1))
  out <- file.path(tempdir(), paste0("rr-test-", as.integer(runif(1, 1, 1e6))))
  scree_result <- plot_scree(fit, reps = 20, seed = 1)
  files <- save_outputs(fit, out, formats = "png", item_plots = TRUE,
                        bootstrap = bs, dimensionality = scree_result)
  expect_true(file.exists(file.path(out, "tables", "item_statistics.csv")))
  expect_true(file.exists(file.path(out, "tables", "person_estimates.csv")))
  expect_true(file.exists(file.path(out, "tables", "dif_anova.csv")))
  expect_true(file.exists(file.path(out, "tables", "bootstrap_fit.csv")))
  expect_true(file.exists(file.path(out, "tables", "bootstrap_person_fit.csv")))
  scree <- utils::read.csv(file.path(out, "tables", "residual_eigenvalues.csv"))
  expect_true(all(c("reference_critical", "parallel_p_adj",
                    "n_reference_requested") %in% names(scree)))
  expect_true(file.exists(file.path(out, "summary.txt")))
  expect_true(file.exists(file.path(out, "plots", "test_information.png")))
  expect_true(file.exists(file.path(out, "plots", "items", "I01_icc.png")))
  expect_gte(length(files), 8 + 8 + 4 * L)
  unlink(out, recursive = TRUE)
})

test_that("save_outputs writes comparative judgement results", {
  set.seed(101)
  dat <- simulate_btl(n_objects = 6, n_judges = 20,
                      reps_per_pair = 2, seed = 101)
  fit <- btl(dat, object_a = "object_a", object_b = "object_b",
             winner = "winner", judge = "judge")
  bs <- suppressWarnings(fit_bootstrap(fit, B = 19, seed = 1))
  out <- tempfile("rasch-btl-output-")
  on.exit(unlink(out, recursive = TRUE), add = TRUE)

  bd <- btl_dimensionality(fit, reps = 20, seed = 1)
  files <- save_outputs(fit, out, formats = "png", item_plots = TRUE,
                        bootstrap = bs, dimensionality = bd)
  expect_true(file.exists(file.path(out, "tables", "object_estimates.csv")))
  expect_true(file.exists(file.path(out, "tables", "pair_fit.csv")))
  expect_true(file.exists(file.path(out, "tables", "bootstrap_object_fit.csv")))
  expect_true(file.exists(file.path(out, "tables", "bootstrap_pair_fit.csv")))
  expect_true(file.exists(file.path(out, "tables", "bootstrap_judge_fit.csv")))
  expect_true(file.exists(file.path(out, "tables", "residual_bimensions.csv")))
  expect_true(file.exists(file.path(out, "plots", "residual_bimensions.png")))
  expect_true(file.exists(file.path(out, "plots", "object_locations.png")))
  expect_true(file.exists(file.path(out, "summary.txt")))
  expect_gt(length(files), 6)
})

test_that("the scree reference is model-simulated and calibrated", {
  # null Rasch data: the observed first eigenvalue must sit AT the reference,
  # not systematically above it (an independent-normal reference would put
  # every null data set 'beyond chance'); planted 2D structure must clear it
  set.seed(1); Np <- 600; L <- 10
  d <- scale(seq(-2, 2, length.out = L), scale = FALSE)[, 1]
  X <- matrix(rbinom(Np * L, 1, plogis(outer(rnorm(Np, 0, 1.3), d, "-"))), Np, L)
  colnames(X) <- sprintf("I%02d", 1:L)
  f1 <- rasch(X)
  pc1 <- residual_pca(f1, 10)
  n_conditional <- 0L
  conditional <- .fit_gen_conditional
  set.seed(101)
  ref1 <- testthat::with_mocked_bindings(
    rasch:::.scree_reference(f1, nrow(pc1$eigen_table), 20),
    .fit_gen_conditional = function(X, tau_list, na_mask) {
      n_conditional <<- n_conditional + 1L
      out <- conditional(X, tau_list, na_mask)
      expect_equal(rowSums(out, na.rm = TRUE), rowSums(X, na.rm = TRUE))
      expect_identical(is.na(out), is.na(X))
      out
    },
    .package = "rasch")
  expect_identical(n_conditional, 20L)
  expect_lt(abs(pc1$eigen_table$eigenvalue[1] /
                  attr(ref1, "mean")[1] - 1), 0.10)
  expect_identical(dim(attr(ref1, "draws")), c(20L, 10L))
  expect_identical(attr(ref1, "n_used"), 20L)
  expect_identical(attr(ref1, "n_requested"), 20L)
  expect_identical(attr(ref1, "n_nonconverged"), 0L)
  expect_identical(attr(ref1, "n_errors"), 0L)
  expect_identical(attr(ref1, "alpha"), 0.05)
  expect_equal(as.numeric(ref1), as.numeric(attr(ref1, "mean")))
  inf1 <- rasch:::.sim_upper_family(
    pc1$eigen_table$eigenvalue, attr(ref1, "draws"), 0.05)
  expect_true(all(inf1$p_adjusted >= inf1$p))
  expect_true(all(inf1$p > 0 & inf1$p <= 1))
  expect_true(all(inf1$p_adjusted > 0 & inf1$p_adjusted <= 1))
  expect_identical(inf1$significant, inf1$p_adjusted <= 0.05)
  expect_identical(inf1$significant,
                   pc1$eigen_table$eigenvalue > inf1$critical)
  expect_true(all(inf1$critical >= inf1$mean))

  set.seed(2); thA <- rnorm(Np, 0, 1.4)
  thB <- 0.3 * thA + sqrt(1 - 0.09) * rnorm(Np, 0, 1.4)
  d2 <- scale(seq(-2, 2, length.out = 16), scale = FALSE)[, 1]
  XA <- matrix(rbinom(Np * 8, 1, plogis(outer(thA, d2[1:8], "-"))), Np, 8)
  XB <- matrix(rbinom(Np * 8, 1, plogis(outer(thB, d2[9:16], "-"))), Np, 8)
  X2 <- cbind(XA, XB); colnames(X2) <- sprintf("D%02d", 1:16)
  f2 <- rasch(X2)
  pc2 <- residual_pca(f2, 10)
  set.seed(102)
  ref2 <- rasch:::.scree_reference(f2, nrow(pc2$eigen_table), 20)
  inf2 <- rasch:::.sim_upper_family(
    pc2$eigen_table$eigenvalue, attr(ref2, "draws"), 0.05)
  expect_gt(pc2$eigen_table$eigenvalue[1], inf2$critical[1])
  pdf(NULL); on.exit(dev.off())
  expect_no_error(et <- plot_scree(f2, reps = 20, seed = 72))
  expect_true(all(c("reference_mean", "reference_critical", "parallel_p",
                    "parallel_p_adj", "parallel_significant",
                    "n_reference", "n_reference_requested",
                    "n_reference_nonconverged", "n_reference_errors",
                    "reference_seed") %in% names(et)))
  expect_s3_class(et, "rasch_scree")
  expect_identical(et$reference_seed, rep(72L, nrow(et)))
  state <- .Random.seed
  invisible(plot_scree(f2, reps = 20, seed = 72))
  expect_identical(.Random.seed, state)
  expect_equal(plot_scree(f2, result = et), et)
  expect_true(all(et$parallel_p > 0 & et$parallel_p <= 1))
  expect_true(all(et$parallel_p_adj >= et$parallel_p))
  expect_identical(et$parallel_significant, et$parallel_p_adj <= 0.05)
  expect_identical(attr(et, "parallel_adjustment"),
                   "single-step leave-one-out maximum-standardised-statistic")
  expect_error(plot_scree(f2, reps = 19), "between 20")
  expect_no_error(plot_scree(f2, parallel = FALSE, reps = 1))
})

test_that("simulated upper-tail decisions use finite familywise probabilities", {
  draws <- cbind(a = seq_len(20), b = rev(seq_len(20)))
  out <- rasch:::.sim_upper_family(c(100, 20), draws, alpha = 0.05)
  expect_equal(unname(out$p), c(1 / 21, 2 / 21))
  expect_true(all(out$p_adjusted >= out$p))
  expect_identical(out$significant, out$p_adjusted <= 0.05)
  expect_identical(out$significant,
                   unname(c(100, 20) > out$critical))
  expect_true(out$significant[1])
  expect_false(out$significant[2])
  manual_max <- vapply(seq_len(nrow(draws)), function(i) {
    training <- draws[-i, , drop = FALSE]
    max((draws[i, ] - colMeans(training)) /
          apply(training, 2, stats::sd))
  }, numeric(1))
  expect_equal(out$max_null, manual_max)
})

test_that("the scree null retains the observed category structure", {
  set.seed(103)
  X <- matrix(
    sample(0:2, 300 * 5, replace = TRUE), 300, 5,
    dimnames = list(NULL, paste0("P", 1:5))
  )
  fit <- rasch(X, model = "PCM")
  shortened <- fit
  shortened$m[1L] <- shortened$m[1L] - 1L
  refused <- tryCatch(testthat::with_mocked_bindings(
    rasch:::.scree_reference(fit, k = 2L, reps = 20L, seed = 1L),
    rasch = function(...) shortened,
    .package = "rasch"
  ), error = identity)
  expect_s3_class(refused, "rasch_fit_bootstrap_refusal")
  expect_match(conditionMessage(refused),
               "only 0 of 20.*retained the fitted response structure")
  expect_identical(refused$B, 20L)
  expect_identical(refused$B_used, 0L)
  expect_identical(refused$B_nonconverged, 0L)
  expect_identical(refused$B_errors, 20L)
})

test_that("the scree null retains enough draws for its finite reference", {
  set.seed(104)
  X <- matrix(
    sample(0:1, 300 * 5, replace = TRUE), 300, 5,
    dimnames = list(NULL, paste0("D", 1:5))
  )
  fit <- rasch(X)
  shortened <- fit
  shortened$m[1L] <- 0L
  calls <- 0L
  refused <- tryCatch(testthat::with_mocked_bindings(
    rasch:::.scree_reference(fit, k = 2L, reps = 20L, seed = 2L),
    rasch = function(...) {
      calls <<- calls + 1L
      if (calls == 1L) shortened else fit
    },
    .package = "rasch"
  ), error = identity)
  expect_s3_class(refused, "rasch_fit_bootstrap_refusal")
  expect_match(conditionMessage(refused), "only 19 of 20")
  expect_identical(refused$B, 20L)
  expect_identical(refused$B_used, 19L)
  expect_identical(refused$B_errors, 1L)
})

test_that("the subset mean difference is reported without a test", {
  # An easy/hard split of unidimensional data: the two subset estimates
  # carry different estimation bias, so a t-test of their mean difference
  # rejects however unidimensional the data are. The difference is
  # descriptive; the inference belongs to the person-level comparisons.
  set.seed(1); N <- 800; L <- 20
  d <- seq(-2.5, 2.5, length.out = L)
  X <- matrix(rbinom(N * L, 1, plogis(outer(rnorm(N), d, "-"))), N, L,
              dimnames = list(NULL, sprintf("I%02d", 1:L)))
  fit <- rasch(X)
  dt <- dimensionality_test(fit, items_positive = sprintf("I%02d", 1:10),
                            items_negative = sprintf("I%02d", 11:20),
                            min_score_points = 2)
  expect_false(dt$multidimensional)
  expect_null(dt$paired_t)
  smd <- dt$subset_mean_difference
  expect_true(is.list(smd))
  expect_true(is.finite(smd$mean_difference) && smd$mean_difference < 0)
  expect_match(smd$note, "no test of this mean")
  expect_false(any(c("t", "df", "p") %in% names(smd)))
  # the description is printed, so it is readable without the object
  out <- capture.output(print(dt))
  expect_match(paste(out, collapse = "\n"),
               "Mean difference between subset estimates: .*descriptive, not a test")
  expect_false(any(grepl("Paired t", out, fixed = TRUE)))
})

test_that("a saved subtest with the superseded paired t-test is refused", {
  # The displays print a paired t of the subset means whenever the field is
  # present, so a result saved before the field was withdrawn must not be
  # restored beside the current verdict.
  set.seed(2); N <- 200; L <- 10
  X <- matrix(rbinom(N * L, 1, plogis(outer(
    rnorm(N), seq(-1.5, 1.5, length.out = L), "-"))), N, L,
    dimnames = list(NULL, sprintf("I%02d", 1:L)))
  fit <- rasch(X)
  dt <- dimensionality_test(fit, items_positive = sprintf("I%02d", 1:5),
                            items_negative = sprintf("I%02d", 6:10),
                            min_score_points = 2)
  expect_identical(dt$algorithm, "person-subset-comparison-1")
  expect_no_error(.validate_dimensionality_test(dt, fit))
  old <- dt
  old$algorithm <- NULL
  old$subset_mean_difference <- NULL
  old$paired_t <- list(mean_difference = -0.5, t = -17.4, df = 290,
                       p = 1.4e-34)
  attr(old, "result_signature") <- NULL
  attr(old, "result_signature") <- .fit_boot_md5(old)
  expect_error(.validate_dimensionality_test(old, fit), "superseded")
})

test_that("a saved project keeps its data when its subtest is superseded", {
  # The stamp must not cost the user the rest of the analysis: the superseded
  # result is dropped with a warning, as for every other superseded tag, and
  # the data, fits and history reopen.
  d <- simulate_rasch(200, 8, seed = 91)
  fit <- rasch(d)
  current <- dimensionality_test(fit, items_positive = sprintf("I%02d", 1:4),
                                 items_negative = sprintf("I%02d", 5:8),
                                 min_score_points = 2)
  expect_identical(current$algorithm, "person-subset-comparison-1")
  old <- current
  old$algorithm <- NULL
  old$subset_mean_difference <- NULL
  old$paired_t <- list(mean_difference = -0.4, t = -6.1, df = 180, p = 5e-09)
  attr(old, "result_signature") <- NULL
  attr(old, "result_signature") <- .fit_boot_md5(old)
  project <- .seal_app_project(list(
    format = "rasch-shiny-project", schema = 2L,
    model_type = "rasch", data = as.data.frame(d), base_fit = fit,
    rasch_steps = list(), btl_steps = list(), kept_fits = list(base = fit),
    settings = list(), results = list(subtest = old,
                                      display = list(item = "I01"))))
  path <- tempfile(fileext = ".rasch")
  on.exit(unlink(path), add = TRUE)
  saveRDS(project, path)
  expect_error(.save_app_project(project, path), "superseded person-subset")
  expect_warning(restored <- .read_app_project(path),
                 "paired t-test of the subset means.*omitted")
  expect_null(restored$results$subtest)
  for (field in c("data", "base_fit", "rasch_steps", "btl_steps", "kept_fits",
                  "settings"))
    expect_identical(restored[[field]], project[[field]])
  expect_identical(restored$results$display, project$results$display)
  expect_match(attr(restored, "rasch_project_legacy_dropped"),
               "person-subset dimensionality test")
  expect_no_error(.validate_app_project(restored))
  expect_no_error(.save_app_project(restored, path))
  expect_no_warning(.read_app_project(path))

  changed <- project
  changed$results$subtest$prop_significant <- 0.9
  saveRDS(changed, path)
  expect_error(.read_app_project(path), "changed since they were saved")

  project$results$subtest <- current
  project <- .seal_app_project(project)
  expect_no_error(.save_app_project(project, path))
  expect_no_warning(restored <- .read_app_project(path))
  expect_identical(restored$results$subtest, current)

  project$schema <- 1L
  project$binding <- NULL
  project$results$subtest <- old
  saveRDS(project, path)
  expect_warning(restored <- .read_app_project(path), "schema-1")
  expect_null(restored$results$subtest)
  expect_identical(restored$base_fit, fit)
  expect_no_error(.validate_app_project(restored))
})
