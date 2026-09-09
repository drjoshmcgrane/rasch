# Phase-B capabilities: response-dependence magnitude, the spread LUB
# screen, multidimensionality magnitude, tailored analysis for guessing,
# traditional statistics, rack/stack reshaping, and the new displays.

test_that("dependence_magnitude recovers a simulated dichotomous d", {
  set.seed(3); N <- 1500; d_true <- 1
  d0 <- seq(-1.5, 1.5, length.out = 8)
  th <- rnorm(N)
  X <- matrix(rbinom(N * 8, 1, plogis(outer(th, d0, "-"))), N, 8)
  # item 5 depends on item 4: its difficulty shifts -d when x4 = 1, +d when 0
  X[, 5] <- rbinom(N, 1, plogis(th - (d0[5] + d_true * (1 - 2 * X[, 4]))))
  colnames(X) <- paste0("I", 1:8)
  fit <- rasch(X)
  dm <- dependence_magnitude(fit, dependent = "I5", independent = "I4")
  expect_lt(abs(dm$d - d_true), 0.35)
  expect_lt(dm$p, 0.001)
  expect_true(is.infinite(dm$df))
  expect_equal(dm$p, 2 * stats::pt(-abs(dm$t), dm$df))
  expect_equal(nrow(dm$thresholds), 1)
  # resolved items replace the originals in the refit
  expect_false(any(c("I4", "I5") %in% dm$refit$items$item))
  expect_true(all(c("I5|I4=0", "I5|I4=1") %in% dm$refit$items$item))
  # an independent pair carries no signal
  dm0 <- dependence_magnitude(fit, dependent = "I7", independent = "I2")
  expect_lt(abs(dm0$d), 0.3)
  expect_gt(dm0$p, 0.01)
  # guards
  expect_error(dependence_magnitude(fit, "I5", "I5"), "different items")
})

test_that("dependence_magnitude se uses the joint covariance of the refit", {
  set.seed(11); N <- 900
  d0 <- seq(-1.2, 1.2, length.out = 10)
  X <- matrix(rbinom(N * 10, 1, plogis(outer(rnorm(N), d0, "-"))), N, 10)
  colnames(X) <- paste0("I", 1:10)
  dm <- dependence_magnitude(rasch(X), dependent = "I6", independent = "I5")
  rf <- dm$refit
  thr <- rf$thresholds
  item_of <- rf$items$item[thr$item]
  lo <- thr[item_of == "I6|I5=0" & thr$k == 1, ]
  hi <- thr[item_of == "I6|I5=1" & thr$k == 1, ]
  cv <- rf$est$cov_tau
  # the reported se is the delta-method contrast from the sandwich
  # covariance, not the independence pooling of A&M eqs. 24.9-24.11
  se_joint <- sqrt((cv[lo$id, lo$id] + cv[hi$id, hi$id] -
                    2 * cv[lo$id, hi$id]) / 4)
  se_pooled <- sqrt((lo$se^2 + hi$se^2) / 4)
  expect_equal(dm$se, se_joint, tolerance = 1e-10)
  # the resolved locations are negatively correlated through the shared
  # comparators, so the joint se must exceed the independence pooling
  expect_gt(dm$se, se_pooled)
  expect_equal(dm$thresholds$se_k, se_joint, tolerance = 1e-10)
})

test_that("dependence_magnitude withholds inference from an invalid covariance", {
  set.seed(111)
  X <- matrix(rbinom(700 * 8, 1, .5), 700, 8,
              dimnames = list(NULL, paste0("I", 1:8)))
  fit <- rasch(X)
  good <- dependence_magnitude(fit, dependent = "I6", independent = "I5")
  bad_refit <- good$refit
  bad_refit$est$cov_tau[1, 1] <- -1e6

  guarded <- testthat::with_mocked_bindings(
    dependence_magnitude(fit, dependent = "I6", independent = "I5"),
    .rasch_refit = function(...) bad_refit,
    .package = "rasch")
  expect_equal(guarded$d, good$d)
  expect_true(is.na(guarded$se))
  expect_true(is.na(guarded$p))
  expect_true(all(is.na(guarded$thresholds$se_k)))
  expect_match(guarded$note, "not positive semidefinite")

  # A clustered resolved covariance must carry its finite reference through
  # the public magnitude object and probability.
  clustered_refit <- good$refit
  clustered_refit$est$cluster_support <- list(
    repeated = TRUE, n = 11L, effective = 10.5)
  clustered_refit$est$cluster_inference <- TRUE
  clustered <- testthat::with_mocked_bindings(
    dependence_magnitude(fit, dependent = "I6", independent = "I5"),
    .rasch_refit = function(...) clustered_refit,
    .package = "rasch")
  expect_equal(clustered$df, 10)
  expect_equal(clustered$p,
               2 * stats::pt(-abs(clustered$t), df = clustered$df))

  unsupported_refit <- good$refit
  unsupported_refit$est$cluster_support <- list(
    repeated = FALSE, n = 16L, effective = 16)
  unsupported_refit$est$cluster_inference <- FALSE
  unsupported <- testthat::with_mocked_bindings(
    dependence_magnitude(fit, dependent = "I6", independent = "I5"),
    .rasch_refit = function(...) unsupported_refit,
    .package = "rasch")
  expect_true(is.na(unsupported$se))
  expect_true(is.na(unsupported$t))
  expect_true(is.na(unsupported$p))
  expect_true(all(is.na(unsupported$thresholds$se_k)))
  expect_match(unsupported$note, "independent.*support")
})

test_that("dependence resolution retains controls and refuses constrained polytomous thresholds", {
  set.seed(12)
  X <- matrix(rbinom(1000 * 8, 1, .5), 1000, 8,
              dimnames = list(NULL, paste0("I", 1:8)))
  f <- rasch(X, n_groups = 7, maxit = 75, tol = 1e-9)
  dm <- dependence_magnitude(f, "I5", "I4")
  expect_equal(dm$refit$refit_spec$n_groups, 7)
  expect_equal(dm$refit$refit_spec$maxit, 75)

  d <- simulate_rasch(500, 7, model = "RSM", n_categories = 4, seed = 13)
  fr <- rasch(d, model = "RSM", id = "id")
  expect_error(dependence_magnitude(fr, "I04", "I03"),
               "unconstrained PCM")
})

test_that("spread_test flags a dependent subtest by the LUB", {
  set.seed(5); N <- 800
  d0 <- rep(c(-0.5, 0, 0.5), 3)[1:9]
  th <- rnorm(N)
  X <- matrix(rbinom(N * 9, 1, plogis(outer(th, d0, "-"))), N, 9)
  # items 2 and 3 strongly follow item 1
  X[, 2] <- ifelse(runif(N) < 0.9, X[, 1], X[, 2])
  X[, 3] <- ifelse(runif(N) < 0.9, X[, 1], X[, 3])
  colnames(X) <- paste0("I", 1:9)
  fit2 <- combine_items(rasch(X), list(c("I1", "I2", "I3"),
                                       c("I4", "I5", "I6"),
                                       c("I7", "I8", "I9")))
  st <- spread_test(fit2)
  expect_equal(nrow(st), 3)
  expect_true(all(st$eligible))
  expect_equal(st$lub, rep(0.55, 3))
  expect_equal(st$below_bound, st$spread < st$lub)
  expect_true(all(is.infinite(st$df)))
  expect_equal(st$p, stats::pt(st$t, df = st$df), tolerance = 1e-12)
  expect_equal(st$p_adj, p.adjust(st$p, method = "holm"))
  expect_equal(st$dependent, st$p_adj < 0.05)
  expect_equal(attr(st, "alpha"), 0.05)
  expect_identical(attr(st, "p_adjust"), "holm")
  dep_row <- grep("I1", st$item)
  ind_rows <- setdiff(seq_len(3), dep_row)
  expect_true(st$dependent[dep_row])
  expect_lt(st$spread[dep_row], 0.3)
  expect_true(all(st$spread[ind_rows] > st$spread[dep_row]))

  # spread_test() reaches the internal fit directly so it can pass person IDs
  pc_zero <- pcml_pc(fit2$X)
  pc_zero$components$spread_se[1] <- 0
  st_zero <- testthat::with_mocked_bindings(
    spread_test(fit2),
    .pcml_pc_fit = function(...) pc_zero,
    .package = "rasch")
  expect_true(is.na(st_zero$t[1]))
  expect_true(is.na(st_zero$p_adj[1]))

  pc_cluster <- pc_zero
  pc_cluster$components$spread_se[1] <- 0.2
  pc_cluster$cluster_support <- list(
    repeated = TRUE, n = 12L, effective = 11.5)
  pc_cluster$cluster_inference <- TRUE
  st_cluster <- testthat::with_mocked_bindings(
    spread_test(fit2),
    .pcml_pc_fit = function(...) pc_cluster,
    .package = "rasch")
  expect_equal(st_cluster$df[st_cluster$eligible],
               rep(11, sum(st_cluster$eligible)))
  expect_equal(st_cluster$p[st_cluster$eligible],
               stats::pt(st_cluster$t[st_cluster$eligible], 11))

  pc_unsupported <- pc_zero
  pc_unsupported$components$spread_se[1] <- 0.2
  pc_unsupported$cluster_support <- list(
    repeated = FALSE, n = 16L, effective = 16)
  pc_unsupported$cluster_inference <- FALSE
  st_unsupported <- testthat::with_mocked_bindings(
    spread_test(fit2),
    .pcml_pc_fit = function(...) pc_unsupported,
    .package = "rasch")
  expect_true(all(is.na(st_unsupported$df[st_unsupported$eligible])))
  expect_true(all(is.na(st_unsupported$se)))
  expect_true(all(is.na(st_unsupported$t)))
  expect_true(all(is.na(st_unsupported$p[st_unsupported$eligible])))
  expect_match(attr(st_unsupported, "note"),
               "independent-person support")

  pcm <- rasch(simulate_rasch(400, 5, model = "PCM", n_categories = 3,
                              seed = 52), id = "id")
  expect_error(spread_test(pcm), "no recorded superitems")
  pcm$subtest_map <- list(I01 = c("source1", "source2"))
  pcm$subtest_binary <- c(I01 = FALSE)
  ps <- spread_test(pcm)
  expect_false(ps$eligible)
  expect_true(is.na(ps$lub) && is.na(ps$below_bound) &&
                is.na(ps$p) && is.na(ps$p_adj) && is.na(ps$dependent))
  expect_error(spread_test(fit2, alpha = 1), "strictly between")
  expect_error(spread_test(fit2, p_adjust = "invalid"), "p.adjust.methods")
})

test_that("dimensionality_magnitude reproduces the Andrich (2016) block", {
  set.seed(7); N <- 800
  common <- rnorm(N); u1 <- rnorm(N); u2 <- rnorm(N); c_true <- 0.8
  d0 <- rep(seq(-1, 1, length.out = 10), 2)
  X <- sapply(1:20, function(i) rbinom(N, 1,
    plogis(common + c_true * (if (i <= 10) u1 else u2) - d0[i])))
  colnames(X) <- paste0("I", 1:20)
  fit <- rasch(X)
  dm <- dimensionality_magnitude(fit, list(paste0("I", 1:10), paste0("I", 11:20)))
  tab <- dm$table
  # subtest reliability drops; c is recovered roughly; rho = 1/(1+c2); A = S/(S+c2)
  expect_true(all(tab$subtest < tab$run1))
  expect_lt(abs(tab$c[tab$index == "PSI"] - c_true), 0.35)
  expect_equal(tab$rho, 1 / (1 + tab$c2), tolerance = 1e-10)
  expect_equal(tab$A, 2 / (2 + tab$c2), tolerance = 1e-10)
  # a unidimensional scale yields c near zero and rho near 1
  X1 <- sapply(1:20, function(i) rbinom(N, 1, plogis(common - d0[i])))
  colnames(X1) <- paste0("I", 1:20)
  dm1 <- dimensionality_magnitude(rasch(X1),
                                  list(paste0("I", 1:10), paste0("I", 11:20)))
  expect_lt(dm1$table$c2[1], 0.25)
  expect_gt(dm1$table$rho[1], 0.8)
  # guard: every item must be assigned
  expect_error(dimensionality_magnitude(fit, list(paste0("I", 1:6))),
               "at least two")
})

test_that("dimensionality magnitude withholds an unusable reliability ratio", {
  fit <- rasch(simulate_rasch(300, 8, seed = 121), id = "id")
  fit$alpha$alpha <- -0.1
  refit <- fit
  refit$alpha$alpha <- 0.5
  out <- testthat::with_mocked_bindings(
    dimensionality_magnitude(
      fit, list(sprintf("I%02d", 1:4), sprintf("I%02d", 5:8))),
    combine_items = function(...) refit,
    .package = "rasch"
  )
  alpha_row <- out$table$index == "alpha"
  expect_true(all(is.na(unlist(
    out$table[alpha_row, c("c2", "c", "rho", "A")], use.names = FALSE))))
  expect_true(all(is.finite(unlist(
    out$table[!alpha_row, c("c2", "c", "rho", "A")], use.names = FALSE))))
})

test_that("tailored_analysis shows the guessing signature", {
  set.seed(11); N <- 900
  d0 <- seq(-2, 2.5, length.out = 10); th <- rnorm(N)
  P <- 0.25 + 0.75 * plogis(outer(th, d0, "-"))   # guessing floor 0.25
  X <- matrix(rbinom(N * 10, 1, P), N, 10)
  colnames(X) <- paste0("I", 1:10)
  ta <- tailored_analysis(rasch(X), chance = 0.25)
  expect_gt(ta$n_removed, 50)
  expect_identical(ta$se_method, "none")
  expect_true(all(is.na(ta$table$se)))
  expect_true(all(is.na(ta$table$p_adj)))
  # the hardest items become harder under tailoring, on the common origin,
  # and clearly more so than the easy items
  hard <- order(ta$table$initial, decreasing = TRUE)[1:2]
  easy <- order(ta$table$initial)[1:2]
  expect_gt(mean(ta$table$shift[hard]), 0.1)
  expect_gt(mean(ta$table$shift[hard]), mean(ta$table$shift[easy]) + 0.1)
  expect_lt(mean(abs(ta$table$shift[easy])), 0.25)
  # step-4 fit: items fixed (se 0), persons re-estimated on original data
  expect_true(all(ta$anchored$thresholds$se == 0))
  expect_equal(nrow(ta$anchored$person), N)
  # no-guessing data keeps the difficult items in place
  X0 <- matrix(rbinom(N * 10, 1, plogis(outer(th, d0, "-"))), N, 10)
  colnames(X0) <- paste0("I", 1:10)
  ta0 <- tailored_analysis(rasch(X0), chance = 0.25)
  expect_lt(mean(ta0$table$shift[order(ta0$table$initial,
                                       decreasing = TRUE)[1:2]]), 0.3)
  expect_error(tailored_analysis(rasch(X), chance = c(0.2, 0.25)),
               "chance.*probability")
  expect_error(tailored_analysis(rasch(X), chance = 1),
               "strictly between")
  expect_error(tailored_analysis(rasch(X), chance = 0.25,
                                 se_method = "bootstrap", boot_reps = 49),
               "whole number")
  expect_error(tailored_analysis(rasch(X), chance = 0.25,
                                 se_method = "bootstrap", seed = 1.5),
               "whole number")
  expect_error(tailored_analysis(
    rasch(X), chance = 0.25,
    anchor_items = matrix(c("I1", "I2"), ncol = 1L)),
    "ordinary vector")
  expect_error(tailored_analysis(
    rasch(X), chance = 0.25, anchor_items = c("I1", "I1")),
    "named more than once")
})

test_that("tailored_analysis bootstrap repeats the complete procedure", {
  # At equality the package's strict p_adj < .05 rule still has no rejection
  # region; 40m, rather than 40m - 1, is the first usable request.
  expect_equal(rasch:::.tailored_boot_floor(79L, 2L), 0.05)
  expect_lt(rasch:::.tailored_boot_floor(80L, 2L), 0.05)

  set.seed(812)
  N <- 220L
  d0 <- seq(-1.5, 1.8, length.out = 6)
  th <- rnorm(N)
  P <- 0.25 + 0.75 * plogis(outer(th, d0, "-"))
  X <- matrix(rbinom(N * 6, 1, P), N, 6,
              dimnames = list(NULL, paste0("I", 1:6)))
  # at 50 replicates and 6 items the Holm-adjusted p floor (2m/(B+1))
  # exceeds 0.05, and tailored_analysis must SAY so -- the warning is part
  # of the contract being tested here, not noise to suppress
  fit <- rasch(X)
  set.seed(99)
  caller_stream <- .Random.seed
  expect_warning(
    ta <- tailored_analysis(fit, chance = 0.25,
                            se_method = "bootstrap", boot_reps = 50,
                            seed = 812),
    "smallest achievable")
  expect_identical(.Random.seed, caller_stream)
  expect_identical(ta$seed, 812L)
  expect_identical(ta$se_method, "bootstrap")
  expect_gte(ta$boot_reps_used, 45L)
  expect_identical(ta$boot_minimum_usable, 45L)
  expect_identical(ta$boot_reps,
                   ta$boot_reps_used + ta$boot_reps_nonconverged +
                     ta$boot_reps_errors)
  expect_true(all(is.finite(ta$table$se)))
  expect_true(all(is.finite(ta$table$ci_low)))
  expect_true(all(is.finite(ta$table$ci_high)))
  expect_true(all(is.finite(ta$table$p_adj)))

  original_rasch <- rasch
  calls <- 0L
  refused <- tryCatch(suppressWarnings(with_mocked_bindings(
    tailored_analysis(fit, chance = 0.25,
                      se_method = "bootstrap", boot_reps = 50),
    rasch = function(...) {
      calls <<- calls + 1L
      if (calls <= 2L) original_rasch(...) else stop("forced inner failure")
    },
    .package = "rasch")), error = identity)
  expect_s3_class(refused, "rasch_fit_bootstrap_refusal")
  expect_identical(refused$B, 50L)
  expect_identical(refused$B_used, 0L)
  expect_identical(refused$B_nonconverged, 0L)
  expect_identical(refused$B_errors, 50L)
})

test_that("tailored bootstrap keeps repeated-person rows together", {
  set.seed(813)
  id <- c("A", "A", "B", "B", "B", "C", "C")
  bs <- rasch:::.tailored_boot_rows(id)
  copies <- split(id[bs$rows], bs$id)
  expect_true(all(vapply(copies, function(x) length(unique(x)) == 1L,
                         logical(1))))
  expect_equal(lengths(copies),
               vapply(copies, function(x) sum(id == x[1]), integer(1)))
})

test_that("ctt_table reports the classical companions", {
  set.seed(13)
  d0 <- seq(-1.5, 1.5, length.out = 8)
  X <- matrix(rbinom(500 * 8, 1, plogis(outer(rnorm(500), d0, "-"))), 500, 8)
  colnames(X) <- paste0("I", 1:8)
  fit <- rasch(X)
  ct <- ctt_table(fit)
  expect_equal(ct$alpha, fit$alpha$alpha, tolerance = 1e-10)
  expect_equal(ct$table$facility, unname(colMeans(X)), tolerance = 1e-10)
  expect_true(all(ct$table$item_total > 0.2))
  expect_true(all(ct$table$item_rest < ct$table$item_total))
  expect_true(all(ct$table$di > 0.1))
  expect_equal(ct$sem, ct$sd * sqrt(1 - ct$alpha), tolerance = 1e-10)
  # facility ordering follows difficulty
  expect_equal(order(ct$table$facility, decreasing = TRUE), order(d0))

  virtual <- fit
  class(virtual) <- c("rasch_mfrm", class(virtual))
  expect_error(ctt_table(virtual), "several frame or facet response cells")
  expect_error(guttman_table(virtual), "several frame or facet response cells")
})

test_that("a negative alpha is not used as reliability for the CTT SEM", {
  set.seed(122)
  N <- 1000L
  X <- matrix(0L, N, 4L,
              dimnames = list(NULL, paste0("I", 1:4)))
  chosen <- sample.int(4L, N, replace = TRUE)
  X[cbind(seq_len(N), chosen)] <- 1L
  noisy <- runif(N) < 0.1
  X[noisy, ] <- matrix(rbinom(sum(noisy) * 4L, 1L, 0.5), sum(noisy), 4L)
  fit <- structure(list(X = X, m = rep(1L, 4L)), class = "rasch")
  out <- ctt_table(fit)
  expect_lt(out$alpha, 0)
  expect_true(is.na(out$sem))
  expect_match(out$note, "SEM withheld.*alpha is negative")
})

test_that("rack_data and stack_data reshape repeated measurements", {
  d <- expand.grid(pid = 1:50, t = 1:3)
  d$Q1 <- rbinom(150, 1, 0.6); d$Q2 <- rbinom(150, 2, 0.5)
  r <- rack_data(d, person = "pid", time = "t", items = c("Q1", "Q2"))
  expect_equal(dim(r), c(50, 1 + 6))
  expect_true(all(c("Q1@1", "Q2@3") %in% names(r)))
  expect_equal(r$`Q1@2`, d$Q1[d$t == 2][match(r$id, d$pid[d$t == 2])])
  s <- stack_data(d, person = "pid", time = "t", items = c("Q1", "Q2"))
  expect_equal(nrow(s), 150)
  expect_true(is.factor(s$time))
  expect_equal(s$Q1, d$Q1)
  expect_identical(s$id, d$pid)
  expect_equal(length(unique(s$row_id)), nrow(s))
  # duplicate person-time rows are an error when racking
  expect_error(rack_data(rbind(d, d[1, ]), "pid", "t", c("Q1", "Q2")),
               "more than one row")
  expect_error(stack_data(rbind(d, d[1, ]), "pid", "t", c("Q1", "Q2")),
               "same time point")
  expect_error(rack_data(d, "pid", "pid", c("Q1", "Q2")), "distinct")
  expect_error(stack_data(d, "pid", "pid", c("Q1", "Q2")), "distinct")
  dn <- d; names(dn)[1L] <- "1"
  expect_error(rack_data(dn, 1, "1", c("Q1", "Q2")), "character string")
  expect_error(stack_data(dn, 1, "1", c("Q1", "Q2")), "character string")
  expect_error(rack_data(d, "pid", "t", factor(c("Q1", "Q2"))),
               "at least one item column")
  expect_error(stack_data(d, "pid", "t", factor(c("Q1", "Q2"))),
               "at least one item column")
  dd <- data.frame(pid = rep(1:2, 2),
                   t = rep(as.Date(c("2020-01-01", "2020-01-02")), each = 2),
                   Q = 1:4)
  expect_true(all(c("Q@2020-01-01", "Q@2020-01-02") %in%
                    names(rack_data(dd, "pid", "t", "Q"))))

  # Separators in source values cannot merge distinct person-time cells,
  # and generated output names cannot silently duplicate one another.
  dc <- data.frame(pid = c("a@b", "a"), t = c("c", "b@c"), Q = 0:1)
  sc <- stack_data(dc, "pid", "t", "Q")
  expect_equal(length(unique(sc$row_id)), 2L)
  dr <- data.frame(pid = 1:2, t = c("2", "1@2"), Q = 0:1, `Q@1` = 1:0,
                   check.names = FALSE)
  expect_error(rack_data(dr, "pid", "t", c("Q", "Q@1")),
               "not unique")
  ds <- data.frame(pid = 1:2, t = 1:2, id = 0:1)
  expect_error(stack_data(ds, "pid", "t", "id"), "reserved")

  # Whitespace is not part of an identifier or occasion label.  The values
  # used to validate the design must also be the values used to align it.
  dw <- data.frame(pid = c(" P1", "P2 ", "P1 ", " P2"),
                   t = c(" T1", "T1 ", " T2", "T2 "), Q = 1:4)
  rw <- rack_data(dw, "pid", "t", "Q")
  expect_identical(rw$id, c("P1", "P2"))
  expect_identical(names(rw), c("id", "Q@T1", "Q@T2"))
  expect_equal(rw$`Q@T1`, 1:2)
  expect_equal(rw$`Q@T2`, 3:4)
  sw <- stack_data(dw, "pid", "t", "Q")
  expect_identical(sw$id, c("P1", "P2", "P1", "P2"))
  expect_identical(levels(sw$time), c("T1", "T2"))
  expect_error(stack_data(transform(dw, pid = c("P1", " P1 ", "P2", "P3"),
                                    t = c("T1", " T1 ", "T2", "T2")),
                          "pid", "t", "Q"), "same time point")
})

test_that("the new displays draw without error", {
  set.seed(17)
  tau <- list(c(-1, 0.2), c(-0.5, 0.6), c(-1.2, 0), c(0, 1), c(-0.8, 0.8),
              c(-0.3, 0.9))
  X <- sapply(tau, function(tt) vapply(rnorm(300), function(b)
    sample(0:2, 1, prob = item_moments(b, tt)$P), 0L))
  colnames(X) <- paste0("Q", 1:6)
  fit <- rasch(X)
  pdf(NULL); on.exit(dev.off())
  expect_no_error(plot_ccc(fit, "Q1", observed = TRUE))
  expect_no_error(plot_threshold_prob(fit, "Q1", observed = TRUE))
  expect_no_error(plot_pcc(fit, person = 1))
  expect_no_error(plot_resid_dist(fit, "items"))
  expect_no_error(plot_resid_dist(fit, "persons", "natural"))
  # paired t-test is part of the dimensionality report
  dt <- dimensionality_test(fit, min_score_points = 2)
  expect_true(is.list(dt$paired_t))
  expect_true(is.finite(dt$paired_t$p))
})

test_that("rasch(pc_components) routes estimation through pcml_pc", {
  set.seed(19)
  tau <- list(c(-1, 0, 1), c(-0.6, 0.2, 1.1), c(-1.3, -0.2, 0.9),
              c(-0.9, 0.1, 1.2))
  X <- sapply(tau, function(tt) vapply(rnorm(400), function(b)
    sample(0:3, 1, prob = item_moments(b, tt)$P), 0L))
  colnames(X) <- paste0("Q", 1:4)
  fit <- rasch(X, pc_components = 2)
  expect_s3_class(fit, "rasch")
  expect_true(!is.null(fit$est$components))
  expect_true(all(is.na(fit$est$components$skewness)))   # rank 2: no skewness
  expect_true(any(grepl("principal component", fit$notes)))
  # equal spread within each item: threshold spacings constant
  for (tl in fit$tau_list) expect_lt(diff(range(diff(tl))), 1e-8)
  # full-rank pc reproduces free estimation exactly through four thresholds
  fit4 <- rasch(X, pc_components = 4)
  free <- rasch(X)
  expect_equal(fit4$thresholds$tau, free$thresholds$tau, tolerance = 1e-5)
  # guards
  expect_error(rasch(X, model = "RSM", pc_components = 2), "PCM only")
  expect_error(rasch(X, pc_components = 2,
                     anchors = data.frame(item = "Q1", k = 1, tau = 0)),
               "anchors")
})

test_that("tailored analysis refuses an externally fixed origin", {
  set.seed(47)
  X <- matrix(rbinom(500 * 8, 1, .5), 500, 8,
              dimnames = list(NULL, paste0("I", 1:8)))
  f <- rasch(X, anchors = data.frame(item = "I1", k = 1, tau = 0))
  expect_error(tailored_analysis(f), "unanchored calibration")
})

test_that("inferential magnitude procedures refuse unconverged calibrations", {
  set.seed(49)
  X <- matrix(rbinom(500 * 8, 1, .5), 500, 8,
              dimnames = list(NULL, paste0("I", 1:8)))
  f <- rasch(X)
  f$est$converged <- FALSE
  expect_error(dependence_magnitude(f, "I2", "I1"), "did not converge")
  expect_error(spread_test(f), "did not converge")
  expect_error(dimensionality_test(f), "did not converge")
  expect_error(dimensionality_magnitude(
    f, list(paste0("I", 1:4), paste0("I", 5:8))), "did not converge")
  expect_error(tailored_analysis(f), "did not converge")
  expect_error(residual_correlations(f), "did not converge")
  expect_error(residual_pca(f), "did not converge")
})

test_that("dependence resolution refuses a generated-name collision", {
  set.seed(50)
  X <- matrix(rbinom(800 * 5, 1, .5), 800, 5,
              dimnames = list(NULL, c("A", "B", "B|A=0", "C", "D")))
  f <- rasch(X)
  expect_error(dependence_magnitude(f, "B", "A"),
               "resolved-item name already exists")
})

test_that("dependence magnitude routes MFRM virtual items through a PCM", {
  set.seed(51)
  d <- expand.grid(person = sprintf("P%03d", 1:220),
                   item = c("A", "B", "C", "D"),
                   rater = c("R1", "R2"), stringsAsFactors = FALSE)
  d$score <- stats::rbinom(nrow(d), 1, .5)
  f <- rasch_mfrm(d, "person", "item", "score", "rater")
  vm <- f$virtual_map$vkey[f$virtual_map$rater == "R1"]
  z <- dependence_magnitude(f, dependent = vm[2], independent = vm[1])
  expect_s3_class(z, "rasch_dependence")
  expect_identical(z$refit$model, "PCM")

  ef <- structure(f, class = c("rasch_efrm", class(f)))
  expect_error(dependence_magnitude(ef, vm[2], vm[1]),
               "mutually exclusive EFRM")
})
