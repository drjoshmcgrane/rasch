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
  expect_equal(st$lub, rep(0.55, 3))
  dep_row <- grep("I1", st$item)
  ind_rows <- setdiff(seq_len(3), dep_row)
  expect_true(st$dependent[dep_row])
  expect_lt(st$spread[dep_row], 0.3)
  expect_true(all(st$spread[ind_rows] > st$spread[dep_row]))
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

test_that("tailored_analysis shows the guessing signature", {
  set.seed(11); N <- 900
  d0 <- seq(-2, 2.5, length.out = 10); th <- rnorm(N)
  P <- 0.25 + 0.75 * plogis(outer(th, d0, "-"))   # guessing floor 0.25
  X <- matrix(rbinom(N * 10, 1, P), N, 10)
  colnames(X) <- paste0("I", 1:10)
  ta <- tailored_analysis(rasch(X), chance = 0.25)
  expect_gt(ta$n_removed, 50)
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
  # duplicate person-time rows are an error when racking
  expect_error(rack_data(rbind(d, d[1, ]), "pid", "t", c("Q1", "Q2")),
               "more than one row")
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
  dt <- dimensionality_test(fit)
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
  # full-rank pc reproduces free estimation exactly at <= 3 thresholds
  fit4 <- rasch(X, pc_components = 4)
  free <- rasch(X)
  expect_equal(fit4$thresholds$tau, free$thresholds$tau, tolerance = 1e-5)
  # guards
  expect_error(rasch(X, model = "RSM", pc_components = 2), "PCM only")
  expect_error(rasch(X, pc_components = 2,
                     anchors = data.frame(item = "Q1", k = 1, tau = 0)),
               "anchors")
})
