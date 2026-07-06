simP <- function(theta, tau) { x <- 0:length(tau); p <- exp(x * theta - c(0, cumsum(tau))); p / sum(p) }

test_that("dimensionality test separates 1D from 2D data", {
  set.seed(1); Np <- 1500; L <- 20
  d <- scale(seq(-2, 2, length.out = L), scale = FALSE)[, 1]; th <- rnorm(Np, 0, 1.4)
  X1 <- matrix(rbinom(Np * L, 1, plogis(outer(th, d, "-"))), Np, L); colnames(X1) <- sprintf("U%02d", 1:L)
  dt1 <- dimensionality_test(rasch(X1, model = "PCM"))
  expect_false(dt1$multidimensional)

  set.seed(2)
  thA <- rnorm(Np, 0, 1.4); thB <- 0.3 * thA + sqrt(1 - 0.3^2) * rnorm(Np, 0, 1.4)
  XA <- matrix(rbinom(Np * 10, 1, plogis(outer(thA, d[1:10], "-"))), Np, 10)
  XB <- matrix(rbinom(Np * 10, 1, plogis(outer(thB, d[11:20], "-"))), Np, 10)
  X2 <- cbind(XA, XB); colnames(X2) <- sprintf("D%02d", 1:20)
  dt2 <- dimensionality_test(rasch(X2, model = "PCM"))
  expect_true(dt2$multidimensional)
  expect_gt(dt2$prop_significant, dt1$prop_significant)
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
  expect_true(fit$power_of_fit %in% c("reasonable", "good", "excellent"))
})

test_that("save_outputs writes the full set of tables and plots", {
  set.seed(10); Np <- 300; L <- 6
  d <- scale(seq(-1.5, 1.5, length.out = L), scale = FALSE)[, 1]
  X <- matrix(rbinom(Np * L, 1, plogis(outer(rnorm(Np), d, "-"))), Np, L)
  colnames(X) <- sprintf("I%02d", 1:L)
  fit <- rasch(data.frame(X, g = rep(c("x", "y"), Np / 2)), factors = "g")
  out <- file.path(tempdir(), paste0("rr-test-", as.integer(runif(1, 1, 1e6))))
  files <- save_outputs(fit, out, formats = "png", item_plots = TRUE)
  expect_true(file.exists(file.path(out, "tables", "item_statistics.csv")))
  expect_true(file.exists(file.path(out, "tables", "person_estimates.csv")))
  expect_true(file.exists(file.path(out, "tables", "dif_anova.csv")))
  expect_true(file.exists(file.path(out, "summary.txt")))
  expect_true(file.exists(file.path(out, "plots", "test_information.png")))
  expect_true(file.exists(file.path(out, "plots", "items", "I01_icc.png")))
  expect_gte(length(files), 8 + 8 + 4 * L)
  unlink(out, recursive = TRUE)
})
