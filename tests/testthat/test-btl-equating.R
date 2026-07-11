# pkgload::load_all() does not regenerate NAMESPACE, so the freshly added
# print method is not yet in the S3 table; register it for the test session
# exactly as the roxygen @export tag will once NAMESPACE is regenerated.
registerS3method("print", "rasch_btl_equate", print.rasch_btl_equate)

test_that("btl_equate recovers the origin shift and flags nothing when no object drifts", {
  set.seed(20)
  objs <- paste0("O", 1:10)
  beta <- setNames(seq(-2, 2, length.out = 10), objs)

  # two panels judging overlapping object sets: 10 objects in all, 7 common.
  # each panel's sum-zero origin is the mean of its OWN set, so the two
  # origins differ and equating must recover the difference.
  set1 <- paste0("O", 1:9)                 # unique: O8, O9
  set2 <- c(paste0("O", 1:7), "O10")       # unique: O10
  common <- intersect(set1, set2)
  expect_equal(length(common), 7L)

  sim_panel <- function(objset, betas, reps = 60) {
    pr <- t(utils::combn(objset, 2))
    d <- data.frame(a = rep(pr[, 1], each = reps),
                    b = rep(pr[, 2], each = reps), stringsAsFactors = FALSE)
    p <- plogis(betas[d$a] - betas[d$b])
    d$win <- ifelse(runif(nrow(d)) < p, d$a, d$b)
    btl(d, "a", "b", "win")
  }

  f1 <- sim_panel(set1, beta)
  f2 <- sim_panel(set2, beta)
  eq <- btl_equate(f1, f2)

  expect_s3_class(eq, "rasch_btl_equate")
  expect_equal(eq$n_common, 7L)

  # sum-zero origins differ by mean(beta over set2) - mean(beta over set1);
  # loc2 + shift lands fit2 on fit1's scale, so shift ~ that difference
  expected_shift <- mean(beta[set2]) - mean(beta[set1])
  expect_lt(abs(eq$shift - expected_shift), 0.35)

  # essentially nothing drifts (at most the ~alpha rate)
  expect_lte(sum(eq$table$drifting), 1L)

  # fit2's equated locations track the truth (a pure shift of it) ~ perfectly
  truth <- beta[eq$equated$object]
  expect_gt(cor(eq$equated$location, truth), 0.97)
})

test_that("btl_equate flags a planted drift and essentially only that object", {
  set.seed(6)
  objs <- paste0("O", 1:10)
  beta <- setNames(seq(-2, 2, length.out = 10), objs)
  set1 <- paste0("O", 1:9)
  set2 <- c(paste0("O", 1:7), "O10")

  sim_panel <- function(objset, betas, reps = 60) {
    pr <- t(utils::combn(objset, 2))
    d <- data.frame(a = rep(pr[, 1], each = reps),
                    b = rep(pr[, 2], each = reps), stringsAsFactors = FALSE)
    p <- plogis(betas[d$a] - betas[d$b])
    d$win <- ifelse(runif(nrow(d)) < p, d$a, d$b)
    btl(d, "a", "b", "win")
  }

  # panel 2 values one common object (O4) 1.2 logits higher than panel 1 does
  beta2 <- beta
  beta2["O4"] <- beta2["O4"] + 1.2

  f1 <- sim_panel(set1, beta)
  f2 <- sim_panel(set2, beta2)
  eq <- btl_equate(f1, f2)

  drift <- setNames(eq$table$drifting, eq$table$object)
  expect_true(drift["O4"])                       # the planted object flags
  expect_equal(sum(drift[names(drift) != "O4"]), 0L)   # and only that one
  # the planted object dominates the drift statistics
  expect_equal(eq$table$object[which.max(abs(eq$table$t))], "O4")
})

test_that("a bank built from fit1's own objects equates to a zero shift with no drift", {
  set.seed(22)
  objs <- paste0("O", 1:8)
  beta <- setNames(seq(-1.8, 1.8, length.out = 8), objs)
  pr <- t(utils::combn(objs, 2))
  d <- data.frame(a = rep(pr[, 1], each = 50),
                  b = rep(pr[, 2], each = 50), stringsAsFactors = FALSE)
  d$win <- ifelse(runif(nrow(d)) < plogis(beta[d$a] - beta[d$b]), d$a, d$b)
  f1 <- btl(d, "a", "b", "win")

  bank <- data.frame(object = f1$objects$object,
                     location = f1$objects$location,
                     se = f1$objects$se, stringsAsFactors = FALSE)
  eq <- btl_equate(f1, bank)

  expect_lt(abs(eq$shift), 1e-8)
  expect_equal(sum(eq$table$drifting), 0L)
  expect_equal(eq$n_common, nrow(f1$objects))
})

test_that("btl_equate guards fewer than three common objects and non-btl input", {
  set.seed(23)
  objs <- paste0("O", 1:6)
  beta <- setNames(seq(-1.5, 1.5, length.out = 6), objs)
  pr <- t(utils::combn(objs, 2))
  d <- data.frame(a = rep(pr[, 1], each = 40),
                  b = rep(pr[, 2], each = 40), stringsAsFactors = FALSE)
  d$win <- ifelse(runif(nrow(d)) < plogis(beta[d$a] - beta[d$b]), d$a, d$b)
  f1 <- btl(d, "a", "b", "win")

  # a bank sharing only two objects with fit1
  bank2 <- data.frame(object = c("O1", "O2", "Z1", "Z2"),
                      location = c(-1, -0.5, 0.5, 1),
                      se = rep(0.2, 4), stringsAsFactors = FALSE)
  expect_error(btl_equate(f1, bank2), "three common")

  # non-btl fit1
  expect_error(btl_equate(42, bank2), "btl")
  # non-btl, non-bank fit2
  expect_error(btl_equate(f1, 42), "btl fit or a bank")
})

test_that("plot_btl_equate draws without error", {
  set.seed(24)
  objs <- paste0("O", 1:10)
  beta <- setNames(seq(-2, 2, length.out = 10), objs)
  set1 <- paste0("O", 1:9)
  set2 <- c(paste0("O", 1:7), "O10")

  sim_panel <- function(objset, betas, reps = 50) {
    pr <- t(utils::combn(objset, 2))
    d <- data.frame(a = rep(pr[, 1], each = reps),
                    b = rep(pr[, 2], each = reps), stringsAsFactors = FALSE)
    p <- plogis(betas[d$a] - betas[d$b])
    d$win <- ifelse(runif(nrow(d)) < p, d$a, d$b)
    btl(d, "a", "b", "win")
  }
  f1 <- sim_panel(set1, beta)
  f2 <- sim_panel(set2, beta)

  pdf(NULL)
  on.exit(dev.off())
  res <- plot_btl_equate(f1, f2)
  expect_s3_class(res, "rasch_btl_equate")

  # print method exercises its formatting path too
  expect_output(print(res), "Common-object equating")
})
