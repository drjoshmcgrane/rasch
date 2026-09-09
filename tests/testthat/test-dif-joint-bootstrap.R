test_that("conditional bootstrap replays joint incomplete-panel adjustment", {
  set.seed(9162)
  n <- 400L
  A <- rep(c("a0", "a1"), each = n / 2)
  B <- ifelse(runif(n) < ifelse(A == "a0", .2, .8), "b1", "b0")
  id <- c(seq_len(n), which(A == "a1"))
  factors <- data.frame(A = A[id], B = B[id],
    occasion = c(rep("pre", n), rep("post", sum(A == "a1"))))
  theta <- rnorm(n)
  X <- sapply(seq(-1, 1, length.out = 6), function(b)
    rbinom(length(id), 1, plogis(theta[id] - b)))
  colnames(X) <- paste0("I", seq_len(ncol(X)))
  fit <- rasch(X, id = id, factors = factors)
  da <- dif_anova(fit, within = "occasion", n_groups = 2)
  expect_match(da$between_covariance, "CR3")
  db <- suppressWarnings(dif_bootstrap(fit, da, B = 3L, workers = 1L,
                                      seed = 9163))
  expect_equal(db$B_used, 3L)
  expect_equal(nrow(db$replicates$p), 3L)
  expect_true(all(is.finite(db$replicates$p)))
  expect_no_error(.validate_dif_bootstrap(db, fit, da))
})
