test_that("unavailable total-fit probability is explicit in the summary table", {
  set.seed(1)
  X <- matrix(rbinom(400L * 5L, 1, 0.5), 400L, 5L)
  colnames(X) <- paste0("I", seq_len(ncol(X)))
  fit <- rasch(X, id = rep(seq_len(200L), 2L))

  tab <- fit_summary_table(fit)
  row <- tab$statistic == "Approximate asymptotic item-trait probability"
  expect_identical(tab$value[row], "unavailable")
  expect_false(any(grepl("p = \\)$", tab$value, perl = TRUE)))
})
