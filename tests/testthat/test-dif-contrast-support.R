test_that("planned DIF contrasts cannot renormalise away sparse weighted cells", {
  set.seed(812)
  design <- expand.grid(group = c("A", "B"), stratum = letters[1:5],
    person = seq_len(40))
  theta <- rnorm(nrow(design))
  prob <- plogis(outer(theta, seq(-1, 1, length.out = 6), "-"))
  prob[, 1] <- plogis(theta + 1 - .9 * (design$group == "B"))
  X <- matrix(rbinom(length(prob), 1, prob), nrow(design), 6,
    dimnames = list(NULL, paste0("I", 1:6)))
  sparse <- design$stratum == "e" & design$person > 5
  X[sparse, 1] <- NA
  factors <- design[c("group", "stratum")]
  fit <- rasch(X, factors = factors)
  cellmap <- unique(data.frame(cell = as.character(.factor_cells(factors)),
    factors))
  all_five <- setNames(ifelse(cellmap$group == "B", .2, -.2), cellmap$cell)
  first_four <- all_five
  first_four[cellmap$stratum == "e"] <- 0
  only_a <- all_five
  only_a[cellmap$stratum != "a"] <- 0
  contrasts <- list(all_five = all_five, first_four = first_four,
                    stratum_a = only_a)
  result <- dif_contrasts(fit, items = c("I1", "I2"), contrasts = contrasts)
  tab <- result$table
  omitted <- tab$item == "I1" & tab$contrast == "all_five"
  expect_equal(sum(omitted), 1)
  expect_true(all(is.na(unlist(tab[omitted,
    c("estimate", "se", "statistic", "p", "p_adj", "lower", "upper")]))))
  expect_true(all(is.finite(tab$estimate[!omitted])))
  expect_true(all(is.finite(tab$p[!omitted])))
  expect_match(paste(result$notes, collapse = " "),
    "I1 \\[all_five\\].*required contrast cell")
  expect_identical(result$family_n, 6L)
  expect_equal(tab$p_adj[!omitted],
    p.adjust(tab$p[!omitted], "holm", n = 6))

  # Zero-weight cells need no support, and their deliberately narrower
  # contrast must agree with the separately requested four-stratum analysis.
  control <- dif_contrasts(fit, items = "I1",
    contrasts = list(first_four = first_four))
  at <- tab$item == "I1" & tab$contrast == "first_four"
  expect_equal(tab$estimate[at], control$table$estimate)
  expect_equal(tab$se[at], control$table$se)
  expect_equal(tab$p[at], control$table$p)

  # Automatic marginal contrasts declare the same five-stratum question.
  auto <- dif_contrasts(fit, items = "I1")
  fam <- .dif_contrast_family(factors, cellmap, character())
  group_contrast <- names(fam$meta)[vapply(fam$meta,
    function(x) identical(x$factors, "group"), logical(1))]
  expect_length(group_contrast, 1L)
  expect_true(is.na(auto$table$estimate[auto$table$contrast == group_contrast]))
  expect_true(is.na(auto$table$p_adj[auto$table$contrast == group_contrast]))
})

test_that("fully supported DIF contrasts retain their declared weighting", {
  set.seed(813)
  factors <- expand.grid(group = c("A", "B"), stratum = letters[1:5],
    person = seq_len(40))[c("group", "stratum")]
  X <- matrix(rbinom(nrow(factors) * 6, 1, .5), nrow(factors), 6,
    dimnames = list(NULL, paste0("I", 1:6)))
  fit <- rasch(X, factors = factors)
  grp <- .factor_cells(factors)
  cells <- unique(data.frame(cell = as.character(grp), factors))
  weights <- setNames(ifelse(cells$group == "B", .2, -.2), cells$cell)
  out <- dif_contrasts(fit, items = "I1", contrasts = list(group = weights))
  resolved <- .dif_resolve(fit, "I1", grp, 20)
  w <- weights[resolved$levs]
  expect_equal(out$table$estimate, sum(w * resolved$loc))
  expect_equal(out$table$se, sqrt(drop(t(w) %*% resolved$vloc %*% w)))
  expect_true(is.finite(out$table$p))
})

test_that("an unweighted middle level is not required by a linear trend", {
  # contr.poly leaves the middle level of an odd-K linear contrast at about
  # -1e-17. That residue must not be read as a cell the trend needs.
  three <- factor(c("1", "2", "3"), ordered = TRUE)
  fam <- .dif_factor_contrasts(three, "band")
  expect_identical(unname(fam[["band: linear"]][["2"]]), 0)
  expect_identical(unname(.dif_leading(three)$weights[["2"]]), 0)

  set.seed(814)
  n <- c(90L, 10L, 100L)
  theta <- rnorm(sum(n))
  X <- matrix(rbinom(length(theta) * 5L, 1,
                     plogis(outer(theta, seq(-1, 1, length.out = 5), "-"))),
              length(theta), 5L, dimnames = list(NULL, paste0("I", 1:5)))
  band <- factor(rep(c("1", "2", "3"), n), ordered = TRUE)
  fit <- rasch(X, factors = data.frame(band = band))
  out <- dif_contrasts(fit, factors = "band", items = c("I1", "I2"),
                       min_n = 20)
  tab <- out$table
  lin <- tab$contrast == "band: linear"
  expect_equal(sum(lin), 2L)
  expect_true(all(is.finite(tab$estimate[lin])))
  expect_true(all(is.finite(tab$p[lin])))
  # the quadratic contrast does weight the thin middle level, so it stays
  # withheld and explained
  expect_true(all(is.na(tab$estimate[!lin])))
  expect_match(paste(out$notes, collapse = " "),
               "quadratic\\]: estimate and inference withheld")
})
