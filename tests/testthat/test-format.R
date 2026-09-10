
test_that("results tables print without scientific notation", {
  # base R prints a p of 4e-83 as "4.00e-83"; a large sample produces
  # probabilities that small routinely, and they are unreadable that way
  d <- data.frame(item = c("I1", "I2"), n = c(10L, 20L),
                  location = c(1.23456, -0.5), p_adj = c(4e-83, 0.42),
                  flagged = c(TRUE, FALSE))
  out <- capture.output(print(.tag_tables(d)))
  expect_false(any(grepl("e[-+][0-9]", out)))
  expect_true(any(grepl("< 0.001", out, fixed = TRUE)))
  expect_true(any(grepl("1.235", out, fixed = TRUE)))   # fixed decimals
  expect_true(any(grepl("\\b10\\b", out)))              # integers stay integers

  # tagging keeps the frame a data frame and survives subsetting
  td <- .tag_tables(d)
  expect_s3_class(td, "data.frame")
  expect_true(is.data.frame(td))
  expect_s3_class(td[, c("item", "p_adj")], "rasch_table")
  expect_identical(as.data.frame(td)$p_adj, d$p_adj)   # values untouched

  # a table that already carries a class of its own keeps its own printing
  keep <- d
  class(keep) <- c("something_else", "data.frame")
  expect_false(inherits(.tag_tables(keep), "rasch_table"))
})

test_that("table printing validates display controls", {
  x <- .tag_tables(data.frame(p = 0.25))
  expect_error(print(x, digits = c(2, 3)), "digits")
  expect_error(print(x, digits = -1), "digits")
  expect_error(print(x, digits = matrix(2)), "digits")
  expect_error(print(x, n = NA_real_), "`n`")
  expect_error(print(x, n = 1.5), "`n`")
  expect_error(print(x, n = matrix(2)), "`n`")
  expect_output(print(x, n = Inf), "0.250")
  expect_output(print(x, n = 0), "more row")
})

test_that("observed and expected proportions are not read as probabilities", {
  # obs_p and est_p are category proportions: a category nobody chose has an
  # observed proportion of exactly zero, which "< 0.001" would misreport
  expect_false(.is_pcol("obs_p"))
  expect_false(.is_pcol("est_p"))
  expect_true(.is_pcol("p"))
  expect_true(.is_pcol("p_adj"))
  expect_true(.is_pcol("total_p"))
  # an empty category beside a populated one: the zero must read as zero
  expect_equal(.fmt_df(data.frame(obs_p = c(0, 0.2273)))$obs_p,
               c("0.000", "0.227"))
  expect_equal(.fmt_df(data.frame(p_adj = c(0, 0.2273)))$p_adj,
               c("< 0.001", "0.227"))
})

test_that("a fitted object's tables print through the shared formatter", {
  set.seed(1)
  N <- 1200; K <- 6
  th <- stats::rnorm(N); dl <- seq(-1.2, 1.2, length.out = K)
  X <- vapply(seq_len(K), function(i)
    stats::rbinom(N, 1, stats::plogis(th - dl[i])), numeric(N))
  colnames(X) <- sprintf("I%02d", seq_len(K))
  f <- rasch(data.frame(id = seq_len(N), X), id = "id")
  expect_s3_class(f$items, "rasch_table")
  expect_s3_class(f$item_trait, "rasch_table")
  expect_false(any(grepl("e[-+][0-9]", capture.output(print(f$items)))))
  expect_false(any(grepl("e[-+][0-9]", capture.output(print(f)))))
})

test_that("a saved table carries plain decimals rather than exponents", {
  p <- tempfile(fileext = ".csv")
  on.exit(unlink(p))
  .write_csv_plain(data.frame(item = "I1", p_adj = 2.3e-17, loc = 1.2345678), p)
  txt <- readLines(p)
  expect_false(any(grepl("e[-+][0-9]", txt)))
  # and the precision is kept, not rounded away
  expect_equal(utils::read.csv(p)$p_adj, 2.3e-17)
  expect_equal(utils::read.csv(p)$loc, 1.2345678)
})

test_that("a deliberate refusal is signalled apart from a fault", {
  skip_on_cran()
  # the application showed considered refusals in the same red as a crash,
  # because a caller had no way to tell them apart
  d <- simulate_efrm(n_per_group = 200, items_per_set = 5, n_sets = 2,
                     n_groups = 2, seed = 4)
  tr <- attr(d, "truth")
  f <- rasch_efrm(d, item_sets = tr$item_sets, groups = "group", id = "id",
                  boot_reps = 0)

  # residual components are undefined on structurally disjoint columns
  expect_error(residual_pca(f), class = "rasch_refusal")
  # a frame-defining factor has no within-frame contrast left to test
  expect_error(dif_anova(f, "group"), class = "rasch_refusal")
  # and the message still reads as before, so existing expectations hold
  expect_error(dif_anova(f, "group"), "define\\(s\\) the EFRM frame structure")

  # a genuine fault must NOT be dressed as a refusal
  e <- tryCatch(residual_pca("nonsense"), error = function(e) e)
  expect_false(inherits(e, "rasch_refusal"))
  expect_s3_class(e, "error")
})

test_that("the ETS categories follow the published rule on the logit scale", {
  # 2.35 delta units to the logit, so A tops out at 0.426 and C starts at
  # 0.638; C additionally tests against the A ceiling rather than against zero
  a_cut <- 1 / 2.35; c_cut <- 1.5 / 2.35
  small <- 0.001                      # a standard error tight enough to decide

  expect_equal(.ets_category(a_cut - 0.01, small, 1e-10), "A")
  expect_equal(.ets_category(a_cut + 0.01, small, 1e-10), "B+")
  expect_equal(.ets_category(c_cut + 0.01, small, 1e-10), "C+")
  expect_equal(.ets_category(-(c_cut + 0.01), small, 1e-10), "C-")

  # not significantly different from zero is A whatever the magnitude
  expect_equal(.ets_category(1.2, 0.9, 0.4), "A")
  # large but not significantly beyond the A ceiling is B, not C
  expect_equal(.ets_category(c_cut + 0.01, 0.5, 0.01), "B+")
  expect_true(is.na(.ets_category(NA_real_, NA_real_, NA_real_)))
  expect_true(is.na(.ets_category(Inf, small, 0.01)))
  expect_true(is.na(.ets_category(1, Inf, 0.01)))
  expect_true(is.na(.ets_p_beyond(1, Inf)))
})

test_that("ETS classifications do not replace unavailable tests with A or B", {
  expect_true(all(is.na(.ets_category(c(.2, .8, -.8), .1,
                                      c(NA_real_, NaN, Inf)))))
  expect_true(all(is.na(.ets_category(.8, c(0, -1), .01))))
  expect_true(all(is.na(.ets_category(c(.8, -.8), .1, .01,
                                      p_beyond = NA_real_))))
  # The interval-null test cannot affect these already determined classes.
  expect_identical(.ets_category(c(.2, .5, .8), .1, c(.01, .01, .2),
                                p_beyond = NA_real_), c("A", "B+", "A"))

  d <- simulate_rasch(300, 6, n_groups = 2, seed = 922817)
  fit <- rasch(d, id = "id", factors = "group")
  original <- split_items
  testthat::local_mocked_bindings(split_items = function(...) {
    out <- original(...)
    # A degenerate resolved covariance must not imply negligible DIF.
    out$est$cov_tau[,] <- 0
    out
  }, .package = "rasch")
  size <- dif_size(fit, fit$items$item[1L], "group")$pairs
  expect_true(all(is.finite(size$difference)))
  expect_true(all(size$se == 0))
  expect_true(all(is.na(size$p_adj)))
  expect_true(all(is.na(size$ets)))
})

test_that("dif_size reports an ETS category beside the magnitude", {
  set.seed(2)
  N <- 2000; K <- 8
  th <- stats::rnorm(N); dl <- seq(-1.2, 1.2, length.out = K)
  g <- rep(c("ref", "foc"), each = N / 2)
  X <- vapply(seq_len(K), function(i)
    stats::rbinom(N, 1, stats::plogis(th - dl[i])), numeric(N))
  X[g == "foc", 3] <- stats::rbinom(sum(g == "foc"), 1,
    stats::plogis(th[g == "foc"] - dl[3] - 0.9))
  colnames(X) <- sprintf("I%02d", seq_len(K))
  f <- rasch(data.frame(id = seq_len(N), X, grp = g), id = "id",
             factors = "grp")

  big <- dif_size(f, "I03", "grp")$pairs
  expect_true("ets" %in% names(big))
  expect_match(big$ets, "^C")               # 0.9 logits is 2.1 delta
  clean <- dif_size(f, "I06", "grp")$pairs
  expect_equal(clean$ets, "A")
})

test_that("a polytomous DIF comparison reports signed area descriptively", {
  skip_on_cran()
  ta <- c(-1, 0, 1) + 0.4
  tb <- c(-1, 0, 1) - 0.2
  numeric_area <- stats::integrate(function(th)
    vapply(th, function(z) item_moments(z, tb)$E -
                     item_moments(z, ta)$E, 0), -20, 20)$value
  expect_equal(numeric_area, sum(ta - tb), tolerance = 1e-6)

  set.seed(5)
  N <- 2000; K <- 6
  th <- stats::rnorm(N, 0, 1.3)
  tau <- lapply(seq_len(K), function(i) c(-1, 0, 1) +
                  seq(-0.5, 0.5, length.out = K)[i])
  g <- rep(c("ref", "foc"), each = N / 2)
  sh <- ifelse(g == "foc", 0.9, 0)
  X <- vapply(seq_len(K), function(i)
    vapply(seq_len(N), function(p) {
      s <- if (i == 2L) sh[p] else 0
      sample(0:3, 1, prob = item_moments(th[p] - s, tau[[i]])$P) }, 0),
    numeric(N))
  colnames(X) <- sprintf("I%02d", seq_len(K))
  f <- rasch(data.frame(id = seq_len(N), X, grp = g), id = "id",
             factors = "grp")
  p <- dif_size(f, "I02", "grp")$pairs
  expect_true(all(f$items$max == 3))          # genuinely polytomous
  expect_true(is.na(p$ets))                   # ETS is dichotomous
  expect_equal(p$signed_area, 3 * p$difference)
  expect_identical(dif_size(f, "I02", "grp")$classification,
                   "PCM signed expected-score area (descriptive)")

  # A group-specific empty category changes the fitted threshold structure;
  # neither a PCM location contrast nor its signed area is then comparable.
  X2 <- X
  X2[g == "foc", 2][X2[g == "foc", 2] == 3] <- 2
  f2 <- rasch(data.frame(id = seq_len(N), X2, grp = g), id = "id",
              factors = "grp")
  d2 <- dif_size(f2, "I02", "grp")
  expect_true(is.na(d2$pairs$difference))
  expect_true(is.na(d2$pairs$signed_area))
  expect_match(d2$notes, "different observed response-category")

  # Equal maxima are not enough: different intermediate observed categories
  # would be collapsed differently by the two resolved calibrations.
  X3 <- X
  X3[g == "ref", 2][X3[g == "ref", 2] == 1] <- 0
  X3[g == "foc", 2][X3[g == "foc", 2] == 2] <- 1
  f3 <- rasch(data.frame(id = seq_len(N), X3, grp = g), id = "id",
              factors = "grp")
  d3 <- dif_size(f3, "I02", "grp")
  expect_true(is.na(d3$pairs$difference))
  expect_true(is.na(d3$pairs$signed_area))
})

test_that("crossed factor cells cannot collapse when labels contain separators", {
  x <- data.frame(A = c("a:b", "a"), B = c("c", "b:c"),
                  check.names = FALSE)
  z <- .factor_cells(x, sep = ":")
  expect_equal(nlevels(z), 2L)
  expect_false(z[1] == z[2])
  expect_true(all(grepl("a:b:c", as.character(z), fixed = TRUE)))
})

test_that("factor keys are stable across subsets and punctuation", {
  x <- data.frame(
    A = c("a:b", "a", "a:b", NA),
    B = c("c", "b:c", "c", "N;"),
    stringsAsFactors = FALSE
  )
  k <- .factor_keys(x)
  expect_identical(k[c(1, 3)], .factor_keys(x[c(1, 3), ]))
  expect_false(k[1] == k[2])
  expect_false(k[4] == .factor_keys(data.frame(A = "N;", B = "N;")))
})

test_that("EFRM information designs do not parse set labels", {
  fit <- structure(list(
    tau_list = rep(list(0), 4),
    virtual_map = data.frame(
      item = paste0("I", 1:4), group = rep("G", 4),
      set = c("A+B", "C", "A", "B+C"), stringsAsFactors = FALSE),
    X = matrix(c(1, 0, NA, NA,
                 NA, NA, 1, 0), nrow = 2, byrow = TRUE)
  ), class = c("rasch_efrm", "rasch"))
  z <- .design_blocks(fit)
  expect_length(z, 2L)
  expect_setequal(unname(z), list(1:2, 3:4))
  expect_equal(length(unique(names(z))), 2L)
})
