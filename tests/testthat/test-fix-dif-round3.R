# Regressions for the round-three DIF fixes: availability of the robust
# uniform test under a rank-deficient nuisance column, and the notes that
# explain every withheld statistic.

# the adjacent-cell design of the report: one item unseen in the (b, y) cell
# leaves a quantile class interval occupied by that cell alone, so the
# site-by-interval column is aliased while every group term stays estimable
sim_aliased_cell <- function() {
  set.seed(12); n <- 800L; L <- 7L
  d <- seq(-1.6, 1.6, length.out = L)
  th <- rnorm(n)
  g <- rep(c("a", "b"), each = n / 2)
  s <- rep(c("x", "y"), length.out = n)
  X <- matrix(rbinom(n * L, 1, plogis(outer(th, d, "-"))), n, L)
  colnames(X) <- paste0("I", seq_len(L))
  X[g == "b" & s == "y", "I1"] <- NA
  rasch(data.frame(X, grp = g, site = s), factors = c("grp", "site"))
}

test_that("an aliased nuisance column does not withhold an estimable term", {
  set.seed(31); n <- 240L
  f1 <- factor(rep(c("a", "b"), each = n / 2))
  f2 <- factor(rep(c("x", "y"), length.out = n))
  ci <- factor(rep(1:3, length.out = n))
  # interval 4 holds (b, y) persons only: f2:ci is aliased there, f1 is not
  by <- which(f1 == "b" & f2 == "y")[1:30]
  ci <- factor(replace(as.character(ci), by, "4"))
  d <- data.frame(z = rnorm(n), f1 = f1, f2 = f2, ci = ci)
  tl <- c("f1", "f2", "ci", "f1:ci", "f2:ci")
  got <- .dif_type2(d, tl, variance = "hc3",
                    robust_terms = c("f1", "f2"))

  # the classical model itself aliases one nuisance column only
  m1 <- stats::lm(z ~ f2 + ci + f2:ci + f1, d)
  expect_equal(names(stats::coef(m1))[is.na(stats::coef(m1))], "f2y:ci4")

  # the robust Wald on the columns lm retained is the value reported
  est <- !is.na(stats::coef(m1))
  X <- stats::model.matrix(m1)[, est, drop = FALSE]
  Xi <- solve(crossprod(X))
  h <- pmin(stats::hatvalues(m1), 1 - 1e-8)
  ae <- stats::residuals(m1) / pmax(1 - h, 1e-8)
  V <- Xi %*% crossprod(X * ae) %*% Xi
  j <- which(colnames(X) == "f1b")
  expect_equal(got$F_value[got$term == "f1"],
               unname(stats::coef(m1)["f1b"]^2 / V[j, j]))
  expect_true(is.finite(got$p[got$term == "f1"]))
  expect_true(is.finite(got$F_value[got$term == "f2"]))
  expect_true(all(is.na(got$unavailable_reason[got$term %in% c("f1", "f2")])))
})

test_that("a term aliased with the retained design is withheld with its reason", {
  set.seed(32); n <- 240L
  f1 <- factor(rep(c("a", "b", "c"), length.out = n))
  ci <- factor(rep(1:3, each = n / 3))
  # every level-c person, and only those, falls in interval 4: the f1
  # contrasts are no longer separable from the class interval
  cc <- which(f1 == "c")
  ci <- factor(replace(as.character(ci), cc, "4"))
  d <- data.frame(z = rnorm(n), f1 = f1, ci = ci)
  got <- .dif_type2(d, c("f1", "ci"), variance = "hc3", robust_terms = "f1")
  row <- got$term == "f1"
  expect_true(any(row))
  expect_true(is.na(got$F_value[row]))
  expect_true("unavailable_reason" %in% names(got))
  expect_match(got$unavailable_reason[row],
               "aliased with the rest of the retained design")
})

test_that("dif_anova tests every complete item when one item aliases a cell", {
  skip_on_cran()
  fit <- sim_aliased_cell()
  da <- dif_anova(fit)
  uni <- da$terms[da$terms$item != "I1" &
                    da$terms$term %in% c("grp", "site"), ]
  expect_equal(nrow(uni), 12L)
  expect_true(all(is.finite(uni$F_value)))
  expect_true(all(is.finite(uni$p)))
  # the item carrying the structural hole is still tested on its own design
  own <- da$terms[da$terms$item == "I1" &
                    da$terms$term %in% c("grp", "site"), ]
  expect_true(all(is.finite(own$F_value)))
  # no uniform hypothesis is left unanswered, so resolve_dif() no longer
  # qualifies its verdict
  r <- resolve_dif(fit, max_splits = 0)
  expect_equal(r$n_untested, 0L)
  expect_false(grepl("not estimable", r$stopped))
})

test_that("the aggregate withheld count matches the per-row notes", {
  skip_on_cran()
  da <- dif_anova(sim_aliased_cell(), effects = "factorial")
  rows <- sum(grepl("unavailable because", da$notes))
  aggregate <- grep("requested item-term test\\(s\\) were not estimable",
                    da$notes, value = TRUE)
  expect_length(aggregate, 1L)
  expect_equal(as.integer(sub("^([0-9]+) .*$", "\\1", aggregate)), rows)
  # the structurally absent interaction is the withheld test, with its reason
  expect_match(paste(da$notes, collapse = " "),
               "I1 [grp:site]: unavailable because the item's own design",
               fixed = TRUE)
})

test_that("p_adjust = 'none' is not described as an adjustment family", {
  skip_on_cran()
  da <- dif_anova(sim_aliased_cell(), effects = "factorial",
                  p_adjust = "none")
  notes <- paste(da$notes, collapse = " ")
  expect_match(notes, "unavailable because")
  expect_false(grepl("none adjustment family", notes, fixed = TRUE))
  expect_match(notes, "carries no multiplicity adjustment")
})

test_that("a category-structure refusal names the item without a stray space", {
  d <- simulate_rasch(600, 6, model = "PCM", n_categories = 4, seed = 703)
  group <- factor(rep(c("A", "B"), each = 300))
  d$I02[group == "A" & d$I02 == 1L] <- 2L
  d$group <- group
  fit <- rasch(d, id = "id", factors = "group",
               items = sprintf("I%02d", 1:6))
  notes <- paste(dif_contrasts(fit, factors = "group", items = "I02")$notes,
                 collapse = " ")
  expect_match(notes,
               "I02: resolved contrasts withheld because groups have different observed response-category structures",
               fixed = TRUE)
  expect_false(grepl("I02 :", notes, fixed = TRUE))
})

test_that("dif_size carries a refused MFRM refit's own reason", {
  skip_on_cran()
  d <- simulate_mfrm(n_persons = 100, n_items = 5, n_raters = 3,
                     n_categories = 2, seed = 31)
  ids <- unique(d$person)
  d$group <- rep(c("A", "B"), length.out = length(ids))[match(d$person, ids)]
  mf <- rasch_mfrm(d, "person", "item", "score", facets = "rater",
                   factors = "group", n_groups = 2)
  real_split <- split_items
  refuse_i3 <- function(fit, items, by, ...) {
    if (any(grepl("^I3", items)))
      stop("refit did not converge for I3 cells", call. = FALSE)
    real_split(fit, items, by = by, ...)
  }
  ds <- testthat::with_mocked_bindings(
    dif_size(mf, "I3", by = "group"),
    split_items = refuse_i3, .package = "rasch")
  notes <- paste(ds$notes, collapse = " ")
  expect_match(notes, "the split refit is unavailable: refit did not converge")
  expect_false(grepl("response-category", notes))
  expect_true(all(is.na(ds$pairs$difference)))
})

test_that("resolve_dif says why a split carries no magnitude", {
  skip_on_cran()
  set.seed(1); n <- 600L
  d <- seq(-2, 2, length.out = 8)
  g <- rep(c("a", "b"), each = n / 2)
  sh <- matrix(0, n, 8); sh[g == "b", 3] <- 1.2
  X <- matrix(rbinom(n * 8, 1, plogis(outer(rnorm(n), d, "-") - sh)), n, 8)
  colnames(X) <- paste0("I", 1:8)
  fit <- rasch(data.frame(X, grp = g), factors = "grp")
  r <- testthat::with_mocked_bindings(
    resolve_dif(fit),
    dif_posthoc = function(...) stop("post-hoc refit failed", call. = FALSE),
    .package = "rasch")
  expect_equal(r$n_splits, 1L)
  expect_true(is.na(r$splits$magnitude))
  expect_match(paste(r$notes, collapse = " "),
               "I3 [grp]: split, but its DIF magnitude is unavailable",
               fixed = TRUE)
  expect_match(paste(r$notes, collapse = " "), "post-hoc refit failed")
})
