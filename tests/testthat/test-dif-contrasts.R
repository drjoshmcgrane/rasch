# Planned DIF contrasts: the auto-derived family, logit-scale estimation on
# resolved locations, and person-level residual tests for stacked
# repeated-measures designs (Maxwell & Delaney 2004; Hagquist & Andrich 2017).

test_that("the auto family follows the factor structure and finds planted DIF", {
  set.seed(42); n <- 800
  gender <- rep(c("m", "f"), n / 2)
  age <- sample(c("20", "30", "40", "50"), n, replace = TRUE)
  d <- seq(-1.5, 1.5, length.out = 8)
  th <- rnorm(n)
  shift <- matrix(0, n, 8)
  shift[gender == "f", 3] <- 0.8
  shift[, 5] <- (as.numeric(age) - 35) / 10 * 0.4
  X <- matrix(rbinom(n * 8, 1, plogis(outer(th, d, "-") - shift)), n, 8)
  colnames(X) <- paste0("I", 1:8)
  fit <- rasch(data.frame(X, gender = gender, age = age),
               factors = c("gender", "age"))
  dc <- dif_contrasts(fit, items = c("I3", "I5", "I7"))

  # the derived family: difference, linear + quadratic trends, interaction
  expect_setequal(dc$family$contrast,
                  c("gender: m - f", "age: linear", "age: quadratic",
                    "gender(m - f) x age(linear)"))
  t <- dc$table
  # planted uniform gender DIF on I3, in the planted direction
  g3 <- t[t$item == "I3" & t$contrast == "gender: m - f", ]
  expect_true(g3$significant && g3$practical && g3$estimate < -0.5)
  # planted linear age drift on I5
  a5 <- t[t$item == "I5" & t$contrast == "age: linear", ]
  expect_true(a5$significant && a5$estimate > 0.5)
  # the clean item stays clean
  expect_false(any(t$significant[t$item == "I7"]))
  # every contrast reads as a difference of weighted averages (cells strings
  # are rounded to 2 dp for display, hence the loose tolerance)
  expect_true(all(abs(vapply(strsplit(dc$family$cells, ", "), function(cc)
    sum(abs(as.numeric(sub("^.* ", "", cc)))), 0) - 2) < 0.05))
  # The independent-row calibration uses the limiting normal reference.
  expect_equal(t$statistic, t$estimate / t$se, tolerance = 1e-10)
  expect_true(all(is.infinite(t$df)))
  expect_equal(t$p, 2 * stats::pt(-abs(t$statistic), df = t$df))
  expect_equal(t$lower,
               t$estimate - stats::qt(0.975, df = t$df) * t$se)
})

test_that("planned DIF contrasts withhold Wald inference from invalid covariance", {
  set.seed(142)
  n <- 500
  X <- matrix(rbinom(n * 6, 1, .5), n, 6,
              dimnames = list(NULL, paste0("I", 1:6)))
  grp <- factor(rep(c("a", "b"), each = n / 2))
  fit <- rasch(data.frame(X, grp = grp), factors = "grp")
  bad_refit <- split_items(fit, "I2", by = fit$factors$grp)
  split_rows <- grep("^I2 \\(", bad_refit$items$item)
  split_ids <- bad_refit$thresholds$id[
    bad_refit$thresholds$item %in% split_rows]
  bad_refit$est$cov_tau[split_ids[1], split_ids[1]] <- -1e6

  guarded <- testthat::with_mocked_bindings(
    dif_contrasts(fit, items = "I2"),
    .rasch_refit = function(...) bad_refit,
    .package = "rasch")
  expect_true(all(is.finite(guarded$table$estimate)))
  expect_true(all(is.na(guarded$table$se)))
  expect_true(all(is.na(guarded$table$p)))
  expect_false(any(guarded$table$significant))
  expect_true(any(grepl("not positive semidefinite", guarded$notes)))

  supported_refit <- split_items(fit, "I2", by = fit$factors$grp)
  supported_refit$est$cluster_support <- list(
    repeated = FALSE, n = 16L, effective = 16)
  supported_refit$est$cluster_inference <- FALSE
  unsupported <- testthat::with_mocked_bindings(
    dif_contrasts(fit, items = "I2"),
    .rasch_refit = function(...) supported_refit,
    .package = "rasch")
  expect_true(all(is.finite(unsupported$table$estimate)))
  expect_true(all(is.na(unsupported$table$se)))
  expect_true(all(is.na(unsupported$table$df)))
  expect_true(all(is.na(unsupported$table$p)))
  expect_match(paste(unsupported$notes, collapse = " "),
               "independent-person support")
})

test_that("stacked designs use person-level scores and detect drift over time", {
  skip_on_cran()
  set.seed(9); n <- 400
  d <- seq(-1.5, 1.5, length.out = 8)
  th <- rnorm(n)
  gender <- rep(c("m", "f"), n / 2)
  gen <- function(shift4) {
    s <- matrix(0, n, 8); s[, 4] <- shift4
    matrix(rbinom(n * 8, 1, plogis(outer(th, d, "-") - s)), n, 8)
  }
  X <- rbind(gen(0), gen(0.6)); colnames(X) <- paste0("I", 1:8)
  dat <- data.frame(X, time = rep(c("1", "2"), each = n),
                    gender = rep(gender, 2))
  id <- rep(sprintf("P%03d", 1:n), 2)
  fit <- rasch(dat, factors = c("time", "gender"), id = id)
  dc <- dif_contrasts(fit, items = c("I4", "I6"), id = id)

  # time varies within id, so it is detected as within-subject
  expect_identical(dc$within, "time")
  expect_true(dc$paired)
  t <- dc$table
  w4 <- t[t$item == "I4" & t$contrast == "time: 2 - 1", ]
  # The time contrast is paired within person and marginalised over gender;
  # its nuisance-cell means use a Welch--Satterthwaite reference.
  expect_gt(w4$df, n - 10)
  expect_lt(w4$df, n)
  expect_true(w4$significant && w4$estimate > 0.4 && w4$statistic > 0)
  expect_true(all(is.na(t$se)))
  expect_true(all(is.na(t$lower)) && all(is.na(t$upper)))
  # clean item, null gender effect, null interaction all stay quiet
  expect_false(any(t$significant[t$item == "I6"]))
  expect_false(t$significant[t$item == "I4" & t$contrast == "gender: m - f"])
  expect_false(t$significant[t$item == "I4" &
                             t$contrast == "time(2 - 1) x gender(m - f)"])
  # The resolved magnitude uses the person-clustered calibration covariance.
  ds <- dif_size(fit, "I4", by = "time")
  expect_true(any(abs(ds$pairs$difference) > 0.3))
  expect_true(all(is.finite(ds$levels$se)))
  expect_true(all(is.finite(ds$pairs$se)))
  expect_type(ds$pairs$significant, "logical")
  # the fitted identifier is used automatically when it is not repeated in
  # the call
  auto <- dif_contrasts(fit, items = "I4", within = "time")
  expect_true(auto$paired)
  expect_identical(auto$within, "time")
  expect_equal(auto$table$df[auto$table$contrast == "time: 2 - 1"], w4$df)
  expect_error(dif_contrasts(fit, items = "I4", within = "time",
                             id = seq_len(10)), "one value per")
})

test_that("planned contrasts validate the declared within-person structure", {
  set.seed(901)
  n <- 80L
  X <- matrix(rbinom(n * 5L, 1L, 0.5), n, 5L,
              dimnames = list(NULL, paste0("I", 1:5)))
  between <- factor(rep(c("A", "B"), each = n / 2L))
  ordinary <- rasch(data.frame(X, group = between), factors = "group")
  expect_error(
    dif_contrasts(ordinary, items = "I1", within = "group"),
    "need repeated person ids")
  expect_error(
    dif_posthoc(ordinary, "I1", "group", within = "group"),
    "need repeated person ids")

  id <- rep(sprintf("P%03d", seq_len(n)), 2L)
  occasion <- factor(rep(c("T1", "T2"), each = n))
  group <- rep(between, 2L)
  stacked <- rasch(
    data.frame(rbind(X, X), occasion = occasion, group = group),
    id = id, factors = c("occasion", "group"))
  expect_error(
    dif_contrasts(stacked, items = "I1", within = "group"),
    "declared within-subject never vary")
  expect_error(
    dif_contrasts(stacked, items = "I1", within = character(0)),
    "vary within persons but are not declared")
})

test_that("within-person follow-ups marginalise nuisance cells consistently", {
  set.seed(91)
  n_a <- 100L; n_b <- 900L; n <- n_a + n_b
  h <- factor(c(rep("a", n_a), rep("b", n_b)))
  time <- rep(c("t1", "t2"), each = n)
  id <- rep(sprintf("P%04d", seq_len(n)), 2)
  theta <- rnorm(n)
  difficulty <- seq(-1, 1, length.out = 5)
  make_wave <- function()
    matrix(rbinom(n * 5, 1, plogis(outer(theta, difficulty, "-"))), n, 5)
  X <- rbind(make_wave(), make_wave())
  colnames(X) <- paste0("I", seq_len(ncol(X)))
  dat <- data.frame(X, time = time, h = rep(h, 2))
  fit <- rasch(dat, factors = c("time", "h"), id = id)

  # The equally weighted time contrast is +1 in nuisance cell a and -1 in b,
  # hence zero after marginalising equally over the two cells. A shortcut that
  # averages people instead would target 0.1(+1) + 0.9(-1) = -0.8.
  eps <- ave(rnorm(n, 0, 0.2), h, FUN = function(x) x - mean(x))
  delta <- ifelse(h == "a", 1, -1) + eps
  fit$residuals[, 1] <- c(-delta / 2, delta / 2)

  dc <- dif_contrasts(fit, items = "I1", within = "time")
  row <- dc$table[dc$table$contrast == "time: t2 - t1", ]
  expect_equal(row$statistic, 0, tolerance = 1e-10)
  expect_equal(row$p, 1, tolerance = 1e-10)

  # dif_posthoc() uses the same marginal comparison.
  ph <- dif_posthoc(fit, "I1", "time", within = "time")
  expect_equal(ph$table$statistic, row$statistic, tolerance = 1e-10)
  expect_equal(ph$table$p, row$p, tolerance = 1e-10)
})

test_that("mixed follow-ups retain equal margins over an imbalanced nuisance factor", {
  set.seed(92)
  n_s1 <- 120L; n_s2 <- 1080L; n <- n_s1 + n_s2
  site <- factor(c(rep("s1", n_s1), rep("s2", n_s2)))
  group <- factor(unlist(lapply(c(n_s1, n_s2), function(nn)
    rep(c("A", "B"), each = nn / 2))))
  id <- rep(sprintf("P%03d", seq_len(n)), 2)
  time <- factor(rep(c("t1", "t2"), each = n))
  theta <- rnorm(n)
  difficulty <- seq(-1, 1, length.out = 4)
  make_wave <- function()
    matrix(rbinom(n * 4, 1, plogis(outer(theta, difficulty, "-"))), n, 4)
  X <- rbind(make_wave(), make_wave())
  colnames(X) <- paste0("I", seq_len(ncol(X)))
  dat <- data.frame(X, time = time, group = rep(group, 2),
                    site = rep(site, 2))
  fit <- rasch(dat, factors = c("time", "group", "site"), id = id)

  # The time-by-group contrast is +1 at site s1 and -1 at site s2. Its
  # equal-site marginal value is therefore zero, although a person-frequency
  # shortcut would target 0.1(+1) + 0.9(-1) = -0.8.
  eps <- ave(rnorm(n, 0, 0.2), interaction(site, group),
             FUN = function(x) x - mean(x))
  site_effect <- ifelse(site == "s1", 1, -1)
  delta <- eps + ifelse(group == "B", site_effect, 0)
  fit$residuals[, 1] <- c(-delta / 2, delta / 2)

  ph <- dif_posthoc(fit, "I1", term = c("time", "group"),
                    within = "time")
  expect_equal(nrow(ph$table), 1L)
  expect_equal(ph$table$statistic, 0, tolerance = 1e-10)
  expect_equal(ph$table$p, 1, tolerance = 1e-10)
})

test_that("custom cell-weight contrasts are accepted and normalised", {
  set.seed(5); n <- 500
  g <- sample(c("a", "b", "c"), n, replace = TRUE)
  d <- seq(-1, 1, length.out = 6)
  sh <- matrix(0, n, 6); sh[g == "c", 2] <- 0.9   # DIF on I2 only
  X <- matrix(rbinom(n * 6, 1, plogis(outer(rnorm(n), d, "-") - sh)), n, 6)
  colnames(X) <- paste0("I", 1:6)
  fit <- rasch(data.frame(X, g = g), factors = "g")
  dc <- dif_contrasts(fit, items = "I2",
                      contrasts = list("c vs rest" = c(a = -1, b = -1, c = 2)))
  r <- dc$table
  expect_equal(nrow(r), 1L)
  expect_true(r$significant && r$estimate > 0.5)
  expect_error(dif_contrasts(fit, items = "I2",
                             contrasts = list(bad = c(x = 1, y = -1))),
               "design cells")
  shaped <- matrix(c(-1, -1, 2), 1L)
  names(shaped) <- c("a", "b", "c")
  expect_error(dif_contrasts(fit, items = "I2",
                             contrasts = list(bad = shaped)),
               "plain vector")
})

test_that("DIF post-hocs give marginal pairs and pure interaction magnitudes", {
  set.seed(31); n <- 1800
  g1 <- factor(rep(c("a", "b", "c"), each = n / 3))
  g2 <- factor(rep(rep(c("x", "y", "z"), each = n / 9), 3))
  d <- seq(-1.3, 1.3, length.out = 6)
  shift <- ifelse(g1 == "c" & g2 == "z", 1.0, 0)
  X <- matrix(rbinom(n * 6, 1,
    plogis(outer(rnorm(n), d, "-") - outer(shift, c(0, 1, 0, 0, 0, 0)))),
    n, 6)
  colnames(X) <- paste0("I", 1:6)
  fit <- rasch(data.frame(X, g1 = g1, g2 = g2),
               factors = c("g1", "g2"))

  main <- dif_posthoc(fit, "I2", "g1")
  expect_s3_class(main, "rasch_dif_posthoc")
  expect_equal(nrow(main$table), choose(3, 2))
  expect_true(all(main$table$p_adj >= main$table$p - 1e-12))

  intr <- dif_posthoc(fit, "I2", "g1:g2")
  expect_equal(nrow(intr$table), choose(3, 2)^2)
  expect_identical(intr$type, "interaction magnitude")

  # The c-a by z-x row is the difference-in-differences of the jointly
  # resolved cell locations, not merely the largest pair of cell means.
  cells <- dif_size(fit, "I2", by = c("g1", "g2"))$levels
  loc <- setNames(cells$location, cells$level)
  manual <- (loc[["c:z"]] - loc[["a:z"]]) -
    (loc[["c:x"]] - loc[["a:x"]])
  row <- intr$table[intr$table$contrast == "c - a x z - x", ]
  expect_equal(row$estimate, manual, tolerance = 1e-8)
  expect_gt(abs(row$estimate), 0.6)
})

test_that("an item whose split refit is refused withholds only its own rows", {
  # In group b, item I3 scores 0 only inside all-zero patterns, so its
  # split copy loses that category during calibration and split_items
  # refuses. The refusal must not abort the other items' contrasts.
  set.seed(11); n <- 200; L <- 6
  d <- seq(-1.5, 1.5, length.out = L); th <- rnorm(n)
  X <- matrix(rbinom(n * L, 1, plogis(outer(th, d, "-"))), n, L)
  colnames(X) <- paste0("I", seq_len(L))
  g <- rep(c("a", "b"), each = n / 2); b <- g == "b"
  X[b, "I3"] <- 1L
  X[which(b)[1:6], ] <- 0L
  fit <- rasch(data.frame(X, grp = g, check.names = FALSE), factors = "grp")
  expect_error(split_items(fit, "I3", by = "grp"),
               "cannot preserve the fitted score structure")

  dc <- dif_contrasts(fit)
  tab <- dc$table
  expect_equal(nrow(tab), L)
  refused <- tab$item == "I3"
  expect_true(all(is.na(unlist(tab[refused,
    c("estimate", "se", "statistic", "p", "p_adj")]))))
  expect_true(all(is.finite(tab$estimate[!refused])))
  expect_true(all(is.finite(tab$p[!refused])))
  expect_match(paste(dc$notes, collapse = " "),
               "I3: resolved contrasts withheld because the split refit is")
  # the withheld item keeps its place in the multiplicity family
  expect_equal(dc$family_n, L)
  expect_equal(tab$p_adj[!refused],
               p.adjust(tab$p[!refused], "holm", n = L))
})

test_that("a refused split refit is not reported as a category mismatch", {
  # An anchored item cannot be split, but both groups answer it on the same
  # 0:1 structure. The withholding note must give the refusal's own reason
  # and must not also blame the response categories.
  set.seed(7); n <- 300L; L <- 6L
  d <- seq(-1.5, 1.5, length.out = L); th <- rnorm(n)
  X <- matrix(rbinom(n * L, 1, plogis(outer(th, d, "-"))), n, L)
  colnames(X) <- c("A1", "A2", paste0("I", 3:L))
  df <- data.frame(X, grp = rep(c("a", "b"), each = n / 2),
                   check.names = FALSE)
  f0 <- rasch(df, factors = "grp")
  th0 <- f0$thresholds
  th0$item <- f0$items$item[th0$item]
  fit <- rasch(df, factors = "grp",
               anchors = th0[th0$item %in% c("A1", "A2"), c("item", "k", "tau")])

  dc <- dif_contrasts(fit)
  notes <- paste(dc$notes, collapse = " ")
  expect_match(notes, "A1: resolved contrasts withheld because the split refit")
  expect_match(notes, "an anchored item cannot be split")
  expect_false(grepl("response-category", notes))
  anchored <- dc$table$item %in% c("A1", "A2")
  expect_true(all(is.na(dc$table$estimate[anchored])))
  expect_true(all(is.finite(dc$table$estimate[!anchored])))
})
