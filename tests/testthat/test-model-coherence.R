test_that("BTL-EFRM likelihood matches its reported fitted probabilities", {
  d <- simulate_btl_efrm(
    n_objects_per_set = 5, n_sets = 3, n_panels = 3,
    n_judges_per_panel = 10, reps_within = 8, reps_cross = 4,
    panel_units = c(0.75, 1, 1.35), set_units = c(1, 1.25, 0.8),
    set_origins = c(0, 0.4, -0.3), seed = 472)
  tr <- attr(d, "truth")
  fit <- btl_efrm(
    d, object_a = "object_a", object_b = "object_b", winner = "winner",
    judge = "judge", panels = "panel", object_sets = tr$object_sets,
    se_method = "conditional", min_link = 5)

  y <- fit$comparisons$response
  p <- fit$comparisons$expected
  ll <- sum(y * log(p) + (1 - y) * log1p(-p))
  expect_equal(fit$loglik, ll, tolerance = 1e-8)
  expect_equal(fit$equal_unit$loglik_frames, ll, tolerance = 1e-8)
})

test_that("DIF contrasts refuse incompatible score structures", {
  d <- simulate_rasch(800, 6, model = "PCM", n_categories = 4, seed = 701)
  group <- factor(rep(c("A", "B"), each = 400))
  d$I02[group == "A" & d$I02 == 1L] <- 2L
  d$group <- group
  fit <- rasch(d, id = "id", factors = "group",
               items = sprintf("I%02d", 1:6))

  dc <- dif_contrasts(fit, factors = "group", items = "I02")
  expect_true(all(is.na(dc$table$estimate)))
  expect_true(all(is.na(dc$table$p_adj)))
  expect_match(paste(dc$notes, collapse = " "),
               "different observed response-category structures")
  expect_error(dif_posthoc(fit, "I02", "group"),
               "different observed response-category structures")
  expect_true(all(is.na(dif_size(fit, "I02", "group")$pairs$difference)))
})

test_that("pooled MFRM contrasts refuse incompatible score structures", {
  d <- simulate_mfrm(120, 4, 3, n_categories = 3, seed = 702)
  group <- setNames(rep(c("A", "B"), each = 60), unique(d$person))
  d$group <- group[d$person]
  d$score[d$item == "I2" & d$group == "A" & d$score == 1L] <- 2L
  fit <- rasch_mfrm(d, person = "person", item = "item", score = "score",
                    facets = "rater", factors = "group")

  dc <- dif_contrasts(fit, factors = "group", items = "I2")
  expect_true(all(is.na(dc$table$estimate)))
  expect_true(all(is.na(dc$table$p_adj)))
  expect_match(paste(dc$notes, collapse = " "),
               "same observed response-category structure")
})

test_that("multifactor DIF magnitudes match the jointly adjusted estimand", {
  set.seed(703)
  cell <- rep(c("AX", "AY", "BX", "BY"), c(600, 200, 200, 600))
  A <- factor(substr(cell, 1, 1)); B <- factor(substr(cell, 2, 2))
  n <- length(cell); theta <- rnorm(n)
  difficulty <- seq(-1.5, 1.5, length.out = 8)
  shift <- (A == "B") * 1 + (B == "Y") * 2
  X <- sapply(seq_along(difficulty), function(j)
    rbinom(n, 1, plogis(theta - difficulty[j] -
                          if (j == 3L) shift else 0)))
  colnames(X) <- paste0("I", seq_len(ncol(X)))
  fit <- rasch(data.frame(X, A = A, B = B), factors = c("A", "B"))
  out <- dif_anova(fit, sizes = TRUE)

  a <- out$sizes[out$sizes$item == "I3" & out$sizes$term == "A", ]
  b <- out$sizes[out$sizes$item == "I3" & out$sizes$term == "B", ]
  expect_equal(out$sizes, out$posthoc)
  expect_lt(abs(abs(a$estimate) - 1), 0.35)
  expect_lt(abs(abs(b$estimate) - 2), 0.35)
})

test_that("BTL DIF magnitudes marginalise all jointly fitted factors", {
  set.seed(704)
  counts <- c(AX = 24, AY = 8, BX = 8, BY = 24)
  cells <- rep(names(counts), counts)
  judges <- sprintf("J%02d", seq_along(cells))
  A <- setNames(substr(cells, 1, 1), judges)
  B <- setNames(substr(cells, 2, 2), judges)
  objects <- paste0("O", 1:6)
  beta <- setNames(seq(-1, 1, length.out = 6), objects)
  pairs <- t(utils::combn(objects, 2))
  rows <- lapply(judges, function(j) {
    z <- data.frame(a = rep(pairs[, 1], each = 2),
                    b = rep(pairs[, 2], each = 2), judge = j)
    shift <- (A[j] == "B") * 1 + (B[j] == "Y") * 2
    eta <- beta[z$a] - beta[z$b] +
      shift * ((z$a == "O3") - (z$b == "O3"))
    z$winner <- ifelse(runif(nrow(z)) < plogis(eta), z$a, z$b)
    z
  })
  fit <- btl(do.call(rbind, rows), "a", "b", "winner", judge = "judge")
  out <- btl_dif(fit, list(A = A, B = B), objects = "O3")

  a <- out$sizes[out$sizes$term == "A", ]
  b <- out$sizes[out$sizes$term == "B", ]
  expect_match(a$level_a, "B:X")
  expect_match(a$level_a, "B:Y")
  expect_lt(abs(abs(a$difference) - 1), 0.35)
  expect_lt(abs(abs(b$difference) - 2), 0.35)
})

test_that("structural refits never renumber the fitted scores", {
  set.seed(705)
  X <- matrix(rbinom(1000, 1, 0.5), 200, 5,
              dimnames = list(NULL, paste0("I", 1:5)))
  X[, 2] <- X[, 1]
  fit <- rasch(X)
  expect_error(combine_items(fit, c("I1", "I2")),
               "cannot preserve the fitted score structure")

  d <- simulate_rasch(600, 6, model = "PCM", n_categories = 3, seed = 706)
  group <- factor(rep(c("A", "B"), each = 300))
  d$I02[group == "A" & d$I02 == 1L] <- 2L
  d$group <- group
  fit2 <- rasch(d, id = "id", factors = "group",
                items = sprintf("I%02d", 1:6))
  expect_error(split_items(fit2, "I02", "group"),
               "without changing its score structure")
})

test_that("resolved EFRM items retain their source score structure", {
  d <- simulate_efrm(n_per_group = 100, items_per_set = 4, n_sets = 1,
                     n_groups = 2, n_categories = 4, seed = 707)
  tr <- attr(d, "truth")
  item <- tr$item_sets[[1]][1]
  group <- sort(unique(d$group))[1]
  d[[item]][d$group == group & d[[item]] == 1L] <- 2L
  fit <- rasch_efrm(d, item_sets = tr$item_sets, groups = "group",
                    id = "id", boot_reps = 0)
  expect_error(resolve_frames(fit, item, boot_reps = 0),
               "cannot preserve the fitted score structure")
})

test_that("ETS decisions use adjusted probabilities throughout", {
  difference <- c(0.8, 0.8, 0.8)
  se <- rep(0.1, 3)
  p_raw <- c(0.02, 0.03, 0.04)
  p_adj <- p.adjust(p_raw, method = "holm")
  beyond_adj <- p.adjust(.ets_p_beyond(difference, se), method = "holm")
  expect_identical(.ets_category(difference, se, p_adj, 0.05, beyond_adj),
                   rep("A", 3))

  # A cluster-based magnitude uses the same finite reference for the ETS
  # interval-null test as for its ordinary Wald probability.
  a_cut <- 1 / ETS_DELTA_PER_LOGIT
  finite_ref <- .ets_p_beyond(0.8, 0.3, df = 9)
  expect_equal(finite_ref,
               stats::pt((-a_cut - 0.8) / 0.3, 9) +
                 stats::pt((a_cut - 0.8) / 0.3, 9))
  expect_gt(finite_ref, .ets_p_beyond(0.8, 0.3))
})

test_that("Rasch anchors name their required columns at the boundary", {
  X <- matrix(rep(c(0L, 1L), 250), 100, 5,
              dimnames = list(NULL, paste0("I", 1:5)))
  expect_error(rasch(X, anchors = data.frame(foo = 1)),
               "missing: item, k, tau")
  expect_error(rasch(X, anchors = data.frame(item = character(),
                                              k = integer(), tau = numeric())),
               "at least one")
  expect_error(pcml(X, anchors = data.frame(item = character(),
                                             k = integer(), tau = numeric())),
               "at least one")
})
