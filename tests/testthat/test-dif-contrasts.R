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
  # z-statistics equal estimate/se in the independent-rows case
  expect_equal(t$statistic, t$estimate / t$se, tolerance = 1e-10)
})

test_that("stacked designs use person-level scores and detect drift over time", {
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
  fit <- rasch(dat, factors = c("time", "gender"))
  id <- rep(sprintf("P%03d", 1:n), 2)
  dc <- dif_contrasts(fit, items = c("I4", "I6"), id = id)

  # time varies within id, so it is detected as within-subject
  expect_identical(dc$within, "time")
  expect_true(dc$paired)
  t <- dc$table
  w4 <- t[t$item == "I4" & t$contrast == "time: 2 - 1", ]
  # paired t (df = persons - 1), significant, positive drift, sign-aligned
  expect_equal(w4$df, n - 1)
  expect_true(w4$significant && w4$estimate > 0.4 && w4$statistic > 0)
  # clean item, null gender effect, null interaction all stay quiet
  expect_false(any(t$significant[t$item == "I6"]))
  expect_false(t$significant[t$item == "I4" & t$contrast == "gender: m - f"])
  expect_false(t$significant[t$item == "I4" &
                             t$contrast == "time(2 - 1) x gender(m - f)"])
  # within requires id
  expect_error(dif_contrasts(fit, items = "I4", within = "time"), "id")
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
})
