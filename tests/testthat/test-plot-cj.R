# plot_cj(): the joint calibration map

cj_plot_sim <- function(seed = 1, shift = 0, rankings = FALSE) {
  set.seed(seed)
  delta <- seq(-1.5, 1.5, length.out = 8)
  names(delta) <- sprintf("I%02d", 1:8)
  theta <- stats::rnorm(300)
  X <- sapply(delta, function(d) as.integer(stats::runif(300) < stats::plogis(theta - d)))
  d_j <- delta; d_j["I03"] <- d_j["I03"] + shift
  pairs <- t(utils::combn(names(delta), 2))[sample(28, 300, replace = TRUE), ]
  p_a <- stats::plogis(0.6 * (d_j[pairs[, 1]] - d_j[pairs[, 2]]))
  cj <- data.frame(a = pairs[, 1], b = pairs[, 2],
                   winner = ifelse(stats::runif(300) < p_a, pairs[, 1], pairs[, 2]),
                   stringsAsFactors = FALSE)
  rk <- NULL
  if (rankings)
    rk <- do.call(rbind, lapply(1:40, function(j) {
      rem <- sample(names(delta), 4); ord <- character(0)
      while (length(rem) > 1) {
        pick <- sample(rem, 1, prob = exp(1.1 * d_j[rem]))
        ord <- c(ord, pick); rem <- setdiff(rem, pick)
      }
      data.frame(ranking = j, item = c(ord, rem), rank = 1:4,
                 stringsAsFactors = FALSE)
    }))
  rasch_cj(X, comparisons = cj, rankings = rk, object_a = "a", object_b = "b",
           winner = "winner")
}

test_that("plot_cj draws the combined and separate calibrations", {
  fit <- cj_plot_sim()
  expect_silent(plot_cj(fit))
  expect_silent(plot_cj(fit, frames = FALSE))
  expect_null(plot_cj(fit))
  # a moving item under three frames takes the highlighted path
  moved <- cj_plot_sim(shift = 2, rankings = TRUE)
  tab <- moved$invariance$items
  expect_true("I03" %in% tab$item[!is.na(tab$p_adj) & tab$p_adj < 0.05])
  expect_silent(plot_cj(moved))
})

test_that("plot_cj refuses what it cannot draw", {
  fit <- cj_plot_sim()
  expect_error(plot_cj(list()), "rasch_cj")
  expect_error(plot_cj(fit, frames = NA), "`frames`")
  expect_error(plot_cj(fit, frames = c(TRUE, FALSE)), "`frames`")
  stalled <- fit; stalled$converged <- FALSE
  expect_error(plot_cj(stalled), "did not converge")
  persons <- fit; persons$mode <- "persons"
  expect_error(plot_cj(persons), "person-mode")
})
