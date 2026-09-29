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

test_that("frame markers align subset and disconnected origins, not differences", {
  set.seed(9525)
  delta <- setNames(seq(-2, 2, length.out = 6), paste0("I", 1:6))
  theta <- rnorm(500)
  X <- sapply(delta, function(d) rbinom(500, 1, plogis(theta - d)))
  make_pairs <- function(objects) {
    pairs <- t(combn(objects, 2))
    pairs <- pairs[sample(nrow(pairs), 600, TRUE), ]
    data.frame(object_a = pairs[, 1], object_b = pairs[, 2],
      winner = ifelse(runif(600) < plogis(delta[pairs[, 1]] - delta[pairs[, 2]]),
                      pairs[, 1], pairs[, 2]))
  }
  high <- make_pairs(names(delta)[4:6])
  for (cmp in list(high, rbind(high, make_pairs(names(delta)[1:3])))) {
    fit <- rasch_cj(X, cmp, units = c(comparisons = 1))
    expect_true(fit$converged)
    resp <- subset(fit$frame_locations, frame == "responses")
    expect_equal(nrow(resp), ncol(X))
    expect_true(all(is.finite(resp$location)))
    expect_equal(mean(resp$location), mean(fit$items$location))
    loc <- subset(fit$frame_locations, frame == "comparisons")
    expect_equal(nrow(loc), length(unique(c(cmp$object_a, cmp$object_b))))
    for (block in unique(loc$block)) {
      at <- loc[loc$block == block, ]
      idx <- match(at$item, fit$items$item)
      expect_equal(mean(at$location), mean(fit$items$location[idx]))
      # Only the arbitrary origin changes; every within-block difference
      # is exactly the separate calibration's original difference.
      expect_equal(diff(at$location), diff(fit$items$location_comparisons[idx]))
      expect_gt(abs(mean(at$location)), .5)
    }
    # The plot uses these aligned locations, not the raw separate column.
    calls <- list()
    local_mocked_bindings(points = function(x, y, ...) {
      calls[[length(calls) + 1L]] <<- x
    }, .package = "rasch")
    plot_cj(fit)
    order_items <- fit$items$item[order(fit$items$location)]
    expect_equal(calls[[3]], loc$location[match(order_items, loc$item)])
    before <- fit$invariance
    plot_cj(fit)
    expect_identical(fit$invariance, before)
    old <- fit; old$frame_locations <- NULL
    expect_warning(plot_cj(old), "separate-frame origins are unavailable")
    expect_silent(plot_cj(old, frames = FALSE))
    stalled <- rasch_cj(X, cmp, units = c(comparisons = 1), maxit = 1)
    expect_false(stalled$converged)
    expect_equal(stalled$invariance$n_contrasts, nrow(loc))
    expect_null(stalled$invariance$items)
  }
})
