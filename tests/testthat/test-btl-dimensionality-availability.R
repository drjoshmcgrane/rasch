check_btl_dimensionality_withheld <- function(result, fit, reps = 20L) {
  expect_true(is.na(result$leading_structured))
  expect_false(result$reference$inference_available)
  expect_true(all(is.finite(result$bimensions$strength)))
  expect_true(all(is.finite(result$residual_matrix)))
  expect_true(all(is.na(unlist(result$reference[
    c("mean", "p95", "p", "p_adj")]))))
  expect_true(all(is.na(result$bimensions$ref_mean)))
  expect_true(all(is.na(result$bimensions$ref_p95)))
  expect_length(result$reference$draws, reps)
  expect_identical(result$reference$n_used, reps)
  expect_no_error(.validate_btl_dimensionality(result, fit))
  printed <- capture.output(print(result))
  expect_true(any(grepl("inference withheld", printed, fixed = TRUE)))
  expect_false(any(grepl("adjusted p =", printed, fixed = TRUE)))
}

test_that("incomplete BTL pair coverage withholds the whole inference display", {
  d <- simulate_btl(5, 15, reps_per_pair = 30, seed = 9681)
  pair <- sort(c(d$object_a[1L], d$object_b[1L]))
  removed <- (d$object_a == pair[1L] & d$object_b == pair[2L]) |
    (d$object_a == pair[2L] & d$object_b == pair[1L])
  incomplete <- d[!removed, ]
  fit <- btl(incomplete, "object_a", "object_b", winner = "winner",
              judge = "judge")
  result <- btl_dimensionality(fit, reps = 20L, seed = 9691,
                              independent_comparisons = TRUE)
  check_btl_dimensionality_withheld(result, fit)

  complete_fit <- btl(d, "object_a", "object_b", winner = "winner",
                       judge = "judge")
  complete <- btl_dimensionality(complete_fit, reps = 20L, seed = 9691,
                                independent_comparisons = TRUE)
  expect_true(complete$reference$inference_available)
  expect_true(all(is.finite(unlist(complete$reference[
    c("mean", "p95", "p", "p_adj")]))))
  expect_identical(complete$leading_structured,
                   complete$reference$p_adj <= 0.05)

  # A correctly signed result from an older version can still contain the
  # withheld reference. Do not let it re-enter reports or an app cache.
  legacy <- result
  legacy$reference$inference_available <- NULL
  legacy$reference$mean <- mean(legacy$reference$draws)
  legacy$reference$p95 <- as.numeric(quantile(legacy$reference$draws, .95))
  legacy$reference$p <- legacy$reference$p_adj <- .01
  legacy$bimensions$ref_mean[1L] <- legacy$reference$mean
  legacy$bimensions$ref_p95[1L] <- legacy$reference$p95
  attr(legacy, "result_signature") <- NULL
  attr(legacy, "result_signature") <- .fit_boot_md5(legacy)
  expect_error(.validate_btl_dimensionality(legacy, fit), "recompute")
  expect_false(any(grepl("adjusted p =", capture.output(print(legacy)),
                         fixed = TRUE)))

  legend <- NULL
  testthat::local_mocked_bindings(
    .rr_legend = function(position, labels, ...) legend <<- labels,
    .package = "rasch")
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off())
  expect_no_error(plot_btl_scree(result))
  expect_identical(legend, "Observed")
  expect_no_error(plot_btl_scree(legacy))
  expect_identical(legend, "Observed")
  expect_no_error(plot_btl_scree(complete))
  expect_identical(legend, c("Observed", "Null reference band"))
})

test_that("shared comparison order cannot retain a dimensionality probability", {
  objects <- LETTERS[1:6]
  beta <- setNames(seq(-1.2, 1.2, length.out = 6L), objects)
  pairs <- t(utils::combn(objects, 2L))
  rows <- list()
  set.seed(1)
  for (j in seq_len(12L)) {
    last <- NA_character_
    for (r in seq_len(nrow(pairs))) {
      a <- pairs[r, 1L]; b <- pairs[r, 2L]
      eta <- beta[a] - beta[b]
      if (!is.na(last) && last == a) eta <- eta + 1.5
      if (!is.na(last) && last == b) eta <- eta - 1.5
      winner <- if (runif(1L) < plogis(eta)) a else b
      last <- winner
      rows[[length(rows) + 1L]] <- data.frame(
        a, b, winner, judge = sprintf("J%02d", j), order = r)
    }
  }
  fit <- btl(do.call(rbind, rows), "a", "b", winner = "winner",
              judge = "judge", order = "order")
  expect_true(.btl_order_variation(fit$dependence_data,
                                   fit$objects$object)$shared)
  result <- btl_dimensionality(fit, reps = 20L, seed = 9691,
                              independent_comparisons = TRUE)
  check_btl_dimensionality_withheld(result, fit)
  expect_match(paste(result$notes, collapse = " "), "confounded")
})

test_that("BTL frame dimensionality applies the incomplete-pair guard", {
  d <- simulate_btl_efrm(4, 2, 5, 2, 10, 10, seed = 71)
  sets <- attr(d, "truth")$object_sets
  missing_pair <- sets[[1L]][1:2]
  removed <- (d$object_a == missing_pair[1L] & d$object_b == missing_pair[2L]) |
    (d$object_a == missing_pair[2L] & d$object_b == missing_pair[1L])
  fit <- btl_efrm(d[!removed, ], "object_a", "object_b", winner = "winner",
    judge = "judge", panels = "panel", object_sets = sets,
    se_method = "conditional", boot_reps = 0)
  result <- btl_dimensionality(fit, reps = 20L, seed = 9691,
                              independent_comparisons = TRUE)
  check_btl_dimensionality_withheld(result, fit)
})

test_that("judge-clustered dimensionality requires an explicit independence assumption", {
  d <- simulate_btl(6, 30, reps_per_pair = 30, seed = 15)
  fit <- btl(d, "object_a", "object_b", "winner", judge = "judge")
  default <- btl_dimensionality(fit, reps = 20L, seed = 8)
  check_btl_dimensionality_withheld(default, fit)
  expect_false(default$reference$independent_comparisons)
  expect_match(paste(default$notes, collapse = " "), "not supply a cluster-robust")

  sensitivity <- btl_dimensionality(fit, reps = 20L, seed = 8,
                                     independent_comparisons = TRUE)
  expect_true(sensitivity$reference$inference_available)
  expect_true(sensitivity$reference$independent_comparisons)
  expect_true(is.finite(sensitivity$reference$p))
  expect_match(paste(sensitivity$notes, collapse = " "), "not cluster-robust")
  expect_equal(default$residual_matrix, sensitivity$residual_matrix)
  expect_equal(default$bimensions$strength, sensitivity$bimensions$strength)
  expect_equal(default$reference$draws, sensitivity$reference$draws)
  expect_no_error(.validate_btl_dimensionality(sensitivity, fit))

  repeated <- d[rep(seq_len(nrow(d)), each = 10L), ]
  fit_repeated <- btl(repeated, "object_a", "object_b", "winner", judge = "judge")
  expect_equal(fit$objects$location, fit_repeated$objects$location, tolerance = 1e-8)
  expect_equal(fit$objects$se, fit_repeated$objects$se, tolerance = 1e-8)
  check_btl_dimensionality_withheld(
    btl_dimensionality(fit_repeated, reps = 20L, seed = 8), fit_repeated)

  legacy <- sensitivity
  legacy$reference$independent_comparisons <- NULL
  attr(legacy, "result_signature") <- NULL
  attr(legacy, "result_signature") <- .fit_boot_md5(legacy)
  expect_error(.validate_btl_dimensionality(legacy, fit), "recompute")
})

test_that("the dimensionality assumption has a strict and reversible API", {
  d <- simulate_btl(5, 15, reps_per_pair = 30, seed = 9681)
  fit <- btl(d, "object_a", "object_b", "winner")
  automatic <- btl_dimensionality(fit, reps = 20L, seed = 3)
  expect_true(automatic$reference$inference_available)
  expect_true(automatic$reference$independent_comparisons)
  check_btl_dimensionality_withheld(btl_dimensionality(fit,
    reps = 20L, seed = 3, independent_comparisons = FALSE), fit)
  for (bad in list(NA, 1, "TRUE", c(TRUE, FALSE), logical()))
    expect_error(btl_dimensionality(fit, reps = 20L,
      independent_comparisons = bad), "must be TRUE, FALSE, or NULL")

  frame_data <- simulate_btl_efrm(4, 2, 5, 2, 10, 10, seed = 71)
  frame <- btl_efrm(frame_data, "object_a", "object_b", "winner", "judge",
    "panel", object_sets = attr(frame_data, "truth")$object_sets,
    se_method = "conditional")
  check_btl_dimensionality_withheld(
    btl_dimensionality(frame, reps = 20L, seed = 3), frame)
  explicit <- btl_dimensionality(frame, reps = 20L, seed = 3,
                                 independent_comparisons = TRUE)
  expect_true(explicit$reference$inference_available)
  expect_no_error(.validate_btl_dimensionality(explicit, frame))
})
