test_that("BTL recovery uses comparisons retained before boundary removal", {
  d <- simulate_btl(5, 10, reps_per_pair = 2, seed = 1)
  fit <- btl(d, "object_a", "object_b", winner = "winner", judge = "judge")

  # O5 is undefeated for this deterministic design and is reported as an
  # extrapolated boundary row.  Its source rows still belong to recovery.
  expect_true(any(fit$objects$extreme %in% TRUE))
  expect_lt(nrow(fit$comparisons), nrow(fit$observed_comparisons))
  rec <- sim_recovery(fit, d)
  expect_s3_class(rec, "rasch_recovery")
  ext <- as.character(fit$objects$object[fit$objects$extreme %in% TRUE])
  expect_false(any(rec$pieces[["object location"]]$label %in% ext))
  expect_equal(nrow(rec$pieces[["object location"]]),
               sum(!(fit$objects$extreme %in% TRUE)))

  # The changed row is one of O5's source comparisons, which is absent from
  # the fitted likelihood.  It must still be rejected by provenance checks.
  boundary_row <- which(d$object_a == "O5" | d$object_b == "O5")[1L]
  changed_boundary <- d
  changed_boundary$winner[boundary_row] <- if (
    changed_boundary$winner[boundary_row] == changed_boundary$object_a[boundary_row])
    changed_boundary$object_b[boundary_row] else
      changed_boundary$object_a[boundary_row]
  expect_error(sim_recovery(fit, changed_boundary),
               "comparison, outcome or judge")

  # A winless object exercises the other boundary direction.  O5 remains
  # undefeated, so both boundary rows must be represented by the same source
  # provenance check.
  dw <- d
  oi <- dw$object_a == "O1" | dw$object_b == "O1"
  dw$winner[oi] <- ifelse(dw$object_a[oi] == "O1",
                          dw$object_b[oi], dw$object_a[oi])
  fw <- btl(dw, "object_a", "object_b", winner = "winner", judge = "judge")
  expect_true(sum(fw$objects$extreme %in% TRUE) >= 2L)
  expect_s3_class(sim_recovery(fw, dw), "rasch_recovery")
})

test_that("BTL recovery refuses a changed source and supports old fits", {
  d <- simulate_btl(5, 10, reps_per_pair = 10, seed = 1)
  fit <- btl(d, "object_a", "object_b", winner = "winner", judge = "judge")
  expect_false(any(fit$objects$extreme %in% TRUE))

  changed <- d
  changed$winner[1L] <- if (changed$winner[1L] == changed$object_a[1L])
    changed$object_b[1L] else changed$object_a[1L]
  expect_error(sim_recovery(fit, changed),
               "comparison, outcome or judge")

  # Saved fits predating observed_comparisons fall back to their fitted rows.
  old_fit <- fit
  old_fit$observed_comparisons <- NULL
  expect_s3_class(sim_recovery(old_fit, d), "rasch_recovery")
})

test_that("BTL recovery handles graded and count comparison representations", {
  graded <- simulate_btl(5, 10, reps_per_pair = 2, model = "graded",
                         n_categories = 4, seed = 1)
  gf <- btl(graded, "object_a", "object_b", response = "response",
            judge = "judge")
  expect_s3_class(sim_recovery(gf, graded), "rasch_recovery")

  # Count compression is an equivalent representation of the original
  # generated records; recovery expands its weights canonically.
  d <- simulate_btl(5, 10, reps_per_pair = 2, seed = 2)
  key <- paste(d$object_a, d$object_b, d$winner, d$judge, sep = "\034")
  first <- !duplicated(key)
  counted <- d[first, c("object_a", "object_b", "winner", "judge"),
               drop = FALSE]
  counted$count <- as.numeric(tabulate(match(key, key[first])))
  cf <- btl(counted, "object_a", "object_b", winner = "winner",
            judge = "judge", count = "count")
  expect_s3_class(sim_recovery(cf, d), "rasch_recovery")
})
