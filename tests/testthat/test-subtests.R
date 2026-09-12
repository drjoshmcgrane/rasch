
test_that("drop_items removes items and refits, for rasch and efrm", {
  d <- simulate_rasch(300, 8, seed = 1)
  f <- rasch(d, id = "id")
  f2 <- drop_items(f, "I03")
  expect_equal(nrow(f2$items), nrow(f$items) - 1L)
  expect_false("I03" %in% f2$items$item)
  expect_match(tail(f2$notes, 1), "dropped item")
  expect_s3_class(f2, "rasch")

  e <- simulate_efrm(n_per_group = 250, items_per_set = 8, n_sets = 2,
                     n_groups = 2, set_unit_ratio = 1.3, seed = 7)
  tr <- attr(e, "truth")
  fe <- rasch_efrm(e, item_sets = tr$item_sets, groups = "group", id = "id",
                   boot_reps = 0)
  fe2 <- drop_items(fe, "S1I03")
  expect_s3_class(fe2, "rasch_efrm")
  expect_false("S1I03" %in% names(fe2$set_of))
  expect_equal(length(names(fe2$set_of)), length(names(fe$set_of)) - 1L)
  expect_setequal(unique(fe2$set_of), unique(fe$set_of))   # both sets survive
  expect_equal(nrow(fe2$person), nrow(fe$person))          # persons preserved

  # guards
  expect_error(drop_items(fe, "NOPE"), "not in the fit")
  expect_error(drop_items(fe, names(fe$set_of)[1:8]), "empty set")
  expect_error(drop_items(f, character()), "at least one non-missing item")
})

test_that("drop_items keeps the standard errors the source fit had", {
  # Structural refits reproduce the requested uncertainty calculation. The
  # successful count is reported separately and may be smaller.
  d <- simulate_efrm(n_per_group = 400, items_per_set = 6, n_sets = 2,
                     n_groups = 1, set_unit_ratio = 1.4, seed = 3)
  tr <- attr(d, "truth")
  drop1 <- tr$item_sets[[1]][1]

  with_se <- rasch_efrm(d, item_sets = tr$item_sets, groups = "group",
                        id = "id")
  expect_true(any(is.finite(with_se$alpha_table$se_log_alpha)))
  expect_identical(with_se$boot_reps_requested, 300L)
  expect_gte(with_se$boot_reps_used, 150L)
  expect_identical(with_se$boot_reps_failed,
                   with_se$boot_reps_requested - with_se$boot_reps_used)
  kept <- drop_items(with_se, drop1)
  expect_true(all(is.finite(kept$alpha_table$se_log_alpha)))
  # The two centred set units are one hypothesis: the first row carries the
  # adjusted probability and the second, which restates it, is withheld.
  kept_adj <- kept$efrm_vs_rasch$unit_tests$p_adj
  expect_true(is.finite(kept_adj[1]))
  expect_true(is.na(kept_adj[2]))

  # a fit asked for no standard errors keeps none, and stays cheap
  without <- rasch_efrm(d, item_sets = tr$item_sets, groups = "group",
                        id = "id", boot_reps = 0)
  expect_true(all(is.na(without$alpha_table$se_log_alpha)))
  expect_identical(without$boot_reps_requested, 0L)
  expect_identical(without$boot_reps_used, 0L)
  expect_identical(without$boot_reps_failed, 0L)
  expect_true(all(is.na(drop_items(without, drop1)$alpha_table$se_log_alpha)))
  # and an explicit zero still wins over the source fit's own standard errors
  expect_true(all(is.na(
    drop_items(with_se, drop1, boot_reps = 0)$alpha_table$se_log_alpha)))
})
