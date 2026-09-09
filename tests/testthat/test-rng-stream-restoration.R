.stream_rng_setup <- function() {
  old_kind <- RNGkind()
  old_seed <- .sim_seed_capture()
  list(kind = old_kind, seed = old_seed)
}

.stream_rng_reset <- function(old) {
  suppressWarnings(do.call(RNGkind, as.list(old$kind)))
  .sim_seed_restore(old$seed)
}

.expect_unchanged_stream <- function(operation, normal, refused = FALSE) {
  RNGkind("Mersenne-Twister", normal, "Rejection")
  set.seed(113)
  invisible(rnorm(1)) # leave one cached normal when Box-Muller is used
  expected <- rnorm(4)
  set.seed(113)
  invisible(rnorm(1))
  before <- .Random.seed
  if (refused) expect_error(operation(), "cannot preserve the Box-Muller")
  else operation()
  expect_identical(.Random.seed, before)
  expect_identical(rnorm(4), expected)
  expect_identical(RNGkind()[2L], normal)
}

test_that("seeded simulators refuse Box-Muller without changing its cache", {
  old <- .stream_rng_setup()
  on.exit(.stream_rng_reset(old), add = TRUE)
  calls <- list(
    function() simulate_rasch(40, 5, seed = 932),
    function() simulate_btl(5, 10, 5, seed = 932),
    function() simulate_mfrm(30, 4, 3, seed = 932),
    function() simulate_efrm(60, 4, 2, 2, seed = 932),
    function() simulate_btl_efrm(4, 2, 10, 2, 5, 5, seed = 932))
  for (operation in calls) {
    .expect_unchanged_stream(operation, "Inversion")
    .expect_unchanged_stream(operation, "Box-Muller", refused = TRUE)
  }
})

test_that("bootstrap refusal precedes draws with and without an explicit seed", {
  old <- .stream_rng_setup()
  on.exit(.stream_rng_reset(old), add = TRUE)
  fit <- rasch(simulate_rasch(180, 6, seed = 932),
               factors = data.frame(cohort = rep(c("A", "B"), 90)))
  da <- dif_anova(fit, n_groups = 2)
  bt <- btl(simulate_btl(5, 12, 10, seed = 419),
            "object_a", "object_b", winner = "winner", judge = "judge")
  for (seed in list(713L, NULL)) {
    operations <- list(
      function() fit_bootstrap(fit, B = 2, workers = 1, seed = seed),
      function() fit_bootstrap(bt, B = 2, workers = 1, seed = seed),
      function() dif_bootstrap(fit, da, B = 2, workers = 1, seed = seed),
      function() .dim_bootstrap(fit, 1:3, 4:6, TRUE, 1, .05,
                                B = 2, workers = 1, seed = seed))
    for (operation in operations)
      .expect_unchanged_stream(operation, "Box-Muller", refused = TRUE)
    .expect_unchanged_stream(function()
      sim_replicate(simulate_rasch, 2, n_persons = 40, n_items = 5,
                     seed = seed), "Box-Muller", refused = TRUE)
  }
  .expect_unchanged_stream(function()
    suppressWarnings(fit_bootstrap(fit, B = 2, workers = 1, seed = 713)),
    "Inversion")
  .expect_unchanged_stream(function()
    expect_error(simulate_rasch(-1, 5, seed = 713), "n_persons"), "Inversion")
})

test_that("unseeded simulations continue to support Box-Muller", {
  old <- .stream_rng_setup()
  on.exit(.stream_rng_reset(old), add = TRUE)
  RNGkind("Mersenne-Twister", "Box-Muller", "Rejection")
  set.seed(122)
  d <- simulate_rasch(40, 5)
  set.seed(122)
  expect_identical(simulate_rasch(40, 5), d)
})

test_that("extended-frame seed paths use the same stream guard", {
  old <- .stream_rng_setup()
  on.exit(.stream_rng_reset(old), add = TRUE)
  d <- simulate_efrm(120, 5, 1, 2, seed = 881)
  sets <- attr(d, "truth")$item_sets
  ef <- rasch_efrm(d, item_sets = sets, groups = "group", id = "id",
                   boot_reps = 0, workers = 1)
  cj <- simulate_btl_efrm(4, 2, 10, 2, 5, 5, seed = 932)
  cj_sets <- attr(cj, "truth")$object_sets
  operations <- list(
    function() rasch_efrm(d, item_sets = sets, groups = "group", id = "id",
                          boot_reps = 0, workers = 1, seed = 123),
    function() btl_efrm(cj, "object_a", "object_b", winner = "winner",
                        judge = "judge", panels = "panel",
                        object_sets = cj_sets, boot_reps = 0, workers = 1,
                        se_method = "conditional", seed = 123),
    function() frame_invariance(ef, se_method = "bootstrap", boot_reps = 30,
                                seed = 123))
  for (operation in operations)
    .expect_unchanged_stream(operation, "Box-Muller", refused = TRUE)
})
