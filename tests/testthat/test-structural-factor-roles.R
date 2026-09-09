make_structural_collision_efrm <- function(factor_name = "I01") {
  d <- simulate_efrm(100, 4, 1, 2, seed = 882)
  X <- as.matrix(d[, 2:5])
  colnames(X) <- sprintf("I%02d", 1:4)
  factors <- data.frame(value = rep(c("A", "B"), length.out = nrow(X)),
                        stringsAsFactors = FALSE)
  names(factors) <- factor_name
  ids <- sprintf("E%03d", seq_len(nrow(X)))
  fit <- rasch_efrm(
    X, item_sets = list(S1 = colnames(X)), groups = d$group,
    id = ids, factors = factors, n_groups = 6,
    se_method = "hybrid", boot_reps = 0, workers = 1, seed = 19,
    maxit = 41, tol = 2e-7)
  list(fit = fit, factors = factors, ids = ids)
}

test_that("ordinary structural refits keep separate item and factor roles", {
  d <- simulate_rasch(250, 6, seed = 935)
  X <- as.matrix(d[, setdiff(names(d), "id")])
  colnames(X)[2L] <- "id"
  ids <- sprintf("A%03d", seq_len(nrow(X)))
  factors <- data.frame(I01 = rep(c("G1", "G2"), length.out = nrow(X)))
  anchors <- data.frame(item = "I01", k = 1L, tau = 0.25)
  fit <- rasch(X, id = ids, factors = factors, anchors = anchors,
               n_groups = 6, maxit = 75, tol = 5e-9)

  dropped <- drop_items(fit, "I06")
  combined <- combine_items(fit, c("I05", "I06"))
  for (refit in list(dropped, combined)) {
    expect_true(refit$est$converged)
    expect_identical(names(refit$factors), "I01")
    expect_identical(as.character(refit$factors$I01),
                     as.character(factors$I01))
    expect_identical(as.character(refit$person$I01),
                     as.character(factors$I01))
    expect_identical(as.character(refit$person$id), ids)
    expect_true("id" %in% refit$items$item)
    expect_identical(as.character(refit$est$anchors$item), "I01")
    expect_equal(refit$refit_spec$n_groups, 6)
    expect_identical(refit$refit_spec$maxit, 75)
    expect_identical(refit$refit_spec$tol, 5e-9)
    expect_false(any(grepl(".rasch_refit_factor_",
                           c(names(refit$factors), names(refit$person)),
                           fixed = TRUE)))
  }
})

test_that("keyed mixed-score refits retain a same-named external factor", {
  set.seed(934)
  n <- 300
  raw <- data.frame(
    M1 = sample(c("A", "B", "C"), n, TRUE),
    M2 = sample(c("A", "B", "C"), n, TRUE),
    I3 = rbinom(n, 1, .5), I4 = rbinom(n, 1, .5),
    check.names = FALSE)
  X <- as.matrix(raw)
  key <- data.frame(
    item = rep(c("M1", "M2"), each = 2L),
    option = rep(c("A", "B"), 2L), score = rep(c(2L, 1L), 2L))
  factors <- data.frame(M1 = rep(c("G1", "G2"), length.out = n))
  ids <- sprintf("K%03d", seq_len(n))
  fit <- rasch(X, id = ids, factors = factors, key = key,
               n_groups = 7, maxit = 80, tol = 1e-9)

  dropped <- drop_items(fit, "I4")
  combined <- combine_items(fit, c("M2", "I3"))
  expect_equal(dropped$X[, c("M1", "M2")], fit$X[, c("M1", "M2")])
  expect_equal(dropped$mc$map, fit$mc$map)
  expect_identical(colnames(combined$mc$raw), "M1")
  expect_equal(combined$mc$map$M1, fit$mc$map$M1)
  for (refit in list(dropped, combined)) {
    expect_identical(names(refit$factors), "M1")
    expect_identical(as.character(refit$factors$M1),
                     as.character(factors$M1))
    expect_identical(as.character(refit$person$M1),
                     as.character(factors$M1))
    expect_identical(as.character(refit$person$id), ids)
    expect_equal(refit$refit_spec$n_groups, 7)
    expect_identical(refit$refit_spec$maxit, 80)
    expect_identical(refit$refit_spec$tol, 1e-9)
    expect_false(any(grepl(".rasch_refit_factor_",
                           c(names(refit$factors), names(refit$person)),
                           fixed = TRUE)))
  }
})

test_that("EFRM structural refits preserve colliding factor and frame metadata", {
  made <- make_structural_collision_efrm()
  fit <- made$fit
  dropped <- drop_items(fit, "I04", boot_reps = 0)
  resolved <- resolve_frames(fit, "I02", boot_reps = 0)

  for (refit in list(dropped, resolved)) {
    expect_true(refit$est$converged)
    expect_identical(names(refit$factors), c("group", "I01"))
    expect_identical(as.character(refit$factors$I01),
                     as.character(made$factors$I01))
    expect_identical(as.character(refit$person$I01),
                     as.character(made$factors$I01))
    expect_identical(refit$frame_group, fit$frame_group)
    expect_identical(as.character(refit$person$id), made$ids)
    expect_identical(refit$refit_spec$groups, "group")
    expect_identical(refit$refit_spec$factors, "I01")
    expect_identical(refit$refit_spec$n_groups, 6L)
    expect_identical(refit$refit_spec$maxit, 41)
    expect_identical(refit$refit_spec$tol, 2e-7)
    expect_identical(refit$refit_spec$se_method, "hybrid")
    expect_identical(refit$refit_spec$boot_reps, 0L)
    expect_identical(refit$refit_spec$workers, 1L)
    expect_identical(refit$refit_spec$seed, 19L)
    expect_false(any(grepl(".rasch_refit_factor_",
                           c(names(refit$factors), names(refit$person),
                             refit$frame_group), fixed = TRUE)))
  }
  expect_identical(colnames(dropped$X),
                   fit$virtual_map$vkey[fit$virtual_map$item != "I04"])
  expect_setequal(grep("I02", names(resolved$set_of), value = TRUE),
                  c("I02 (g1)", "I02 (g2)"))
})

test_that("crossed punctuated EFRM frame roles survive a structural replay", {
  d <- simulate_efrm(100, 5, 1, 2, seed = 913)
  truth <- attr(d, "truth")
  items <- unlist(truth$item_sets, use.names = FALSE)
  d[["age:band"]] <- rep(c("North", "South"), length.out = nrow(d))
  d$g2 <- as.character(d$group)
  d$cohort <- rep(c("C1", "C2"), each = nrow(d) / 2L)
  fit <- rasch_efrm(
    d, item_sets = truth$item_sets, groups = c("age:band", "g2"), id = "id",
    factors = "cohort", items = items, boot_reps = 0,
    maxit = 44, tol = 3e-7)

  refit <- drop_items(fit, tail(items, 1L), boot_reps = 0)
  expect_true(refit$est$converged)
  expect_identical(refit$frame_group,
                   c("age:band:g2", "age:band", "g2"))
  expect_identical(names(refit$factors),
                   c("age:band:g2", "age:band", "g2", "cohort"))
  expect_identical(refit$refit_spec$groups, c("age:band", "g2"))
  expect_identical(refit$refit_spec$factors, "cohort")
  expect_false(is.null(refit$phi_factorial))
  expect_false(is.null(refit$phi_factorial_tests))
  expect_true(any(grepl("`age:band`", refit$phi_factorial$term,
                        fixed = TRUE)))
  expect_false(any(grepl("``age:band``", refit$phi_factorial$term,
                         fixed = TRUE)))
  expect_false(any(grepl(".rasch_refit_factor_",
                         c(names(refit$factors), names(refit$person),
                           refit$frame_group, refit$phi_factorial$term,
                           refit$phi_factorial_tests$term), fixed = TRUE)))
})

test_that("factorial alias restoration is one-pass and prefix-safe", {
  generated <- .structural_factor_aliases(
    data.frame(I1 = "A", ".rasch_refit_factor_1__suffix" = "B",
               check.names = FALSE), reserved = "I1")
  expect_false(any(vapply(generated$original, function(x)
    grepl(generated$replay[1L], x, fixed = TRUE), logical(1))))

  aliases <- list(
    original = c("age:band", "region", "cohort"),
    replay = c(".rasch_refit_factor_1", ".rasch_refit_factor_11", "cohort"))
  fake <- list(
    factors = data.frame(
      .rasch_refit_factor_1 = "A", .rasch_refit_factor_11 = "B",
      cohort = "C", check.names = FALSE),
    person = data.frame(
      id = 1, .rasch_refit_factor_1 = "A",
      .rasch_refit_factor_11 = "B", cohort = "C", check.names = FALSE),
    frame_group = c("cell", ".rasch_refit_factor_1",
                    ".rasch_refit_factor_11"),
    phi_factorial = data.frame(term = c(
      ".rasch_refit_factor_1A",
      ".rasch_refit_factor_11B",
      ".rasch_refit_factor_1A:.rasch_refit_factor_11B")),
    phi_factorial_tests = data.frame(term = c(
      ".rasch_refit_factor_1", ".rasch_refit_factor_11")))
  restored <- .restore_structural_factor_names(
    fake, aliases, c("cell", "age:band", "region"))
  expect_identical(names(restored$factors), c("age:band", "region", "cohort"))
  expect_identical(names(restored$person),
                   c("id", "age:band", "region", "cohort"))
  expect_identical(restored$phi_factorial$term,
    c("`age:band`A", "regionB", "`age:band`A:regionB"))
  expect_identical(restored$phi_factorial_tests$term,
                   c("`age:band`", "region"))

  unchanged <- list(
    factors = data.frame("age:band" = "A", check.names = FALSE),
    person = data.frame(id = 1, "age:band" = "A", check.names = FALSE),
    frame_group = c("cell", "age:band"),
    phi_factorial = data.frame(term = "`age:band`A"),
    phi_factorial_tests = data.frame(term = "`age:band`"))
  same <- .restore_structural_factor_names(unchanged,
    list(original = "age:band", replay = "age:band"),
    c("cell", "age:band"))
  expect_identical(same$phi_factorial$term, "`age:band`A")
  expect_identical(same$phi_factorial_tests$term, "`age:band`")
})

test_that("frame bootstrap inference is invariant to an item-named factor", {
  collision <- make_structural_collision_efrm("I01")$fit
  renamed <- make_structural_collision_efrm("cohort")$fit
  one <- frame_invariance(collision, se_method = "bootstrap",
                          boot_reps = 30, seed = 27)
  two <- frame_invariance(renamed, se_method = "bootstrap",
                          boot_reps = 30, seed = 27)

  expect_identical(one$boot_reps_used, 30L)
  expect_identical(one$boot_reps_errors, 0L)
  expect_identical(two$boot_reps_used, 30L)
  expect_identical(two$boot_reps_errors, 0L)
  expect_equal(one$locations, two$locations, tolerance = 1e-10)
  expect_equal(one$discrimination, two$discrimination, tolerance = 1e-10)
  expect_equal(one$summary, two$summary, tolerance = 1e-10)
})
