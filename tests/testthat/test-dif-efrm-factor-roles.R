.make_efrm_factor_role_fit <- function(factor_name, seed = 881L) {
  d <- simulate_efrm(n_per_group = 100, items_per_set = 4,
                     n_sets = 1, n_groups = 2, seed = seed)
  truth <- attr(d, "truth")
  names(d)[names(d) == "group"] <- "frame"
  d[[factor_name]] <- factor(rep(c("A", "B"), length.out = nrow(d)))
  fit <- rasch_efrm(
    d, item_sets = truth$item_sets, groups = "frame", id = "id",
    factors = factor_name, boot_reps = 0, workers = 1)
  list(fit = fit, factor_name = factor_name)
}

test_that("EFRM DIF refits preserve named frame and factor roles", {
  fits <- lapply(c("group", "cohort"), .make_efrm_factor_role_fit)
  refits <- lapply(fits, function(z) {
    da <- dif_anova(z$fit, factors = z$factor_name, n_groups = 2)
    rasch:::.dif_boot_refit_efrm(z$fit$X, z$fit, da$bootstrap_design)
  })

  expect_true(all(vapply(refits, function(z)
    identical(rasch:::.fit_boot_status(z), "ok"), logical(1))))
  expect_true(all(vapply(refits, function(z)
    inherits(z$dif, "rasch_dif"), logical(1))))
  expect_identical(refits[[1L]]$dif$factor_names, "group")
  expect_identical(refits[[2L]]$dif$factor_names, "cohort")
  expect_equal(refits[[1L]]$dif$terms$F_value,
               refits[[2L]]$dif$terms$F_value, tolerance = 1e-12)
  expect_equal(refits[[1L]]$dif$terms$p,
               refits[[2L]]$dif$terms$p, tolerance = 1e-12)
  expect_equal(fits[[1L]]$fit$person$theta,
               fits[[2L]]$fit$person$theta, tolerance = 1e-12)

  # The public path must retain the complete factor family as well.
  da <- dif_anova(fits[[1L]]$fit,
                  factors = c("frame", "group"), n_groups = 2)
  expect_identical(da$factor_names, "group")
  expect_true(any(grepl("excluded", da$notes, fixed = TRUE)))
  db <- suppressWarnings(dif_bootstrap(fits[[1L]]$fit, da, B = 2,
                                       workers = 1, seed = 113))
  expect_identical(db$B_used, 2L)
  da2 <- dif_anova(fits[[2L]]$fit, factors = "cohort", n_groups = 2)
  db2 <- suppressWarnings(dif_bootstrap(fits[[2L]]$fit, da2, B = 2,
                                        workers = 1, seed = 113))
  expect_equal(unname(db$replicates$F), unname(db2$replicates$F),
               tolerance = 1e-12)
  expect_equal(unname(db$replicates$p), unname(db2$replicates$p),
               tolerance = 1e-12)
})

test_that("EFRM DIF refits preserve crossed frame names and collision factors", {
  d <- simulate_efrm(n_per_group = 80, items_per_set = 3,
                     n_sets = 1, n_groups = 2, seed = 882)
  truth <- attr(d, "truth")
  d$frame_a <- rep(c("A", "B"), each = nrow(d) / 2L)
  d$frame_b <- rep(rep(c("X", "Y"), each = nrow(d) / 4L), 2L)
  d$group <- NULL
  d[["frame_a:frame_b"]] <- factor(
    rep(c("L", "H"), length.out = nrow(d)))
  fit <- rasch_efrm(
    d, item_sets = truth$item_sets,
    groups = c("frame_a", "frame_b"), id = "id",
    factors = "frame_a:frame_b", boot_reps = 0, workers = 1)
  expect_setequal(fit$frame_group,
                  c(".frame_group", "frame_a", "frame_b"))
  da <- dif_anova(fit, factors = "frame_a:frame_b", n_groups = 2)
  expect_identical(da$factor_names, "frame_a:frame_b")

  refit <- rasch:::.dif_boot_refit_efrm(fit$X, fit, da$bootstrap_design)
  expect_identical(rasch:::.fit_boot_status(refit), "ok")
  expect_identical(refit$dif$factor_names, "frame_a:frame_b")
  db <- suppressWarnings(dif_bootstrap(fit, da, B = 1,
                                       workers = 1, seed = 114))
  expect_identical(db$B_used, 1L)
})

test_that("EFRM DIF refits allow an external factor to share an item name", {
  d <- simulate_efrm(n_per_group = 100, items_per_set = 4,
                     n_sets = 1, n_groups = 2, seed = 883)
  truth <- attr(d, "truth")
  X <- as.matrix(d[truth$item_sets[[1L]]])
  colnames(X) <- sprintf("I%02d", seq_len(ncol(X)))
  item_sets <- list(set1 = colnames(X))
  factors <- data.frame(I01 = factor(rep(c("A", "B"),
                                         length.out = nrow(X))))
  fit <- rasch_efrm(X, item_sets = item_sets, groups = d$group,
                    id = d$id, factors = factors,
                    boot_reps = 0, workers = 1)
  da <- dif_anova(fit, factors = "I01", n_groups = 2)
  refit <- rasch:::.dif_boot_refit_efrm(fit$X, fit,
                                        da$bootstrap_design)
  expect_identical(rasch:::.fit_boot_status(refit), "ok")
  expect_identical(refit$dif$factor_names, "I01")
  db <- suppressWarnings(dif_bootstrap(fit, da, B = 1,
                                       workers = 1, seed = 115))
  expect_identical(db$B_used, 1L)
})
