test_that("frame helpers support crossed factors and punctuated labels", {
  skip_on_cran()
  set.seed(13)
  N <- 800; K <- 8
  a <- rep(c("A+B", "C.D"), each = N / 2)
  b <- rep(rep(c("X", "Y"), each = N / 4), 2)
  theta <- rnorm(N, 0, 1.3)
  delta <- seq(-1.4, 1.4, length.out = K)
  X <- vapply(delta, function(d) rbinom(N, 1, plogis(theta - d)), numeric(N))
  colnames(X) <- paste0("item:", seq_len(K))
  d <- data.frame(id = seq_len(N), X, a = a, b = b, check.names = FALSE)
  f <- suppressWarnings(rasch_efrm(
    d, item_sets = list(core = colnames(X)), groups = c("a", "b"),
    id = "id", items = colnames(X), boot_reps = 0))

  expect_s3_class(frame_invariance(f), "rasch_frame_invariance")
  expect_s3_class(drop_items(f, "item:8", boot_reps = 0), "rasch_efrm")
  expect_s3_class(resolve_frames(f, "item:8", boot_reps = 0), "rasch_efrm")
})

test_that("frame invariance is unchanged by the first item's missingness and order", {
  d <- simulate_efrm(n_per_group = 350, items_per_set = 8, n_sets = 1,
                     n_groups = 2, seed = 19)
  tr <- attr(d, "truth")
  items <- unlist(tr$item_sets, use.names = FALSE)
  set.seed(20)
  d[sample(which(d$group == levels(factor(d$group))[1]), 45), items[1]] <- NA
  f1 <- rasch_efrm(d, item_sets = tr$item_sets, groups = "group", id = "id",
                   boot_reps = 0)
  inv1 <- frame_invariance(f1)

  rev_items <- rev(items)
  f2 <- rasch_efrm(d[, c("id", rev_items, "group")],
                   item_sets = list(set1 = rev_items), groups = "group",
                   id = "id", boot_reps = 0)
  inv2 <- frame_invariance(f2)
  a <- inv1$locations[order(inv1$locations$item), ]
  b <- inv2$locations[order(inv2$locations$item), ]
  expect_equal(a$difference, b$difference, tolerance = 1e-8)
  expect_equal(a$se, b$se, tolerance = 1e-8)
})

test_that("polytomous frame differences use the threshold-weighted origin", {
  d <- simulate_efrm(n_per_group = 450, items_per_set = 7, n_sets = 1,
                     n_groups = 2, n_categories = 4, seed = 37)
  tr <- attr(d, "truth")
  items <- unlist(tr$item_sets, use.names = FALSE)
  f <- rasch_efrm(d, item_sets = tr$item_sets, groups = "group", id = "id",
                  boot_reps = 0)
  z <- frame_invariance(f)
  vm <- f$virtual_map
  m <- vapply(z$locations$item, function(it) {
    j <- which(vm$item == it & vm$group == z$locations$frame_1[1])[1]
    f$m[match(vm$vkey[j], colnames(f$X))]
  }, 0)
  expect_equal(sum(m * z$locations$difference), 0, tolerance = 1e-8)

  g2 <- d$group == levels(factor(d$group))[2]
  d[g2 & d[[items[1]]] == 3, items[1]] <- 2
  f2 <- rasch_efrm(d, item_sets = tr$item_sets, groups = "group", id = "id",
                   boot_reps = 0)
  z2 <- frame_invariance(f2)
  expect_true(items[1] %in% z2$excluded$item)
  expect_equal(z2$locations$p_adj,
               p.adjust(z2$locations$p, "holm",
                        n = nrow(z2$locations) + nrow(z2$excluded)))
})

test_that("frame invariance withholds boundary slope tests without shrinking the family", {
  cmp <- data.frame(p = c(.01, .03))
  dsc <- data.frame(p = c(.02, .04), statistic = c(2.4, 2.1),
                    disc_boundary = c(TRUE, FALSE))
  excluded <- data.frame(item = "unavailable")
  z <- .frame_invariance_probabilities(
    cmp, dsc, excluded, "bootstrap", alpha = .05, adjust = "holm")
  # Two tested and one unavailable item contribute a location and a
  # discrimination comparison each: the predeclared family has six members.
  expect_identical(z$family_n, 6L)
  expect_true(is.na(z$discrimination$p[1]))
  expect_true(is.na(z$discrimination$statistic[1]))
  expect_true(is.na(z$discrimination$p_adj[1]))
  expect_equal(z$locations$p_adj,
               p.adjust(c(.01, .03, .04), "holm", n = 6L)[1:2])
  expect_equal(z$discrimination$p_adj[2],
               p.adjust(c(.01, .03, .04), "holm", n = 6L)[3])
})

test_that("frame invariance withholds Wald probabilities at zero uncertainty", {
  z <- .frame_invariance_wald(c(0.2, 0, 0.3), c(0, 0, 0.1))
  expect_true(all(is.na(z$statistic[1:2])))
  expect_true(all(is.na(z$p[1:2])))
  expect_equal(z$statistic[3], 3)
  expect_equal(z$p[3], 2 * pnorm(-3))
})

test_that("frame bootstrap comparisons keep the observed centring family", {
  observed <- data.frame(
    set = "S", frame_1 = "A", frame_2 = "B", item = c("I1", "I2"),
    difference = c(.2, -.2), disc_ratio = c(1.1, .9))
  key <- .factor_keys(observed[c("set", "frame_1", "frame_2", "item")])
  expect_equal(
    .frame_invariance_boot_vector(observed[2:1, ], key, "difference"),
    observed$difference)
  expect_equal(
    .frame_invariance_boot_vector(observed[2:1, ], key, "disc_ratio", log),
    log(observed$disc_ratio))
  expect_null(.frame_invariance_boot_vector(
    observed[1, ], key, "difference"))
  expect_null(.frame_invariance_boot_vector(
    rbind(observed, transform(observed[1, ], item = "I3")),
    key, "difference"))
  expect_null(.frame_invariance_boot_vector(
    rbind(observed[1, ], observed[1, ]), key, "difference"))
})

test_that("a sparse polytomous bootstrap cannot enlarge the frame family", {
  d <- simulate_efrm(n_per_group = 80, items_per_set = 6, n_sets = 1,
                     n_groups = 2, n_categories = 3, seed = 2026)
  truth <- attr(d, "truth")
  fit <- rasch_efrm(d, item_sets = truth$item_sets, groups = "group",
                    id = "id", boot_reps = 0)
  observed <- .frame_invariance_conditional(fit)
  recovered <- observed
  recovered$locations <- rbind(
    recovered$locations,
    transform(recovered$locations[1L, ], item = "sparse recovered item"))
  recovered$discrimination <- rbind(
    recovered$discrimination,
    transform(recovered$discrimination[1L, ],
              item = "sparse recovered item"))
  calls <- 0L
  conditional <- function(...) {
    calls <<- calls + 1L
    if (calls == 2L) recovered else observed
  }
  result <- testthat::with_mocked_bindings(
    frame_invariance(
      fit, se_method = "bootstrap", boot_reps = 40, seed = 11),
    .frame_invariance_conditional = conditional,
    .efrm_refit = function(...) fit,
    .package = "rasch")
  expect_identical(result$boot_reps_used, 39L)
  expect_identical(result$boot_reps_nonconverged, 0L)
  expect_identical(result$boot_reps_errors, 1L)
  expect_identical(result$algorithm,
                   "frame-invariance-complete-family-1")
  expect_no_error(.validate_frame_invariance(result, fit))
})

test_that("frame-invariance bootstrap refits the units and controls one family", {
  skip_on_cran()
  d <- simulate_efrm(n_per_group = 160, items_per_set = 6, n_sets = 1,
                     n_groups = 2, seed = 41)
  tr <- attr(d, "truth")
  f <- rasch_efrm(d, item_sets = tr$item_sets, groups = "group", id = "id",
                  boot_reps = 0)
  set.seed(42)
  before <- .Random.seed
  z <- frame_invariance(f, se_method = "bootstrap", boot_reps = 30, seed = 43)
  expect_identical(.Random.seed, before)
  expect_equal(z$boot_reps_used, 30)
  expect_equal(z$boot_reps_nonconverged, 0)
  expect_equal(z$boot_reps_errors, 0)
  expect_equal(z$boot_minimum_usable, 30)
  expect_no_error(.validate_frame_invariance(z, f))
  expect_output(print(z),
                "30/30 usable; 0 non-converged; 0 other failures",
                fixed = TRUE)
  expect_true(all(is.finite(z$locations$se)))
  p_all <- p.adjust(c(z$locations$p, z$discrimination$p), "holm")
  expect_equal(z$locations$p_adj, head(p_all, nrow(z$locations)))
  expect_equal(z$discrimination$p_adj,
               tail(p_all, nrow(z$discrimination)))
})

test_that("frame invariance requires adequate support in every frame", {
  set.seed(4401)
  n <- 60; L <- 6
  grp <- rep(c("A", "B"), each = n / 2)
  X <- matrix(rbinom(n * L, 1,
    plogis(outer(rnorm(n), seq(-1.2, 1.2, length.out = L), "-"))), n, L)
  colnames(X) <- paste0("I", seq_len(L))
  fit <- rasch_efrm(data.frame(X, group = grp),
                    item_sets = list(all = colnames(X)), groups = "group",
                    boot_reps = 0)
  expect_error(frame_invariance(fit), "at least 50 persons")
})

test_that("frame invariance does not select comparisons by calibration success", {
  skip_on_cran()
  d <- simulate_efrm(n_per_group = 120, items_per_set = 6, n_sets = 1,
                     n_groups = 3, seed = 9183)
  tr <- attr(d, "truth")
  items <- unlist(tr$item_sets, use.names = FALSE)
  g3 <- which(d$group == levels(factor(d$group))[3])
  h <- floor(length(g3) / 2)
  d[g3[seq_len(h)], items[4:6]] <- NA
  d[g3[h + seq_len(length(g3) - h)], items[1:3]] <- NA

  fit <- rasch_efrm(d, item_sets = tr$item_sets, groups = "group", id = "id",
                    boot_reps = 0)
  expect_true(fit$est$converged)
  expect_error(
    frame_invariance(fit),
    "usable separate calibration.*set1/g3.*cannot be selected")
  expect_null(.frame_invariance_conditional(fit, strict = FALSE))
})

test_that("frame items dropped by a separate calibration remain in the family", {
  d <- simulate_efrm(n_per_group = 150, items_per_set = 6, n_sets = 1,
                     n_groups = 2, n_categories = 3, seed = 771)
  tr <- attr(d, "truth")
  items <- unlist(tr$item_sets, use.names = FALSE)
  g1 <- levels(factor(d$group))[1]
  d[d$group == g1, items[1]] <- 0L
  fit <- rasch_efrm(d, item_sets = tr$item_sets, groups = "group", id = "id",
                    boot_reps = 0)
  expect_true(fit$est$converged)
  z <- frame_invariance(fit)
  expect_identical(z$excluded$item, items[1])
  expect_match(z$excluded$reason, "dropped or rescored")
  expect_equal(nrow(z$locations), 5L)
  expect_equal(z$locations$p_adj,
               p.adjust(z$locations$p, "holm", n = 6L))
  expect_identical(z$family_n, 6L)
})

test_that("a lone comparable frame item remains in the unavailable family", {
  d <- simulate_efrm(n_per_group = 180, items_per_set = 6, n_sets = 2,
                     n_groups = 2, n_categories = 3, seed = 887)
  tr <- attr(d, "truth")
  s1 <- tr$item_sets[[1]]
  focal <- d$group == levels(d$group)[2]
  for (item in s1[-1]) d[focal & d[[item]] == 1L, item] <- 2L
  fit <- rasch_efrm(d, item_sets = tr$item_sets, groups = "group", id = "id",
                    boot_reps = 0)
  z <- frame_invariance(fit)
  expect_setequal(z$excluded$item, s1)
  expect_match(z$excluded$reason[z$excluded$item == s1[1]],
               "frame origin")
  expect_identical(z$family_n, 12L)
  expect_equal(z$locations$p_adj,
               p.adjust(z$locations$p, "holm", n = 12L))
})

test_that("frame invariance refuses when no category structure is comparable", {
  set.seed(490)
  n <- 300
  group <- rep(c("A", "B"), each = n / 2)
  theta <- rnorm(n)
  X <- sapply(seq(-1, 1, length.out = 6), function(delta) {
    p <- plogis(theta - delta)
    ifelse(group == "A", rbinom(n, 1, p), rbinom(n, 2, p))
  })
  colnames(X) <- paste0("I", 1:6)
  fit <- rasch_efrm(data.frame(X, group = group),
                    item_sets = list(all = colnames(X)), groups = "group",
                    boot_reps = 0)
  expect_true(fit$est$converged)
  expect_error(frame_invariance(fit),
               "no frame pair retains at least two items")
})

test_that("drop_items preserves anchors and principal-component PCM", {
  d <- simulate_rasch(400, 8, seed = 4)
  anchors <- data.frame(item = "I01", k = 1, tau = 1)
  f <- rasch(d, id = "id", anchors = anchors, n_groups = 7,
             maxit = 80, tol = 1e-9)
  f2 <- drop_items(f, "I08")
  expect_equal(as.character(f2$est$anchors$item), "I01")
  expect_equal(f2$refit_spec$n_groups, 7)
  expect_error(drop_items(f, "I01"),
               "externally anchored item.*cannot be dropped")
  expect_error(drop_items(f, c("I01", "I08")),
               "externally anchored item.*cannot be dropped")

  p <- simulate_rasch(400, 8, model = "PCM", n_categories = 4, seed = 8)
  fp <- rasch(p, id = "id", pc_components = 2)
  fp2 <- drop_items(fp, "I08")
  expect_equal(fp2$est$n_components, 2)
  expect_equal(fp2$refit_spec$pc_components, 2)
})

test_that("drop_items preserves polytomous option scoring", {
  set.seed(31)
  d <- data.frame(
    M1 = sample(c("A", "B", "C"), 600, TRUE),
    M2 = sample(c("A", "B", "C"), 600, TRUE),
    I3 = rbinom(600, 1, 0.5), I4 = rbinom(600, 1, 0.5),
    check.names = FALSE)
  key <- data.frame(
    item = rep(c("M1", "M2"), each = 2),
    option = rep(c("A", "B"), 2),
    score = rep(c(2L, 1L), 2))
  f <- rasch(d, key = key)
  f2 <- drop_items(f, "I4")
  expect_equal(f2$mc$map, f$mc$map)
  expect_equal(f2$X[, c("M1", "M2")], f$X[, c("M1", "M2")])
  expect_true(is.data.frame(f2$refit_spec$key))
})

test_that("subtests and DIF splits retain the active Rasch specification", {
  set.seed(32)
  n <- 600
  grp <- rep(c("A", "B"), each = n / 2)
  d <- data.frame(
    M1 = sample(c("A", "B", "C"), n, TRUE),
    M2 = sample(c("A", "B", "C"), n, TRUE),
    I3 = rbinom(n, 1, .5), I4 = rbinom(n, 1, .5), grp = grp,
    check.names = FALSE)
  key <- data.frame(
    item = rep(c("M1", "M2"), each = 2),
    option = rep(c("A", "B"), 2), score = rep(c(2L, 1L), 2))
  f <- rasch(d, factors = "grp", key = key, n_groups = 7,
             maxit = 80, tol = 1e-9)

  sp <- split_items(f, "M1", by = "grp")
  expect_setequal(colnames(sp$mc$raw), c("M2", "M1 (A)", "M1 (B)"))
  expect_equal(sp$refit_spec$n_groups, 7)
  expect_equal(sp$refit_spec$maxit, 80)

  su <- combine_items(f, list(c("M1", "I3")))
  expect_identical(colnames(su$mc$raw), "M2")
  expect_true("M1+I3" %in% su$items$item)
  expect_equal(su$refit_spec$n_groups, 7)
  expect_equal(su$subtest_map[["M1+I3"]], c("M1", "I3"))
  expect_false(su$subtest_binary[["M1+I3"]])

  X <- f$X
  colnames(X)[1] <- "M1 (original)"
  fp <- rasch(data.frame(X, grp = grp, check.names = FALSE), factors = "grp")
  spp <- split_items(fp, "M1 (original)", by = "grp")
  expect_true(all(spp$split_map[grep("M1", names(spp$split_map))] ==
                    "M1 (original)"))
})

test_that("split provenance survives later ordinary structural refits", {
  d <- simulate_rasch(500, 8, n_groups = 2, seed = 351)
  f <- rasch(d, id = "id", factors = "group")
  sp <- split_items(f, "I03", by = "group")
  expect_equal(.n_unsplit_sources(.split_source_map(sp)), 7L)

  dr <- drop_items(sp, "I08")
  expect_identical(names(dr$split_map), colnames(dr$X))
  expect_true(all(dr$split_map[grep("I03", names(dr$split_map))] == "I03"))
  expect_equal(.n_unsplit_sources(.split_source_map(dr)), 6L)

  co <- combine_items(sp, c("I07", "I08"))
  expect_identical(names(co$split_map), colnames(co$X))
  expect_true(all(co$split_map[grep("I03", names(co$split_map))] == "I03"))
  expect_identical(unname(co$split_map[["I07+I08"]]), "I07+I08")
  expect_equal(.n_unsplit_sources(.split_source_map(co)), 6L)

  split_copy <- names(sp$split_map)[names(sp$split_map) != sp$split_map][1L]
  expect_false(is.na(split_copy))
  expect_error(combine_items(sp, c(split_copy, "I07")),
               "group-specific split item.*cannot be combined")

  map <- c("I1 (A)" = "I1", "I1 (B)" = "I1", I2 = "I2")
  expect_setequal(.split_source_items(c("I1 (A)", "I1 (A)", "I2"), map),
                  c("I1", "I2"))
})

test_that("structural refits do not silently transform external anchors", {
  d <- simulate_rasch(500, 8, seed = 35)
  a <- data.frame(item = "I01", k = 1, tau = 1)
  f <- rasch(d, id = "id", anchors = a)
  expect_error(split_items(f, "I01", by = rep(c("A", "B"), each = 250)),
               "anchored item cannot be split")
  expect_error(combine_items(f, list(c("I01", "I02"))),
               "anchored item cannot be combined")
  expect_equal(as.character(
    split_items(f, "I02", by = rep(c("A", "B"), each = 250))$
      est$anchors$item), "I01")
  expect_equal(as.character(combine_items(f, list(c("I02", "I03")))$
                              est$anchors$item), "I01")
})

test_that("frame resolution needs set origins as well as global unit links", {
  d <- simulate_efrm(n_per_group = 350, items_per_set = 6, n_sets = 2,
                     n_groups = 2, seed = 23)
  tr <- attr(d, "truth")
  f <- rasch_efrm(d, item_sets = tr$item_sets, groups = "group", id = "id",
                  boot_reps = 0)
  # The second set can identify phi, but cannot identify separate origins
  # for two disjoint groups of item versions in the first set. Resolving
  # every item there leaves a threshold-only flat likelihood direction.
  H <- NULL
  original <- .likelihood_curvature_ok
  testthat::local_mocked_bindings(.likelihood_curvature_ok = function(h, ...) {
    H <<- h
    original(h, ...)
  }, .package = "rasch")
  expect_error(resolve_frames(f, tr$item_sets[[1]], boot_reps = 0),
               "did not reach an identified local maximum")
  ev <- eigen(H, symmetric = TRUE, only.values = TRUE)$values
  expect_lt(min(abs(ev)) / max(abs(ev)), 1e-10)
  expect_lte(max(ev), 1e-8 * max(abs(ev)))

  # One common item retains the origin link. The second set still supplies
  # the unit link, so two common items are not required in every set.
  resolved <- resolve_frames(f, tr$item_sets[[1]][-1L], boot_reps = 0)
  expect_s3_class(resolved, "rasch_efrm")
  expect_true(resolved$est$converged)
  expect_true(all(is.finite(resolved$phi_table$phi)))
})

test_that("ETS classification uses its interval-null rule", {
  # This helper uses the supplied probabilities; public callers adjust them
  # over their planned comparison family before classification.
  expect_equal(.ets_category(0.8, 0.1, 0.01), "C+")
  expect_equal(.ets_category(0.8, 0.1, 0.20), "A")
  expect_equal(.ets_category(0.6, 0.2, 0.01), "B+")
})

test_that("DIF resolution returns its final residual-DIF table", {
  d <- simulate_rasch(500, 8, n_groups = 2, seed = 41)
  f <- rasch(d, id = "id", factors = "group")
  rr <- resolve_dif(f, max_splits = 0)
  expect_named(rr, c("algorithm", "fit", "splits", "n_splits", "stopped",
                     "dif", "notes", "effects", "n_remaining_dif",
                     "n_nonuniform"))
  expect_equal(rr$n_remaining_dif, if (is.null(rr$dif)) 0L else
    length(.split_source_items(rr$dif$item, .split_source_map(rr$fit))))
  expect_error(resolve_dif(f, min_anchors = ncol(f$X)), "min_anchors")

  bad <- f; bad$est$converged <- FALSE
  expect_error(dif_anova(bad), "did not converge")
  expect_error(dif_size(bad, "I01", "group"), "did not converge")
  expect_error(dif_contrasts(bad), "did not converge")
  expect_error(equate_tests(bad, f, independent = TRUE), "did not converge")
})

test_that("structural refits require valid item and group labels", {
  f <- rasch(simulate_rasch(240, 6, seed = 411), id = "id")
  expect_error(drop_items(f, NA_character_), "non-missing item name")
  expect_error(split_items(f, "   ", rep(c("A", "B"), each = 120)),
               "non-missing item name")
  expect_error(combine_items(f, list(c("I01", NA_character_))),
               "non-missing item names")

  g <- rep(c("A", " B "), each = 120)
  s <- split_items(f, "I01", g)
  expect_true(all(c("I01 (A)", "I01 (B)") %in% colnames(s$X)))
  expect_false(any(grepl("\\( B \\)", colnames(s$X))))
})

test_that("frame invariance compares exact observed category structures", {
  d <- simulate_efrm(n_per_group = 250, items_per_set = 6, n_sets = 1,
                     n_groups = 2, seed = 44)
  tr <- attr(d, "truth")
  f <- rasch_efrm(d, item_sets = tr$item_sets, groups = "group", id = "id",
                  boot_reps = 0)
  expect_error(testthat::with_mocked_bindings(
    frame_invariance(f),
    .item_location_covariance = function(fit) {
      V <- diag(nrow(fit$items))
      V[1L, 1L] <- -1
      V
    },
    .package = "rasch"),
    "not positive semidefinite")
  vm <- f$virtual_map
  g2 <- unique(vm$group)[2]
  item <- tr$item_sets[[1]][1]
  v <- vm$vkey[vm$item == item & vm$group == g2]
  f$X[f$X[, v] == 1, v] <- 2L
  z <- frame_invariance(f)
  expect_true(item %in% z$excluded$item)
  expect_true(all(grepl("category structure", z$excluded$reason)))
})

test_that("crossed EFRM cells remain distinct with colon-bearing levels", {
  d <- simulate_efrm(n_per_group = 220, items_per_set = 6, n_sets = 1,
                     n_groups = 2, seed = 45)
  tr <- attr(d, "truth")
  d$g1 <- ifelse(d$group == unique(d$group)[1], "a:b", "a")
  d$g2 <- ifelse(d$group == unique(d$group)[1], "c", "b:c")
  d[["g1:g2"]] <- rep(c("x", "y"), length.out = nrow(d))
  f <- rasch_efrm(d, item_sets = tr$item_sets, groups = c("g1", "g2"),
                  id = "id", factors = c("g1", "g2", "g1:g2"),
                  boot_reps = 0)
  expect_equal(nrow(f$phi_table), 2L)
  expect_false(anyDuplicated(f$phi_table$group) > 0L)
  expect_false(anyDuplicated(f$virtual_map$vkey) > 0L)
  expect_false(anyDuplicated(names(f$factors)) > 0L)
  expect_true("g1:g2" %in% names(f$factors))
  expect_true(f$frame_group[1] != "g1:g2")
})

test_that("frame resolution refuses a generated item-name collision", {
  d <- simulate_efrm(n_per_group = 180, items_per_set = 6, n_sets = 1,
                     n_groups = 2, seed = 46)
  tr <- attr(d, "truth")
  it <- tr$item_sets[[1]][1]
  other <- tr$item_sets[[1]][2]
  clash <- paste0(it, " (", levels(factor(d$group))[1], ")")
  names(d)[names(d) == other] <- clash
  sets <- tr$item_sets
  sets[[1]][sets[[1]] == other] <- clash
  f <- rasch_efrm(d, item_sets = sets, groups = "group", id = "id",
                  boot_reps = 0)
  expect_error(resolve_frames(f, it, boot_reps = 0),
               "resolved-item name already exists")
  f$est$converged <- FALSE
  expect_error(frame_invariance(f), "did not converge")
})
