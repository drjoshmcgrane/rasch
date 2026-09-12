.partial_information_fixture <- function(kind) {
  m <- c(1L, 2L, 3L, 1L)
  vm <- data.frame(vkey = paste0("v", 1:4), item = paste0("I", 1:4),
                   set = rep("core", 4), group = rep("G", 4),
                   rater = rep("A", 4))
  structure(list(est = list(converged = TRUE),
    alpha_table = data.frame(set = "core", alpha = 1), m = m,
    tau_list = lapply(m, function(k) rep(0, k)),
    disc = c(0.5, 1, 1.5, 2), virtual_map = vm, facet_spec = "rater",
    items = data.frame(item = vm$vkey),
    X = matrix(c(0, 1, NA, NA, NA, 1, 2, NA, 1, 0, 2, 1,
                 NA, NA, NA, NA, 1, 2, NA, NA), 5, byrow = TRUE,
               dimnames = list(NULL, vm$vkey))),
    class = c(paste0("rasch_", kind), "rasch"))
}

test_that("structural information uses exact observed item patterns", {
  for (kind in c("efrm", "mfrm")) {
    fit <- .partial_information_fixture(kind)
    blocks <- .design_blocks(fit)
    expect_setequal(unname(blocks), list(1:2, 2:3, 1:4))
    expect_length(unique(names(blocks)), 3L)
    expect_error(test_information(fit, grid = 0, items = 1e100), "lie in")
    expect_error(test_information(fit, grid = 0, items = 1 + 1i), "whole numbers")
    ti <- test_information(fit, grid = 0)
    # Equal thresholds at zero give a uniform score on 0:m, with variance
    # m(m+2)/12. This reference does not call the package's moment code.
    contribution <- fit$disc^2 * fit$m * (fit$m + 2) / 12
    expected <- vapply(blocks, function(ii) sum(contribution[ii]), 0)
    expect_equal(ti$info, unname(expected), tolerance = 1e-12)
    expect_equal(ti$sem, unname(1 / sqrt(expected)), tolerance = 1e-12)
    subset <- test_information(fit, grid = 0, items = c("v1", "v2"))
    # blocks that differ only outside the selection are one design once
    # restricted to it, so the curve is reported once
    restricted <- unique(lapply(blocks, intersect, 1:2))
    expect_equal(subset$info, vapply(restricted, function(ii)
      sum(contribution[ii]), 0), tolerance = 1e-12)

    # The same correction must reach the plotted expected totals.
    curves <- list()
    grDevices::pdf(NULL)
    tryCatch(with_mocked_bindings(plot_tcc(fit, grid = c(0, 0.1)),
      lines = function(x, y, ...) curves[[length(curves) + 1L]] <<- y,
      .package = "rasch"), finally = grDevices::dev.off())
    expect_equal(vapply(curves, `[`, 0, 1L),
      unname(vapply(blocks, function(ii) sum(fit$m[ii]) / 2, 0)),
      tolerance = 1e-12)
  }
})

test_that("partial administration labels cannot merge distinct designs", {
  for (kind in c("efrm", "mfrm")) {
    fit <- .partial_information_fixture(kind)
    fit$virtual_map$item <- c("a+b", "c", "a", "b+c")
    fit$X <- rbind(c(0, 1, NA, NA), c(NA, NA, 2, 1))
    blocks <- .design_blocks(fit)
    expect_setequal(unname(blocks), list(1:2, 3:4))
    expect_length(unique(names(blocks)), 2L)
    expect_length(unique(test_information(fit, grid = 0)$design), 2L)
  }
})

test_that("EFRM information agrees with fitted partial-pattern score curves", {
  d <- simulate_efrm(n_per_group = 120, items_per_set = 6, n_sets = 1,
                     n_groups = 2, seed = 91731)
  sets <- attr(d, "truth")$item_sets
  for (g in unique(d$group))
    d[head(which(d$group == g), 30), sets[[1]][1:2]] <- NA
  fit <- rasch_efrm(d, item_sets = sets, groups = "group", id = "id",
                    boot_reps = 0, workers = 1)
  ti <- test_information(fit, grid = 0)
  sc <- fit$score_curves[abs(fit$score_curves$theta) < 1e-8, ]
  expect_equal(nrow(ti), 4L)
  expect_equal(sort(ti$sem), sort(sc$sem), tolerance = 1e-8)
  # One enumeration and one labeller: a score curve, its information and the
  # curve plots cannot name the same administration differently.
  expect_identical(sort(sc$design), sort(ti$design))
  expect_identical(sort(unique(fit$score_curves$design)),
                   sort(names(.design_blocks(fit))))
  # A label identifies one administration: as many distinct labels as there
  # are distinct observed patterns, counted without the shared enumeration.
  expect_identical(length(unique(sc$design)),
                   nrow(unique(!is.na(fit$X))))
  expect_identical(sum(sc$n_persons), nrow(fit$X))
})

test_that("person-item information excludes designs outside the selected subgroup", {
  for (kind in c("efrm", "mfrm")) {
    fit <- .partial_information_fixture(kind)
    fit$person <- data.frame(theta = seq(-1, 1, length.out = 5))
    fit$factors <- data.frame(subgroup = c("short", "other", "full", NA, "short"),
                              frame = rep("G", 5))
    fit$frame_group <- "frame"
    fit$thresholds <- data.frame(item = rep(1:4, fit$m),
                                 tau = unlist(fit$tau_list))
    original <- fit
    captured <- NULL
    information_fun <- test_information
    grDevices::pdf(NULL)
    tryCatch(with_mocked_bindings(
      plot_pimap(fit, group = "subgroup: short", information = TRUE),
      test_information = function(fit, grid, items = NULL) {
        captured <<- information_fun(fit, grid, items)
        captured
      }, .package = "rasch"), finally = grDevices::dev.off())
    expect_length(unique(captured$design), 1L)
    expected <- vapply(captured$theta, function(th)
      sum(vapply(1:2, function(i) fit$disc[i]^2 *
        item_moments(th, fit$tau_list[[i]], fit$disc[i])$V, 0)), 0)
    expect_equal(captured$info, expected, tolerance = 1e-12)
    expect_identical(fit, original)
  }
})

# A fit whose only variation is which cells were answered.
.nonresponse_fixture <- function(kind, X) {
  k <- ncol(X)
  vm <- data.frame(vkey = colnames(X), item = paste0("I", seq_len(k)),
                   set = rep("core", k), group = rep("G", k),
                   rater = rep("A", k), stringsAsFactors = FALSE)
  structure(list(est = list(converged = TRUE),
    alpha_table = data.frame(set = "core", alpha = 1),
    m = rep(1L, k), tau_list = rep(list(0), k), disc = rep(1, k),
    virtual_map = vm, facet_spec = "rater",
    items = data.frame(item = vm$vkey), X = X),
    class = c(paste0("rasch_", kind), "rasch"))
}

test_that("item nonresponse leaves each person the items they answered", {
  gaps <- utils::combn(10L, 2L)[, 1:40]     # 40 persons, 40 distinct gaps
  X <- matrix(1L, 40L, 10L, dimnames = list(NULL, paste0("v", 1:10)))
  X[cbind(1:40, gaps[1, ])] <- NA
  X[cbind(1:40, gaps[2, ])] <- NA
  for (kind in c("efrm", "mfrm")) {
    fit <- .nonresponse_fixture(kind, X)
    # nobody answered all ten, so no ten-item curve may be drawn, and the
    # eight a person did answer are all their own curve may claim
    expect_no_warning(blocks <- .design_blocks(fit))
    expect_length(blocks, 40L)
    expect_true(all(lengths(blocks) == 8L))
    ti <- test_information(fit, grid = 0)
    expect_equal(nrow(ti), 40L)
    expect_equal(max(ti$info), 8 * 0.25, tolerance = 1e-12)
  }
  # twenty persons who did answer all ten add that design and change no other
  X <- rbind(X, matrix(1L, 20L, 10L, dimnames = list(NULL, colnames(X))))
  for (kind in c("efrm", "mfrm")) {
    fit <- .nonresponse_fixture(kind, X)
    expect_no_warning(blocks <- .design_blocks(fit))
    expect_length(blocks, 41L)
    expect_equal(sum(lengths(blocks) == 10L), 1L)
  }
})

test_that("a design taken by few persons is not absorbed into a fuller one", {
  X <- matrix(1L, 216L, 10L, dimnames = list(NULL, paste0("v", 1:10)))
  X[201:206, 6:10] <- NA                    # a short form, given to 6 of 216
  X[cbind(207:216, 1:10)] <- NA             # ten persons, one stray gap each
  for (kind in c("efrm", "mfrm")) {
    fit <- .nonresponse_fixture(kind, X)
    blocks <- .design_blocks(fit)
    # the short form keeps its own curve however few persons took it, and so
    # does every one-gap pattern: 1:10, 1:5, and ten of nine items
    expect_length(blocks, 12L)
    expect_true(any(vapply(blocks, function(ii) setequal(ii, 1:5), NA)))
    ti <- test_information(fit, grid = 0)
    # the short form is 5 of the 10 items, so it cannot carry the whole
    # test's information
    expect_equal(min(ti$info), max(ti$info) / 2, tolerance = 1e-12)
  }
})

test_that("a gap is never filled in from the form it sits in", {
  X <- matrix(NA_integer_, 75L, 10L, dimnames = list(NULL, paste0("v", 1:10)))
  X[1:42, 1:6] <- 1L                        # form A, 30 persons plus 12 with
  X[cbind(31:42, rep(1:6, 2))] <- NA        # a gap inside it
  X[43:75, 5:10] <- 1L                      # form B, 30 persons plus 3 with
  X[cbind(73:75, 5:7)] <- NA                # a gap inside it
  obs <- !is.na(X)
  for (kind in c("efrm", "mfrm")) {
    fit <- .nonresponse_fixture(kind, X)
    blocks <- .design_blocks(fit)
    expect_true(all(vapply(blocks, function(ii)
      any(rowSums(obs[, ii, drop = FALSE]) == length(ii)), NA)))
    ti <- test_information(fit, grid = 0)
    expect_equal(max(ti$info), 6 * 0.25, tolerance = 1e-12)
  }
})

test_that("how a pattern is read does not depend on the group size", {
  gaps <- utils::combn(8L, 2L)
  for (n in c(15L, 60L)) {
    X <- matrix(1L, n, 8L, dimnames = list(NULL, paste0("v", 1:8)))
    take <- gaps[, rep_len(seq_len(ncol(gaps)), n)]
    X[cbind(seq_len(n), take[1, ])] <- NA
    X[cbind(seq_len(n), take[2, ])] <- NA
    # the same gaps, left by 1 person or shared by 2 or 3 of them, are read
    # the same way: one curve over the six items those persons answered
    seen <- nrow(unique(!is.na(X)))
    for (kind in c("efrm", "mfrm")) {
      fit <- .nonresponse_fixture(kind, X)
      expect_no_warning(blocks <- .design_blocks(fit))
      expect_length(blocks, seen)
      expect_true(all(lengths(blocks) == 6L))
    }
  }
})

test_that("a restricted design is labelled by the items that produced it", {
  for (kind in c("efrm", "mfrm")) {
    fit <- .partial_information_fixture(kind)
    ti <- test_information(fit, grid = 0, items = c("v1", "v2"))
    expect_false(any(grepl("I3", ti$design) | grepl("I4", ti$design)))
    expect_true(all(grepl("I2", ti$design)))
    # 1:2 and 1:4 are the same design once restricted to the selection
    expect_equal(nrow(ti), 2L)
    expect_equal(anyDuplicated(ti$design), 0L)
  }
})

test_that("a design block is never wider than someone's administration", {
  set.seed(4)
  X <- matrix(NA_integer_, 200L, 30L, dimnames = list(NULL, paste0("v", 1:30)))
  X[, 1:10] <- 1L                           # core, taken by everyone
  half <- rep(1:2, length.out = 200L)
  X[half == 1, 11:20] <- 1L                 # one extension, to half of them
  X[half == 2, 21:30] <- 1L                 # the other, to the other half
  X[!is.na(X) & runif(length(X)) < 0.3] <- NA
  obs <- !is.na(X)
  for (kind in c("efrm", "mfrm")) {
    fit <- .nonresponse_fixture(kind, X)
    blocks <- .design_blocks(fit)
    # at this much nonresponse hardly anyone answers a whole form, so most
    # patterns are nobody but their own person's -- and none of them may be
    # merged with another, which would sum two mutually exclusive extensions
    # into a form nobody sat
    expect_true(all(vapply(blocks, function(ii)
      any(rowSums(obs[, ii, drop = FALSE]) == length(ii)), NA)))
    expect_false(any(vapply(blocks, function(ii)
      any(ii %in% 11:20) && any(ii %in% 21:30), NA)))
    ti <- test_information(fit, grid = 0)
    expect_equal(max(ti$info), max(rowSums(obs)) / 4, tolerance = 1e-12)
  }
})

test_that("mutually exclusive forms are never summed into one design", {
  # three 10-item forms of 90 persons; 60 of each form leave one item of
  # their own form unanswered, so most patterns are shared by six persons
  X <- matrix(NA_integer_, 270L, 30L, dimnames = list(NULL, paste0("v", 1:30)))
  for (f in 0:2) {
    rows <- f * 90L + seq_len(90L)
    X[rows, f * 10L + 1:10] <- 1L
    X[cbind(rows[31:90], f * 10L + rep(1:10, each = 6L))] <- NA
  }
  owner <- rep(1:3, each = 10L)
  for (kind in c("efrm", "mfrm")) {
    fit <- .nonresponse_fixture(kind, X)
    blocks <- .design_blocks(fit)
    expect_length(blocks, 33L)              # 3 whole forms, 30 with a gap
    expect_true(all(vapply(blocks, function(ii)
      length(unique(owner[ii])) == 1L, NA)))
    ti <- test_information(fit, grid = 0)
    expect_equal(max(ti$info), 10 * 0.25, tolerance = 1e-12)
  }
})
