
test_that("the item summary reports a discrimination that orders items", {
  # the index is biased upward and its level is not interpretable, so the
  # claim under test is the ordering: a planted steep item and a planted flat
  # one should come back at the extremes of the column
  set.seed(4)
  N <- 2000; K <- 6
  th <- stats::rnorm(N); dl <- seq(-1.2, 1.2, length.out = K)
  a <- rep(1, K); a[2] <- 1.8; a[5] <- 0.5
  X <- vapply(seq_len(K), function(i)
    stats::rbinom(N, 1, stats::plogis(a[i] * (th - dl[i]))), numeric(N))
  colnames(X) <- sprintf("I%d", seq_len(K))
  f <- rasch(data.frame(id = seq_len(N), X), id = "id")

  expect_true("disc" %in% names(f$items))
  expect_true(all(is.finite(f$items$disc)))
  expect_identical(which.max(f$items$disc), 2L)
  expect_identical(which.min(f$items$disc), 5L)
})

test_that("the discrimination index covers polytomous items", {
  # the earlier implementation fitted a binomial glm and returned NA for
  # anything with more than two categories
  skip_on_cran()
  set.seed(9)
  N <- 2000; K <- 5
  th <- stats::rnorm(N, 0, 1.3)
  tau <- lapply(seq_len(K), function(i)
    c(-1, 0, 1) + seq(-0.6, 0.6, length.out = K)[i])
  a <- rep(1, K); a[2] <- 1.9; a[4] <- 0.5
  X <- vapply(seq_len(K), function(i)
    vapply(th, function(t)
      sample(0:3, 1, prob = item_moments(t, tau[[i]], disc = a[i])$P), 0),
    numeric(N))
  colnames(X) <- sprintf("I%d", seq_len(K))
  f <- rasch(data.frame(id = seq_len(N), X), id = "id")

  expect_true(all(f$items$max == 3))
  expect_true(all(is.finite(f$items$disc)))
  expect_identical(which.max(f$items$disc), 2L)
  expect_identical(which.min(f$items$disc), 4L)
})

test_that("an item with no variation has no discrimination", {
  # a wholly constant item never reaches the item table -- the fit drops it
  # and says so -- so the guard is exercised where it actually fires: an item
  # that varies overall but not among the persons the index is computed on
  d <- simulate_rasch(300, 6, seed = 2)
  d$I02 <- 1L
  f <- suppressWarnings(rasch(d, id = "id"))
  expect_false("I02" %in% f$items$item)
  expect_match(paste(f$notes, collapse = " "), "constant item")

  tau <- list(0, 0.5)
  X <- cbind(c(0, 1, 0, 1, 1, 0), c(1, 1, 1, 1, 1, 1))
  disc <- .item_discrim(seq(-1, 1, length.out = 6), X, tau,
                        extreme = rep(FALSE, 6))
  expect_true(is.finite(disc[1]))
  expect_true(is.na(disc[2]))
})

test_that("frame_invariance results print without scientific notation", {
  skip_on_cran()
  # every other analysis entry point tags its tables; this one did not, so a
  # user reading inv$locations saw the exponents print.rasch_table removes
  d <- simulate_efrm(n_per_group = 300, items_per_set = 6, n_sets = 2,
                     n_groups = 2, set_unit_ratio = 1.3, seed = 7)
  tr <- attr(d, "truth")
  f <- rasch_efrm(d, item_sets = tr$item_sets, groups = "group", id = "id",
                  boot_reps = 0)
  inv <- frame_invariance(f)
  expect_s3_class(inv$locations, "rasch_table")
  expect_s3_class(inv$discrimination, "rasch_table")
  expect_s3_class(inv$summary, "rasch_table")
  expect_false(any(grepl("e[-+][0-9]", capture.output(print(inv$locations)))))
  # the returned elements match what the documentation claims
  expect_setequal(names(inv), c("locations", "discrimination", "summary",
                                "excluded", "algorithm", "alpha", "adjust",
                                "se_method",
                                "family_n", "boot_reps", "boot_reps_used",
                                "boot_reps_nonconverged", "boot_reps_errors",
                                "boot_minimum_usable", "bootstrap_stratified",
                                "seed", "fit_signature", "result_signature"))
})

test_that("a frame model's virtual item name resolves to its source item", {
  skip_on_cran()
  # the application selects an item from fit$items, which for a frame model
  # names it by the frame that took it; drop_items() works on the source item
  d <- simulate_efrm(n_per_group = 300, items_per_set = 6, n_sets = 2,
                     n_groups = 2, set_unit_ratio = 1.3, seed = 7)
  tr <- attr(d, "truth")
  f <- rasch_efrm(d, item_sets = tr$item_sets, groups = "group", id = "id",
                  boot_reps = 0)
  virtual <- f$items$item[1]
  src <- names(f$set_of)
  expect_false(virtual %in% src)               # the mismatch the app hit
  bare <- sub(":[^:]*$", "", virtual)
  expect_true(bare %in% src)
  expect_error(drop_items(f, virtual), "not in the fit")
  expect_s3_class(drop_items(f, bare, boot_reps = 0), "rasch_efrm")
})

test_that("the report and saved outputs carry the invariance test", {
  # a report that shows frame units without the test of the assumption behind
  # them invites the reader to trust them further than the analysis warrants
  skip_on_cran()
  set.seed(11)
  N <- 500; K <- 8
  phi <- c(0.85, 1.18); delta <- seq(-1.5, 1.5, length.out = K)
  mk <- function(g, dsc, shift) {
    th <- stats::rnorm(N, 0, 1.3)
    X <- vapply(seq_len(K), function(i)
      stats::rbinom(N, 1, stats::plogis(
        phi[g] * dsc[i] * (th - delta[i] - shift[i]))), numeric(N))
    colnames(X) <- sprintf("I%02d", seq_len(K)); X
  }
  dsc <- rep(1, K); dsc[3] <- 1.8
  sh <- rep(0, K); sh[6] <- 0.6
  d <- data.frame(id = sprintf("P%04d", seq_len(2 * N)),
                  rbind(mk(1, rep(1, K), rep(0, K)), mk(2, dsc, sh)),
                  group = rep(c("g1", "g2"), each = N), check.names = FALSE)
  f <- rasch_efrm(d, item_sets = list(set1 = sprintf("I%02d", seq_len(K))),
                  groups = "group", id = "id", boot_reps = 0)

  html <- tempfile(fileext = ".html")
  on.exit(unlink(html), add = TRUE)
  report_warnings <- character(0)
  withCallingHandlers(
    report_html(f, html),
    warning = function(w) {
      report_warnings <<- c(report_warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  expect_gt(length(report_warnings), 0L)
  expect_true(all(grepl("residual PCA is undefined", report_warnings)))
  x <- readLines(html, warn = FALSE)
  expect_true(any(grepl("Item invariance across frames", x)))
  # the section reports whichever way the test came out; which branch fires
  # is a question about power, tested elsewhere, not about the report
  expect_true(any(grepl("Locations differing across frames", x)) ||
                any(grepl("No available item-location comparison differs", x)))
  expect_true(any(grepl("Descriptive discrimination comparisons", x)))
  expect_true(any(grepl("bootstrap", x, ignore.case = TRUE)))
  expect_true(any(grepl("<th>rmsd</th>", x, fixed = TRUE))) # summary is there
  expect_false(any(grepl(">[0-9.]+e[-+][0-9]+<", x)))

  dir <- tempfile(); dir.create(dir)
  on.exit(unlink(dir, recursive = TRUE), add = TRUE)
  save_warnings <- character(0)
  withCallingHandlers(
    save_outputs(f, dir, formats = "png", item_plots = FALSE),
    warning = function(w) {
      save_warnings <<- c(save_warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  expect_gt(length(save_warnings), 0L)
  expect_true(all(grepl("residual PCA is undefined", save_warnings)))
  saved <- list.files(dir, pattern = "invariance", recursive = TRUE,
                      full.names = TRUE)
  expect_length(saved, 3L)
  loc <- utils::read.csv(grep("locations", saved, value = TRUE))
  # the full table, not the curated screen version
  expect_true(all(c("location_1", "location_2") %in% names(loc)))
})
