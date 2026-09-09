test_that("automatic DIF follow-ups keep externally supplied factors", {
  set.seed(22)
  n <- 600L
  group <- factor(rep(c("A", "B"), each = n / 2L))
  theta <- rnorm(n)
  difficulty <- seq(-1.2, 1.2, length.out = 6L)
  shift <- matrix(0, n, length(difficulty))
  shift[group == "B", 3L] <- 1.2
  X <- matrix(rbinom(n * ncol(shift), 1,
                     plogis(outer(theta, difficulty, "-") - shift)),
              n, ncol(shift))
  colnames(X) <- paste0("I", seq_len(ncol(X)))

  # The fit has no stored factor.  The nominated factor is wholly external,
  # so automatic follow-ups must receive the normalized data frame rather
  # than trying to look its name up in fit$factors.
  fit <- rasch(X)
  supplied <- data.frame(group = group)
  auto <- dif_anova(fit, factors = supplied, n_groups = 3L, sizes = TRUE)
  direct <- dif_posthoc(fit, "I3", term = "group", factors = supplied)

  got <- auto$posthoc[auto$posthoc$item == "I3" &
                        auto$posthoc$term == "group", , drop = FALSE]
  expect_equal(nrow(got), nrow(direct$table))
  expect_equal(got$estimate, direct$table$estimate, tolerance = 1e-12)
  expect_true(any(got$significant))
  expect_identical(auto$followup_algorithm, "normalized-design-1")
})

test_that("authenticated obsolete automatic follow-ups cannot reach reports", {
  set.seed(2201)
  n <- 500L
  group <- factor(rep(c("A", "B"), each = n / 2L))
  stored_group <- factor(sample(as.character(group)))
  theta <- rnorm(n)
  difficulty <- seq(-1.2, 1.2, length.out = 6L)
  shift <- matrix(0, n, length(difficulty))
  shift[group == "B", 3L] <- 1.8
  X <- matrix(rbinom(n * ncol(shift), 1,
                     plogis(outer(theta, difficulty, "-") - shift)),
              n, ncol(shift))
  colnames(X) <- paste0("I", seq_len(ncol(X)))
  fit <- rasch(data.frame(X, group = stored_group), factors = "group")
  supplied <- data.frame(group = group)

  current <- dif_anova(fit, factors = supplied, n_groups = 3L, sizes = TRUE)
  direct <- dif_posthoc(fit, "I3", term = "group", factors = supplied)
  current_i3 <- current$posthoc[current$posthoc$item == "I3" &
                                  current$posthoc$term == "group", ,
                                drop = FALSE]
  expect_equal(current_i3$estimate, direct$table$estimate, tolerance = 1e-12)

  # Recreate the pre-fix constructor in-process: the omnibus uses `supplied`,
  # but its automatic post-hoc call resolves only that data frame's names
  # against the stale factors stored in the fit.
  posthoc_impl <- dif_posthoc
  obsolete <- testthat::with_mocked_bindings(
    dif_anova(fit, factors = supplied, n_groups = 3L, sizes = TRUE),
    dif_posthoc = function(fit, item, term = NULL, factors = NULL, ...) {
      posthoc_impl(fit, item, term = term, factors = names(factors), ...)
    },
    .package = "rasch")
  obsolete_i3 <- obsolete$posthoc[obsolete$posthoc$item == "I3" &
                                    obsolete$posthoc$term == "group", ,
                                  drop = FALSE]
  stale <- dif_posthoc(fit, "I3", term = "group", factors = "group")
  expect_equal(obsolete_i3$estimate, stale$table$estimate, tolerance = 1e-12)
  expect_gt(abs(obsolete_i3$estimate - current_i3$estimate), 0.5)

  # Refusal occurs only after the existing result hash has authenticated the
  # object. Removing the stamp without resealing is a generic integrity error;
  # an authentically shaped pre-fix result reaches the explicit compatibility
  # refusal even though its omnibus algorithm is already joint-between-1.
  tampered <- obsolete
  tampered$followup_algorithm <- NULL
  expect_error(.validate_dif_result(tampered, fit),
               "incomplete or internally inconsistent")
  obsolete <- tampered
  obsolete$result_signature <- NULL
  obsolete$result_signature <- .fit_boot_md5(unclass(obsolete))
  expect_identical(obsolete$algorithm, "joint-between-1")
  expect_false(.dif_followups_current(obsolete))
  expect_error(.validate_dif_result(obsolete, fit),
               "may use earlier factor values")

  report <- tempfile(fileext = ".html")
  on.exit(unlink(report), add = TRUE)
  expect_error(report_html(fit, report, dif = obsolete),
               "may use earlier factor values")
  expect_false(file.exists(report))

  # A pre-stamp analysis with no stored follow-up rows remains valid: there is
  # no obsolete estimate that a report could accidentally publish.
  legacy_primary <- current
  legacy_primary[c("sizes", "posthoc", "posthoc_family_n",
                   "followup_algorithm")] <- NULL
  legacy_primary$result_signature <- NULL
  legacy_primary$result_signature <- .fit_boot_md5(unclass(legacy_primary))
  expect_true(.dif_followups_current(legacy_primary))
  expect_no_error(.validate_dif_result(legacy_primary, fit))

  malformed <- legacy_primary
  malformed$sizes <- list()
  expect_false(.dif_followups_current(malformed))
})

test_that("automatic DIF follow-ups preserve replacement factors in incomplete panels", {
  set.seed(11)
  n <- 400L
  group <- factor(rep(c("A", "B"), each = n / 2L))
  old_group <- factor(sample(as.character(group)))
  id <- rep(sprintf("P%03d", seq_len(n)), 2L)
  occasion <- factor(rep(c("t1", "t2"), each = n))
  theta <- rnorm(n)
  difficulty <- seq(-1.2, 1.2, length.out = 5L)
  # Keep the occasion-specific item effect explicit rather than relying on a
  # mutable closure state; this also makes the incomplete-panel construction
  # easy to read.
  wave <- function(occasion_shift) {
    shift <- matrix(0, n, length(difficulty))
    shift[group == "B", 2L] <- 1
    shift[group == "B", 3L] <- 1.2
    shift[, 3L] <- shift[, 3L] + occasion_shift
    matrix(rbinom(n * ncol(shift), 1,
                  plogis(outer(theta, difficulty, "-") - shift)),
           n, ncol(shift))
  }
  X <- rbind(wave(0), wave(0.5))
  colnames(X) <- paste0("I", seq_len(ncol(X)))
  # Item I2 has an incomplete panel: group A is missing at t2 and group B at
  # t1.  The follow-up is therefore routed through the mixed-design branch.
  X[n + which(group == "A"), 2L] <- NA
  X[which(group == "B"), 2L] <- NA
  # The item whose follow-up is compared also has incomplete panels, but
  # retains group support at both occasions so its effect is estimable.
  X[n + seq_len(50L), 3L] <- NA
  X[n + n / 2L + seq_len(20L), 3L] <- NA

  fit <- rasch(data.frame(X, group = rep(old_group, 2L), occasion = occasion),
               id = id, factors = c("group", "occasion"))
  supplied <- data.frame(group = rep(group, 2L), occasion = occasion)
  auto <- dif_anova(fit, factors = supplied, within = "occasion",
                    n_groups = 3L, sizes = TRUE)
  direct <- dif_posthoc(fit, "I3", term = "group", factors = supplied,
                        within = "occasion")
  stale <- dif_posthoc(fit, "I3", term = "group",
                       factors = c("group", "occasion"), within = "occasion")

  got <- auto$posthoc[auto$posthoc$item == "I3" &
                        auto$posthoc$term == "group", , drop = FALSE]
  expect_equal(got$estimate, direct$table$estimate, tolerance = 1e-12)
  expect_gt(abs(got$estimate - stale$table$estimate), 0.5)
  expect_true(any(got$significant))
  expect_true(any(grepl("Incomplete panels", auto$notes, fixed = TRUE)))
})
