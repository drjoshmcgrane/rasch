make_dif_boot_fit <- function(seed = 5101, N = 140L, L = 5L,
                              polytomous = FALSE) {
  set.seed(seed)
  group <- factor(rep(c("A", "B"), length.out = N))
  theta <- rnorm(N) + ifelse(group == "B", 0.25, 0)
  if (!polytomous) {
    difficulty <- seq(-1.2, 1.2, length.out = L)
    X <- matrix(rbinom(N * L, 1, plogis(outer(theta, difficulty, "-"))),
                N, L)
  } else {
    tau <- seq(-1, 1, length.out = 3L)
    X <- vapply(seq_len(L), function(j) {
      eta <- vapply(0:3, function(k)
        k * (theta - seq(-.8, .8, length.out = L)[j]) -
          c(0, cumsum(tau))[k + 1L], numeric(N))
      apply(eta, 1L, function(z) sample.int(4L, 1L,
        prob = exp(z - max(z)))) - 1L
    }, integer(N))
  }
  colnames(X) <- paste0("I", seq_len(L))
  rasch(X, model = "PCM", factors = data.frame(group = group), n_groups = 2)
}

test_that("conditional DIF bootstrap repeats the complete declared family", {
  fit <- make_dif_boot_fit()
  da <- dif_anova(fit, n_groups = 2, p_adjust = "holm", alpha = .1)
  set.seed(77)
  before <- .Random.seed
  db <- suppressWarnings(dif_bootstrap(fit, da, B = 3, workers = 1,
                                       seed = 912))

  expect_s3_class(db, "rasch_dif_bootstrap")
  expect_identical(.Random.seed, before)
  expect_equal(db$B_used, 3L)
  expect_equal(db$family_n,
    sum(!da$term_ids %in% c("Residuals", "ci")))
  expect_identical(dim(db$replicates$F), c(3L, db$family_n))
  expect_identical(dim(db$replicates$p), c(3L, db$family_n))
  expect_identical(db$algorithm, "reference-minp-1")
  family <- !da$term_ids %in% c("Residuals", "ci")
  expected_marginal <- vapply(seq_len(db$family_n), function(j)
    (1 + sum(db$replicates$p[, j] <= da$terms$p[family][j])) /
      (db$B_used + 1), numeric(1))
  expect_equal(db$terms$p_boot[family], expected_marginal)
  expect_true(all(is.finite(db$terms$p_boot[db$terms$n_boot > 0L])))
  expect_true(all(db$terms$p_boot_adj[db$terms$n_boot > 0L] >=
                  db$terms$p_boot[db$terms$n_boot > 0L]))
  expect_true(all(db$summary$n_boot_uniform == 3L))
  expect_true(all(db$summary$n_boot_nonuniform == 3L))
  expect_identical(db$alpha, .1)
  expect_match(paste(capture.output(print(db)), collapse = "\n"),
               "Global-null familywise bootstrap flags", fixed = TRUE)
})

test_that("conditional DIF bootstrap is seed reproducible and term-safe", {
  fit <- make_dif_boot_fit(seed = 5102)
  names(fit$factors) <- "ci"
  da <- dif_anova(fit, factors = "ci", n_groups = 2)
  a <- suppressWarnings(dif_bootstrap(fit, da, B = 2, workers = 1, seed = 5))
  b <- suppressWarnings(dif_bootstrap(fit, da, B = 2, workers = 1, seed = 5))
  expect_identical(a$replicates, b$replicates)
  expect_true(all(is.finite(a$summary$p_uniform_boot_adj)))
  expect_true(all(is.finite(a$summary$p_nonuniform_boot_adj)))
  expect_no_error(rasch:::.validate_dif_bootstrap(a, fit, da))

  changed <- fit
  changed$person$id[1] <- "another person"
  expect_error(rasch:::.validate_dif_bootstrap(a, changed),
               "different fitted model")
})

test_that("conditional generation preserves score and missingness in PCM data", {
  fit <- make_dif_boot_fit(seed = 5103, N = 120, L = 4,
                           polytomous = TRUE)
  fit$X[seq(1, nrow(fit$X), by = 7), 2] <- NA
  fit <- rasch(fit$X, model = "PCM", factors = fit$factors, n_groups = 2)
  set.seed(3)
  xb <- rasch:::.fit_gen_conditional(
    fit$X, fit$tau_list, is.na(fit$X))
  expect_identical(is.na(xb), is.na(fit$X))
  expect_equal(rowSums(xb, na.rm = TRUE), rowSums(fit$X, na.rm = TRUE))
  da <- dif_anova(fit, n_groups = 2)
  expect_s3_class(suppressWarnings(
    dif_bootstrap(fit, da, B = 1, workers = 1, seed = 4)),
    "rasch_dif_bootstrap")
})

test_that("conditional DIF bootstrap validates scope and provenance", {
  fit <- make_dif_boot_fit(seed = 5104)
  da <- dif_anova(fit, n_groups = 2)
  old <- da
  old$bootstrap_design <- NULL
  expect_error(dif_bootstrap(fit, old, B = 1, workers = 1),
               "predates conditional-bootstrap")
  expect_error(dif_bootstrap(fit, da, B = 0), "whole number")
  expect_error(dif_bootstrap(fit, da, workers = 0), "whole number")
  expect_error(dif_bootstrap(fit, da, seed = -1), "whole number")

  expect_error(dif_bootstrap(list(), B = 1), "must be a fitted model")
})

test_that("explanatory and Multiple Ratings DIF retain score sufficiency", {
  set.seed(5110)
  N <- 120L
  group <- factor(rep(c("A", "B"), each = N / 2L))
  theta <- rnorm(N)
  X <- sapply(seq(-1, 1, length.out = 5L), function(d)
    rbinom(N, 1L, plogis(theta - d)))
  colnames(X) <- paste0("I", seq_len(ncol(X)))
  q <- data.frame(item = colnames(X), feature = seq_len(ncol(X)))
  ef <- rasch_explanatory(X, q, ~ feature,
                          factors = data.frame(group = group), n_groups = 2)
  ed <- dif_anova(ef, n_groups = 2)
  eb <- suppressWarnings(dif_bootstrap(ef, ed, B = 1, workers = 1,
                                       seed = 5111))
  expect_identical(eb$model_kind, "explanatory")
  expect_identical(eb$null_method, "person-score conditional")
  relaxed <- relax_explanatory(ef, "I3", "location")
  rd <- dif_anova(relaxed, n_groups = 2)
  rb <- suppressWarnings(dif_bootstrap(relaxed, rd, B = 1, workers = 1,
                                       seed = 5119))
  expect_equal(nrow(relaxed$explanatory$relaxations), 1L)
  expect_equal(rb$B_used, 1L)

  set.seed(5128)
  group_p <- factor(rep(c("A", "B"), each = N / 2L))
  x <- seq(-1, 1, length.out = 5L)
  theta_p <- rnorm(N)
  Xp <- vapply(x, function(xi) vapply(theta_p, function(th)
    sample.int(3L, 1L,
      prob = item_moments(th, c(-.7, .7) + .3 * xi)$P) - 1L,
    integer(1)), integer(N))
  colnames(Xp) <- paste0("P", seq_len(ncol(Xp)))
  ep <- rasch_explanatory(
    Xp, data.frame(item = colnames(Xp), x = x), ~ x + threshold,
    factors = data.frame(group = group_p), n_groups = 2)
  epd <- dif_anova(ep, n_groups = 2)
  epb <- suppressWarnings(dif_bootstrap(ep, epd, B = 1, workers = 1,
                                        seed = 5129))
  expect_identical(ep$explanatory_model, "LPCM")
  expect_equal(epb$B_used, 1L)

  d <- simulate_mfrm(n_persons = 90, n_items = 4, n_raters = 3,
                     n_categories = 3,
                     interaction = list(rater = "R2", item = "I2", bias = .4),
                     seed = 5112)
  ids <- unique(d$person)
  d$group <- rep(c("A", "B"), length.out = length(ids))[
    match(d$person, ids)]
  mf <- rasch_mfrm(d, "person", "item", "score", facets = "rater",
                   interaction = "rater", factors = "group", n_groups = 2)
  md <- dif_anova(mf, n_groups = 2)
  mb <- suppressWarnings(dif_bootstrap(mf, md, B = 1, workers = 1,
                                       seed = 5113))
  expect_identical(mb$model_kind, "mfrm")
  expect_identical(mb$null_method, "person-score conditional")
  expect_identical(mf$interaction, "rater")
  expect_equal(mb$B_used, 1L)
})

test_that("Multiple Ratings bootstrap keeps internal roles distinct from factors", {
  d <- simulate_mfrm(n_persons = 90, n_items = 4, n_raters = 3,
                     n_categories = 3, seed = 5126)
  ids <- unique(d$person)
  d$..person <- rep(c("A", "B"), length.out = length(ids))[
    match(d$person, ids)]
  fit <- rasch_mfrm(d, "person", "item", "score", facets = "rater",
                    factors = "..person", n_groups = 2)
  da <- dif_anova(fit, n_groups = 2)
  db <- suppressWarnings(dif_bootstrap(fit, da, B = 1, workers = 1,
                                       seed = 5127))
  expect_equal(db$B_used, 1L)
})

test_that("Extended Frames DIF conditions within item sets", {
  d <- simulate_efrm(n_per_group = 90, items_per_set = 4, n_sets = 2,
                     n_groups = 2, seed = 5114)
  tr <- attr(d, "truth")
  d$site <- rep(c("A", "B"), length.out = nrow(d))
  fit <- rasch_efrm(d, item_sets = tr$item_sets, groups = "group", id = "id",
                    factors = "site", boot_reps = 0, n_groups = 2)
  set.seed(5115)
  xb <- rasch:::.dif_gen_efrm(fit)
  block <- interaction(fit$virtual_map$set, fit$virtual_map$group,
                       drop = TRUE)
  for (cols in split(seq_len(ncol(fit$X)), block))
    expect_equal(rowSums(xb[, cols, drop = FALSE], na.rm = TRUE),
                 rowSums(fit$X[, cols, drop = FALSE], na.rm = TRUE))
  da <- dif_anova(fit, n_groups = 2)
  db <- suppressWarnings(dif_bootstrap(fit, da, B = 1, workers = 1,
                                       seed = 1))
  expect_identical(db$model_kind, "efrm")
  expect_identical(db$null_method, "item-set-score conditional")
})

test_that("Extended Frames DIF retains crossed frame cells", {
  d <- simulate_efrm(n_per_group = 140, items_per_set = 6, n_sets = 1,
                     n_groups = 2, seed = 5120)
  tr <- attr(d, "truth")
  names(d)[names(d) == "group"] <- "efrm_factor_0002"
  d[["school year"]] <- rep(c("N", "S"), length.out = nrow(d))
  d$site <- rep(rep(c("A", "B"), each = 2L), length.out = nrow(d))
  fit <- rasch_efrm(
    d, item_sets = tr$item_sets,
    groups = c("efrm_factor_0002", "school year"), id = "id",
    factors = "site", boot_reps = 0, n_groups = 2)
  da <- dif_anova(fit, n_groups = 2)
  db <- suppressWarnings(dif_bootstrap(fit, da, B = 1, workers = 1,
                                       seed = 5121))
  expect_equal(nrow(fit$phi_table), 4L)
  expect_true(all(c("efrm_factor_0002", "school year") %in%
                    names(fit$factors)))
  expect_false("school.year" %in% names(fit$factors))
  expect_setequal(fit$frame_group,
                  c("efrm_factor_0002:school year",
                    "efrm_factor_0002", "school year"))
  expect_identical(da$factor_names, "site")
  expect_true(all(c("efrm_factor_0002", "school year",
                    "efrm_factor_0002:school year") %in%
                  fit$phi_factorial_tests$term))
  expect_equal(db$B_used, 1L)
})

test_that("a structural one-level frame does not block another DIF factor", {
  set.seed(5122)
  n <- 120L
  X <- matrix(rbinom(n * 6L, 1L, 0.5), n, 6L,
              dimnames = list(NULL, paste0("I", 1:6)))
  d <- data.frame(id = sprintf("P%03d", seq_len(n)), X,
                  frame = "all",
                  site = rep(c("A", "B"), each = n / 2L),
                  check.names = FALSE)
  fit <- rasch_efrm(
    d, item_sets = list(core = colnames(X)), groups = "frame", id = "id",
    factors = "site", boot_reps = 0L, n_groups = 2L)
  expect_no_error(da <- dif_anova(fit, n_groups = 2L))
  expect_identical(da$factor_names, "site")
  expect_match(paste(da$notes, collapse = " "),
               "frame factor.*excluded")
})

test_that("Comparative Judgement DIF uses the fitted-outcome null", {
  d <- simulate_btl(n_objects = 5, n_judges = 20, reps_per_pair = 40,
                    seed = 5117)
  fit <- btl(d, "object_a", "object_b", winner = "winner", judge = "judge")
  source_settings <- list(model_type = "btl", bt_a = "object_a",
                           bt_b = "object_b", bt_win = "winner",
                           bt_judge = "judge")
  attr(fit, "rasch_app_source") <- list(
    data = as.data.frame(d), settings = source_settings,
    resources = list(), simulation = list())
  judges <- unique(fit$comparisons$judge)
  group <- setNames(rep(c("A", "B"), length.out = length(judges)), judges)
  da <- btl_dif(fit, group, objects = "O3", min_n = 10)
  db <- suppressWarnings(dif_bootstrap(fit, da, B = 1, workers = 1,
                                       seed = 5118))
  expect_identical(db$model_kind, "btl")
  expect_identical(db$null_method, "fitted-outcome")
  expect_true(all(c("object", "p_uniform_boot_adj",
                    "p_nonuniform_boot_adj") %in% names(db$summary)))
  expect_no_error(rasch:::.validate_dif_bootstrap(db, fit, da))
  changed_da <- da
  changed_da$summary$p_uniform_adj[1L] <- .777
  expect_error(rasch:::.validate_btl_dif_result(changed_da, fit),
               "incomplete or internally inconsistent")
  expect_error(dif_bootstrap(fit, B = 1), "explicit btl_dif")

  project <- .seal_app_project(list(
    format = "rasch-shiny-project", schema = 2L,
    data = d, model_type = "btl", base_fit = fit,
    rasch_steps = list(), btl_steps = list(), settings = source_settings,
    results = list(btl_dif = da,
      btl_dif_meta = list(judge_col = "judge", fitted_judge_col = "judge"),
      dif_bootstrap = list(db = db, B = 1L, seed = 5118L, kind = "bdif"))))
  expect_no_error(.validate_app_project(project))

  d2 <- d
  d2$winner[seq_len(20L)] <- ifelse(
    d2$winner[seq_len(20L)] == d2$object_a[seq_len(20L)],
    d2$object_b[seq_len(20L)], d2$object_a[seq_len(20L)])
  fit2 <- btl(d2, "object_a", "object_b", winner = "winner",
              judge = "judge")
  primary_only <- .seal_app_project(list(
    format = "rasch-shiny-project", schema = 2L,
    data = d2, model_type = "btl", base_fit = fit2,
    rasch_steps = list(), btl_steps = list(), settings = list(),
    results = list(btl_dif = da)))
  expect_error(.validate_app_project(primary_only), "primary DIF analysis")

  project$results$btl_dif <- btl_dif(
    fit, group, objects = "O3", min_n = 10, alpha = .1)
  project <- .seal_app_project(project)
  expect_error(.validate_app_project(project), "different DIF analysis")

  boundary <- fit
  boundary$objects$extreme[1] <- TRUE
  expect_error(dif_bootstrap(boundary, da, B = 1),
               "undefeated or winless")
})

test_that("Comparative Judgement DIF retains category and history models", {
  d <- simulate_btl(5, 24, reps_per_pair = 12, model = "polytomous",
                    n_categories = 4, seed = 5122)
  d$response <- ordered(d$response, levels = 0:3)
  fit <- btl(d, "object_a", "object_b", response = "response",
             judge = "judge", thresholds = "pc")
  judges <- unique(fit$comparisons$judge)
  group <- setNames(rep(c("A", "B"), length.out = length(judges)), judges)
  da <- btl_dif(fit, group, objects = "O3", min_n = 10)
  db <- suppressWarnings(dif_bootstrap(fit, da, B = 1, workers = 1,
                                       seed = 5123))
  expect_equal(fit$m, 3L)
  expect_equal(db$B_used, 1L)

  h <- simulate_btl(
    5, 30, reps_per_pair = 12,
    dependence = list(exposure = .2, carry_over = .1), seed = 5124)
  hf <- btl(h, "object_a", "object_b", winner = "winner", judge = "judge",
            order = "order", position = TRUE)
  judges <- unique(hf$comparisons$judge)
  group <- setNames(rep(c("A", "B"), length.out = length(judges)), judges)
  hd <- btl_dif(hf, group, objects = "O3", min_n = 10)
  hb <- suppressWarnings(dif_bootstrap(hf, hd, B = 1, workers = 1,
                                       seed = 5125))
  expect_setequal(hf$dependence$effect,
                  c("exposure", "carry_over", "position"))
  expect_equal(hb$B_used, 1L)
})

test_that("DIF bootstrap travels only with its DIF analysis and active fit", {
  fit <- make_dif_boot_fit(seed = 5105, N = 120)
  da <- dif_anova(fit, n_groups = 2)
  db <- suppressWarnings(dif_bootstrap(fit, da, B = 1, workers = 1,
                                       seed = 8))

  expect_error(report_html(fit, tempfile(fileext = ".html"),
                           dif_bootstrap = db), "must accompany")
  html <- tempfile(fileext = ".html")
  on.exit(unlink(html), add = TRUE)
  expect_identical(suppressWarnings(report_html(
    fit, html, dif = da, dif_bootstrap = db)), html)
  txt <- paste(readLines(html, warn = FALSE), collapse = "\n")
  expect_match(txt, "Bootstrap sensitivity analysis", fixed = TRUE)

  out <- tempfile("dif-bootstrap-export-")
  on.exit(unlink(out, recursive = TRUE), add = TRUE)
  files <- suppressWarnings(save_outputs(
    fit, out, formats = "png", item_plots = FALSE,
    dif = da, dif_bootstrap = db, dpi = 72))
  expect_true(any(grepl("dif_conditional_bootstrap[.]csv$", files)))
  expect_true(any(grepl("dif_conditional_bootstrap_accounting[.]csv$", files)))

  project <- .seal_app_project(list(
    format = "rasch-shiny-project", schema = 2L,
    data = as.data.frame(fit$X), model_type = "rasch", base_fit = fit,
    rasch_steps = list(), btl_steps = list(),
    settings = list(dif_effects = da$effects, dif_alpha = da$alpha),
    results = list(dif = da,
      dif_bootstrap = list(db = db, B = 1L, seed = 8L))))
  expect_no_error(.validate_app_project(project))

  old_db <- db
  old_db$algorithm <- NULL
  unsigned_db <- unclass(old_db)
  unsigned_db$result_signature <- NULL
  old_db$result_signature <- .fit_boot_md5(unsigned_db)
  old_project <- project
  old_project$results$dif_bootstrap$db <- old_db
  old_project <- .seal_app_project(old_project)
  old_path <- tempfile(fileext = ".rasch")
  on.exit(unlink(old_path), add = TRUE)
  saveRDS(old_project, old_path)
  expect_warning(restored_old <- .read_app_project(old_path),
                 "earlier raw-F marginal")
  expect_null(restored_old$results$dif_bootstrap)
  expect_match(attr(restored_old, "rasch_project_legacy_dropped"),
               "raw-F marginal")
  expect_no_error(.validate_app_project(restored_old))

  missing_primary <- project
  missing_primary$results$dif <- NULL
  missing_primary <- .seal_app_project(missing_primary)
  expect_error(.validate_app_project(missing_primary),
               "no accompanying primary DIF")

  other_dif <- dif_anova(fit, effects = "factorial", alpha = .1)
  mismatched <- project
  mismatched$results$dif <- other_dif
  mismatched$settings <- list(dif_effects = other_dif$effects,
                              dif_alpha = other_dif$alpha)
  mismatched <- .seal_app_project(mismatched)
  expect_error(.validate_app_project(mismatched),
               "different DIF analysis")

  wrong_settings <- project
  wrong_settings$settings$dif_alpha <- .1
  wrong_settings <- .seal_app_project(wrong_settings)
  expect_error(.validate_app_project(wrong_settings),
               "does not match the restored DIF settings")

  other <- make_dif_boot_fit(seed = 5106, N = 120)
  primary_only <- project
  primary_only$results$dif_bootstrap <- NULL
  primary_only$base_fit <- other
  primary_only <- .seal_app_project(primary_only)
  expect_error(.validate_app_project(primary_only), "primary DIF analysis")

  project$base_fit <- other
  project <- .seal_app_project(project)
  expect_error(.validate_app_project(project),
               "does not belong to the active fit")
})

test_that("DIF bootstrap integrity is checked before restore or export", {
  fit <- make_dif_boot_fit(seed = 5130, N = 100)
  da <- dif_anova(fit, n_groups = 2)
  db <- suppressWarnings(dif_bootstrap(fit, da, B = 2, workers = 1,
                                       seed = 5131))

  changed_da <- da
  changed_da$terms$p_adj[1L] <- -1
  expect_error(rasch:::.validate_dif_result(changed_da, fit),
               "incomplete or internally inconsistent")
  changed_da <- da
  changed_da$summary$p_uniform_adj[1L] <- .777
  expect_error(rasch:::.validate_dif_result(changed_da, fit),
               "incomplete or internally inconsistent")

  broken <- list()
  broken$missing_accounting <- db
  broken$missing_accounting[c(
    "B", "B_used", "B_failed", "B_nonconverged", "B_errors",
    "family_n")] <- NULL
  broken$bad_accounting <- db
  broken$bad_accounting$B_errors <- broken$bad_accounting$B_errors + 1L
  broken$bad_replicate <- db
  colnames(broken$bad_replicate$replicates$F)[1] <- "wrong family member"
  broken$bad_minimum <- db
  broken$bad_minimum$replicates$min_p[1] <- 1
  broken$bad_terms <- db
  broken$bad_terms$terms$p_boot[1] <- .123
  broken$bad_summary <- db
  broken$bad_summary$summary$p_uniform_boot_adj[1] <- .123
  broken$bad_kind <- db
  broken$bad_kind$model_kind <- "btl"
  broken$bad_adjustment <- db
  broken$bad_adjustment$adjustment <- "Holm"
  broken$old_algorithm <- db
  broken$old_algorithm$algorithm <- NULL
  broken$bad_excluded_flag <- db
  excluded <- da$term_ids %in% c("Residuals", "ci")
  broken$bad_excluded_flag$terms$significant_boot[which(excluded)[1L]] <- NA

  for (nm in names(broken)) {
    x <- broken[[nm]]
    expect_error(rasch:::.validate_dif_bootstrap(x, fit, da),
                 "incomplete or internally inconsistent", info = nm)
  }

  out <- tempfile("bad-dif-bootstrap-export-")
  expect_false(dir.exists(out))
  expect_error(save_outputs(
    fit, out, formats = "png", item_plots = FALSE, dif = da,
    dif_bootstrap = broken$missing_accounting, dpi = 72),
    "incomplete or internally inconsistent")
  expect_false(dir.exists(out))

  report <- tempfile(fileext = ".html")
  expect_error(report_html(
    fit, report, dif = da,
    dif_bootstrap = broken$missing_accounting),
    "incomplete or internally inconsistent")
  expect_false(file.exists(report))

  project <- .seal_app_project(list(
    format = "rasch-shiny-project", schema = 2L,
    data = as.data.frame(fit$X), model_type = "rasch", base_fit = fit,
    rasch_steps = list(), btl_steps = list(),
    settings = list(dif_effects = da$effects, dif_alpha = da$alpha),
    results = list(dif = da, dif_bootstrap = list(
      db = broken$missing_accounting, B = 2L, seed = 5131L))))
  expect_error(.validate_app_project(project),
               "incomplete or internally inconsistent")
  project_file <- tempfile(fileext = ".rasch")
  expect_error(.save_app_project(project, project_file),
               "incomplete or internally inconsistent")
  expect_false(file.exists(project_file))
})
