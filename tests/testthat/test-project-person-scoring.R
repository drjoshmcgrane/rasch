.scoring_project <- function(fit) .seal_app_project(list(
  format = "rasch-shiny-project", schema = 2L, package_version = "test",
  created = "2026-09-09", data = as.data.frame(fit$X),
  model_type = .app_fit_family(fit), base_fit = fit,
  rasch_steps = list(), btl_steps = list(), rcode = "fit <- rasch(data)",
  kept_fits = list(), kept_fit_code = list(),
  settings = list(model_type = .app_fit_family(fit), item_cols = colnames(fit$X)),
  resources = list(), simulation = list(), results = list()))

.first_root_wle <- function(tau_list, disc = 1) {
  max_score <- sum(lengths(tau_list))
  theta <- se <- setNames(rep(NA_real_, max_score + 1L), 0:max_score)
  for (r in 0:max_score) {
    g <- function(th) {
      mo <- lapply(tau_list, item_moments, theta = th, disc = disc)
      v <- sum(vapply(mo, `[[`, 0, "V"))
      r - sum(vapply(mo, `[[`, 0, "E")) +
        sum(vapply(mo, `[[`, 0, "mu3")) / (2 * v)
    }
    root <- uniroot(g, .person_root_interval(tau_list, disc),
                    tol = 1e-9 / disc)$root
    theta[as.character(r)] <- root
    se[as.character(r)] <- 1 / (disc * sqrt(sum(vapply(tau_list,
      function(tt) item_moments(root, tt, disc = disc)$V, 0))))
  }
  list(theta = theta, se = se)
}

test_that("affected old scoring is refused throughout saved projects", {
  set.seed(1)
  X <- matrix(sample(0:1, 500 * 7, replace = TRUE), 500, 7,
              dimnames = list(NULL, paste0("I", 1:7)))
  anchors <- data.frame(item = c("I1", "I2", "I3", "I5", "I7"),
                         k = 1L, tau = c(-6, -6, -6, 6, 6))
  current <- rasch(X, model = "PCM", anchors = anchors)
  expect_identical(current$person_scoring_algorithm, "max-wle-1")
  old <- local({
    local_mocked_bindings(person_wle = .first_root_wle)
    .assemble_fit(current$model, current$X, current$est, current$person$id,
                  current$factors, current$n_groups, current$notes)
  })
  old$person_scoring_algorithm <- NULL
  # Every dependent field was assembled at the old roots, rather than
  # creating an inconsistent fixture by replacing person theta alone.
  expect_gt(max(abs(old$person$theta - current$person$theta)), 6)
  expect_gt(abs(old$psi$PSI - current$psi$PSI), .1)
  project <- .scoring_project(current)
  entry <- list(type = "test", label = "Prior analysis", details = list(),
                code = "fit <- rasch(data)", created = "2026-09-09", fit = old)
  path <- tempfile(fileext = ".rasch")
  on.exit(unlink(path))
  for (where in c("base", "history", "kept")) {
    p <- project
    if (where == "base") p$base_fit <- old
    if (where == "history") p$rasch_steps <- list(entry)
    if (where == "kept") p$kept_fits <- list(previous = old)
    p <- .seal_app_project(p)
    saveRDS(p, path)
    before <- tools::md5sum(path)
    expect_error(.read_app_project(path), "superseded or unverified person scoring")
    expect_identical(tools::md5sum(path), before)
  }
  p <- .scoring_project(old); p$schema <- 1L; p$binding <- NULL
  saveRDS(p, path)
  expect_error(.read_app_project(path), "refit this analysis")
})

test_that("compatible unversioned scores remain readable without resealing results", {
  fit <- rasch(simulate_rasch(100, 6, seed = 8621))
  fit$person_scoring_algorithm <- NULL
  p <- .scoring_project(fit)
  path <- tempfile(fileext = ".rasch")
  on.exit(unlink(path))
  saveRDS(p, path)
  expect_no_warning(actual <- .read_app_project(path))
  expect_identical(actual, p)
  expect_no_error(.validate_app_person_scoring(fit, "fit"))
  fit$person_scoring_algorithm <- "unknown"
  expect_error(.validate_app_person_scoring(fit, "fit"), "unverified person scoring")

  # Authenticate before diagnosing changed fields in an older bundle.
  p$base_fit$person$theta[1L] <- p$base_fit$person$theta[1L] + 1
  saveRDS(p, path)
  expect_error(.read_app_project(path), "changed since they were saved")
})

test_that("old-scoring checks include unused conversion scores and unequal frame units", {
  fit <- rasch(simulate_rasch(100, 6, seed = 8621))
  fit$person_scoring_algorithm <- NULL
  fit$score_table$theta[1L] <- fit$score_table$theta[1L] + 1
  expect_error(.validate_app_person_scoring(fit, "fit"), "refit")

  d <- simulate_efrm(n_per_group = 70, items_per_set = 5,
                     n_groups = 2, n_sets = 1, seed = 8651)
  tr <- attr(d, "truth")
  ef <- rasch_efrm(d, item_sets = tr$item_sets, groups = "group", id = "id",
                   boot_reps = 0, workers = 1)
  expect_identical(ef$person_scoring_algorithm, "max-wle-1")
  expect_gt(length(unique(ef$disc)), 1L)
  ef$person_scoring_algorithm <- NULL
  expect_no_error(.validate_app_person_scoring(ef, "EFRM fit"))
  ef$person$theta[1L] <- ef$person$theta[1L] + 1
  expect_error(.validate_app_person_scoring(ef, "EFRM fit"), "refit")
})

test_that("saved kurtosis calibrations require the corrected polynomial", {
  d <- simulate_rasch(250, 6, model = "PCM", n_categories = 5, seed = 8652)
  fit <- rasch(d, pc_components = 4)
  expect_identical(fit$est$pc_algorithm, "guttman-four-1")
  expect_no_error(.validate_app_person_scoring(fit, "PC fit"))
  fit$est$pc_algorithm <- NULL
  expect_error(.validate_app_person_scoring(fit, "PC fit"), "corrected principal-component kurtosis")
  # Lower-order PC fits never used the changed polynomial.
  reduced <- rasch(d, pc_components = 3)
  reduced$est$pc_algorithm <- NULL
  expect_no_error(.validate_app_person_scoring(reduced, "PC fit"))
})
