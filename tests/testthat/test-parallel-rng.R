.rng_socket_support <- function() {
  skip_on_cran()
  skip_if_not(.rasch_namespace_is_installed(),
              "parallel RNG integration needs the installed current package")
  probe <- try(parallel::makePSOCKcluster(2L), silent = TRUE)
  skip_if(inherits(probe, "try-error"), "local socket clusters unavailable")
  parallel::stopCluster(probe)
}

test_that("socket bootstrap jobs inherit all RNG kinds", {
  .rng_socket_support()
  old_kind <- RNGkind(); old_seed <- .sim_seed_capture()
  on.exit({
    suppressWarnings(do.call(RNGkind, as.list(old_kind)))
    .sim_seed_restore(old_seed)
  }, add = TRUE)
  # Each job seeds itself, exactly as the fitted-outcome and conditional
  # bootstrap generators do. Odd normal-vector lengths exercise Box-Muller
  # as well as the uniform and discrete samplers.
  job <- function(i) {
    set.seed(713L + i)
    list(pid = Sys.getpid(), kind = RNGkind(),
         draw = c(runif(5), rnorm(5), sample.int(1000, 5)))
  }
  for (kind in list(c("Mersenne-Twister", "Inversion", "Rejection"),
                   c("L'Ecuyer-CMRG", "Inversion", "Rejection"),
                   c("Mersenne-Twister", "Box-Muller", "Rejection"),
                   c("Mersenne-Twister", "Inversion", "Rounding"))) {
    suppressWarnings(do.call(RNGkind, as.list(kind)))
    serial <- .rasch_boot_apply(4L, job, workers = 1L)
    parallel <- suppressWarnings(.rasch_boot_apply(4L, job, workers = 2L))
    # A serial fallback must not make this integration test falsely pass.
    expect_true(all(vapply(parallel, function(x) x$pid != Sys.getpid(), TRUE)))
    expect_identical(lapply(parallel, `[[`, "kind"),
                     lapply(serial, `[[`, "kind"))
    expect_identical(lapply(parallel, `[[`, "draw"),
                     lapply(serial, `[[`, "draw"))
  }
})

test_that("socket bootstrap workers load the coordinator's installation", {
  .rng_socket_support()
  job <- function(i) list(
    pid = Sys.getpid(),
    path = normalizePath(getNamespaceInfo(asNamespace("rasch"), "path")),
    conditional_generator = exists(".fit_gen_conditional",
      asNamespace("rasch"), inherits = FALSE))
  out <- .rasch_boot_apply(2L, job, workers = 2L)
  expected <- normalizePath(getNamespaceInfo(asNamespace("rasch"), "path"))
  expect_true(all(vapply(out, function(z) z$pid != Sys.getpid(), TRUE)))
  expect_true(all(vapply(out, function(z) identical(z$path, expected), TRUE)))
  expect_true(all(vapply(out, `[[`, TRUE, "conditional_generator")))
})

test_that("seeded item-fit bootstraps agree under non-default RNGs", {
  .rng_socket_support()
  old_workers <- options(rasch.max_workers = 2L,
                         rasch.efrm.max_workers = NULL)
  on.exit(options(old_workers), add = TRUE)
  skip_if(.rasch_available_workers() < 2L,
          "system limits allow only serial public bootstraps")
  old_kind <- RNGkind(); old_seed <- .sim_seed_capture()
  on.exit({
    suppressWarnings(do.call(RNGkind, as.list(old_kind)))
    .sim_seed_restore(old_seed)
  }, add = TRUE)
  set.seed(419)
  theta <- rnorm(200)
  X <- vapply(seq(-1, 1, length.out = 5), function(d)
    rbinom(200, 1, plogis(theta - d)), numeric(200))
  colnames(X) <- paste0("I", 1:5)
  fit <- rasch(X)
  RNGkind("L'Ecuyer-CMRG")
  set.seed(821)
  before <- .sim_seed_capture()
  serial <- suppressWarnings(fit_bootstrap(fit, B = 5, workers = 1, seed = 713))
  expect_identical(.sim_seed_capture(), before)
  parallel <- suppressWarnings(fit_bootstrap(fit, B = 5, workers = 2, seed = 713))
  expect_identical(.sim_seed_capture(), before)
  expect_identical(RNGkind()[1L], "L'Ecuyer-CMRG")
  expect_identical(parallel$B_used, 5L)
  expect_identical(parallel$replicates, serial$replicates)
  expect_identical(parallel$items, serial$items)
  expect_identical(parallel$persons, serial$persons)
})
