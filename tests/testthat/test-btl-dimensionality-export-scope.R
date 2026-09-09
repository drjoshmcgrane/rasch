test_that("CJ archives retain the dimensionality reference assumption", {
  d <- simulate_btl(6, 20, reps_per_pair = 10, seed = 9874)
  fit <- btl(d, "object_a", "object_b", "winner", judge = "judge")
  out <- tempfile("btl-reference-scope-")
  on.exit(unlink(out, recursive = TRUE), add = TRUE)
  local_mocked_bindings(.rr_save_plot = function(...) character(0),
                         .package = "rasch")
  for (independent in c(FALSE, TRUE)) {
    destination <- file.path(out, if (independent) "conditional" else
                               "descriptive")
    result <- btl_dimensionality(fit, reps = 20L, seed = 4,
      independent_comparisons = independent)
    save_outputs(fit, destination, formats = "png", item_plots = FALSE,
                   dimensionality = result)
    tab <- utils::read.csv(file.path(destination, "tables",
      "residual_dimensionality_reference.csv"))
    expect_identical(tab$independent_comparisons, independent)
    expect_identical(tab$inference_available, independent)
    expect_identical(is.finite(tab$p), independent)
    expect_match(tab$note, if (independent) "not cluster-robust" else
      "inference is withheld")
  }
})
