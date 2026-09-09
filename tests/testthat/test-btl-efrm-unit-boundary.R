boundary_frame_data <- function(b_within_wins = 12L, reversed = TRUE) {
  block <- function(a, b, wins) data.frame(
    a = a, b = b, win = c(rep(a, wins), rep(b, 100L - wins)))
  d <- do.call(rbind, list(
    block("A1", "A2", 12L), block("B1", "B2", b_within_wins),
    block("A1", "B1", if (reversed) 12L else 50L),
    block("A2", "B1", if (reversed) 50L else 88L),
    block("A1", "B2", if (reversed) 50L else 12L),
    block("A2", "B2", if (reversed) 88L else 50L)))
  d$judge <- paste0("J", rep(1:20, length.out = nrow(d)))
  d$panel <- "P"
  d
}

fit_boundary_frame <- function(d, method = "conditional") {
  btl_efrm(d, "a", "b", "win", "judge", "panel",
    object_sets = list(A = c("A1", "A2"), B = c("B1", "B2")),
    se_method = method, boot_reps = 30L, workers = 1L, seed = 812)
}

test_that("zero linking units cannot be replaced by unit one", {
  d <- boundary_frame_data()
  for (method in c("conditional", "bootstrap", "judge_bootstrap"))
    expect_error(fit_boundary_frame(d, method), "zero boundary: B")
})

test_that("genuinely flat within-set locations retain the placement convention", {
  f <- fit_boundary_frame(boundary_frame_data(b_within_wins = 50L))
  expect_true(f$converged)
  expect_true(is.na(f$alpha_table$alpha[f$alpha_table$set == "B"]))
  expect_equal(f$objects$beta_set[f$objects$set == "B"], c(0, 0))
  expect_true(all(is.finite(f$objects$location)))
  expect_match(paste(f$notes, collapse = " "), "unit\\(s\\) unidentified")
})

test_that("both bootstrap paths reject zero-unit refits", {
  original <- .btlef_stage2
  d <- boundary_frame_data()
  cross <- substr(d$a, 1L, 1L) != substr(d$b, 1L, 1L)
  dc <- d[cross, ]
  rejected <- original(dc$a, dc$b, as.integer(dc$win == dc$a),
    rep(1, nrow(dc)), substr(dc$a, 1L, 1L), substr(dc$b, 1L, 1L),
    c(A1 = -1, A2 = 1, B1 = -1, B2 = 1), c("A", "B"), 300L, 1e-8)
  expect_true(rejected$alpha_boundary[["B"]])
  expect_false(rejected$rank_ok)
  expect_false(rejected$converged)
  calls <- 0L
  testthat::local_mocked_bindings(.btlef_stage2 = function(...) {
    calls <<- calls + 1L
    if (calls == 1L) original(...) else rejected
  }, .package = "rasch")
  for (method in c("bootstrap", "judge_bootstrap")) {
    calls <- 0L
    expect_error(fit_boundary_frame(boundary_frame_data(reversed = FALSE),
                                    method), "only 0 usable fits")
    expect_equal(calls, 31L)
  }
})
